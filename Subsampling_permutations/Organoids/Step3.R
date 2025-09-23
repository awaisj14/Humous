#!/usr/bin/env Rscript

# run_ordinals_serial_updated.R
# Nested ordinal regression pipeline (AGE & DIFF) + landscape ordinals
# - Auto-discovers IntO/IntHO_* files
# - Removes IntH2 cells
# - Label hygiene for age_ek/diff_ek
# - Reproducible sampling
# - Robust checks for integrated@scale.data
# - Mean coefficients across permutations + comparison (optional if refs exist)
# - Landscape overlap metrics

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(plyr)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(DescTools)
})

# 0) Utilities (expects: training_ordi_full, custom_red_and_pred, costM, ordi_normalize)
source("/home/users/j/javed/Humous_reviews/lib_misc.R")

# ---------- PARAMETERS ----------
set.seed(123)            # reproducible sampling
input_dir  <- "IntO"     # where IntHO_*.rds live
models_dir <- "models"
land_dir   <- "lands"
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(land_dir,   recursive = TRUE, showWarnings = FALSE)

# Label maps (ensure these match your metadata values)
age_map  <- c("1-2 m"=1, "2-3 m"=2, "3 m"=3, "3.5-5 m"=4, "6-7 m"=5)
diff_map <- c("RG"=1, "IPC"=2, "N"=3)

# Modeling params
n_lvl2   <- 5      # number of level-2 fits per permutation (averaged)
train_n  <- 50     # target cells per age×diff group for training; auto-downscales to smallest group
max_iter <- 1000
ngenesselect <- 25
lambda_full  <- 10
lambda_red   <- 10
epsilon      <- 1e-7
nfolds       <- 20

# ---------- HELPERS ----------
# Safe sample without replacement; if group smaller than k, take all
sample_safe <- function(v, k) {
  v <- as.character(v)
  if (length(v) <= k) v else sample(v, k, replace = FALSE)
}

# Extract zero-padded permutation ID from filename (fallback to seq index)
get_perm_label <- function(fpath, idx) {
  # matches ...IntHO_full.list_f_perm_0X.rds (X can be 1 or 2 digits)
  m <- regexec("IntHO_full\\.list_f_perm_0(\\d+)\\.rds$", basename(fpath))
  g <- regmatches(basename(fpath), m)
  if (length(g) == 1 && length(g[[1]]) == 2) {
    sprintf("%02d", as.integer(g[[1]][2]))
  } else {
    sprintf("%02d", idx) # fallback: running index
  }
}

# Ensure integrated@scale.data is present (compute if missing/empty)
ensure_scaled <- function(seu) {
  stopifnot("integrated" %in% Assays(seu))
  sd_mat <- seu@assays$integrated@scale.data
  if (is.null(sd_mat) || nrow(sd_mat) == 0 || ncol(sd_mat) == 0) {
    DefaultAssay(seu) <- "integrated"
    seu <- ScaleData(seu, assay = "integrated", vars.to.regress = "nFeature_RNA", verbose = FALSE)
  }
  seu
}

# Label hygiene: trim, normalize dashes/spaces
clean_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\u2013|\u2014", "-", x)  # en/em dash -> hyphen
  x <- gsub("\\s+", " ", x)           # collapse multiple spaces
  x
}

# ---------- DISCOVER INPUT FILES ----------
# Expecting files like: IntO/IntHO_full.list_f_perm_0X.rds
int_files <- list.files(input_dir, pattern = "^IntHO_full\\.list_f_perm_0\\d+\\.rds$", full.names = TRUE)
if (length(int_files) == 0L) {
  stop("No input files found in '", input_dir, "'. Expected pattern: IntHO_full.list_f_perm_0X.rds")
}

# Natural sort by numeric suffix
perm_nums <- as.integer(sub("^.*_perm_0(\\d+)\\.rds$", "\\1", basename(int_files)))
ord <- order(perm_nums, na.last = TRUE)
int_files <- int_files[ord]
perm_nums <- perm_nums[ord]

cat("Discovered input files:\n")
print(basename(int_files))

n_perms <- length(int_files)

# ---------- STORAGE FOR COEFS ----------
all_fullAGE <- vector("list", n_perms)

# ---------- MAIN LOOP ----------
for (i in seq_len(n_perms)) {
  perm_label <- get_perm_label(int_files[i], i)  # "01", "02", ...
  cat(sprintf("\n=== Permutation %s ===\n", perm_label))
  
  # Load integrated object
  IntH <- readRDS(int_files[i])
  
  # Remove IntH2 cells if column present; otherwise warn and continue
  if ("dataset_proto" %in% colnames(IntH@meta.data)) {
    IntH <- subset(IntH, subset = dataset_proto != "IntH2")
  } else {
    warning(sprintf("Perm %s: 'dataset_proto' not in metadata; skipping IntH2 removal", perm_label))
  }
  
  # Ensure integrated assay scaled
  if (!("integrated" %in% Assays(IntH))) stop("Perm ", perm_label, ": missing 'integrated' assay")
  IntH <- ensure_scaled(IntH)
  sd_mat <- IntH@assays$integrated@scale.data
  if (is.null(sd_mat) || nrow(sd_mat) == 0 || ncol(sd_mat) == 0) stop("Perm ", perm_label, ": empty integrated@scale.data after ScaleData")
  
  # Clean labels and map to integers
  IntH$age_ek  <- clean_label(IntH$age_ek)
  IntH$diff_ek <- clean_label(IntH$diff_ek)
  IntH$y_age   <- as.integer(age_map[ IntH$age_ek ])
  IntH$y_diff  <- as.integer(diff_map[ IntH$diff_ek ])
  
  if (anyNA(IntH$y_age)) {
    miss <- unique(IntH$age_ek[is.na(IntH$y_age)])
    stop("Perm ", perm_label, ": unmapped age_ek levels: ", paste(miss, collapse = ", "))
  }
  if (anyNA(IntH$y_diff)) {
    miss <- unique(IntH$diff_ek[is.na(IntH$y_diff)])
    stop("Perm ", perm_label, ": unmapped diff_ek levels: ", paste(miss, collapse = ", "))
  }
  
  # Stratified train set by age×diff
  IntH$group <- paste0(IntH$y_age, "_", IntH$y_diff)
  grp_counts <- table(IntH$group)
  sample_n1  <- min(train_n, as.integer(min(grp_counts)))
  if (sample_n1 < train_n) warning(sprintf("Perm %s: train_n reduced to %d due to small groups", perm_label, sample_n1))
  
  train_cells1 <- unlist(lapply(names(grp_counts), function(g) {
    cs <- colnames(IntH)[IntH$group == g]
    sample_safe(cs, sample_n1)
  }), use.names = FALSE)
  IntH_train <- subset(IntH, cells = train_cells1)
  if (length(colnames(IntH_train)) != length(train_cells1)) stop("Perm ", perm_label, ": mismatch in IntH_train cell count")
  
  # ----- Level-2 fits for AGE (average coefficients) -----
  lvl2_coefs <- vector("list", n_lvl2)
  for (j in seq_len(n_lvl2)) {
    # resample balanced cells for each lvl2 fit
    cells_j <- unlist(lapply(names(grp_counts), function(g) {
      cs <- colnames(IntH)[IntH$group == g]
      sample_safe(cs, sample_n1)
    }), use.names = FALSE)
    tmp <- subset(IntH, cells = cells_j)
    
    # Ensure scale.data still aligned
    X <- t(tmp@assays$integrated@scale.data)
    y <- tmp$y_age
    
    m2 <- training_ordi_full(
      datamatrix   = X,
      target       = y,
      cost         = costM(X = X, y = y),
      lambda_full  = lambda_full,
      epsilon_full = epsilon,
      maxiter      = max_iter
    )
    if (is.list(m2) && !is.null(m2$coefficients))      c2 <- m2$coefficients
    else if (is.atomic(m2))                             c2 <- m2
    else                                                c2 <- suppressWarnings(coef(m2))
    
    # Align coefficients to genes
    gene_names <- rownames(sd_mat)
    if (is.null(names(c2))) {
      names(c2) <- gene_names
    } else {
      # reorder/intersect if needed
      keep <- intersect(names(c2), gene_names)
      c2 <- c2[keep]
      # pad missing with 0 to keep vector lengths consistent
      if (length(keep) < length(gene_names)) {
        z <- setNames(numeric(length(gene_names)), gene_names)
        z[names(c2)] <- c2
        c2 <- z
      }
    }
    lvl2_coefs[[j]] <- c2
  }
  mean_coef         <- Reduce(`+`, lvl2_coefs) / n_lvl2
  all_fullAGE[[i]]  <- mean_coef
  write_csv(
    data.frame(gene = names(mean_coef), coeff = as.numeric(mean_coef)),
    file.path(models_dir, sprintf("mean_fullAGE_coeffs_perm_%s.csv", perm_label))
  )
  
  # ----- Reduced AGE model & predictions -----
  redA <- custom_red_and_pred(
    fullmodel    = lvl2_coefs[[1]],
    xtrain       = t(IntH_train@assays$integrated@scale.data),
    target       = IntH_train$y_age,
    xtest        = t(IntH@assays$integrated@scale.data),
    ngenesselect = ngenesselect,
    lambda_red   = lambda_red,
    epsilon_red  = epsilon,
    maxiter      = max_iter,
    nfolds       = nfolds
  )
  predA <- redA$pred
  if (!is.numeric(predA)) stop("Perm ", perm_label, ": AGE pred not numeric")
  if (is.null(names(predA))) names(predA) <- colnames(IntH)
  IntH$ordi_age <- scales::rescale(predA[colnames(IntH)], to = c(0, 1))
  if (anyNA(IntH$ordi_age)) stop("Perm ", perm_label, ": NA in ordi_age")
  saveRDS(redA$redmodel, file = file.path(models_dir, sprintf("redAGE_perm_%s.rds", perm_label)))
  
  # ----- Reduced DIFF model & predictions -----
  Xtr <- t(IntH_train@assays$integrated@scale.data)
  ytr <- IntH_train$y_diff
  fullD <- training_ordi_full(
    datamatrix   = Xtr,
    target       = ytr,
    cost         = costM(X = Xtr, y = ytr),
    lambda_full  = lambda_full,
    epsilon_full = epsilon,
    maxiter      = max_iter
  )
  redD <- custom_red_and_pred(
    fullmodel    = fullD,
    xtrain       = t(IntH_train@assays$integrated@scale.data),
    target       = IntH_train$y_diff,
    xtest        = t(IntH@assays$integrated@scale.data),
    ngenesselect = ngenesselect,
    lambda_red   = lambda_red,
    epsilon_red  = epsilon,
    maxiter      = max_iter,
    nfolds       = nfolds
  )
  predD <- redD$pred
  if (!is.numeric(predD)) stop("Perm ", perm_label, ": DIFF pred not numeric")
  if (is.null(names(predD))) names(predD) <- colnames(IntH)
  IntH$ordi_diff <- scales::rescale(predD[colnames(IntH)], to = c(0, 1))
  if (anyNA(IntH$ordi_diff)) stop("Perm ", perm_label, ": NA in ordi_diff")
  saveRDS(redD$redmodel, file = file.path(models_dir, sprintf("redDIFF_perm_%s.rds", perm_label)))
  
  # ----- Normalized ordinals + save landscape -----
  IntH$ordi_diff_norm <- ordi_normalize(
    IntH$ordi_diff, IntH$y_diff, IntH$y_age,
    limits1 = c(0, 0.24), limits2 = c(0.21, 0.47), limits3 = c(0.4, 1),
    applyWinsor = TRUE
  )
  IntH$ordi_age_norm <- ordi_normalize(
    IntH$ordi_age, IntH$y_age, IntH$y_diff,
    applyWinsor = TRUE
  )
  
  saveRDS(IntH, file = file.path(land_dir, sprintf("IntH_L_perm_%s.rds", perm_label)))
  cat(sprintf("Permutation %s landscape saved\n", perm_label))
}

# ---------- Aggregate fullAGE coefficients ----------
stopifnot(length(all_fullAGE) > 0)
gene_ref <- names(all_fullAGE[[1]])
coef_df  <- data.frame(gene = gene_ref, stringsAsFactors = FALSE)

for (i in seq_len(n_perms)) {
  lab <- get_perm_label(int_files[i], i)
  coef_df[[sprintf("perm%s", lab)]] <- as.numeric(all_fullAGE[[i]][gene_ref])
}
perm_cols <- grep("^perm", names(coef_df), value = TRUE)
coef_df$mean_perm <- rowMeans(coef_df[perm_cols], na.rm = TRUE)

write_csv(coef_df, file = file.path(models_dir, "coef_fullAGE_perms.csv"))

# ---------- Compare with reference (optional; only if files found) ----------
ref_IntH_path <- "/home/users/j/javed/Humous_reviews/Org/IntO_OG/IntO.rds"
ref_mod_path  <- "/home/users/j/javed/Humous_reviews/Org/ordiO_OG/model_fullAGE_O.rds"

if (file.exists(ref_IntH_path) && file.exists(ref_mod_path)) {
  ref_IntH <- readRDS(ref_IntH_path)
  ref_mod  <- readRDS(ref_mod_path)
  
  if (is.list(ref_mod) && !is.null(ref_mod$coefficients)) {
    coef_ref <- ref_mod$coefficients
  } else if (is.atomic(ref_mod)) {
    coef_ref <- ref_mod
  } else {
    coef_ref <- suppressWarnings(coef(ref_mod))
  }
  
  # Align names to ref genes
  ref_genes <- rownames(ref_IntH@assays$integrated@scale.data)
  if (is.null(names(coef_ref))) names(coef_ref) <- ref_genes
  
  ref_df  <- data.frame(gene = names(coef_ref), ref = as.numeric(coef_ref), stringsAsFactors = FALSE)
  plot_df <- inner_join(ref_df, coef_df[, c("gene", "mean_perm")], by = "gene")
  
  p <- ggplot(plot_df, aes(x = ref, y = mean_perm)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(
      x = "Original fullAGE coeff",
      y = "Mean perm fullAGE coeff",
      title = "Original vs Mean Perm Coeffs"
    ) +
    annotate(
      "text",
      x = min(plot_df$ref, na.rm = TRUE),
      y = max(plot_df$mean_perm, na.rm = TRUE),
      label = sprintf("rho=%.2f", cor(plot_df$ref, plot_df$mean_perm, use = "complete.obs")),
      hjust = 0, vjust = 1
    )
  
  pdf(file.path(models_dir, "coef_fullAGE_comparison.pdf"), width = 6, height = 6)
  print(p); dev.off()
} else {
  message("Reference files not found; skipping comparison plot.")
}

# ---------- Landscape overlap metrics ----------
perm_files <- list.files(land_dir, pattern = "^IntH_L_perm_\\d+\\.rds$|^IntH_L_perm_\\d{2}\\.rds$", full.names = TRUE)
if (length(perm_files) >= 2) {
  cell_lists <- lapply(perm_files, function(f) colnames(readRDS(f)))
  n <- length(cell_lists)
  jaccard_mat <- matrix(NA_real_, n, n, dimnames = list(basename(perm_files), basename(perm_files)))
  
  for (a in seq_len(n)) for (b in seq_len(n)) {
    inter_ab <- length(intersect(cell_lists[[a]], cell_lists[[b]]))
    union_ab <- length(union(cell_lists[[a]], cell_lists[[b]]))
    jaccard_mat[a, b] <- if (union_ab == 0) NA_real_ else inter_ab / union_ab
  }
  
  jac_vals <- jaccard_mat[lower.tri(jaccard_mat)]
  metrics <- data.frame(
    min_jaccard  = min(jac_vals, na.rm = TRUE),
    mean_jaccard = mean(jac_vals, na.rm = TRUE),
    max_jaccard  = max(jac_vals, na.rm = TRUE),
    common_cells = length(Reduce(intersect, cell_lists))
  )
  print(metrics)
  write_csv(metrics, file = file.path(land_dir, "landscape_overlap_metrics.csv"))
} else {
  message("Not enough landscapes to compute overlap metrics (need >= 2).")
}

cat("\nAll permutations and landscape ordinals done.\nResults in:", models_dir, "and", land_dir, "\n")
