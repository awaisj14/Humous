#!/usr/bin/env Rscript

# run_ordinals_serial.R (updated)
# Nested ordinal regression pipeline (AGE & DIFF) + landscape ordinals
# Level1 permutations + Level2 subsampling for robust fullAGE estimation
# Now: full landscape for all perms + single subset called landscapeSubset
# Plus: metrics summary of overlap among IntH_L_perm_* objects

# 1) Load utilities & libraries
source("lib_misc.R")   # Provides training_ordi_full, custom_red_and_pred, costM, ordi_normalize
library(Seurat)
library(dplyr)         # use dplyr:: explicitly below
library(readr)         # use readr:: explicitly below
library(plyr)
library(ggplot2)
library(scales)
library(purrr)
library(DescTools)

# 2) Prepare output directories
models_dir <- "models"
land_dir   <- "lands"
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(land_dir,   recursive = TRUE, showWarnings = FALSE)

# 3) Parameters
n_perms  <- 10    # Number of Level1 permutations
n_lvl2   <- 5     # Number of Level2 subsamples per permutation
train_n  <- 150   # Desired training cells per (age × diff) group
land_n   <- 250   # Desired landscape training cells per group
max_iter <- 1000  # Max iterations for ordinal SVM

# 4) Storage for averaged fullAGE coefficients
all_fullAGE <- vector("list", n_perms)

# 5) Level1 permutation loop
for (i in seq_len(n_perms)) {
  cat(sprintf("\n=== Permutation %02d ===\n", i))
  
  # 5.1) Load integrated Seurat object
  int_file <- sprintf("IntH_perm_%02d.rds", i)
  if (!file.exists(int_file)) stop("Missing file: ", int_file)
  IntH <- readRDS(int_file)
  
  # 5.2) Sanity checks on assay
  if (!"integrated" %in% Assays(IntH)) stop("Perm ", i, ": missing 'integrated' assay")
  sd_mat <- IntH@assays$integrated@scale.data
  if (is.null(sd_mat) || nrow(sd_mat)==0 || ncol(sd_mat)==0) stop("Perm ", i, ": empty integrated@scale.data")
  
  # 5.3) Recode diff and bin ages
  IntH$diff_ek[IntH$diff_ek=="Neuron"] <- "N"
  IntH$age_ek <- plyr::revalue(
    as.character(IntH$age_ek),
    c("12"="12","15"="15-16","16"="15-16","17"="17-18","18"="17-18",
      "20"="20-21","21"="20-21","23"="23-24","24"="23-24")
  )
  
  # 5.4) Numeric encoding & check
  age_map  <- c("12"=1,"15-16"=2,"17-18"=3,"20-21"=4,"23-24"=5)
  diff_map <- c("RG"=1,"IPC"=2,"N"=3)
  IntH$y_age  <- as.integer(age_map[IntH$age_ek])
  IntH$y_diff <- as.integer(diff_map[IntH$diff_ek])
  if (anyNA(IntH$y_age) || anyNA(IntH$y_diff)) stop("Perm ", i, ": NA in y_age/y_diff labels")
  
  # 5.5) Balanced sampling for Level1 training
  IntH$group <- paste0(IntH$y_age, "_", IntH$y_diff)
  grp_counts  <- table(IntH$group)
  sample_n1   <- min(train_n, min(grp_counts))
  if (sample_n1 < train_n) warning(sprintf("Perm %02d: train_n reduced to %d", i, sample_n1))
  train_cells1 <- unlist(lapply(names(grp_counts), function(g) {
    cs <- colnames(IntH)[IntH$group == g]
    sample(cs, sample_n1)
  }))
  IntH_train <- subset(IntH, cells = train_cells1)
  if (length(colnames(IntH_train)) != length(train_cells1))
    stop("Perm ", i, ": mismatch in IntH_train cell count")
  
  # 5.6) Level2 subsampling & fullAGE estimation
  lvl2_coefs <- vector("list", n_lvl2)
  for (j in seq_len(n_lvl2)) {
    cells_j <- unlist(lapply(names(grp_counts), function(g) {
      cs <- colnames(IntH)[IntH$group == g]
      sample(cs, sample_n1)
    }))
    tmp <- subset(IntH, cells = cells_j)
    m2 <- training_ordi_full(
      datamatrix   = t(tmp@assays$integrated@scale.data),
      target       = tmp$y_age,
      cost         = costM(X = t(tmp@assays$integrated@scale.data), y = tmp$y_age),
      lambda_full  = 0.05, epsilon_full = 1e-7, maxiter = max_iter
    )
    # extract & name coefficients
    if      (is.list(m2)   && !is.null(m2$coefficients)) c2 <- m2$coefficients
    else if (is.atomic(m2))                           c2 <- m2
    else                                              c2 <- suppressWarnings(coef(m2))
    if (is.null(names(c2))) names(c2) <- rownames(sd_mat)
    lvl2_coefs[[j]] <- c2
  }
  mean_coef       <- Reduce(`+`, lvl2_coefs) / n_lvl2
  all_fullAGE[[i]] <- mean_coef
  readr::write_csv(
    data.frame(gene = names(mean_coef), coeff = mean_coef),
    file.path(models_dir, sprintf("mean_fullAGE_coeffs_perm_%02d.csv", i))
  )
  
  # 5.7) Reduced AGE modeling (use first Level2 model)
  redA <- custom_red_and_pred(
    fullmodel    = lvl2_coefs[[1]],
    xtrain       = t(IntH_train@assays$integrated@scale.data),
    target       = IntH_train$y_age,
    xtest        = t(IntH@assays$integrated@scale.data),
    ngenesselect = 25, lambda_red = 0.01, epsilon_red = 1e-7,
    maxiter      = max_iter, nfolds = 20
  )
  predA <- redA$pred
  if (!is.numeric(predA)) stop("Perm ", i, ": AGE pred not numeric")
  if (is.null(names(predA))) names(predA) <- colnames(IntH)
  if (!all(colnames(IntH) %in% names(predA))) stop("Perm ", i, ": missing AGE pred names")
  IntH$ordi_age <- scales::rescale(predA[colnames(IntH)], to = c(0,1))
  if (anyNA(IntH$ordi_age)) stop("Perm ", i, ": NA in ordi_age after rescaling")
  saveRDS(redA$redmodel, file = file.path(models_dir, sprintf("redAGE_perm_%02d.rds", i)))
  
  metricsAGE <- IntH@meta.data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(min = min(ordi_age),
                     max = max(ordi_age),
                     median = median(ordi_age),
                     sd = sd(ordi_age)) %>%
    as.data.frame()
  readr::write_csv(metricsAGE, file = file.path(models_dir, sprintf("metricsAGE_perm_%02d.csv", i)))
  
  # 5.8) Full & reduced DIFF modeling
  redD <- custom_red_and_pred(
    fullmodel    = training_ordi_full(
      datamatrix   = t(IntH_train@assays$integrated@scale.data),
      target       = IntH_train$y_diff,
      cost         = costM(X = t(IntH_train@assays$integrated@scale.data), y = IntH_train$y_diff),
      lambda_full  = 0.1, epsilon_full = 1e-7, maxiter = max_iter
    ),
    xtrain       = t(IntH_train@assays$integrated@scale.data),
    target       = IntH_train$y_diff,
    xtest        = t(IntH@assays$integrated@scale.data),
    ngenesselect = 25, lambda_red = 0.01, epsilon_red = 1e-7,
    maxiter      = max_iter, nfolds = 20
  )
  predD <- redD$pred
  if (!is.numeric(predD)) stop("Perm ", i, ": DIFF pred not numeric")
  if (is.null(names(predD))) names(predD) <- colnames(IntH)
  if (!all(colnames(IntH) %in% names(predD))) stop("Perm ", i, ": missing DIFF pred names")
  IntH$ordi_diff <- scales::rescale(predD[colnames(IntH)], to = c(0,1))
  if (anyNA(IntH$ordi_diff)) stop("Perm ", i, ": NA in ordi_diff after rescaling")
  saveRDS(redD$redmodel, file = file.path(models_dir, sprintf("redDIFF_perm_%02d.rds", i)))
  
  metricsDIFF <- IntH@meta.data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(min = min(ordi_diff),
                     max = max(ordi_diff),
                     median = median(ordi_diff),
                     sd = sd(ordi_diff)) %>%
    as.data.frame()
  readr::write_csv(metricsDIFF, file = file.path(models_dir, sprintf("metricsDIFF_perm_%02d.csv", i)))
  
  # 5.9) Landscape ordinals preparation & sanity
  diff_meta <- readr::read_csv("DiffTypes_Hm_Int.csv")
  IntH$diff_ek <- diff_meta$DiffTypes_final[match(colnames(IntH), diff_meta$cells)]
  
  # subsample training set for landscape ordinals
  grpL     <- table(IntH$group)
  sample_nL<- min(land_n, min(grpL))
  if (sample_nL < land_n) warning(sprintf("Perm %02d: land_n reduced to %d", i, sample_nL))
  cellsL   <- unlist(lapply(names(grpL), function(g) {
    cs <- colnames(IntH)[IntH$group == g]
    sample(cs, sample_nL)
  }))
  IntH_trainL <- subset(IntH, cells = cellsL)
  
  # For all permutations: use full IntH as landscape
  IntH_L <- IntH
  
  # For permutation 1 only: also save the CSV‐defined subset
  if (i == 1) {
    lands_sel <- readr::read_csv("SelectCells_Hm_Int_forLandscapes_400.csv")
    IntH_landscapeSubset <- subset(IntH, cells = lands_sel$x)
    saveRDS(IntH_landscapeSubset,
            file = file.path(land_dir, sprintf("landscapeSubset_perm_%02d.rds", i)))
  }
  
  # AGE on landscape
  redA_L <- custom_red_and_pred(
    fullmodel    = training_ordi_full(
      datamatrix   = t(IntH_trainL@assays$integrated@scale.data),
      target       = IntH_trainL$y_age,
      cost         = costM(X = t(IntH_trainL@assays$integrated@scale.data),
                           y = IntH_trainL$y_age),
      lambda_full  = 0.05, epsilon_full = 1e-7, maxiter = max_iter
    ),
    xtrain       = t(IntH_trainL@assays$integrated@scale.data),
    target       = IntH_trainL$y_age,
    xtest        = t(IntH_L@assays$integrated@scale.data),
    ngenesselect = 25, lambda_red = 0.01, epsilon_red = 1e-7,
    maxiter      = max_iter, nfolds = 20
  )
  predA_L <- redA_L$pred
  if (is.null(names(predA_L))) names(predA_L) <- colnames(IntH_L)
  IntH_L$ordi_age <- scales::rescale(predA_L[colnames(IntH_L)], to = c(0,1))
  if (anyNA(IntH_L$ordi_age)) stop("Perm ", i, ": NA in landscape ordi_age")
  
  # DIFF on landscape
  redD_L <- custom_red_and_pred(
    fullmodel    = training_ordi_full(
      datamatrix   = t(IntH_trainL@assays$integrated@scale.data),
      target       = IntH_trainL$y_diff,
      cost         = costM(X = t(IntH_trainL@assays$integrated@scale.data),
                           y = IntH_trainL$y_diff),
      lambda_full  = 0.1, epsilon_full = 1e-7, maxiter = max_iter
    ),
    xtrain       = t(IntH_trainL@assays$integrated@scale.data),
    target       = IntH_trainL$y_diff,
    xtest        = t(IntH_L@assays$integrated@scale.data),
    ngenesselect = 25, lambda_red = 0.01, epsilon_red = 1e-7,
    maxiter      = max_iter, nfolds = 20
  )
  predD_L <- redD_L$pred
  if (is.null(names(predD_L))) names(predD_L) <- colnames(IntH_L)
  IntH_L$ordi_diff <- scales::rescale(predD_L[colnames(IntH_L)], to = c(0,1))
  if (anyNA(IntH_L$ordi_diff)) stop("Perm ", i, ": NA in landscape ordi_diff")
  
  # Normalize landscape ordinals
  IntH_L$ordi_diff_norm <- ordi_normalize(
    ordiscore        = IntH_L$ordi_diff,
    ordigroups       = IntH_L$y_diff,
    ordigroups_other = IntH_L$y_age,
    limits1          = c(0, 0.24),
    limits2          = c(0.21, 0.47),
    limits3          = c(0.4, 1),
    applyWinsor      = TRUE
  )
  IntH_L$ordi_age_norm <- ordi_normalize(
    ordiscore        = IntH_L$ordi_age,
    ordigroups       = IntH_L$y_age,
    ordigroups_other = IntH_L$y_diff,
    applyWinsor      = TRUE
  )
  if (anyNA(IntH_L$ordi_age_norm) || anyNA(IntH_L$ordi_diff_norm))
    stop("Perm ", i, ": NA after normalizing landscape ordinals")
  
  saveRDS(IntH_L,
          file = file.path(land_dir, sprintf("IntH_L_perm_%02d.rds", i)))
  cat(sprintf("Permutation %02d landscape saved\n", i))
}

# 6) Aggregate & save fullAGE coefficients across Level1 permutations
coef_df <- data.frame(gene = names(all_fullAGE[[1]]), stringsAsFactors = FALSE)
for (i in seq_len(n_perms)) {
  coef_df[[sprintf("perm%02d", i)]] <- all_fullAGE[[i]]
}
coef_df$mean_perm <- rowMeans(coef_df[grep("^perm", names(coef_df))])
readr::write_csv(coef_df, file = file.path(models_dir, "coef_fullAGE_perms.csv"))

# 7) Compare to OG_humous reference & plot
ref_mod <- readRDS(file.path("OG_humous", "model_fullAGE_H.rds"))
if (is.list(ref_mod) && !is.null(ref_mod$coefficients)) {
  coef_ref <- ref_mod$coefficients
} else if (is.atomic(ref_mod)) {
  coef_ref <- ref_mod
} else {
  coef_ref <- suppressWarnings(coef(ref_mod))
}
ref_IntH  <- readRDS("OG_humous/IntH.rds")
genes_ref <- rownames(ref_IntH@assays$integrated@scale.data)
if (length(coef_ref) != length(genes_ref) || is.null(names(coef_ref))) {
  names(coef_ref) <- genes_ref
}
ref_df <- data.frame(gene = names(coef_ref), ref = as.numeric(coef_ref),
                     stringsAsFactors = FALSE)
plot_df <- dplyr::inner_join(ref_df, coef_df, by = "gene")
p <- ggplot(plot_df, aes(x = ref, y = mean_perm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6) + theme_bw() +
  labs(x = "Original fullAGE coeff", y = "Mean perm fullAGE coeff",
       title = "Original vs Mean Perm Coeffs") +
  annotate("text",
           x = min(plot_df$ref, na.rm = TRUE),
           y = max(plot_df$mean_perm, na.rm = TRUE),
           label = sprintf("rho=%.2f",
                           cor(plot_df$ref, plot_df$mean_perm,
                               use = "complete.obs")),
           hjust = 0, vjust = 1)
pdf(file.path(models_dir, "coef_fullAGE_comparison.pdf"), width = 6, height = 6)
print(p)
dev.off()

cat("\nAll permutations and landscape ordinals done.\nResults in:", 
    models_dir, "and", land_dir, "\n")

# 8) Metrics: overlap of final IntH_L_perm objects
perm_files <- list.files(land_dir,
                         pattern = "^IntH_L_perm_\\d+\\.rds$",
                         full.names = TRUE)
cell_lists <- lapply(perm_files, function(f) colnames(readRDS(f)))
names(cell_lists) <- basename(perm_files)
n <- length(cell_lists)
jaccard_mat <- matrix(NA, n, n,
                      dimnames = list(names(cell_lists),
                                      names(cell_lists)))
for (a in seq_len(n)) {
  for (b in seq_len(n)) {
    inter_ab <- length(intersect(cell_lists[[a]], cell_lists[[b]]))
    union_ab <- length(union(cell_lists[[a]], cell_lists[[b]]))
    jaccard_mat[a, b] <- inter_ab / union_ab
  }
}
# summarize
jac_vals <- jaccard_mat[lower.tri(jaccard_mat)]
metrics <- data.frame(
  min_jaccard  = min(jac_vals),
  mean_jaccard = mean(jac_vals),
  max_jaccard  = max(jac_vals),
  common_cells = length(Reduce(intersect, cell_lists))
)
print(metrics)
readr::write_csv(metrics, file = file.path(land_dir, "landscape_overlap_metrics.csv"))
