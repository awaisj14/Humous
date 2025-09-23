#!/usr/bin/env Rscript

# run_ordinals_serial_full.R
# Nested ordinal regression pipeline (AGE & DIFF) + landscape ordinals
# IntO inputs, removal of IntH2 cluster, new age/diff maps, 4 perms, lambda=10
# Saves: per-permutation mean full-model coefficients (AGE & DIFF),
#        reduced models (AGE & DIFF), landscapes per perm,
#        aggregated coefficient tables, and comparison plots.

# 1) Load utilities & libraries
source("/home/users/j/javed/humous_reviews/lib_misc.R")   # training_ordi_full, custom_red_and_pred, costM, ordi_normalize
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

# 2) Prepare output directories
models_dir <- "models"
land_dir   <- "lands"
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(land_dir,   recursive = TRUE, showWarnings = FALSE)

# 3) Parameters
n_perms   <- 4
n_lvl2    <- 5
train_n   <- 50        # cells per (age,diff) group for reduced models
max_iter  <- 1000
lambda_f  <- 10
lambda_r  <- 10
ngenes_r  <- 25
nfolds_r  <- 20
land_n    <- 100       # not used directly here but kept for compatibility
base_seed <- 777

age_map  <- c("1-2 m"=1, "2-3 m"=2, "3-5 m"=3, "6-7 m"=4)
diff_map <- c("RG"=1, "IPC"=2, "N"=3)

# 4) Storage for full-model coefficients across permutations
all_fullAGE  <- vector("list", n_perms)
all_fullDIFF <- vector("list", n_perms)

# 5) Main loop over permutations
for (i in seq_len(n_perms)) {
  cat(sprintf("\n=== Permutation %02d ===\n", i))
  set.seed(base_seed + i)
  
  int_file <- file.path("IntO", sprintf("IntHO_perm_%d.rds", i))
  if (!file.exists(int_file)) stop("Missing file: ", int_file)
  IntH <- readRDS(int_file)
  IntH <- subset(IntH, subset = dataset_proto != "IntH2")
  
  if (!"integrated" %in% Assays(IntH)) stop("Perm ", i, ": missing 'integrated' assay")
  sd_mat <- IntH@assays$integrated@scale.data
  if (is.null(sd_mat) || nrow(sd_mat)==0 || ncol(sd_mat)==0) stop("Perm ", i, ": empty integrated@scale.data")
  
  IntH$y_age  <- as.integer(age_map[ as.character(IntH$age_ek) ])
  IntH$y_diff <- as.integer(diff_map[ as.character(IntH$diff_ek) ])
  if (anyNA(IntH$y_age) || anyNA(IntH$y_diff)) stop("Perm ", i, ": NA in y_age/y_diff labels")
  
  # Balanced sampling across (age,diff) groups for reduced models
  IntH$group   <- paste0(IntH$y_age, "_", IntH$y_diff)
  grp_counts   <- table(IntH$group)
  sample_n1    <- min(train_n, min(grp_counts))
  if (sample_n1 < train_n) warning(sprintf("Perm %02d: train_n reduced to %d (min group size)", i, sample_n1))
  
  train_cells1 <- unlist(lapply(names(grp_counts), function(g) {
    cs <- colnames(IntH)[IntH$group == g]
    set.seed(base_seed + i + as.integer(gsub("_","",g)))
    sample(cs, sample_n1)
  }))
  IntH_train <- subset(IntH, cells = train_cells1)
  if (length(colnames(IntH_train)) != length(train_cells1)) stop("Perm ", i, ": mismatch in IntH_train cell count")
  
  # --------- Level-2 bootstrap / stability averaging: AGE ---------
  lvl2_age_coefs <- vector("list", n_lvl2)
  for (j in seq_len(n_lvl2)) {
    set.seed(base_seed + i*100 + j)
    cells_j <- unlist(lapply(names(grp_counts), function(g) {
      cs <- colnames(IntH)[IntH$group == g]
      sample(cs, sample_n1)
    }))
    tmpA <- subset(IntH, cells = cells_j)
    mA <- training_ordi_full(
      datamatrix   = t(tmpA@assays$integrated@scale.data),
      target       = tmpA$y_age,
      cost         = costM(X = t(tmpA@assays$integrated@scale.data), y = tmpA$y_age),
      lambda_full  = lambda_f, epsilon_full = 1e-7, maxiter = max_iter
    )
    if (is.list(mA) && !is.null(mA$coefficients)) cA <- mA$coefficients
    else if (is.atomic(mA))                        cA <- mA
    else                                           cA <- suppressWarnings(coef(mA))
    if (is.null(names(cA))) names(cA) <- rownames(sd_mat)
    lvl2_age_coefs[[j]] <- cA
  }
  mean_coefAGE        <- Reduce(`+`, lvl2_age_coefs) / n_lvl2
  all_fullAGE[[i]]    <- mean_coefAGE
  write_csv(data.frame(gene = names(mean_coefAGE), coeff = mean_coefAGE),
            file.path(models_dir, sprintf("mean_fullAGE_coeffs_perm_%02d.csv", i)))
  
  # Reduced AGE model + predictions
  redA <- custom_red_and_pred(
    fullmodel    = lvl2_age_coefs[[1]],
    xtrain       = t(IntH_train@assays$integrated@scale.data),
    target       = IntH_train$y_age,
    xtest        = t(IntH@assays$integrated@scale.data),
    ngenesselect = ngenes_r, lambda_red = lambda_r, epsilon_red = 1e-7,
    maxiter      = max_iter, nfolds = nfolds_r
  )
  predA <- redA$pred
  if (!is.numeric(predA)) stop("Perm ", i, ": AGE pred not numeric")
  if (is.null(names(predA))) names(predA) <- colnames(IntH)
  IntH$ordi_age <- scales::rescale(predA[colnames(IntH)], to = c(0,1))
  if (anyNA(IntH$ordi_age)) stop("Perm ", i, ": NA in ordi_age")
  saveRDS(redA$redmodel, file = file.path(models_dir, sprintf("redAGE_perm_%02d.rds", i)))
  
  # --------- Level-2 bootstrap / stability averaging: DIFF ---------
  lvl2_diff_coefs <- vector("list", n_lvl2)
  for (j in seq_len(n_lvl2)) {
    set.seed(base_seed + i*200 + j)
    cells_j <- unlist(lapply(names(grp_counts), function(g) {
      cs <- colnames(IntH)[IntH$group == g]
      sample(cs, sample_n1)
    }))
    tmpD <- subset(IntH, cells = cells_j)
    mD <- training_ordi_full(
      datamatrix   = t(tmpD@assays$integrated@scale.data),
      target       = tmpD$y_diff,
      cost         = costM(X = t(tmpD@assays$integrated@scale.data), y = tmpD$y_diff),
      lambda_full  = lambda_f, epsilon_full = 1e-7, maxiter = max_iter
    )
    if (is.list(mD) && !is.null(mD$coefficients)) cD <- mD$coefficients
    else if (is.atomic(mD))                        cD <- mD
    else                                           cD <- suppressWarnings(coef(mD))
    if (is.null(names(cD))) names(cD) <- rownames(sd_mat)
    lvl2_diff_coefs[[j]] <- cD
  }
  mean_coefDIFF        <- Reduce(`+`, lvl2_diff_coefs) / n_lvl2
  all_fullDIFF[[i]]    <- mean_coefDIFF
  write_csv(data.frame(gene = names(mean_coefDIFF), coeff = mean_coefDIFF),
            file.path(models_dir, sprintf("mean_fullDIFF_coeffs_perm_%02d.csv", i)))
  
  # Reduced DIFF model + predictions
  redD <- custom_red_and_pred(
    fullmodel    = lvl2_diff_coefs[[1]],
    xtrain       = t(IntH_train@assays$integrated@scale.data),
    target       = IntH_train$y_diff,
    xtest        = t(IntH@assays$integrated@scale.data),
    ngenesselect = ngenes_r, lambda_red = lambda_r, epsilon_red = 1e-7,
    maxiter      = max_iter, nfolds = nfolds_r
  )
  predD <- redD$pred
  if (!is.numeric(predD)) stop("Perm ", i, ": DIFF pred not numeric")
  if (is.null(names(predD))) names(predD) <- colnames(IntH)
  IntH$ordi_diff <- scales::rescale(predD[colnames(IntH)], to = c(0,1))
  if (anyNA(IntH$ordi_diff)) stop("Perm ", i, ": NA in ordi_diff")
  saveRDS(redD$redmodel, file = file.path(models_dir, sprintf("redDIFF_perm_%02d.rds", i)))
  
  # Normalizations + save landscape
  IntH$ordi_diff_norm <- ordi_normalize(
    IntH$ordi_diff, IntH$y_diff, IntH$y_age,
    limits1 = c(0,0.24), limits2 = c(0.21,0.47), limits3 = c(0.4,1), applyWinsor = TRUE
  )
  IntH$ordi_age_norm  <- ordi_normalize(
    IntH$ordi_age,  IntH$y_age,  IntH$y_diff, applyWinsor = TRUE
  )
  saveRDS(IntH, file = file.path(land_dir, sprintf("IntH_L_perm_%02d.rds", i)))
  cat(sprintf("Permutation %02d landscape saved\n", i))
}

# 6) Aggregate fullAGE coefficients across permutations
coef_df_AGE <- data.frame(gene = names(all_fullAGE[[1]]), stringsAsFactors = FALSE)
for (i in seq_len(n_perms)) {
  coef_df_AGE[[sprintf("perm%02d", i)]] <- all_fullAGE[[i]]
}
coef_df_AGE$mean_perm <- rowMeans(coef_df_AGE[grep("^perm", names(coef_df_AGE))])
write_csv(coef_df_AGE, file = file.path(models_dir, "coef_fullAGE_perms.csv"))

# 6b) Aggregate fullDIFF coefficients across permutations
coef_df_DIFF <- data.frame(gene = names(all_fullDIFF[[1]]), stringsAsFactors = FALSE)
for (i in seq_len(n_perms)) {
  coef_df_DIFF[[sprintf("perm%02d", i)]] <- all_fullDIFF[[i]]
}
coef_df_DIFF$mean_perm <- rowMeans(coef_df_DIFF[grep("^perm", names(coef_df_DIFF))])
write_csv(coef_df_DIFF, file = file.path(models_dir, "coef_fullDIFF_perms.csv"))

# 7) Compare with reference (AGE)
ref_IntH <- readRDS("/home/users/j/javed/humous_reviews/Org/IntO_OG/IntO.rds")
ref_modA <- readRDS("/home/users/j/javed/humous_reviews/Org/ordiO_OG/model_fullAGE_O.rds")
coef_refA <- if (is.list(ref_modA) && !is.null(ref_modA$coefficients)) ref_modA$coefficients
else if (is.atomic(ref_modA)) ref_modA else suppressWarnings(coef(ref_modA))
if (is.null(names(coef_refA))) names(coef_refA) <- rownames(ref_IntH@assays$integrated@scale.data)

ref_dfA  <- data.frame(gene = names(coef_refA), ref = as.numeric(coef_refA), stringsAsFactors = FALSE)
plot_dfA <- dplyr::inner_join(ref_dfA, coef_df_AGE, by = "gene")

pA <- ggplot(plot_dfA, aes(x = ref, y = mean_perm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(
    x = "Original fullAGE coeff",
    y = "Mean perm fullAGE coeff",
    title = "Original vs Mean Perm Coeffs (fullAGE)"
  ) +
  annotate(
    "text",
    x = min(plot_dfA$ref, na.rm = TRUE),
    y = max(plot_dfA$mean_perm, na.rm = TRUE),
    label = sprintf("rho=%.2f", cor(plot_dfA$ref, plot_dfA$mean_perm, use = "complete.obs")),
    hjust = 0, vjust = 1
  )
pdf(file.path(models_dir, "coef_fullAGE_comparison.pdf"), width = 6, height = 6); print(pA); dev.off()

# 7b) Compare with reference (DIFF)
ref_modD <- readRDS("/home/users/j/javed/humous_reviews/Org/ordiO_OG/model_fullDIFF_O.rds")
coef_refD <- if (is.list(ref_modD) && !is.null(ref_modD$coefficients)) ref_modD$coefficients
else if (is.atomic(ref_modD)) ref_modD else suppressWarnings(coef(ref_modD))
if (is.null(names(coef_refD))) names(coef_refD) <- rownames(ref_IntH@assays$integrated@scale.data)

ref_dfD  <- data.frame(gene = names(coef_refD), ref = as.numeric(coef_refD), stringsAsFactors = FALSE)
plot_dfD <- dplyr::inner_join(ref_dfD, coef_df_DIFF, by = "gene")

pD <- ggplot(plot_dfD, aes(x = ref, y = mean_perm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(
    x = "Original fullDIFF coeff",
    y = "Mean perm fullDIFF coeff",
    title = "Original vs Mean Perm Coeffs (fullDIFF)"
  ) +
  annotate(
    "text",
    x = min(plot_dfD$ref, na.rm = TRUE),
    y = max(plot_dfD$mean_perm, na.rm = TRUE),
    label = sprintf("rho=%.2f", cor(plot_dfD$ref, plot_dfD$mean_perm, use = "complete.obs")),
    hjust = 0, vjust = 1
  )
pdf(file.path(models_dir, "coef_fullDIFF_comparison.pdf"), width = 6, height = 6); print(pD); dev.off()

# 8) Overlap metrics across saved landscapes
perm_files <- list.files(land_dir, pattern = "^IntH_L_perm_\\d+\\.rds$", full.names = TRUE)
cell_lists <- lapply(perm_files, function(f) colnames(readRDS(f)))
n <- length(cell_lists)
jaccard_mat <- matrix(NA, n, n, dimnames = list(basename(perm_files), basename(perm_files)))
for (a in seq_len(n)) {
  for (b in seq_len(n)) {
    inter_ab <- length(intersect(cell_lists[[a]], cell_lists[[b]]))
    union_ab <- length(union(cell_lists[[a]], cell_lists[[b]]))
    jaccard_mat[a, b] <- inter_ab / union_ab
  }
}
jac_vals <- jaccard_mat[lower.tri(jaccard_mat)]
metrics <- data.frame(
  min_jaccard  = min(jac_vals),
  mean_jaccard = mean(jac_vals),
  max_jaccard  = max(jac_vals),
  common_cells = length(Reduce(intersect, cell_lists))
)
print(metrics)
write_csv(metrics, file = file.path(land_dir, "landscape_overlap_metrics.csv"))

cat("\nAll permutations and landscape ordinals done.\nResults in:", models_dir, "and", land_dir, "\n")
