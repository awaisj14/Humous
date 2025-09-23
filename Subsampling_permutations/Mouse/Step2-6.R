#!/usr/bin/env Rscript

# run_ordinals_from_draws.R
# SLURM-compatible pipeline continuing from precomputed draw_cells_perm_XX.rds
# Uses pre-sampled cells and runs integration, modeling, normalization,
# landscape generation, and optional reference comparison.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(parallel)
  library(matrixStats)
  library(entropy)
  library(DescTools)
  library(abind)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i == length(args)) return(TRUE)
  args[i + 1]
}

# ----------------------------
# CLI Arguments
# ----------------------------
lib_misc   <- get_arg("--lib", "lib_misc.R")
lib_land   <- get_arg("--landlib", "lib_lands.R")
base_out   <- get_arg("--out_root", "permutations_all")
n_perms    <- as.integer(get_arg("--n_perms", "5"))
cores      <- as.integer(get_arg("--cores", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
ref_lmr    <- get_arg("--ref_lmr")
merged_rds <- get_arg("--merged", "merged_all.rds")
perm_i     <- as.integer(get_arg("--perm", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))

source(lib_misc)
source(lib_land)

options(mc.cores = cores)

message("🚀 Running permutation ", perm_i, " / ", n_perms)
i <- perm_i
i_padded <- sprintf("%02d", i)

# ----------------------------
# Create output directory
# ----------------------------
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Reference landscape loading
# ----------------------------
L_ref <- if (!is.na(ref_lmr) && file.exists(ref_lmr)) {
  x <- readRDS(ref_lmr)
  if (length(dim(x)) == 3) {
    if (inherits(x, "DelayedArray")) as.array(x) else x
  } else stop("Invalid reference landscape.")
} else NULL

# ----------------------------
# Utility functions
# ----------------------------
flatten_sd <- function(L3) matrixStats::rowSds(apply(L3, 1, as.vector) |> t(), na.rm = TRUE)
calc_entropy_from_L3 <- function(L3) apply(L3, 1, function(m) entropy::entropy(as.vector(m), method = "CS"))
per_gene_corr <- function(L3a, L3b) {
  g <- intersect(rownames(L3a), rownames(L3b))
  tibble(gene = g, corr = sapply(g, function(x) {
    a <- as.vector(L3a[x,,]); b <- as.vector(L3b[x,,])
    if (sd(a) > 0 && sd(b) > 0) suppressWarnings(cor(a, b)) else NA_real_
  }))
}

# ----------------------------
# Checkpoints
# ----------------------------
checkpoints <- list(
  merged_f = file.path(base_out, paste0("merged_f_perm", i_padded, ".rds")),
  AGE      = file.path(base_out, paste0("fullAGE_perm", i_padded, ".rds")),
  DIFF     = file.path(base_out, paste0("fullDIFF_perm", i_padded, ".rds")),
  grid     = file.path(base_out, paste0("gridlist_perm", i_padded, ".rds")),
  IntM     = file.path(base_out, paste0("IntM_perm", i_padded, ".rds")),
  LMR      = file.path(base_out, paste0("L_MR_M_perm", i_padded, ".rds")),
  cells    = file.path(base_out, paste0("draw_cells_perm_", i_padded, ".rds")),
  compare  = if (!is.null(L_ref)) file.path(base_out, paste0("perm_vs_ref_per_gene_perm", i_padded, ".csv")) else NULL
)

# ----------------------------
# Skip if this permutation already completed
# ----------------------------
if (file.exists(checkpoints$LMR)) {
  message("✔ Permutation ", i, " already completed. Skipping.")
  quit(status = 0)
}

# ----------------------------
# Load merged object + selected cells
# ----------------------------
merged_all <- readRDS(merged_rds)
if (!file.exists(checkpoints$cells)) {
  stop("❌ Missing precomputed draw file: ", checkpoints$cells)
}
selected <- readRDS(checkpoints$cells)

# Subset merged object to selected cells
merged_f <- subset(merged_all, cells = selected)
saveRDS(merged_f, checkpoints$merged_f)

# ----------------------------
# Integration
# ----------------------------
message("🔗 Running integration...")
full.list.f <- SplitObject(merged_f, split.by = "dataset")
IntM <- dataset_integration(full.list.f, nfeats = 10000,
                            FindIntDims = 1:15,
                            FindInt_kScore = 20,
                            IntKweight = 20)
DefaultAssay(IntM) <- "integrated"
IntM <- ScaleData(IntM, vars.to.regress = "nFeature_RNA") |>
  RunPCA() |>
  RunUMAP(dims = 1:25, return.model = TRUE)

# ----------------------------
# Train AGE + DIFF models
# ----------------------------
IntM$y_age <- as.integer(plyr::revalue(IntM$age_ek,
                                       c("12"="1","13"="2","14"="3","15"="4","16-17"="5")))
IntM$y_diff <- as.integer(plyr::revalue(IntM$diff_ek,
                                        c("RG"="1","IPC"="2","N"="3")))
IntM$grouping_ordi <- paste0(IntM$y_age, "_", IntM$y_diff)
IntM_train <- cell.selector(IntM, colnames(IntM), IntM$grouping_ordi, n = 50)

message("🧠 Training AGE model...")
Xtrain <- t(as.matrix(IntM_train@assays$integrated@scale.data))
ytrain <- IntM_train$y_age
AGEmod <- training_ordi_full(Xtrain, ytrain, costM(Xtrain, ytrain),
                             lambda_full = 0.05,
                             epsilon_full = 1e-7,
                             maxiter = 1000)
saveRDS(AGEmod, checkpoints$AGE)

message("🧠 Training DIFF model...")
ytrain_diff <- IntM_train$y_diff
DIFFmod <- training_ordi_full(Xtrain, ytrain_diff, costM(Xtrain, ytrain_diff),
                              lambda_full = 0.05,
                              epsilon_full = 1e-7,
                              maxiter = 1000)
saveRDS(DIFFmod, checkpoints$DIFF)

# ----------------------------
# Normalize + predict
# ----------------------------
message("🌄 Generating grid and landscape...")
DefaultAssay(IntM) <- "RNA"
IntM <- NormalizeData(IntM) |> ScaleData(vars.to.regress = "nFeature_RNA")
sc_data <- IntM@assays$RNA@scale.data
if (min(sc_data) < 0) IntM@assays$RNA@scale.data <- sc_data + abs(min(sc_data))

IntM$ordi_age <- as.numeric(t(as.matrix(IntM@assays$integrated@scale.data)) %*% AGEmod)
IntM$ordi_diff <- as.numeric(t(as.matrix(IntM@assays$integrated@scale.data)) %*% DIFFmod)
IntM$ordi_age_norm <- ordi_normalize(IntM$ordi_age, IntM$y_age, IntM$y_diff, applyWinsor = TRUE)
IntM$ordi_diff_norm <- ordi_normalize(IntM$ordi_diff, IntM$y_diff, IntM$y_age,
                                      limits1 = c(0,0.26),
                                      limits2 = c(0.21,0.5),
                                      limits3 = c(0.4,1),
                                      applyWinsor = TRUE)
saveRDS(IntM, checkpoints$IntM)

# ----------------------------
# Grid + landscape
# ----------------------------
gridlist <- list(grid_medR = knn_array_medres(IntM$ordi_age_norm,
                                              IntM$ordi_diff_norm,
                                              k = 200L))
saveRDS(gridlist, checkpoints$grid)

L_MR_M <- local({
  all_genes <- rownames(IntM)
  m <- IntM@assays$RNA@scale.data[all_genes, ]
  landscapes_MR_M <- knn_rowMeans(m, gridlist$grid_medR)
  rownames(landscapes_MR_M) <- all_genes
  landscapes_MR_M
})
saveRDS(L_MR_M, checkpoints$LMR)

# ----------------------------
# Compare to reference landscape
# ----------------------------
if (!is.null(L_ref)) {
  message("🔍 Comparing to reference...")
  sd_ref <- flatten_sd(L_ref)
  sd_perm <- flatten_sd(L_MR_M)
  ent_ref <- calc_entropy_from_L3(L_ref)
  ent_perm <- calc_entropy_from_L3(L_MR_M)
  corr_df <- per_gene_corr(L_MR_M, L_ref)
  df <- corr_df |>
    mutate(sd_ref = sd_ref[gene],
           sd_perm = sd_perm[gene],
           ent_ref = ent_ref[gene],
           ent_perm = ent_perm[gene])
  write_csv(df, checkpoints$compare)
}

message("✅ Permutation ", i, " complete.")

# ----------------------------
# Pairwise Jaccard overlap — only run on last job
# ----------------------------
if (i == n_perms) {
  message("🔎 Calculating between-permutation overlaps...")
  sel_list <- lapply(seq_len(n_perms), function(k) {
    readRDS(file.path(base_out, sprintf("draw_cells_perm_%02d.rds", k)))
  })
  names(sel_list) <- sprintf("perm_%02d", seq_len(n_perms))
  
  pm <- length(sel_list)
  jaccard_matrix <- matrix(NA_real_, pm, pm,
                           dimnames = list(names(sel_list), names(sel_list)))
  for (a in seq_len(pm)) for (b in seq_len(pm)) {
    inter <- length(intersect(sel_list[[a]], sel_list[[b]]))
    uni   <- length(union(sel_list[[a]], sel_list[[b]]))
    jaccard_matrix[a,b] <- if (uni == 0) NA_real else inter/uni
  }
  
  write.csv(jaccard_matrix,
            file = file.path(base_out, "jaccard_between_perms.csv"),
            row.names = TRUE)
  message("✅ Pairwise Jaccard overlap saved.")
}
