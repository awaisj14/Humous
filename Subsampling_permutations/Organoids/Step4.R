#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript Step4_landscapes.R <perm_id> [lands_dir out_obj_dir]")
}
perm_id <- as.integer(args[1])
if (is.na(perm_id)) stop("First arg must be a permutation ID (integer).")

lands_dir   <- if (length(args) >= 2) args[2] else "lands"
out_obj_dir <- if (length(args) >= 3) args[3] else "LandsH"

# ---- dependencies you already have ----
if (!file.exists("lib_misc.R"))  stop("lib_misc.R not found in CWD.")
if (!file.exists("lib_lands.R")) stop("lib_lands.R not found in CWD.")
source("lib_misc.R")   # for gene_selector
source("lib_lands.R")  # for knn_array_* and knn_rowMeans

# ---- helpers ----
ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pick_grouping <- function(meta) {
  # Prefer 'dataset' (as per your earlier scripts), otherwise try 'dataset_proto' then error
  if ("dataset" %in% colnames(meta)) return(meta$dataset)
  if ("dataset_proto" %in% colnames(meta)) return(meta$dataset_proto)
  stop("$dataset (or $dataset_proto) missing in metadata; required by gene_selector().")
}

pick_covariate <- function(seu, assay_to_use) {
  cands <- c(
    paste0("nFeature_", assay_to_use),
    paste0("nCount_",   assay_to_use),
    "nFeature_RNA", "nCount_RNA"
  )
  cands[cands %in% colnames(seu@meta.data)][1]
}

# ---- main ----
infile <- file.path(lands_dir, sprintf("IntH_L_perm_%02d.rds", perm_id))
if (!file.exists(infile)) stop("Landscape RDS not found: ", infile)

IntH_ordi <- readRDS(infile)
if (!inherits(IntH_ordi, "Seurat")) stop("Loaded object is not Seurat: ", infile)

# Assay handling (Seurat v5 safe): prefer RNA if present
assay_names  <- names(IntH_ordi@assays)
if (length(assay_names) == 0) stop("No assays present in: ", infile)
assay_to_use <- if ("RNA" %in% assay_names) "RNA" else assay_names[1]
DefaultAssay(IntH_ordi) <- assay_to_use
# Drop integrated assay if present (we'll work in a single assay space)
if ("integrated" %in% assay_names) IntH_ordi[["integrated"]] <- NULL

# Gene selection by dataset grouping
grouping <- pick_grouping(IntH_ordi@meta.data)
genes_keep <- gene_selector(
  seuratobj          = IntH_ordi,
  grouping           = grouping,
  threshold_grouping = 0,
  threshold_ncells   = 9
)
if (length(genes_keep) == 0) stop("gene_selector() returned 0 genes.")

LandsS_H <- subset(IntH_ordi, features = genes_keep)

# Normalize + scale (regress the best available nFeature_/nCount_ covariate)
covariate <- pick_covariate(LandsS_H, assay_to_use)
LandsS_H  <- NormalizeData(LandsS_H, assay = assay_to_use, verbose = FALSE)

if (length(covariate)) {
  LandsS_H <- ScaleData(
    LandsS_H,
    assay = assay_to_use,
    vars.to.regress = covariate,
    features = rownames(LandsS_H),
    verbose = FALSE
  )
} else {
  message("[WARN] No matching nFeature_/nCount_ covariate; scaling without regression.")
  LandsS_H <- ScaleData(
    LandsS_H,
    assay = assay_to_use,
    features = rownames(LandsS_H),
    verbose = FALSE
  )
}

# Shift scale.data to be non-negative (required for downstream grids)
sdat <- GetAssayData(LandsS_H, slot = "scale.data", assay = assay_to_use)
smin <- suppressWarnings(min(sdat))
if (is.finite(smin) && smin < 0) {
  LandsS_H@assays[[assay_to_use]]@scale.data <- sdat + abs(smin)
}

# Ordinal sanity check (should exist from Step 3)
if (is.null(LandsS_H$ordi_age_norm) || is.null(LandsS_H$ordi_diff_norm)) {
  stop("ordi_age_norm/ordi_diff_norm missing in metadata. Was the ordinal Step 3 completed?")
}

# ---- outputs ----
ensure_dir(out_obj_dir)
outfile_obj <- file.path(out_obj_dir, sprintf("LandsS_H_perm_%02d.rds", perm_id))
saveRDS(LandsS_H, outfile_obj)
message("Saved LandsS_H → ", outfile_obj)

# Grids (kNN smoothing on the normalized ordinal axes)
gridlistH <- list(
  grid_medR  = knn_array_medres (LandsS_H$ordi_age_norm, LandsS_H$ordi_diff_norm, k = 100L),
  grid_lowR  = knn_array_lowres (LandsS_H$ordi_age_norm, LandsS_H$ordi_diff_norm, k = 100L),
  grid_highR = knn_array_highres(LandsS_H$ordi_age_norm, LandsS_H$ordi_diff_norm, k = 100L)
)
saveRDS(gridlistH, file.path(lands_dir, sprintf("grids_perm_%02d.rds", perm_id)))

# Medium-res matrix (row-means of scaled expression over medium-res grid bins)
m  <- GetAssayData(LandsS_H, slot = "scale.data", assay = assay_to_use)
cm <- knn_rowMeans(m, gridlistH$grid_medR)
rownames(cm) <- rownames(m)
outf <- file.path(lands_dir, sprintf("L_MR_H_perm_%02d.rds", perm_id))
saveRDS(cm, outf)

cat("Permutation", sprintf("%02d", perm_id), "→ saved", outf, "\n")
