#!/usr/bin/env Rscript

# ── SLURM ARRAY SETUP ─────────────────────────────────────────────────────────
# Usage: sbatch --array=1-100 --cpus-per-task=1 integration_script.R
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  perm_id <- args[1]
} else {
  idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  stopifnot(!is.na(idx))
  perm_files <- list.files("subsamples", pattern = "^perm_.*\\.rds$", full.names = TRUE)
  stopifnot(idx <= length(perm_files))
  perm_id <- tools::file_path_sans_ext(basename(perm_files[idx]))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(digest)
  library(Matrix)
})
options(Seurat.object.assay.version = "v3")

# ── PATHS ─────────────────────────────────────────────────────────────────────
IntH_path <- "/home/users/j/javed/humous_reviews/OG_humous/IntH.rds"
out_dir   <- "IntO"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
perm_path <- file.path("subsamples", paste0(perm_id, ".rds"))
message("=== Processing ", perm_id, " ===")

# ── LOAD HELPERS ──────────────────────────────────────────────────────────────
source("/home/users/j/javed/humous_reviews/lib_misc.R")
stopifnot(exists("double_dataset_integration", mode = "function"))

# ── STEP 1: Load and build IntH2 ──────────────────────────────────────────────
IntH <- readRDS(IntH_path)
IntH[["integrated"]] <- as(IntH[["integrated"]], "Assay")
mat <- as(GetAssayData(IntH, assay = "integrated", slot = "data"), "dgCMatrix")
IntH2 <- CreateSeuratObject(counts = mat, meta.data = IntH@meta.data)
IntH2@assays$RNA@data <- IntH2@assays$RNA@counts
IntH2[["RNA"]] <- as(IntH2[["RNA"]], "Assay")
IntH2$dataset_proto <- "IntH2"
DefaultAssay(IntH2) <- "RNA"

# ── STEP 2: Load permutation-specific organoid ─────────────────────────────────
mergedO_f <- readRDS(perm_path)
mergedO_f@assays <- map(mergedO_f@assays, ~ as(.x, "Assay"))

# ── STEP 3: Merge, split, normalize & find features ───────────────────────────
merged_list <- merge(IntH2, mergedO_f)
full.list.int <- SplitObject(merged_list, split.by = "dataset_proto")
full.list.int <- imap(full.list.int, function(x, nm) {
  if (nm != "IntH2") x <- NormalizeData(x, verbose = FALSE)
  FindVariableFeatures(x, selection.method = "mvp", nfeatures = 10000, verbose = FALSE)
})

# Extract and cap HVGs
hvg_lists <- map(full.list.int, VariableFeatures)
common <- reduce(hvg_lists, intersect)
if (length(common) == 0) {
  combined <- reduce(hvg_lists, union)
  features <- if (length(combined) > 10000) sample(combined, 10000) else combined
  message("⚠️  Using union-derived features: ", length(features))
} else if (length(common) > 10000) {
  features <- sample(common, 10000)
  message("🔀 Sampled 10000 intersecting features")
} else {
  features <- common
}

# ── STEP 4: Double dataset integration ────────────────────────────────────────
# double_dataset_integration expects `full.list` in scope
full.list <- full.list.int
IntHO <- double_dataset_integration(
  list_seuratobjs    = full.list.int,
  nfeats             = 10000,
  FindIntDims        = 1:15,
  FindInt_kScore     = 20,
  IntKweight         = 20,
  reference_dataset  = "IntH2",
  features           = features
)

# ── SAVE POST-INTEGRATION OBJECT ─────────────────────────────────────────────
post_int_path <- file.path(out_dir, paste0("IntHO_", perm_id, "_postIntegration.rds"))
saveRDS(IntHO, post_int_path)
message("💾 Saved post-integration: ", post_int_path)

# ── STEP 5: Scaling, PCA, UMAP & save final ──────────────────────────────────
DefaultAssay(IntHO) <- "integrated"
IntHO <- IntHO %>%
  ScaleData(vars.to.regress   = "nFeature_RNA") %>%
  RunPCA() %>%
  RunUMAP(dims               = 1:25,
          return.model       = TRUE,
          min.dist           = 0.5,
          n.neighbors        = 50,
          local.connectivity = 6)

final_path <- file.path(out_dir, paste0("IntHO_", perm_id, ".rds"))
saveRDS(IntHO, final_path)
message("✅ Saved final IntHO: ", final_path)

# ── DIAGNOSTIC PLOTS ─────────────────────────────────────────────────────────
plots_pdf <- file.path(out_dir, paste0("IntHO_", perm_id, "_plots.pdf"))
pdf(plots_pdf, width = 10, height = 8)
DimPlot(IntHO, group.by = "dataset_proto", pt.size = 0.5, split.by = "dataset")
DimPlot(IntHO, group.by = "age_ek",        pt.size = 0.5, split.by = "dataset_proto")
DimPlot(IntHO, group.by = "diff_ek",       pt.size = 0.5, split.by = "dataset_proto")
dev.off()
message("🔍 Saved diagnostic plots to ", plots_pdf)
