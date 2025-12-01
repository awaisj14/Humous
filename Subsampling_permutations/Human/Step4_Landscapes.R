#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
perm_id <- as.integer(args[1])
if (is.na(perm_id)) stop("Please provide permutation ID as first argument")

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(png)

source("lib_misc.R")
source("lib_lands.R")

# 1) load the permuted object
infile <- sprintf("lands/IntH_L_perm_%02d.rds", perm_id)
IntH_ordi <- readRDS(infile)

# 2) switch to RNA assay & drop integrated
DefaultAssay(IntH_ordi) <- "RNA"
IntH_ordi[["integrated"]] <- NULL

# 3) gene selection
genestokeep_h <- gene_selector(
  seuratobj       = IntH_ordi,
  grouping        = IntH_ordi$dataset,
  threshold_grouping = 0,
  threshold_ncells   = 9
)
LandsS_H <- subset(IntH_ordi, features = genestokeep_h)

# 4) normalize & scale (make all positives)
LandsS_H <- LandsS_H %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = "nFeature_RNA",
            features = rownames(LandsS_H))
# shift so min = 0
LandsS_H@assays$RNA@scale.data <-
  LandsS_H@assays$RNA@scale.data +
  abs(min(LandsS_H@assays$RNA@scale.data))

# 5) build grids
gridlistH <- list(
  grid_medR  = knn_array_medres (LandsS_H$ordi_age_norm, LandsS_H$ordi_diff_norm, k = 100L),
  grid_lowR  = knn_array_lowres(LandsS_H$ordi_age_norm, LandsS_H$ordi_diff_norm, k = 100L),
  grid_highR = knn_array_highres(LandsS_H$ordi_age_norm, LandsS_H$ordi_diff_norm, k = 100L)
)

# 6) medium‐resolution, all genes
L_MR_H <- local({
  m <- LandsS_H@assays$RNA@scale.data
  cm <- knn_rowMeans(m, gridlistH$grid_medR)
  rownames(cm) <- rownames(m)
  cm
})

# 7) save result
outf <- sprintf("lands/L_MR_H_perm_%02d.rds", perm_id)
saveRDS(L_MR_H, outf)
cat("Permutation", perm_id, "→ saved", outf, "\n")
