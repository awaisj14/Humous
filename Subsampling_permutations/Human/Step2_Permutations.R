#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
source("lib_misc.R")    # dataset_integration

# read SLURM_ARRAY_TASK_ID (1–10)
perm_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(perm_id)) stop("SLURM_ARRAY_TASK_ID not set")

# load the appropriate prefiltered list
full.list.f <- readRDS(sprintf("full.list_f_%02d.rds", perm_id))

cat(sprintf("[Perm %02d] Running integration…\n", perm_id))
IntH <- dataset_integration(
  list_seuratobjs = full.list.f,
  nfeats          = 10000,
  FindIntDims     = 1:15,
  FindInt_kScore  = 20,
  IntKweight      = 20
)
DefaultAssay(IntH) <- "integrated"
IntH <- IntH %>%
  ScaleData(vars.to.regress="nFeature_RNA") %>%
  RunPCA() %>%
  RunUMAP(
    dims               = 1:15,
    return.model       = TRUE,
    min.dist           = 0.8,
    n.neighbors        = 100,
    local.connectivity = 6
  )
IntH$diff_ek <- ifelse(IntH$diff_ek=="Neuron","N", as.character(IntH$diff_ek))
cat(sprintf("[Perm %02d] Done. Saving IntH_perm_%02d.rds\n", perm_id, perm_id))
saveRDS(IntH, file = sprintf("IntH_perm_%02d.rds", perm_id))

