#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
perm <- as.integer(args[1])
if (is.na(perm)) stop("Usage: Rscript make_all_perm_pngs.R <perm_id>")

library(parallel)
library(ggplot2)
library(png)
library(raster)

# adjust these if you put your scripts elsewhere
source("lib_misc.R")
source("lib_lands.R")

# 1) load the precomputed medium-res landscape for this perm
L_MR_H <- readRDS(sprintf("lands/L_MR_H_perm_%02d.rds", perm))

# 2) make sure output dirs exist (one subfolder per perm)
web_dir <- sprintf("pngLs/web_pngLs/O/perm_%02d", perm)
ml_dir  <- sprintf("pngLs/ml_pngLs/O/perm_%02d",  perm)
dir.create(web_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(ml_dir,  recursive=TRUE, showWarnings=FALSE)

genes <- rownames(L_MR_H)

# Helper to open/close PNG safely and catch errors
safe_png <- function(file, width, height, expr, label) {
  png(file, width=width, height=height)
  on.exit(dev.off(), add=TRUE)
  tryCatch(
    expr,
    error = function(e) {
      message(sprintf("[perm %02d][%s] Error plotting %s: %s", 
                      perm, label, basename(file), e$message))
    }
  )
}

# 3) WEBSITE PNGs (500×500) - safe version
mclapply(mc.cores=16, genes, function(g) {
  mat <- as.array(L_MR_H[g, , ])
  if (max(mat) == min(mat)) mat[] <- 0
  mat <- scales::rescale(mat, to=c(0,1))
  r   <- raster(mat)
  
  out <- sprintf("%s/%s_H.png", web_dir, g)
  safe_png(
    file   = out,
    width  = 500,
    height = 500,
    expr   = print(r %>% website_ggraster() + ggtitle(paste(g, "(Org)"))),
    label  = "WEB"
  )
  NULL
})

# 4) DEEP-LEARNING PNGs (224×224) - safe version
mclapply(mc.cores=16, genes, function(g) {
  mat <- as.array(L_MR_H[g, , ])
  if (max(mat) == min(mat)) mat[] <- 0
  mat <- scales::rescale(mat, to=c(0,1))
  r   <- raster(mat)
  
  out <- sprintf("%s/%s_O.png", ml_dir, g)
  safe_png(
    file   = out,
    width  = 224,
    height = 224,
    expr   = print(r %>% ML_ggraster()),
    label  = "ML"
  )
  NULL
})

cat("Done permutation", perm, "\n")
