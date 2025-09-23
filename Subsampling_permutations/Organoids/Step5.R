#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(parallel)
  library(ggplot2)
  library(png)
  library(raster)
  library(scales)
})

# ---- deps ----
if (!file.exists("lib_misc.R"))  stop("lib_misc.R not found in CWD.")
if (!file.exists("lib_lands.R")) stop("lib_lands.R not found in CWD.")
source("lib_misc.R")   # for website_ggraster()
source("lib_lands.R")  # for ML_ggraster(), etc.

# ---- args ----
# Usage: Rscript Step5_make_pngs.R <perm_id> [lands_dir]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript Step5_make_pngs.R <perm_id> [lands_dir]")
perm <- as.integer(args[1]); if (is.na(perm)) stop("perm_id must be integer")
lands_dir <- if (length(args) >= 2) args[2] else "lands"

# ---- helpers ----
ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Cairo device (seam-free) with error guard
safe_png <- function(file, width, height, expr, label, res = 96) {
  png(filename = file, width = width, height = height, units = "px",
      type = "cairo-png", bg = "white", res = res, antialias = "subpixel")
  on.exit(dev.off(), add = TRUE)
  tryCatch(expr, error = function(e) {
    message(sprintf("[perm %02d][%s] Error plotting %s: %s",
                    perm, label, basename(file), e$message))
  })
}

# parallel policy: override via env (e.g., SLURM sets STEP5_CORES=16)
ncores <- suppressWarnings(as.integer(Sys.getenv("STEP5_CORES", "1")))
if (is.na(ncores) || ncores < 1) ncores <- 1
message(sprintf("[INFO] Step5 using %d core(s) for plotting", ncores))
pllapply <- function(X, FUN) {
  if (ncores > 1) mclapply(X, FUN, mc.cores = ncores, mc.preschedule = FALSE) else lapply(X, FUN)
}

# ---- inputs / outputs ----
infile <- file.path(lands_dir, sprintf("L_MR_H_perm_%02d.rds", perm))
if (!file.exists(infile)) stop(sprintf("No L_MR_H for perm %02d at %s", perm, infile))

L_MR_H <- readRDS(infile)
if (length(dim(L_MR_H)) != 3L) stop("L_MR_H is not a 3D array: ", infile)

web_dir <- file.path(lands_dir, "web_pngLs", "H", sprintf("perm_%02d", perm))
ml_dir  <- file.path(lands_dir, "ml_pngLs",  "H", sprintf("perm_%02d", perm))
ensure_dir(web_dir); ensure_dir(ml_dir)

genes <- rownames(L_MR_H)
if (is.null(genes)) stop("L_MR_H rownames (genes) are NULL.")
total <- length(genes)

# Process in chunks to show progress
chunk_size <- 200L
chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))

# ---- WEBSITE PNGs (500×500) ----
processed <- 0L
for (ch in chunks) {
  pllapply(ch, function(g) {
    mat <- as.array(L_MR_H[g, , , drop = TRUE])
    if (!is.numeric(mat)) mat <- suppressWarnings(matrix(as.numeric(mat), nrow = nrow(mat)))
    if (max(mat, na.rm = TRUE) == min(mat, na.rm = TRUE)) mat[] <- 0
    mat <- scales::rescale(mat, to = c(0, 1))
    r   <- raster::raster(mat, xmn = 0, xmx = ncol(mat), ymn = 0, ymx = nrow(mat))
    out <- file.path(web_dir, sprintf("%s_O.png", g))
    safe_png(out, 500, 500, {
      print(website_ggraster(r) + ggtitle(paste0(g, " (Org)")))
    }, "WEB", res = 96)
    NULL
  })
  processed <- processed + length(ch)
  message(sprintf("[perm %02d] %d/%d (WEB)", perm, processed, total))
}

# ---- ML PNGs (224×224) ----
processed <- 0L
for (ch in chunks) {
  pllapply(ch, function(g) {
    mat <- as.array(L_MR_H[g, , , drop = TRUE])
    if (!is.numeric(mat)) mat <- suppressWarnings(matrix(as.numeric(mat), nrow = nrow(mat)))
    if (max(mat, na.rm = TRUE) == min(mat, na.rm = TRUE)) mat[] <- 0
    mat <- scales::rescale(mat, to = c(0, 1))
    r   <- raster::raster(mat, xmn = 0, xmx = ncol(mat), ymn = 0, ymx = nrow(mat))
    out <- file.path(ml_dir, sprintf("%s_O.png", g))
    safe_png(out, 224, 224, {
      print(ML_ggraster(r))
    }, "ML", res = 96)
    NULL
  })
  processed <- processed + length(ch)
  message(sprintf("[perm %02d] %d/%d (ML)", perm, processed, total))
}

cat(sprintf("[DONE] perm %02d PNGs written to:\n  %s\n  %s\n",
            perm,
            normalizePath(web_dir, winslash = "/"),
            normalizePath(ml_dir,  winslash = "/")))
