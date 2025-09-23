#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(matrixStats)
  library(ggplot2)
  library(entropy)
  library(dplyr)
  library(readr)
  library(methods)
})

options(stringsAsFactors = FALSE, warn = 1) # show warnings promptly

# ---------------- CONFIG ----------------
get_env <- function(key, default = NULL) {
  val <- Sys.getenv(key, unset = NA_character_)
  if (is.na(val) || val == "") default else val
}

BASE         <- get_env("BASE")
PERM_IDS     <- as.integer(strsplit(get_env("PERM_IDS", "1,2,3,4,5"), ",")[[1]])
USE_AVG      <- as.logical(get_env("USE_AVG", "FALSE"))
REF_LMR_PATH <- get_env("REF_LMR_PATH", "landsO_OG/L_MR_O.rds")

lands_dir <- file.path(BASE, "lands")
out_dir   <- file.path(BASE, "step8", "compare_to_reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] run_dir: ", BASE)
message("[INFO] lands_dir: ", lands_dir)
message("[INFO] out_dir: ", out_dir)
message("[INFO] comparing to ref: ", REF_LMR_PATH)

# ---------------- HELPERS ----------------

fail <- function(...) stop(sprintf(...), call. = FALSE)

# Try to find/load a 3D array from various container types.
load_L3 <- function(path) {
  if (!file.exists(path)) fail("Reference landscape not found: %s", path)
  x <- readRDS(path)
  
  # 1) Already a 3D array
  if (!is.null(dim(x)) && length(dim(x)) == 3) return(x)
  
  # 2) DelayedArray/HDF5-backed
  if (inherits(x, "DelayedArray")) {
    x2 <- as.array(x)
    if (length(dim(x2)) == 3) return(x2)
  }
  
  # 3) List containing a 3D array
  if (is.list(x)) {
    cand <- Filter(function(y) !is.null(dim(y)) && length(dim(y)) == 3, x)
    if (length(cand) >= 1) return(cand[[1]])
  }
  
  # 4) S4 object with a slot containing a 3D array
  if (methods::is(x, "S4")) {
    for (s in methods::slotNames(x)) {
      y <- methods::slot(x, s)
      if (!is.null(dim(y)) && length(dim(y)) == 3) return(y)
    }
  }
  
  fail("Could not find a 3D array in: %s\n[DEBUG] class: %s%s",
       path,
       paste(class(x), collapse = ", "),
       if (!is.null(dim(x))) paste0(" | dim: ", paste(dim(x), collapse = "x")) else "")
}

# Ensure rownames exist; if missing, synthesize consistent names.
ensure_rownames <- function(L3, prefix = "gene") {
  if (is.null(rownames(L3))) {
    rownames(L3) <- paste0(prefix, seq_len(dim(L3)[1]))
  }
  L3
}

# Flatten per-gene SD across [ny, nx] grid
flatten_sd <- function(L3) {
  dims <- dim(L3)
  if (length(dims) != 3L) fail("flatten_sd expects a 3D array")
  flat <- matrix(NA_real_, nrow = dims[1], ncol = dims[2] * dims[3])
  for (i in seq_len(dims[1])) flat[i, ] <- as.vector(L3[i, , ])
  out <- rowSds(flat, na.rm = TRUE)
  names(out) <- rownames(L3)
  out
}

# Per-gene entropy over the ny x nx matrix
calc_entropy_from_L3 <- function(L3) {
  if (length(dim(L3)) != 3L) fail("calc_entropy_from_L3 expects a 3D array")
  vals <- vapply(seq_len(dim(L3)[1]), function(i) {
    m <- as.matrix(L3[i, , ])
    tryCatch(entropy::entropy(m, method = "CS"), error = function(e) NA_real_)
  }, numeric(1))
  names(vals) <- rownames(L3)
  vals
}

# Gene-wise Pearson r between two landscapes (same ny,nx)
per_gene_corr <- function(L3_a, L3_b) {
  if (!identical(dim(L3_a)[2:3], dim(L3_b)[2:3])) {
    fail("Grid mismatch between landscapes: %s vs %s",
         paste(dim(L3_a)[2:3], collapse = "x"),
         paste(dim(L3_b)[2:3], collapse = "x"))
  }
  g <- intersect(rownames(L3_a), rownames(L3_b))
  if (length(g) == 0L) warning("[WARN] No overlapping genes between perm and reference")
  
  res <- data.frame(gene = g, corr = NA_real_)
  for (i in seq_along(g)) {
    a <- as.vector(L3_a[g[i], , ])
    b <- as.vector(L3_b[g[i], , ])
    if (all(is.finite(a)) && all(is.finite(b)) && sd(a) > 0 && sd(b) > 0) {
      res$corr[i] <- suppressWarnings(cor(a, b, method = "pearson"))
    }
  }
  res
}

# Robust loader for perm/avg landscape
load_perm_landscape <- function(tag, lands_dir, perm_ids = NULL) {
  if (identical(tag, "avg")) {
    files <- list.files(lands_dir, pattern = "^L_MR_H_perm_avg.*\\.rds$", full.names = TRUE)
    if (!length(files)) return(NA_character_)
    if (!is.null(perm_ids) && length(files) > 1) {
      pid_str <- paste(sprintf("%02d", perm_ids), collapse = "-")
      files <- grep(pid_str, files, value = TRUE)
    }
    if (length(files) == 0) return(NA_character_)
    if (length(files) > 1) warning("[WARN] Multiple avg files matched, using first: ", files[1])
    return(files[1])
  } else {
    cand <- c(
      file.path(lands_dir, sprintf("L_MR_H_perm_%02d.rds", as.integer(tag))),
      file.path(lands_dir, sprintf("L_MR_H_perm_%d.rds",  as.integer(tag)))
    )
    f <- cand[file.exists(cand)]
    if (length(f)) return(f[1]) else return(NA_character_)
  }
}

# ---------------- LOAD REFERENCE ----------------
L_ref <- load_L3(REF_LMR_PATH)
L_ref <- ensure_rownames(L_ref, prefix = "gene")
if (length(dim(L_ref)) != 3L) fail("Reference is not a 3D array after loading.")
message("[INFO] Reference dims: ", paste(dim(L_ref), collapse = " x "))

sd_ref      <- flatten_sd(L_ref)
entropy_ref <- calc_entropy_from_L3(L_ref)

# ---------------- WHICH PERMS ----------------
targets <- if (isTRUE(USE_AVG)) "avg" else sprintf("%02d", PERM_IDS)
all_summaries <- list()
any_outputs <- FALSE

# ---------------- PROCESS EACH PERM ----------------
for (tag in targets) {
  L_file <- load_perm_landscape(tag, lands_dir, if (identical(tag, "avg")) PERM_IDS else NULL)
  if (is.na(L_file) || !file.exists(L_file)) {
    warning(sprintf("[SKIP] missing: %s",
                    ifelse(is.na(L_file), paste0("avg file in ", lands_dir), L_file)))
    next
  }
  
  obj <- readRDS(L_file)
  if (!is.null(dim(obj)) && length(dim(obj)) == 3) {
    L_perm <- obj
  } else if (inherits(obj, "DelayedArray")) {
    L_perm <- as.array(obj)
  } else if (is.list(obj)) {
    cand <- Filter(function(y) !is.null(dim(y)) && length(dim(y)) == 3, obj)
    if (length(cand) >= 1) {
      L_perm <- cand[[1]]
    } else {
      warning("[SKIP] no 3D array found in ", L_file); next
    }
  } else if (methods::is(obj, "S4")) {
    L_perm <- NULL
    for (s in methods::slotNames(obj)) {
      y <- methods::slot(obj, s)
      if (!is.null(dim(y)) && length(dim(y)) == 3) { L_perm <- y; break }
    }
    if (is.null(L_perm)) { warning("[SKIP] no 3D array slot found in ", L_file); next }
  } else {
    warning("[SKIP] unsupported object type in ", L_file, " (class: ", paste(class(obj), collapse=", "), ")")
    next
  }
  
  L_perm <- ensure_rownames(L_perm, prefix = "gene")
  if (length(dim(L_perm)) != 3L) { warning("[SKIP] perm is not 3D in ", L_file); next }
  
  message("[INFO] Perm ", tag, " dims: ", paste(dim(L_perm), collapse = " x "))
  
  if (!identical(dim(L_perm)[2:3], dim(L_ref)[2:3])) {
    warning("[SKIP] grid mismatch for tag ", tag, ": perm ",
            paste(dim(L_perm)[2:3], collapse="x"), " vs ref ",
            paste(dim(L_ref)[2:3], collapse="x"))
    next
  }
  
  corr_df <- per_gene_corr(L_perm, L_ref)
  sd_perm      <- flatten_sd(L_perm)
  entropy_perm <- calc_entropy_from_L3(L_perm)
  
  df <- corr_df %>%
    mutate(
      sd_perm  = unname(sd_perm[gene]),
      sd_ref   = unname(sd_ref[gene]),
      ent_perm = unname(entropy_perm[gene]),
      ent_ref  = unname(entropy_ref[gene])
    )
  
  csv_path <- file.path(out_dir, sprintf("perm_%s_vs_ref_per_gene.csv", tag))
  write_csv(df, csv_path)
  message("[OK] per-gene table -> ", csv_path)
  
  finite_corr <- df$corr[is.finite(df$corr)]
  if (length(finite_corr) >= 1L) {
    p <- ggplot(df, aes(corr)) +
      geom_histogram(bins = 50, color = "grey30", fill = "steelblue") +
      theme_bw() +
      labs(title = paste("Perm", tag, "vs Reference: per-gene Pearson r"),
           x = "Pearson r", y = "Genes")
    pdf_path <- file.path(out_dir, sprintf("perm_%s_vs_ref_corr_hist.pdf", tag))
    pdf(pdf_path, width = 6, height = 4)
    print(p)
    invisible(dev.off())
    message("[OK] histogram -> ", pdf_path)
  } else {
    warning("[NOTE] No finite correlations for tag ", tag, "; histogram skipped.")
  }
  
  s <- data.frame(
    target = tag,
    n_genes = sum(!is.na(df$corr)),
    corr_median = suppressWarnings(median(df$corr, na.rm = TRUE)),
    corr_mean   = suppressWarnings(mean(df$corr, na.rm = TRUE)),
    corr_q05    = suppressWarnings(quantile(df$corr, 0.05, na.rm = TRUE, names = FALSE)),
    corr_q95    = suppressWarnings(quantile(df$corr, 0.95, na.rm = TRUE, names = FALSE)),
    sd_perm_median  = suppressWarnings(median(df$sd_perm,  na.rm = TRUE)),
    sd_ref_median   = suppressWarnings(median(df$sd_ref,   na.rm = TRUE)),
    ent_perm_median = suppressWarnings(median(df$ent_perm, na.rm = TRUE)),
    ent_ref_median  = suppressWarnings(median(df$ent_ref,  na.rm = TRUE))
  )
  all_summaries[[length(all_summaries) + 1]] <- s
  any_outputs <- TRUE
}

# ---------------- FINAL SUMMARY ----------------
if (length(all_summaries)) {
  summ <- bind_rows(all_summaries)
  summ_path <- file.path(out_dir, "summary_perm_vs_ref.csv")
  write_csv(summ, summ_path)
  message("[OK] summary -> ", summ_path)
  suppressWarnings(print(summ))
} else if (isTRUE(any_outputs)) {
  message("[WARN] Per-gene outputs existed but no summaries were produced.")
} else {
  message("[WARN] No comparisons were produced (no valid perm files matched or grids mismatched).")
}

# ---------------- SESSION INFO ----------------
message("[INFO] Session info (key pkgs):")
message("matrixStats: ", as.character(packageVersion("matrixStats")),
        " | entropy: ",   as.character(packageVersion("entropy")),
        " | ggplot2: ",   as.character(packageVersion("ggplot2")),
        " | dplyr: ",     as.character(packageVersion("dplyr")),
        " | readr: ",     as.character(packageVersion("readr")))
