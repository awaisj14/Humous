#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ----------------------------
# Settings
# ----------------------------
base_out   <- "Mouse_models"                         
n_perms    <- 5
models_dir <- file.path(base_out, "compare_out")     
dir.create(models_dir, showWarnings = FALSE)

# Input files
orig_age_f   <- file.path(base_out, "Original_humous", "model_fullAGE_M.rds")
merged_rds   <- file.path(base_out, "merged_all.rds")   # <-- merged object you used for training

# ----------------------------
# Load gene names
# ----------------------------
merged_all <- readRDS(merged_rds)
gene_names <- rownames(merged_all@assays$integrated@scale.data)

# ----------------------------
# Helper to extract coefficients
# ----------------------------
extract_coef <- function(mod, genes) {
  if (is.list(mod) && !is.null(mod$coefficients)) {
    coef <- mod$coefficients
  } else if (is.atomic(mod)) {
    coef <- mod
  } else {
    coef <- suppressWarnings(coef(mod))
  }
  
  # Always enforce gene names
  names(coef) <- genes[seq_along(coef)]
  coef
}

# ----------------------------
# Collect coefficients from permutations
# ----------------------------
collect_perm_coefs <- function(prefix, genes) {
  all_coefs <- list()
  for (i in seq_len(n_perms)) {
    i_padded <- sprintf("%02d", i)
    f <- file.path(base_out, paste0(prefix, "_perm", i_padded, ".rds"))
    if (!file.exists(f)) stop("Missing file: ", f)
    mod <- readRDS(f)
    coef <- extract_coef(mod, genes)
    all_coefs[[i]] <- coef
  }
  
  mat <- do.call(cbind, lapply(all_coefs, function(x) x[genes]))
  colnames(mat) <- sprintf("perm_%02d", seq_len(n_perms))
  
  coef_df <- data.frame(
    gene = genes,
    mean_perm = rowMeans(mat, na.rm = TRUE),
    sd_perm   = apply(mat, 1, sd, na.rm = TRUE),
    mat,
    stringsAsFactors = FALSE
  )
  coef_df
}

# ----------------------------
# Compare AGE permutations vs original
# ----------------------------
coef_df_age <- collect_perm_coefs("fullAGE", gene_names)

# load reference
ref_mod_age <- readRDS(orig_age_f)
coef_ref_age <- extract_coef(ref_mod_age, gene_names)

ref_df <- data.frame(
  gene = names(coef_ref_age),
  ref  = as.numeric(coef_ref_age),
  stringsAsFactors = FALSE
)

plot_df <- inner_join(ref_df, coef_df_age, by = "gene")

# correlations
pearson_r  <- cor(plot_df$ref, plot_df$mean_perm, use = "complete.obs", method = "pearson")
spearman_r <- cor(plot_df$ref, plot_df$mean_perm, use = "complete.obs", method = "spearman")

# save CSV
out_csv <- file.path(models_dir, "coef_fullAGE_perms.csv")
write_csv(plot_df, out_csv)

# scatterplot
p <- ggplot(plot_df, aes(x = ref, y = mean_perm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(
    x = "Original fullAGE coeff",
    y = "Mean perm fullAGE coeff",
    title = "Original vs Mean Perm Coeffs (AGE)"
  ) +
  annotate(
    "text",
    x = min(plot_df$ref, na.rm = TRUE),
    y = max(plot_df$mean_perm, na.rm = TRUE),
    label = sprintf("Pearson=%.2f\nSpearman=%.2f", pearson_r, spearman_r),
    hjust = 0,
    vjust = 1
  )

out_pdf <- file.path(models_dir, "coef_fullAGE_comparison.pdf")
pdf(out_pdf, width = 6, height = 6)
print(p)
dev.off()

message("✅ Comparison complete with real gene names. CSV + PDF saved in: ", models_dir)



# Function to extract and rank gene weights from bmrm::nrbm models
extract_gene_weights <- function(model, top_n = 50) {
  # Extract coefficients and gene names
  coefs <- as.numeric(model)
  genes <- attr(model, "gradient")
  
  # Make sure lengths match
  if (length(coefs) != length(genes)) {
    stop("Mismatch between number of coefficients and genes")
  }
  
  # Build dataframe
  df <- data.frame(
    gene = genes,
    weight = coefs,
    abs_weight = abs(coefs)
  )
  
  # Rank by absolute weight (importance)
  df <- df[order(-df$abs_weight), ]
  
  # Optionally return only top_n
  if (!is.null(top_n)) {
    df <- head(df, top_n)
  }
  
  return(df)
}

# Example usage
model <- readRDS("Mouse_models/fullDIFF_perm05.rds")
gene_weights <- extract_gene_weights(model, top_n = 100)

# View the top weighted genes
head(gene_weights, 20)
