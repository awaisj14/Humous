######## LISI analysis ########
IntH <- readRDS("~/Humous_reviews/IntH.rds")

library(lisi)
library(kBET)
library(dplyr)
library(Seurat)

embeddings <- Embeddings(IntH, reduction = "umap")  # or "pca"
meta <- IntH@meta.data
table(IntH$dataset)
  # 2. Compute LISI
# iLISI -> integration quality (batch mixing)
# cLISI -> conservation quality (cell type preservation)
lisi_res <- compute_lisi(embeddings, meta, c("dataset", "diff_ek"))

# Add results back to Seurat metadata
IntH$lisi_batch   <- lisi_res$dataset
IntH$lisi_celltype <- lisi_res$diff_ek

summary(IntH$lisi_batch)
summary(IntH$lisi_celltype)

# ----------------------------
# 3. Compute kBET
# ----------------------------
# kBET expects numeric matrix of embeddings + batch vector
set.seed(123)
kbet_res_H <- kBET(
  as.matrix(embeddings),
  batch = meta$dataset,
  do.pca = FALSE,
  k0 = 50   # number of neighbors
)

# Mean rejection rate (lower is better)
mean_reject_rate_H <- mean(kbet_res_H$summary$kBET.observed)
mean_reject_rate_H

# ----------------------------
# 4. Summarize
# ----------------------------
cat("iLISI (higher is better):\n")
print(summary(IntH$lisi_batch))

cat("\ncLISI (lower is better):\n")
print(summary(IntH$lisi_celltype))

cat("\nMean kBET rejection rate (lower is better): ", mean_reject_rate, "\n")




IntM <- readRDS("~/Humous_reviews/IntM.rds")

library(lisi)
library(kBET)
library(dplyr)
library(Seurat)

embeddings <- Embeddings(IntM, reduction = "umap")  # or "pca"
meta <- IntM@meta.data
table(IntM$dataset)
# 2. Compute LISI
# iLISI -> integration quality (batch mixing)
# cLISI -> conservation quality (cell type preservation)
lisi_res <- compute_lisi(embeddings, meta, c("dataset", "diff_ek"))

# Add results back to Seurat metadata
IntM$lisi_batch   <- lisi_res$dataset
IntM$lisi_celltype <- lisi_res$diff_ek

summary(IntM$lisi_batch)
summary(IntM$lisi_celltype)

# ----------------------------
# 3. Compute kBET
# ----------------------------
# kBET expects numeric matrix of embeddings + batch vector
set.seed(123)
kbet_res_M <- kBET(
  as.matrix(embeddings),
  batch = meta$dataset,
  do.pca = FALSE,
  k0 = 50   # number of neighbors
)

# Mean rejection rate (lower is better)
mean_reject_rate_M <- mean(kbet_res_M$summary$kBET.observed)
mean_reject_rate_M

# ----------------------------
# 4. Summarize
# ----------------------------
cat("iLISI (higher is better):\n")
print(summary(IntM$lisi_batch))

cat("\ncLISI (lower is better):\n")
print(summary(IntM$lisi_celltype))

cat("\nMean kBET rejection rate (lower is better): ", mean_reject_rate, "\n")



IntO <- readRDS("~/Humous_reviews/IntO.rds")

library(lisi)
library(kBET)
library(dplyr)
library(Seurat)

embeddings <- Embeddings(IntO, reduction = "umap")  # or "pca"
meta <- IntO@meta.data
table(IntO$dataset)
# 2. Compute LISI
# iLISI -> integration quality (batch mixing)
# cLISI -> conservation quality (cell type preservation)
lisi_res <- compute_lisi(embeddings, meta, c("dataset", "diff_ek"))

# Add results back to Seurat metadata
IntO$lisi_batch   <- lisi_res$dataset
IntO$lisi_celltype <- lisi_res$diff_ek

summary(IntO$lisi_batch)
summary(IntO$lisi_celltype)

# ----------------------------
# 3. Compute kBET
# ----------------------------
# kBET expects numeric matrix of embeddings + batch vector
set.seed(123)
kbet_res_O <- kBET(
  as.matrix(embeddings),
  batch = meta$dataset,
  do.pca = FALSE,
  k0 = 50   # number of neighbors
)

# Mean rejection rate (lower is better)
mean_reject_rate_O <- mean(kbet_res_O$summary$kBET.observed)
mean_reject_rate_O

# ----------------------------
# 4. Summarize
# ----------------------------
cat("iLISI (higher is better):\n")
print(summary(IntO$lisi_batch))

cat("\ncLISI (lower is better):\n")
print(summary(IntO$lisi_celltype))

cat("\nMean kBET rejection rate (lower is better): ", mean_reject_rate, "\n")



# ----------------------------
# Export LISI values
# ----------------------------
# For each Seurat object, pull out the relevant metadata
lisi_H <- IntH@meta.data[, c("lisi_batch", "lisi_celltype")]
lisi_O <- IntO@meta.data[, c("lisi_batch", "lisi_celltype")]
lisi_M <- IntM@meta.data[, c("lisi_batch", "lisi_celltype")]

# Write to CSV
write.csv(lisi_H, "IntH_LISI.csv", row.names = TRUE)
write.csv(lisi_O, "IntO_LISI.csv", row.names = TRUE)
write.csv(lisi_M, "IntM_LISI.csv", row.names = TRUE)

# ----------------------------
# Export kBET results
# ----------------------------
# The main useful table is usually in $summary
write.csv(kbet_res$summary,   "kBET_summary_H.csv", row.names = FALSE)
write.csv(kbet_res_O$summary, "kBET_summary_O.csv", row.names = FALSE)
write.csv(kbet_res_M$summary, "kBET_summary_M.csv", row.names = FALSE)

# If you also want the full results (not just summary):
write.csv(kbet_res$results,   "kBET_results_H.csv", row.names = FALSE)
write.csv(kbet_res_O$results, "kBET_results_O.csv", row.names = FALSE)
write.csv(kbet_res_M$results, "kBET_results_M.csv", row.names = FALSE)

IntO
DimPlot(IntO, group.by = "dataset")
