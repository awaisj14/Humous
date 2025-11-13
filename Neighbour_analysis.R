library(Seurat)
library(RANN)
library(dplyr)

# ----------------------------
# 1. Filter IntM by age_ek
# ----------------------------
IntM_filtered <- subset(IntM, !(age_ek %in% c("12")))
meta <- IntM_filtered@meta.data
meta$cond <- trimws(meta$cond)

# ----------------------------
# 2. Extract PCA embeddings
# ----------------------------
embeddings <- Embeddings(IntM_filtered, reduction = "pca")[, 1:30]

# ----------------------------
# 3. Compute nearest neighbors
# ----------------------------
k <- 50
nn <- nn2(embeddings, k = k)
rownames(nn$nn.idx) <- rownames(embeddings)

# ----------------------------
# 4. Count Control vs Hu_Ctl neighbors
# ----------------------------
control_types <- c("Control", "Hu_Ctl")

count_neighbors_by_type <- function(cell) {
  neighbors <- nn$nn.idx[cell, ]
  neighbor_conds <- meta[neighbors, "cond"]
  sapply(control_types, function(ct) sum(neighbor_conds == ct))
}

neighbor_counts_split <- t(sapply(rownames(meta), count_neighbors_by_type))
neighbor_counts_split <- as.data.frame(neighbor_counts_split)
colnames(neighbor_counts_split) <- c("neighbors_Control", "neighbors_Hu_Ctl")
neighbor_counts_split$cell_id <- rownames(neighbor_counts_split)

# ----------------------------
# 5. Merge into metadata
# ----------------------------
meta$cell_id <- rownames(meta)
meta <- meta[, !grepl("^neighbors_", colnames(meta))]  # drop old neighbor cols
meta <- left_join(meta, neighbor_counts_split, by = "cell_id")
rownames(meta) <- meta$cell_id
meta$cell_id <- NULL
IntM_filtered@meta.data <- meta

# ----------------------------
# 6. Per-cell results
# ----------------------------
per_cell_neighbors <- meta %>%
  dplyr::select(cond, neighbors_Control, neighbors_Hu_Ctl) %>%
  dplyr::mutate(
    total_neighbors = neighbors_Control + neighbors_Hu_Ctl,
    prop_Control = neighbors_Control / k,
    prop_Hu_Ctl  = neighbors_Hu_Ctl / k,
    closer_to = ifelse(prop_Control > prop_Hu_Ctl, "Control", "Hu_Ctl"),
    prop_Control_vs_HuCtl = ifelse(
      (neighbors_Control + neighbors_Hu_Ctl) > 0,
      neighbors_Control / (neighbors_Control + neighbors_Hu_Ctl),
      NA_real_
    )
  )

# Export per-cell CSV
write.csv(per_cell_neighbors, "PCA_noE12_neighbors_per_cell.csv", row.names = TRUE)

# ----------------------------
# 7. Summaries
# ----------------------------
groups_of_interest <- c("Control", "hEGR1", "hIRF1", "hIRF1_hEGR1", "Junb")

## (a) Global proportions
summary_table_global <- per_cell_neighbors %>%
  dplyr::filter(cond %in% groups_of_interest) %>%
  dplyr::group_by(cond) %>%
  dplyr::summarise(
    mean_prop_Control = mean(prop_Control),
    mean_prop_Hu_Ctl  = mean(prop_Hu_Ctl),
    n_CloserToControl = sum(closer_to == "Control"),
    n_CloserToHuCtl   = sum(closer_to == "Hu_Ctl"),
    .groups = "drop"
  )

print(summary_table_global)

## (b) Normalized Control-vs-HuCtl only
summary_table_vs <- per_cell_neighbors %>%
  dplyr::filter(cond %in% groups_of_interest) %>%
  dplyr::group_by(cond) %>%
  dplyr::summarise(
    mean_prop_Control_vs_HuCtl = mean(prop_Control_vs_HuCtl, na.rm = TRUE),
    n_CloserToControl = sum(prop_Control_vs_HuCtl > 0.5, na.rm = TRUE),
    n_CloserToHuCtl   = sum(prop_Control_vs_HuCtl < 0.5, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table_vs)

# ----------------------------
# 8. Optional: visualize PCA
# ----------------------------
DimPlot(IntM_filtered, group.by = "cond", reduction = "pca")

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggforce) # for stat_ellipse with better control

# Get PCA embeddings and metadata
pca_df <- as.data.frame(Embeddings(IntM_filtered, reduction = "pca")[, 1:2])
pca_df$cond <- IntM_filtered@meta.data$cond

centroids <- pca_df %>%
  dplyr::group_by(cond) %>%
  dplyr::summarise(PC1 = mean(PC_1), PC2 = mean(PC_2))


ggplot(pca_df, aes(x = PC_1, y = PC_2, color = cond)) +
  geom_point(alpha = 0.5, size = 1) +
  stat_ellipse(aes(group = cond), type = "norm", level = 0.68, size = 1) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = cond), 
             shape = 8, size = 3, stroke = 1.5) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA Plot with Group Ellipses", 
       x = "PC1", 
       y = "PC2") +
  theme(legend.title = element_blank())



