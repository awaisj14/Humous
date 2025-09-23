library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)
library(scDblFinder)
library(tidyr)
library(tidyverse)
library(qs)
library(SingleCellExperiment)
library(tidyverse)

c25 <- c(
  "gray", "black", # red
  
  "red",
  "pink", # lt pink
  "green",
  "blue", # lt purple
  "magenta", # lt orange
  "yellow", "darkturquoise",
  "maroon", "orchid", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
color.predict = c("gray", 
                  "black"
) 




color.predict_spatial = c("green4", 
                          "palegreen3",
                          "#20b2aa"
) 

packages <- c("Seurat", "ggplot2", "tidyr", "readr"
              , "plyr", "stringr", "harmony"
              , "cowplot", "reshape2","ggpubr"
              , "gsubfn", "tibble", "gplots"
              , "Matrix", "dplyr", "pbapply", "schex"
              , "UpSetR", "extraDistr", "ape"
              , "stats"
)
lapply(packages, library, character.only = TRUE)

options(Seurat.object.assay.version = 'v3')






####### dataset #######
data <- Read10X_h5('/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/Control/cellbender_output_50epo_filtered_seurat.h5')
metadata <- read.csv(
  file = "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/Control/cellbender_output_50epo_metrics.csv",
  header = TRUE
)
LF1 <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = "Control",meta.data = metadata
)

data <- Read10X_h5('/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/hEGR1/cellbender_output_50epo_filtered_seurat.h5')
metadata <- read.csv(
  file = "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/hEGR1/cellbender_output_50epo_metrics.csv",
  header = TRUE
)
LF2 <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = "hEGR1",meta.data = metadata
)

data <- Read10X_h5('/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/hIRF1/cellbender_output_50epo_filtered_seurat.h5')
metadata <- read.csv(
  file = "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/hIRF1/cellbender_output_50epo_metrics.csv",
  header = TRUE
)
LF3 <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = "hIRF1",meta.data = metadata
)

data <- Read10X_h5('/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/hIRF1_hEGR1/cellbender_output_50epo_filtered_seurat.h5')
metadata <- read.csv(
  file = "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/hIRF1_hEGR1/cellbender_output_50epo_metrics.csv",
  header = TRUE
)
LF4 <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = "hIRF1_hEGR1",meta.data = metadata
)
pbmc <- merge(LF1, y= c(LF2,LF3, LF4))

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident", ncol = 3)
pbmc <- subset(x = pbmc,subset = nCount_RNA > 500 & nCount_RNA < 20000 & nFeature_RNA >1000)
pbmc <- subset(x = pbmc,subset = percent.mt < 5)
pbmc <- SCTransform(pbmc)

DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc$clusters <- pbmc$seurat_clusters
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters")
DimPlot(pbmc, reduction = "umap", group.by = "clusters", split.by="orig.ident")
FeaturePlot(pbmc, features = "Junb")
FeaturePlot(pbmc, features = "Neurod2")
FeaturePlot(pbmc, features = "Eomes")
VlnPlot(pbmc, features = "Junb", group.by = "orig.ident")
tenX_matrix <- "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/matrix"
out <- "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/out_dblfinder"
print(paste0("Using the following counts directory: ", tenX_matrix))
## make sure the directory exists ###
dir.create(out, recursive = TRUE)

## Read in Data ##

sce <- as.SingleCellExperiment(pbmc)
doublet_ratio <- ncol(sce)/1000*0.008
sce <- scDblFinder(sce, dbr=doublet_ratio)

### Make a dataframe of the results ###
results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)


write_delim(results, path = paste0(out,"/scDblFinder_doublets_singlets.tsv"), delim = "\t")

### Calculate number of doublets and singlets ###
summary <- as.data.frame(table(results$scDblFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
write_delim(summary, paste0(out,"/scDblFinder_doublet_summary.tsv"), "\t")
pbmc <- AddMetaData(pbmc, results$scDblFinder_DropletType, col.name = 'scDblFinder_DropletType')
table(pbmc$scDblFinder_DropletType)
Idents(pbmc) <- "scDblFinder_DropletType"
pbmc <- subset(pbmc, ident = "singlet")
DimPlot(pbmc, reduction = "umap", group.by = "scDblFinder_DropletType")

saveRDS(pbmc, "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/merged_doublet.rds")
pbmc <- readRDS("/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/merged_doublet.rds")


Idents(pbmc) <- "seurat_clusters"
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# Filter markers with avg_log2FC > 1 and take top 100 per cluster
top10 <- pbmc.markers %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()


# Reorder genes accordingly
top10 <- top10 %>% arrange(cluster, desc(avg_log2FC))

DimPlot(pbmc, group.by = "seurat_clusters", cols = c25, label = T)
FeaturePlot(pbmc, features = "Eomes")
FeaturePlot(pbmc, features = "Neurog2")

Idents(pbmc) <- "seurat_clusters"
pbmc <- RenameIdents(
  object = pbmc,
  c('0' = 'RG','1' = 'IP','2' = 'IP','3' = 'RG',
    '4' = 'IP', "5" = "N","6" = "N","7" = "N","8" = "IP","9" = "IP","10" = "IP","11" = "RG",
    "12" = "RG","13" = "RG","14" = "Microglia","15" = "Vasculature","16" = "N","17" = "N","18" = "N"
    
  ))

pbmc$celltype <- pbmc@active.ident
DimPlot(pbmc, group.by = "celltype", cols = c25)
FeaturePlot(pbmc, features = "hEGR1")
FeaturePlot(pbmc, features = "Junb")
VlnPlot(pbmc, split.by = "celltype", group.by = "orig.ident",features = "Junb")

DefaultAssay(pbmc) <- "SCT"
pbmc@reductions$umap_rna <- pbmc@reductions$umap
pbmc@reductions$pca_rna <- pbmc@reductions$pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, group.by = "celltype", cols = c25)
table(pbmc$celltype)
saveRDS(pbmc, "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/merged_all_annot.rds")
pbmc <- readRDS("/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/merged_all_annot.rds")
Idents(pbmc) <- "celltype"
pbmc <- subset(pbmc, ident = c("Microglia","Vasculature"), invert = T)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- RunUMAP(pbmc, dims = 1:5, n.neighbors = 200)
pbmc <- FindClusters(pbmc, resolution = 0.3)

DimPlot(pbmc, group.by = "seurat_clusters", cols = c25, label = T)
FeaturePlot(pbmc, features = "Neurod2")

Idents(pbmc) <- "seurat_clusters"
pbmc <- RenameIdents(
  object = pbmc,
  c('0' = 'RG','1' = 'IP','2' = 'N','3' = 'IP',
    '4' = 'IP', "5" = "RG",'6' = 'N','7' = 'RG','8' = 'N','9' = 'IP',
    '10' = 'RG'
    
  ))
pbmc$celltype <- pbmc@active.ident

VlnPlot(pbmc, features = "nCount_RNA")
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", cols = c25)
DimPlot(pbmc, reduction = "umap", group.by = "celltype", cols = c25)



####### Plotting #######

library(Seurat)
library(dplyr)
library(ggplot2)

# Get metadata
meta <- pbmc@meta.data
# Set desired order for celltype
meta$celltype <- factor(meta$celltype, levels = c("RG", "IP", "N"))
# Summarize cell counts by celltype and orig.ident
cell_counts <- meta %>%
  dplyr::group_by(orig.ident, celltype) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()

# Optional: Calculate proportion within each orig.ident
cell_props <- cell_counts %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup()

# Define custom colors
celltype_colors <- c(
  "RG" = "#e2c683",
  "IP" = "#86c9c2",
  "N"  = "#60348a"
)

# Plot with custom colors
ggplot(cell_props, aes(x = orig.ident, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  ylab("Proportion of Cells") +
  xlab("Sample (orig.ident)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Cell Type")

library(gprofiler2)
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
pbmc <- CellCycleScoring(pbmc, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = TRUE)


library(dplyr)
library(ggplot2)

# Subset metadata to only RG cells
meta <- pbmc@meta.data
rg_meta <- meta %>%
  dplyr::filter(celltype == "RG") %>%
  dplyr::mutate(Phase = factor(Phase, levels = c("G1", "S", "G2M")))  # Optional: set order

# Count RG cells per orig.ident and Phase
rg_counts <- rg_meta %>%
  dplyr::group_by(orig.ident, Phase) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()

# Calculate proportions within each orig.ident
rg_props <- rg_counts %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup()

# Plot: stacked barplot
ggplot(rg_props, aes(x = orig.ident, y = proportion, fill = Phase)) +
  geom_bar(stat = "identity") +
  ylab("Proportion of RG Cells") +
  xlab("Sample (orig.ident)") +
  scale_fill_brewer(palette = "Set2") +  # Or use your own color vector
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Cell Cycle Phase", title = "RG Cell Cycle Phase Distribution by Sample")


# Create a combined celltype + orig.ident label
pbmc$celltype_orig <- paste(pbmc$celltype, pbmc$orig.ident, sep = "_")

# DotPlot across combined groups
DotPlot(pbmc, features = c("Junb","Cebpd","Glis3","Irf1"), group.by = "celltype_orig") +
  scale_color_gradient(low = "white", high = "#0e518f") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Junb Expression across Cell Type and Sample")




######## Plotting 2 #########
library(Seurat)
library(dplyr)
library(ggplot2)

# Create a new column in pbmc@meta.data with the combined label
pbmc$custom_label <- with(pbmc@meta.data, case_when(
  celltype == "RG" & Phase == "G1"  ~ "RG_G1",
  celltype == "RG" & Phase == "S"   ~ "RG_S",
  celltype == "RG" & Phase == "G2M" ~ "RG_G2M",
  celltype == "RG"                  ~ "RG_other",
  celltype == "IP" & Phase == "G1"  ~ "IP_G1",
  celltype == "IP" & Phase %in% c("S", "G2M") ~ "IP_SG2M",
  celltype == "IP"                  ~ "IP_other",
  celltype == "N"                   ~ "N",
  TRUE                              ~ "Other"  # fallback, should be rare
))

# Set custom order if you like (optional)
pbmc$custom_label <- factor(pbmc$custom_label, levels = c(
  "RG_G1", "RG_S", "RG_G2M", "RG_other",
  "IP_G1", "IP_SG2M", "IP_other",
  "N"
))

# Define the colors
custom_colors <- c(
  "RG_G1"     = "#cba350",
  "RG_S"      = "#fbd487",
  "RG_G2M"    = "#f9bb4d",
  "RG_other"  = "#e5c373",
  "IP_G1"     = "#59a299",
  "IP_SG2M"   = "#76cbc5",
  "IP_other"  = "#76cbc5",
  "N"         = "#b899c8"
)

# Plot
DimPlot(pbmc, group.by = "custom_label", cols = custom_colors, pt.size = 1)+NoLegend()
table(pbmc$custom_label)
FeaturePlot(pbmc, features = c("Neurog2"))
DimPlot(pbmc, group.by = "seurat_clusters")
ref <- readRDS("/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Final_030923_Junb.rds")

library(Seurat)
library(ggplot2)

# Define your custom colors: Control = pale gray, others = muted reds
custom_colors <- c(
  "Control" = "#e5e8e6",       # pale gray
  "hEGR1" = "#f4a6a6",         # light muted red
  "hIRF1" = "#e06666",         # medium muted red
  "hIRF1_hEGR1" = "#990000"    # deep muted red (darkest)
)

# Ensure orig.ident is ordered correctly
pbmc$orig.ident <- factor(pbmc$orig.ident,
                          levels = c("Control", "hEGR1", "hIRF1", "hIRF1_hEGR1"))

# Generate UMAP plot
DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", shuffle = T, pt.size = 1) +
  scale_color_manual(values = custom_colors) + NoLegend()



######## RG specific plotting ########
saveRDS(pbmc, "/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/merged_all_annot_final.rds")
pbmc <- readRDS("/Users/javed/Documents/Humous_indexing/Analysis_reviewers/scRNA/Custom_genome/merged_all_annot_final.rds")


Idents(pbmc) <- "celltype"
pbmc <- subset(pbmc, ident = "RG")

DefaultAssay(pbmc) <- "RNA"
DotPlot(pbmc, features = c("Hopx"), group.by = "celltype_orig") +
  scale_color_gradient(low = "white", high = "#0e518f") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Junb Expression across Cell Type and Sample")


VlnPlot(pbmc, features = "Hopx", group.by = "orig.ident")


########################################
# Libraries
########################################
library(Seurat)
library(Matrix)
library(dplyr)
library(readr)
library(VennDiagram)
library(grid)  # required for grid.draw()

########################################
# Load data
########################################
junbOE_df <- read_csv("DEGs_RG_markers_Junb_ek.csv")
pbmc$celltype
Idents(pbmc) <- "celltype"
pbmc <- subset(pbmc, ident = "RG")
Idents(pbmc) <- "orig.ident"
EGR1 <- FindMarkers(pbmc, ident.1 = "hEGR1",ident.2 = "Control")
IRF1 <- FindMarkers(pbmc, ident.1 = "hIRF1",ident.2 = "Control")
IRF1EGR1 <- FindMarkers(pbmc, ident.1 = "hIRF1_hEGR1",ident.2 = "Control")
# Make sure gene names are a proper column
IRF1$gene <- rownames(IRF1)
EGR1$gene <- rownames(EGR1)
IRF1EGR1$gene <- rownames(IRF1EGR1)
library(VennDiagram)
library(grid)
library(gridExtra)

# --- Define Junb subsets (JunbOE) ---
junb_repressed <- unique(subset(junbOE_df, avg_log2FC >=  0.25)$gene)  # higher in Control (repressed in Junb)
junb_activated <- unique(subset(junbOE_df, avg_log2FC <= -0.25)$gene)  # higher in Junb   (activated in Junb)

# tiny helper to turn VennDiagram's return into a grob we can arrange
venn_grob <- function(x_named_list, fills) {
  v <- venn.diagram(
    x = x_named_list,
    filename = NULL,                  # return grob(s), don't write file
    fill = fills,
    alpha = 0.6,
    cex = 1.2,
    cat.cex = 1.2,
    cat.pos = 0,
    cat.dist = 0.05,
    margin = 0.1
  )
  grid::grobTree(v)
}

compare_with_junb <- function(df, label, out_prefix) {
  # safety: enforce expected columns + character vectors
  stopifnot(all(c("gene","avg_log2FC") %in% colnames(df)))
  df$gene <- as.character(df$gene)
  
  df_pos <- unique(subset(df, avg_log2FC >  0.5)$gene)
  df_neg <- unique(subset(df, avg_log2FC < -0.5)$gene)
  
  # overlaps
  overlap_pos_with_junb_activated <- intersect(df_pos, junb_activated)
  overlap_neg_with_junb_repressed <- intersect(df_neg, junb_repressed)
  
  cat("\n====", label, "====\n")
  cat("Positive genes:", length(df_pos), "\n")
  cat("Negative genes:", length(df_neg), "\n")
  cat("Overlap (positive ∩ Junb-activated):", length(overlap_pos_with_junb_activated), "\n")
  cat("Overlap (negative ∩ Junb-repressed):", length(overlap_neg_with_junb_repressed), "\n")
  
  # save overlap gene lists
  write.csv(data.frame(gene = overlap_pos_with_junb_activated),
            paste0(out_prefix, "_pos_vs_Junb_activated.csv"), row.names = FALSE)
  write.csv(data.frame(gene = overlap_neg_with_junb_repressed),
            paste0(out_prefix, "_neg_vs_Junb_repressed.csv"), row.names = FALSE)
  
  # Colors (your scheme): light gray, blackish
  venn_colors <- c("#E6E7E8", "#231F1F")
  
  # build **named** lists with setNames()
  x_pos <- setNames(
    list(df_pos, junb_activated),
    c(paste0(label, "_positive"), "Junb_activated")
  )
  x_neg <- setNames(
    list(df_neg, junb_repressed),
    c(paste0(label, "_negative"), "Junb_repressed")
  )
  
  venn1 <- venn_grob(x_pos, fills = venn_colors)
  venn2 <- venn_grob(x_neg, fills = venn_colors)
  
  # combined PDF per TF
  pdf(paste0(out_prefix, "_venn.pdf"), width = 10, height = 5)
  gridExtra::grid.arrange(venn1, venn2, ncol = 2,
                          top = paste("Overlap of", label, "with Junb-regulated genes"))
  dev.off()
  cat("Saved Venns →", paste0(out_prefix, "_venn.pdf"), "\n")
  
  invisible(list(
    pos_genes = df_pos,
    neg_genes = df_neg,
    overlap_pos_with_junb_activated = overlap_pos_with_junb_activated,
    overlap_neg_with_junb_repressed = overlap_neg_with_junb_repressed
  ))
}

# Run comparisons
out_irf1      <- compare_with_junb(IRF1,      "IRF1",      "IRF1")
out_egr1      <- compare_with_junb(EGR1,      "EGR1",      "EGR1")
out_irf1egr1  <- compare_with_junb(IRF1EGR1,  "IRF1EGR1",  "IRF1EGR1")
