library(Seurat)
library(Signac)
library(monocle3)
library(cicero)
library(SeuratDisk)
library(SeuratWrappers)
setwd("~/Humous_reviews")
wang <- readRDS("wang_subset_clean.rds")
wang_multi <- readRDS("snMultiome_atlas_Seurat_object.rds")
wang_only <- subset(wang_multi, cells = colnames(wang))

wang_only@assays$RNA <- NULL
wang_only@assays$SCT <- NULL
wang_only@assays$integrated <- NULL
DefaultAssay(wang_only) <- "ATAC"
saveRDS(wang_only, "wang_peaks.rds")


DefaultAssay(wang_only) <- "ATAC"
# convert to CellDataSet format and make the cicero object
input_cds <- as.cell_data_set(x = wang_only)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

chromosome_length <- read.table("hg38.chrom.sizes")

conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run

#Save results (Optional)
saveRDS(conns, "wang_cicero_connections.Rds")

# Check results
head(conns)
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = "wang_all_peaks.csv")
write.csv(x = conns, file = "wang_cicero_connections.csv")
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(wang_only) <- links

saveRDS(wang_only, 'wang_atac_links.rds')
pbmc <- readRDS("~/Humous_reviews/integrated_all.rds")
DefaultAssay(pbmc) <- "RNA"

VlnPlot(pbmc, features = "IRF1", group.by = "diff_ek")

pbmc$diff_ek <- as.character(pbmc$diff_ek)
pbmc$age_ek <- as.character(pbmc$age_ek)
pbmc$dataset <- as.character(pbmc$dataset)
pbmc$species <- as.character(pbmc$species)
table(pbmc$dataset)
table(pbmc$species)
table(pbmc$age_ek)
Idents(pbmc) <- "species"

pbmc <- subset(pbmc, ident= "H")

pbmc@assays$SCT <- NULL
pbmc@assays$integrated <- NULL
SaveH5Seurat(pbmc, filename = "wang_vb_humous_RNA.h5Seurat")
Convert("wang_vb_humous_RNA.h5Seurat", dest = "h5ad")
