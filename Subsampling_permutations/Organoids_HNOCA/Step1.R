library(Seurat)
library(reticulate)
library(anndata)
library(plyr)
data <- read_h5ad("HNOCA_subset.h5ad") # Obtained from https://cellxgene.cziscience.com/collections/de379e5f-52d0-498c-9801-0f850823c847 1.7M cells non-extended set
obs <- data$obs
counts <- data$X

pbmc <- CreateSeuratObject(counts = t(counts), meta.data = obs)
table(pbmc$age_ek)

table(pbmc$annot_level_2)

pbmc$diff_ek <- pbmc$annot_level_2

pbmc@meta.data$diff_ek <- mapvalues(
  pbmc@meta.data$diff_ek,
  from = c(
    "Dorsal Telencephalic IP",
    "Dorsal Telencephalic NPC",
    "Dorsal Telencephalic Neuron"
  ),
  to   = c("IPC", "RG", "N")
)

pbmc@meta.data$diff_ek <- factor(
  pbmc@meta.data$diff_ek,
  levels = c("IPC", "RG", "N", setdiff(levels(pbmc@meta.data$diff_ek), c("IPC","RG","N")))
)

table(pbmc$diff_ek)

pbmc$assay_differentiation
Idents(pbmc) <- "assay_differentiation"
pbmc <- subset(pbmc, idents = c("Pellegrini, 2020 (doi: 10.1126/science.aaz5626); hChPO","Quadrato, 2023 (doi: no_doi)","Sawada, 2020 (doi: 10.1038/s41380-020-0844-z)",
                                "Xiang, 2019 (doi: 10.1016/j.stem.2018.12.015)","Qian et al., 2016 (doi: 10.1016/j.cell.2016.04.032)","Jo et al., 2016 (doi: 10.1016/j.stem.2016.07.005)",
                                "Huang, 2021 (doi: 10.1016/j.stem.2021.04.006)","Fiorenzano, 2021 (doi: 10.1038/s41467-021-27464-5); standard","Fiorenzano, 2021 (doi: 10.1038/s41467-021-27464-5); silk",
                                "Fiorenzano, 2021 (doi:10.1038/s41467-021-27464-5); silk+laminin","Birey, 2017 (doi: 10.1038/nature22330)","Andersen, 2020 (doi: 10.1016/j.cell.2020.11.017); without DAPT",
                                "Andersen, 2020 (doi: 10.1016/j.cell.2020.11.017)"),invert = T)

pbmc$dataset_proto <- pbmc$assay_differentiation

pbmc@meta.data$dataset_proto <- mapvalues(
  pbmc@meta.data$dataset_proto,
  from = c(
    "Bhaduri, 2020 (doi: 10.1038/s41586-020-1962-0); directed",
    "Bhaduri, 2020 (doi: 10.1038/s41586-020-1962-0); most directed",
    "Esk, 2020 (doi: 10.1126/science.abb5390)",
    "Lancaster, 2014 (doi: 10.1038/nprot.2014.158)",
    "Miura, 2020 (doi: 10.1038/s41587-020-00763-w)",
    "Pasca, 2015 (doi: 10.1038/nmeth.3415)",
    "Pellegrini, 2020 (doi: 10.1126/science.aaz5626); hCO",
    "Qian, 2020 (doi: 10.1016/j.stem.2020.02.002)",
    "Quadrato, 2017 (doi: 10.1038/protex.2017.049)",
    "Trujillo, 2019 (doi: 10.1016/j.stem.2019.08.002)",
    "Velasco, 2019 (doi: 10.1038/s41586-019-1289-x)",
    "Watanabe, 2017 (doi: 10.1016/j.celrep.2017.09.047)",
    "Yoon, 2019 (doi: 10.1038/s41592-018-0255-0)"
    
  ),
  to   = c("Kb20_1", "Kb20_2", "Ke20","Kl14","Pm20","Pp15","Lp20","Mq20","Aq17","Mt19","Av19","Nw17","Py19")
)

pbmc@meta.data$dataset_proto <- factor(
  pbmc@meta.data$dataset_proto,
  levels = c("Kb20_1", "Kb20_2", "Ke20","Kl14","Pm20","Pp15","Lp20","Mq20","Aq17","Mt19","Av19","Nw17","Py19", setdiff(levels(pbmc@meta.data$diff_ek), c("Kb20_1", "Kb20_2", "Ke20","Kl14","Pm20","Pp15","Lp20","Mq20","Aq17","Mt19","Av19","Nw17","Py19")))
)
pbmc$dataset_proto <- droplevels(pbmc$dataset_proto)


pbmc@meta.data$diff_ek <- mapvalues(
  pbmc@meta.data$diff_ek,
  from = c(
    "Dorsal Telencephalic IP",
    "Dorsal Telencephalic NPC",
    "Dorsal Telencephalic Neuron"
  ),
  to   = c("IPC", "RG", "N")
)

pbmc@meta.data$diff_ek <- factor(
  pbmc@meta.data$diff_ek,
  levels = c("IPC", "RG", "N", setdiff(levels(pbmc@meta.data$diff_ek), c("IPC","RG","N")))
)


library(dplyr)
library(Seurat)

# 1. Extract metadata and compute per‐group sizes
meta <- pbmc@meta.data %>% 
  mutate(cell = rownames(.))
counts <- meta %>% count(dataset_proto)
min_n <- min(counts$n)
message("Each subsample will draw ", min_n, " cells per (dataset) group.")

# 2. Set up
n_perm <- 5
subsamples_cells <- vector("list", n_perm)
names(subsamples_cells) <- paste0("perm_", seq_len(n_perm))

# 3. For each diff_ek × age_ek “bin”…
group_list <- split(meta$cell, list(meta$dataset_proto), drop = TRUE)

table(pbmc$dataset_proto)

for (bin_cells in group_list) {
  bin_size   <- length(bin_cells)
  total_need <- n_perm * min_n
  
  # a) If we have at least n_perm*min_n cells, sample that many once…
  if (bin_size >= total_need) {
    pool <- sample(bin_cells, total_need, replace = FALSE)
    
    # b) Otherwise, take all cells + draw the remainder with replacement
  } else {
    extra_needed <- total_need - bin_size
    pool <- c(
      bin_cells,
      sample(bin_cells, extra_needed, replace = TRUE)
    )
    pool <- sample(pool)  # shuffle
  }
  
  # 4. Slice that “pool” into n_perm chunks of size = min_n
  for (i in seq_len(n_perm)) {
    idx_start <- (i - 1) * min_n + 1
    idx_end   <- i * min_n
    subsamples_cells[[i]] <- c(
      subsamples_cells[[i]],
      pool[idx_start:idx_end]
    )
  }
}

# 5. Build Seurat objects
subsamples <- lapply(subsamples_cells, function(cells) subset(pbmc, cells = cells))
names(subsamples) <- names(subsamples_cells)

# 6. Compute overlap matrix
overlap <- matrix(
  0, n_perm, n_perm,
  dimnames = list(names(subsamples), names(subsamples))
)
for (i in seq_len(n_perm)) {
  for (j in seq_len(n_perm)) {
    overlap[i, j] <- length(
      intersect(subsamples_cells[[i]], subsamples_cells[[j]])
    )
  }
}
print(overlap)

# 7. Save each subsampled object
output_dir <- "subsamples"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (name in names(subsamples)) {
  file_path <- file.path(output_dir, paste0(name, ".rds"))
  saveRDS(subsamples[[name]], file = file_path)
  message("Saved ", name, " to ", file_path)
}
