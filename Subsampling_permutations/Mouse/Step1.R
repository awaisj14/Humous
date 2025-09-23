source("lib_misc.R")
library(Seurat) ; library(dplyr); library(readr) ; library(ggplot2) ; library(plyr) ; library(dplyr) ; library(scales) ;library(SummarizedExperiment)
; library(Matrix)
options(Seurat.object.assay.version = "v3")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
Ad21=readRDS("Ad21.rds") #smb://nasac-m2.isis.unige.ch/m-GJabaudon/GJabaudon/Esther-Lucia-Julien/humous_Lcopy/data/Ad21.rds
My17=readRDS("My17.rds")
LMl21=readRDS("ForebrainDorsal_SumExp_LaManno.RDS")
Jt19=readRDS("Jt19.rds")
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
counts <- Ad21@assays$RNA@counts
rownames(counts) <- toupper(rownames(counts))
meta <- Ad21@meta.data
Ad21 <- CreateSeuratObject(counts, meta.data = meta, project = "Ad21")

counts <- My17@assays$RNA@counts
rownames(counts) <- toupper(rownames(counts))
meta <- My17@meta.data
My17 <- CreateSeuratObject(counts, meta.data = meta, project = "My17")



cts   <- assay(LMl21) 
genes <- as.character(rowData(LMl21)$Gene)
counts_by_gene <- rowsum(cts, group = genes, reorder = FALSE)
counts_by_gene
rownames(counts_by_gene) <- toupper(rownames(counts_by_gene))
meta <- LMl21@colData
LMl21 <- CreateSeuratObject(counts_by_gene, meta.data = as.data.frame(meta), project = "LMl21")

counts <- Jt19@assays$RNA@counts
meta <- Jt19@meta.data
Jt19 <- CreateSeuratObject(counts, meta.data = meta, project = "Jt19")

rownames(Ad21)
rownames(LMl21)
rownames(My17)
rownames(Jt19)


# load annotations prepared by esther, filter and attach age and diff vectors

# Ad21
annot_Ad21 <- read_csv("DiffTypes_Ad21.csv") ; table(annot_Ad21$age,annot_Ad21$DiffTypes,useNA="always") ; dim(annot_Ad21)
Ad21 <- subset(Ad21,cells=annot_Ad21$cells[annot_Ad21$DiffTypes != "CR"]) # subset
Ad21$age_ek <- annot_Ad21$age[match(colnames(Ad21),annot_Ad21$cells)]
Ad21$diff_ek <- annot_Ad21$DiffTypes[match(colnames(Ad21),annot_Ad21$cells)]
table(Ad21$age_ek,Ad21$diff_ek,useNA = "always")
Ad21$dataset <- "Ad21"

# My17
annot_My17 <- read_csv("DiffTypes_My17.csv") ; table(annot_My17$age,annot_My17$DiffTypes,useNA="always") ; dim(annot_My17)
My17 <- subset(My17,cells=annot_My17$cells[ ! annot_My17$DiffTypes %in% c("CR","SP","IPCearly") & annot_My17$age!="11" ] ) # subset
My17$age_ek <- annot_My17$age[match(colnames(My17),annot_My17$cells)]
My17$diff_ek <- annot_My17$DiffTypes[match(colnames(My17),annot_My17$cells)]
table(My17$age_ek,My17$diff_ek,useNA = "always") 
My17$dataset <- "My17"

# LMl21
annot_LMl21 <- read_csv("DiffTypes_LMl21.csv")  ; table(annot_LMl21$age,annot_LMl21$DiffTypes,useNA="always") ; dim(annot_LMl21)
LMl21 <- subset(LMl21,cells=annot_LMl21$cells[ ! annot_LMl21$DiffTypes %in% c("CR")]) # subset
LMl21$age_ek <- annot_LMl21$age[match(colnames(LMl21),annot_LMl21$cells)]
LMl21$diff_ek <- annot_LMl21$DiffTypes[match(colnames(LMl21),annot_LMl21$cells)]
table(LMl21$age_ek,LMl21$diff_ek,useNA = "always") # all good
LMl21$dataset <- "LMl21"

# Jt19
annot_Jt19 <- read_csv("DiffTypes_Jt19.csv")  ; table(annot_Jt19$age,annot_Jt19$DiffTypes,useNA="always") ; dim(annot_Jt19)
# all good, no need to subset
Jt19$age_ek <- annot_Jt19$age[match(colnames(Jt19),annot_Jt19$cells)] 
table(Jt19$age_ek,useNA = "always") # there are some NA cells, remove them
Jt19$age_ek <- ifelse(is.na(Jt19$age_ek),"toremove",as.character(Jt19$age_ek))
Jt19 <- subset(Jt19,age_ek!="toremove")
Jt19$diff_ek <- annot_Jt19$DiffTypes[match(colnames(Jt19),annot_Jt19$cells)]
table(Jt19$age_ek,Jt19$diff_ek,useNA = "always") 
Jt19$dataset <- "Jt19"

# Step 1: Get shared features across all four datasets
shared_features_all <- Reduce(intersect, list(
  rownames(Ad21),
  rownames(Jt19),
  rownames(My17),
  rownames(LMl21)
))


# Step 2: Subset each object to the shared features
Ad21_sub <- subset(Ad21, features = shared_features_all)
Jt19_sub <- subset(Jt19, features = shared_features_all)
My17_sub <- subset(My17, features = shared_features_all)
LMl21_sub <- subset(LMl21, features = shared_features_all)

# Step 3: Merge the subsetted objects
merged_all <- merge(Ad21_sub, y = list(Jt19_sub, My17_sub, LMl21_sub), 
                    add.cell.ids = c("Ad21", "Jt19", "My17", "LMl21"), 
                    project = "AllMerged")

VlnPlot(merged_all, features = "NEUROD2", group.by = "dataset")


# Step 1: Get shared features between Ad21 and LMl21
shared_features_ad_lml <- intersect(rownames(Ad21), rownames(LMl21))

# Step 2: Subset to shared features
Ad21_sub2 <- subset(Ad21, features = shared_features_ad_lml)
LMl21_sub2 <- subset(LMl21, features = shared_features_ad_lml)

# Step 3: Merge the two
merged_ad_lml <- merge(Ad21_sub2, y = LMl21_sub2, 
                       add.cell.ids = c("Ad21", "LMl21"), 
                       project = "Ad21_LMl21_Merged")

table(merged_all$age_ek ,useNA = "always") # all good
table(merged_ad_lml$age_ek ,useNA = "always") # all good

merged_ad_lml$age_ek <- as.character(merged_ad_lml$age_ek)
merged_all$age_ek <- as.character(merged_all$age_ek)

# re-define age groups
merged_all$age_ek <- plyr::revalue(merged_all$age_ek,c("12"="12","13"="13","14"="14","15"="15","16"="16-17","17"="16-17"))
merged_ad_lml$age_ek <- plyr::revalue(merged_ad_lml$age_ek,c("12"="12","13"="13","14"="14","15"="15","16"="16-17","17"="16-17"))

table(merged_all$age_ek)
table(merged_ad_lml$age_ek)

saveRDS(merged_all, "merged_all.rds")
saveRDS(merged_ad_lml, "merged_ad_lml.rds")


# Load necessary libraries
library(data.table)       # for fast CSV reading
library(DelayedArray)     # assuming L_ref is a DelayedArray

# Step 1: Read landinfo.csv
landinfo <- fread("landinfo.csv")

# Step 2: Filter for type == "P"
genes_PM <- landinfo[type == "P" & cond == "M", gene]  # Replace 'gene' if your gene column is named differently
genes_PM <- toupper(genes_PM)
# Step 3: Subset L_ref by gene names (rows)
# Make sure gene names are rownames of L_ref
L_ref_subset <- L_ref[rownames(L_ref) %in% genes_P,,]

saveRDS(L_ref_subset, "landsM_OG/L_MR_M_patterned.rds")













#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# -----------------------------
# CONFIGURATION
# -----------------------------
merged_rds <- "merged_all.rds"
base_out   <- "permutations_all"
n_perms    <- 5
cells_per_ds <- 4000    # max cells per dataset
base_seed  <- 101

# -----------------------------
# FUNCTION: global random sampling with dataset cap
# -----------------------------
global_cap_selector <- function(merged, n_cap = 4000) {
  meta <- merged@meta.data
  all_cells <- rownames(meta)
  all_ds <- as.character(meta$dataset)  # ensure character
  
  # Shuffle globally
  shuffled <- sample(all_cells, length(all_cells))
  
  # Initialize per-dataset counters
  unique_ds <- unique(all_ds)
  picked_per_ds <- setNames(rep(0L, length(unique_ds)), unique_ds)
  
  picked <- character(0)
  
  for (cell in shuffled) {
    # Use numeric index, not name-based indexing
    idx <- match(cell, all_cells)
    ds <- all_ds[idx]
    
    # Safety check: skip cells without dataset assignment
    if (is.na(ds)) next
    
    # Skip if dataset cap reached
    if (picked_per_ds[ds] >= n_cap) next
    
    picked <- c(picked, cell)
    picked_per_ds[ds] <- picked_per_ds[ds] + 1
  }
  
  picked
}

# -----------------------------
# MAIN LOGIC
# -----------------------------
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
merged_all <- readRDS(merged_rds)

all_selected <- list()
for (i in seq_len(n_perms)) {
  i_padded <- sprintf("%02d", i)
  out_file <- file.path(base_out, paste0("draw_cells_perm_", i_padded, ".rds"))
  
  # Skip existing permutations
  if (file.exists(out_file)) {
    message("✔ Perm ", i, " already exists — skipping")
    all_selected[[i]] <- readRDS(out_file)
    next
  }
  
  # Attempt unique selection
  seed_attempt <- 0
  repeat {
    set.seed(base_seed + i + seed_attempt * 997 + sample.int(1e5, 1))
    selected <- global_cap_selector(merged_all, n_cap = cells_per_ds)
    
    identical_found <- FALSE
    for (prev_idx in seq_along(all_selected)) {
      if (setequal(all_selected[[prev_idx]], selected)) {
        identical_found <- TRUE
        seed_attempt <- seed_attempt + 1
        message("⚠️ Identical selection detected with perm ", prev_idx,
                " — retrying (attempt ", seed_attempt, ")")
        break
      }
    }
    if (!identical_found) break
  }
  
  all_selected[[i]] <- selected
  saveRDS(selected, out_file)
  message("✅ Saved unique selection for perm ", i)
}

# -----------------------------
# JACCARD OVERLAP REPORT
# -----------------------------
names(all_selected) <- sprintf("perm_%02d", seq_len(n_perms))

pm <- length(all_selected)
jaccard_matrix <- matrix(NA_real_, pm, pm,
                         dimnames = list(names(all_selected), names(all_selected)))

for (a in seq_len(pm)) {
  for (b in seq_len(pm)) {
    inter <- length(intersect(all_selected[[a]], all_selected[[b]]))
    uni   <- length(union(all_selected[[a]], all_selected[[b]]))
    jaccard_matrix[a,b] <- if (uni == 0) NA_real else inter / uni
  }
}

write.csv(jaccard_matrix,
          file = file.path(base_out, "jaccard_between_perms.csv"),
          row.names = TRUE)

message("📊 Pairwise Jaccard overlap saved to: ",
        file.path(base_out, "jaccard_between_perms.csv"))








#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# -----------------------------
# CONFIGURATION
# -----------------------------
merged_rds <- "merged_ad_lml.rds"
base_out   <- "permutations_adlml"
n_perms    <- 5
cells_per_ds <- 6000    # max cells per dataset
base_seed  <- 101

# -----------------------------
# FUNCTION: global random sampling with dataset cap
# -----------------------------
global_cap_selector <- function(merged, n_cap = 4000) {
  meta <- merged@meta.data
  all_cells <- rownames(meta)
  all_ds <- as.character(meta$dataset)  # ensure character
  
  # Shuffle globally
  shuffled <- sample(all_cells, length(all_cells))
  
  # Initialize per-dataset counters
  unique_ds <- unique(all_ds)
  picked_per_ds <- setNames(rep(0L, length(unique_ds)), unique_ds)
  
  picked <- character(0)
  
  for (cell in shuffled) {
    # Use numeric index, not name-based indexing
    idx <- match(cell, all_cells)
    ds <- all_ds[idx]
    
    # Safety check: skip cells without dataset assignment
    if (is.na(ds)) next
    
    # Skip if dataset cap reached
    if (picked_per_ds[ds] >= n_cap) next
    
    picked <- c(picked, cell)
    picked_per_ds[ds] <- picked_per_ds[ds] + 1
  }
  
  picked
}

# -----------------------------
# MAIN LOGIC
# -----------------------------
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
merged_all <- readRDS(merged_rds)

all_selected <- list()
for (i in seq_len(n_perms)) {
  i_padded <- sprintf("%02d", i)
  out_file <- file.path(base_out, paste0("draw_cells_perm_", i_padded, ".rds"))
  
  # Skip existing permutations
  if (file.exists(out_file)) {
    message("✔ Perm ", i, " already exists — skipping")
    all_selected[[i]] <- readRDS(out_file)
    next
  }
  
  # Attempt unique selection
  seed_attempt <- 0
  repeat {
    set.seed(base_seed + i + seed_attempt * 997 + sample.int(1e5, 1))
    selected <- global_cap_selector(merged_all, n_cap = cells_per_ds)
    
    identical_found <- FALSE
    for (prev_idx in seq_along(all_selected)) {
      if (setequal(all_selected[[prev_idx]], selected)) {
        identical_found <- TRUE
        seed_attempt <- seed_attempt + 1
        message("⚠️ Identical selection detected with perm ", prev_idx,
                " — retrying (attempt ", seed_attempt, ")")
        break
      }
    }
    if (!identical_found) break
  }
  
  all_selected[[i]] <- selected
  saveRDS(selected, out_file)
  message("✅ Saved unique selection for perm ", i)
}

# -----------------------------
# JACCARD OVERLAP REPORT
# -----------------------------
names(all_selected) <- sprintf("perm_%02d", seq_len(n_perms))

pm <- length(all_selected)
jaccard_matrix <- matrix(NA_real_, pm, pm,
                         dimnames = list(names(all_selected), names(all_selected)))

for (a in seq_len(pm)) {
  for (b in seq_len(pm)) {
    inter <- length(intersect(all_selected[[a]], all_selected[[b]]))
    uni   <- length(union(all_selected[[a]], all_selected[[b]]))
    jaccard_matrix[a,b] <- if (uni == 0) NA_real else inter / uni
  }
}

write.csv(jaccard_matrix,
          file = file.path(base_out, "jaccard_between_perms.csv"),
          row.names = TRUE)

message("📊 Pairwise Jaccard overlap saved to: ",
        file.path(base_out, "jaccard_between_perms.csv"))
