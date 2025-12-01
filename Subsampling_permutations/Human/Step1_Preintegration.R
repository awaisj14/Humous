# 0) load your helper library
source("lib_misc.R")   # defines dataset_integration, etc.
library(Seurat); library(dplyr); library(readr); library(plyr)

# 1) load raw and prefilter by Esther’s list i.e. RG IP N celltypes from the datasets
full.list <- list(
  Kb21 = readRDS("Hm/Kb21.rds"),
  Kn17 = readRDS("Hm/Kn17.rds"),
  Pt21 = readRDS("Hm/Pt21.rds")
)
Sel_all_Hm <- read_csv("SelectCells_all_Hm.csv")
full.list_f <- lapply(names(full.list), function(ds) {
  obj <- full.list[[ds]]
  keep <- intersect(Sel_all_Hm$x, colnames(obj))
  so   <- subset(obj, cells = keep)
  so$dataset <- ds
  so
})
names(full.list_f) <- names(full.list)

# 2) Annotate with age_ek/diff_ek
annot_paths <- list(
  Kb21 = "Hm/DiffTypes_Kb21.csv",
  Kn17 = "Hm/DiffTypes_Kn17.csv",
  Pt21 = "Hm/DiffTypes_Pt21.csv"
)
for(ds in names(full.list_f)) {
  annot <- read_csv(annot_paths[[ds]])
  cols  <- colnames(full.list_f[[ds]])
  idx   <- match(cols, annot$cells)
  full.list_f[[ds]]$age_ek  <- annot$age[idx]
  types  <- annot$DiffTypes[idx]
  full.list_f[[ds]]$diff_ek <- ifelse(types=="Neuron","N",types)
}

# 3) Merge and re‐bin ages
merged <- merge(full.list_f$Kb21, y=list(full.list_f$Kn17, full.list_f$Pt21))
merged$age_ek <- plyr::revalue(
  as.character(merged$age_ek),
  c("12"="12","15"="15-16","16"="15-16",
    "17"="17-18","18"="17-18",
    "20"="20-21","21"="20-21",
    "23"="23-24","24"="23-24")
)

# 4) Save the *base* merged object immediately before any subsampling
saveRDS(merged, file="merged_prefiltered.rds")


# 5) Now generate 10 random‐draw lists and corresponding full.list_f_i objects
set.seed(2025)
n_perms      <- 10
sample_ages  <- c("15-16","17-18","20-21","23-24")

for(i in seq_len(n_perms)) {
  # pick cells
  sel <- unlist(lapply(sample_ages, function(age) {
    cells_age <- colnames(merged)[merged$age_ek==age]
    sample(cells_age, min(6000, length(cells_age)))
  }))
  # always include week12
  sel <- c(sel, colnames(merged)[merged$age_ek=="12"])
  saveRDS(sel,      file=sprintf("draw_cells_perm_%02d.rds", i))
  
  # build and save full.list_f_i
  perm_list <- lapply(names(full.list_f), function(ds) {
    subset(full.list_f[[ds]], cells=intersect(sel, colnames(full.list_f[[ds]])))
  })
  names(perm_list) <- names(full.list_f)
  saveRDS(perm_list, file=sprintf("full.list_f_%02d.rds", i))
}


