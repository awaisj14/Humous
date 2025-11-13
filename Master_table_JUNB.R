library(dplyr)
library(readr)

# -------------------------
# 1. Read input files
# -------------------------
rg_grn   <- read_csv("RG_GRN_Final.csv")
landinfo <- read_csv("landinfo.csv")

# Remove commented lines ("# ...") when reading
cau_top  <- read_delim("JUNB_CAU_TOP.txt", 
                       delim = "\t", 
                       col_names = FALSE, 
                       comment = "#")

junb_only <- read_delim("JUNB_only_genes.txt", 
                        delim = "\t", 
                        col_names = FALSE, 
                        comment = "#")

ctrl_only <- read_delim("Control_only_genes.txt", 
                        delim = "\t", 
                        col_names = FALSE, 
                        comment = "#")

degs     <- read_csv("DEGs_RG_markers_Junb_ek.csv")

# -------------------------
# 2. Prepare each dataset
# -------------------------

# GRN: keep JUNB sources only
rg_grn_sub <- rg_grn %>%
  dplyr::filter(source == "JUNB") %>%
  dplyr::mutate(GRN_regulation = ifelse(coef_mean > 0, "positive", "negative")) %>%
  dplyr::select(gene = target, GRN_regulation)

# Landinfo: split H vs M
land_H <- landinfo %>%
  dplyr::filter(cond == "H") %>%
  dplyr::select(gene, H_Lands = AnnotTypes)

land_M <- landinfo %>%
  dplyr::filter(cond == "M") %>%
  dplyr::select(gene, M_Lands = AnnotTypes)
land_M$gene <- toupper(land_M$gene)

# CAU_TOP: mark binding
cau_top <- cau_top %>%
  transmute(gene = X1, binding = TRUE)

# JUNB_only vs Control_only
junb_only <- junb_only %>%
  transmute(gene = X1, H3K27ac = "GAIN")
junb_only$gene <- toupper(junb_only$gene)

ctrl_only <- ctrl_only %>%
  transmute(gene = X1, H3K27ac = "LOSS")
ctrl_only$gene <- toupper(ctrl_only$gene)

h3k27ac <- bind_rows(junb_only, ctrl_only)
h3k27ac$gene <- toupper(h3k27ac$gene)

# DEGs: keep gene, log2FC, p_val
degs_sub <- degs %>%
  dplyr::select(gene, 
         Log2FC_snRNAseq = avg_log2FC, 
         p_val_snRNAseq = p_val)
degs_sub$gene <- toupper(degs_sub$gene)

# -------------------------
# 3. Master table assembly
# -------------------------
master <- rg_grn_sub %>%
  full_join(land_H, by = "gene") %>%
  full_join(land_M, by = "gene") %>%
  full_join(cau_top, by = "gene") %>%
  full_join(h3k27ac, by = "gene") %>%
  full_join(degs_sub, by = "gene") %>%
  dplyr::mutate(binding = ifelse(is.na(binding), FALSE, binding))

# -------------------------
# 4. Save output
# -------------------------
write_csv(master, "Master_JUNB_Table.csv")





##### Sources ########
library(dplyr)
library(readr)

# -------------------------
# 1. Read input files
# -------------------------
rg_grn   <- read_csv("RG_GRN_Final.csv")
landinfo <- read_csv("landinfo.csv")

# -------------------------
# 2. Filter for JUNB target
# -------------------------
junb_sources <- rg_grn %>%
  dplyr::filter(target == "JUNB")

# -------------------------
# 3. Prepare landinfo annotations for H, M, and O
# -------------------------
land_H <- landinfo %>%
  dplyr::filter(cond == "H") %>%
  dplyr::select(gene, Annot_H = AnnotTypes)

land_M <- landinfo %>%
  dplyr::filter(cond == "M") %>%
  dplyr::select(gene, Annot_M = AnnotTypes)

land_O <- landinfo %>%
  dplyr::filter(cond == "O") %>%
  dplyr::select(gene, Annot_O = AnnotTypes)

# Ensure consistent casing
land_H$gene <- toupper(land_H$gene)
land_M$gene <- toupper(land_M$gene)
land_O$gene <- toupper(land_O$gene)
junb_sources$source <- toupper(junb_sources$source)

# -------------------------
# 4. Join annotation info to JUNB source genes
# -------------------------
junb_sources_annot <- junb_sources %>%
  left_join(land_H, by = c("source" = "gene")) %>%
  left_join(land_M, by = c("source" = "gene")) %>%
  left_join(land_O, by = c("source" = "gene"))

# -------------------------
# 5. Save output
# -------------------------
write_csv(junb_sources_annot, "JUNB_Source_Genes_with_Annot.csv")

# Load the data
df <- read.csv("Master_JUNB_Table.csv", sep = ",", header = TRUE)

# 1. How many are GRN_regulation positive
pos_count <- sum(df$GRN_regulation == "positive", na.rm = TRUE)

# 2. Of those, how many have H_Lands == "RG" or "RG_IP"
pos_RG <- subset(df, GRN_regulation == "positive" & H_Lands %in% c("RG"))
pos_RG_count <- nrow(pos_RG)

# 3. Of those, how many have binding not TRUE
pos_RG_notbinding <- subset(pos_RG, binding != TRUE)
pos_RG_notbinding_count <- nrow(pos_RG_notbinding)

# 4. Among GRN positive genes, how many have H3K27ac == "GAIN" or "LOSS"
pos_gain <- sum(pos_RG$GRN_regulation == "positive" & pos_RG$H3K27ac == "GAIN", na.rm = TRUE)
pos_loss <- sum(df$GRN_regulation == "positive" & df$H3K27ac == "LOSS", na.rm = TRUE)

# Print results
cat("1. Total GRN positive genes:", pos_count, "\n")
cat("2. GRN positive with H_Lands RG or RG_IP:", pos_RG_count, "\n")
cat("3. Of those, binding not TRUE:", pos_RG_notbinding_count, "\n")
cat("4. GRN positive with H3K27ac GAIN:", pos_gain, "\n")
cat("   GRN positive with H3K27ac LOSS:", pos_loss, "\n")
