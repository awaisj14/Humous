#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# File paths
junb_file <- "bed_files/JunbK27_top0.01.peaks.fixed.bed"
ctrl_file <- "bed_files/ControlK27_top0.01.peaks.fixed.bed"
#!/usr/bin/env Rscript

# Import BEDs
junb <- import(junb_file, format = "BED")
ctrl <- import(ctrl_file, format = "BED")

# Keep canonical chromosomes
canonical_chrs <- paste0("chr", c(1:19, "X", "Y", "M"))
junb <- keepSeqlevels(junb, canonical_chrs, pruning.mode = "coarse")
ctrl <- keepSeqlevels(ctrl, canonical_chrs, pruning.mode = "coarse")

# Find overlaps with at least 50 bp
hits <- findOverlaps(junb, ctrl, minoverlap = 50)

# Indices of overlapping peaks
common_junb_idx <- unique(queryHits(hits))
common_ctrl_idx <- unique(subjectHits(hits))

# Subsets
common_junb <- junb[common_junb_idx]
junb_only   <- junb[-common_junb_idx]
ctrl_only   <- ctrl[-common_ctrl_idx]

# --- Export full BEDs (all columns preserved) ---
export(common_junb, "Junb_vs_Ctrl_common.bed", format = "BED")
export(junb_only,  "Junb_vs_Ctrl_JunbOnly.bed", format = "BED")
export(ctrl_only,  "Junb_vs_Ctrl_CtrlOnly.bed", format = "BED")

# --- Export GREAT-ready BED3 (chr start end only) ---
toBED3 <- function(gr) {
  data.frame(
    seqnames = as.character(seqnames(gr)),
    start = start(gr) - 1,   # 0-based start
    end = end(gr)
  )
}

common_bed3   <- toBED3(common_junb)
junb_only_bed3 <- toBED3(junb_only)
ctrl_only_bed3 <- toBED3(ctrl_only)

write.table(common_bed3, "Junb_vs_Ctrl_common_forGREAT_H3K27ac.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(junb_only_bed3, "Junb_vs_Ctrl_JunbOnly_forGREAT_H3K27ac.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ctrl_only_bed3, "Junb_vs_Ctrl_CtrlOnly_forGREAT_H3K27ac.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- Summary ---
cat("Peak comparison (≥50 bp overlap, canonical chromosomes only):\n")
cat("  Common peaks:   ", length(common_junb), "\n")
cat("  Junb-only peaks:", length(junb_only), "\n")
cat("  Ctrl-only peaks:", length(ctrl_only), "\n")
cat("\nFiles written:\n")
cat("  Full BEDs:   Junb_vs_Ctrl_common.bed, Junb_vs_Ctrl_JunbOnly.bed, Junb_vs_Ctrl_CtrlOnly.bed\n")
cat("  GREAT BED3:  Junb_vs_Ctrl_common_forGREAT.bed, Junb_vs_Ctrl_JunbOnly_forGREAT.bed, Junb_vs_Ctrl_CtrlOnly_forGREAT.bed\n")





# Files - These are coming from the GREAT analysis of the specific bed files above
control_file <- "Gene_lists/JUNB_OE_CONTROL_ONLY_H3K27ac_peaks.txt"
junb_file    <- "Gene_lists/JUNB_OE_JUNB_ONLY_H3K27ac_peaks.txt"

# Function to read GREAT output and keep only gene names
read_great_genes <- function(file) {
  lines <- read.delim(file, header = FALSE, comment.char = "#")
  genes <- lines[[1]]   # first column is gene name
  unique(genes)
}

# Read lists
genes_control <- read_great_genes(control_file)
genes_junb    <- read_great_genes(junb_file)

# Find overlap
common_genes <- intersect(genes_control, genes_junb)

# Save results
write.table(common_genes, "Common_genes_CONTROLvsJUNB.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Print summary
cat("Number of genes in Control-only list: ", length(genes_control), "\n")
cat("Number of genes in Junb-only list:    ", length(genes_junb), "\n")
cat("Number of overlapping genes:          ", length(common_genes), "\n")
cat("Overlapping genes written to Common_genes_CONTROLvsJUNB.txt\n")


# Function to read GREAT output and keep only gene names
read_great_genes <- function(file) {
  lines <- read.delim(file, header = FALSE, comment.char = "#")
  genes <- lines[[1]]   # first column is gene name
  unique(genes)
}

# Read lists
genes_control <- read_great_genes(control_file)
genes_junb    <- read_great_genes(junb_file)

# Overlap and exclusive sets
common_genes   <- intersect(genes_control, genes_junb)
control_only   <- setdiff(genes_control, genes_junb)
junb_only      <- setdiff(genes_junb, genes_control)

# Save results
write.table(common_genes, "Common_genes_CONTROLvsJUNB.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(control_only, "CONTROL_only_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(junb_only, "JUNB_only_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Print summary
cat("Genes in Control-only list: ", length(genes_control), "\n")
cat("Genes in Junb-only list:    ", length(genes_junb), "\n")
cat("Shared genes:               ", length(common_genes), "\n")
cat("Unique to Control:          ", length(control_only), "\n")
cat("Unique to Junb:             ", length(junb_only), "\n")
cat("\nFiles written:\n")
cat("  Common_genes_CONTROLvsJUNB.txt\n")
cat("  CONTROL_only_genes.txt\n")
cat("  JUNB_only_genes.txt\n")

# Function: GO enrichment (Biological Process) with clusterProfiler
run_go <- function(gene_list, label) {
  ids <- bitr(gene_list, fromType = "SYMBOL",
              toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  if (nrow(ids) == 0) {
    cat("⚠️ No mappable genes for", label, "\n")
    return(NULL)
  }
  
  ego <- enrichGO(gene = ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  if (is.null(ego) || nrow(ego) == 0) {
    cat("⚠️ No enriched GO terms for", label, "\n")
    return(NULL)
  }
  
  # Save table
  out_file <- paste0(label, "_GO_enrichment.tsv")
  write.table(as.data.frame(ego), out_file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save plots
  pdf(paste0(label, "_GO_barplot.pdf"), width = 8, height = 6)
  print(barplot(ego, showCategory = 20, title = paste(label, "GO BP")))
  dev.off()
  
  pdf(paste0(label, "_GO_dotplot.pdf"), width = 8, height = 6)
  print(dotplot(ego, showCategory = 20, title = paste(label, "GO BP")))
  dev.off()
  
  cat("✅ GO enrichment + plots saved for", label, "\n")
  return(ego)
}

# Run enrichment
ego_control <- run_go(control_only, "CONTROL_only")
ego_junb    <- run_go(junb_only, "JUNB_only")

# Print summary
cat("\nSummary:\n")
cat("Genes in Control-only list: ", length(genes_control), "\n")
cat("Genes in Junb-only list:    ", length(genes_junb), "\n")
cat("Shared genes:               ", length(common_genes), "\n")
cat("Unique to Control:          ", length(control_only), "\n")
cat("Unique to Junb:             ", length(junb_only), "\n")




library(readr)
library(dplyr)
library(eulerr)

# --- DEG sets ---
junbOE_df <- read_csv("DEGs_RG_markers_Junb_ek.csv")

upregulated_genes <- junbOE_df %>%
  filter(avg_log2FC > 0.25) %>%
  pull(gene) %>%
  toupper() %>%
  unique()

downregulated_genes <- junbOE_df %>%
  filter(avg_log2FC < -0.25) %>%
  pull(gene) %>%
  toupper() %>%
  unique()

# --- GRN JUNB targets ---
grn_df <- read_csv("RG_GRN_Final.csv")

junb_targets_pos <- grn_df %>%
  filter(source == "JUNB", coef_mean > 0) %>%
  pull(target) %>%
  toupper() %>%
  unique()

junb_targets_neg <- grn_df %>%
  filter(source == "JUNB", coef_mean < 0) %>%
  pull(target) %>%
  toupper() %>%
  unique()

# --- H3K27ac sets ---
junb_only <- scan("JUNB_only_genes.txt", what = character(), quiet = TRUE) %>%
  toupper() %>%
  unique()

control_only <- scan("CONTROL_only_genes.txt", what = character(), quiet = TRUE) %>%
  toupper() %>%
  unique()

# ========================
# Euler A: Upregulation
# ========================
venn_up <- euler(
  list(
    DEGs_up      = upregulated_genes,
    JUNB_targets_pos = junb_targets_pos,
    H3K27ac_gain = junb_only
  )
)

pdf("Euler_Upregulated.pdf", width = 7, height = 7)
plot(
  venn_up,
  fills = c("tomato", "goldenrod1", "skyblue"),
  alpha = 0.5,
  labels = list(cex = 1.1),
  quantities = list(cex = 1.3),
  main = "Upregulated: DEGs ∩ JUNB_targets(+) ∩ H3K27ac Gain"
)
dev.off()

# ========================
# Euler B: Downregulation
# ========================
venn_down <- euler(
  list(
    DEGs_down     = downregulated_genes,
    JUNB_targets_neg = junb_targets_neg,
    H3K27ac_loss  = control_only
  )
)

pdf("Euler_Downregulated.pdf", width = 7, height = 7)
plot(
  venn_down,
  fills = c("firebrick3", "darkolivegreen3", "grey60"),
  alpha = 0.5,
  labels = list(cex = 1.1),
  quantities = list(cex = 1.3),
  main = "Downregulated: DEGs ∩ JUNB_targets(-) ∩ H3K27ac Loss"
)
dev.off()










library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# --- Prepare input as a named list ---
gene_sets <- list(
  Activated   = junb_overlap_activated,
  Repressed = control_overlap_repressed
)

# --- Map symbols to Entrez IDs ---
gene_sets_entrez <- lapply(gene_sets, function(gs) {
  ids <- bitr(gs, fromType = "SYMBOL",
              toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ids$ENTREZID
})

# --- Run compareCluster with GO ---
cc_go <- compareCluster(
  geneCluster = gene_sets_entrez,
  fun         = "enrichGO",
  OrgDb       = org.Mm.eg.db,
  keyType     = "ENTREZID",
  ont         = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable    = TRUE
)

# --- Save results ---
write.table(as.data.frame(cc_go), "CompareCluster_GO_JUNBvsCONTROL.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Plots ---
pdf("CompareCluster_GO_dotplot.pdf", width = 10, height = 6)
dotplot(cc_go, showCategory = 10) + ggplot2::ggtitle("GO comparison: JUNB vs CONTROL")
dev.off()


library(dplyr)

# Extract results table
cc_res <- as.data.frame(cc_go)

# Order by GeneRatio (highest first) within each gene set
cc_res <- cc_res %>%
  group_by(Cluster) %>%
  arrange(desc(GeneRatio), .by_group = TRUE)

# Convert back to compareClusterResult
cc_go@compareClusterResult <- cc_res

# Now plot
pdf("CompareCluster_GO_dotplot_byGeneRatio.pdf", width = 10, height = 6)
dotplot(cc_go, showCategory = 10) +
  ggplot2::ggtitle("GO comparison: JUNB vs CONTROL (ordered by GeneRatio)")
dev.off()

pdf("CompareCluster_GO_barplot_byGeneRatio.pdf", width = 10, height = 6)





