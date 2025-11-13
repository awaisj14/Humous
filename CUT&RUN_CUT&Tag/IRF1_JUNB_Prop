## ---- PACKAGES ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)     # import/export + liftOver
  library(stringr)
  # For motif analysis:
  # install.packages("BiocManager")
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Ptroglodytes.UCSC.panTro2)
  library(JASPAR2022)
  library(TFBSTools)
  library(Biostrings)
  library(motifmatchr)
})
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19",
                       "BSgenome.Ptroglodytes.UCSC.panTro2",
                       "JASPAR2022","TFBSTools","Biostrings","motifmatchr"))


HARs <- read.csv("/Users/javed/Documents/Humous_indexing/Analysis_reviewers/HARs/GSE140983_COUNTS_ANNOTATED.txt.gz", sep = "\t")
HARs



## ---- CONFIG ----
# Paths to TF peaks (hg38)
junb_bed <- "/USers/javed/Documents/Humous_indexing/Analysis_reviewers/CUTTAG/Org/New/SEACR_GREAT/JunbCAU_top0.01.peaks.stringent_cleaned.bed"
irf1_bed <- "/USers/javed/Documents/Humous_indexing/Analysis_reviewers/CUTTAG/CR/SEACR_GREAT/Irf1PCAU_top0.01.stringent_cleaned.bed"
chain_hg38_to_hg19 <- "hg38ToHg19.over.chain"

min_total_counts <- 50  # threshold for "high-count" HARs

#--------------------------------------------------------------------
# 1. PREP HARs BY SPECIES
#--------------------------------------------------------------------
make_gr <- function(chr, start, end) GRanges(seqnames=chr, ranges=IRanges(start=start, end=end))

parse_alignment_coords <- function(aln_vec) {
  x <- str_extract(aln_vec, "Chr[0-9XYM]+:[0-9]+-[0-9]+")
  tibble(raw = x) |>
    separate(raw, into=c("Chr","Range"), sep=":") |>
    separate(Range, into=c("Start","Stop"), sep="-", convert=TRUE) |>
    mutate(chr = str_replace(Chr, "^Chr", "chr")) |>
    select(chr, Start, Stop)
}

HARs <- HARs %>% mutate(total_counts = Inert_Counts + Comp_Counts)

# hg19 subset
HARs_hg19_df <- HARs %>% filter(Species == "hg19")
hg19_coords <- parse_alignment_coords(HARs_hg19_df$Alignment)
HARs_hg19 <- make_gr(hg19_coords$chr, hg19_coords$Start, hg19_coords$Stop)
mcols(HARs_hg19) <- HARs_hg19_df

# panTro2 subset
HARs_pt2_df <- HARs %>% filter(Species == "PanTro2")
chr_pt2 <- ifelse(str_detect(HARs_pt2_df$Chr,"^chr"), HARs_pt2_df$Chr, paste0("chr", HARs_pt2_df$Chr))
HARs_pt2 <- make_gr(chr_pt2, HARs_pt2_df$Start, HARs_pt2_df$Stop)
mcols(HARs_pt2) <- HARs_pt2_df

#--------------------------------------------------------------------
# 2. LIFTOVER JUNB / IRF1 peaks FROM hg38 TO hg19
#--------------------------------------------------------------------
JUNB_hg38 <- import(junb_bed)
IRF1_hg38 <- import(irf1_bed)

chain_38to19 <- import.chain(chain_hg38_to_hg19)

JUNB_hg19 <- unlist(liftOver(JUNB_hg38, chain_38to19))
IRF1_hg19 <- unlist(liftOver(IRF1_hg38, chain_38to19))

#--------------------------------------------------------------------
# 3. OVERLAPS WITH hg19 HARs
#--------------------------------------------------------------------
ov_JUNB <- findOverlaps(HARs_hg19, JUNB_hg19, ignore.strand=TRUE)
ov_IRF1 <- findOverlaps(HARs_hg19, IRF1_hg19, ignore.strand=TRUE)

har_idx_JUNB <- unique(queryHits(ov_JUNB))
har_idx_IRF1 <- unique(queryHits(ov_IRF1))

summary_all <- tibble(
  set = c("JUNB_overlaps","IRF1_overlaps","Either","Both"),
  n = c(length(har_idx_JUNB),
        length(har_idx_IRF1),
        length(union(har_idx_JUNB, har_idx_IRF1)),
        length(intersect(har_idx_JUNB, har_idx_IRF1)))
)
print(summary_all)

#--------------------------------------------------------------------
# 4. FILTER HIGH-COUNT HARs AND RECOUNT
#--------------------------------------------------------------------
HARs_hg19_high <- HARs_hg19[HARs_hg19$total_counts >= min_total_counts]
ov_JUNB_high <- findOverlaps(HARs_hg19_high, JUNB_hg19, ignore.strand=TRUE)
ov_IRF1_high <- findOverlaps(HARs_hg19_high, IRF1_hg19, ignore.strand=TRUE)

summary_high <- tibble(
  set = c("JUNB_highCount","IRF1_highCount","Either_high","Both_high"),
  n = c(length(unique(queryHits(ov_JUNB_high))),
        length(unique(queryHits(ov_IRF1_high))),
        length(union(unique(queryHits(ov_JUNB_high)),
                     unique(queryHits(ov_IRF1_high)))),
        length(intersect(unique(queryHits(ov_JUNB_high)),
                         unique(queryHits(ov_IRF1_high)))))
)
print(summary_high)

#--------------------------------------------------------------------
# 5. MOTIF ENRICHMENT FOR panTro2 HARs
#--------------------------------------------------------------------
# get JASPAR motifs for JUNB and IRF1
pfms <- getMatrixSet(JASPAR2022, opts=list(collection="CORE", name=c("JUNB","IRF1"), all_versions=FALSE))
pwms <- lapply(pfms, toPWM)

# get sequences for panTro2 HARs
seqs_pt2 <- getSeq(BSgenome.Ptroglodytes.UCSC.panTro2, HARs_pt2)

# motif matching
motif_matches <- matchMotifs(pwms, seqs_pt2, genome=BSgenome.Ptroglodytes.UCSC.panTro2, out="matches")

motif_summary <- tibble(
  motif = names(pwms),
  n_with = colSums(assay(motif_matches)),
  n_total = length(seqs_pt2),
  fraction = n_with / n_total
)
print(motif_summary)



library(GenomicRanges)

# Option: set a flank/buffer around HARs (bp). Use 0 for exact (≥1bp) overlaps.
har_flank <- 0
HARs_hg19_win <- if (har_flank > 0) flank(HARs_hg19, width = har_flank, both = TRUE) else HARs_hg19

frac_overlap <- function(peaks, regions) {
  ov <- findOverlaps(peaks, regions, ignore.strand = TRUE)
  n_peaks <- length(peaks)
  n_peaks_in_regions <- length(unique(queryHits(ov)))   # count each peak once
  data.frame(
    n_peaks = n_peaks,
    n_in_HARs = n_peaks_in_regions,
    fraction = n_peaks_in_regions / n_peaks
  )
}

prop_JUNB <- frac_overlap(JUNB_hg19, HARs_hg19_win)
prop_IRF1 <- frac_overlap(IRF1_hg19, HARs_hg19_win)

prop_JUNB
prop_IRF1

library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(patchwork)  # optional, to arrange the two pies side by side

# Find overlaps
ov_JUNB <- findOverlaps(JUNB_hg19, HARs_hg19, ignore.strand = TRUE)
ov_IRF1 <- findOverlaps(IRF1_hg19, HARs_hg19, ignore.strand = TRUE)

# Count unique peaks overlapping HARs
JUNB_in  <- length(unique(queryHits(ov_JUNB)))
JUNB_out <- length(JUNB_hg19) - JUNB_in

IRF1_in  <- length(unique(queryHits(ov_IRF1)))
IRF1_out <- length(IRF1_hg19) - IRF1_in

# Create a tidy dataframe for plotting
pie_data <- tibble(
  TF = c(rep("JUNB", 2), rep("IRF1", 2)),
  Category = rep(c("In_HARs", "Not_in_HARs"), 2),
  Count = c(JUNB_in, JUNB_out, IRF1_in, IRF1_out)
) %>%
  group_by(TF) %>%
  mutate(Fraction = Count / sum(Count),
         Label = paste0(Category, "\n", round(Fraction * 100, 2), "%"))

pie_data


ov_JUNB_merged <- findOverlaps(HARs_hg19_merged, JUNB_hg19, ignore.strand = TRUE)
ov_IRF1_merged <- findOverlaps(HARs_hg19_merged, IRF1_hg19, ignore.strand = TRUE)

HARs_JUNB_idx <- unique(queryHits(ov_JUNB_merged))
HARs_IRF1_idx <- unique(queryHits(ov_IRF1_merged))

length(intersect(HARs_JUNB_idx, HARs_IRF1_idx))  # shared HARs

# get HAR sequences that overlap JUNB vs IRF1
seqs_JUNB_HARs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, HARs_hg19_merged[HARs_JUNB_idx])
seqs_IRF1_HARs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, HARs_hg19_merged[HARs_IRF1_idx])

# then run matchMotifs or motifmatchr on both sets to see motif enrichment patterns


# JUNB and IRF1 overlap objects (from your merged HARs)
ov_JUNB_merged <- findOverlaps(HARs_hg19_merged, JUNB_hg19, ignore.strand = TRUE)
ov_IRF1_merged <- findOverlaps(HARs_hg19_merged, IRF1_hg19, ignore.strand = TRUE)

# HAR indices with overlaps
HARs_JUNB_idx <- unique(queryHits(ov_JUNB_merged))
HARs_IRF1_idx <- unique(queryHits(ov_IRF1_merged))

# Extract those HARs
HARs_JUNB_bound <- HARs_hg19_merged[HARs_JUNB_idx]
HARs_IRF1_bound <- HARs_hg19_merged[HARs_IRF1_idx]

# If you want HARs bound by both:
HARs_both_bound <- HARs_hg19_merged[intersect(HARs_JUNB_idx, HARs_IRF1_idx)]

# Convert to a readable data frame
HARs_JUNB_df <- as.data.frame(HARs_JUNB_bound)[, c("seqnames", "start", "end")]
HARs_IRF1_df <- as.data.frame(HARs_IRF1_bound)[, c("seqnames", "start", "end")]
HARs_both_df <- as.data.frame(HARs_both_bound)[, c("seqnames", "start", "end")]

# Print
HARs_JUNB_df
HARs_IRF1_df
HARs_both_df
