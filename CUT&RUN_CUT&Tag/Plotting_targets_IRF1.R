#install.packages(c("data.table","dplyr","ggplot2","ggVennDiagram"), dependencies = TRUE)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)

# --- Helper to read GREAT gene lists (first col = gene symbol) ---
read_great_genes <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]          # drop GREAT header lines
  genes <- sub("\\t.*$", "", lines)                # keep text before first TAB
  genes <- toupper(trimws(genes))
  genes[nzchar(genes)] |> unique()
}

# --- Load GREAT gene list ---
great_cau <- read_great_genes("IrfP_CAU.txt")

# --- Load GRN and extract IRF1 targets ---
grn <- fread("RG_GRN_final.csv", data.table = FALSE)
ln <- tolower(names(grn))
src_col <- which(ln == "source")[1]; tgt_col <- which(ln == "target")[1]
stopifnot(!is.na(src_col), !is.na(tgt_col))
grn <- grn[, c(src_col, tgt_col)]; names(grn) <- c("source","target")

IRF1_targets <- grn |>
  filter(toupper(source) == "IRF1") |>
  pull(target) |>
  toupper() |>
  unique()

# --- Sets for Venn ---
sets <- list(
  CAU_binding = great_cau,
  IRF1_GRN    = IRF1_targets
)

# --- Plot Venn ---
p <- ggVennDiagram(sets, label_alpha = 0) +
  labs(title = "IRF1 GRN targets vs IRF1 binding (GREAT): CAU") +
  theme(legend.position = "none")
print(p)

# --- Useful overlap counts ---
IRF1_cau <- length(intersect(IRF1_targets, great_cau))
only_cau <- length(setdiff(intersect(great_cau, IRF1_targets), character(0)))

summary_df <- tibble::tibble(
  metric = c("IRF1∩CAU", "IRF1∩CAU only",
             "Total IRF1 GRN", "Total CAU binding"),
  n      = c(IRF1_cau, only_cau,
             length(IRF1_targets), length(great_cau))
)
print(summary_df)
