library(tidyverse)

# --- load perm files ---
# match both styles: "perm_01_vs_ref_per_gene.csv" and "perm_vs_ref_per_gene_perm05.csv"
files <- list.files(pattern = "perm.*vs_ref_per_gene.*\\.csv")

df_all <- files %>%
  map_df(~ read.csv(.x) %>% mutate(file = .x))

# average corr per gene
df_avg <- df_all %>%
  group_by(gene) %>%
  summarise(mean_corr = mean(corr, na.rm = TRUE))

# --- load landinfo ---
landinfo <- read.csv("landinfo.csv")

# filter for cond == H and type == P
land_subset <- landinfo %>%
  filter(cond == "M", type == "P")
land_subset$gene <- toupper(land_subset$gene)
# keep only those genes
df_filtered <- df_avg %>%
  semi_join(land_subset, by = "gene")

# categorize
df_filtered <- df_filtered %>%
  mutate(category = case_when(
    mean_corr > 0.5 ~ "High",
    mean_corr >= 0.25 & mean_corr <= 0.5 ~ "Mid",
    mean_corr >= 0 & mean_corr < 0.25 ~ "Low",
    mean_corr < 0 ~ "Anti"
  ))

# count + percentage
counts <- df_filtered %>%
  count(category) %>%
  mutate(percentage = round(100 * n / sum(n), 1),
         label = paste0(category, "\n", n, " (", percentage, "%)"))

# --- Pie chart with counts + percentages ---
ggplot(counts, aes(x = "", y = n, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(
    "High" = "forestgreen",
    "Mid"  = "steelblue",
    "Low"  = "orange",
    "Anti" = "firebrick"
  )) +
  labs(title = "H + P subset gene categories",
       fill = "Category") +
  theme_void()

# --- Optional: Bar chart ---
ggplot(counts, aes(x = category, y = n, fill = category)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", percentage, "%)")),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    "High" = "forestgreen",
    "Mid"  = "steelblue",
    "Low"  = "orange",
    "Anti" = "firebrick"
  )) +
  labs(title = "Counts of H + P subset gene categories",
       x = "Category", y = "Number of Genes") +
  theme_minimal()
