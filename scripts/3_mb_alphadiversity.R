# BARIA: Associations between muscle mass and gut microbiota
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(MetBrewer)
library(grid)
library(ggthemes)
library(ggpubr)
library(patchwork)

# Theme
renoir_cols_20 <- met.brewer("Renoir", n = 20)
fill_cols_2 <- scale_fill_manual(values = c("yes" = renoir_cols_20[18], "no" = renoir_cols_20[5]))

theme_minimal_custom <- function(base_size = 14, base_family = "sans") {

  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(0.8), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(0.8)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(colour = "black"),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      panel.grid.major = element_line(colour = "#dddddd", linewidth = 0.4, linetype = "22"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "bottom",
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )

}

# Data
baria_muscle_ab <- readRDS("data/20260722_BARIA_muscle_clinical.RDS") # metadata/clinical data
baria_mb <- readRDS("data/ps.BARIA.metaphlan.706.2548.RDS")

# Filter out participants taking antibiotics
baria_muscle <- baria_muscle_ab |> 
  filter(abx_v0 == "no")

# qc
sample_sums(baria_mb) |>
  summary() # adds up to 100

# Subset; keep only samples with one run or first run of samples with duplicates
run1_mb <- prune_samples(
  sample_data(baria_mb)$Extra_data == "NA" | sample_data(baria_mb)$Extra_data == "rep1",
  baria_mb
)

# Diversity metrics
# Convert OTU table to matrix
matrix_mb <- as(otu_table(run1_mb), "matrix")

# vegan requires a matrix with samples as rows and taxa as cols
if (taxa_are_rows(run1_mb)) { 
  matrix_mb <- t(matrix_mb) 
}

# Shannon, Simpson, Richness 
# note to self: calculated per sample, no filtering per visit needed at this stage
shannon <- vegan::diversity(matrix_mb, index = "shannon")
simpson <- vegan::diversity(matrix_mb, index = "simpson")
richness <- vegan::specnumber(matrix_mb)

# Create a df with Shannon, Simpson, Richness
alpha <- tibble(
  Sample = names(shannon),
  shannon = shannon,
  simpson = simpson,
  richness = richness
)

# Merge with metadata
baria_mb_df <- as(sample_data(run1_mb), "data.frame")

alpha_meta <- baria_mb_df |> 
  tibble::rownames_to_column(var = "Sample") |> 
  mutate(
    visit = str_extract(Time_Point, "\\d"),
    visit = if_else(visit == "1", "0", visit),
    visit = as.factor(visit),
    id = Subject_ID
  ) |> 
  select(Sample, id, visit) |> 
  left_join(alpha, by = "Sample") |>  
  inner_join(baria_muscle, by = "id") |> 
  relocate(id, .before = Sample)

#### Baseline Plots ####
alpha_v0 <- alpha_meta |>
  filter(visit == "0") |>
  arrange(id, Sample) |>
  distinct(id, .keep_all = TRUE)

## Shannon ##
# Boxplot
shannon_ffmi_v0 <- alpha_v0 |> 
  mutate(low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes", after = 0L)) |> 
  ggplot(aes(x = low_ffmi_v0, y = shannon, fill = low_ffmi_v0)) +
  geom_boxplot() +
  #geom_jitter(position = position_dodge(0.75)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "Shannon index", y = "Shannon index", x = "Low baseline FFMI") +
  fill_cols_2 +
  labs(fill = "Low baseline FFMI") + 
  theme_minimal_custom()
ggsave(shannon_ffmi_v0, filename = "graphs/alphadiversity/shannon_ffmi_v0.pdf", width = 7, height = 5)

# Violin
shannon_violin_ffmi_v0 <- alpha_v0 |> 
  mutate(low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes", after = 0L)) |>  
  ggplot(aes(x = low_ffmi_v0, y = shannon)) +
  geom_violin(aes(fill = low_ffmi_v0), trim = FALSE) +
  geom_boxplot(fill = "white", width = 0.1) +
  labs(x = "", y = "Shannon index", title = "Shannon index", fill = "Baseline FFMI status") +
  stat_compare_means( 
    tip.length = 0, 
    hide.ns = TRUE, 
    label.x = 1.5,
    method = "wilcox.test",
    label = "p.signif"
  ) + 
  fill_cols_2 +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  scale_x_discrete(labels = c("yes" = "Low FFMI", "no" = "Moderate/high FFMI")) +
  theme_minimal_custom() +
  theme(legend.position = "none")
ggsave(shannon_violin_ffmi_v0, filename = "graphs/alphadiversity/shannon_violin_ffmi_v0.pdf", width = 6, height = 5)

## Simpson ##
# Boxplot
simpson_ffmi_v0 <- alpha_v0 |> 
  mutate(low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes", after = 0L)) |> 
  ggplot(aes(x = low_ffmi_v0, y = simpson, fill = low_ffmi_v0)) +
  geom_boxplot() +
  #geom_jitter(position = position_dodge(0.75)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "Simpson index", y = "Simpson index", x = "Low baseline FFMI") +
  fill_cols_2 +
  labs(fill = "Low baseline FFMI") + 
  theme_minimal_custom()
ggsave(simpson_ffmi_v0, filename = "graphs/alphadiversity/simpson_ffmi_v0.pdf", width = 7, height = 5)

# Violin
simpson_violin_ffmi_v0 <- alpha_v0 |> 
  mutate(low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes", after = 0L)) |>  
  ggplot(aes(x = low_ffmi_v0, y = simpson)) +
  geom_violin(aes(fill = low_ffmi_v0), trim = FALSE) +
  geom_boxplot(fill = "white", width = 0.1) +
  labs(x = "", y = "Simpson index", title = "Simpson index", fill = "Baseline FFMI status") +
  stat_compare_means( 
    tip.length = 0, 
    hide.ns = TRUE, 
    label.x = 1.5,
    method = "wilcox.test",
    label = "p.signif"
  ) + 
  fill_cols_2 +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  scale_x_discrete(labels = c("yes" = "Low FFMI", "no" = "Moderate/high FFMI")) +
  theme_minimal_custom() +
  theme(legend.position = "none")
ggsave(simpson_violin_ffmi_v0, filename = "graphs/alphadiversity/simpson_violin_ffmi_v0.pdf", width = 6, height = 5)

## Richness ##
# Boxplot
richness_ffmi_v0 <- alpha_v0 |> 
  mutate(low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes", after = 0L)) |> 
  ggplot(aes(x = low_ffmi_v0, y = richness, fill = low_ffmi_v0)) +
  geom_boxplot() +
  #geom_jitter(position = position_dodge(0.75)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "Richness", y = "Richness", x = "Low baseline FFMI") +
  fill_cols_2 +
  labs(fill = "Low baseline FFMI") + 
  theme_minimal_custom()
ggsave(richness_ffmi_v0, filename = "graphs/alphadiversity/richness_ffmi_v0.pdf", width = 7, height = 5)

# Violin
richness_violin_ffmi_v0 <- alpha_v0 |> 
  mutate(low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes", after = 0L)) |>  
  ggplot(aes(x = low_ffmi_v0, y = richness)) +
  geom_violin(aes(fill = low_ffmi_v0), trim = FALSE) +
  geom_boxplot(fill = "white", width = 0.1) +
  labs(x = "", y = "Richness", title = "Richness", fill = "Baseline FFMI status") +
  stat_compare_means( 
    tip.length = 0, 
    hide.ns = TRUE, 
    label.x = 1.5,
    method = "wilcox.test",
    label = "p.signif"
  ) + 
  fill_cols_2 +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  scale_x_discrete(labels = c("yes" = "Low FFMI", "no" = "Moderate/high FFMI")) +
  theme_minimal_custom() +
  theme(legend.position = "none")
ggsave(richness_violin_ffmi_v0, filename = "graphs/alphadiversity/richness_violin_ffmi_v0.pdf", width = 6, height = 5)

## Combine into a panel 
# Horizontal
# Baseline (Shannon, Simpson, Richness)
# Violins
alpha_panel_ffmi_v0 <- 
  (shannon_violin_ffmi_v0 + simpson_violin_ffmi_v0 + richness_violin_ffmi_v0) +
  plot_layout(guides = "collect") &
  # theme(legend.position = "bottom") &
  theme(aspect.ratio = 0.8)
ggsave(alpha_panel_ffmi_v0, filename = "graphs/alphadiversity/alpha_panel_ffmi_v0.pdf", width = 12, height = 8)

# Vertical (for pptx)
alpha_panel_ffmi_v0_vertical <-
  shannon_violin_ffmi_v0 /
  simpson_violin_ffmi_v0 /
  richness_violin_ffmi_v0 +
  plot_layout(guides = "collect") &
  theme(legend.position = "none", plot.title = element_blank())
ggsave(alpha_panel_ffmi_v0_vertical, filename = "graphs/alphadiversity/alpha_panel_ffmi_v0_vertical.pdf", width = 5, height = 8)
