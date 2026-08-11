# Alpha diversity
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

renoir_cols_20 <- met.brewer("Renoir", n = 20)
fill_cols_2 <- scale_fill_manual(values = c("yes" = renoir_cols_20[18], "no" = renoir_cols_20[5]))
x_axis_labels_ffmi <- scale_x_discrete(labels = c("yes" = "Low FFMI", "no" = "Moderate/high FFMI"))

# Data
baria_mb <- readRDS("data/processed_data/260811_BARIA_mb_clean.RDS") # already contains necessary metadata for the grouping

# Diversity metrics
# Convert OTU table to matrix and transpose 
matrix_mb <- as(otu_table(baria_mb), "matrix") |> 
  t() # vegan requires a matrix with samples as rows and taxa as cols

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

# Metadata
alpha_meta <- as(sample_data(baria_mb), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(Sample, id, visit, low_ffmi_v0) |>
  left_join(alpha, by = "Sample")

#### Baseline Plots ####
alpha_v0 <- alpha_meta |>
  filter(visit == "v0") |>
  arrange(id, Sample) |>
  distinct(id, .keep_all = TRUE) |> 
  mutate(low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes"))

## Shannon ##
# Boxplot
shannon_ffmi_v0 <- alpha_v0 |> 
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
  x_axis_labels_ffmi + 
  theme_minimal_custom() +
  theme(legend.position = "none")
ggsave(shannon_violin_ffmi_v0, filename = "graphs/alphadiversity/shannon_violin_ffmi_v0.pdf", width = 6, height = 5)

## Simpson ##
# Boxplot
simpson_ffmi_v0 <- alpha_v0 |> 
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
  x_axis_labels_ffmi +
  theme_minimal_custom() +
  theme(legend.position = "none")
ggsave(simpson_violin_ffmi_v0, filename = "graphs/alphadiversity/simpson_violin_ffmi_v0.pdf", width = 6, height = 5)

## Richness ##
# Boxplot
richness_ffmi_v0 <- alpha_v0 |> 
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
  x_axis_labels_ffmi +
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
ggsave(alpha_panel_ffmi_v0_vertical, filename = "graphs/alphadiversity/alpha_panel_ffmi_v0_vertical.pdf", width = 3.8, height = 8)
