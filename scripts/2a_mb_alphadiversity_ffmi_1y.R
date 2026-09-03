# Alpha diversity - baseline microbiome by 1-year FFMI loss group
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
      panel.grid.major = element_line(
        colour = "#dddddd",
        linewidth = 0.4,
        linetype = "22"
      ),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "bottom",
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )
}

renoir_cols_20 <- met.brewer("Renoir", n = 20)

fill_cols_2 <- scale_fill_manual(
  name = "FFMI loss at 1 year",
  values = c(
    "high" = renoir_cols_20[18],
    "modest/low" = renoir_cols_20[5]
  )
)

x_axis_labels_loss <- scale_x_discrete(
  labels = c(
    "high" = "High",
    "modest/low" = "Modest/low"
  )
)

# Data
baria_mb <- readRDS(
  "data/processed_data/BARIA_mb_clean.RDS"
)

# Output folder
dir.create("results/graphs/alphadiversity", recursive = TRUE, showWarnings = FALSE)

# Diversity metrics
# Convert OTU table to matrix with samples as rows and taxa as columns
matrix_mb <- as(otu_table(baria_mb), "matrix")

if (taxa_are_rows(baria_mb)) {
  matrix_mb <- t(matrix_mb)
}

# Shannon, Simpson and richness
shannon <- vegan::diversity(matrix_mb, index = "shannon")
simpson <- vegan::diversity(matrix_mb, index = "simpson")
richness <- vegan::specnumber(matrix_mb)

# Create diversity dataframe
alpha <- tibble(
  Sample = names(shannon),
  shannon = shannon,
  simpson = simpson,
  richness = richness
)

# Add metadata
alpha_meta <- as(sample_data(baria_mb), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(
    Sample,
    id,
    visit,
    perc_change_ffmi_v4_group
  ) |>
  left_join(alpha, by = "Sample")

#### Baseline ####

alpha_v0 <- alpha_meta |>
  filter(
    visit == "v0",
    !is.na(perc_change_ffmi_v4_group)
  ) |>
  arrange(id, Sample) |>
  distinct(id, .keep_all = TRUE) |>
  mutate(
    perc_change_ffmi_v4_group = factor(
      perc_change_ffmi_v4_group,
      levels = c("high", "modest/low")
    )
  )

# Check sample numbers
alpha_v0 |>
  count(perc_change_ffmi_v4_group)


#### Shannon ####

# Boxplot
shannon_ffmi_loss_v0 <- alpha_v0 |>
  ggplot(
    aes(
      x = perc_change_ffmi_v4_group,
      y = shannon,
      fill = perc_change_ffmi_v4_group
    )
  ) +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  labs(
    title = "Shannon index",
    y = "Shannon index",
    x = "1-year FFMI loss",
    fill = "1-year FFMI loss"
  ) +
  fill_cols_2 +
  x_axis_labels_loss +
  theme_minimal_custom()

ggsave(
  filename = "results/graphs/alphadiversity/shannon_by_ffmi_loss_v4.pdf",
  plot = shannon_ffmi_loss_v0,
  width = 7,
  height = 5
)

# Violin
shannon_violin_ffmi_loss_v0 <- alpha_v0 |>
  ggplot(
    aes(
      x = perc_change_ffmi_v4_group,
      y = shannon
    )
  ) +
  geom_violin(
    aes(fill = perc_change_ffmi_v4_group),
    trim = FALSE
  ) +
  geom_boxplot(
    fill = "white",
    width = 0.1
  ) +
  stat_compare_means(
    tip.length = 0,
    hide.ns = TRUE,
    label.x = 1.5,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(
    x = NULL,
    y = "Shannon index",
    title = "Shannon index",
    fill = "1-year FFMI loss group"
  ) +
  fill_cols_2 +
  x_axis_labels_loss +
  theme_minimal_custom() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "results/graphs/alphadiversity/shannon_violin_by_ffmi_loss_v4.pdf",
  plot = shannon_violin_ffmi_loss_v0,
  width = 6,
  height = 5
)


#### Simpson ####

# Boxplot
simpson_ffmi_loss_v0 <- alpha_v0 |>
  ggplot(
    aes(
      x = perc_change_ffmi_v4_group,
      y = simpson,
      fill = perc_change_ffmi_v4_group
    )
  ) +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  labs(
    title = "Simpson index",
    y = "Simpson index",
    x = "1-year FFMI loss",
    fill = "1-year FFMI loss"
  ) +
  fill_cols_2 +
  x_axis_labels_loss +
  theme_minimal_custom()

ggsave(
  filename = "results/graphs/alphadiversity/simpson_by_ffmi_loss_v4.pdf",
  plot = simpson_ffmi_loss_v0,
  width = 7,
  height = 5
)

# Violin
simpson_violin_ffmi_loss_v0 <- alpha_v0 |>
  ggplot(
    aes(
      x = perc_change_ffmi_v4_group,
      y = simpson
    )
  ) +
  geom_violin(
    aes(fill = perc_change_ffmi_v4_group),
    trim = FALSE
  ) +
  geom_boxplot(
    fill = "white",
    width = 0.1
  ) +
  stat_compare_means(
    tip.length = 0,
    hide.ns = TRUE,
    label.x = 1.5,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(
    x = NULL,
    y = "Simpson index",
    title = "Simpson index",
    fill = "1-year FFMI loss group"
  ) +
  fill_cols_2 +
  x_axis_labels_loss +
  theme_minimal_custom() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "results/graphs/alphadiversity/simpson_violin_by_ffmi_loss_v4.pdf",
  plot = simpson_violin_ffmi_loss_v0,
  width = 6,
  height = 5
)


#### Richness ####

# Boxplot
richness_ffmi_loss_v0 <- alpha_v0 |>
  ggplot(
    aes(
      x = perc_change_ffmi_v4_group,
      y = richness,
      fill = perc_change_ffmi_v4_group
    )
  ) +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  labs(
    title = "Richness",
    y = "Richness",
    x = "1-year FFMI loss",
    fill = "1-year FFMI loss"
  ) +
  fill_cols_2 +
  x_axis_labels_loss +
  theme_minimal_custom()

ggsave(
  filename = "results/graphs/alphadiversity/richness_by_ffmi_loss_v4.pdf",
  plot = richness_ffmi_loss_v0,
  width = 7,
  height = 5
)

# Violin
richness_violin_ffmi_loss_v0 <- alpha_v0 |>
  ggplot(
    aes(
      x = perc_change_ffmi_v4_group,
      y = richness
    )
  ) +
  geom_violin(
    aes(fill = perc_change_ffmi_v4_group),
    trim = FALSE
  ) +
  geom_boxplot(
    fill = "white",
    width = 0.1
  ) +
  stat_compare_means(
    tip.length = 0,
    hide.ns = TRUE,
    label.x = 1.5,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(
    x = NULL,
    y = "Richness",
    title = "Richness",
    fill = "1-year FFMI loss group"
  ) +
  fill_cols_2 +
  x_axis_labels_loss +
  theme_minimal_custom() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "results/graphs/alphadiversity/richness_violin_by_ffmi_loss_v4.pdf",
  plot = richness_violin_ffmi_loss_v0,
  width = 6,
  height = 5
)


#### Combined panels ####

# Horizontal
alpha_panel_ffmi_loss_v0 <-
  (
    shannon_violin_ffmi_loss_v0 +
      simpson_violin_ffmi_loss_v0 +
      richness_violin_ffmi_loss_v0
  ) +
  plot_layout(guides = "collect") &
  theme(aspect.ratio = 0.8)

ggsave(
  filename = "results/graphs/alphadiversity/alpha_panel_by_ffmi_loss_v4.pdf",
  plot = alpha_panel_ffmi_loss_v0,
  width = 12,
  height = 8
)

# Vertical
alpha_panel_ffmi_loss_v0_vertical <-
  shannon_violin_ffmi_loss_v0 /
  simpson_violin_ffmi_loss_v0 /
  richness_violin_ffmi_loss_v0 +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    plot.title = element_blank()
  )

ggsave(
  filename = "results/graphs/alphadiversity/alpha_panel_by_ffmi_loss_v4_vertical.pdf",
  plot = alpha_panel_ffmi_loss_v0_vertical,
  width = 3.8,
  height = 8
)