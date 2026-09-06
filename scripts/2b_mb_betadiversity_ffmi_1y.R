# Beta diversity — Baseline microbiome stratified by 1-year FFMI loss group
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)
library(MetBrewer)
library(grid)
library(ggthemes)
library(ggpubr)
library(ggsci)

# Theme
source("scripts/assets/functions.R")
renoir_cols_20 <- met.brewer("Renoir", n = 20)
fill_cols_2 <- scale_fill_manual(
  name = "FFMI loss at 1 year",
  values = c("high" = renoir_cols_20[18], "modest/low" = renoir_cols_20[5]),
  labels = c("high" = "High", "modest/low" = "Modest/low")
)
color_cols_2 <- scale_color_manual(
  name = "FFMI loss at 1 year",
  values = c("high" = renoir_cols_20[18],"modest/low" = renoir_cols_20[5]),
  labels = c("high" = "High", "modest/low" = "Modest/low")
)

# Data
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")

# Output folder
dir.create("results/graphs/betadiversity", recursive = TRUE, showWarnings = FALSE)

# Convert OTU table to matrix
matrix_mb <- as(otu_table(baria_mb), "matrix") |> 
  t() # vegan requires samples as rows and taxa as columns

# Metadata
mb_meta_v0 <- as(sample_data(baria_mb), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(Sample, id, visit, perc_change_ffmi_v4_group) |>
  filter(visit == "v0", !is.na(perc_change_ffmi_v4_group)) |>
  arrange(id, Sample) |>
  distinct(id, .keep_all = TRUE) |>
  mutate(perc_change_ffmi_v4_group = factor(perc_change_ffmi_v4_group, levels = c("high", "modest/low")))

# Check group sizes
mb_meta_v0 |>
  count(perc_change_ffmi_v4_group)

# Keep baseline samples and order matrix to match metadata
mb_v0 <- matrix_mb[mb_meta_v0$Sample, , drop = FALSE]

# Check exact alignment
stopifnot(all(rownames(mb_v0) == mb_meta_v0$Sample))

#### Bray-Curtis distance ####
bray_v0 <- vegan::vegdist(mb_v0, method = "bray")

#### PCoA ####
pcoord_v0 <- ape::pcoa(bray_v0, correction = "cailliez")
expl_variance_bray_v0 <- pcoord_v0$values$Rel_corr_eig * 100 # Percentage of variation explained
expl_variance_bray_v0[1:2] # Check first two axes

# Extract first two PCoA axes and add metadata
bray_2_v0 <- as.data.frame(pcoord_v0$vectors[, c("Axis.1", "Axis.2")]) |>
  rownames_to_column("Sample") |>
  left_join(mb_meta_v0, by = "Sample")

#### PERMANOVA ####
# Test whether overall microbiome composition differs by 1-year FFMI loss group
adonis_v0 <- vegan::adonis2(bray_v0 ~ perc_change_ffmi_v4_group,data = mb_meta_v0, permutations = 999)
p_adonis_v0 <- adonis_v0$`Pr(>F)`[1]
r2_adonis_v0 <- adonis_v0$R2[1]

#### Homogeneity of dispersion ####
disp_v0 <- vegan::betadisper(bray_v0,group = mb_meta_v0$perc_change_ffmi_v4_group)
disp_anova_v0 <- anova(disp_v0)
p_disp_anova_v0 <- disp_anova_v0$`Pr(>F)`[1]

#### Plot ####
bray_ffmi_loss_1y <- bray_2_v0 |>
  ggplot(aes(x = Axis.1, y = Axis.2)) +
  stat_ellipse(
    geom = "polygon", aes(color = perc_change_ffmi_v4_group, fill = perc_change_ffmi_v4_group),
    type = "norm", alpha = 0.13, linewidth = 1, show.legend = FALSE
  ) +
  geom_point(aes(color = perc_change_ffmi_v4_group), size = 2.5, alpha = 0.8) +
  labs(
    title = "PCoA Bray-Curtis distance",
    x = paste0("PCoA1 (", round(expl_variance_bray_v0[1], 1), "%)"),
    y = paste0("PCoA2 (", round(expl_variance_bray_v0[2], 1), "%)")
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0(
      "PERMANOVA",
      "\nR² = ", formatC(r2_adonis_v0, format = "f", digits = 3), ", p = ", formatC(p_adonis_v0, format = "f", digits = 3),
      "\nDispersion p = ", formatC(p_disp_anova_v0, format = "f", digits = 3)
    ),
    hjust = 1.05, vjust = 1.2, size = 3
  ) +
  fill_cols_2 +
  color_cols_2 +
  theme_minimal_custom() +
  theme(legend.position = "left")
ggsave(filename = "results/graphs/betadiversity/bray_ffmi_loss_1y.pdf", plot = bray_ffmi_loss_1y, width = 10, height = 7)
