# Beta diversity
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(MetBrewer)
library(grid)
library(ggthemes)
library(ggpubr)
library(ggsci)

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
fill_cols_2 <- scale_fill_manual(
  values = c("yes" = renoir_cols_20[18], "no" = renoir_cols_20[5]),
  labels = c("yes" = "Low FFMI", "no" = "Moderate/high FFMI")
)
color_cols_2 <- scale_color_manual(
  values = c("yes" = renoir_cols_20[18], "no" = renoir_cols_20[5]),
  labels = c("yes" = "Low FFMI", "no" = "Moderate/high FFMI")
)

# Data
baria_mb <- readRDS("data/processed_data/260811_BARIA_mb_clean.RDS") # already contains necessary metadata for the grouping

# Convert OTU table to matrix and transpose 
matrix_mb <- as(otu_table(baria_mb), "matrix") |> 
  t() # vegan requires a matrix with samples as rows and taxa as cols

# Metadata (baseline only)
mb_meta_v0 <- as(sample_data(baria_mb), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(Sample, id, visit, low_ffmi_v0) |> 
  filter(visit == "v0")

# Keep baseline samples and order matrix to match metadata
mb_v0 <- matrix_mb[mb_meta_v0$Sample, , drop = FALSE] # keeps the result as a matrix

# Check exact alignment
stopifnot(all(rownames(mb_v0) == mb_meta_v0$Sample))

### Compute Bray-Curtis distance ###
bray_v0 <- vegan::vegdist(mb_v0, method = "bray")

# PCoA
pcoord_v0 <- ape::pcoa(bray_v0, correction = "cailliez") # corrects for negative eigenvalues

# Percentage of variation explained by each PCoA axis
expl_variance_bray_v0 <- pcoord_v0$values$Rel_corr_eig * 100

# Check how much variation is explained by PCoA1 and 2
expl_variance_bray_v0[1:2]

# Extract first 2 PCoA axes and add Sample names
bray_2_v0 <- as.data.frame(pcoord_v0$vectors[, c("Axis.1", "Axis.2")]) |>
  rownames_to_column("Sample") |>
  left_join(mb_meta_v0, by = "Sample") # merge with ordered metadata 

## PERMANOVA
# Test whether overall mb composition differs by baseline low FFMI group
adonis_v0 <- vegan::adonis2(bray_v0 ~ low_ffmi_v0, data = mb_meta_v0, permutations = 999, by = "terms")
p_adonis_v0 <- adonis_v0$`Pr(>F)`[1]
r2_adonis_v0 <- adonis_v0$R2[1]

# Check homogeneity of dispersion (variances)
disp_v0 <- vegan::betadisper(bray_v0, group = mb_meta_v0$low_ffmi_v0)
disp_anova_v0 <- anova(disp_v0)
p_disp_anova_v0 <- disp_anova_v0$`Pr(>F)`[1]

## Plots ##
bray_ffmi_v0 <- bray_2_v0 |> 
  mutate(low_ffmi_v0 = factor(low_ffmi_v0, levels = c("yes", "no"))) |> 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(color = low_ffmi_v0), size = 2.5, alpha = 0.8) +
  xlab(paste0("PCoA1 (", round(expl_variance_bray_v0[1], digits = 1),"%)")) +
  ylab(paste0("PCoA2 (", round(expl_variance_bray_v0[2], digits = 1),"%)")) +
  labs(color = "", fill = "", title = "PCoA Bray-Curtis Distance") +
  stat_ellipse(
    geom = "polygon", 
    aes(color = low_ffmi_v0, fill = low_ffmi_v0),
    type = "norm", alpha = 0.13, linewidth = 1,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0(
      "PERMANOVA",
      "\nR² = ", formatC(r2_adonis_v0, format = "f", digits = 3),
      ", p = ", formatC(p_adonis_v0, format = "f", digits = 3),
      "\nDispersion p = ", formatC(p_disp_anova_v0, format = "f", digits = 3)
  ),
  hjust = 1.0, vjust = 1.2,
  size = 3
  ) +
  fill_cols_2 +
  color_cols_2 +
  theme_minimal_custom() +
  theme(legend.position = "left") 
ggsave(bray_ffmi_v0, filename = "graphs/betadiversity/bray_ffmi_v0.pdf", width = 10, height = 7)