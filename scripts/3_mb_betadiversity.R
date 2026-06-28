# Baria muscle loss & microbiota project: Beta diversity
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
library(ggsci)

# Theme
manet_cols_30 <- met.brewer("Manet", n = 30)
fill_cols_2 <- scale_fill_manual(values = c("yes" = manet_cols_30[14], "no" = manet_cols_30[30]))
color_cols_2 <- scale_color_manual(values = c("yes" = manet_cols_30[14], "no" = manet_cols_30[30]))

theme_Publication <- function(base_size=14, base_family="sans") {
  
  (theme_foundation(base_size = base_size, base_family = base_family) + 
    theme(
      plot.title = element_text(face = "bold", size = rel(0.8), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(0.8)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(), 
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.key.size = unit(0.2, "cm"),
      legend.spacing = unit(0, "cm"),
      plot.margin = unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
} 

# Data
baria_muscle <- read_rds("data/20260624_BARIA_muscle_clinical.RDS") # clinical data
baria_mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")

# qc
sample_sums(baria_mb) |>
  summary() # adds up to 100

# Subset; keep only samples with one run or first run of samples with duplicates
run1_mb <- prune_samples(
  sample_data(baria_mb)$Extra_data == "NA" | sample_data(baria_mb)$Extra_data == "rep1",
  baria_mb
)

# Convert OTU table to matrix
matrix_mb <- as(otu_table(run1_mb), "matrix")

# vegan requires a matrix with samples as rows and taxa as cols
if (taxa_are_rows(run1_mb)) { 
  matrix_mb <- t(matrix_mb) 
}

# checks
dim(matrix_mb) # should be n samples (rows) x n taxa (cols)
rownames(matrix_mb)[1:5] # Samples
colnames(matrix_mb)[1:5] # Taxa
summary(rowSums(matrix_mb)) # should be ~100

# Convert sample data into a df
mb_df <- sample_data(run1_mb) |> 
  data.frame() |> 
  tibble::rownames_to_column(var = "Sample")

# Merge with metadata
mb_meta <- mb_df |> 
  mutate(
    visit = str_extract(Time_Point, "\\d"),
    visit = if_else(visit == "1", "0", visit),
    visit = as.factor(visit),
    id = Subject_ID
  ) |> 
  select(Sample, id, visit) |> 
  inner_join(baria_muscle, by = "id") |> 
  relocate(id, .before = Sample)

# Keep only baseline samples for baseline beta diversity analysis
mb_meta_v0 <- mb_meta |> 
  filter(visit == "0")

# Check whether all samples from metadata are in the (baseline) mb matrix
all(mb_meta_v0$Sample %in% rownames(matrix_mb)) # returns TRUE

# Keep only samples that also have matched metadata; order to match mb_meta exactly for downstream analyses
mb_v0 <- matrix_mb[mb_meta_v0$Sample, , drop = FALSE]

# Final check
all(rownames(mb_v0) == mb_meta_v0$Sample)

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
  tibble::rownames_to_column("Sample") |>
  left_join(mb_meta_v0, by = "Sample") # merge with ordered metadata 

## PERMANOVA
# Test whether overall mb composition differs by baseline low FFMI group
adonis_v0 <- vegan::adonis2(bray_v0 ~ low_ffmi_v0, data = mb_meta_v0, permutations = 999, by = "terms")
p_adonis_v0 <- adonis_v0$`Pr(>F)`[1]

# Check homogeneity of dispertion (variances)
disp_v0 <- vegan::betadisper(bray_v0, group = mb_meta_v0$low_ffmi_v0)
anova(disp_v0)

## Plots ##
# Baseline
centroids_v0 <- bray_2_v0 |> 
  group_by(low_ffmi_v0) |> 
  summarize(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE),
    .groups = "drop"
  )

bray_ffmi_v0 <- bray_2_v0 |> 
  mutate(low_ffmi_v0 = factor(low_ffmi_v0, levels = c("yes", "no"))) |> 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(color = low_ffmi_v0), size = 3, alpha = 0.9) +
  xlab(paste0("PCo1 (", round(expl_variance_bray_v0[1], digits = 1),"%)")) +
  ylab(paste0("PCo2 (", round(expl_variance_bray_v0[2], digits = 1),"%)")) +
  labs(color = "", fill = "", title = "PCoA Bray-Curtis Distance") +
  stat_ellipse(
    geom = "polygon", 
    aes(color = low_ffmi_v0, fill = low_ffmi_v0),
    type = "norm", alpha = 0.13, linewidth = 1) +
  geom_point(
    data = centroids_v0, 
    aes(x = Axis.1, y = Axis.2),
    size = 6, shape = 21
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("PERMANOVA p = ", round(p_adonis_v0, 3)),
    hjust = 1.0, vjust = 1.2,
    size = 3
  ) +
  fill_cols_2 +
  color_cols_2 +
  theme_Publication() +
  theme(legend.position = "top") 
ggsave(bray_ffmi_v0, filename = "graphs/betadiversity/bray_ffmi_v0.pdf", width = 10, height = 7)