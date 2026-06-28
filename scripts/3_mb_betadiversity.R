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
manet_cols <- met.brewer("Manet", n = 20)
fill_cols_asm1y <- scale_fill_manual(values = c("high" = manet_cols[10], "low/modest" = manet_cols[20]))
color_cols_asm1y <- scale_color_manual(values = c("high" = manet_cols[10], "low/modest" = manet_cols[20]))

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
bray_v0_2 <- as.data.frame(pcoord_v0$vectors[, c("Axis.1", "Axis.2")]) |>
  tibble::rownames_to_column("Sample") |>
  left_join(mb_meta_v0, by = "Sample") # merge with ordered metadata 






### Bray-Curtis per visit ###
# set up function
bray_per_v <- function(v, df = bray_meta, abd_tab = bray_wide) {
  set.seed(1312)

  meta_v <- df |> 
    filter(visit == v) |> 
    distinct(Sample, .keep_all = TRUE) # drop all other cols

  abd_tab_v <- abd_tab[meta_v$Sample, , drop = FALSE] # keep only Samples from that visit 

  bray_v <- vegan::vegdist(abd_tab_v, method = "bray", na.rm = TRUE)

  vegan::adonis2(bray_v ~ asm_change_v4_group, data = meta_v, permutations = 999, by = "terms")
}

# run it for each visit
visits <- sort(unique(bray_meta$visit)) 
results_per_v <- lapply(visits, bray_per_v)
names(results_per_v) <- visits

## Plots ##
# Baseline
bray_asm1y_v0 <- bray2 |> 
  filter(visit == 0) |> 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(color = asm_change_v4_group), size = 3, alpha = 0.9) +
  xlab(paste0("PCo1 (", round(expl_variance_bray[1], digits = 1),"%)")) +
  ylab(paste0("PCo2 (", round(expl_variance_bray[2], digits = 1),"%)")) +
  labs(color = "", fill = "", title = "PCoA Bray-Curtis Distance") +
  stat_ellipse(
    geom = "polygon", 
    aes(color = asm_change_v4_group, fill = asm_change_v4_group),
    type = "norm", alpha = 0.13, linewidth = 1.2) +
  fill_cols_asm1y +
  color_cols_asm1y +
  theme_Publication() +
  theme(legend.position = "top") +
  geom_text(
    data = as.data.frame(results_per_v$`0`) |> slice(1),
    aes(x = Inf, y = Inf,
      label = paste0("p = ", round(`Pr(>F)`, 3))),
      hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE)
ggsave(bray_asm1y_v0, filename = "graphs/betadiversity/bray_asm1y_v0.pdf", width = 10, height = 7)
