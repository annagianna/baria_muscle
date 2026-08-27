# Differential abundance analysis with LinDA
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(MicrobiomeStat)
library(phyloseq)
library(grid)
library(MetBrewer)

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

# Data
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
baria_muscle_long <- readRDS("data/processed_data/BARIA_muscle_long.RDS")

### Prepare data for LinDA ###
# Filter out poorly annotated ("GGB"-containing) taxa
baria_mb_species <- prune_taxa(
  str_detect(rownames(otu_table(baria_mb)), "GGB\\d+", negate = TRUE),
  baria_mb
)

# Extract abundance matrix & metadata
mb_otu <- as(otu_table(baria_mb_species), "matrix") # taxa are rows
mb_sample_data <- as(sample_data(baria_mb_species), "data.frame")

# Join metadata from long dataset
linda_meta_long <- mb_sample_data |>
  filter(visit %in% c("v0", "v4", "v5")) |>
  rownames_to_column("Sample") |>
  left_join(
    baria_muscle_long |> 
      filter(visit %in% c("v0", "v4", "v5")) |> 
      select(id, visit, age_v0, ffmi, fmi, perc_change_weight_kg, t2d_v0, t2d_labs, dm_meds),
    by = c("visit", "id")
  ) |> 
  filter(!is.na(ffmi), !is.na(perc_change_weight_kg)) |> 
  group_by(id) |> 
  mutate(
    # Create between and within vars for LinDA
    ffmi_between = mean(ffmi),
    ffmi_within = ffmi - ffmi_between
  ) |> 
  ungroup() |> 
  mutate(
    visit = factor(visit, levels = c("v0", "v4", "v5")),
    id = factor(id)
  )

# Match abundance matrix to long metadata
mb_otu_long <- mb_otu[ , linda_meta_long$Sample, drop = FALSE]

# Compute prevalence and mean abundance
species_prevalence <- rowMeans(mb_otu_long > 0)
species_abundance <- rowMeans(mb_otu_long)

## Filter
# Keep species detected in at least 50% of samples, with mean relative abundance >= 0.01%
species_keep <- tibble(
  species = names(species_prevalence),
  prevalence = species_prevalence,
  mean_abundance = species_abundance
) |>
  filter(
    prevalence >= 0.50,
    mean_abundance >= 0.01
  ) |>
  pull(species)

length(species_keep)

# Filter abundance matrix
mb_otu_keep <- mb_otu_long[species_keep, , drop = FALSE]

#### Run LinDA ####
linda_ffmi <- linda(
  feature.dat = mb_otu_keep,
  meta.dat = linda_meta_long,
  formula = "~ ffmi_within + ffmi_between + visit + age_v0 + sex + perc_change_weight_kg + (1 | id)",
  feature.dat.type = "proportion",
  prev.filter = 0, mean.abund.filter = 0, # already set filters above
  p.adj.method = "BH",
  alpha = 0.05
)
