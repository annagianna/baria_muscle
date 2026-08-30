# Functional pathways (HUMAnN)
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(grid)
library(lmerTest)
library(broom.mixed)
library(phyloseq)
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
renoir_15 <- met.brewer("Renoir", n = 15)

# Data
baria_humann <- readRDS("data/raw_data/BARIA.humann4.profiles.2026.581.910.RDS")
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
linda_species_signif <- readRDS("data/processed_data/LinDA_significant_species.RDS")

# HUMAnN cleaning
humann_clean <- baria_humann |>
  rownames_to_column(var = "pathway") |>
  mutate(
    pathway_id = str_extract(pathway, "^[^:]+"),
    pathway_name = str_remove(pathway, "^[^:]+:\\s*"),
    pathway_name = str_replace(pathway_name, "^.", toupper),
    pathway_type = case_when(pathway %in% c("UNMAPPED", "UNINTEGRATED") ~ "non-pathway", TRUE ~ "pathway")
  )

# Long dataset
humann_long <- humann_clean |>
  filter(pathway_type == "pathway") |>
  pivot_longer(
    cols = -c(pathway, pathway_id, pathway_name, pathway_type),
    names_to = "Sample_humann",
    values_to = "pathway_abundance"
  ) |>
  mutate(
    Sample = Sample_humann |>
      str_remove("_Abundance$") |>
      str_replace("\\.V\\.1\\.", ".V-1.") # HUMAnN data baseline V1, shotgun V-1
  ) |>
  inner_join(
    as(sample_data(baria_mb), "data.frame") |>
      rownames_to_column("Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  ) |>
  filter(visit %in% c("v0", "v4", "v5")) |>
  select(Sample, Sample_humann, id, visit, pathway, pathway_id, pathway_name, pathway_abundance)

# Extract significant species from LinDA
linda_species_long <- as(otu_table(baria_mb), "matrix") |>
  as.data.frame() |>
  rownames_to_column(var = "species") |>
  filter(species %in% linda_species_signif) |>
  pivot_longer(
    cols = -species,
    names_to = "Sample",
    values_to = "species_abundance"
  ) |>
  left_join(
    as(sample_data(baria_mb), "data.frame") |>
      rownames_to_column("Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  )



head(

