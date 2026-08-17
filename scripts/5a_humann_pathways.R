# Functional pathways (HUMAnN)
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(grid)
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

# Data
humann <- readRDS("data/raw_data/BARIA.humann4.profiles.2026.581.910.RDS")
baria_mb <- readRDS("data/processed_data/260817_BARIA_mb_clean.RDS")
baria_muscle_long <- readRDS("data/processed_data/260817_BARIA_muscle_long.RDS")

# Clean up
humann_long <- humann |> 
  rownames_to_column(var = "pathway") |> 
  pivot_longer(
    cols = -pathway,
    names_to = "Sample",
    values_to = "pathway_abundance"
  ) |> 
  mutate(
    Sample = Sample |> 
      str_remove("_Abundance$") |> 
      str_replace("V\\.1", "V-1")
  ) |> 
  filter(
    str_detect(Sample, "\\.Fecal\\."),
    Sample %in% sample_names(baria_mb),
    !pathway %in% c("UNMAPPED", "UNINTEGRATED")
  ) |> 
  mutate(
    pathway_id = str_extract(pathway, "[A-Z0-9-]+(?=:)"),
    pathway_name = str_extract(pathway, "(?<=:).*") 
      |> trimws()
  ) |> 
  left_join(
    as(sample_data(baria_mb), "data.frame") |> 
      rownames_to_column(var = "Sample") |> 
      select(Sample, id, visit),
    by = "Sample"
  ) |> 
  mutate(pathway_abundance = replace_na(pathway_abundance, 0)) |> # Treat undetected pathways as zero abundance
  group_by(pathway_id) |> 
  mutate(
    # Keep pathways with at least > 5 CPM in at least 50% of baseline samples
    prevalence_filter = mean(pathway_abundance[visit == "v0"] > 5, na.rm = TRUE) >= 0.5 
  ) |> 
  ungroup() |> 
  mutate(log10_pathway_abundance = log10(pathway_abundance + 1))


