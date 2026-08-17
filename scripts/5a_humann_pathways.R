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
renoir_15 <- met.brewer("Renoir", n = 15)

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
  filter(mean(pathway_abundance[visit == "v0"] > 5) >= 0.5) |> # Keep pathways with >5 CPM in at least 50% of baseline samples
  ungroup() |> 
  mutate(log10_pathway_abundance = log10(pathway_abundance + 1))

#### Baseline pathway - FFMI associations ####
# Prepare baseline data
humann_v0 <- humann_long |>
  filter(visit == "v0") |>
  left_join(
    baria_muscle_long |>
      filter(visit == "v0") |>
      select(id, age_v0, sex, t2d_v0, dm_meds_v0, statins_v0, ffmi, fmi) |>
      rename(ffmi_v0 = ffmi, fmi_v0 = fmi),
    by = "id"
  )

# Nest participant-level data separately for each pathway
model_data_ffmi_v0_nested <- humann_v0 |>
  group_by(pathway_id, pathway_name) |>
  nest()

# Function to run pathway LMs
run_pathway_lm <- function(model_data, model_formula) {
  
  model_data |>
    mutate(
      model = map(data, \(x) lm(model_formula, data = x)),
      results = map(model, broom::tidy, conf.int = TRUE)
    ) |>
    ungroup() |> # Remove pathway grouping before FDR correction across all pathways
    select(pathway_id, pathway_name, results) |>
    unnest(results) |>
    filter(term == "log10_pathway_abundance") |>
    mutate(
      p_fdr = p.adjust(p.value, method = "BH"),
      signif = if_else(p_fdr < 0.05, "significant", "not significant")
    )
}

# Model formulas
model1 <- ffmi_v0 ~ log10_pathway_abundance + age_v0 + sex
model2 <- ffmi_v0 ~ log10_pathway_abundance + age_v0 + sex + fmi_v0
model3 <- ffmi_v0 ~ log10_pathway_abundance + age_v0 + sex + fmi_v0 + t2d_v0 + dm_meds_v0 + statins_v0

# Run pathway LM for each model formula
model_pathway_ffmi_v0_1_results <- run_pathway_lm(model_data_ffmi_v0_nested, model1)
model_pathway_ffmi_v0_2_results <- run_pathway_lm(model_data_ffmi_v0_nested, model2)
model_pathway_ffmi_v0_3_results <- run_pathway_lm(model_data_ffmi_v0_nested, model3)

# Pathways significantly associated with baseline FFMI after FDR correction
model_pathway_ffmi_v0_1_signif <- model_pathway_ffmi_v0_1_results |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

model_pathway_ffmi_v0_2_signif <- model_pathway_ffmi_v0_2_results |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

model_pathway_ffmi_v0_3_signif <- model_pathway_ffmi_v0_3_results |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

# No pathways were significantly associated with baseline FFMI after FDR correction

#### Baseline pathway abundance & FFMI trajectories LMMs ####
# Prepare baseline pathway abundance
humann_v0 <- humann_long |>
  left_join(
    baria_muscle_long |>
      filter(visit == "v0") |>
      select(id, age_v0, sex, t2d_v0, dm_meds_v0, statins_v0, ffmi, fmi) |>
      rename(ffmi_v0 = ffmi, fmi_v0 = fmi),
    by = "id"
  ) 