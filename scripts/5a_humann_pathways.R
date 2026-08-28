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
humann_long <- readRDS("data/processed_data/BARIA_humann_pathways_long.RDS")
baria_muscle_long <- readRDS("data/processed_data/BARIA_muscle_long.RDS")

# Feature selection & transform for this analysis
humann_long <- humann_long |>
  group_by(pathway_id) |>
  filter(mean(pathway_abundance[visit == "v0"] > 5) >= 0.2) |>
  ungroup() |>
  mutate(log10_pathway_abundance = log10(pathway_abundance + 1))

#### Baseline pathway - FFMI associations ####
# Prepare baseline data
humann_v0 <- humann_long |>
  filter(visit == "v0") |>
  left_join(
    baria_muscle_long |>
      filter(visit == "v0") |>
      select(id, age_v0, sex, t2d_v0, dm_meds, statins, ffmi, fmi) |>
      rename(dm_meds_v0 = dm_meds, statins_v0 = statins, ffmi_v0 = ffmi, fmi_v0 = fmi),
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
    select(pathway_id, pathway_name, results) |>
    unnest(results) |>
    filter(term == "log10_pathway_abundance") |>
    ungroup() |> # Remove pathway grouping before FDR correction across all pathways
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
# Merge baseline pathway abundance with longitudinal FFMI data
# Prepare LMM data
model_data_pathway_lmm <- humann_v0 |>
  select(id, pathway_id, pathway_name, baseline_pathway_abundance = pathway_abundance, log10_baseline_pathway_abundance = log10_pathway_abundance) |>
  inner_join(
    baria_muscle_long |>
      filter(visit %in% c("v0", "v4", "v5"), !is.na(ffmi)) |>
      mutate(
        visit = factor(visit, levels = c("v0", "v4", "v5")),
        age_centered_v0 = age_v0 - mean(age_v0, na.rm = TRUE)
      ) |>
      select(id, visit, ffmi, age_centered_v0, sex, perc_change_weight_kg),
    by = "id",
    relationship = "many-to-many" # expected: each id fans out across pathway x visit
  )

# Function for pathway LMMs
run_pathway_lmm <- function(lmm_data, follow_up) {

  lmm_data |>
    filter(visit %in% c("v0", follow_up)) |>
    droplevels() |>
    group_by(pathway_id, pathway_name) |>
    nest() |>
    mutate(
      model = map(
        data,
        ~ lmerTest::lmer(ffmi ~ log10_baseline_pathway_abundance * visit + age_centered_v0 + sex + perc_change_weight_kg + (1 | id),
        data = .x, REML = FALSE
        )
      ),
      results = map(model, ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE))
    ) |>
    select(pathway_id, pathway_name, results) |>
    unnest(results) |>
    filter(str_detect(term, "log10_baseline_pathway_abundance:visit")) |>
    ungroup() |> # Remove pathway grouping before FDR correction across all pathways
    mutate(
      p_fdr = p.adjust(p.value, method = "BH"),
      signif = if_else(p_fdr < 0.05, "significant", "not significant"),
      visit = follow_up
    )
}

# Run pathway LMM for each follow-up visit
lmm_pathway_v4_results <- run_pathway_lmm(model_data_pathway_lmm, "v4")
lmm_pathway_v5_results <- run_pathway_lmm(model_data_pathway_lmm, "v5")

# Significant pathways
lmm_pathway_v4_signif <- lmm_pathway_v4_results |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

lmm_pathway_v5_signif <- lmm_pathway_v5_results |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

c(
  n_pathways_tested = n_distinct(model_data_pathway_lmm$pathway_id),
  v4_results = nrow(lmm_pathway_v4_results),
  v5_results = nrow(lmm_pathway_v5_results),
  v4_signif = nrow(lmm_pathway_v4_signif),
  v5_signif = nrow(lmm_pathway_v5_signif)
)

pathway_lmm_nominal_overlap <- lmm_pathway_v4_results |>
  filter(p.value < 0.05) |>
  select(
    pathway_id,
    pathway_name,
    estimate_v4 = estimate,
    p_v4 = p.value,
    p_fdr_v4 = p_fdr
  ) |>
  inner_join(
    lmm_pathway_v5_results |>
      filter(p.value < 0.05) |>
      select(
        pathway_id,
        estimate_v5 = estimate,
        p_v5 = p.value,
        p_fdr_v5 = p_fdr
      ),
    by = "pathway_id"
  ) |>
  mutate(
    same_direction = sign(estimate_v4) == sign(estimate_v5)
  )

nrow(pathway_lmm_nominal_overlap)
pathway_lmm_nominal_overlap
