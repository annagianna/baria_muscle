# Baria muscle mass & gut microbiota project: Associations between baseline species & FFMI trajectories
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(ggdag)
library(lmerTest)
library(broom.mixed)
library(ggrepel)
library(ggthemes)
library(MetBrewer)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)

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
baria_muscle_ab <- read_rds("data/20260722_BARIA_muscle_clinical.RDS") # clinical data
baria_mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")

# Filter out participants taking antibiotics
baria_muscle <- baria_muscle_ab |> 
  filter(abx_v0 == "no")

# Keep only samples with one run or first run of samples with duplicates (same approach as in previous scripts)
run1_mb <- prune_samples(
  sample_data(baria_mb)$Extra_data == "NA" | sample_data(baria_mb)$Extra_data == "rep1",
  baria_mb
)

# Filter out poorly annotated ("GGB"-containing) taxa
run1_mb_clean <- prune_taxa(
  str_detect(rownames(otu_table(run1_mb)), "GGB\\d+", negate = TRUE),
  run1_mb
)

# Extract abundance matrix
matrix_mb <- as(otu_table(run1_mb_clean), "matrix")

# Transponse samples as rows, taxa as cols
if (taxa_are_rows(run1_mb_clean)) {
  matrix_mb <- t(matrix_mb)
}

# Sample data as df
sample_data_df <- as(sample_data(run1_mb_clean), "data.frame") |> 
  rownames_to_column(var = "Sample") |> 
  mutate(
    visit = case_when(
      str_detect(Time_Point, "V-1") ~ "0",
      str_detect(Time_Point, "V4") ~ "4",
      str_detect(Time_Point, "V5") ~ "5",
      TRUE ~ NA_character_
    ),
    visit = as.factor(visit),
    id = Subject_ID
  )

# Join matrix_mb with sample metadata and filter baseline samples
mb_v0 <- matrix_mb |>
  as.data.frame() |>
  rownames_to_column(var = "Sample") |>
  left_join(
    sample_data_df |>
      select(Sample, id, visit),
    by = "Sample"
  ) |> 
  relocate(id, visit, .before = everything()) |> 
  filter(
    visit == "0",
    id %in% baria_muscle$id # already filtered for !is.na(ffmi_v0) in cleaning script + abx use above
  )

species_v0 <- mb_v0 |> 
  select(all_of(colnames(matrix_mb)))

# Compute revalence and mean abundance
species_v0_prevalence <- colMeans(species_v0 > 0, na.rm = FALSE) # (= proportion of patients/samples where species is detected)
species_v0_abundance <- colMeans(species_v0, na.rm = FALSE) # mean relative abundance per species

## Filter
# Keep species detected in at least 20% of baseline samples, with mean relative abundance >= 0.01%
species_v0_keep <- tibble(
  species = names(species_v0_prevalence),
  prevalence_v0 = species_v0_prevalence,
  mean_abundance_v0 = species_v0_abundance
) |>
  filter(
    prevalence_v0 >= 0.20,
    mean_abundance_v0 >= 0.01
  ) |>
  pull(species)

# Pivot clinical data longer
baria_muscle_long <- baria_muscle |>
  select(id, age_v0, sex, ffmi_v0, ffmi_v4, ffmi_v5, fmi_v0, t2d_v0, dm_meds_v0, statins_v0) |>
  #filter(!is.na(ffmi_v4), !is.na(ffmi_v5)) |> 
  pivot_longer(
    cols = c(ffmi_v0, ffmi_v4, ffmi_v5),
    names_to = "visit",
    names_prefix = "ffmi_v",
    values_to = "ffmi"
  ) |>
  mutate(visit = factor(visit, levels = c("0", "4", "5"))) |> 
  filter(!is.na(ffmi))

# Model data (long format) (merge long metadata with mb + prevalence + abundance matrices)
lmm_data_ffmi <- mb_v0 |>
  pivot_longer(
    cols = all_of(species_v0_keep),
    names_to = "species",
    values_to = "baseline_abundance"
  ) |> 
  select(-visit, -Sample) |> 
  inner_join(baria_muscle_long, by = "id") |> 
  filter(!is.na(species))

# Log10-transform abundance and calculate a species-specific pseudocount
lmm_data_ffmi_log10 <- lmm_data_ffmi |>
  group_by(species) |> 
  mutate(
    min_baseline_abundance = min(baseline_abundance[baseline_abundance > 0], na.rm = TRUE),
    log10_baseline_abundance = log10(baseline_abundance + min_baseline_abundance / 2)
  ) |> 
  ungroup()

# Create consistent species labels for use across all models
species_labels <- tibble(species = sort(unique(lmm_data_ffmi_log10$species))) |>
  mutate(
    species_label = str_extract(species, "s__[^|]+"),
    species_label = str_remove(species_label, "^s__"),
    #species_label = str_remove(species_label, "_SGB\\d+$"),
    species_label = str_replace_all(species_label, "_", " "),

    sgb = str_extract(species, "t__SGB\\d+"),
    sgb = str_remove(sgb, "^t__")
  )

# Nest by species
lmm_data_ffmi_nested <- lmm_data_ffmi_log10 |>
  group_by(species) |>
  nest()


#### Model 1: Age, sex ####
# Fit one mixed model per species
lmm1_ffmi <- lmm_data_ffmi_nested |>
  mutate(
    model = map(
      data,
      ~ lmerTest::lmer(ffmi ~ log10_baseline_abundance * visit + age_v0 + sex + (1 | id), data = .x, REML = FALSE)
    )
  )

# extract coeff tables for each model
lmm1_ffmi_tidy <- lmm1_ffmi |> 
  mutate(results = map(model, broom.mixed::tidy, effects = "fixed", conf.int = TRUE))


lmm1_ffmi_results <- lmm1_ffmi_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(str_detect(term, "log10_baseline_abundance:visit")) |> 
  group_by(term) |> # separately for v4 and v5
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values
    signif = if_else(p_fdr < 0.05, "significant", "not significant") # for plots
  ) |> 
  ungroup() |> 
  mutate(follow_up = case_when(term == "log10_baseline_abundance:visit4" ~ "1 year", term == "log10_baseline_abundance:visit5" ~ "2 years")) |> 
  left_join(species_labels, by = "species")

# Signif. results for 1y (v4)
lmm1_ffmi_v4_signif <- lmm1_ffmi_results |> 
  filter(term == "log10_baseline_abundance:visit4", signif == "significant")

# Signif. results for 2y (v5)
lmm1_ffmi_v5_signif <- lmm1_ffmi_results |> 
  filter(term == "log10_baseline_abundance:visit5", signif == "significant")