# Cross-sectional associations between species abundance and FFMI
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(broom)
library(MetBrewer)
library(grid)

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
shape_manual_signif <- scale_shape_manual(
    values = c("significant" = 16, "not significant" = 1), 
    labels = c("significant" = "< 0.05", "not significant" = "\u2265 0.05")
  )
color_manual_models123 <- scale_color_manual(
    values = c(
      "Age and sex" = renoir_15[14],
      "Age, sex and FMI" = renoir_15[6],
      "Age, sex, FMI, T2D, antidiabetic medication and statin use" = renoir_15[2]
  ))

# Data
baria_muscle_wide <- readRDS("data/processed_data/260817_BARIA_muscle_wide.RDS")
baria_mb <- readRDS("data/processed_data/260817_BARIA_mb_clean.RDS")

# Filter out poorly annotated ("GGB"-containing) taxa for species-level associations
baria_mb_species <- prune_taxa(
  str_detect(rownames(otu_table(baria_mb)), "GGB\\d+", negate = TRUE),
  baria_mb
)

# Extract abundance matrix
matrix_mb <- as(otu_table(baria_mb_species), "matrix") |> 
  t()

# Build baseline model df
mb_v0 <- as(sample_data(baria_mb_species), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(Sample, id, visit) |>
  filter(visit == "v0") |>
  left_join(
    baria_muscle_wide |>
      select(id, age_v0, sex, t2d_v0, dm_meds_v0, statins_v0, ffmi_v0, fmi_v0, low_ffmi_v0),
    by = "id"
  ) |>
  left_join(
    matrix_mb |>
      as.data.frame() |>
      rownames_to_column(var = "Sample"),
    by = "Sample"
  )

## Calculate prevalence and mean relative abundance per species
# Use only mb species cols from matrix_mb
species_v0 <- mb_v0 |>
  select(all_of(colnames(matrix_mb))) 

# Prevalence (= proportion of patients/samples where species is detected)
# nts: relative prevalence, not absolute count
species_v0_prevalence <- colMeans(species_v0 > 0, na.rm = FALSE)

# Mean relative abundance per species
species_v0_abundance <- colMeans(species_v0, na.rm = FALSE)

# Combine prevalence and abundance into one filtering overview table
species_v0_fltr_check <- tibble(
  species = names(species_v0_prevalence),
  prevalence = species_v0_prevalence,
  mean_abundance = species_v0_abundance
)

### Filter
# Keep species detected in at least 20% of baseline samples, with mean relative abundance >= 0.01%
species_v0_keep <- species_v0_fltr_check |>
  filter(
    prevalence >= 0.20,
    mean_abundance >= 0.01
  ) |>
  pull(species)

length(species_v0_keep) # check

# Filter baseline species
species_v0_fltr <- species_v0 |>
  select(all_of(species_v0_keep))

# Check
rowSums(species_v0_fltr, na.rm = TRUE) |>
  summary()

##### FFMI #####
## Associations with baseline FFMI ##
# Prepare data for linear models
model_data_ffmi_v0 <- mb_v0 |> 
  select(id, age_v0, sex, t2d_v0, dm_meds_v0, statins_v0, ffmi_v0, fmi_v0, low_ffmi_v0, all_of(species_v0_keep))

# Pivot longer -> reshape to one row per participant-species combination
model_data_ffmi_v0_long <- model_data_ffmi_v0 |> 
  pivot_longer(
    cols = all_of(species_v0_keep),
    names_to = "species",
    values_to = "abundance"
  ) |> 
  arrange(species)

# Log10-transform abundance using a species-specific pseudocount
model_data_ffmi_v0_long_log <- model_data_ffmi_v0_long |>
  group_by(species) |>
  mutate(
    min_abundance = min(abundance[abundance > 0], na.rm = TRUE),
    log10_abundance = log10(abundance + min_abundance / 2)
  ) |>
  ungroup()

# Create consistent cleaned species labels for use across all models
species_labels <- tibble(
  species = sort(unique(model_data_ffmi_v0_long_log$species))
) |>
  mutate(
    species_label = str_extract(species, "s__[^|]+"),
    species_label = str_remove(species_label, "^s__"),
    #species_label = str_remove(species_label, "_SGB\\d+$"),
    species_label = str_replace_all(species_label, "_", " "),

    sgb = str_extract(species, "t__SGB\\d+"),
    sgb = str_remove(sgb, "^t__")
  )

# Nest participant-level data separately for each species
model_data_ffmi_v0_nested <- model_data_ffmi_v0_long_log |> 
  group_by(species) |> 
  nest()

#### Model 1: Adjusted for age and sex #####
# Fit one lm per species
model_ffmi_v0_1 <- model_data_ffmi_v0_nested |> 
  mutate( # add model as a new col to the (nested) df
    model = map( # applies lm to each nested species-specific df
      data,
      \(x) lm(ffmi_v0 ~ log10_abundance + age_v0 + sex, data = x)
    )
  )

# Extract coefficient tables for all species models
model_ffmi_v0_1_tidy <- model_ffmi_v0_1 |> 
  mutate(results = map(model, broom::tidy, conf.int = TRUE))

# Unnest coefficient tables and keep only abundance coefficient, add FDR-adjusted p-values
model_ffmi_v0_1_results <- model_ffmi_v0_1_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(term == "log10_abundance") |> # keep only abundance term
  ungroup() |> # Remove species grouping before FDR correction across all species
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values
    signif = if_else(p_fdr < 0.05, "significant", "not significant") # for plots
  ) |> 
  left_join(species_labels, by = "species") |> 
  mutate(model = "model1")

# Species significantly associated with baseline FFMI after FDR correction
model_ffmi_v0_1_signif <- model_ffmi_v0_1_results |> 
  filter(p_fdr < 0.05) |> 
  arrange(p_fdr) |>
  relocate(species_label, .after = species)

#### Model  2: Adjusted for age, sex and adiposity (FMI) ####
# Fit one lm per species
model_ffmi_v0_2 <- model_data_ffmi_v0_nested |> 
  mutate( # add model as a new col to the (nested) df
    model = map( # applies lm to each nested species-specific df
      data,
      \(x) lm(ffmi_v0 ~ log10_abundance + age_v0 + sex + fmi_v0, data = x)
    )
  )

# Extract coefficient tables for all species models
model_ffmi_v0_2_tidy <- model_ffmi_v0_2 |> 
  mutate(results = map(model, broom::tidy, conf.int = TRUE))

# Unnest coefficient tables and keep only abundance coefficient, add FDR-adjusted p-values
model_ffmi_v0_2_results <- model_ffmi_v0_2_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(term == "log10_abundance") |> # keep only abundance term
  ungroup() |> 
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values
    signif = if_else(p_fdr < 0.05, "significant", "not significant") # for plots
  ) |> 
  left_join(species_labels, by = "species") |> 
  mutate(model = "model2")

# Species significantly associated with baseline FFMI after FDR correction
model_ffmi_v0_2_signif <- model_ffmi_v0_2_results |> 
  filter(p_fdr < 0.05) |> 
  arrange(p_fdr) |>
  relocate(species_label, .after = species)

#### Model 3: Extensive/Sensitivity model; adjusted for age, sex, adiposity/FMI, T2D, antidiabetic medication, statins
# Fit one lm per species
model_ffmi_v0_3 <- model_data_ffmi_v0_nested |> 
  mutate(
    model = map(
      data,
      \(x) lm(ffmi_v0 ~ log10_abundance + age_v0 + sex + fmi_v0 + t2d_v0 + dm_meds_v0 + statins_v0, data = x)
    )
  )

# Extract coefficient tables for all species models
model_ffmi_v0_3_tidy <- model_ffmi_v0_3 |> 
  mutate(results = map(model, broom::tidy, conf.int = TRUE))

# Keep abundance coefficient and calculate FDR-adjusted p-values
model_ffmi_v0_3_results <- model_ffmi_v0_3_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(term == "log10_abundance") |> 
  ungroup() |> 
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"),
    signif = if_else(
      p_fdr < 0.05,
      "significant",
      "not significant"
    )
  ) |> 
  left_join(species_labels, by = "species") |> 
  mutate(model = "model3")

# Species significantly associated with baseline FFMI after FDR correction
model_ffmi_v0_3_signif <- model_ffmi_v0_3_results |> 
  filter(p_fdr < 0.05) |> 
  arrange(p_fdr) |> 
  relocate(species_label, .after = species)

# Number of FDR-significant species per model
c(
  model1 = nrow(model_ffmi_v0_1_signif),
  model2 = nrow(model_ffmi_v0_2_signif),
  model3 = nrow(model_ffmi_v0_3_signif)
)
