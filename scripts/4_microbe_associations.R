# Baria muscle mass & gut microbiota project: Cross-sectional per-microbe associations
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(dagitty)
library(ggdag)
library(broom)
library(ggrepel)
library(ggthemes)
library(MetBrewer)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)

theme_Publication <- function(base_size = 14, base_family = "sans") {
  
  (theme_foundation(base_size = base_size, base_family = base_family) + 
    theme(
      plot.title = element_text(face = "bold", size = rel(0.8), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA, fill = NA),
      plot.background = element_rect(colour = NA, fill = NA),
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
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}

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
baria_muscle_ab <- read_rds("data/20260720_BARIA_muscle_clinical.RDS") # clinical data
baria_mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")

# qc
sample_sums(baria_mb) |>
  summary() # adds up to 100

# Filter out participants taking antibiotics
baria_muscle <- baria_muscle_ab |> 
  filter(abx_v0 == "no")

# Keep only samples with one run or first run of samples with duplicates (same approach as in previous scripts)
run1_mb <- prune_samples(
  sample_data(baria_mb)$Extra_data == "NA" | sample_data(baria_mb)$Extra_data == "rep1",
  baria_mb
)

# Filter out poorly annotated ("GGB"-containing) taxa for species-level associations
run1_mb_clean <- prune_taxa(
  str_detect(rownames(otu_table(run1_mb)), "GGB\\d+", negate = TRUE),
  run1_mb
)

# Extract abundance matrix
matrix_mb <- as(otu_table(run1_mb_clean), "matrix")

# vegan requires a matrix with samples as rows, taxa as cols
if (taxa_are_rows(run1_mb_clean)) {
  matrix_mb <- t(matrix_mb)
}

# Check that each sample sums to ~100% to confirm that table sums rel. abundances
rowSums(matrix_mb, na.rm = TRUE) |> 
  summary()

# Create analysis table: abundance matrix + mb sample data + clinical metadata
# Sample data as df
run1_mb_clean_df <- as(sample_data(run1_mb_clean), "data.frame") |> 
  tibble::rownames_to_column(var = "Sample") |> 
  mutate(
    visit = case_when(
      str_detect(Time_Point, "V-1") ~ "0",
      str_detect(Time_Point, "V4") ~ "4",
      str_detect(Time_Point, "V5") ~ "5",
      TRUE ~ NA_character_
    ),
    visit = as.factor(visit),
    id = Subject_ID
  ) |> 
  select(-Time_Point, -Subject_ID)

# Join all together
mb_df <- matrix_mb |> # abundance matrix
  as.data.frame() |> 
  tibble::rownames_to_column(var = "Sample") |> 
  left_join(
    run1_mb_clean_df |> # mb sample data
      select(Sample, visit, id),
    by = "Sample"
  ) |> 
  left_join(baria_muscle, by = "id") |> 
  relocate(id, visit, .before = everything())

# Keep baseline samples only for cross-sectional associations (baseline species abundance with FFMI)
mb_v0 <- mb_df |> 
  filter(visit == "0")

# Check baseline sample size and availability of FFMI vars
mb_v0 |> 
  summarize(
    n_samples = n(),
    n_ffmi_v0 = sum(!is.na(ffmi_v0)),
    n_low_ffmi_v0 = sum(!is.na(low_ffmi_v0))
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
### Draw DAG ###
# Model 1: Primary DAG
dag_ffmi_v0_1 <- dagitty::dagitty("
dag {
microbiome [exposure]
FFMI [outcome]

age -> microbiome
age -> FFMI

sex -> microbiome
sex -> FFMI

microbiome -> metabolites
metabolites -> FFMI
}")

dagitty::adjustmentSets(
  dag_ffmi_v0_1,
  exposure = "microbiome",
  outcome = "FFMI",
  type = "minimal",
  effect = "total"
)

dagitty::paths(
  dag_ffmi_v0_1,
  from = "microbiome",
  to = "FFMI"
)

# Model 2: Sensitivity DAG (with adjustment for adiposity/FMI)
dag_ffmi_v0_2s <- dagitty::dagitty("
dag {
microbiome [exposure]
FFMI [outcome]

age -> microbiome
age -> FFMI

sex -> microbiome
sex -> FFMI

adiposity -> microbiome
adiposity -> FFMI

microbiome -> metabolites
metabolites -> FFMI
}")

dagitty::adjustmentSets(
  dag_ffmi_v0_2s,
  exposure = "microbiome",
  outcome = "FFMI",
  type = "minimal",
  effect = "total"
)

# Model 3: Extensive / sensitivity DAG 
# Additional adjustment for T2D and antidiabetic medication
dag_ffmi_v0_3s <- dagitty::dagitty("
dag {
microbiome [exposure]
FFMI [outcome]

age -> microbiome
age -> FFMI
age -> T2D

sex -> microbiome
sex -> FFMI
sex -> T2D

adiposity -> microbiome
adiposity -> FFMI
adiposity -> T2D
adiposity -> inflammation
adiposity -> insulin_resistance

T2D -> microbiome
T2D -> FFMI
T2D -> metformin

metformin -> microbiome
metformin -> FFMI

microbiome -> inflammation
inflammation -> FFMI

microbiome -> FFMI

FFMI -> insulin_resistance
}")

dagitty::adjustmentSets(
  dag_ffmi_v0_3s,
  exposure = "microbiome",
  outcome = "FFMI",
  type = "minimal",
  effect = "total"
)

# Prepare data for linear models
model_data_ffmi_v0 <- mb_v0 |> 
  select(id, age_v0, sex, t2d_v0, dm_meds_v0, metformin_v0, ffmi_v0, fmi_v0, low_ffmi_v0, all_of(species_v0_keep))

# Pivot longer -> reshape to one row per participant-species combination
model_data_ffmi_v0_long <- model_data_ffmi_v0 |> 
  pivot_longer(
    cols = all_of(species_v0_keep),
    names_to = "species",
    values_to = "abundance"
  ) |> 
  arrange(species) # ensure consistent order before make.unique()

# Find the smallest non-zero relative abundance
min_abundance <- model_data_ffmi_v0_long |> 
  filter(abundance > 0) |> 
  summarize(min_abundance = min(abundance)) |> 
  pull(min_abundance)

# Log10-transform abundance using half the minimum non-zero abundance as pseudocount
model_data_ffmi_v0_long_log <- model_data_ffmi_v0_long |> 
  mutate(log10_abundance = log10(abundance + (min_abundance/2)))

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

# Create consistent cleaned species labels for use across all models
species_labels <- tibble(
  species = sort(unique(model_data_ffmi_v0_long_log$species))
) |>
  mutate(
    species_label = str_extract(species, "s__[^|]+"),
    species_label = str_remove(species_label, "^s__"),
    species_label = str_remove(species_label, "_SGB\\d+$"),
    species_label = str_replace_all(species_label, "_", " ")
  )

# Unnest coefficient tables and keep only abundance coefficient, add FDR-adjusted p-values
model_ffmi_v0_1_results <- model_ffmi_v0_1_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(term == "log10_abundance") |> # keep only abundance term
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
  relocate(species_label, .after = species) |> 
  ungroup()

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
  relocate(species_label, .after = species) |> 
  ungroup()

# Prepare data and labels to plot models 1 & 2 together
# Species significant after FDR correction in at least one model (1 and/or 2)
models_ffmi_1_2_signif <- union(model_ffmi_v0_1_signif$species, model_ffmi_v0_2_signif$species)

# Create unique labels only among (significant) species included in the combined figure
models_ffmi_1_2_species_labels <- species_labels |> 
  filter(species %in% models_ffmi_1_2_signif) |> 
  arrange(species_label, species) |> 
  mutate(
    species_label_unique = make.unique(species_label)
  ) |> 
  select(species, species_label_unique)

# Combine estimates from both models for species significant in at least one model
models_ffmi_v0_1_2_results <- bind_rows(model_ffmi_v0_1_results, model_ffmi_v0_2_results) |> 
  filter(species %in% models_ffmi_1_2_signif) |> 
  select(-species_label) |> 
  left_join(models_ffmi_1_2_species_labels, by = "species") |> # use unique labels
  mutate(model = factor(model, levels = c("model1", "model2"), labels = c("Age and sex", "Age, sex and FMI")))

# Order species according to model 1 estimates
species_order_model1 <- models_ffmi_v0_1_2_results |> 
  filter(model == "Age and sex") |> 
  arrange(estimate) |> 
  pull(species_label_unique)

### Combined forest plot for models 1 and 2 
forest_model_ffmi_1_2 <- models_ffmi_v0_1_2_results |> 
  mutate(
    species_label_unique = factor(species_label_unique, levels = species_order_model1),
    signif = factor(signif, levels = c("significant", "not significant"))
  ) |> 
  ggplot(aes(x = estimate, y = species_label_unique, color = model, shape = signif, group = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high), 
    height = 0.5, 
    position = position_dodge(width = 0.6),
    show.legend = FALSE
  ) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("Age and sex" = renoir_15[14], "Age, sex and FMI" = renoir_15[6])) +
  scale_shape_manual(values = c("significant" = 16, "not significant" = 1), labels = c("significant" = "< 0.05", "not significant" = "≥ 0.05")) +
  labs(
    x = "kg/m² change in baseline FFMI per 1-unit increase in log10 abundance",
    y = NULL,
    colour = NULL,
    shape = "FDR p-value"
  ) +
  theme_minimal_custom() +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "italic"))
ggsave(plot = forest_model_ffmi_1_2, filename = "graphs/microbe_associations/forest_model_ffmi_1_2.pdf", width = 8, height = 6)

#### Model 3: Extensive/Sensitivity model; adjusted for age, sex adiposity, dm, antidiabetic medication
# Check numbers for T2D and antidiabetic medication
model_data_ffmi_v0 |> 
  count(t2d_v0, dm_meds_v0)

#### Model 3: Adjusted for age, sex, adiposity, T2D and antidiabetic medication ####
# Fit one lm per species
model_ffmi_v0_3 <- model_data_ffmi_v0_nested |> 
  mutate(
    model = map(
      data,
      \(x) lm(ffmi_v0 ~ log10_abundance + age_v0 + sex + fmi_v0 + t2d_v0 + metformin_v0, data = x)
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
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"),
    signif = if_else(
      p_fdr < 0.05,
      "significant",
      "not significant"
    )
  ) |> 
  left_join(species_labels, by = "species") #|> 
 # mutate(posneg = if_else(estimate > 0, "positive", "negative"))

# Species significantly associated with baseline FFMI after FDR correction
model_ffmi_v0_3_signif <- model_ffmi_v0_3_results |> 
  filter(p_fdr < 0.05) |> 
  arrange(p_fdr) |> 
  relocate(species_label, .after = species) |> 
  ungroup()

# Create unique labels among species significant in Model 3
model_ffmi_v0_3_species_labels <- model_ffmi_v0_3_signif |> 
  select(species, species_label) |> 
  arrange(species_label, species) |> 
  mutate(species_label_unique = make.unique(species_label)) |> 
  select(species, species_label_unique)


# Forest plot model 3
forest_model_ffmi_v0_3 <- model_ffmi_v0_3_signif |>
  select(-species_label) |> 
  left_join(model_ffmi_v0_3_species_labels, by = "species") |> 
  arrange(estimate) |> 
  mutate(species_label_unique = factor(species_label_unique, levels = species_label_unique)) |> 
  ggplot(aes(x = estimate, y = species_label_unique)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = renoir_15[8]), height = 0.5) +
  geom_point(size = 2.5, color = renoir_15[8]) +
  labs(
    x = "kg/m² change in baseline FFMI per 1-unit increase in log10 abundance",
    y = NULL,
    title = "Species associated with baseline FFMI",
    subtitle = "Adjusted for age, sex, FMI, T2D and metformin; FDR < 0.05"
  ) +
  theme_minimal_custom() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic")
  )
ggsave(plot = forest_model_ffmi_v0_3, filename = "graphs/microbe_associations/forest_model_ffmi_v0_3.pdf", width = 8, height = 6)

