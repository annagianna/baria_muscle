# Baria muscle mass & gut microbiota project: Per-microbe associations
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(broom)

# Data
baria_muscle <- read_rds("data/20260624_BARIA_muscle_clinical.RDS") # clinical data
baria_mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")

# qc
sample_sums(baria_mb) |>
  summary() # adds up to 100

# Keep only samples with one run or first run of samples with duplicates (same approach as in previous scripts)
run1_mb <- prune_samples(
  sample_data(baria_mb)$Extra_data == "NA" | sample_data(baria_mb)$Extra_data == "rep1",
  baria_mb
)

# Extract abundance matrix
matrix_mb <- as(otu_table(run1_mb), "matrix")

# vegan requires a matrix with samples as rows, taxa as cols
if (taxa_are_rows(run1_mb)) {
  matrix_mb <- t(matrix_mb)
}

# Check that each sample sums to ~100% to confirm that table sums rel. abundances
rowSums(matrix_mb, na.rm = TRUE) |> 
  summary()

# Create analysis table: abundance matrix + mb sample data + clinical metadata
# Sample data as df
run1_mb_df <- as(sample_data(run1_mb), "data.frame") |> 
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
    run1_mb_df |> # mb sample data
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
# Use only mb species cols
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
# Keep species detected in at least 20% of baseline samples, with mean rel. abundance >= 0.01%
species_v0_keep <- species_v0_fltr_check |>
  filter(
    prevalence >= 0.20,
    mean_abundance >= 0.01
  ) |>
  pull(species)

length(species_v0_keep)

# Filter baseline species
species_v0_fltr <- species_v0 |>
  select(all_of(species_v0_keep))

# Check
rowSums(species_v0_fltr, na.rm = TRUE) |>
  summary()

#### Associations
# Prepare data for linear models
model_data_v0 <- mb_v0 |> 
  select(id, age_v0, sex, bmi_v0, ffmi_v0, low_ffmi_v0, all_of(species_v0_keep))

# Pivot longer
model_data_v0_long <- model_data_v0 |> 
  pivot_longer(
    cols = all_of(species_v0_keep),
    names_to = "species",
    values_to = "abundance"
  )

## Test for one species
species1_data <- model_data_ffmi_v0_long |> 
  filter(species == species_v0_keep[1]) |> 
  select(id, ffmi_v0, age_v0, sex, bmi_v0, species, abundance)

# Run first model (for one species)
species1_model <- lm(
  ffmi_v0 ~ abundance + age_v0 + sex + bmi_v0,
  data = species1_data
)

summary(species1_model)

## All species
# Create one nested df per species to fit the same model repeatedly
model_data_v0_nested <- model_data_v0_long |> 
  group_by(species) |> 
  nest()

# Fit one lm per species
model_v0 <- model_data_v0_nested |> 
  mutate( # add model as a new col to the (nested) df
    model = map( # applies lm to each nested species-specific df
      data,
      \(x) lm(ffmi_v0 ~ abundance + age_v0 + sex + bmi_v0, data = x)
    )
  )







### Low FFMI 


