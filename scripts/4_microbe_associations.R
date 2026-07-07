# Baria muscle mass & gut microbiota project: Per-microbe associations
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(broom)
library(ggrepel)
library(ggthemes)
library(MetBrewer)

# Theme
manet_15 <- met.brewer("Manet", n = 15)

theme_Publication <- function(base_size = 14, base_family = "sans") {
  
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
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}

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

##### FFMI #####
#### Associations with baseline FFMI ####
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
species1_data <- model_data_v0_long |> 
  filter(species == species_v0_keep[1]) |> 
  select(id, ffmi_v0, age_v0, sex, bmi_v0, species, abundance)

# Run first model (for one species)
species1_model <- lm(
  ffmi_v0 ~ abundance + age_v0 + sex + bmi_v0,
  data = species1_data
)

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

# Extract coefficient tables for all species models
model_v0_tidy <- model_v0 |> 
  mutate(results = map(model, broom::tidy, conf.int = TRUE))

# Unnest coefficient tables and keep only abundance coefficient, add FDR-adjusted p-values
model_v0_results <- model_v0_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(term == "abundance") |> # keep only abundance term
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values

    signif = if_else(p_fdr < 0.05, "significant", "not significant"),
    posneg = if_else(estimate > 0, "positive", "negative"),

    species_label = str_extract(species, "s__[^|]+"), # extract only species part
    species_label = str_remove(species_label, "^s__"), # remove prefix

    strain_label = str_extract(species, "t__[^|]+"),
    strain_label = str_remove(strain_label, "^t__"),

    species_strain_label = paste(species_label, strain_label, sep = " ")
  ) 

# Top 10 significant species
model_v0_top10 <- model_v0_results |> 
  arrange(p_fdr) |> # sort by (ascending) FDR-adjusted p-value
  select(species, estimate, conf.low, conf.high, p.value, p_fdr) |> 
  ungroup() |> 
  slice_head(n = 10) # keep only the top 10 significant species

# Species significantly associated with baseline FFMI after FDR correction
model_v0_signif <- model_v0_results |> 
  filter(p_fdr < 0.05) |> 
  arrange(p_fdr) |>
  relocate(c(species_label, strain_label, species_strain_label), .after = species) |> 
  ungroup()
model_ffm_v4_signif$species_strain_label

### Forest plot
forest_model_v0 <- model_v0_signif |> 
  mutate(species_strain_label = fct_reorder(species_strain_label, estimate)) |> 
  ggplot(aes(x = estimate, y = species_strain_label)) +
  geom_point(aes(color = posneg), size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = posneg), height = 0.2) +
  labs(
    x = "Beta-coefficients (95% CI)",
    y = NULL,
    title = "Species associated with baseline FFMI",
    subtitle = "Linear models adjusted for age, sex, and BMI; FDR < 0.05"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("positive" = manet_15[4], "negative" = manet_15[9]))
ggsave(plot = forest_model_v0, filename = "graphs/microbe_associations/forest_model_v0.pdf", width = 8, height = 6)

### Bar plot
col_model_v0 <- model_v0_signif |> 
  mutate(species_strain_label = fct_reorder(species_strain_label, estimate)) |> 
  ggplot(aes(x = estimate, y = species_strain_label, fill = posneg)) +
  geom_col(width = 0.8) +
  labs(x = "Beta-coefficients") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("positive" = manet_15[4], "negative" = manet_15[9]))
ggsave(plot = col_model_v0, filename = "graphs/microbe_associations/col_model_v0.pdf", width = 8, height = 6)

### Volcano plot
volcano_model_v0 <-  model_v0_results |> 
  ggplot(aes(x = estimate, y = -log10(p_fdr))) +
  geom_point(aes(color = signif), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_label_repel(
    data = filter(model_v0_results, signif == "significant"),
    aes(label = species_strain_label)
  ) +
  labs(x = "Beta-coefficient") +
  scale_color_manual(values = c("significant" = manet_15[9], "not significant" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(plot = volcano_model_v0, filename = "graphs/microbe_associations/volcano_model_v0.pdf", width = 8, height = 6)




##### FFM #######

#### Associations with baseline FFM ####
# Prepare data for linear models
model_data_ffm_v0 <- mb_v0 |> 
  select(id, age_v0, sex, bmi_v0, ffm_kg_v0, all_of(species_v0_keep))

# Pivot longer
model_data_ffm_v0_long <- model_data_ffm_v0 |> 
  pivot_longer(
    cols = all_of(species_v0_keep),
    names_to = "species",
    values_to = "abundance"
  )

# Create one nested df per species to fit the same model repeatedly
model_data_ffm_v0_nested <- model_data_ffm_v0_long |> 
  group_by(species) |> 
  nest()

# Fit one lm per species
model_ffm_v0 <- model_data_ffm_v0_nested |> 
  mutate( # add model as a new col to the (nested) df
    model = map( # applies lm to each nested species-specific df
      data,
      \(x) lm(ffm_kg_v0 ~ abundance + age_v0 + sex + bmi_v0, data = x)
    )
  )

# Extract coefficient tables for all species models
model_ffm_v0_tidy <- model_ffm_v0 |> 
  mutate(results = map(model, broom::tidy, conf.int = TRUE))

# Unnest coefficient tables and keep only abundance coefficient, add FDR-adjusted p-values
model_ffm_v0_results <- model_ffm_v0_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(term == "abundance") |> # keep only abundance term
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values

    signif = if_else(p_fdr < 0.05, "significant", "not significant"),
    posneg = if_else(estimate > 0, "positive", "negative"),

    species_label = str_extract(species, "s__[^|]+"), # extract only species part
    species_label = str_remove(species_label, "^s__"), # remove prefix

    strain_label = str_extract(species, "t__[^|]+"),
    strain_label = str_remove(strain_label, "^t__"),

    species_strain_label = paste(species_label, strain_label, sep = " ")
  ) 

# Top 10 significant species
model_ffm_v0_top10 <- model_ffm_v0_results |> 
  arrange(p_fdr) |> # sort by (ascending) FDR-adjusted p-value
  select(species, estimate, conf.low, conf.high, p.value, p_fdr) |> 
  ungroup() |> 
  slice_head(n = 10) # keep only the top 10 significant species

# Species significantly associated with FFMI after FDR correction
model_ffm_v0_signif <- model_ffm_v0_results |> 
  filter(p_fdr < 0.05) |> 
  arrange(p_fdr) |>
  relocate(c(species_label, strain_label, species_strain_label), .after = species) |> 
  ungroup()

### Forest plot
forest_model_ffm_v0 <- model_ffm_v0_signif |> 
  mutate(species_strain_label = fct_reorder(species_strain_label, estimate)) |> 
  ggplot(aes(x = estimate, y = species_strain_label)) +
  geom_point(aes(color = posneg), size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = posneg), height = 0.2) +
  labs(
    x = "Beta-coefficients (95% CI)",
    y = NULL,
    title = "Species associated with baseline FFM",
    subtitle = "Linear models adjusted for age, sex, and BMI; FDR < 0.05"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("positive" = manet_15[15], "negative" = manet_15[8]))
ggsave(plot = forest_model_ffm_v0, filename = "graphs/microbe_associations/forest_model_ffm_v0.pdf", width = 8, height = 6)

### Bar plot
col_model_ffm_v0 <- model_ffm_v0_signif |> 
  mutate(species_strain_label = fct_reorder(species_strain_label, estimate)) |> 
  ggplot(aes(x = estimate, y = species_strain_label, fill = posneg)) +
  geom_col(width = 0.8) +
  labs(x = "Beta-coefficients") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("positive" = manet_15[15], "negative" = manet_15[8]))
ggsave(plot = col_model_ffm_v0, filename = "graphs/microbe_associations/col_model_ffm_v0.pdf", width = 8, height = 6)

### Volcano plot
volcano_model_ffm_v0 <-  model_ffm_v0_results |> 
  ggplot(aes(x = estimate, y = -log10(p_fdr))) +
  geom_point(aes(color = signif), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_label_repel(
    data = filter(model_ffm_v0_results, signif == "significant"),
    aes(label = species_strain_label)
  ) +
  labs(x = "Beta-coefficient") +
  scale_color_manual(values = c("significant" = manet_15[9], "not significant" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(plot = volcano_model_ffm_v0, filename = "graphs/microbe_associations/volcano_model_ffm_v0.pdf", width = 8, height = 6)