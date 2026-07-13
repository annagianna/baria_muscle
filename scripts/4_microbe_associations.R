# Baria muscle mass & gut microbiota project: Per-microbe associations
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
# Keep species detected in at least 10% of baseline samples, with mean relative abundance >= 0.01%
species_v0_keep <- species_v0_fltr_check |>
  filter(
    prevalence >= 0.10,
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
#### Associations with baseline FFMI ####
### Draw DAG ###
# Primary DAG
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

# Sensitivity DAG (with adjustment for adiposity/FMI)
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

# Prepare data for linear models
model_data_ffmi_v0 <- mb_v0 |> 
  select(id, age_v0, sex, ffmi_v0, fmi_v0, low_ffmi_v0, all_of(species_v0_keep))

# Pivot longer
model_data_ffmi_v0_long <- model_data_ffmi_v0 |> 
  pivot_longer(
    cols = all_of(species_v0_keep),
    names_to = "species",
    values_to = "abundance"
  ) |> 
  arrange(species) # to make the order deterministic before make.unique()

# Create one nested df per species to fit the same model repeatedly
model_data_ffmi_v0_nested <- model_data_ffmi_v0_long |> 
  group_by(species) |> 
  nest()

# First model: adjusted for age and sex
# Fit one lm per species
model_ffmi_v0_1 <- model_data_ffmi_v0_nested |> 
  mutate( # add model as a new col to the (nested) df
    model = map( # applies lm to each nested species-specific df
      data,
      \(x) lm(ffmi_v0 ~ abundance + age_v0 + sex, data = x)
    )
  )

# Extract coefficient tables for all species models
model_ffmi_v0_1_tidy <- model_ffmi_v0_1 |> 
  mutate(results = map(model, broom::tidy, conf.int = TRUE))

# Unnest coefficient tables and keep only abundance coefficient, add FDR-adjusted p-values
model_ffmi_v0_1_results <- model_ffmi_v0_1_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(term == "abundance") |> # keep only abundance term
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values

    signif = if_else(p_fdr < 0.05, "significant", "not significant"), # for plots
    posneg = if_else(estimate > 0, "positive", "negative"),

    species_label = str_extract(species, "s__[^|]+"), # extract only species part
    species_label = str_remove(species_label, "^s__"), # remove prefix
    species_label = str_remove(species_label, "_SGB\\d+$"), # remove suffix / SGB identifier
    species_label = str_replace_all(species_label, "_", " "),
    species_label = make.unique(species_label)
  )

# Species significantly associated with baseline FFMI after FDR correction
model_ffmi_v0_1_signif <- model_ffmi_v0_1_results |> 
  filter(p_fdr < 0.05) |> 
  arrange(p_fdr) |>
  relocate(species_label, .after = species) |> 
  ungroup()

### Forest plot
forest_model_ffmi_v0_1 <- model_ffmi_v0_1_signif |> 
  mutate(species_label = fct_reorder(species_label, estimate)) |> 
  ggplot(aes(x = estimate, y = species_label)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = posneg), size = 2.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = posneg), height = 0.5) +
  labs(
    x = "Estimate (95% CI)",
    y = NULL,
    title = "Species associated with baseline FFMI",
    subtitle = "Linear models adjusted for age and sex; FDR < 0.05"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic")
  ) +
  scale_color_manual(values = c("positive" = renoir_15[14], "negative" = renoir_15[6]))
ggsave(plot = forest_model_ffmi_v0_1, filename = "graphs/microbe_associations/forest_model_ffmi_v0_1.pdf", width = 8, height = 6)

### Bar plot
col_model_ffmi_v0_1 <- model_ffmi_v0_1_signif |> 
  mutate(species_label = fct_reorder(species_label, estimate)) |> 
  ggplot(aes(x = estimate, y = species_label, fill = posneg)) +
  geom_col(width = 0.8) +
  labs(x = "Estimate") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic")
  ) +
  scale_fill_manual(values = c("positive" = renoir_15[14], "negative" = renoir_15[6]))
ggsave(plot = col_model_ffmi_v0_1, filename = "graphs/microbe_associations/col_model_ffmi_v0_1.pdf", width = 8, height = 6)

### Volcano plot
volcano_model_ffmi_v0_1 <-  model_ffmi_v0_1_results |> 
  ggplot(aes(x = estimate, y = -log10(p_fdr))) +
  geom_point(aes(color = signif), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_label_repel(
    data = filter(model_ffmi_v0_1_results, signif == "significant"),
    aes(label = species_label),
    fontface = "italic"
  ) +
  labs(x = "Beta-coefficient") +
  scale_color_manual(values = c("significant" = renoir_15[6], "not significant" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(plot = volcano_model_ffmi_v0_1, filename = "graphs/microbe_associations/volcano_model_ffmi_v0_1.pdf", width = 8, height = 6)
