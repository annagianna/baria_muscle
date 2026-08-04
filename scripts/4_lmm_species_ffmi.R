# Baria muscle mass & gut microbiota project: Associations between baseline species & FFMI trajectories
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(broom.mixed)
library(ggrepel)
library(ggthemes)
library(ggpubr)
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
baria_muscle <- read_rds("data/20260731_BARIA_muscle_clinical.RDS") # clinical data
baria_mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")

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

# Transpose samples as rows, taxa as cols
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

# Compute prevalence and mean abundance
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
  mutate(age_centered_v0 = age_v0 - mean(age_v0, na.rm = TRUE)) |> 
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

# 1y model data
lmm_data_ffmi_v4 <- lmm_data_ffmi_log10 |> 
  filter(visit %in% c("0", "4")) |> 
  droplevels()

# 2y model data
lmm_data_ffmi_v5 <- lmm_data_ffmi_log10 |> 
  filter(visit %in% c("0", "5")) |> 
  droplevels()

# Nest by species
# 1y
lmm_data_ffmi_v4_nested <- lmm_data_ffmi_v4 |>
  group_by(species) |>
  nest()

# 2y
lmm_data_ffmi_v5_nested <- lmm_data_ffmi_v5 |>
  group_by(species) |>
  nest()

#### Model 1: Age, sex ####
### Baseline -> 1y (v0-v4) ###
lmm1_ffmi_v4 <- lmm_data_ffmi_v4_nested |>
  mutate(
    model = map(
      data,
      ~ lmerTest::lmer(ffmi ~ log10_baseline_abundance * visit + age_centered_v0 + sex + (1 | id), data = .x, REML = FALSE)
    )
  )

# extract coeff tables for each model
lmm1_ffmi_v4_tidy <- lmm1_ffmi_v4 |> 
  mutate(results = map(model, broom.mixed::tidy, effects = "fixed", conf.int = TRUE))

lmm1_ffmi_v4_results <- lmm1_ffmi_v4_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(str_detect(term, "log10_baseline_abundance:visit")) |> 
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values
    signif = if_else(p_fdr < 0.05, "significant", "not significant") # for plots
  ) |>  
  left_join(species_labels, by = "species")

### Baseline -> 2y (v0-v5) ###
lmm1_ffmi_v5 <- lmm_data_ffmi_v5_nested |>
  mutate(
    model = map(
      data,
      ~ lmerTest::lmer(ffmi ~ log10_baseline_abundance * visit + age_centered_v0 + sex + (1 | id), data = .x, REML = FALSE)
    )
  )

# extract coeff tables for each model
lmm1_ffmi_v5_tidy <- lmm1_ffmi_v5 |> 
  mutate(results = map(model, broom.mixed::tidy, effects = "fixed", conf.int = TRUE))

lmm1_ffmi_v5_results <- lmm1_ffmi_v5_tidy |> 
  select(species, results) |> 
  unnest(results) |> 
  filter(str_detect(term, "log10_baseline_abundance:visit")) |> 
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"), # add FDR-adjusted p-values
    signif = if_else(p_fdr < 0.05, "significant", "not significant") # for plots
  ) |>
  left_join(species_labels, by = "species")

# Combine v4 and v5 results
lmm1_ffmi_results <- bind_rows(lmm1_ffmi_v4_results, lmm1_ffmi_v5_results) |> 
  mutate(
    follow_up = case_when(
      term == "log10_baseline_abundance:visit4" ~ "Baseline -> 1 year", 
      term == "log10_baseline_abundance:visit5" ~ "Baseline -> 2 years"
    )
  )

# Significant species in model 1
lmm1_ffmi_signif_species <- lmm1_ffmi_results |> 
  filter(signif == "significant") |> 
  distinct(species) |> 
  pull(species)

# Create unique species labels for significant species in model 1
lmm1_ffmi_species_labels <- species_labels |> 
  filter(species %in% lmm1_ffmi_signif_species) |> 
  arrange(species_label, species) |> 
  group_by(species_label) |> 
  mutate(
    species_label_unique = if(n() > 1) {
      paste(species_label, sgb, sep = " ")
    } else {
      species_label
    }
  ) |> 
  ungroup() |> 
  select(species, species_label_unique)

# Add labels to signif results for each interval
# Species significantly associated with FFMI trajectory from baseline to year 1
lmm1_ffmi_v4_signif <- lmm1_ffmi_v4_results |> 
  filter(signif == "significant") |> 
  arrange(p_fdr)
  
# Species significantly associated with FFMI trajectory from baseline to year 2
lmm1_ffmi_v5_signif <- lmm1_ffmi_v5_results |> 
  filter(signif == "significant") |> 
  arrange(p_fdr)

# Overlap of significant species between y1 & y2 in model1 (age, sex)
lmm1_ffmi_overlap <- lmm1_ffmi_v4_signif |> 
  select(
    species, 
    term_v4 = term, 
    estimate_v4 = estimate,
    p.value_v4 = p.value,
    p_fdr_v4 = p_fdr
  ) |> 
  inner_join(
    lmm1_ffmi_v5_signif |> 
      select(
        species, 
        term_v5 = term, 
        estimate_v5 = estimate,
        p.value_v5 = p.value,
        p_fdr_v5 = p_fdr
      ),
      by = "species"
  ) |> 
  left_join(lmm1_ffmi_species_labels, by = "species") |> 
  relocate(species_label_unique, .after = species) |> 
  mutate(
    posneg_v4 = if_else(estimate_v4 > 0, "positive", "negative"),
    posneg_v5 = if_else(estimate_v5 > 0, "positive", "negative"),
    consistent_direction = posneg_v4 == posneg_v5
  ) |> 
  arrange(p_fdr_v4, p_fdr_v5)

lmm1_ffmi_overlap_species <- lmm1_ffmi_overlap |> 
  pull(species)

# Species-specific abundance percentiles (for plotting)
lmm1_ffmi_percentiles <- lmm_data_ffmi_log10 |> 
  filter(species %in% lmm1_ffmi_overlap_species) |> 
  distinct(species, id, .keep_all = TRUE) |> 
  group_by(species) |> 
  summarize(
    log10_baseline_abundance_p10 = quantile(log10_baseline_abundance, probs = 0.1, na.rm = TRUE),
    log10_baseline_abundance_p90 = quantile(log10_baseline_abundance, probs = 0.9, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  left_join(
    lmm1_ffmi_overlap |> 
      select(species, species_label_unique),
    by = "species"
  )


#### Plot observed FFMI trajectories ####
plot_ffmi_observed <- baria_muscle_long |> 
  select(id, visit, ffmi) |>
  mutate(visit = factor(visit, levels = c("0", "4", "5"), labels = c("Baseline", "1 year", "2 years"))) |> 
  ggplot(aes(x = visit, y = ffmi, group = id)) +

  # individual observed trajectories
  geom_line(color = "grey70", alpha = 0.2, linewidth = 0.6) +
  geom_point(color = "grey70", alpha = 0.2, size = 0.6) +

  # group mean trajectories
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1.2) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", color = "black", size = 1.5) +
  labs(
    title = "Observed FFMI trajectories",
    x = NULL,
    y = "FFMI (kg/m²)"
  ) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.01))) +
  theme_minimal_custom()
ggsave(plot = plot_ffmi_observed, filename = "graphs/lmm_species_ffmi/plot_ffmi_observed.pdf", width = 10, height = 8)


#### Plot function for significant species in 1y & 2y models #####
plot_lmm1_ffmi <- function(species_name, model_df, model_data) {
  
  # Extract the model for one species
  species_model <- model_df |> 
    filter(species == species_name) |> 
    pull(model) |> 
    first()

  # Identify whether supplied data are v4 or v5
  follow_up <- model_data |> 
    filter(visit != "0") |> 
    distinct(visit) |> 
    pull() |> 
    as.character()

  follow_up_label <- case_when(
    follow_up == "4" ~ "1 year",
    follow_up == "5" ~ "2 years",
    TRUE ~ NA_character_
  )

  # Interaction terms for individual species plot titles
  interaction_term <- paste0("log10_baseline_abundance:visit", follow_up)
  interaction_beta <- broom.mixed::tidy(species_model, effects = "fixed") |> 
    filter(term == interaction_term) |> 
    pull(estimate)

  # Observed species data
  species_data <- model_data |> 
    filter(species == species_name) |>
    left_join(lmm1_ffmi_percentiles, by = "species")

  # Population-averaged predictions for P10 & P90
  species_predicted_p10_p90 <- species_data |> 
    select(id, visit, age_centered_v0, sex, log10_baseline_abundance_p10, log10_baseline_abundance_p90) |> # retain only cols needed for model
    pivot_longer(
      cols = c(log10_baseline_abundance_p10, log10_baseline_abundance_p90),
      names_to = "percentile",
      values_to = "log10_baseline_abundance"
    ) |> 
    mutate(
      percentile = str_extract(percentile, "p\\d+"),
      predicted_ffmi = predict(species_model, newdata = pick(everything()), re.form = NA)
    ) |> 
    group_by(visit, percentile) |> 
    summarize(predicted_ffmi = mean(predicted_ffmi, na.rm = TRUE), .groups = "drop") |> 
    mutate(visit = factor(visit, levels = c("0", follow_up), labels = c("Baseline", follow_up_label)
    )
  )

  species_data |> 
    mutate(visit = factor(visit, levels = c("0", follow_up), labels = c("Baseline", follow_up_label))) |> 
    ggplot(aes(x = visit, y = ffmi, group = id)) +

    # observed data (individual participant trajectories)
    geom_line(color = "grey70", alpha = 0.2, linewidth = 0.6) +
    geom_point(color = "grey70", alpha = 0.2, size = 0.6) +

    # predicted data for P10 & P90
    geom_line(
     data = species_predicted_p10_p90,
      aes(x = visit, y = predicted_ffmi, group = percentile, color = percentile),
      linewidth = 1.4
    ) +
    scale_colour_manual(
      values = c("p10" = renoir_15[6], "p90" = renoir_15[14]),
      labels = c("p10" = "P10", "p90" = "P90"),
      name = "Baseline abundance"
    ) +
    
    # Observed data (population average)
    # stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1, color = "black", linetype = "13", alpha = 0.8) +
    # stat_summary(aes(group = 1), fun = mean, geom = "point", size = 0.6, color = "black") +
    
    labs(
      title = unique(species_data$species_label_unique),
      x = NULL,
      y = "FFMI (kg/m²)"
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("βinteraction = ", round(interaction_beta, 2)), 
      hjust = 1, vjust = 1.4, size = 3
    ) +
    scale_x_discrete(expand = expansion(mult = c(0, 0.01))) +
    theme_minimal_custom() +
    theme(plot.title = element_text(face = "italic", hjust = 0.5))
  
}

# Call plot function for each one of the overlap species
# 1 year
plot_lmm1_ffmi_v4 <- map(lmm1_ffmi_overlap_species, ~ plot_lmm1_ffmi(species_name = .x, model_df = lmm1_ffmi_v4, model_data = lmm_data_ffmi_v4))

# 2 years
plot_lmm1_ffmi_v5 <- map(lmm1_ffmi_overlap_species, ~ plot_lmm1_ffmi(species_name = .x, model_df = lmm1_ffmi_v5, model_data = lmm_data_ffmi_v5))

# Panels
# 1 year
panel_lmm1_ffmi_v4 <- ggarrange(
  plotlist = plot_lmm1_ffmi_v4,
  ncol = 3, nrow = 2,
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)
ggsave(plot = panel_lmm1_ffmi_v4, filename = "graphs/lmm_species_ffmi/panel_lmm1_ffmi_v4.pdf", width = 14, height = 8)

# 2 years
panel_lmm1_ffmi_v5 <- ggarrange(
  plotlist = plot_lmm1_ffmi_v5,
  ncol = 3, nrow = 2,
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)
ggsave(plot = panel_lmm1_ffmi_v5, filename = "graphs/lmm_species_ffmi/panel_lmm1_ffmi_v5.pdf", width = 14, height = 8)

#### Species significant in only one time interval #####
# 1 year
lmm1_ffmi_v4_only <- lmm1_ffmi_v4_signif |>
  anti_join(lmm1_ffmi_v5_signif, by = "species") |> 
  arrange(-estimate)

llm1_ffmi_v4_top_pos <- lmm1_ffmi_v4_signif |> 
  slice_head(n = 3)

llm1_ffmi_v4_top_neg <- lmm1_ffmi_v4_signif |> 
  slice_tail(n = 3)

# 2 years
lmm1_ffmi_v5_only <- lmm1_ffmi_v5_signif |> 
  anti_join(lmm1_ffmi_v4_signif, by = "species") |> 
  arrange(-estimate)

llm1_ffmi_v5_top_pos <- lmm1_ffmi_v5_signif |> 
  slice_head(n = 3)

llm1_ffmi_v5_top_neg <- lmm1_ffmi_v5_signif |> 
  slice_tail(n = 3)


