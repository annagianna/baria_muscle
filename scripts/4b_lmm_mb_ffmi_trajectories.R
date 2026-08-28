# Associations between baseline species abundance & FFMI trajectories
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(broom.mixed)
library(ggrepel)
library(ggthemes)
library(ggpubr)
library(grid)
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
color_manual_abundance <- scale_color_manual(
  values = c("Low" = renoir_15[6], "High" = renoir_15[14]),
  labels = c("Low" = "Low", "High" = "High")
)
fill_manual_abundance <- scale_fill_manual(
  values = c("Low" = renoir_15[6], "High" = renoir_15[14]),
  labels = c("Low" = "Low", "High" = "High")
)

# Data
baria_muscle_long <- readRDS("data/processed_data/BARIA_muscle_long.RDS")
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")

# Filter out poorly annotated ("GGB"-containing) taxa
baria_mb_species <- prune_taxa(str_detect(rownames(otu_table(baria_mb)), "GGB\\d+", negate = TRUE), baria_mb)

# Extract abundance matrix
matrix_mb <- as(otu_table(baria_mb_species), "matrix") |> 
  t()

# Create baseline mb data
mb_v0 <- matrix_mb |>
  as.data.frame() |>
  rownames_to_column(var = "Sample") |>
  left_join(
    as(sample_data(baria_mb_species), "data.frame") |>
      rownames_to_column(var = "Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  ) |>
  filter(visit == "v0") |>
  relocate(id, visit, .before = everything())

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
    prevalence_v0 >= 0.50,
    mean_abundance_v0 >= 0.1
  ) |>
  pull(species)

length(species_v0_keep)

# Prepare subset needed for lmm
baria_muscle_lmm <- baria_muscle_long |>
  filter(visit %in% c("v0", "v4", "v5"), !is.na(ffmi)) |>
  mutate(
    visit = factor(visit, levels = c("v0", "v4", "v5")),
    age_centered_v0 = age_v0 - mean(age_v0, na.rm = TRUE)
  )

# Model data (long format) (merge long metadata with mb + prevalence + abundance matrices)
lmm_data_ffmi <- mb_v0 |>
  pivot_longer(
    cols = all_of(species_v0_keep),
    names_to = "species",
    values_to = "baseline_abundance"
  ) |> 
  select(-visit, -Sample) |> 
  inner_join(
    baria_muscle_lmm |> 
      select(id, visit, ffmi, perc_change_weight_kg, perc_change_ffmi, delta_ffmi, age_centered_v0, sex),
    by = "id"
  ) |> 
  group_by(species) |> 
  mutate(
    # Log10-transform abundance and calculate a species-specific pseudocount
    min_baseline_abundance = min(baseline_abundance[baseline_abundance > 0], na.rm = TRUE),
    log10_baseline_abundance = log10(baseline_abundance + min_baseline_abundance / 2)
  ) |> 
  ungroup()

# Create consistent species labels for use across all models
species_labels <- tibble(species = sort(unique(lmm_data_ffmi$species))) |>
  mutate(
    species_label = str_extract(species, "s__[^|]+"),
    species_label = str_remove(species_label, "^s__"),
    #species_label = str_remove(species_label, "_SGB\\d+$"),
    species_label = str_replace_all(species_label, "_", " "),

    sgb = str_extract(species, "t__SGB\\d+"),
    sgb = str_remove(sgb, "^t__")
  )

#### Models ####
##### Model function ####
run_lmm_ffmi <- function(lmm_data, follow_up, species_labels) {

  lmm_data |>
    filter(visit %in% c("v0", follow_up)) |>
    droplevels() |>
    group_by(species) |>
    nest() |>
    mutate(
      model = map(
        data, 
        ~ lmerTest::lmer(ffmi ~ log10_baseline_abundance * visit + age_centered_v0 + sex + perc_change_weight_kg + (1 | id), 
        data = .x, REML = FALSE
        )),
      results = map(model, ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE))
    ) |>
    select(species, results) |>
    unnest(results) |>
    filter(str_detect(term, "log10_baseline_abundance:visit")) |>
    ungroup() |> # Remove species grouping before FDR correction across all species
    mutate(
      p_fdr = p.adjust(p.value, method = "BH"),
      signif = if_else(p_fdr < 0.05, "significant", "not significant"),
      visit = follow_up
    ) |>
    left_join(species_labels, by = "species")

}

# Run the function for each follow-up visit (v4 & v5)
lmm_ffmi_v4_results <- run_lmm_ffmi(lmm_data = lmm_data_ffmi, follow_up = "v4", species_labels = species_labels)
lmm_ffmi_v5_results <- run_lmm_ffmi(lmm_data = lmm_data_ffmi, follow_up = "v5", species_labels = species_labels)

# Combine results for both intervals (v4 & v5)
lmm_ffmi_results <- bind_rows(lmm_ffmi_v4_results, lmm_ffmi_v5_results)

# Species significantly associated with FFMI trajectory at 1 year
lmm_ffmi_v4_signif <- lmm_ffmi_v4_results |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

# Species significantly associated with FFMI trajectory at 2 years
lmm_ffmi_v5_signif <- lmm_ffmi_v5_results |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

# Combine significant species across follow-up intervals
lmm_ffmi_signif <- bind_rows(lmm_ffmi_v4_signif, lmm_ffmi_v5_signif)
  
#### Plots ####
#### Plot observed FFMI trajectories ####
plot_ffmi_observed <- baria_muscle_long |> 
  filter(visit %in% c("v0", "v4", "v5")) |>
  select(id, visit, ffmi) |>
  mutate(visit = factor(visit, levels = c("v0", "v4", "v5"), labels = c("Baseline", "1 year", "2 years"))) |> 
  ggplot(aes(x = visit, y = ffmi, group = id)) +

  # Individual trajectories
  geom_line(color = "grey70", alpha = 0.2, linewidth = 0.6) +
  geom_point(color = "grey70", alpha = 0.2, size = 0.6) +

  # Group mean trajectories
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1.2) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", color = "black", size = 1.5) +
  labs(
    title = "Observed FFMI trajectories",
    x = "Baseline log10-abundance",
    y = "FFMI (kg/m²)"
  ) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.01))) +
  theme_minimal_custom()
ggsave(plot = plot_ffmi_observed, filename = "graphs/lmm_species_ffmi/plot_ffmi_observed.pdf", width = 10, height = 8)

# Prepare plotting data for significant v4 species
lmm_v4_plot_data <- lmm_data_ffmi |>
  filter(
    species %in% lmm_ffmi_v4_signif$species,
    visit %in% c("v0", "v4")
  ) |>
  group_by(species) |>
  mutate(
    abundance_tertile_low = quantile(baseline_abundance[visit == "v0"], probs = 1/3, na.rm = TRUE),
    abundance_tertile_high = quantile(baseline_abundance[visit == "v0"], probs = 2/3, na.rm = TRUE),
    abundance_group = case_when(
      baseline_abundance <= abundance_tertile_low ~ "Low",
      baseline_abundance >= abundance_tertile_high ~ "High"
    )
  ) |>
  ungroup() |>
  filter(!is.na(abundance_group)) |>
  left_join(
    lmm_ffmi_v4_signif |>
      select(species, species_label),
    by = "species"
  )

# Plot %FFMI change by baseline abundance group
plot_perc_ffmi_change_v4 <- lmm_v4_plot_data |>
  filter(visit == "v4") |>
  ggplot(aes(x = abundance_group, y = perc_change_ffmi, color = abundance_group, fill = abundance_group)) +
  stat_boxplot(geom = "errorbar", width = 0.15) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 1) +
  labs(
    title = unique(lmm_v4_plot_data$species_label),
    x = "Baseline log10-abundance",
    y = "FFMI change from baseline (%)",
    color = "Baseline log10-abundance",
    fill = "Baseline log10-abundance"
  ) +
  color_manual_abundance +
  fill_manual_abundance +
  theme_minimal_custom() +
  theme(plot.title = element_text(face = "italic"), legend.position = "none")
ggsave(plot = plot_perc_ffmi_change_v4, filename = "graphs/lmm_species_ffmi/plot_perc_ffmi_change_v4_signif.pdf", width = 10, height = 7)
