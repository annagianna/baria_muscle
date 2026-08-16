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
color_manual_overlap_species <- scale_color_manual(values = c("1 year" = renoir_15[5], "2 years" = renoir_15[3]))
color_manual_abundance <- scale_color_manual(values = c("Low baseline abundance" = renoir_15[6], "High baseline abundance" = renoir_15[14]))
fill_manual_abundance <- scale_fill_manual(values = c("Low baseline abundance" = renoir_15[6], "High baseline abundance" = renoir_15[14]))

# Data
baria_muscle_long <- readRDS("data/processed_data/260816_BARIA_muscle_long.RDS")
baria_mb <- readRDS("data/processed_data/260816_BARIA_mb_clean.RDS")
lm_ffmi_signif_v0 <- readRDS("tables/260816_model_ffmi_v0_1_signif.RDS")

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
    prevalence_v0 >= 0.20,
    mean_abundance_v0 >= 0.01
  ) |>
  pull(species)

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
      select(id, visit, ffmi, perc_change_weight_kg, age_centered_v0, sex),
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

#### Model: Age, sex, %weight change ####
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

## Significant species
lmm_ffmi_signif <- lmm_ffmi_results |>
  filter(signif == "significant") |> 
  group_by(species_label) |>
  mutate(
    # add sgb identifier only for species appearing more than once in the significant results in both intervals
    species_label_unique = case_when(n_distinct(species) > 1 ~ paste(species_label, sgb), TRUE ~ species_label)
  ) |>
  ungroup()

# Species significant at both v4 and v5
lmm_ffmi_overlap <- lmm_ffmi_signif |>
  group_by(species) |>
  filter(n_distinct(visit) == 2) |>
  ungroup() |> 
  mutate(
    visit = factor(visit, levels = c("v4", "v5"), labels = c("1 year", "2 years")),
    signif = factor(signif, levels = c("significant", "not significant"))
  ) |> 
  group_by(species_label_unique) |> 
  mutate(estimate_v4 = estimate[visit == "1 year"]) |> 
  ungroup() |> 
  mutate(species_label_unique = fct_reorder(species_label_unique, estimate_v4))

# Species significant only at v4
lmm_ffmi_v4_only <- lmm_ffmi_signif |>
  group_by(species) |>
  filter(n_distinct(visit) == 1, visit == "v4") |>
  ungroup()

# Species significant only at v5
lmm_ffmi_v5_only <- lmm_ffmi_signif |>
  group_by(species) |>
  filter(n_distinct(visit) == 1, visit == "v5") |>
  ungroup()
  
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
    x = NULL,
    y = "FFMI (kg/m²)"
  ) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.01))) +
  theme_minimal_custom()
ggsave(plot = plot_ffmi_observed, filename = "graphs/lmm_species_ffmi/plot_ffmi_observed.pdf", width = 10, height = 8)

#### Forest plots ####
# Species significant in both v4 and v5
lmm_ffmi_overlap_forest <- lmm_ffmi_overlap |> 
  ggplot(aes(x = estimate, y = species_label_unique, group = visit, color = visit)) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.5, 
    position = position_dodge(width = 0.8, reverse = TRUE),
    show.legend = FALSE
  ) +
  labs(
    x = "Difference in FFMI change (kg/m²) per 1-unit increase in baseline log10 abundance",
    y = NULL,
    colour = NULL
  ) +
  geom_point(size = 2.5, position = position_dodge(width = 0.8, reverse = TRUE)) + 
  color_manual_overlap_species +
  theme_minimal_custom() +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "italic"))
ggsave(plot = lmm_ffmi_overlap_forest, filename = "graphs/lmm_species_ffmi/lmm_ffmi_overlap_forest.pdf", width = 13, height = 8)


#### Species overlapping with baseline associations ####
# Reshape longitudinal overlap results to wide format
lm_lmm_ffmi_overlap <- lmm_ffmi_overlap |> 
  mutate(visit = case_when(visit == "1 year" ~ "v4", visit == "2 years" ~ "v5")) |> 
  select(species, species_label_unique, visit, estimate, p_fdr) |> 
  pivot_wider(
    names_from = visit,
    values_from = c(estimate, p_fdr),
    names_glue = "{.value}_{visit}"
  ) |> 
  inner_join(
    lm_ffmi_signif_v0 |> 
      select(
        species,
        estimate_v0 = estimate,
        p_fdr_v0 = p_fdr),
    by = "species"
  ) |> 
  select(contains("species"), contains("p_fdr"), contains("estimate"))

lm_lmm_ffmi_overlap_table <- lm_lmm_ffmi_overlap |>
  relocate(contains("v0"), .before = contains("v4")) |> 
  relocate(contains("v4"), .before = contains("v5")) |> 
  rename(
    Species = species_label_unique,
    `Baseline β` = estimate_v0,
    `Baseline FDR p` = p_fdr_v0,
    `1-year β` = estimate_v4,
    `1-year FDR p` = p_fdr_v4,
    `2-year β` = estimate_v5,
    `2-year FDR p` = p_fdr_v5
  ) |>
  select(-species)
write.csv(lm_lmm_ffmi_overlap_table, "tables/lm_lmm_ffmi_overlap_table.csv", row.names = FALSE)

### Plot trajectories for overlap species between cross-sectional and trajecotry associations ###
lm_lmm_plot_data <- lmm_data_ffmi |> 
  filter(species %in% lm_lmm_ffmi_overlap$species) |> 
  group_by(species) |> 
  filter(any(visit %in% c("v4", "v5"))) |>
  mutate(
    abundance_tertile_low = quantile(baseline_abundance[visit == "v0"], probs = 1/3, na.rm = TRUE),
    abundance_tertile_high = quantile(baseline_abundance[visit == "v0"], probs = 2/3, na.rm = TRUE),
    abundance_group = case_when(
      baseline_abundance <= abundance_tertile_low ~ "Low baseline abundance",
      baseline_abundance >= abundance_tertile_high ~ "High baseline abundance"
    )
  ) |> 
  ungroup() |> 
  filter(!is.na(abundance_group)) |> 
  left_join(
    lm_lmm_ffmi_overlap |> 
      select(species, species_label_unique),
    by = "species"
  ) |> 
  mutate(
    visit = factor(visit, levels = c("v0", "v4", "v5"), labels = c("Baseline", "1 year", "2 years")),
    abundance_group = factor(abundance_group, levels = c("Low baseline abundance", "High baseline abundance"))
  )

# Plot function
plot_ffmi_trajectory <- function(species_name, plot_data) {

  species_data <- plot_data |> 
    filter(species == species_name)

  jitter_position <- position_jitter(width = 0.1, height = 0, seed = 2026)

  ggplot(species_data, aes(x = visit, y = ffmi, color = abundance_group, fill = abundance_group)) +

    # Individual trajectories
    geom_line(aes(group = id), color = "grey70", alpha = 0.4, linewidth = 0.4, position = jitter_position) +

    # FFMI distribution
    stat_boxplot(geom = "errorbar", width = 0.1) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.4) +

    # Individual observations
    geom_point(aes(color = abundance_group), width = 0.2, alpha = 0.8, size = 1, position = jitter_position) +
    
    facet_wrap(~ abundance_group) +
    
    labs(
      title = unique(species_data$species_label_unique),
      x = NULL, y = "FFMI (kg/m²)",
      color = "Baseline abundance",
      fill = "Baseline abundance"
    ) +
    color_manual_abundance +
    fill_manual_abundance +
    theme_minimal_custom() +
    theme(
      plot.title = element_text(face = "italic"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
  )

}

# Run function for each species
trajectory_plots <- map(lm_lmm_ffmi_overlap$species, ~ plot_ffmi_trajectory(species_name = .x, plot_data = lm_lmm_plot_data))

# Save plots
ggsave(plot = trajectory_plots[[1]], "graphs/lmm_species_ffmi/trajectory_species_1.pdf", width = 12, height = 7)
ggsave(plot = trajectory_plots[[2]], "graphs/lmm_species_ffmi/trajectory_species_1.pdf", width = 12, height = 7)
ggsave(plot = trajectory_plots[[3]], "graphs/lmm_species_ffmi/trajectory_species_1.pdf", width = 12, height = 7)

