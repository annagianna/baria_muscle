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

# Data
baria_muscle_long <- readRDS("data/processed_data/260815_BARIA_muscle_long.RDS")
baria_mb <- readRDS("data/processed_data/260815_BARIA_mb_clean.RDS")
lm_ffmi_signif_v0 <- readRDS("tables/260811_model_ffmi_v0_1_signif.RDS")

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

# Combine results for both intervals
lmm_ffmi_results <- bind_rows(lmm_ffmi_v4_results, lmm_ffmi_v5_results)

lmm_ffmi_results |> 
  count(visit, signif)

lmm_ffmi_results |>
  filter(signif == "significant") |>
  arrange(visit, desc(abs(estimate))) |>
  select(
    visit,
    species,
    species_label,
    estimate,
    conf.low,
    conf.high,
    p.value,
    p_fdr
  )  

# Unique species labels for significant species in each model
  lmm_ffmi_species_labels <- species_labels |> 
    filter(species %in% lmm_ffmi_signif_species) |> 
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

# Combine v4 and v5 results
lmm1_ffmi_results <- bind_rows(lmm1_ffmi_v4_results, lmm1_ffmi_v5_results)
  
# Significant species
lmm_ffmi_signif <- lmm_ffmi_results |>
  filter(signif == "significant")

# Species significant at both v4 and v5
lmm_ffmi_overlap <- lmm_ffmi_signif |>
  group_by(species) |>
  filter(n_distinct(visit) == 2) |>
  ungroup()

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






#### Species overlapping with baseline associations ####
lm_lmm_ffmi_overlap <- lmm1_ffmi_overlap |> 
  inner_join(
    lm_ffmi_signif_v0 |> 
      select(-species_label), 
    by = "species") |> 
  rename(
    estimate_v0 = estimate,
    p_fdr_v0 = p_fdr
  ) |> 
  select(species, species_label_unique, contains("p_fdr_"), contains("estimate")) |> 
  relocate(contains("v0"), .before = contains("v4")) |> 
  relocate(contains("v4"), .before = contains("v5"))

lm_lmm_ffmi_overlap_table <- lm_lmm_ffmi_overlap |>
  rename(
    Species = species_label_unique,
    `Baseline β` = estimate_v0,
    `Baseline FDR` = p_fdr_v0,
    `1-year β` = estimate_v4,
    `1-year FDR` = p_fdr_v4,
    `2-year β` = estimate_v5,
    `2-year FDR` = p_fdr_v5
  ) |>
  select(-species)

write.csv(
  lm_lmm_ffmi_overlap_table,
  "tables/lm_lmm_ffmi_overlap_table.csv",
  row.names = FALSE
)

