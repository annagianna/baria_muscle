# Differential abundance analysis with LinDA
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(MicrobiomeStat)
library(phyloseq)
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
model_color_manual <- scale_color_manual(values = c(renoir_15[6], renoir_15[14]))

# Data
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
baria_muscle_long <- readRDS("data/processed_data/BARIA_muscle_long.RDS")

### Prepare data for LinDA ###
# Filter out poorly annotated ("GGB"-containing) taxa
baria_mb_species <- prune_taxa(
  str_detect(rownames(otu_table(baria_mb)), "GGB\\d+", negate = TRUE),
  baria_mb
)

# Extract abundance matrix & metadata
mb_otu <- as(otu_table(baria_mb_species), "matrix") # taxa are rows
mb_sample_data <- as(sample_data(baria_mb_species), "data.frame")

# Join metadata from long dataset
linda_meta_long <- mb_sample_data |>
  filter(visit %in% c("v0", "v4", "v5")) |>
  rownames_to_column("Sample") |>
  left_join(
    baria_muscle_long |>
      filter(visit %in% c("v0", "v4", "v5")) |> 
      select(id, visit, age_v0, ffmi, fmi, perc_change_weight_kg, t2d_v0, t2d_labs, dm_meds),
    by = c("visit", "id")
  ) |> 
  filter(!is.na(ffmi), !is.na(perc_change_weight_kg)) |> 
  group_by(id) |> 
  mutate(
    ffmi_between = mean(ffmi, na.rm = TRUE),
    ffmi_within = ffmi - ffmi_between,
    fmi_between = mean(fmi, na.rm = TRUE),
    fmi_within = fmi - fmi_between
  ) |> 
  ungroup() |> 
  left_join(
    baria_muscle_long |> 
      filter(visit == "v0") |> 
      distinct(id, .keep_all = TRUE) |> 
      select(id, dm_meds_v0 = dm_meds),
    by = "id"
  ) |> 
  mutate(
    visit = factor(visit, levels = c("v0", "v4", "v5")),
    id = factor(id)
  )

# Match abundance matrix to long metadata
mb_otu_long <- mb_otu[ , linda_meta_long$Sample, drop = FALSE]

# Compute prevalence and mean abundance
species_prevalence <- rowMeans(mb_otu_long > 0)
species_abundance <- rowMeans(mb_otu_long)

## Filter
# Keep species detected in at least 20% of samples, with mean relative abundance >= 0.01%
species_keep <- tibble(
  species = names(species_prevalence),
  prevalence = species_prevalence,
  mean_abundance = species_abundance
) |>
  filter(
    prevalence >= 0.20,
    mean_abundance >= 0.01
  ) |>
  pull(species)

length(species_keep)

# Create species labels
species_labels <- tibble(species = species_keep) |>
  mutate(
    species_label = str_extract(species, "(?<=s__)[^|]+"),
    species_label = str_replace_all(species_label, "_", " ")
  )

# Filter abundance matrix
mb_otu_keep <- mb_otu_long[species_keep, , drop = FALSE]

#### LinDA ####
# Define formulas
linda_formulas <- list(
  model1 = "~ ffmi_within + ffmi_between + visit + age_v0 + sex + perc_change_weight_kg + (1 | id)",
  model2 = "~ ffmi_within + ffmi_between + visit + age_v0 + sex + perc_change_weight_kg + t2d_v0 + dm_meds_v0 + (1 | id)"
)

# LinDA function
run_linda_ffmi <- function(formula, model_name) {

  linda <- MicrobiomeStat::linda(
    feature.dat = mb_otu_keep, meta.dat = linda_meta_long, formula = formula, feature.dat.type = "proportion",
    prev.filter = 0, mean.abund.filter = 0, p.adj.method = "BH", alpha = 0.05
  )

  linda$output$ffmi_within |>
    rownames_to_column("species") |>
    left_join(species_labels, by = "species") |>
    mutate(
      model = model_name,
      signif = padj < 0.05,
      ci_lower = log2FoldChange - 1.96 * lfcSE,
      ci_upper = log2FoldChange + 1.96 * lfcSE
    )
  
}

# Run LinDA function for each formula
linda_ffmi_results <- imap_dfr(linda_formulas, ~ run_linda_ffmi(formula = .x, model_name = .y))

# Extract significant species with clean labels
species_linda_signif <- linda_ffmi_results |>
  filter(signif) |>
  distinct(species, species_label)

### Plot ###
linda_plot_data <- linda_ffmi_results |>
  filter(species %in% species_linda_signif$species) |> 
  mutate(
    model = factor(model, levels = c("model1", "model2"), labels = c("Age + sex + % weight change from baseline", "+ baseline T2D + diabetes medication"))
  )

## Save signif results
# Vector of significant species in at least one model
saveRDS(species_linda_signif,"data/processed_data/LinDA_significant_species.RDS")

# Save signif results table
species_linda_signif_summary <- linda_ffmi_results |>
  filter(species %in% species_linda_signif$species) |>
  select(model, species, species_label, log2FoldChange, lfcSE, ci_lower, ci_upper, pvalue, padj, signif) |>
  arrange(species_label, model)
write_csv(species_linda_signif_summary, "results/tables/LinDA_significant_species.csv")

### Forest plot of signif species ###
linda_forest_plot <- linda_plot_data |> 
  ggplot(aes(x = log2FoldChange, y = species_label, color = model, group = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = signif), size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") + 
  labs(x = expression("Change in log"[2]*" abundance per 1 kg/m"^2*" higher FFMI"), y = NULL, color = "Model") +
  guides(color = guide_legend(ncol = 1)) +
  model_color_manual +
  theme_minimal_custom() +
  theme(axis.text.y = element_text(face = "italic"))
ggsave("graphs/LinDA/LinDA_ffmi_forest_plot.pdf", linda_forest_plot, width = 12, height = 8)

