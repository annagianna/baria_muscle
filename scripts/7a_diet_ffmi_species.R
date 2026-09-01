# Dietary intake and FFMI-associated microbiome features
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(grid)
library(phyloseq)
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
      panel.grid.major = element_line(
        colour = "#dddddd",
        linewidth = 0.4,
        linetype = "22"
      ),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "bottom",
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )
}

renoir_15 <- met.brewer("Renoir", n = 15)

# Data
baria_diet <- readRDS("data/raw_data/251002_BARIA_macronutrients.RDS")
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
baria_humann <- readRDS("data/raw_data/BARIA.humann4.profiles.2026.581.910.RDS")
linda_species_signif <- readRDS("data/processed_data/LinDA_significant_species.RDS")

# Clean dietary data
baria_diet_clean <- baria_diet |> 
  ungroup() |> 
  rename(id = ID, visit = Visit) |> 
  mutate(
    id = str_remove(id, "^BARIA_"),
    visit = recode(as.character(visit), "V1" = "v0")
  )

# Join baseline dietary data with sample data
diet_mb <- as(sample_data(baria_mb), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(Sample, id, visit) |>
  filter(visit == "v0") |>
  inner_join(baria_diet_clean, by = c("id", "visit"))

# Dietary vars of interest
diet_vars <- c("TotalCal_kcal", "Carbs_g", "Protein_g", "Fat_g", "SatFat_g", "Fibers_g", "Alcohol_g", "Water_ml")

# Extract significant LinDA species and join dietary data
diet_species_long <- as(otu_table(baria_mb), "matrix") |>
  as.data.frame() |>
  rownames_to_column(var = "species") |>
  filter(species %in% linda_species_signif$species) |>
  pivot_longer(cols = -species, names_to = "Sample", values_to = "species_abundance") |>
  left_join(linda_species_signif, by = "species") |>
  inner_join(
    diet_mb |>
      select(Sample, id, visit, all_of(diet_vars)),
    by = "Sample"
  ) |> 
   pivot_longer(
    cols = all_of(diet_vars),
    names_to = "diet_feature",
    values_to = "diet_value"
  )

### Correlations ###
# Spearman correlations between dietary intake and significant LinDA species
diet_species_cor <- diet_species_long |>
  group_by(species, species_label, diet_feature) |>
  summarise(
    n = sum(complete.cases(species_abundance, diet_value)),
    test = list(
      cor.test(species_abundance, diet_value, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) |>
  mutate(
    rho = map_dbl(test, ~ unname(.x$estimate)),
    p = map_dbl(test, ~ .x$p.value)
  ) |>
  select(-test) |>
  ungroup() |>
  mutate(padj = p.adjust(p, method = "BH"))

# Plot
# Plot data
diet_species_plot_data <- diet_species_cor |>
  mutate(
    signif = padj < 0.05,
    species_label = factor(species_label, levels = linda_species_signif$species_label),
    diet_feature = dplyr::recode(
      diet_feature, TotalCal_kcal = "Energy (kcal)", Carbs_g = "Carbohydrates (g)", Protein_g = "Protein (g)", Fat_g = "Fat (g)",
      SatFat_g = "Saturated fat (g)", Fibers_g = "Fiber (g)", Alcohol_g = "Alcohol (g)", Water_ml = "Water (mL)"
    ),
    diet_feature = factor(
      diet_feature,
      levels = c("Energy (kcal)", "Carbohydrates (g)", "Protein (g)", "Fat (g)", "Saturated fat (g)", "Fiber (g)", "Alcohol (g)", "Water (mL)")
    )
  )

# Symmetric colour scale
rho_max_diet <- max(abs(diet_species_plot_data$rho), na.rm = TRUE)
fill_cor_diet <- scale_fill_gradient2(
  low = renoir_15[14],
  mid = "white",
  high = renoir_15[6],
  midpoint = 0,
  limits = c(-rho_max_diet, rho_max_diet),
  name = "Spearman\nrho"
)

# Heatmap
diet_species_heatmap <- ggplot(
  diet_species_plot_data,
  aes(x = diet_feature, y = species_label, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = diet_species_plot_data |> filter(signif), shape = 8, size = 2.3) +
  fill_cor_diet +
  labs(x = NULL, y = NULL) +
  theme_minimal_custom(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave("results/figures/Diet_LinDA_species_heatmap.pdf", diet_species_heatmap, width = 8, height = 4)
