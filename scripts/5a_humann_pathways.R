# Functional pathways (HUMAnN)
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
baria_humann <- readRDS("data/raw_data/BARIA.humann4.profiles.2026.581.910.RDS")
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
linda_species_signif <- readRDS("data/processed_data/LinDA_significant_species.RDS")

# HUMAnN cleaning
humann_clean <- baria_humann |>
  rownames_to_column(var = "pathway") |>
  mutate(
    pathway_id = str_extract(pathway, "^[^:]+"),
    pathway_name = str_remove(pathway, "^[^:]+:\\s*"),
    pathway_name = str_replace(pathway_name, "^.", toupper),
    pathway_type = case_when(pathway %in% c("UNMAPPED", "UNINTEGRATED") ~ "non-pathway", TRUE ~ "pathway"), 
    across(where(is.numeric), ~ replace_na(.x, 0))
  )

# Long dataset
humann_long <- humann_clean |>
  filter(pathway_type == "pathway") |>
  pivot_longer(
    cols = -c(pathway, pathway_id, pathway_name, pathway_type),
    names_to = "Sample_humann",
    values_to = "pathway_abundance"
  ) |>
  mutate(
    Sample = Sample_humann |>
      str_remove("_Abundance$") |>
      str_replace("\\.V\\.1\\.", ".V-1.") # HUMAnN data baseline V1, shotgun V-1
  ) |>
  inner_join(
    as(sample_data(baria_mb), "data.frame") |>
      rownames_to_column("Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  ) |>
  filter(visit %in% c("v0", "v4", "v5")) |>
  select(Sample, Sample_humann, id, visit, pathway, pathway_id, pathway_name, pathway_abundance)

# Extract significant species from LinDA
linda_species_long <- as(otu_table(baria_mb), "matrix") |>
  as.data.frame() |>
  rownames_to_column(var = "species") |>
  filter(species %in% linda_species_signif$species) |>
  pivot_longer(cols = -species, names_to = "Sample", values_to = "species_abundance") |>
  left_join(linda_species_signif, by = "species") |>
  left_join(
    as(sample_data(baria_mb), "data.frame") |>
      rownames_to_column("Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  )

# Filter HUMAnN pathways based on baseline prevalence (>5 CPM in at least 50% of samples)
humann_keep_v0 <- humann_long |> 
  filter(visit == "v0") |> 
  group_by(pathway_id, pathway_name) |> 
  summarize(prevalence = mean(pathway_abundance > 5), .groups = "drop") |> 
  filter(prevalence >= 0.50)

# Keep filtered baseline pathways across all visits (v0, v4, v5) & join significant LinDA species
humann_linda <- humann_long |>
  filter(pathway_id %in% humann_keep_v0$pathway_id) |> 
  inner_join(linda_species_long, by = c("Sample", "id", "visit"), relationship = "many-to-many")

### Correlations ###
# Species-pathway Spearman correlations per visit
humann_linda_cor <- humann_linda |>
  group_by(visit, species, species_label, pathway_id, pathway_name) |>
  summarise(
    n = sum(complete.cases(species_abundance, pathway_abundance)),
    test = list(cor.test(species_abundance, pathway_abundance, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) |>
  mutate(
    rho = map_dbl(test, ~ unname(.x$estimate)),
    p = map_dbl(test, ~ .x$p.value)
  ) |>
  select(-test) |>
  group_by(visit, species) |>
  mutate(padj = p.adjust(p, method = "BH")) |>
  ungroup()

# Strongest significant baseline pathway associations per species
humann_linda_cor_top10 <- humann_linda_cor |>
  filter(visit == "v0",padj < 0.05) |>
  mutate(direction = if_else(rho > 0, "positive", "negative")) |>
  group_by(species_label, direction) |> 
  arrange(desc(abs(rho)), .by_group = TRUE) |> 
  slice_head(n = 10) |> 
  mutate(rank_direction = row_number()) |> 
  ungroup() |> 
  group_by(species_label) |> 
  mutate(# Select up to 10 strongest signif baseline pathway associations per species (up to 5 from each direction where available)
    priority = if_else(rank_direction <= 5, 1L, 2L)
  ) |>
  arrange(priority, desc(abs(rho)), .by_group = TRUE) |> 
  slice_head(n = 10) |> 
  ungroup() |> 
  select(species_label, pathway_id, pathway_name, rho, padj) |> 
  print(n = 40)

# Pull union of pathways
pathways_top <- humann_linda_cor_top10 |>
  distinct(pathway_id) |>
  pull(pathway_id)

### Plot ###
# Plot data
humann_linda_plot_data <- humann_linda_cor |>
  filter(visit == "v0", pathway_id %in% pathways_top) |>
  mutate(
    signif = padj < 0.05,
    species_label = factor(species_label, levels = linda_species_signif$species_label
    )
  )

rho_max <- max(abs(humann_linda_plot_data$rho), na.rm = TRUE)
fill_cor_manual <- scale_fill_gradient2(low = renoir_15[14], mid = "white", high = renoir_15[6], midpoint = 0, limits = c(-rho_max, rho_max), name = "Spearman\nrho")

## Heatmap
# Species x axis
humann_linda_heatmap_x <- humann_linda_plot_data |> 
  ggplot(aes(x = species_label, y = pathway_name, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(
    data = humann_linda_plot_data |> filter(signif),
    shape = 8,
    size = 2.3
  ) +
  fill_cor_manual +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal_custom(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave("results/figures/HUMAnN_LinDA_pathway_heatmap_x.pdf", humann_linda_heatmap_x, width = 8, height = 10)

# Species y axis (for pptx)
humann_linda_heatmap_y <- humann_linda_plot_data |> 
  ggplot(aes(x = pathway_name, y = species_label, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(
    data = humann_linda_plot_data |> filter(signif),
    shape = 8,
    size = 2.3
  ) +
  fill_cor_manual +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal_custom(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(face = "italic"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave("results/figures/HUMAnN_LinDA_pathway_heatmap_y.pdf", humann_linda_heatmap_y, width = 12, height = 5)

## Bubble plot
humann_linda_bubble <- humann_linda_plot_data |> 
  ggplot(aes(x = species_label, y = pathway_name)) +
  geom_point(aes(size = abs(rho), fill = rho), shape = 21, colour = "grey30", stroke = 0.3) +
  geom_text(aes(label = if_else(signif, "*", "")), size = 3) +
  scale_size_continuous(range = c(1, 7), guide = "none") +
  fill_cor_manual +
  labs(x = NULL, y = NULL) +
  theme_minimal_custom(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave("results/figures/HUMAnN_LinDA_pathway_bubble.pdf", humann_linda_bubble, width = 8, height = 10)
