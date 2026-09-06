# Baseline microbial composition stratified by 1-year %FFMI change group
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(grid)
library(MetBrewer)

# Theme
source("scripts/assets/functions.R")
renoir_20 <- met.brewer("Renoir", n = 20)

# Data
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")

# Output folder
dir.create("results/graphs/composition", recursive = TRUE, showWarnings = FALSE)

# Melt phyloseq object into dataframe
melted_mb <- psmelt(baria_mb)

#### Composition plots stratified by 1-year % FMI change group ####
# Keep baseline samples with known 1-year FFMI
melted_mb_v0 <- melted_mb |>
  filter(
    visit == "v0",
    !is.na(perc_change_ffmi_v4_group)
  ) |>
  mutate(perc_change_ffmi_v4_group = factor(perc_change_ffmi_v4_group, levels = c("high", "modest/low")))

#### Species-level composition ####
# Identify top 20 species at baseline
top20_species_v0 <- melted_mb_v0 |>
  group_by(Sample, Tax) |>
  summarize(Abundance = sum(Abundance),.groups = "drop") |>
  group_by(Tax) |>
  summarize(Abundance = mean(Abundance), .groups = "drop") |>
  arrange(desc(Abundance)) |>
  slice_head(n = 20)

# Collapse all other species into "Other species"
species_ffmi_change_v0 <- melted_mb_v0 |>
  mutate(Species2 = if_else(Tax %in% top20_species_v0$Tax, Tax,"Other species")) |>
  group_by(Species2, Sample, perc_change_ffmi_v4_group) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
  group_by(Species2, perc_change_ffmi_v4_group) |>
  summarize(Abundance = mean(Abundance), .groups = "drop")

# Create species order based on abundance in high %FFMI-change group
species_order_high_ffmi_change_v0 <- species_ffmi_change_v0 |>
  filter(
    perc_change_ffmi_v4_group == "high",
    Species2 != "Other species"
  ) |>
  arrange(desc(Abundance)) |>
  pull(Species2)

# Apply species order
species_ffmi_change_v0 <- species_ffmi_change_v0 |>
  mutate(Species2 = factor(Species2, levels = c(species_order_high_ffmi_change_v0, "Other species")))

# Assign colours according to abundance order in high-change group
set.seed(13)
top20_species_v0_high_colours <- c("Other species" = "grey75", setNames(sample(renoir_20), species_order_high_ffmi_change_v0))

# Composition plot
species_comp_ffmi_change_v0 <- species_ffmi_change_v0 |>
  ggplot(aes(x = perc_change_ffmi_v4_group, y = Abundance, fill = Species2)) +
  geom_bar(stat = "identity", colour = "black", width = 0.65, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = top20_species_v0_high_colours) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
  scale_x_discrete(labels = c("high" = "High", "modest/low" = "Modest/low")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance (%)", title = "Species", fill = NULL) +
  theme_minimal_composition() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = rel(0.9))
  )
ggsave(filename = "results/graphs/composition/species_comp_by_ffmi_change_v4.pdf", plot = species_comp_ffmi_change_v0, width = 14, height = 8)
