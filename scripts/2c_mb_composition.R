# Baseline microbial composition
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(grid)
library(MetBrewer)

# Theme
theme_minimal_composition <- function(base_size = 14, base_family = "sans") {

  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      panel.grid.major.y = element_line(colour = "#dddddd", linewidth = 0.4, linetype = "22"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = rel(0.7)),
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )
}
renoir_20 <- met.brewer("Renoir", n = 20)

# Data
baria_mb <- readRDS("data/processed_data/260818_BARIA_mb_clean.RDS") # already contains necessary metadata for the grouping

# Melt into df
melted_mb <- psmelt(baria_mb)

#### Composition Plots ####
#### Statified by FFMI status at baseline ####
### Species level ###
# Summarize per group and identify top 20 species at baseline
top20_species_v0 <- melted_mb |> 
  filter(visit == "v0") |>
  group_by(Sample, Species) |> 
  summarize(Abundance = sum(Abundance), .groups = "drop") |> 
  group_by(Species) |> 
  summarize(Abundance = mean(Abundance)) |> 
  arrange(-Abundance) |> 
  slice_head(n = 20)

### Baseline composition ###
species_ffmi_v0 <- melted_mb |>
  filter(visit == "v0", !is.na(low_ffmi_v0)) |>
  mutate(
    Species2 = if_else(
      Species %in% top20_species_v0$Species,
      Species,
      "Other species" # Collapse other species
    )
  ) |>
  group_by(Species2, Sample, low_ffmi_v0) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
  group_by(Species2, low_ffmi_v0) |> # species per muscle mass group
  summarize(Abundance = mean(Abundance), .groups = "drop") # avg abundance per species per FFMI group

# check sum
species_ffmi_v0 |> # check
  group_by(low_ffmi_v0) |> 
  summarize(sum_Abundance = sum(Abundance)) # adds up to 100

# Create order based on the abundance in the first (low FFMI) group
species_order_low_ffmi_v0 <- species_ffmi_v0 |> 
  filter(
    low_ffmi_v0 == "yes",
    Species2 != "Other species"
  ) |> 
  arrange(-Abundance) |> 
  pull(Species2)

species_ffmi_v0 <- species_ffmi_v0 |>
  mutate(Species2 = factor(Species2, levels = c(species_order_low_ffmi_v0, "Other species")))
  
# Assign colors according to abundance order in the low FFMI group
set.seed(13)
top20_species_v0_low_colours <- c("Other species" = "grey75", setNames(sample(renoir_20), species_order_low_ffmi_v0))

# Composition plot
species_comp_ffmi_v0 <- species_ffmi_v0 |> 
  mutate(
    low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes", after = 0L) # low baseline FFMI first
  ) |> 
  ggplot(aes(x = low_ffmi_v0, y = Abundance, fill = Species2)) +
  geom_bar(
    stat = "identity",
    color = "black",
    width = 0.65,
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(values = top20_species_v0_low_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "Relative abundance (%)", title = "Species", fill = "") +
  scale_x_discrete(labels = c("yes" = "Low FFMI", "no" = "Moderate/high FFMI")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal_composition() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = rel(0.9))
  )
ggsave(plot = species_comp_ffmi_v0, "graphs/composition/species_comp_ffmi_v0.pdf", width = 14, height = 8)

