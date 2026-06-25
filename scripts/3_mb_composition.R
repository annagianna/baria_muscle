# Baria: Associations between muscle mass loss and gut microbiota (Baseline)
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(MetBrewer)
library(ggthemes)
library(patchwork)

# Theme
manet_5 <- met.brewer("Manet", n = 5)
manet_20 <- met.brewer("Manet", n = 20)

theme_composition <- function(base_size = 14, base_family = "sans") {
  
  (theme_foundation(base_size=base_size, base_family=base_family) + 
    theme(
      plot.title = element_text(face = "bold",size = rel(1.0), hjust = 0.5), 
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle = 90,vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      #axis.text.x =  element_text(angle = 45, hjust = 1),
      axis.text = element_text(), 
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "right",
      legend.key.size= unit(0.4, "cm"),
      legend.spacing  = unit(0, "cm"),
      legend.text = element_text(size = rel(0.7)),
      plot.margin = unit(c(10,5,5,5),"mm"),
      strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
  ))
}

theme_Publication <- function(base_size=14, base_family="sans") {
  
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
      strip.background=element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
} 

# Data
baria_muscle <- read_rds("data/20260624_BARIA_muscle_clinical.RDS") # metadata/clinical data CHECK FOR MOST RECENT VERSION
mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")

# QC
# Check that the raw metaphlan rel. abundances sum to 100% per sample
melted_mb |> 
  group_by(Sample) |> 
  summarize(total_abundance = sum(Abundance)) |> 
  summarize(
    min = min(total_abundance),
    mean = mean(total_abundance),
    max = max(total_abundance)
  ) # returns 100 for all

# Identify samples with multiple runs (per sample per visit)
run2_samples <- melted_mb |> 
  distinct(Sample, Time_Point, Subject_ID) |> 
  count(Subject_ID, Time_Point) |> 
  filter(n > 1)

# Merge with metadata
mb <- melted_mb |> 
  filter(!is.na(Extra_data) | Extra_data == "rep1") |> 
  select(-c(Study, Sample_Type, Data_Type)) |> 
  mutate(
    visit = str_extract(Time_Point, "\\d"),
    visit = if_else(visit == "1", "0", visit),
    id = Subject_ID
  ) |> 
  tibble::rownames_to_column(var = "otu") |> 
  select(-Subject_ID, -Time_Point, -OTU) |> 
  left_join(baria_muscle, by = join_by(id)) |> 
  group_by(id, visit, Species) |>
  ungroup() |> 
  relocate(visit, .before = otu) |> 
  relocate(id, .before = visit)

# QC sums
mb |> 
  group_by(Sample) |> 
  summarise(total_abundance = sum(Abundance, na.rm = TRUE)) |> 
  summarise(
    min = min(total_abundance),
    mean = mean(total_abundance),
    max = max(total_abundance)
  ) # returns 100 for all

#### Composition Plots ####
#### Statified by FFMI status at baseline ####
### Species level ###
# Summarize per group and identify top 20 species (baseline = V1) 
top20_species <- mb |> 
  # filter(visit == 0) |> # not filtered for longitudinal approach
  group_by(Sample, Species) |> 
  dplyr::summarize(Abundance = sum(Abundance), .groups = "drop") |> 
  group_by(Species) |> 
  dplyr::summarize(Abundance = mean(Abundance)) |> 
  arrange(-Abundance) |> 
  select(Species) |> 
  head(20)

# assign fixed color to each one of the top20 species
top20_species_vector <- top20_species$Species # top 20 species names as vector
set.seed(13)
top20_species_colours <- c("Other species" = "grey63", setNames(sample(manet_20), top20_species_vector))

### Baseline composition ###
species_ffmi_v0 <- mb |> 
  filter(
    visit == 0, # baseline composition (in clinical data v0 for shotgun data v1!)
    !is.na(low_ffmi_v0)
  ) |> 
  mutate(
    Species2 = if_else(
        Species %in% top20_species$Species,
        Species,
        "Other species" # collapse other species
        ),
    Species2 = as.factor(Species2)
  ) |> 
  group_by(Species2, Sample, low_ffmi_v0) |> # species abundance per sample
  dplyr::summarize(Abundance = sum(Abundance)) |> 
  group_by(Species2, low_ffmi_v0) |> # species per muscle mass group
  dplyr::summarize(Abundance = mean(Abundance)) |> # avg abundance per species per FFMI group
  ungroup() |>
  mutate(
    Species2 = fct_reorder(Species2, Abundance),
    Species2 = fct_relevel(Species2, "Other species", after = 0L) # move other species to the front
  )

species_ffmi_v0 |> # check
  group_by(low_ffmi_v0) |> 
  summarize(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
#species_comp_asm1y_v0 <- mb_species_asm1y_v0 |> 
  #mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  #ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Species2)) +
  #geom_bar(stat = "identity", color = "black", width = 0.9) +
  #scale_fill_manual(values = top20_species_colours) +
  #guides(fill = guide_legend(ncol = 1)) +
  #labs(y = "Relative abundance (%)", x = "%ASM change at 1y", title = "Species", fill = "") +
  #scale_y_continuous(expand = c(0, 0)) +
  #theme_composition() +
  #theme(
   # axis.title.x = element_blank(),
   # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  #)
#ggsave("graphs/composition/species_comp_asm1y_v0.pdf", width = 12, height = 10)