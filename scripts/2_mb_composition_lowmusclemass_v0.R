# Baria: associations between low muscle mass and gut microbiota (Baseline)
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(stringr)
library(MetBrewer)
library(grid)
library(ggthemes)

# Theme
manet <- met.brewer("Manet", n = 20)
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
baria_muscle <- read_rds("data/251221_BARIA_muscle_clinical.RDS") # metadata/clinical data
baria_mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")
sample_sums(baria_mb) # adds up to 100

# convert sample_data of baria_mb to a data frame to merge with selected cols from baria_muscle
meta_baria_mb <- meta(baria_mb) |> 
  tibble::rownames_to_column(var = "rownames")
class(meta_baria_mb) # is df

# Add low muscle variables to sample_data() of phyloseq object
new_sample_data <- baria_muscle |> 
  select(id, sex, contains("low_")) |> # to be used for group comparisons (low muscle mass)
  mutate(id = as.character(id)) |> 
  inner_join(meta_baria_mb, join_by(id == Subject_ID)) |> # only keep ids with BIA and microbiome data
  mutate(
    visit = as.integer(str_extract(Time_Point, "\\d+")),
    Extra_data = na_if(Extra_data, "NA")
  ) |> 
  select(-(Study:Time_Point)) |> 
  group_by(id, visit) |> 
  slice(1) |> # keeps first occurence per group (some ids had second runs)
  relocate(rownames, .before = id) |> 
  relocate(visit, .before = id) |> 
  relocate(Extra_data, .after = id) |> 
  tibble::column_to_rownames(var = "rownames")

# update sample_data() within baria_mb
sample_data(baria_mb) <- sample_data(new_sample_data)

# melt phyloseq object into a df
baria_mb_df <- psmelt(baria_mb)
baria_mb_df |> # check
  group_by(Sample) |> 
  summarise(sum_abundance = sum(Abundance)) # adds up to 100

# Compositional plots at species level
# Summarize per group and identify top 20 species (baseline = V1)
top20_species <- baria_mb_df |> 
  # filter(visit == 1) |> # not filtered for longitudinal approach
  group_by(Sample, Species) |> 
  summarize(Abundance = sum(Abundance), .groups = "drop") |> 
  group_by(Species) |> 
  summarize(Abundance = mean(Abundance)) |> 
  arrange(-Abundance) |> 
  select(Species) |> 
  head(20) |> 
  print()

# low ASM
# per group (low ASM at baseline)
baria_mb_low_asm_v0 <- baria_mb_df |> 
  filter(visit == 1) |> # baseline composition (in clinical data v0 for shotgun data v1!)
  mutate(
    Species2 = if_else(
        Species %in% top20_species$Species,
        Species,
        "Other species" # collapse other species
        ),
    Species2 = as.factor(Species2)
  ) |> 
  group_by(Species2, Sample, low_asm_v0) |> # species abundance per sample
  summarize(Abundance = sum(Abundance)) |> 
  group_by(Species2, low_asm_v0) |> # species per muscle mass group
  summarize(Abundance = mean(Abundance)) |> # avg abundance per species per muscle group
  ungroup() |> 
  mutate(
    Species2 = fct_reorder(Species2, Abundance),
    Species2 = fct_relevel(Species2, "Other species", after = 0L) # move other spevies to the front
  ) |> 
  print()

baria_mb_low_asm_v0 |> # check
  group_by(low_asm_v0) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

levels <- levels(baria_mb_low_asm_v0$Species2)

# Composition plots
set.seed(13)
species_comp_low_asm_v0 <- baria_mb_low_asm_v0 |> 
  mutate(low_asm_v0 = fct_relevel(low_asm_v0, "yes", after = 0L)) |> # low asm first
  ggplot(aes(x = low_asm_v0, y = Abundance, fill = Species2)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = rev(c(sample(manet), "grey90")), labels = lev) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y="Composition (relative abundances)", x = "Low baseline ASM", title = "Microbiota composition", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition()

ggsave("graphs/species_comp_low_asm_v0.pdf", width = 12, height = 10)

# low SMM/W
# per group (low SMM/W at baseline)
baria_mb_low_smm_by_weight_v0 <- baria_mb_df |> 
  filter(visit == 1) |> # baseline composition (in clinical data v0 for shotgun data v1!)
  mutate(
    Species2 = if_else(
        Species %in% top20_species$Species,
        Species,
        "Other species" # collapse other species
        ),
    Species2 = as.factor(Species2)
  ) |> 
  group_by(Species2, Sample, low_smm_by_weight_v0) |> # species abundance per sample
  summarize(Abundance = sum(Abundance)) |> 
  group_by(Species2, low_smm_by_weight_v0) |> # species per muscle mass group
  summarize(Abundance = mean(Abundance)) |> # avg abundance per species per muscle group
  ungroup() |> 
  mutate(
    Species2 = fct_reorder(Species2, Abundance),
    Species2 = fct_relevel(Species2, "Other species", after = 0L) # move other spevies to the front
  ) |> 
  print()

baria_mb_low_smm_by_weight_v0 |> # check
  group_by(low_smm_by_weight_v0) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

levels <- levels(baria_mb_low_smm_by_weight_v0$Species2)

# Composition plots
set.seed(13)
species_comp_low_smm_by_weight <- baria_mb_low_smm_by_weight_v0 |> 
  mutate(low_asm_v0 = fct_relevel(low_smm_by_weight_v0, "yes", after = 0L)) |> # low asm first
  ggplot(aes(x = low_asm_v0, y = Abundance, fill = Species2)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = rev(c(sample(manet), "grey90")), labels = lev) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y="Composition (relative abundances)", x = "Low baseline SMM/Weight", title = "Microbiota composition", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition()

ggsave("graphs/species_comp_low_smm_by_weight_v0.pdf", width = 12, height = 10)
