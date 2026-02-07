# Baria: associations between low muscle mass and gut microbiota (Baseline)
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(stringr)
library(MetBrewer)
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
baria_muscle <- read_rds("data/260206_BARIA_muscle_clinical.RDS") # metadata/clinical data CHECK FOR MOST RECENT VERSION
baria_mb <- read_rds("data/ps.BARIA.metaphlan.706.2548.RDS")
sample_sums(baria_mb) # adds up to ~100

# convert sample_data of baria_mb phyloseq obj. to a data frame to merge with selected cols from baria_muscle
meta_baria_mb <- meta(baria_mb) |> 
  tibble::rownames_to_column(var = "rownames") |> 
  mutate(
    visit = as.integer(parse_number(Time_Point)),
    visit = if_else(visit == -1, 0, visit),
    Extra_data = na_if(Extra_data, "NA")
  )
class(meta_baria_mb) # is df
colnames(baria_muscle)
# Add low muscle variables to sample_data() of phyloseq object for grouping
new_sample_data <- baria_muscle |> 
  select(id, sex, asm_change_v4_group) |> 
  mutate(id = as.character(id)) |> 
  inner_join(meta_baria_mb, join_by(id == Subject_ID)) |> # only keep ids with BIA and microbiome data
  select(-(Study:Time_Point)) |> 
  group_by(id, visit) |>
  slice(1) |> # keeps first occurence per group (some ids had second runs)
  relocate(rownames, .before = id) |> 
  relocate(visit, .after = id) |> 
  relocate(Extra_data, .after = visit) |> 
  tibble::column_to_rownames(var = "rownames")

# update sample_data() within baria_mb
sample_data(baria_mb) <- sample_data(new_sample_data)

# melt phyloseq object into a df
baria_mb_df <- psmelt(baria_mb)
baria_mb_df |> # check
  group_by(Sample) |> 
  summarise(sum_abundance = sum(Abundance)) # adds up to 100

# Summarize per group and identify top 20 species (baseline = V1) 
top20_species <- baria_mb_df |> 
  # filter(visit == 0) |> # not filtered for longitudinal approach
  group_by(Sample, Species) |> 
  dplyr::summarize(Abundance = sum(Abundance), .groups = "drop") |> 
  group_by(Species) |> 
  dplyr::summarize(Abundance = mean(Abundance)) |> 
  arrange(-Abundance) |> 
  select(Species) |> 
  head(20) |> 
  print()

# assign fixed color to each one of the top20 species an
top20_species_vector <- top20_species$Species # top 20 species names as vector
set.seed(13)
top20_species_colours <- c(
  "Other species" = "grey90",
  setNames(sample(manet), top20_species_vector)
)


#### Statified by %ASM change group at 1y ####
### Baseline composition ###
mb_asm1y_v0 <- baria_mb_df |> 
  filter(
    visit == 0, # baseline composition (in clinical data v0 for shotgun data v1!)
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(
    Species2 = if_else(
        Species %in% top20_species$Species,
        Species,
        "Other species" # collapse other species
        ),
    Species2 = as.factor(Species2)
  ) |> 
  group_by(Species2, Sample, asm_change_v4_group) |> # species abundance per sample
  dplyr::summarize(Abundance = sum(Abundance)) |> 
  group_by(Species2, asm_change_v4_group) |> # species per muscle mass group
  dplyr::summarize(Abundance = mean(Abundance)) |> # avg abundance per species per muscle group
  ungroup() |>
  mutate(
    Species2 = fct_reorder(Species2, Abundance),
    Species2 = fct_relevel(Species2, "Other species", after = 0L) # move other species to the front
  )

mb_asm1y_v0 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
species_comp_asm1y_v0 <- mb_asm1y_v0 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Species2)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = top20_species_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance", x = "%ASM change at 1y", title = "Microbiota composition", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition()
ggsave("graphs/composition/species_comp_asm1y_v0.pdf", width = 12, height = 10)

### 1y composition ###
baria_mb_df$visit
mb_asm1y_v4 <- baria_mb_df |> 
  filter(
    visit == 4, # 1y
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(
    Species2 = if_else(
        Species %in% top20_species$Species,
        Species,
        "Other species" # collapse other species
        ),
    Species2 = as.factor(Species2)
  ) |> 
  group_by(Species2, Sample, asm_change_v4_group) |> # species abundance per sample
  dplyr::summarize(Abundance = sum(Abundance)) |> 
  group_by(Species2, asm_change_v4_group) |> # species per muscle mass group
  dplyr::summarize(Abundance = mean(Abundance)) |> # avg abundance per species per muscle group
  ungroup() |>
  mutate(
    Species2 = fct_reorder(Species2, Abundance),
    Species2 = fct_relevel(Species2, "Other species", after = 0L) # move other species to the front
  )

mb_asm1y_v4 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
species_comp_asm1y_v4 <- mb_asm1y_v4 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Species2)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = top20_species_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance", x = "%ASM change at 1y", title = "Microbiota composition at 1y", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition()
ggsave("graphs/composition/species_comp_asm1y_v4.pdf", width = 12, height = 10)
