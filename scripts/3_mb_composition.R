# Baria: associations between muscle mass loss and gut microbiota (Baseline)
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(microbiome)
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
baria_muscle <- read_rds("data/260217_BARIA_muscle_clinical.RDS") # metadata/clinical data CHECK FOR MOST RECENT VERSION
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

# Composition Plots
#### Statified by %ASM change group at 1y ####
### Phylum level ###
top5_phyla <- baria_mb_df |> 
  # filter(visit == 0) |>
  group_by(Sample, Phylum) |> 
  dplyr::summarize(Abundance = sum(Abundance), .groups = "drop") |> 
  group_by(Phylum) |> 
  dplyr::summarize(Abundance = mean(Abundance)) |> 
  arrange(-Abundance) |> 
  select(Phylum) |> 
  head(5) |> 
  print()

# assign fixed color to each one of the top 5 phyla
top5_phyla_vector <- top5_phyla$Phylum # top 5 phyla names as vector
set.seed(23)
top5_phyla_colours <- c("Other phyla" = "grey63", setNames(sample(manet_5), top5_phyla_vector))

### Baseline composition ###
mb_phyla_asm1y_v0 <- baria_mb_df |> 
  filter(
    visit == 0, # baseline composition (in clinical data v0 for shotgun data v1!)
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(
    Phylum2 = if_else(
        Phylum %in% top5_phyla$Phylum,
        Phylum,
        "Other phyla" # collapse other phyla
        ),
    Phylum2 = as.factor(Phylum2)
  ) |> 
  group_by(Phylum2, Sample, asm_change_v4_group) |> # abundance per sample
  dplyr::summarize(Abundance = sum(Abundance)) |> 
  group_by(Phylum2, asm_change_v4_group) |> # per %asm change group at 1y
  dplyr::summarize(Abundance = mean(Abundance)) |> # avg abundance per phylum pergroup
  ungroup() |>
  mutate(
    Phylum2 = fct_reorder(Phylum2, Abundance),
    Phylum2 = fct_relevel(Phylum2, "Other phyla", after = 0L) # move other species to the front
  )

mb_phyla_asm1y_v0 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
phyla_comp_asm1y_v0 <- mb_phyla_asm1y_v0 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Phylum2)) +
  geom_bar(stat = "identity", color = "black", width = 0.9) +
  scale_fill_manual(values = top5_phyla_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance (%)", x = "%ASM change at 1y", title = "Phylum", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("graphs/composition/phyla_comp_asm1y_v0.pdf", width = 12, height = 10)

### 1y composition ###
mb_phyla_asm1y_v4 <- baria_mb_df |> 
  filter(
    visit == 4, # 1y
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(
    Phylum2 = if_else(
        Phylum %in% top5_phyla$Phylum,
        Phylum,
        "Other phyla"
        ),
    Phylum2 = as.factor(Phylum2)
  ) |> 
  group_by(Phylum2, Sample, asm_change_v4_group) |> # phyla abundance per sample
  dplyr::summarize(Abundance = sum(Abundance)) |> 
  group_by(Phylum2, asm_change_v4_group) |> # phyla per %asm change at 1 year group
  dplyr::summarize(Abundance = mean(Abundance)) |> # avg abundance per phylum per group
  ungroup() |>
  mutate(
    Phylum2 = fct_reorder(Phylum2, Abundance),
    Phylum2 = fct_relevel(Phylum2, "Other phyla", after = 0L) # move other phyla to the front
  )

mb_phyla_asm1y_v4 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
phyla_comp_asm1y_v4 <- mb_phyla_asm1y_v4 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # high asm loss first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Phylum2)) +
  geom_bar(stat = "identity", color = "black", width = 0.9) +
  scale_fill_manual(values = top5_phyla_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance (%)", x = "%ASM change at 1y", title = "Phylum", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("graphs/composition/phyla_comp_asm1y_v4.pdf", width = 12, height = 10)

### Genus level ###
top20_genera <- baria_mb_df |> 
  # filter(visit == 0) |>
  group_by(Sample, Genus) |> 
  dplyr::summarize(Abundance = sum(Abundance), .groups = "drop") |> 
  group_by(Genus) |> 
  dplyr::summarize(Abundance = mean(Abundance)) |> 
  arrange(-Abundance) |> 
  select(Genus) |> 
  head(20) |> 
  print()

# assign fixed color to each one of the top 20 genera
top20_genera_vector <- top20_genera$Genus # names of top 20 genera as vector
set.seed(18)
top20_genera_colours <- c("Other genera" = "grey63", setNames(sample(manet_20), top20_genera_vector))

### Baseline composition ###
mb_genera_asm1y_v0 <- baria_mb_df |> 
  filter(
    visit == 0,
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(
    Genus2 = if_else(
        Genus %in% top20_genera$Genus,
        Genus,
        "Other genera" # collapse other genera
        ),
    Genus2 = as.factor(Genus2)
  ) |> 
  group_by(Genus2, Sample, asm_change_v4_group) |> # abundance per sample
  dplyr::summarize(Abundance = sum(Abundance)) |> 
  group_by(Genus2, asm_change_v4_group) |> # per %asm change group at 1y
  dplyr::summarize(Abundance = mean(Abundance)) |> # avg abundance per genuys per group
  ungroup() |>
  mutate(
    Genus2 = fct_reorder(Genus2, Abundance),
    Genus2 = fct_relevel(Genus2, "Other genera", after = 0L) # move other genera to the front
  )

mb_genera_asm1y_v0 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
genera_comp_asm1y_v0 <- mb_genera_asm1y_v0 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Genus2)) +
  geom_bar(stat = "identity", color = "black", width = 0.9) +
  scale_fill_manual(values = top20_genera_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance (%)", x = "%ASM change at 1y", title = "Genus", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("graphs/composition/genera_comp_asm1y_v0.pdf", width = 12, height = 10)

### 1y composition ###
mb_genera_asm1y_v4 <- baria_mb_df |> 
  filter(
    visit == 4, # 1y
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(
    Genus2 = if_else(
        Genus %in% top20_genera$Genus,
        Genus,
        "Other genera"
        ),
    Genus2 = as.factor(Genus2)
  ) |> 
  group_by(Genus2, Sample, asm_change_v4_group) |> # abundance per sample
  dplyr::summarize(Abundance = sum(Abundance)) |> 
  group_by(Genus2, asm_change_v4_group) |> # genera per %asm change at 1 year
  dplyr::summarize(Abundance = mean(Abundance)) |> # avg abundance per genus per group
  ungroup() |>
  mutate(
    Genus2 = fct_reorder(Genus2, Abundance),
    Genus2 = fct_relevel(Genus2, "Other genera", after = 0L)
  )

mb_genera_asm1y_v4 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
genera_comp_asm1y_v4 <- mb_genera_asm1y_v4 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Genus2)) +
  geom_bar(stat = "identity", color = "black", width = 0.9) +
  scale_fill_manual(values = top20_genera_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance (%)", x = "%ASM change at 1y", title = "Genus", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("graphs/composition/genera_comp_asm1y_v4.pdf", width = 12, height = 10)

### Species level ###
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

# assign fixed color to each one of the top20 species
top20_species_vector <- top20_species$Species # top 20 species names as vector
set.seed(13)
top20_species_colours <- c("Other species" = "grey63", setNames(sample(manet_20), top20_species_vector))

### Baseline composition ###
mb_species_asm1y_v0 <- baria_mb_df |> 
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

mb_species_asm1y_v0 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
species_comp_asm1y_v0 <- mb_species_asm1y_v0 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Species2)) +
  geom_bar(stat = "identity", color = "black", width = 0.9) +
  scale_fill_manual(values = top20_species_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance (%)", x = "%ASM change at 1y", title = "Species", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("graphs/composition/species_comp_asm1y_v0.pdf", width = 12, height = 10)

### 1y composition ###
mb_species_asm1y_v4 <- baria_mb_df |> 
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

mb_species_asm1y_v4 |> # check
  group_by(asm_change_v4_group) |> 
  summarise(sum_Abundance = sum(Abundance)) # adds up to 100

# Composition plot
species_comp_asm1y_v4 <- mb_species_asm1y_v4 |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> # low asm/height2 first
  ggplot(aes(x = asm_change_v4_group, y = Abundance, fill = Species2)) +
  geom_bar(stat = "identity", color = "black", width = 0.9) +
  scale_fill_manual(values = top20_species_colours) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = "Relative abundance (%)", x = "%ASM change at 1y", title = "Species", fill = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_composition() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("graphs/composition/species_comp_asm1y_v4.pdf", width = 12, height = 10)

# Combine phylum + genus + species plots into one panel
# Dummy plot /shared x axis 
shared_x_asm1y <- ggplot() +
  labs(x = "%ASM change at 1y") +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = -20, r = 5, b = 5, l = 10)
  )

# Baseline composition panel
comp_pgs_v0 <- 
  (phyla_comp_asm1y_v0 + genera_comp_asm1y_v0 + species_comp_asm1y_v0) /
  shared_x_asm1y +
  plot_layout(heights = c(1, 0.005))
ggsave("graphs/composition/comp_pgs_v0.pdf", width = 12, height = 7)

# 1 year composition panel
comp_pgs_v4 <- 
  (phyla_comp_asm1y_v4 + genera_comp_asm1y_v4 + species_comp_asm1y_v4) /
  shared_x_asm1y +
  plot_layout(heights = c(1, 0.005))
ggsave("graphs/composition/comp_pgs_v4.pdf", width = 12, height = 7)
