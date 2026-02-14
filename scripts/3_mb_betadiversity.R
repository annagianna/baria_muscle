# Baria muscle loss & microbiota project: Beta diversity
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(stringr)
library(phyloseq)
library(microbiome)
library(vegan)
library(MetBrewer)
library(grid)
library(ggthemes)
library(ggpubr)
library(patchwork)
library(ggsci)

# Theme
manet_cols <- met.brewer("Manet", n = 20)
fill_cols_asm1y <- scale_fill_manual(values = c("high" = manet_cols[10], "low/modest" = manet_cols[20]))
color_cols_asm1y <- scale_color_manual(values = c("high" = manet_cols[10], "low/modest" = manet_cols[20]))

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
sample_sums(baria_mb) # adds up to 100

# convert sample_data of baria_mb to a data frame to merge with selected cols from baria_muscle
meta_baria_mb <- meta(baria_mb) |> 
  tibble::rownames_to_column(var = "rownames") |> 
  mutate(
    visit = as.integer(parse_number(Time_Point)),
    visit = if_else(visit == -1, 0, visit),
    Extra_data = na_if(Extra_data, "NA")
  )
class(meta_baria_mb) # is df

# Add %ASM change groups to sample_data() of phyloseq object
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

## Calculate Bray-Curtis Distance
# select only needed cols for bray-curtis calculation
bray_cols <- baria_mb_df |> 
  select(Sample, OTU, Abundance)

# Pivot wider for bray-curtis calculation
bray_wide <- bray_cols |> 
  pivot_wider(
    names_from = OTU,
    values_from = Abundance,
    values_fill = 0 # if Abundance 0
  ) |> 
  tibble::column_to_rownames(var = "Sample")

# Make metadata for bray
bray_meta <- baria_mb_df |> 
  select(Sample, id, visit, sex, asm_change_v4_group) |> 
  filter(!is.na(asm_change_v4_group)) |> 
  distinct(Sample, .keep_all = TRUE)

### Compute Bray-Curtis distance ###
bray <- vegan::vegdist(bray_wide[rownames(bray_wide) %in% bray_meta$Sample, ], method = "bray")
pcoord <- ape::pcoa(bray, correction = "cailliez")
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100

bray2 <- as.data.frame(pcoord$vectors[, c("Axis.1", "Axis.2")]) |>
  tibble::rownames_to_column("Sample") |>
  left_join(bray_meta, by = "Sample") # add metadata

### Bray-Curtis per visit ###
# set up function
bray_per_v <- function(v, df = bray_meta, abd_tab = bray_wide) {
  set.seed(1312)

  meta_v <- df |> 
    filter(visit == v) |> 
    distinct(Sample, .keep_all = TRUE) # drop all other cols

  abd_tab_v <- abd_tab[meta_v$Sample, , drop = FALSE] # keep only Samples from that visit 

  bray_v <- vegan::vegdist(abd_tab_v, method = "bray", na.rm = TRUE)

  vegan::adonis2(bray_v ~ asm_change_v4_group, data = meta_v, permutations = 999, by = "terms")
}

# run it for each visit
visits <- sort(unique(bray_meta$visit)) 
results_per_v <- lapply(visits, bray_per_v)
names(results_per_v) <- visits

## Plots ##
# Baseline
bray_asm1y_v0 <- bray2 |> 
  filter(visit == 0) |> 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(color = asm_change_v4_group), size = 3, alpha = 0.9) +
  xlab(paste0("PCo1 (", round(expl_variance_bray[1], digits = 1),"%)")) +
  ylab(paste0("PCo2 (", round(expl_variance_bray[2], digits = 1),"%)")) +
  labs(color = "", fill = "", title = "PCoA Bray-Curtis Distance") +
  stat_ellipse(
    geom = "polygon", 
    aes(color = asm_change_v4_group, fill = asm_change_v4_group),
    type = "norm", alpha = 0.13, linewidth = 1.2) +
  fill_cols_asm1y +
  color_cols_asm1y +
  theme_Publication() +
  theme(legend.position = "top") +
  geom_text(
    data = as.data.frame(results_per_v$`0`),
    aes(x = Inf, y = Inf,
      label = paste0("p = ", round(`Pr(>F)`, 3))),
      hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE)
ggsave(bray_asm1y_v0, filename = "graphs/betadiversity/bray_asm1y_v0.pdf", width = 8, height = 8)

# 1y
bray_asm1y_v4 <- bray2 |> 
  filter(visit == 4) |> 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(color = asm_change_v4_group), size = 3, alpha = 0.9) +
  xlab(paste0("PCo1 (", round(expl_variance_bray[1], digits = 1),"%)")) +
  ylab(paste0("PCo2 (", round(expl_variance_bray[2], digits = 1),"%)")) +
  labs(color = "", fill = "", title = "PCoA Bray-Curtis Distance") +
  stat_ellipse(
    geom = "polygon", 
    aes(color = asm_change_v4_group, fill = asm_change_v4_group),
    type = "norm", alpha = 0.13, linewidth = 1.2) +
  fill_cols_asm1y +
  color_cols_asm1y +
  theme_Publication() +
  theme(legend.position = "top") +
  geom_text(
    data = as.data.frame(results_per_v$`0`),
    aes(x = Inf, y = Inf,
      label = paste0("p = ", round(`Pr(>F)`, 3))),
      hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE)
ggsave(bray_asm1y_v4, filename = "graphs/betadiversity/bray_asm1y_v4.pdf", width = 8, height = 8)
