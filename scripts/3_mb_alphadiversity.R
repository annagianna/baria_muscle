# BARIA: Associations between muscle mass and gut microbiota
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(MetBrewer)
library(grid)
library(ggthemes)
library(ggpubr)
library(patchwork)

# Theme
manet_cols <- met.brewer("Manet", n = 20)
fill_cols_2 <- scale_fill_manual(values = c("high" = manet_cols[10], "low/modest" = manet_cols[20]))

theme_Publication <- function(base_size = 14, base_family = "sans") {
  
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
baria_muscle <- readRDS("data/20260613_BARIA_muscle_clinical.RDS") # metadata/clinical data
baria_mb <- readRDS("data/ps.BARIA.metaphlan.706.2548.RDS")
sample_sums(baria_mb) # adds up to 100

# Diversity metrics
# Convert OTU table to matrix
matrix_mb <- as(otu_table(baria_mb), "matrix")

# vegan requires a matrix with samples as rows and taxa as cols
if (taxa_are_rows(baria_mb)) { 
  matrix_mb <- t(matrix_mb) 
}

# Shannon, Simpson, Richness
shannon <- vegan::diversity(matrix_mb, index = "shannon")
simpson <- vegan::diversity(matrix_mb, index = "simpson")
richness <- vegan::specnumber(matrix_mb)

# Merge with metadata
alpha <- tibble(
  Sample = names(shannon),
  shannon = shannon,
  simpson = simpson,
  richness = richness
)|> 
  left_join(
    as(sample_data(baria_mb), "data.frame") |> 
      tibble::rownames_to_column(var = "Sample") |> 
      mutate(
        visit = str_extract(Time_Point, "\\d"), # visit number without the "V"
        visit = if_else(visit == "1", "0", visit) # Baseline as v0 to match metadata
      ),
      select(Sample, id, visit),
      by = "Sample"
  ) |> 
  inner_join(baria_muscle, by = join_by(Subject_ID == id)) |> 
  mutate(id = Subject_ID) |> 
  relocate(id, .before = Sample)

nrow(alpha)
count(alpha, visit)
table(alpha$low_ffmi_v0)

#### Baseline Plots ####
## Shannon ##
# Boxplot
shannon_ffmi_v0 <- shannon |> 
  filter(
    visit == 0,
    !is.na(asm_change_v4_group)
  ) |> 
  left_join(baria_muscle, ) |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> 
  ggplot(aes(x = asm_change_v4_group, y = shannon, fill = asm_change_v4_group)) +
  geom_boxplot() +
  #geom_jitter(position = position_dodge(0.75)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "Shannon index", y = "Shannon index", x = "%ASM change at 1y") +
  fill_cols_asm1y +
  labs(fill = "%ASM change at 1y") + 
  theme_Publication()
ggsave(shannon_asm1y_box_v0, filename = "graphs/alphadiversity/shannon_asm1y_box_v0.pdf", width = 7, height = 5)

# Violin
shannon_asm1y_violin_v0 <- baria_mb_alpha |> 
  filter(
    visit == 0,
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> 
  ggplot(aes(x = asm_change_v4_group, y = shannon)) +
  geom_violin(aes(fill = asm_change_v4_group)) +
  geom_boxplot(fill = "white", width = 0.1) +
  labs(x = "", y = "Shannon index", title = "Shannon index", fill = "%ASM change at 1y") +
  stat_compare_means( 
    tip.length = 0, 
    hide.ns = TRUE, 
    label.x = 1.5,
    method = "wilcox.test",
    label = "p.signif",
    label.y = max(baria_mb_alpha$shannon, na.rm = TRUE) * 0.99
  ) + 
  fill_cols_asm1y +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_Publication()
ggsave(shannon_asm1y_violin_v0 , filename = "graphs/alphadiversity/shannon_asm1y_violin_v0.pdf", width = 6, height = 5)

## Simpson ##

## Richness ##

## Combine into a panel
#### OLD CODE ###  CHANGE FROM HERE
# Baseline (Shannon, Simpson, Richness)
# Violins
alpha_panel_asm1y_v0 <- 
  (shannon_asm1y_violin_v0 + simpson_asm1y_violin_v0 + richness_asm1y_violin_v0) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  theme(aspect.ratio = 0.8)
ggsave(alpha_panel_asm1y_v0, filename = "graphs/alphadiversity/alpha_panel_asm1y_v0.pdf", width = 12, height = 8)
