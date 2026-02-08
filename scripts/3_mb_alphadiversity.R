# Baria: associations between muscle loss and gut microbiota (Baseline)
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(stringr)
library(MetBrewer)
library(grid)
library(ggthemes)
library(ggpubr)

# Theme
manet_cols <- met.brewer("Manet", n = 20)
fill_cols_asm1y <- scale_fill_manual(values = c("high" = manet_cols[10], "low/modest" = manet_cols[20]))

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

# Add low muscle variables to sample_data() of phyloseq object
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

# Diversity metrics
# Shannon
shannon <- microbiome::alpha(baria_mb, index = "shannon")
shannon_df <- shannon |> 
  tibble::rownames_to_column(var = "rownames") |>
  rename(Sample = rownames)

# join shannon with baria_mb_df
baria_mb_shannon <- baria_mb_df |> 
  inner_join(shannon_df, by = join_by(Sample)) |> 
  rename(shannon = diversity_shannon)

# Simpson
simpson <- microbiome::alpha(baria_mb, index = "simpson")
simpson_df <- simpson |> 
  mutate(simpson = (1 - dominance_simpson)) |> 
  tibble::rownames_to_column(var = "rownames") |> 
  rename(Sample = rownames)

# join simpson with baria_mb_shannon
baria_mb_shannon_simpson <- baria_mb_shannon |> 
  inner_join(simpson_df, by = join_by(Sample))

# Richness
richness <- microbiome::alpha(baria_mb, index = "observed")
richness_df <- richness |> 
  tibble::rownames_to_column(var = "rownames") |>
  rename(Sample = rownames)

# join richness with baria_mb_shannon_simpson
baria_mb_alpha <-  baria_mb_shannon_simpson |> 
  inner_join(richness_df, by = join_by(Sample)) |> 
  rename(richness = observed)

#### Plots ####
# Shannon boxplots 
# Baseline
shannon_asm1y_box_v0 <- baria_mb_alpha |> 
  filter(
    visit == 0,
    !is.na(asm_change_v4_group)
  ) |> 
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

# 1y
shannon_asm1y_box_v4 <- baria_mb_shannon |> 
  filter(
    visit == 4,
    !is.na(asm_change_v4_group)
  ) |> 
  mutate(asm_change_v4_group = fct_relevel(asm_change_v4_group, "high", after = 0L)) |> 
  ggplot(aes(x = asm_change_v4_group, y = shannon, fill = asm_change_v4_group)) +
  geom_boxplot() +
  #geom_jitter(position = position_dodge(0.75)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "Shannon index", y = "Shannon index", x = "%ASM change at 1y") +
  fill_cols_asm1y +
  labs(fill = "%ASM change at 1y") + 
  theme_Publication()
ggsave(shannon_asm1y_box_v4, filename = "graphs/alphadiversity/shannon_asm1y_box_v4.pdf", width = 7, height = 5)

# Shannon violin Plots
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
    label = "p.signif"
  ) +
  fill_cols_asm1y +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  theme_Publication()
ggsave(shannon_asm1y_violin_v0 , filename = "graphs/alphadiversity/shannon_asm1y_violin_v0.pdf", width = 6, height = 5)

# 1y
shannon_asm1y_violin_v4 <- baria_mb_alpha |> 
  filter(
    visit == 4,
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
    label = "p.signif"
  ) +
  fill_cols_asm1y +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  theme_Publication()
ggsave(shannon_asm1y_violin_v4 , filename = "graphs/alphadiversity/shannon_asm1y_violin_v0.pdf", width = 6, height = 5)


## continue here ... Simpson and richness ... 


# Combine into a panel