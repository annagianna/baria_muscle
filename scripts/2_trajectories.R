# Weight and fat-free mass trajectories
# Anna Giannakogeorgou

# Packages
library(dplyr)
library(tidyverse)
library(ggpubr)
library(scales)
library(MetBrewer)
library(wesanderson)
library(patchwork)
library(gt)

# Data
baria_muscle <- read_rds("data/260131_BARIA_muscle_clinical.RDS")

# convert anthropometric data to long format
weight_long <- baria_muscle |> 
  pivot_longer(contains("weight_kg_v"), names_to = "visit", values_to = "weight_kg") |> 
  relocate(visit, weight_kg, .after = id) |> 
  mutate(visit = as.integer(parse_number(visit)))

ffm_long <- baria_muscle |> 
  select(id, contains("ffm_kg_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "ffm_kg") |> 
  mutate(visit = as.integer(parse_number(visit)))

fm_long <- baria_muscle |> 
  select(id, starts_with("fm_kg_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "fm_kg") |> 
  mutate(visit = as.integer(parse_number(visit)))

asm_long <- baria_muscle |> 
  select(id, contains("asm_kg_v"), -asm_kg_v0_tertile) |> 
  pivot_longer(-id, names_to = "visit", values_to = "asm_kg") |> 
  mutate(visit = as.integer(parse_number(visit)))

smm_long <- baria_muscle |> 
  select(id, contains("smm_kg_v"), -smm_kg_v0_tertile) |> 
  pivot_longer(-id, names_to = "visit", values_to = "smm_kg") |> 
  mutate(visit = as.integer(parse_number(visit)))

age_long <- baria_muscle |> 
  select(id, contains("age_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "age") |> 
  mutate(visit = as.integer(parse_number(visit)))

homa_ir_long <- baria_muscle |> 
  select(id, contains("homa_ir_v")) |> 
  pivot_longer(c(homa_ir_v0, homa_ir_v4, homa_ir_v5), names_to = "visit", values_to = "homa_ir") |> 
  mutate(visit = as.integer(parse_number(visit)))

homa_b_long <- baria_muscle |> 
  select(id, contains("homa_b")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "homa_b") |> 
  mutate(visit = as.integer(parse_number(visit)))

hba1c_mmolmol_long <- baria_muscle |> 
  select(id, contains("hba1c_mmolmol_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "hba1c_mmolmol") |> 
  mutate(visit = as.integer(stringr::str_extract(visit, "\\d+$"))) # parse_number() extracts the "1" in hba1c

hba1c_percent_long <- baria_muscle |> 
  select(id, contains("hba1c_percent_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "hba1c_percent") |> 
  mutate(visit = as.integer(stringr::str_extract(visit, "\\d+$"))) # same here

crp_long <- baria_muscle |> 
  select(id, contains("crp_mgl_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "crp_mgl") |> 
  mutate(visit = as.integer(parse_number(visit)))

# create time column based on difftime for trajectories
time <- baria_muscle |> 
  select(id, contains("numeric_y")) |> 
  rename(
    v2 = v0_to_v2_numeric_y,
    v3 = v0_to_v3_numeric_y,
    v4 = v0_to_v4_numeric_y,
    v5 = v0_to_v5_numeric_y,
    v6 = v0_to_v6_numeric_y,
    v7 = v0_to_v7_numeric_y) |>
  mutate(v0 = 0) |> 
  relocate(v0, .before = v2) |> 
  pivot_longer(
    v0:v7,
    names_to = "visit",
    values_to = "time_y"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# Check
time |> 
  group_by(visit) |> 
  summarize(mean = mean(time_y, na.rm = TRUE)) # aligned with the visits based on the baria protocol

# stepwise left joint
add_ffm <- left_join(weight_long, ffm_long, by = c("id", "visit"))
add_fm <- left_join(add_ffm, fm_long, by = c("id", "visit"))
add_asm <- left_join(add_fm, asm_long, by = c("id", "visit"))
add_smm <- left_join(add_asm, smm_long, by = c("id", "visit"))
add_age <- left_join(add_smm, age_long, by = c("id", "visit"))
add_homa_ir <- left_join(add_age, homa_ir_long, by = c("id", "visit"))
add_homa_b <- left_join(add_homa_ir, homa_b_long, by = c("id", "visit"))
add_hba1c_mmolmol <- left_join(add_homa_b, hba1c_mmolmol_long, by = c("id", "visit"))
add_hba1c_percent <- left_join(add_hba1c_mmolmol, hba1c_percent_long, by = c("id", "visit"))
add_crp_mgl <- left_join(add_hba1c_percent, crp_long, by = c("id", "visit"))
add_time <- left_join(add_crp_mgl, time, by = c("id", "visit"))

# clean up + calculate mean weight change from baseline
baria_muscle_long <- add_time |>
  select(
    id, visit, time_y, age, sex, t2d_v0, height_cm, low_asm_v0, low_asm_height2_v0, low_smm_v0, asm_change_v4_group, asm_change_v5_group,
    weight_kg, ffm_kg, fm_kg, asm_kg, smm_kg, homa_ir, homa_b, 
    crp_mgl, hba1c_mmolmol, hba1c_percent
  )|>
  arrange(visit) |>
  mutate(visit = case_when(visit == 0 ~ "Baseline", visit == 4 ~ "1y", visit == 5 ~ "2y", visit == 6 ~ "5y"))

# Plot trajectories
# Themes
gb1_colors <- scale_color_manual(values = wes_palette("GrandBudapest1"))
gb2_colors <- scale_color_manual(values = wes_palette("GrandBudapest2"))
mr3_colors <- scale_color_manual(values = wes_palette("Moonrise3"))
isle1_colors <- scale_color_manual(values = wes_palette("IsleofDogs1"))

### Weight trajectories ###
weight <- baria_muscle_long |>
  filter(!is.na(asm_change_v4_group)) |> 
  ggplot(aes(x = time_y, y = weight_kg, color = asm_change_v4_group)) +
  geom_smooth(se = FALSE, show.legend = FALSE) +
  labs(title = "Body weight", color = "%ASM loss at 1y", x = "Time (years)", y = "Weight (kg)") +
  theme_minimal()

ffm <- baria_muscle_long |>
  filter(!is.na(asm_change_v4_group)) |> 
  ggplot(aes(x = time_y, y = ffm_kg, color = asm_change_v4_group)) +
  geom_smooth(se = FALSE) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = "FFM", color = "%ASM loss at 1y", x = "Time (y)", y = "FFM (kg)") +
  theme_minimal()

fm <- baria_muscle_long |>
  filter(!is.na(asm_change_v4_group)) |> 
  ggplot(aes(x = time_y, y = fm_kg, color = asm_change_v4_group)) +
  geom_smooth(se = FALSE) +
  labs(title = "FM", color = "%ASM loss at 1y", x = "Time (y)", y = "FM (kg)") +
  theme_minimal()

weight / (ffm | fm ) +
  plot_layout(
    guides = "collect"
  ) &
  gb2_colors &
  theme(legend.position = "bottom")
ggsave(plot = asm_bw_1y_plot, filename = "asm_bw_1y.png", path = "graphs/trajectories")

### Scatterplot of %ASM loss vs. %BW loss ####
asm_vs_bw_1y_data <- baria_muscle_clean |> 
  filter(!is.na(asm_change_v4_group)) |> 
  mutate(
    perc_asm_change_1y = perc_asm_change_v4,
    perc_bw_change_1y  = 100 * (weight_kg_v4 - weight_kg_v0) / weight_kg_v0
  ) 

# perform correlation test
cor_asm_bw_1y <- cor.test(asm_vs_bw_1y_data$perc_asm_change_1y, asm_vs_bw_1y_data$perc_bw_change_1y, method = "spearman")
cor_asm_bw_1y_label <- sprintf("ρ = %.3f\np = %.3g", cor_asm_bw_1y$estimate, cor_asm_bw_1y$p.value)

# plot
asm_bw_1y_plot <- ggplot(asm_vs_bw_1y_data, aes(x = perc_asm_change_1y, y = perc_bw_change_1y)) +
  geom_point(aes(color = asm_change_v4_group)) +
  geom_smooth(aes(color = asm_change_v4_group), method = "lm", se = FALSE, color = "grey60", alpha = 0.5, linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", alpha = 0.5, linewidth = 0.8) + # reference 1:1 line
  gb2_colors +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "%ASM change at 1 year",
    y = "%BW change at 1 year"
  ) +
  annotate("text", x = Inf, y = -Inf, label = cor_asm_bw_1y_label, hjust = 1.1, vjust = -0.5, size = 4)
ggsave(plot = asm_bw_1y_plot, filename = "asm_bw_1y.png", path = "graphs/trajectories")

### Proportion area plot ffm/fm/weight ###
# pivot longer to split weight into its components (ffm/fm)
weight_comps <- baria_muscle_long |> 
  select(visit, time_y, weight_kg, ffm_kg, fm_kg) |> 
  filter(!is.na(visit)) |> 
  filter(!is.na(time_y)) |> 
  pivot_longer(
    cols = c(ffm_kg, fm_kg),
    names_to = "weight_component",
    values_to = "mass_kg"
  ) |> 
  print()

# plot 
weight_comps |> 
  mutate(
    weight_component = as.factor(weight_component),
    weight_component = fct_relevel(weight_component, "fm_kg", after = 0L)
  ) |> 
  ggplot(aes(x = time_y, y = mass_kg, fill = weight_component)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(
    values = wes_palette("GrandBudapest2", n = 2),
    labels = c(ffm_kg = "Fat-free mass", fm_kg  = "Fat mass")
  ) +
  labs(
    x = "Time (years)",
    y = "Body weight (kg)",
    fill = "Weight component"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
### need to smoothen this over 
