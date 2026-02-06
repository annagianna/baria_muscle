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
baria_muscle <- read_rds("data/260206_BARIA_muscle_clinical.RDS")

# convert body composition data to long format
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

perc_weight_change_long <- baria_muscle |> 
  select(id, contains("perc_weight_change_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "perc_weight_change") |> 
  mutate(visit = as.integer(parse_number(visit)))

perc_ffm_change_long <- baria_muscle |> 
  select(id, contains("perc_ffm_change_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "perc_ffm_change") |> 
  mutate(visit = as.integer(parse_number(visit)))

perc_fm_change_long <- baria_muscle |> 
  select(id, contains("perc_fm_change_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "perc_fm_change") |> 
  mutate(visit = as.integer(parse_number(visit)))

perc_asm_change_long <- baria_muscle |> 
  select(id, contains("perc_asm_change_v")) |> 
  pivot_longer(-id, names_to = "visit", values_to = "perc_asm_change") |> 
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
  dplyr::summarize(mean = mean(time_y, na.rm = TRUE)) # aligned with the visits based on the baria protocol

# stepwise left joint
add_ffm <- left_join(weight_long, ffm_long, by = c("id", "visit"))
add_fm <- left_join(add_ffm, fm_long, by = c("id", "visit"))
add_asm <- left_join(add_fm, asm_long, by = c("id", "visit"))
add_smm <- left_join(add_asm, smm_long, by = c("id", "visit"))
add_perc_weight_change_long <- left_join(add_smm, perc_weight_change_long, by = c("id", "visit"))
add_perc_ffm_change_long <- left_join(add_perc_weight_change_long, perc_ffm_change_long, by = c("id", "visit"))
add_perc_fm_change_long <- left_join(add_perc_ffm_change_long, perc_fm_change_long, by = c("id", "visit"))
add_perc_asm_change_long <- left_join(add_perc_fm_change_long, perc_asm_change_long, by = c("id", "visit"))
add_age <- left_join(add_perc_asm_change_long, age_long, by = c("id", "visit"))
add_homa_ir <- left_join(add_age, homa_ir_long, by = c("id", "visit"))
add_homa_b <- left_join(add_homa_ir, homa_b_long, by = c("id", "visit"))
add_hba1c_mmolmol <- left_join(add_homa_b, hba1c_mmolmol_long, by = c("id", "visit"))
add_hba1c_percent <- left_join(add_hba1c_mmolmol, hba1c_percent_long, by = c("id", "visit"))
add_crp_mgl <- left_join(add_hba1c_percent, crp_long, by = c("id", "visit"))
add_time <- left_join(add_crp_mgl, time, by = c("id", "visit"))

# clean up
baria_muscle_long <- add_time |>
  select(
    id, visit, time_y, age, sex, t2d_v0, height_cm, low_asm_v0, low_asm_height2_v0, low_smm_v0, asm_change_v4_group, asm_change_v5_group,
    weight_kg, ffm_kg, fm_kg, asm_kg, smm_kg, homa_ir, homa_b, crp_mgl, hba1c_mmolmol, hba1c_percent,
    perc_weight_change, perc_ffm_change, perc_fm_change, perc_asm_change
  )|>
  arrange(visit) |>
  mutate(visit = case_when(visit == 0 ~ "Baseline", visit == 4 ~ "1y", visit == 5 ~ "2y", visit == 6 ~ "5y"))

# Plot trajectories
### Themes ###
manet_cols <- met.brewer("Manet", n = 20)

### Trajectories ###
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

trajectories_asm_1y <-
  weight / (ffm | fm ) +
  plot_layout(
    guides = "collect"
  ) &
  scale_color_manual(values = c("high" = manet_cols[5], "low/modest" = manet_cols[20])) &
  theme(legend.position = "bottom")

ggsave(plot = trajectories_asm_1y, filename = "trajectories_asm_1y.png", path = "graphs/trajectories")

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
  geom_smooth(aes(color = asm_change_v4_group), method = "lm", se = FALSE, color = "grey60", alpha = 0.5, linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", alpha = 0.5, linewidth = 0.5) + # reference 1:1 line
  scale_color_manual(values = c("high" = manet_cols[5], "low/modest" = manet_cols[20])) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "%ASM change at 1 year",
    y = "%BW change at 1 year"
  ) +
  annotate("text", x = Inf, y = -Inf, label = cor_asm_bw_1y_label, hjust = 1.1, vjust = -0.5, size = 4)
ggsave(plot = asm_bw_1y_plot, filename = "asm_bw_1y.png", path = "graphs/trajectories")


### % change BW, FFM, FM ###
# %BW change stratified by %ASM loss at 1 year
plot_stat_data <- baria_muscle_long |> 
  filter(!is.na(visit)) |> 
  filter(!is.na(asm_change_v4_group)) |>
  mutate(visit = factor(visit, levels = c("Baseline", "1y", "2y", "5y"))) |>
  filter(visit != "Baseline")

perc_bw_asm1y <- baria_muscle_long |> 
  filter(!is.na(visit)) |> 
  filter(!is.na(asm_change_v4_group)) |>
  mutate(visit = factor(visit, levels = c("Baseline", "1y", "2y", "5y"))) |>
  filter(visit != "Baseline") |> 
  group_by(visit, asm_change_v4_group) |> 
  dplyr::summarize(
    mean = mean(perc_weight_change, na.rm = TRUE),
    se = sd(perc_weight_change, na.rm = TRUE) / sqrt(sum(!is.na(perc_weight_change))),
    .groups = "drop"
  ) |> 
  ggplot(aes(x = visit, y = mean, group = asm_change_v4_group)) +
  geom_col(aes(fill = asm_change_v4_group), position = position_dodge(width = 0.5), width = 0.5, color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(width = 0.5), width = 0.2) +
  stat_compare_means(
    data = plot_stat_data,
    aes(x = visit, y = perc_weight_change, group = asm_change_v4_group),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  scale_fill_manual(
    values = c("high" = manet_cols[5], "low/modest" = manet_cols[20]),
    name = "%ASM change at 1y"
  ) +
  scale_x_discrete() +
  labs(x = "Time", y = "Mean %BW change") +
  theme_minimal()

# %FFM change stratified by %ASM loss at 1 year
perc_ffm_asm1y <- baria_muscle_long |> 
  filter(!is.na(visit)) |> 
  filter(!is.na(asm_change_v4_group)) |>
  mutate(visit = factor(visit, levels = c("Baseline", "1y", "2y", "5y"))) |>
  filter(visit != "Baseline") |> 
  group_by(visit, asm_change_v4_group) |> 
  dplyr::summarize(
    mean = mean(perc_ffm_change, na.rm = TRUE),
    se = sd(perc_ffm_change, na.rm = TRUE) / sqrt(sum(!is.na(perc_ffm_change))),
    .groups = "drop"
  ) |> 
  ggplot(aes(x = visit, y = mean, group = asm_change_v4_group)) +
  geom_col(aes(fill = asm_change_v4_group), position = position_dodge(width = 0.5), width = 0.5, color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(width = 0.5), width = 0.2) +
  stat_compare_means(
    data = plot_stat_data,
    aes(x = visit, y = perc_ffm_change, group = asm_change_v4_group),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  scale_fill_manual(
    values = c("high" = manet_cols[5], "low/modest" = manet_cols[20]),
    name = "%ASM change at 1y"
  ) +
  scale_x_discrete() +
  labs(x = "Time", y = "Mean %FFM change") +
  theme_minimal() +
  coord_cartesian(ylim = c(-20, 0))

# %FM change stratified by %ASM loss at 1 year
perc_fm_asm1y <- baria_muscle_long |> 
  filter(!is.na(visit)) |> 
  filter(!is.na(asm_change_v4_group)) |>
  mutate(visit = factor(visit, levels = c("Baseline", "1y", "2y", "5y"))) |>
  filter(visit != "Baseline") |> 
  group_by(visit, asm_change_v4_group) |> 
  dplyr::summarize(
    mean = mean(perc_fm_change, na.rm = TRUE),
    se = sd(perc_fm_change, na.rm = TRUE) / sqrt(sum(!is.na(perc_fm_change))),
    .groups = "drop"
  ) |> 
  ggplot(aes(x = visit, y = mean, group = asm_change_v4_group)) +
  geom_col(aes(fill = asm_change_v4_group), position = position_dodge(width = 0.5), width = 0.5, color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(width = 0.5), width = 0.2) +
  stat_compare_means(
    data = plot_stat_data,
    aes(x = visit, y = perc_fm_change, group = asm_change_v4_group),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  scale_fill_manual(
    values = c("high" = manet_cols[5], "low/modest" = manet_cols[20]),
    name = "%ASM change at 1y"
  ) +
  scale_x_discrete() +
  labs(x = "Time", y = "Mean %FM change") +
  theme_minimal()
