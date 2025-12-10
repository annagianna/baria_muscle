# Weight and fat-free mass trajectories
# Anna Giannakogeorgou

# Packages
library(dplyr)
library(tidyverse)
library(lubridate)
library(ggpubr)
library(scales)
library(MetBrewer)
library(wesanderson)
library(patchwork)

# Data
baria_muscle <- read_rds("data/251209_BARIA_muscle_clinical.RDS")
View(baria_muscle)

# Trajectories
# convert data to long format for plots, starting with body weight
weight_long <- baria_muscle |> 
  pivot_longer(
    contains("weight_kg_v"),
    names_to = "visit",
    values_to = "weight_kg"
  ) |> 
  relocate(visit, weight_kg, .after = id) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot longer ffm 
ffm_long <- baria_muscle |> 
  select(id, contains("ffm_kg_v"), -ffm_kg_v0_quintile, -low_ffm_v0) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "ffm_kg"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot longer fm
fm_long <- baria_muscle |> 
  select(id, contains("fm_kg_v")) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "fm_kg"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot longer asm
asm_long <- baria_muscle |> 
  select(id, contains("asm_kg_v"), -asm_kg_v0_quintile, -low_asm_v0) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "asm_kg"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot longer smm
smm_long <- baria_muscle |> 
  select(id, contains("smm_kg_v"), -smm_by_weight_v0_quintile, -low_smm_v0) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "smm_kg"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot longer smm/w
smm_by_weight_long <- baria_muscle |> 
  select(id, contains("smm_by_weight_v"), -smm_by_weight_v0_quintile, -low_smm_by_weight_v0) |>
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "smm_by_weight"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot age longer
age_long <- baria_muscle |> 
  select(id, contains("age_v")) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "age"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot HOMA-IR longer
homa_ir_long <- baria_muscle |> 
  select(id, contains("homa_ir_v")) |> 
  pivot_longer(
    c(homa_ir_v0, homa_ir_v4, homa_ir_v5),
    names_to = "visit",
    values_to = "homa_ir"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot HOMA-B longer
homa_b_long <- baria_muscle |> 
  select(id, contains("homa_b")) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "homa_b"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot hba1c longer (both units)
hba1c_mmolmol_long <- baria_muscle |> 
  select(id, contains("hba1c_mmolmol_v")) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "hba1c_mmolmol"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

hba1c_percent_long <- baria_muscle |> 
  select(id, contains("hba1c_percent_v")) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "hba1c_percent"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# pivot longer crp
crp_long <- baria_muscle |> 
  select(id, contains("crp_mgl_v")) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "crp_mgl"
  ) |> 
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
    v7 = v0_to_v7_numeric_y
  ) |>
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
# weight and ffm
add_ffm <- left_join(weight_long, ffm_long, by = c("id", "visit"))

# add fm
add_fm <- left_join(add_ffm, fm_long, by = c("id", "visit"))

# add asm
add_asm <- left_join(add_fm, asm_long, by = c("id", "visit"))

# add smm
add_smm <- left_join(add_asm, smm_long, by = c("id", "visit"))

# add smm/weight
add_smm_by_weight <- left_join(add_smm, smm_by_weight_long, by = c("id", "visit"))

# add age
add_age <- left_join(add_smm_by_weight, age_long, by = c("id", "visit"))

# add HOMA-IR
add_homa_ir <- left_join(add_age, homa_ir_long, by = c("id", "visit"))

# add HOMA-B
add_homa_b <- left_join(add_homa_ir, homa_b_long, by = c("id", "visit"))

# add hba1c (mmol/mol and %)
add_hba1c_mmolmol <- left_join(add_homa_b, hba1c_mmolmol_long, by = c("id", "visit"))
add_hba1c_percent <- left_join(add_hba1c_mmolmol, hba1c_percent_long, by = c("id", "visit"))

# add crp
add_crp_mgl <- left_join(add_hba1c_percent, crp_long, by = c("id", "visit"))

# add time
add_time <- left_join(add_crp_mgl, time, by = c("id", "visit"))

# clean up
baria_muscle_trajectories <- add_time |>
  select(
    id, visit, time_y, age, sex, t2d_v0, height_cm, low_asm_v0, low_ffm_v0, low_smm_v0, low_smm_by_weight_v0, 
    weight_kg, ffm_kg, fm_kg, asm_kg, smm_kg, smm_by_weight, homa_ir, homa_b, crp_mgl, hba1c_mmolmol, hba1c_percent
  )|>
  arrange(visit) |> 
  filter(!is.na(low_smm_by_weight_v0)) |>  
  filter(time_y >= 0 & time_y <= 5.5) |> # there was a negative value at date v3 (probably incorrect data entry) and one over 900 for v2
  print()

View(baria_muscle_trajectories)

# Trajectories
# Theme
gb1_colors <- scale_color_manual(values = wes_palette("GrandBudapest1"))
gb2_colors <- scale_color_manual(values = wes_palette("GrandBudapest2"))
mr3_colors <- scale_color_manual(values = wes_palette("Moonrise3"))
isle1_colors <- scale_color_manual(values = wes_palette("IsleofDogs1"))

# Weight trajectories
weight <- ggplot(baria_muscle_trajectories, aes(x = time_y, y = weight_kg, color = low_smm_by_weight_v0)) +
  geom_smooth(se = FALSE, show.legend = FALSE) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = "Body weight", color = "Low SMM/W", x = "Time (years)", y = "Weight (kg)")

ffm <- ggplot(baria_muscle_trajectories, aes(x = time_y, y = ffm_kg, color = low_smm_by_weight_v0)) +
  geom_smooth(se = FALSE) +
  labs(title = "FFM", color = "Low SMM/W", x = "Time (y)", y = "FFM (kg)")

fm <- ggplot(baria_muscle_trajectories, aes(x = time_y, y = fm_kg, color = low_smm_by_weight_v0)) +
  geom_smooth(se = FALSE) +
  labs(title = "FM", color = "Low SMM/W", x = "Time (y)", y = "FM (kg)")

weight / (ffm | fm ) +
  plot_layout(
    guides = "collect"
  ) &
  gb2_colors &
  theme(legend.position = "bottom")

# FFM trajectories

# FM trajectories

# Composition plot

