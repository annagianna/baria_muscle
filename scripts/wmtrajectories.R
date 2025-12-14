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
baria_muscle <- read_rds("data/251214_BARIA_muscle_clinical.RDS")

# convert anthropometric data to long format, starting with body weight
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
  select(id, starts_with("fm_kg_v")) |> 
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
  select(id, contains("smm_kg_v"), -smm_kg_v0_quintile) |> 
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
  mutate(visit = as.integer(stringr::str_extract(visit, "\\d+$"))) # parse_number() extracts the "1" in hba1c

hba1c_percent_long <- baria_muscle |> 
  select(id, contains("hba1c_percent_v")) |> 
  pivot_longer(
    -id,
    names_to = "visit",
    values_to = "hba1c_percent"
  ) |> 
  mutate(visit = as.integer(stringr::str_extract(visit, "\\d+$"))) # same here

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

# clean up + calculate mean weight change from baseline
baria_muscle_trajectories <- add_time |>
  select(
    id, visit, time_y, age, sex, t2d_v0, height_cm, low_asm_v0, low_ffm_v0, low_smm_v0, low_smm_by_weight_v0, 
    weight_kg, ffm_kg, fm_kg, asm_kg, smm_kg, smm_by_weight, homa_ir, homa_b, crp_mgl, hba1c_mmolmol, hba1c_percent
  )|>
  arrange(visit) |> 
  filter(!is.na(low_smm_by_weight_v0)) |>  
  filter(time_y >= 0 & time_y <= 5.5) |> # there was a negative value at date v3 (probably incorrect data entry) and one over 900 for v2
  group_by(id) |> 
  mutate(
    delta_weight_kg = weight_kg - weight_kg[visit == 0],
    delta_ffm_kg = ffm_kg - ffm_kg[visit == 0],
    delta_fm_kg = fm_kg - fm_kg[visit == 0],
    delta_smm_kg = smm_kg - smm_kg[visit == 0],
    delta_asm_kg = asm_kg - asm_kg[visit == 0],
  ) |> 
  ungroup() |> 
   mutate(
    visit = case_when(
      visit == 0 ~ "Baseline",
      visit == 4 ~ "1y",
      visit == 5 ~ "2y",
      visit == 6 ~ "5y"
    )
  ) |> 
  print()

# summarize body composition parameters (means + SD) per visit & changes from baseline 
body_comp_summary <- baria_muscle_trajectories |> 
  filter(visit %in% c("Baseline", "1y", "2y", "5y")) |> # visits with available bia measurement
  group_by(visit, low_smm_by_weight_v0) |> 
  summarize(
    # n of available values per variable
    n_weight = sum(!is.na(weight_kg)),
    n_ffm    = sum(!is.na(ffm_kg)),
    n_fm     = sum(!is.na(fm_kg)),
    n_smm    = sum(!is.na(smm_kg)),
    n_asm    = sum(!is.na(asm_kg)),

    # weight
    mean_weight = mean(weight_kg, na.rm = TRUE),
    sd_weight = sd(weight_kg, na.rm = TRUE),
    mean_delta_weight = mean(delta_weight_kg, na.rm = TRUE),
    sd_delta_weight = sd(delta_weight_kg, na.rm = TRUE),

    # FFM
    mean_ffm = mean(ffm_kg, na.rm = TRUE),
    sd_ffm = sd(ffm_kg, na.rm = TRUE),
    mean_delta_ffm = mean(delta_ffm_kg, na.rm = TRUE),
    sd_delta_ffm = sd(delta_ffm_kg, na.rm = TRUE),

    # FM
    mean_fm = mean(fm_kg, na.rm = TRUE),
    sd_fm = sd(fm_kg, na.rm = TRUE),
    mean_delta_fm = mean(delta_fm_kg, na.rm = TRUE),
    sd_delta_fm = sd(delta_fm_kg, na.rm = TRUE),

    # SMM
    mean_smm = mean(smm_kg, na.rm = TRUE),
    sd_smm = sd(smm_kg, na.rm = TRUE),
    mean_delta_smm = mean(delta_smm_kg, na.rm = TRUE),
    sd_delta_smm = sd(delta_smm_kg, na.rm = TRUE),

    # ASM
    mean_asm = mean(asm_kg, na.rm = TRUE),
    sd_asm = sd(asm_kg, na.rm = TRUE),
    mean_delta_asm = mean(delta_asm_kg, na.rm = TRUE),
    sd_delta_asm = sd(delta_asm_kg, na.rm = TRUE)
  ) |> 
  ungroup()

# format as table
# means per visit
body_comp_wide_means <- body_comp_summary |> 
  mutate(
    weight = paste0(round(mean_weight, 1), " ± ", round(sd_weight, 1)),
    ffm = paste0(round(mean_ffm, 1), " ± ", round(sd_ffm, 1)),
    fm = paste0(round(mean_fm, 1), " ± ", round(sd_fm, 1)),
    smm = paste0(round(mean_smm, 1), " ± ", round(sd_smm, 1)),
    asm = paste0(round(mean_asm, 1), " ± ", round(sd_asm, 1)),
  ) |>
  select(visit, low_smm_by_weight_v0, weight, ffm, fm, smm, asm) |> 
  pivot_wider(
    id_cols = visit,
    names_from = low_smm_by_weight_v0,
    values_from = c(weight, ffm, fm, smm, asm),
    names_prefix = "low_"
  ) |> 
  print()

# for deltas
body_comp_wide_delta <- body_comp_summary |> 
  mutate(
    delta_weight = paste0(round(mean_delta_weight, 1), " ± ", round(sd_delta_weight, 1)),
    delta_ffm = paste0(round(mean_delta_ffm, 1), " ± ", round(sd_delta_ffm, 1)),
    delta_fm = paste0(round(mean_delta_fm, 1), " ± ", round(sd_delta_fm, 1)),
    delta_smm = paste0(round(mean_delta_smm, 1), " ± ", round(sd_delta_smm, 1)),
    delta_asm = paste0(round(mean_delta_asm, 1), " ± ", round(sd_delta_asm, 1))
  ) |> # format as table and export
  select(visit, low_smm_by_weight_v0, delta_weight, delta_ffm, delta_fm, delta_smm, delta_asm) |> 
  pivot_wider(
    id_cols = visit, # identifier
    names_from = low_smm_by_weight_v0,
    values_from = c(delta_weight, delta_ffm, delta_fm, delta_smm, delta_asm),
    names_prefix = "low_"
  ) |> 
  print()

# Compute statistics to add to table
# per visit
p_weight <- baria_muscle_trajectories |> 
  filter(visit %in% c("Baseline", "1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_weight = t.test(weight_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

p_ffm <- baria_muscle_trajectories |> 
  filter(visit %in% c("Baseline", "1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_ffm = t.test(ffm_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

p_fm <- baria_muscle_trajectories %>%
  filter(visit %in% c("Baseline", "1y", "2y", "5y")) %>%
  group_by(visit) %>%
  summarise(
    p_fm = t.test(fm_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

p_asm <- baria_muscle_trajectories |> 
  filter(visit %in% c("Baseline", "1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_asm = round(t.test(asm_kg ~ low_smm_by_weight_v0)$p.value, 5),
    .groups = "drop"
  )

p_smm <- baria_muscle_trajectories |> 
  filter(visit %in% c("Baseline", "1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_smm = t.test(smm_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

# add to wide means dataset
body_comp_wide_means_p <- body_comp_wide_means |>
  left_join(p_weight, by = "visit") |>
  left_join(p_ffm, by = "visit") |>
  left_join(p_fm, by = "visit") |>
  left_join(p_smm, by = "visit") |>
  left_join(p_asm, by = "visit")

# deltas
p_delta_weight <- baria_muscle_trajectories |> 
  filter(visit %in% c("1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_delta_weight = t.test(delta_weight_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

p_delta_ffm <- baria_muscle_trajectories |> 
  filter(visit %in% c("1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_delta_ffm = t.test(delta_ffm_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

p_delta_fm <- baria_muscle_trajectories %>%
  filter(visit %in% c("1y", "2y", "5y")) %>%
  group_by(visit) %>%
  summarise(
    p_delta_fm = t.test(delta_fm_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

p_delta_smm <- baria_muscle_trajectories |> 
  filter(visit %in% c("1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_delta_smm = t.test(smm_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

p_delta_asm <- baria_muscle_trajectories |> 
  filter(visit %in% c("1y", "2y", "5y")) |> 
  group_by(visit) |> 
  summarise(
    p_delta_asm = t.test(delta_asm_kg ~ low_smm_by_weight_v0)$p.value,
    .groups = "drop"
  )

# add to wide deltas 
body_comp_wide_delta_p <- body_comp_wide_delta |>
  left_join(p_delta_weight, by = "visit") |>
  left_join(p_delta_ffm, by = "visit") |>
  left_join(p_delta_fm, by = "visit") |>
  left_join(p_delta_smm, by = "visit") |>
  left_join(p_delta_asm, by = "visit")

# create gt table for means per visit
body_comp_gt_means_tbl <- body_comp_wide_means_p |>
  arrange(match(visit, c("Baseline", "1y", "2y", "5y"))) |> 
  relocate(p_weight, .after = weight_low_yes) |> 
  relocate(p_ffm, .after = ffm_low_yes) |> 
  relocate(p_fm, .after = fm_low_yes) |> 
  relocate(p_smm, .after = smm_low_yes) |> 
  gt() |>
  tab_header(title = "Body composition parameters per visit stratified by low muscle mass at baseline") |>
  
  # spanners
  tab_spanner(label = "Weight (kg)", columns = c(weight_low_yes, weight_low_no)) |>
  tab_spanner(label = "FFM (kg)", columns = c(ffm_low_yes, ffm_low_no)) |>
  tab_spanner(label = "FM (kg)", columns = c(fm_low_yes, fm_low_no)) |>
  tab_spanner(label = "SMM (kg)", columns = c(smm_low_yes, smm_low_no)) |>
  tab_spanner(label = "ASM (kg)", columns = c(asm_low_yes, asm_low_no)) |>
  
  # rename subcolumns to low and high mm
  cols_label(
    weight_low_yes = "Low SMM/W", weight_low_no = "High SMM/W",
    ffm_low_yes = "Low SMM/W", ffm_low_no = "High SMM/W", 
    fm_low_yes = "Low SMM/W", fm_low_no = "High SMM/W",
    smm_low_yes = "Low SMM/W", smm_low_no = "High SMM/W",
    asm_low_yes = "Low SMM/W", asm_low_no = "High SMM/W",
    p_weight = "p", p_ffm = "p", p_fm = "p", p_asm = "p", p_smm = "p",
  ) |> 

  # adjust font size for sub-columns (low and high) and for cell body
  tab_style(
    style = cell_text(size = px(10)),
    locations = cells_column_labels(columns = everything())
  ) |>
  tab_style(
    style = cell_text(size = px(10)),
    locations = cells_body(columns = everything())
  ) |> 
  fmt_number(
    columns = c(p_weight, p_ffm, p_fm, p_smm, p_asm),
    decimals = 3
  )

# create gt table for deltas
body_comp_gt_delta_tbl <- body_comp_wide_delta_p |>
  filter((visit != "Baseline")) |> 
  relocate(delta_weight_low_yes, .before = delta_weight_low_no) |> 
  relocate(delta_ffm_low_yes, .before = delta_ffm_low_no) |> 
  relocate(delta_fm_low_yes, .before = delta_fm_low_no) |>
  relocate(delta_smm_low_yes, .before = delta_smm_low_no) |>
  relocate(delta_asm_low_yes, .before = delta_asm_low_no) |>
  relocate(p_delta_weight, .after = delta_weight_low_no) |> 
  relocate(p_delta_ffm, .after = delta_ffm_low_no) |> 
  relocate(p_delta_fm, .after = delta_fm_low_no) |> 
  relocate(p_delta_smm, .after = delta_smm_low_no) |> 
  relocate(p_delta_asm, .after = delta_asm_low_no) |> 
  gt() |>
  tab_header(title = "Differences from baseline") |>
  
  # spanners
  tab_spanner(label = "ΔWeight (kg)", columns = c(delta_weight_low_yes, delta_weight_low_no)) |>
  tab_spanner(label = "ΔFFM (kg)", columns = c(delta_ffm_low_yes, delta_ffm_low_no)) |>
  tab_spanner(label = "ΔFM (kg)", columns = c(delta_fm_low_yes, delta_fm_low_no)) |>
  tab_spanner(label = "ΔSMM (kg)", columns = c(delta_smm_low_yes, delta_smm_low_no)) |>
  tab_spanner(label = "ΔASM (kg)", columns = c(delta_asm_low_yes, delta_asm_low_no)) |>
  
  # rename subcolumns to low and high mm
  cols_label(
    delta_weight_low_yes = "Low SMM/W", delta_weight_low_no = "High SMM/W",
    delta_ffm_low_yes = "Low SMM/W", delta_ffm_low_no = "High SMM/W", 
    delta_fm_low_yes = "Low SMM/W", delta_fm_low_no = "High SMM/W",
    delta_smm_low_yes = "Low SMM/W", delta_smm_low_no = "High SMM/W",
    delta_asm_low_yes = "Low SMM/W", delta_asm_low_no = "High SMM/W",
    p_delta_weight = "p", p_delta_ffm = "p", p_delta_fm = "p", p_delta_smm = "p", p_delta_asm = "p"
  ) |> 

  # adjust font size for sub-columns (low and high) and for cell body
  tab_style(
    style = cell_text(size = px(10)),
    locations = cells_column_labels(columns = everything())
  ) |>
  tab_style(
    style = cell_text(size = px(10)),
    locations = cells_body(columns = everything())
) |>
  fmt_number(
    columns = c(p_delta_weight, p_delta_ffm, p_delta_fm, p_delta_smm, p_delta_asm),
    decimals = 4
  )
body_comp_gt_delta_tbl

# Plot trajectories
# Themes
gb1_colors <- scale_color_manual(values = wes_palette("GrandBudapest1"))
gb2_colors <- scale_color_manual(values = wes_palette("GrandBudapest2"))
mr3_colors <- scale_color_manual(values = wes_palette("Moonrise3"))
isle1_colors <- scale_color_manual(values = wes_palette("IsleofDogs1"))

# Weight trajectories EDA
weight <- ggplot(baria_muscle_trajectories, aes(x = time_y, y = weight_kg, color = low_smm_by_weight_v0)) +
  geom_smooth(se = FALSE, show.legend = FALSE) +
  labs(title = "Body weight", color = "Low SMM/W", x = "Time (years)", y = "Weight (kg)")

ffm <- ggplot(baria_muscle_trajectories, aes(x = time_y, y = ffm_kg, color = low_smm_by_weight_v0)) +
  geom_smooth(se = FALSE) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = "Fat-free mass", color = "Low SMM/W", x = "Time (y)", y = "FFM (kg)")

fm <- ggplot(baria_muscle_trajectories, aes(x = time_y, y = fm_kg, color = low_smm_by_weight_v0)) +
  geom_smooth(se = FALSE) +
  labs(title = "Fat mass", color = "Low SMM/W", x = "Time (y)", y = "FM (kg)")

weight / (ffm | fm ) +
  plot_layout(
    guides = "collect"
  ) &
  isle1_colors &
  theme(legend.position = "bottom")

## better weight loss graphs (check the tirzepatide paper for reference)
