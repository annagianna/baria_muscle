# Weight and fat-free mass trajectories
# Anna Giannakogeorgou

# Definitions

# Packages
library(dplyr)
library(tidyverse)
library(lubridate)
library(ggpubr)
library(scales)

# Data

baria_muscle_clean <- read.csv("data/251208_BARIA_muscle_clinical.csv")
baria_muscle_clean <- baria_muscle_clean |> 
  mutate(
    across(where(is.character), as.factor),
    across(
      starts_with("date"), as.Date)
    )

View(baria_muscle_clean)

# Low muscle mass = lowest quintile on the recommended SMM/W for BIA by the ESPN/EASO consensus
# first calculate lowest quintile

baria_muscle_low <- baria_muscle_clean |> 
  filter( # filtering for device's (Maltron BioScan 920) resolution range excluding high and low outliers
    bia_resistance_50khz_v0 >= 100 & bia_resistance_50khz_v0 <= 1100,
    ffm_percent_v0 + fm_percent_v0 == 100,
    bia_resistance_50khz_v4 >= 100 & bia_resistance_50khz_v4 <= 1100,
    ffm_percent_v4 + fm_percent_v4 == 100,
    bia_resistance_50khz_v5 >= 100 & bia_resistance_50khz_v5 <= 1100,
    ffm_percent_v5 + fm_percent_v5 == 100,
    bia_resistance_50khz_v6 >= 100 & bia_resistance_50khz_v6 <= 1100,
    ffm_percent_v6 + fm_percent_v6 == 100,
  ) |> 
  mutate(
    smm_by_weight_v0 = smm_kg_v0 / weight_kg_v0,
    smm_by_weight_v4 = smm_kg_v4 / weight_kg_v4,
    smm_by_weight_v5 = smm_kg_v5 / weight_kg_v5,
    smm_by_weight_v6 = smm_kg_v6 / weight_kg_v6,
    id = as.factor(id)
  ) |> 
  group_by(sex) |> # different cut-offs for female and male participants
  mutate(
    # Calculate the quintiles
    ffm_kg_v0_quintile = quantile(ffm_kg_v0, probs = 0.20, na.rm = TRUE), 
    smm_kg_v0_quintile = quantile(smm_kg_v0, probs = 0.20, na.rm = TRUE), 
    asm_kg_v0_quintile = quantile(asm_kg_v0, probs = 0.20, na.rm = TRUE),
    smm_by_weight_v0_quintile = quantile(smm_by_weight_v0, probs = 0.20, na.rm = TRUE), # recommended by ESPEN/EASO for SO

    low_ffm_v0 = if_else(ffm_kg_v0 <= ffm_kg_v0_quintile, 1, 0),
    low_smm_v0 = if_else(smm_kg_v0 <= smm_kg_v0_quintile, 1, 0),
    low_asm_v0 = if_else(asm_kg_v0 <= asm_kg_v0_quintile, 1, 0),
    low_smm_by_weight_v0 = if_else(smm_by_weight_v0 <= smm_by_weight_v0_quintile, 1, 0),

    low_ffm_v0 = case_when(low_ffm_v0 == 1 ~ "yes", low_ffm_v0 == 0 ~ "no"),
    low_smm_v0 = case_when(low_smm_v0 == 1 ~ "yes", low_smm_v0 == 0 ~ "no"),
    low_asm_v0 = case_when(low_asm_v0 == 1 ~ "yes", low_asm_v0 == 0 ~ "no"), 
    low_smm_by_weight_v0 = case_when(low_smm_by_weight_v0 == 1 ~ "yes", low_smm_by_weight_v0 == 0 ~ "no"),
  ) |> 
  ungroup()

View(baria_muscle_low)

# eda shows greater delta weight in those with low_smm_by_weight_v0 at 1, 2 and 5 y

# Trajectories
# convert data to long format for plots

baria_muscle_long_weight <- baria_muscle_low |> 
  pivot_longer(
    c(weight_kg_v0, weight_kg_v2, weight_kg_v3, weight_kg_v4, weight_kg_v5, weight_kg_v6),
    names_to = "visit",
    values_to = "weight_kg"
  ) |> 
  mutate(visit = as.integer(parse_number(visit)))

# create time column based on difftime for trajectories

baria_muscle_time <- baria_muscle_low |> 
  select(id, v0_to_v2_numeric_y, v0_to_v3_numeric_y, v0_to_v4_numeric_y, v0_to_v5_numeric_y, v0_to_v6_numeric_y, v0_to_v7_numeric_y) |> 
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
baria_muscle_time |> 
  group_by(visit) |> 
  summarize(mean = mean(time_y, na.rm = TRUE))

# Pivot longer ffm 

# Pivot longer fm

# Pivot longer smm

# Pivot longer asm


# merge with baria_muscle_long

baria_muscle_trajectories <- left_join(baria_muscle_long_weight, baria_muscle_time, by = c("id", "visit"))
View(baria_muscle_overtime)

#######

# Weight trajectories

ggplot(baria_muscle_overtime, aes(x = time_y, y = weight_kg, color = low_asm_v0)) +
  geom_smooth(se = FALSE) +
  scale_color_brewer(palette = "Paired")

# FFM trajectories

# FM trajectories

# Composition plot

