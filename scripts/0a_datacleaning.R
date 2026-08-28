# Baria project: Muscle mass trajectories and gut microbiota following bariatric surgery - Data cleaning
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)
library(phyloseq)

# Open Data and see properties
baria_clinical_data_raw <- readRDS("./data/raw_data/BARIA.clinical.2024-12-09.723.2043.RDS")
baria_mb <- readRDS("data/raw_data/ps.BARIA.metaphlan.706.2548.RDS")

#### Clinical Data ####
# Longitudinal vars
mmt_timepoints <- c(0, 10, 20, 30, 60, 90, 120)

long_vars <- c(
  "date",

  # Body composition
  "bmi", "weight_kg", "wc_cm", "fm_kg", "fm_percent", "ffm_kg", "ffm_percent", "bia_resistance_50khz", "upperleg_cm",

  # Cardiometabolic/lab parameters
  str_c("glucose_mmoll_mmt_", mmt_timepoints), str_c("insulin_pmoll_mmt_", mmt_timepoints), str_c("cpeptide_nmoll_mmt_", mmt_timepoints),
  "glucagon_ngl_mmt_0", "hba1c", "hba1c_mmolmol", "crp_mgl",
  "systolic_bp_mmhg", "diastolic_bp_mmhg",
  "gammagt_ul", "asat_ul", "alat_ul", "triglycerides_mmoll", "tsh_miul", "ft4_pmoll",

  # Medication
  "medication_list", "medication_binary",

  # Nexfin
  "nexfin_hr", "nexfin_dpdt", "nexfin_sv", "nexfin_svi", "nexfin_co", "nexfin_ci", "nexfin_svr", "nexfin_svri"
)

long_vars_pattern <- str_c("^(", str_c(long_vars, collapse = "|"), ")_(v\\d+)$")

# Clinical data
baria_muscle_vars <- baria_clinical_data_raw |>
  select(
    # Baseline & static vars
    id = Subject_ID, date_v0 = date, sg_type = type_surgery, age_v0 = Age, sex, t2d_v0 = dm, 
    bmi_v0 = bmi, weight_kg_v0 = weight, height_cm = height, wc_cm_v0 = taille, upperleg_cm_v0 = upperleg,
    fm_kg_v0 = tbf, fm_percent_v0 = tbf_percent, ffm_kg_v0 = ffm, ffm_percent_v0 = ffm_percent, bia_resistance_50khz_v0 = rawdata_50Khz_Resistance,
    medication_binary_v0 = meds, medication_list_v0 = medication_v0_freetext, sport_v0 = scre_sport,
    systolic_bp_mmhg_v0 = systolic_pressure_v0, diastolic_bp_mmhg_v0 = diastolic_pressure_v0, aht = hypertension,
    glucose_mmoll_v0 = glucose, glucose_mmoll_mmt_0_v0 = min0gluc, insulin_pmoll_mmt_0_v0 = min0insulin,
    glucagon_ngl_mmt_0_v0 = min0glucagon, cpeptide_nmoll_mmt_0_v0 = min0cpept, hba1c_v0 = hba1c, hba1c_mmolmol_v0 = hba1c__IFCC_mmolmol,
    matches("^min(10|20|30|60|90|120)(gluc|insulin|cpept)$"), # MMT vars
    gammagt_ul_v0 = ggt, alat_ul_v0 = alat, asat_ul_v0 = asat, triglycerides_mmoll_v0 = triglycerides,
    crp_mgl_v0 = crp, tsh_miul_v0 = tsh, ft4_pmoll_v0 = ft4,

    # v2-v5 (Repeated vars)
    matches("^V[2-5]_(date|bmi|weight|taille|tbf|tbf_percent|ffm|ffm_percent|rawdata_50Khz_Resistance|upperleg|min0gluc|min0insulin|min0glucagon|min0cpept|hba1c|hba1c__IFCC_mmolmol|crp|ggt|asat|alat|triglycerides|tsh|ft4)$"),
    matches("^V[2-5]_min(10|20|30|60|90|120)(gluc|insulin|cpept)$"),
    matches("^(systolic|diastolic)_pressure_v[45]$"),


    # v4-v5 medication free text
    medication_list_v4 = Medication_v4_freetext,
    medication_list_v5 = medication_v5_freetext,

    # v6 (5 years)
    date_v6 = V6_date, bmi_v6 = V6_bmi_1, weight_kg_v6 = V6_weight_1, wc_cm_v6 = V6_taille_1, fm_kg_v6 = V6_tbf_1, fm_percent_v6 = V6_tbf_percent_1,
    ffm_kg_v6 = V6_ffm_1, ffm_percent_v6 = V6_ffm_percent_1, bia_resistance_50khz_v6 = V6_rawdata_1_50Khz_Resistance, upperleg_cm_v6 = V6_upperleg_1,
    hba1c_v6 = V6_hba1c_1, hba1c_mmolmol_v6 = V6_hba1c_IFCC_mmolmol, crp_mgl_v6 = V6_crp_1,

    # v7 (10 years)
    date_v7 = V7_date, bmi_v7 = V7_BMI, weight_kg_v7 = V7_weight, wc_cm_v7 = V7_taille, upperleg_cm_v7 = V7_upperleg, hba1c_v7 = V7_hba1c,
    hba1c_mmolmol_v7 = V7_hba1c_IFCC_mmolmol, crp_mgl_v7 = V7_CRP,

    # Nexfin data
    matches("^nexfin_(hr|dpdt|sv|svi|co|ci|svr|svri)_", ignore.case = TRUE),
    -nexfin_HR_v0
  ) |>
  rename_with(~ str_replace(.x, "^V([2-5])_(.+)$", "\\2_v\\1"), matches("^V[2-5]_")) |> 
  rename_with(~ str_c(.x, "_v0"), matches("^min(10|20|30|60|90|120)(gluc|insulin|cpept)$")) |>
  rename_with(
    ~ .x |>
      str_replace("^weight_", "weight_kg_") |>
      str_replace("^taille_", "wc_cm_") |>
      str_replace("^tbf_percent_", "fm_percent_") |>
      str_replace("^tbf_", "fm_kg_") |>
      str_replace("^ffm_v", "ffm_kg_v") |>
      str_replace("^rawdata_50Khz_Resistance_", "bia_resistance_50khz_") |>
      str_replace("^upperleg_", "upperleg_cm_") |>
      str_replace("^min(0|10|20|30|60|90|120)gluc_", "glucose_mmoll_mmt_\\1_") |>
      str_replace("^min(0|10|20|30|60|90|120)insulin_", "insulin_pmoll_mmt_\\1_") |>
      str_replace("^min(0|10|20|30|60|90|120)cpept_", "cpeptide_nmoll_mmt_\\1_") |>
      str_replace("^min0glucagon_", "glucagon_ngl_mmt_0_") |>
      str_replace("^hba1c__IFCC_mmolmol_", "hba1c_mmolmol_") |>
      str_replace("^crp_", "crp_mgl_") |>
      str_replace("^ggt_", "gammagt_ul_") |>
      str_replace("^asat_", "asat_ul_") |>
      str_replace("^alat_", "alat_ul_") |>
      str_replace("^triglycerides_", "triglycerides_mmoll_") |>
      str_replace("^tsh_", "tsh_miul_") |>
      str_replace("^ft4_", "ft4_pmoll_") |>
      str_replace("^systolic_pressure_", "systolic_bp_mmhg_") |>
      str_replace("^diastolic_pressure_", "diastolic_bp_mmhg_"),
    matches("_v[2-5]$")
) |> 
  rename_with(
    ~ .x |> 
      str_to_lower() |> # inconsistencies in Nexfin data nomenclature
      str_replace("_v6_1$", "_v6"),
    matches("^nexfin_(hr|dpdt|sv|svi|co|ci|svr|svri)_", ignore.case = TRUE)
  ) |> 
  mutate(
    across(everything(), ~ replace(.x, .x %in% c(-99, -98, -97), NA)), # these values are NA for different reasons in the Baria dataset
    across(.cols = starts_with("date"), ~ if_else(.x %in% c("01-01-2999", "01-01-2997", "01-01-2995"), NA_character_, .x)),
    sex = case_when(sex == "1" ~ "male", sex == "2" ~ "female"),
    t2d_v0 = case_when(t2d_v0 == "1" ~ "yes", t2d_v0 == "2" ~ "no"),
    aht = case_when(aht == "1" ~ "yes", aht == "2" ~ "no"),
    medication_binary_v0 = case_when(medication_binary_v0 == "1" ~ "yes", medication_binary_v0 == "2" ~ "no"),
    sport_v0 = case_when(sport_v0 == "1" ~ "yes", sport_v0 == "2" ~ "no"),
    sg_type = case_when(sg_type == "1" ~ "rygb", sg_type == "2" ~ "omegaloop", sg_type == "3" ~ "sg")
  ) |>
  pivot_longer(
    cols = matches(long_vars_pattern),
    names_to = c(".value", "visit"), 
    names_pattern = "(.+)_(v\\d+)$"
  ) |> 
  mutate(# manual corrections for incorrect BIA data entries
    id = as.character(id),
    fm_kg = case_when(
      id == "212" & visit == "v0" ~ fm_kg / 100,
      id == "600" & visit == "v0" ~ fm_kg * 10,
      id == "30004" & visit == "v4" ~ fm_kg / 100,
      id == "567" & visit == "v4" ~ fm_kg / 10,
      id == "335" & visit == "v4" ~ fm_kg * 10,
      id == "612" & visit == "v0" ~ 52.0,
      id == "23" & visit == "v5" ~ 13.8,
      TRUE ~ fm_kg
    ),
    ffm_kg = case_when(
      id == "238" & visit == "v0" ~ ffm_kg * 10,
      id == "30016" & visit == "v0" ~ ffm_kg * 10,
      id == "50" & visit == "v5" ~ 79.16,
      id == "484" & visit == "v4" ~ 51.3,
      id == "520" & visit == "v4" ~ 57.5,
      id == "451" & visit == "v4" ~ 39.7,
      TRUE ~ ffm_kg
    ),
    fm_percent = case_when(
      id == "551" & visit == "v0" ~ fm_percent / 100,
      id == "572" & visit == "v4" ~ fm_percent / 100,
      id == "120" & visit == "v5" ~ fm_percent / 100,
      TRUE ~ fm_percent
    ),
    ffm_percent = case_when(
      id == "590" & visit == "v0" ~ ffm_percent * 10,
      id == "95"  & visit == "v5" ~ ffm_percent / 100,
      id == "130" & visit == "v5" ~ ffm_percent / 100,
      id == "496" & visit == "v5" ~ ffm_percent / 100,
      id == "364" & visit == "v5" ~ ffm_percent / 10,
      id == "217" & visit == "v5" ~ 67.0,
      id == "257" & visit == "v0" ~ 54.4,
      TRUE ~ ffm_percent
    ),
    # BIA quality & consistency checks
    bia_perc_diff = abs((ffm_percent + fm_percent) - 100),
    bia_kg_diff = abs((ffm_kg + fm_kg) - weight_kg),
    bia_resistance_valid = case_when(
      is.na(bia_resistance_50khz) ~ NA,
      between(bia_resistance_50khz, 110, 1000) ~ TRUE, # Device's resistance range (Maltron BioScan 920)
      TRUE ~ FALSE
    ),
    bia_valid =
      !is.na(ffm_kg) &
      !is.na(fm_kg) &
      !is.na(weight_kg) &
      !is.na(bia_perc_diff) &
      !is.na(bia_kg_diff) &
      bia_perc_diff <= 5 &
      bia_kg_diff <= 5 &
      bia_resistance_valid %in% TRUE &
      !(id == "30023" & visit == "v0") &
      !(id == "30036" & visit == "v0") &
      !(id == "494" & visit == "v4") &
      !(id == "73" & visit == "v4") 
)

## Medication textbox cleaning
baria_muscle_vars_meds <- baria_muscle_vars |>
  mutate(
    medication_list = str_to_lower(medication_list),
    medication_list = case_when(medication_list %in% c("geen", "geen medicatie", "none") ~ "none", TRUE ~ medication_list)
  )

baria_muscle_vars_meds_notypos <- baria_muscle_vars_meds |>
  mutate(medication_list = medication_list |> 
    str_replace_all("\\bamitriptiline\\b", "amitriptyline") |>
    str_replace_all("\\bduloxatine\\b", "duloxetine") |>
    str_replace_all("\\bariprazol\\b", "aripiprazol") |>
    str_replace_all("\\bnifidipine\\b", "nifedipine") |>
    str_replace_all("\\bevothyroxine\\b", "levothyroxine") |>
    str_replace_all("\\bglicliazide\\b", "gliclazide") |>
    str_replace_all("\\bglicazide\\b", "gliclazide") |>
    str_replace_all("\\bamlopdipine\\b", "amlodipine") |>
    str_replace_all("\\bzitalopram\\b", "citalopram") |>
    str_replace_all("\\bperinopril\\b", "perindopril") |>
    str_replace_all("\\bvenlaflaxine\\b", "venlafaxine") |>
    str_replace_all("\\bglimerpiride\\b", "glimepiride") |>
    str_replace_all("\\bhydrchloorthiazide\\b", "hct") |>
    str_replace_all("\\bhydrochloorthiazide\\b", "hct") |>
    str_replace_all("\\bpnatoprazol\\b", "pantoprazole") |>
    str_replace_all("\\bomeprazol\\b", "omeprazole") |>
    str_replace_all("\\besomeprazol\\b", "esomeprazole") |>
    str_replace_all("\\bparvastatine\\b", "pravastatine") |>
    str_replace_all("\\bcandarsetan\\b", "candesartan") |>
    str_replace_all("\\bezatimibe\\b", "ezetimibe")
  )

# Define medication categories
# Diabetes medication
dm_meds <- list(
  metformin = c("metformine", "glucophage"),
  sus = c("gliclazide", "diamicron", "glibenclamide", "daonil", "glimepiride", "amaryl", "glipizide", "minodiab", "tolbutamide"), # Sulfonylureas
  dpp4is = c("sitagliptine", "januvia", "vildagliptine", "galvus", "saxagliptine", "onglyza", "linagliptine", "trajenta", "alogliptine", "vipidia"), # DPP-4 inhibitors
  glp1ras = c("liraglutide", "victoza", "semaglutide", "ozempic", "rybelsus", 
  "exenatide", "byetta", "bydureon", "dulaglutide", "trulicity", "lyxumia", "lixisenatide"), # GLP-1RAs
  sglt2is = c("dapagliflozine", "forxiga", "empagliflozine", "jardiance",
  "canagliflozine", "invokana", "ertugliflozine", "steglatro", "steeglatro"), # SGLT2-inhibitors
  tzds = c("pioglitazon", "pioglitazone"), # thiazolidinediones
  insulin = c("insuline", "lantus", "levemir", "novorapid", "apidra", "toujeo", "tresiba", "degludec",
  "humalog", "novomix", "fiasp", "actrapid", "isofaan", "insulatard", "glargine", "aspart") # insulins (short-, medium- and long-acting)
)

# Create regex patterns for each class of diabetes medication
dm_meds_patterns <- lapply(dm_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"), ")\\b")
})
dm_meds_any_pattern <- str_c("\\b(", str_c((unlist(dm_meds)), collapse = "|"), ")\\b")

# Antihypertensive medication
aht_meds <- list(
  ace_inhibitors = c("perindopril", "candesartan", "enalapril", "ramipril", "lisinopril", "captopril", "fosinopril", "quinapril"),
  arbs = c("losartan", "valsartan", "olmesartan", "irbesartan", "telmisartan", "eprosartan"), # Angiotensin receptor blockers  
  ccbs = c("amlodipine", "lercandipine", "lerdip", "barnidipine", "cyress", "nifedipine", "adalat", "verapamil",
  "felodipine", "nicardipine", "isradipine", "diltiazem"), # Calcium channel blockers
  bblockers = c("metoprolol", "metoprololtart", "metoprololsucc", "selokeen", "carvedilol", "labetalol", "bisoprolol",
  "atenolol", "propanonol", "nebivolol", "solatol"), # Beta-blockers
  a2_agonists = c("clonidine", "moxonidine", "methyldopa"), # Central a2 agonists
  diuretics = c("furosemide", "hct", "hydrochloorthiazide", "spironolacton", "triamteren", "amiloride",
  "dytenzide", "bumetadine", "chloortalidon", "indapamide", "eplerenone"), # diuretics
  combi_aht_meds = c("losartan/hydrochloorthiazide", "lodoz", "preterax", "moduretic") # Combination/others
)

aht_meds_patterns <- lapply(aht_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"),")\\b")
})
aht_meds_any_pattern <- str_c("\\b(", str_c(unlist(aht_meds), collapse = "|"),")\\b")

# Lipid lowering
lipidlowering_meds <- list(
  statins = c("simvastatine", "atorvastatine", "pravastatine", "rosuvastatine", "fluvastatine", 
  "crestor", "lipitor", "zocor", "pravachol", "lescol"), # Statins (HMG-CoA reductase inhibitors)
  ezetimibe = c("ezetimib", "ezetrol"), # Cholesterol absorption inhibitors
  pcsk9is = c("alirocumab", "praluent", "evolocumab", "repatha"), # PCSK9 inhibitors
  inclisiran = c("inclisiran", "leqvio"), # siRNA against PCSK9
  bempedoic_acid = c("bempedoïnezuur", "nexletol", "nustendi"), # ATP citrate lyase inhibitors
  bile_acid_sequestrants = c("colestyramine", "questran"), # Bile acid sequestrants
  fibrates = c("bezalip", "fibraat", "fenofibraat", "lipanthyl", "tricor", "gemfibrozil", "lopid"), # Fibrates (PPAR-α agonists)
  omega3 = c("omega-3-vetzuren", "omacor", "epanova"), # Omega-3 fatty acids
  niacin = c("nicotinezuur", "niaspan") # Nicotinic acid
)

# create regex patterns for each class of lipid-lowering medication
lipidlowering_meds_patterns <- lapply(lipidlowering_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"), ")\\b")
})
lipidlowering_meds_any_pattern <- str_c("\\b(", str_c((unlist(lipidlowering_meds)), collapse = "|"), ")\\b")

# Thyroid medication
thyroid_meds <- c("levothyroxine", "thyrax", "euthyrox")
thyroid_meds_pattern <- str_c("\\b(", str_c(thyroid_meds, collapse = "|"),")\\b")

# Psychiatric medication (that may cause weight loss/gain)
psychiatric_meds <- list(
  ssris = c("escitalopram", "citalopram", "cipramil", "sertraline", "paroxetine", "fluoxetine", "fluvoxamine"), # SSRIS
  tcas = c("amitriptyline", "nortriptyline", "clomipramine"), # Tricyclic antidepressants
  snris = c("venlafaxine", "efexor", "duloxetine"), # SNRIs
  ndris = c("bupropion"), # NDRIs
  antipsychotics = c("aripiprazol", "quetiapine", "olanzapine"), # Atypical antipsychotics
  moods = c("lamotrigine", "pregabalin"), # Mood stabilizers (anticonvulsants)
  adhd_meds = c("methylfenidaat", "concerta", "elvanse"), # ADHD medication/stimulants
  hypnotics = c("zopiclone", "zolpidem") # Sleep medication / hypnotics
)

psychiatric_meds_patterns <- lapply(psychiatric_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"),")\\b")
})
psychiatric_meds_any_pattern <- str_c("\\b(", str_c(unlist(psychiatric_meds), collapse = "|"), ")\\b")

# PPIs (effects on gut microbiota)
ppi <- c("omeprazole", "pantoprazole", "esomeprazole")
ppi_pattern <- str_c("\\b(", str_c(ppi, collapse = "|"), ")\\b")

# Antibiotics
antibiotics <- c(
  #Penicillins
  "amoxicilline", "amoxicillin", "augmentin", "co-amoxiclav", "flucloxacilline", "flucloxacillin", "penicilline", "penicillin",
  # Tetracyclines
  "doxycycline", "minocycline",
  # Macrolides
  "azitromycine", "azithromycin", "claritromycine", "clarithromycin", "erytromycine", "erythromycin",
  # Fluoroquinolones
  "ciprofloxacine", "ciprofloxacin", "levofloxacine", "levofloxacin", "moxifloxacine", "moxifloxacin",
  # Lincosamides
  "clindamycine", "clindamycin",
  # Nitroimidazoles
  "metronidazol", "metronidazole",
  # Urinary tract antibiotics
  "nitrofurantoine", "nitrofurantoin", "fosfomycine", "fosfomycin",
  # Sulfonamides
  "co-trimoxazol", "cotrimoxazol", "trimethoprim", "sulfamethoxazol", "sulfamethoxazole"
)
antibiotics_pattern <- str_c("\\b(", str_c(antibiotics, collapse = "|"), ")\\b")

# Medication class functions
medication_classes <- list(
  metformin = ~ if_else(str_detect(.x, dm_meds_patterns$metformin), "yes", "no"),
  sus = ~ if_else(str_detect(.x, dm_meds_patterns$sus), "yes", "no"),
  dpp4is = ~ if_else(str_detect(.x, dm_meds_patterns$dpp4is), "yes", "no"),
  glp1ras = ~ if_else(str_detect(.x, dm_meds_patterns$glp1ras), "yes", "no"),
  sglt2is = ~ if_else(str_detect(.x, dm_meds_patterns$sglt2is), "yes", "no"),
  tzds = ~ if_else(str_detect(.x, dm_meds_patterns$tzds), "yes", "no"),
  insulin = ~ if_else(str_detect(.x, dm_meds_patterns$insulin), "yes", "no"),
  dm_meds = ~ if_else(str_detect(.x, dm_meds_any_pattern), "yes", "no"),
  ace_inhibitors = ~ if_else(str_detect(.x, aht_meds_patterns$ace_inhibitors), "yes", "no"),
  arbs = ~ if_else(str_detect(.x, aht_meds_patterns$arbs), "yes", "no"),
  ccbs = ~ if_else(str_detect(.x, aht_meds_patterns$ccbs), "yes", "no"),
  bblockers = ~ if_else(str_detect(.x, aht_meds_patterns$bblockers), "yes", "no"),
  a2_agonists = ~ if_else(str_detect(.x, aht_meds_patterns$a2_agonists), "yes", "no"),
  diuretics = ~ if_else(str_detect(.x, aht_meds_patterns$diuretics), "yes", "no"),
  combi_aht_meds = ~ if_else(str_detect(.x, aht_meds_patterns$combi_aht_meds), "yes", "no"),
  aht_meds = ~ if_else(str_detect(.x, aht_meds_any_pattern), "yes", "no"),
  statins = ~ if_else(str_detect(.x, lipidlowering_meds_patterns$statins), "yes", "no"),
  lipidlowering_meds = ~ if_else(str_detect(.x, lipidlowering_meds_any_pattern), "yes", "no"),
  thyroid_meds = ~ if_else(str_detect(.x, thyroid_meds_pattern), "yes", "no"),
  ssris = ~ if_else(str_detect(.x, psychiatric_meds_patterns$ssris), "yes", "no"),
  tcas = ~ if_else(str_detect(.x, psychiatric_meds_patterns$tcas), "yes", "no"),
  snris = ~ if_else(str_detect(.x, psychiatric_meds_patterns$snris), "yes", "no"),
  ndris = ~ if_else(str_detect(.x, psychiatric_meds_patterns$ndris), "yes", "no"),
  antipsychotics = ~ if_else(str_detect(.x, psychiatric_meds_patterns$antipsychotics), "yes", "no"),
  moods = ~ if_else(str_detect(.x, psychiatric_meds_patterns$moods), "yes", "no"),
  adhd_meds = ~ if_else(str_detect(.x, psychiatric_meds_patterns$adhd_meds), "yes", "no"),
  hypnotics = ~ if_else(str_detect(.x, psychiatric_meds_patterns$hypnotics), "yes", "no"),
  psychiatric_meds = ~ if_else(str_detect(.x, psychiatric_meds_any_pattern), "yes", "no"),
  ppi = ~ if_else(str_detect(.x, ppi_pattern), "yes", "no"),
  abx = ~ if_else(str_detect(.x, antibiotics_pattern), "yes", "no")
)

# Create medication categories and subcategories
baria_muscle_vars_meds <- baria_muscle_vars_meds_notypos |>
  mutate(
    across(medication_list, medication_classes, .names = "{.fn}"),
    medication_binary = case_when(visit == "v0" ~ medication_binary, medication_list == "none" ~ "no", !is.na(medication_list) ~ "yes", TRUE ~ NA_character_)
  )

### Define final analysis cohort: valid baseline BIA & available baseline shotgun data & not taking antibiotics
# 1. Valid baseline BIA & 2. no antibiotics
bia_abx_v0_ids <- baria_muscle_vars_meds |>
  filter(
    # Baseline BIA QC
    visit == "v0",
    bia_valid,

    # No antibiotics
    abx == "no"
  ) |>
  pull(id) |> 
  unique()

# 3. Available & valid baseline microbiome/shotgun data
# Keep single runs and first run of duplicated samples (after comparing 1st and 2nd runs)
run1_mb <- prune_samples(
  (sample_data(baria_mb)$Extra_data == "NA" | sample_data(baria_mb)$Extra_data == "rep1") & !str_detect(sample_data(baria_mb)$Time_Point, "^V\\d+re$"),
  baria_mb
)

mb_v0_ids <- sample_data(run1_mb) |> 
  data.frame() |> 
  filter(Time_Point == "V-1") |> 
  pull(Subject_ID) |> 
  as.character() |>
  unique()

bia_abx_mb_ids <- intersect(bia_abx_v0_ids, mb_v0_ids)

## Long dataset ##
baria_muscle_long_all <- baria_muscle_vars_meds |> 
  filter(id %in% bia_abx_mb_ids) |> 
  mutate(date = dmy(date)) |> 
  group_by(id) |> 
  mutate(
    date_baseline = date[visit == "v0"],
    n_years_from_v0 = as.numeric(date - date_baseline) / 365.25,
    age = age_v0 + n_years_from_v0
  ) |>
  mutate(
    hba1c_percent = if_else(hba1c < 15, hba1c, hba1c * 0.0915 + 2.15),
    hba1c_mmolmol = if_else(is.na(hba1c_mmolmol) == FALSE, hba1c_mmolmol, 10.93 * hba1c_percent - 23.5),

    # HOMA-IR & HOMA-2B (insulin unit conversion from pmol/l to uU/ml)
    homa_ir = (insulin_pmoll_mmt_0 / 6.945) * glucose_mmoll_mmt_0 / 22.5,
    homa_b = (20 * (insulin_pmoll_mmt_0 / 6.945)) / (glucose_mmoll_mmt_0 - 3.5),

    # T2D incidence based on lab values at follow-up
    t2d_labs = case_when(
      is.na(hba1c_percent) & is.na(glucose_mmoll_mmt_0) ~ NA_character_,
      hba1c_percent >= 6.5 | glucose_mmoll_mmt_0 >= 7.0 ~ "yes",
      TRUE ~ "no"
    ),

    # Prediabetes based on lab values
    prediab_labs = case_when(
      t2d_labs == "yes" ~ "no",
      is.na(hba1c_percent) & is.na(glucose_mmoll_mmt_0) ~ NA_character_,
      (hba1c_percent >= 5.7 & hba1c_percent <= 6.4) | (glucose_mmoll_mmt_0 >= 5.6 & glucose_mmoll_mmt_0 <= 6.9) ~ "yes",
      TRUE ~ "no"
    ),
    ffmi = if_else(bia_valid, ffm_kg / ((height_cm / 100)^2), NA_real_),
    fmi = if_else(bia_valid, fm_kg / ((height_cm / 100)^2), NA_real_),
    smm_kg = if_else(bia_valid, ((height_cm^2) / bia_resistance_50khz * 0.401) + (age * -0.071) + 5.102 + if_else(sex == "male", 3.825, 0), NA_real_),
    smm_by_weight = smm_kg / weight_kg,

    # Changes from baseline
    perc_change_weight_kg = (weight_kg - weight_kg[visit == "v0"]) / weight_kg[visit == "v0"] * 100,
    perc_change_ffm_kg = if_else(bia_valid, (ffm_kg - ffm_kg[visit == "v0"]) / ffm_kg[visit == "v0"] * 100, NA_real_),
    perc_change_fm_kg = if_else(bia_valid, (fm_kg - fm_kg[visit == "v0"]) / fm_kg[visit == "v0"] * 100, NA_real_),
    delta_ffmi = ffmi - ffmi[visit == "v0"],
    perc_change_ffmi = (ffmi - ffmi[visit == "v0"]) / ffmi[visit == "v0"] * 100
  ) |>
  ungroup() |>
  group_by(sex, visit) |> 
  mutate( # sex-specific tertiles
    ffmi_tertile = quantile(ffmi, probs = 1/3, na.rm = TRUE),
    smm_kg_tertile = quantile(smm_kg, probs = 1/3, na.rm = TRUE),
    smm_by_weight_tertile = quantile(smm_by_weight, probs = 1/3, na.rm = TRUE),
    low_ffmi = if_else(ffmi <= ffmi_tertile, "yes", "no"),
    low_smm = if_else(smm_kg <= smm_kg_tertile, "yes", "no"),
    low_smm_by_weight = if_else(smm_by_weight <= smm_by_weight_tertile, "yes", "no")
  ) |> 
  ungroup()

## Wide dataset ##
baria_muscle_wide <- baria_muscle_long_all |> 
  select(-age_v0) |> # to avoid collision
  pivot_wider(
    names_from = visit,
    values_from =  c(
      all_of(long_vars), all_of(names(medication_classes)), "n_years_from_v0", "age", "hba1c_percent", "homa_ir", "homa_b", "t2d_labs", "prediab_labs",
      contains("ffmi"), contains("fmi"), contains("smm"), contains("bia_"), contains("perc_change_"), contains("delta_")
    ),
    names_glue = "{.value}_{visit}"
  ) |> 
  mutate(
    # New prediabetes occurence at follow-up in participants with NGT at baseline
    across(
      c(prediab_labs_v4, prediab_labs_v5),
      ~ case_when(
        t2d_v0 == "no" & t2d_labs_v0 == "no" & prediab_labs_v0 == "no" & .x == "yes" ~ "yes",
        t2d_v0 == "no" & t2d_labs_v0 == "no" & prediab_labs_v0 == "no" & .x == "no" ~ "no",
        TRUE ~ NA_character_
      ),
      .names = "new_{.col}"
    ),
    # New T2D occurence at follow-up in participants with NGT at baseline
    across(
      c(t2d_labs_v4, t2d_labs_v5),
      ~ case_when(
        t2d_v0 == "no" & t2d_labs_v0 == "no" & prediab_labs_v0 == "no" & .x == "yes" ~ "yes",
        t2d_v0 == "no" & t2d_labs_v0 == "no" & prediab_labs_v0 == "no" & .x == "no" ~ "no",
        TRUE ~ NA_character_
      ),
      .names = "new_{.col}"
    )
  )|> 
  group_by(sex) |> 
  mutate(
    across(starts_with("perc_change_ffmi_v"), ~ quantile(.x, probs = 1/3, na.rm = TRUE), .names = "{.col}_tertile"),
    across(
      starts_with("perc_change_ffmi_v") & !ends_with("_tertile"), 
      ~ if_else(
        .x <= quantile(.x, probs = 1/3, na.rm = TRUE),
        "high", "modest/low"
      ),
      .names = "{.col}_group"
    )
  )  |> 
  ungroup() |> 
  mutate(across(where(is.character) & !matches("^id$"), as.factor))

# Long dataset cleaned up
baria_muscle_long <- baria_muscle_long_all |> 
  select(
    id,visit,

    # Static/baseline covariates intentionally retained
    age_v0, sex, t2d_v0, height_cm, sg_type,

    # Longitudinal variables
    all_of(long_vars), all_of(names(medication_classes)),
    date_baseline, n_years_from_v0, age,
    hba1c_percent, homa_ir, homa_b, t2d_labs, prediab_labs,
    contains("ffmi"), contains("fmi"), contains("smm"), contains("bia_"), contains("perc_change_"), contains("delta_"), starts_with("low_")
  )

#### Microbiome Data Cleaning ####
# Clean microbiome metadata
sample_data(run1_mb)$visit <- case_when(
  sample_data(run1_mb)$Time_Point == "V-1" ~ "v0",
  sample_data(run1_mb)$Time_Point == "V4" ~ "v4",
  sample_data(run1_mb)$Time_Point == "V5" ~ "v5",
  TRUE ~ NA_character_
)

# Clean sample IDs in compositional table
sample_data(run1_mb)$id <- as.character(sample_data(run1_mb)$Subject_ID)
sample_names(run1_mb) <- str_c("BARIA_", str_remove(str_remove(str_remove(sample_names(run1_mb), "BARIA.Metagenome."), ".Fecal"), ".NA"))
sample_names(run1_mb) <- str_replace(str_replace(str_replace(sample_names(run1_mb), "V-1", "v0"), "V4", "v4"), "V5", "v5")
sample_names(run1_mb) <- str_replace(sample_names(run1_mb), "\\.", "_")
sample_names(run1_mb) <- str_replace(sample_names(run1_mb), ".rep\\d", "") # to remove .rep1

all(sample_names(run1_mb) == str_c("BARIA_", sample_data(run1_mb)$id, "_", sample_data(run1_mb)$visit))
# so we couldve also done sample_names(run1_mb) <- str_c("BARIA_", sample_data(run1_mb)$id, "_", sample_data(run1_mb)$visit)
# but that feels less safe

# Restrict to final Baria muscle cohort
baria_mb_clean <- prune_samples(sample_data(run1_mb)$id %in% bia_abx_mb_ids, run1_mb)

# Add relevant metadata to mb
baria_mb_metadata <- as(sample_data(baria_mb_clean), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(-any_of(c("sex", "ffmi_v0", "low_ffmi_v0", "perc_change_ffmi_v4_group", "perc_change_ffmi_v5_group"))) |>
  left_join(
    baria_muscle_wide |>
      select(id, sex, ffmi_v0, low_ffmi_v0, perc_change_ffmi_v4_group, perc_change_ffmi_v5_group),
    by = "id"
  ) |>
  column_to_rownames("Sample")

sample_data(baria_mb_clean) <- sample_data(baria_mb_metadata)

# Add a clean species label to the tax table
tax_table(baria_mb_clean) <- cbind(
  tax_table(baria_mb_clean),
  Tax = tax_table(baria_mb_clean)[, "Species"] |>
    str_remove("^s__") |>
    str_replace_all("_", " ")
)

# Shorten taxa names from the full taxonomy string to "species_SGB####"
species_label <- taxa_names(baria_mb_clean) |>
  str_extract("s__[^|]+") |>
  str_remove("^s__")
sgb <- taxa_names(baria_mb_clean) |>
  str_extract("t__.*$") |>
  str_remove("^t__")
taxa_names(baria_mb_clean) <- if_else(
  str_ends(species_label, sgb), species_label, str_c(species_label, "_", sgb)
)

# Unfiltered (incl. Eukaryota & GGB-labelled genera) - for diversity metrics,
# which should reflect the full profiled community rather than just the
# taxonomically well-annotated subset used for species-level models
baria_mb_unfiltered <- baria_mb_clean

# Drop eukaryotes
baria_mb_clean <- subset_taxa(baria_mb_clean, Kingdom != "k__Eukaryota")
# Drop GGB genera (without proper taxonomy)
baria_mb_clean <- subset_taxa(baria_mb_clean, !str_detect(Genus, "^g__GGB"))

baria_mb_baseline <- prune_samples(baria_mb_clean@sam_data$visit == "v0", baria_mb_clean)
baria_mb_baseline

#### Pathway data cleaning ####
humann <- readRDS("data/raw_data/BARIA.humann4.profiles.2026.581.910.RDS")

humann_long <- humann |>
  rownames_to_column(var = "pathway") |>
  pivot_longer(
    cols = -pathway,
    names_to = "Sample",
    values_to = "pathway_abundance"
  ) |>
  mutate(
    # Raw column names look like "BARIA.Metagenome.<id>.Fecal.V.1.NA_Abundance";
    # rebuild them into the same "BARIA_<id>_<visit>" ids used for the mb
    # sample_names (see the run1_mb renaming above)
    fecal_sample = str_detect(Sample, "\\.Fecal\\."),
    Sample = Sample |>
      str_remove("_Abundance$") |>
      str_replace("V\\.1", "V-1") |>
      str_remove("BARIA.Metagenome.") |>
      str_remove(".Fecal") |>
      str_remove(".NA"),
    Sample = str_c("BARIA_", Sample),
    Sample = Sample |>
      str_replace("V-1", "v0") |>
      str_replace("V4", "v4") |>
      str_replace("V5", "v5") |>
      str_replace("\\.", "_") |>
      str_replace(".rep\\d", "")
  ) |>
  filter(
    fecal_sample,
    Sample %in% sample_names(baria_mb_clean),
    !pathway %in% c("UNMAPPED", "UNINTEGRATED")
  ) |>
  select(-fecal_sample) |>
  mutate(
    pathway_id = str_extract(pathway, "[A-Z0-9-]+(?=:)"),
    pathway_name = str_extract(pathway, "(?<=:).*")
      |> trimws()
  ) |>
  left_join(
    as(sample_data(baria_mb_clean), "data.frame") |>
      rownames_to_column(var = "Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  ) |>
  mutate(pathway_abundance = replace_na(pathway_abundance, 0)) # Treat undetected pathways as zero abundance

# Save clinical data as both RDS and csv files
# Long clinical data
write.csv(baria_muscle_long, "data/processed_data/BARIA_muscle_long.csv")
saveRDS(baria_muscle_long, "data/processed_data/BARIA_muscle_long.RDS")

# Wide clinical data
write.csv(baria_muscle_wide, "data/processed_data/BARIA_muscle_wide.csv")
saveRDS(baria_muscle_wide, "data/processed_data/BARIA_muscle_wide.RDS")

# Save mb data
saveRDS(baria_mb_clean, "data/processed_data/BARIA_mb_clean.RDS")
saveRDS(baria_mb_unfiltered, "data/processed_data/BARIA_mb_clean_unfiltered.RDS")
saveRDS(baria_mb_baseline, "data/processed_data/BARIA_mb_baseline.RDS")

# Save pathway data
saveRDS(humann_long, "data/processed_data/BARIA_humann_pathways_long.RDS")