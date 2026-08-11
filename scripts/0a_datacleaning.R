# Baria project: Muscle mass trajectories and gut microbiota following bariatric surgery - Data cleaning
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)
library(phyloseq)

# Open Data and see properties
baria_clinical_data_raw <- readRDS("./data/raw_data/BARIA.clinical.2024-12-09.723.2043.RDS")
baria_mb <- readRDS("data/raw_data/ps.BARIA.metaphlan.706.2548.RDS")

## Formulas used in this script
# Hba1c(%) = (0,0915 * HbA1c (mmol/mol) + 2,15 (from diabetesfonds.nl)
# Hba1c(mmol/mol) = (10,93 x Hba1c (%)) - 23,5

## T2D Definition
# i. A1C ≥ 6.5% (≥ 48 mmol/mol) OR 
# ii. FPG ≥ 126 mg/dL (≥ 7.0 mmol/L) 
# (iii. 2h OGTT or random plasma glucose (not measured))

## Prediabetes Definition (ADA SOC 2026)
# i. A1C 5.7–6.4% (39–47 mmol/mol) OR
# ii. FPG 100 mg/dL (5.6 mmol/L) to 125 mg/dL (6.9 mmol/L) (IFG) OR
# (iii. 2h gluc (OGTT not available))

# DOI: 10.2337/diacare.21.12.2191 for HOMA calculations
# HOMA-IR = (Fasting insulin (µU/mL) × Fasting glucose (mmol/L)) / 22.5 (for glucose in mmol/L)
# HOMA-B = (20 × Fasting insulin (µU/mL)) / (Fasting glucose (mmol/L) − 3.5)

# Skeletal muscle mass (SMM) by Janssen: SMM [kg] = (height^2 [cm] / BIA-resistance [Ohms] X 0.401) + (gender x 3.825) + (age [years] x - 0.071)] + 5.102 (men = 1; women = 0)

#### Clinical Data ####
# Longitudinal vars
long_vars <- c("date",
  # Body composition
  "bmi", "weight_kg", "wc_cm", "fm_kg", "fm_percent", "ffm_kg", "ffm_percent", "bia_resistance_50khz", "upperleg_cm",

  # Lab/glycemia parameters
  "fasting_glucose_mmoll_mmt", "fasting_insulin_pmoll_mmt", "fasting_glucagon_ngl_mmt", "fasting_cpeptide_nmoll_mmt",
  "hba1c", "hba1c_mmolmol", "crp_mgl",

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
    glucose_mmoll_v0 = glucose, fasting_glucose_mmoll_mmt_v0 = min0gluc, fasting_insulin_pmoll_mmt_v0 = min0insulin, fasting_glucagon_ngl_mmt_v0 = min0glucagon, 
    fasting_cpeptide_nmoll_mmt_v0 = min0cpept, hba1c_v0 = hba1c, hba1c_mmolmol_v0 = hba1c__IFCC_mmolmol,
    gammagt_ul_v0 = ggt, alat_ul_v0 = alat, asat_ul_v0 = asat, triglycerides_mmoll_v0 = triglycerides,
    crp_mgl_v0 = crp, tsh_miul_v0 = tsh, ft4_pmoll_v0 = ft4,

    # v2-v5 (Repeated vars)
    matches("^V[2-5]_(date|bmi|weight|taille|tbf|tbf_percent|ffm|ffm_percent|rawdata_50Khz_Resistance|upperleg|min0gluc|min0insulin|min0glucagon|min0cpept|hba1c|hba1c__IFCC_mmolmol|crp)$"),

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
  rename_with(
    ~ .x |>
      str_replace("^weight_", "weight_kg_") |>
      str_replace("^taille_", "wc_cm_") |>
      str_replace("^tbf_percent_", "fm_percent_") |>
      str_replace("^tbf_", "fm_kg_") |>
      str_replace("^ffm_v", "ffm_kg_v") |>
      str_replace("^rawdata_50Khz_Resistance_", "bia_resistance_50khz_") |>
      str_replace("^upperleg_", "upperleg_cm_") |>
      str_replace("^min0gluc_", "fasting_glucose_mmoll_mmt_") |>
      str_replace("^min0insulin_", "fasting_insulin_pmoll_mmt_") |>
      str_replace("^min0glucagon_", "fasting_glucagon_ngl_mmt_") |>
      str_replace("^min0cpept_", "fasting_cpeptide_nmoll_mmt_") |>
      str_replace("^hba1c__IFCC_mmolmol_", "hba1c_mmolmol_") |>
      str_replace("^crp_", "crp_mgl_"),
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
  )
)

## Medication textbox cleaning
baria_muscle_vars_meds_notypos <- baria_muscle_vars |>
  mutate(medication_list_v0 = str_to_lower(medication_list_v0)) |> 
  mutate(
    medication_list_v0 = if_else(
      medication_list_v0 %in% c("geen", "geen medicatie"),
      "none",
      medication_list_v0 # replaced "geen" and "geen medicatie" with english "none"
      ), # fix typos I noticed in the medication list
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bamitriptiline\\b", "amitriptyline"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bduloxatine\\b", "duloxetine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bariprazol\\b", "aripiprazol"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bnifidipine\\b", "nifedipine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bevothyroxine\\b", "levothyroxine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bglicliazide\\b", "gliclazide"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bglicazide\\b", "gliclazide"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bamlopdipine\\b", "amlodipine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bzitalopram\\b", "citalopram"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bperinopril\\b", "perindopril"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bvenlaflaxine\\b", "venlafaxine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bglimerpiride\\b", "glimepiride"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bhydrchloorthiazide\\b", "hct"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bhydrochloorthiazide\\b", "hct"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bpnatoprazol\\b", "pantoprazole"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bomeprazol\\b", "omeprazole"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\besomeprazol\\b", "esomeprazole"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bparvastatine\\b", "pravastatine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bcandarsetan\\b", "candesartan"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bezatimibe\\b", "ezetimibe")
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

# Use regex patterns to create clean medication columns for all medication categories and subcategories
baria_muscle_vars_meds <- baria_muscle_vars_meds_notypos |> 
  mutate(
    # Diabetes medication
    metformin_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$metformin), "yes", "no"),
    sus_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$sus), "yes", "no"),
    dpp4is_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$dpp4is), "yes", "no"),
    glp1ras_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$glp1ras), "yes", "no"),
    sglt2is_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$sglt2is), "yes", "no"),
    tzds_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$tzds), "yes", "no"),
    insulin_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$insulin), "yes", "no"),
    dm_meds_v0 = if_else(str_detect(medication_list_v0, dm_meds_any_pattern), "yes", "no"), # Use of any diabetes medication

    # Antihypertensive medication
    ace_inhibitors_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$ace_inhibitors), "yes", "no"),
    arbs_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$arbs), "yes", "no"),
    ccbs_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$ccbs), "yes", "no"),
    bblockers_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$bblockers), "yes", "no"),
    a2_agonists_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$a2_agonists), "yes", "no"),
    diuretics_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$diuretics), "yes", "no"),
    combi_aht_meds_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$combi_aht_meds), "yes", "no"),
    aht_meds_v0 = if_else(str_detect(medication_list_v0, aht_meds_any_pattern), "yes", "no"), # general antihypertensive medication

    # Lipid-lowering medication
    statins_v0 = if_else(str_detect(medication_list_v0, lipidlowering_meds_patterns$statins), "yes", "no"),
    lipidlowering_meds_v0 = if_else(str_detect(medication_list_v0, lipidlowering_meds_any_pattern), "yes", "no"), # Use of any lipid-lowering medication

    # thyroid medication
    thyroid_meds_v0 = if_else(str_detect(medication_list_v0, thyroid_meds_pattern), "yes", "no"),

    # psychiatric medication
    ssris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$ssris), "yes", "no"),
    tcas_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$tcas), "yes", "no"),
    snris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$snris), "yes", "no"),
    ndris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$ndris), "yes", "no"),
    antipsychotics_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$antipsychotics), "yes", "no"),
    moods_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$moods), "yes", "no"),
    adhd_meds_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$adhd_meds), "yes", "no"),
    hypnotics_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$hypnotics), "yes", "no"),
    psychiatric_meds_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_any_pattern), "yes", "no"),
    ppi_v0 = if_else(str_detect(medication_list_v0, ppi_pattern), "yes", "no"),
    abx_v0 = if_else(str_detect(medication_list_v0, antibiotics_pattern), "yes", "no")
  ) |> 
  select(-medication_list_v0)

### Define final analysis cohort: valid baseline BIA & available baseline shotgun data & not taking antibiotics
# 1. Valid baseline BIA & 2. no antibiotics
bia_abx_v0_ids <- baria_muscle_vars_meds |>
  filter(
    # Baseline BIA QC
    visit == "v0",
    !is.na(ffm_kg),
    !is.na(fm_kg),
    bia_perc_diff <= 5,
    bia_kg_diff <= 5,
    bia_resistance_valid == TRUE,

    # No antibiotics
    abx_v0 == "no"
  ) |>
  pull(id) |> 
  unique()

# 3. Available & valid baseline microbiome/shotgun data
# Keep single runs and first run of duplicated samples (after comparing 1st and 2nd runs)
run1_mb <- prune_samples(
  sample_data(baria_mb)$Extra_data == "NA" |sample_data(baria_mb)$Extra_data == "rep1",
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
baria_muscle_long <- baria_muscle_vars_meds |> 
  filter(id %in% bia_abx_mb_ids) |> 
  mutate(date = dmy(date)) |> 
  group_by(id) |> 
  mutate(
    date_baseline = date[visit == "v0"],
    n_years_from_v0 = as.numeric(date - date_baseline) / 365.25,
    age = age_v0 + n_years_from_v0
  ) |> 
  ungroup() |> 
  mutate(
    hba1c_percent = if_else(hba1c < 15, hba1c, hba1c * 0.0915 + 2.15),
    hba1c_mmolmol = if_else(is.na(hba1c_mmolmol) == FALSE, hba1c_mmolmol, 10.93 * hba1c_percent - 23.5),

    # HOMA-IR & HOMA-2B (insulin unit conversion from pmol/l to uU/ml)
    homa_ir = (fasting_insulin_pmoll_mmt / 6.945) * fasting_glucose_mmoll_mmt / 22.5,
    homa_b = (20 * (fasting_insulin_pmoll_mmt / 6.945)) / (fasting_glucose_mmoll_mmt - 3.5),

    # T2D incidence based on lab values at follow-up
    t2d_labs = case_when(
      is.na(hba1c_percent) & is.na(fasting_glucose_mmoll_mmt) ~ NA_character_,
      hba1c_percent >= 6.5 | fasting_glucose_mmoll_mmt >= 7.0 ~ "yes",
      TRUE ~ "no"
    ),

    # Prediabetes based on lab values
    prediab_labs = case_when(
      t2d_labs == "yes" ~ "no",
      is.na(hba1c_percent) & is.na(fasting_glucose_mmoll_mmt) ~ NA_character_,
      (hba1c_percent >= 5.7 & hba1c_percent <= 6.4) | (fasting_glucose_mmoll_mmt >= 5.6 & fasting_glucose_mmoll_mmt <= 6.9) ~ "yes",
      TRUE ~ "no"
    ),
    ffmi = ffm_kg / ((height_cm / 100)^2),
    fmi = fm_kg / ((height_cm / 100)^2),
    smm_kg = if_else(
      bia_resistance_valid == TRUE, 
      ((height_cm^2) / bia_resistance_50khz * 0.401) + (age * -0.071) + 5.102 + if_else(sex == "male", 3.825, 0), NA_real_),
    smm_by_weight = smm_kg / weight_kg
  ) |>
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
baria_muscle_wide <- baria_muscle_long |> 
  select(-age_v0) |> # to avoid collision
  pivot_wider(
    names_from = visit,
    values_from =  c(
      all_of(long_vars), "n_years_from_v0", "age", "hba1c_percent", "homa_ir", "homa_b", "t2d_labs", "prediab_labs",
      contains("ffmi"), contains("fmi"), contains("smm"), contains("bia_")
    ),
    names_glue = "{.value}_{visit}"
  ) |> 
  mutate(
    across(c(weight_kg_v4, weight_kg_v5, weight_kg_v6), ~ (.x - weight_kg_v0) / weight_kg_v0 * 100, .names = "perc_change_{.col}"),
    across(c(ffm_kg_v4, ffm_kg_v5), ~ (.x - ffm_kg_v0) / ffm_kg_v0 * 100, .names = "perc_change_{.col}"),
    across(c(fm_kg_v4, fm_kg_v5), ~ (.x - fm_kg_v0) / fm_kg_v0 * 100, .names = "perc_change_{.col}"),
    across(c(ffmi_v4, ffmi_v5), ~ (.x - ffmi_v0), .names = "delta_{.col}"),
    across(c(ffmi_v4, ffmi_v5), ~ (.x - ffmi_v0) / ffmi_v0 * 100, .names = "perc_change_{.col}"),

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
        "high", "low/modest"
      ),
      .names = "{.col}_group"
    )
  )  |> 
  ungroup() |> 
  mutate(across(where(is.character) & !matches("^id$"), as.factor))

#### Microbiome Data Cleaning ####
# Clean microbiome metadata
sample_data(run1_mb)$visit <- case_when(
  sample_data(run1_mb)$Time_Point == "V-1" ~ "v0",
  sample_data(run1_mb)$Time_Point == "V4" ~ "v4",
  sample_data(run1_mb)$Time_Point == "V5" ~ "v5",
  TRUE ~ NA_character_
)

sample_data(run1_mb)$id <- as.character(sample_data(run1_mb)$Subject_ID)

# Restrict to final Baria muscle cohort
baria_mb_clean <- prune_samples(sample_data(run1_mb)$id %in% bia_abx_mb_ids, run1_mb)

# Add relevant metadata to mb
baria_mb_metadata <- as(sample_data(baria_mb_clean), "data.frame") |> 
  rownames_to_column(var = "Sample") |>
  left_join(
    baria_muscle_wide |> 
      select(id, sex, ffmi_v0, low_ffmi_v0),
    by = "id"
  ) |> 
  column_to_rownames("Sample")

sample_data(baria_mb_clean) <- sample_data(baria_mb_metadata)

# Save clinical data as both RDS and csv files
# Long clinical data
write.csv(baria_muscle_long, "data/processed_data/260810_BARIA_muscle_long.csv")
saveRDS(baria_muscle_long, "data/processed_data/260810_BARIA_muscle_long.RDS")

# Wide clinical data
write.csv(baria_muscle_wide, "data/processed_data/260810_BARIA_muscle_wide.csv")
saveRDS(baria_muscle_wide, "data/processed_data/260810_BARIA_muscle_wide.RDS")

# Save mb data
saveRDS(baria_mb_clean, "data/processed_data/260811_BARIA_mb_clean.RDS")
