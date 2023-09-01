#### Unemployment-UI-EQ Project Analysis Script 2 ####
#Author - Kieran Blaikie
#Date - 22 Aug 2023

#Overview - This script aims to perform the TWFE DiD approach introduced by
#           Imai, Kim & Wang. To do so, we:
#             1) Re-formats variables as-needed, truncates wealth, and restrict to the necessary set of variables and analysis observations
#                 - We do so separately for the 'complete case' dataset and each MI dataset
#             2) For the complete case and first MI dataset, we:
#                 a) Create matched set objects from the 'PanelMatch' package using:
#                     - PSM and Mahalanobis Distance for refinement
#                     - Restricting to the 5 and 10 most-similar non-exposed (0,0) observations per set
#                 b) Summarise matching sets and evaluate covariate balance/parallel trend assumptions across refinements
#             3) Using the optimum refinement method and refined set size, we then perform TWFE/DiD estimating the:
#                 - ATT Overall
#                 - CATT conditioning by:
#                     - Baseline EQ quartile, sex, or race and ethnicity
#                     - Length unemployed (if applicable)
#                     - UI recipiency (if applicable)
#                         + Baseline EQ quartile, sex, or race and ethnicity
#                         + Length unemployed (if applicable)
#             4) For MI analyses, we apply Rubin's Rules to construct point estimates and 95% CI

#Loading needed libraries
library(tidyverse) #Data-wrangling & manipulation
library(qs)

#Setting directory
directory <- "R:/Project/precarityR01/PSID/kieran_ui_project/Data/"

#Loading datasets
#Complete Case
cc_data <- qread(paste0(directory, "ui_formatted_complete_case_060823.qs"))

#Multiply Imputed
for (mi in 1:35) {
  assign(paste0("mi", mi, "_data"), qread(paste0(directory, "MI/ui_project_mi_", mi, "_060823.qs")))
  rm(mi)
}

#### Step 0 - Restrict to waves in the same state and occupational sector over the relevant period and create needed variables ####
#NOTE - Where immediate pre-exposure occupation is unknown, I sequentially search for immediate post-exposure occ, then more distant pre-exosure occs
for (dataset in ls(pattern = "_data")) {
  print(paste0("Dataset: ", dataset))
  temp <- eval(parse(text = dataset))
  
  library(data.table)
  temp <- as.data.table(temp)
  temp[, ':=' (occupation_broad = as.character(occupation_broad),
               occupation_broad_lag = as.character(occupation_broad_lag))]
  temp[, base_occupation := fifelse(!is.na(occupation_broad_lag) & occupation_broad_lag != "UnemployedorNILF", occupation_broad_lag,
                                    fifelse(!is.na(occupation_broad) & occupation_broad != "UnemployedorNILF", occupation_broad,
                                            fifelse(!is.na(lag(occupation_broad,2)) & lag(occupation_broad,2) != "UnemployedorNILF", lag(occupation_broad,2),
                                                    fifelse(!is.na(lag(occupation_broad,3)) & lag(occupation_broad,3) != "UnemployedorNILF", lag(occupation_broad,3),
                                                            fifelse(!is.na(lag(occupation_broad,4)) & lag(occupation_broad,4) != "UnemployedorNILF", lag(occupation_broad,4), NA_character_))))), by = .(unique_id)]
  temp[, same_state := fifelse(NROW(unique(state)) == 1, 1, NA_integer_), by = .(unique_id)]
  temp[NROW(unique(state)) >1, same_state := fifelse(is.na(lag(year,1)) == T | is.na(lag(year,2)) == T, 1,
                                                     fifelse(lag(year,1) != (year-2) | lag(year,2) != (year-4), NA_integer_,
                                                             fifelse(lag(state,1) != state | lag(state,2) != state, 0, 1))), by = .(unique_id)]
  temp[, same_occ := fifelse(NROW(unique(base_occupation)) == 1, 1, NA_integer_), by = .(unique_id)]
  temp[NROW(unique(base_occupation)) >1, same_occ := fifelse(is.na(lag(year,1)) == T | is.na(lag(year,2)) == T, 1,
                                                             fifelse(lag(year,1) != (year-2) | lag(year,2) != (year-4), NA_integer_,
                                                                     fifelse(lag(base_occupation,1) != base_occupation | lag(base_occupation,2) != base_occupation, 0, 1))), by = .(unique_id)]
  temp <- as.data.frame(temp)
  detach("package:data.table", unload = T)
  temp <- temp[(temp$same_state == 1 | is.na(temp$same_state)), ]
  temp <- temp[!is.na(temp$unique_id),]
  
  assign(dataset, temp)
  rm(temp)
}

#### Step 1 - Formatting, renaming, restricting variables and restricting observations ####
#NOTE - Here we: 
#         1) Convert covariates to factors
#         2) Re-parameterise 'year' to be sequential integers = (year - 1997)/2
#         3) Truncate family wealth at the 2nd and 98th percentile (family income has already been truncated and standardized)
#         4) Dichotomize or make ordinal covariates (PanelMatch only supports numerical covariates)
#            Note 'exact matching' can be done on categorical covariates (e.g. state, occupation)
#         5) Restricting respondent observations to those without missing exposure or outcome values and 
#            restricting respondent observation sets to the maximum with observations still consecutive (e.g. T=1/2/3, not 1/3/4)
for (dataset in ls(pattern = "_data")) {
  print(paste0("Dataset: ", dataset))
  temp <- eval(parse(text = dataset))
  
  #1) Convert covariates to factors
  factors <- c("state", "region", "gender", "race", "ethnicity", "nativity", 
               "parents_poor", "marital_status", "education", "SRH", 
               "base_occupation", "SRH_lag", "marital_status_lag", 
               "education_lag")
  temp[, factors] <- lapply(temp[, factors], factor)
  rm(factors)
  
  #2) Re-parameterize year as sequential integers
  temp$year <- as.integer(temp$year/2)
  
  #3) Truncate income and wealth variables
  temp$family_wealth_no_home_equity <- ifelse(temp$family_wealth_no_home_equity < quantile(temp$family_wealth_no_home_equity, c(0.02), na.rm=T),
                                              quantile(temp$family_wealth_no_home_equity, c(0.02), na.rm=T), temp$family_wealth_no_home_equity)
  temp$family_wealth_no_home_equity <- ifelse(temp$family_wealth_no_home_equity > quantile(temp$family_wealth_no_home_equity, c(0.98), na.rm=T),
                                              quantile(temp$family_wealth_no_home_equity, c(0.98), na.rm=T), temp$family_wealth_no_home_equity)
  temp$family_wealth_no_home_equity_lag <- ifelse(temp$family_wealth_no_home_equity_lag < quantile(temp$family_wealth_no_home_equity_lag, c(0.02), na.rm=T),
                                                  quantile(temp$family_wealth_no_home_equity_lag, c(0.02), na.rm=T), temp$family_wealth_no_home_equity_lag)
  temp$family_wealth_no_home_equity_lag <- ifelse(temp$family_wealth_no_home_equity_lag > quantile(temp$family_wealth_no_home_equity_lag, c(0.98), na.rm=T),
                                                  quantile(temp$family_wealth_no_home_equity_lag, c(0.98), na.rm=T), temp$family_wealth_no_home_equity_lag)
  
  #4) Dichotomize or make ordinal categorical covariates
  temp %>% mutate(male = ifelse(gender == "Male", 1, 0),
                  black = ifelse(is.na(race), NA_integer_, ifelse(race == "Black", 1, 0)),
                  other = ifelse(is.na(race), NA_integer_, ifelse(race == "Other", 1, 0)),
                  hispanic = ifelse(is.na(ethnicity), NA_integer_, ifelse(ethnicity == "Hispanic", 1, 0)),
                  non_native = ifelse(is.na(nativity), NA_integer_, ifelse(nativity == "NotUS", 1, 0)),
                  childhood_ses_poor = ifelse(is.na(parents_poor), NA_integer_, ifelse(parents_poor == "Poor", 1, 0)),
                  married = ifelse(is.na(marital_status), NA_integer_, ifelse(marital_status == "MarriedCohabiting", 1, 0)),
                  married_lag = ifelse(is.na(marital_status_lag), NA_integer_, ifelse(marital_status_lag == "MarriedCohabiting", 1, 0)),
                  lesshs = ifelse(is.na(education), NA_integer_, ifelse(education == "LessHS", 1, 0)),
                  college = ifelse(is.na(education), NA_integer_, ifelse(education == "College", 1, 0)),
                  lesshs_lag = ifelse(is.na(education_lag), NA_integer_, ifelse(education_lag == "LessHS", 1, 0)),
                  college_lag = ifelse(is.na(education_lag), NA_integer_, ifelse(education_lag == "College", 1, 0)),
                  srh_vgood_exc = ifelse(is.na(SRH), NA_integer_, ifelse(SRH %in% c("Verygood", "Excellent"), 1, 0)),
                  srh_vgood_exc_lag = ifelse(is.na(SRH_lag), NA_integer_, ifelse(SRH_lag %in% c("Verygood", "Excellent"), 1, 0)),
                  base_occ_farmforfish = ifelse(is.na(base_occupation), NA_integer_, ifelse(base_occupation == "Farmingforestryandfishing", 1, 0)),
                  base_occ_managerial = ifelse(is.na(base_occupation), NA_integer_, ifelse(base_occupation == "Managerial", 1, 0)),
                  base_occ_military = ifelse(is.na(base_occupation), NA_integer_, ifelse(base_occupation == "Military", 1, 0)),
                  base_occ_opfablab = ifelse(is.na(base_occupation), NA_integer_, ifelse(base_occupation == "Operatorsfabricatorsandlaborers", 1, 0)),
                  base_occ_precision = ifelse(is.na(base_occupation), NA_integer_, ifelse(base_occupation == "Precisionproductioncraftandrepair", 1, 0)),
                  base_occ_profspec = ifelse(is.na(base_occupation), NA_integer_, ifelse(base_occupation == "Professionalspecialty", 1, 0)),
                  base_occ_services = ifelse(is.na(base_occupation), NA_integer_, ifelse(base_occupation == "Services", 1, 0))) -> temp
  
  #Restricting to necessary variables
  temp <- temp %>% select(unique_id, year, unemp_status_lead, 
                          ui_receipt_6mo_post_unemp, max_week_duration, max_weekly_ben_dep, 
                          Y_lead, Y_lead_year, Y_base, Y_base_year, 
                          state, age, male, black, other, hispanic, non_native, childhood_ses_poor, 
                          married, married_lag, lesshs, lesshs_lag, college, college_lag, 
                          disabl_limits_work, srh_vgood_exc, srh_vgood_exc_lag, family_income,
                          family_income_lag, family_wealth_no_home_equity, family_wealth_no_home_equity_lag, 
                          base_occ_farmforfish, base_occ_managerial, base_occ_military, 
                          base_occ_opfablab, base_occ_precision, base_occ_profspec, base_occ_services,
                          prepost, mo_unemp_pre_int_min, mo_unemp_until_emp, eligible, include, gsp_per_cap, unemp_rate_quart)
  
  #Renaming variables
  names(temp)[c(3:7)] <- c("exposure", "moderator_receipt", "moderator_max_weeks", "moderator_max_amount", "outcome")
  
  #Standardising post-exposure EQ outcome
  temp$outcome <- (temp$outcome - mean(temp$outcome, na.rm = T)) / sd(temp$outcome, na.rm = T)
  
  #Restricting to observations where individuals are provisionally eligible and waves are consecutive
  temp <- temp[temp$include == 1 & !is.na(temp$include), ]
  temp$include <- NULL
  
  #Restricting to the maximum set per person where analysis observations (i,t) or (i',t) are not missing exposure or outcome status
  temp <- temp %>% filter(!is.na(exposure))
  temp <- temp %>% filter(!is.na(outcome)) #Missing lagged outcome status is considered fine.
  
  #Restricting to the maximum set per person where analysis observations per i are still consecutive and there are at least 3 observations
  #Creating needed vars in temp
  temp %>% mutate(max_period = NA_real_, include = NA_real_) -> temp
  
  #Identifying all included individuals
  individuals <- unique(temp$unique_id)
  
  #Creating a loop where, for each person, we:
  # - Create a running counter for the maximum number of consecutive eligible waves
  # - Identify the earliest and latest year relevant to this set of eligible observations
  # - Fill in the 'max_period' and 'include' vars for that person
  for (i in 1:NROW(individuals)) {
    ind <- individuals[i]
    print(paste0("Dataset: ", dataset, ". Person: ", i, ". ", round((i/NROW(individuals)*100), 2), "%"))
    max_consecutive_count <- 0
    current_consecutive_count <- 0
    eligible_start_year <- NA
    eligible_end_year <- NA
    for (j in 1:nrow(temp[temp$unique_id == ind, ])) {
      if (temp[temp$unique_id == ind, ]$eligible[j] == 1) {
        current_consecutive_count <- current_consecutive_count + 1
        if (current_consecutive_count > max_consecutive_count) {
          max_consecutive_count <- current_consecutive_count
          eligible_start_year <- (temp[temp$unique_id == ind, ]$year[j - (current_consecutive_count - 1)]) - 2
          eligible_end_year <- temp[temp$unique_id == ind, ]$year[j]
        }
      } else {
        current_consecutive_count <- 0
      }
    }
    temp$max_period[temp$unique_id == ind] <- max_consecutive_count
    temp[temp$unique_id == ind, ]$include <- ifelse(temp[temp$unique_id == ind, ]$year >= eligible_start_year & temp[temp$unique_id == ind, ]$year <= eligible_end_year, 1, 0)
  }
  rm(individuals, ind, max_consecutive_count, current_consecutive_count, eligible_start_year, eligible_end_year, i, j)
  
  temp <- temp[temp$include == 1 & !is.na(temp$include), ]
  
  #Saving updated 'dataset' after determining waves to keep for analyses
  if (dataset == "cc_data") {
    qsave(temp, paste0(directory, "ui_formatted_complete_case_analysis_220823.qs"))
  } else {
    mi <- dataset
    mi <- str_replace(mi, "mi", "")
    mi <- str_replace(mi, "_data", "")
    qsave(temp, paste0(directory, "MI/ui_project_mi_", mi, "_analysis_220823.qs"))
    rm(mi)
  }
  assign(dataset, temp)
  rm(temp)
}

#### Step 2 - Creating matched set objects ####
#Loading PanelMatch
library(PanelMatch)

## Part A - Creating matched set objects
for (data in c("cc_data", "mi1_data")) {
  print(paste0("Dataset: ", data))
  temp <- eval(parse(text = data))
  
  #Creating the initial matched set without weighting
  assign(paste0(data, "_matched_unrefined"), PanelMatch(lag = 2, time.id = "year", unit.id = "unique_id",
                                                        treatment = "exposure", refinement.method = "none",
                                                        data = temp, match.missing = FALSE,
                                                        size.match = 10, qoi = "att", outcome.var = "outcome", exact.match.variables = c("state"),
                                                        lead = 0, forbid.treatment.reversal = FALSE,
                                                        use.diagonal.variance.matrix = TRUE))
  
  #Creating matched sets using Mahalanobis and PSM refinement with 5 and 10 i' observation set sizes
  for (refinement in c("mahalanobis", "ps.match")) {
    if (refinement == "mahalanobis") { ref_name <- "mahalanobis"}
    if (refinement == "ps.match") { ref_name <- "psm"}
    for (size in c(5, 10)) {
      assign(paste0(data, "_matched_", ref_name, "_", size), PanelMatch(lag = 2, time.id = "year", unit.id = "unique_id",
                                                                        treatment = "exposure", refinement.method = refinement,
                                                                        data = temp, match.missing = FALSE,
                                                                        covs.formula = ~ Y_base + Y_base_year + age + male + black + other + 
                                                                          hispanic + non_native + childhood_ses_poor + married + lesshs + 
                                                                          college + srh_vgood_exc + family_income + family_wealth_no_home_equity + 
                                                                          married_lag + lesshs_lag + college_lag + srh_vgood_exc_lag + family_income_lag + 
                                                                          family_wealth_no_home_equity_lag + base_occ_farmforfish + base_occ_managerial +
                                                                          base_occ_military + base_occ_opfablab + base_occ_precision + base_occ_profspec + 
                                                                          base_occ_services + mo_unemp_pre_int_min + gsp_per_cap + unemp_rate_quart,
                                                                        size.match = size, qoi = "att", outcome.var = "outcome", exact.match.variables = c("state"),
                                                                        lead = 0, forbid.treatment.reversal = FALSE,
                                                                        use.diagonal.variance.matrix = TRUE))
      rm(size)
    }
    rm(refinement, ref_name)
  }
  rm(temp)
}
rm(data)

## Part B - Summarise matching sets and evaluate covariate balance/parallel trend assumptions across refinements
# Complete Case Comparisons
#Extracting matched.set objects
cc_data_mset_unrefined <- cc_data_matched_unrefined$att
cc_data_mset_mahalanobis_5 <- cc_data_matched_mahalanobis_5$att
cc_data_mset_mahalanobis_10 <- cc_data_matched_mahalanobis_10$att
cc_data_mset_psm_5 <- cc_data_matched_psm_5$att
cc_data_mset_psm_10 <- cc_data_matched_psm_10$att

#Creating Balance Scatterplot
balance_scatter(matched_set_list = list(cc_data_mset_mahalanobis_5, cc_data_mset_mahalanobis_10, cc_data_mset_psm_5, cc_data_mset_psm_10),
                covariates = c("Y_base", "Y_base_year", "age", "male", "black", "other", 
                               "hispanic", "non_native", "childhood_ses_poor", "married", "lesshs", 
                               "college", "srh_vgood_exc", "family_income", "family_wealth_no_home_equity", 
                               "married_lag", "lesshs_lag", "college_lag", "srh_vgood_exc_lag", "family_income_lag", 
                               "family_wealth_no_home_equity_lag", "base_occ_farmforfish", "base_occ_managerial",
                               "base_occ_military", "base_occ_opfablab", "base_occ_precision", "base_occ_profspec", 
                               "base_occ_services", "mo_unemp_pre_int_min", "gsp_per_cap", "unemp_rate_quart"),
                data = cc_data, pchs = c(1,2,15,16), xlim = c(0,0.4), ylim = c(0,0.4))
legend(x = 0, y = 0.4, legend = c("Mahalanobis-5", "Mahalanobis-10", "PSM-5", "PSM-10"),
       pch = c(1,2,15,16), pt.cex = 1, bty = "n", ncol = 1, cex = 1, bg = "white")

# First MI Dataset Comparisons
#Extracting matched.set objects
mi1_data_mset_unrefined <- mi1_data_matched_unrefined$att
mi1_data_mset_mahalanobis_5 <- mi1_data_matched_mahalanobis_5$att
mi1_data_mset_mahalanobis_10 <- mi1_data_matched_mahalanobis_10$att
mi1_data_mset_psm_5 <- mi1_data_matched_psm_5$att
mi1_data_mset_psm_10 <- mi1_data_matched_psm_10$att

#Creating Balance Scatterplot
balance_scatter(matched_set_list = list(mi1_data_mset_mahalanobis_5, mi1_data_mset_mahalanobis_10, mi1_data_mset_psm_5, mi1_data_mset_psm_10),
                covariates = c("Y_base", "Y_base_year", "age", "male", "black", "other", 
                               "hispanic", "non_native", "childhood_ses_poor", "married", "lesshs", 
                               "college", "srh_vgood_exc", "family_income", "family_wealth_no_home_equity", 
                               "married_lag", "lesshs_lag", "college_lag", "srh_vgood_exc_lag", "family_income_lag", 
                               "family_wealth_no_home_equity_lag", "base_occ_farmforfish", "base_occ_managerial",
                               "base_occ_military", "base_occ_opfablab", "base_occ_precision", "base_occ_profspec", 
                               "base_occ_services", "mo_unemp_pre_int_min", "gsp_per_cap", "unemp_rate_quart"),
                data = mi1_data, pchs = c(1,2,15,16), xlim = c(0,0.4), ylim = c(0,0.4))
legend(x = 0, y = 0.4, legend = c("Mahalanobis-5", "Mahalanobis-10", "PSM-5", "PSM-10"),
       pch = c(1,2,15,16), pt.cex = 1, bty = "n", ncol = 1, cex = 1, bg = "white")

#### Step 3 - Using the optimum refinement method and refined set size, we then perform TWFE/DiD ####
#NOTE - Based on Step 2, propensity score matching with up to 10 i' observations per set appears optimum
#For the complete case dataset and each MI dataset We estimate the following quantities:
#  - ATT Overall
#  - CATT conditioning by:
#    - Baseline EQ quartile, sex, or race and ethnicity
#    - Length unemployed (if applicable)
#    - UI recipiency (if applicable)
#      + UI maximum duration (< Median, >= Median)
#      + UI maximum generosity (< Median, >= Median)
#      + Baseline EQ quartile, sex, or race and ethnicity
#      + Length unemployed (if applicable)

#Using UI state-level data to determine state median maximum generosity and duration over the 2001-2017 period
ui_state_data <- read.csv(paste0(directory, "sarah_state_ui_dataset_050123.csv"))
ui_state_data <- ui_state_data[ui_state_data$state != "ok", ] #Excluding 2016 record with state=="ok" as a duplicate of the 2016 record with state=="OK"
ui_state_data$year <- as.numeric(substr(ui_state_data$date, start = 5, stop = 8)) #Creating year var 
ui_state_data %>% filter(year %in% c(2001:2017)) -> ui_state_data

med_max_generosity <- quantile(ui_state_data$max_weekly_ben_dep, c(0.5), na.rm=T)
med_max_duration <- quantile(ui_state_data$max_week_duration, c(0.5), na.rm=T)

datasets <- ls()[endsWith(ls(), "_data") == T]
for (dataset in datasets) {
  print(paste0("Dataset: ", dataset))
  temp <- eval(parse(text = dataset))
  
  #Creating stratifying variables
  temp %>% mutate(baseline_eq_quartile = case_when(is.na(Y_base) ~ NA_real_,
                                                   Y_base < quantile(Y_base, c(0.25), na.rm=T) ~ 1,
                                                   Y_base < quantile(Y_base, c(0.5), na.rm=T) ~ 2,
                                                   Y_base < quantile(Y_base, c(0.75), na.rm=T) ~ 3,
                                                   TRUE ~ 4),
                  race_eth = case_when(hispanic == 1 ~ "Hispanic",
                                       is.na(hispanic) ~ NA_character_,
                                       black == 1 ~ "NHBlack",
                                       other == 1 ~ "NHOther",
                                       black == 0 & other == 0 ~ "NHWhite",
                                       TRUE ~ NA_character_),
                  lessmed_length_unemp = case_when(exposure == 0 ~ NA_real_,
                                                   exposure == 1 & mo_unemp_until_emp < quantile(mo_unemp_until_emp[exposure == 1], c(0.5), na.rm=T) ~ 1,
                                                   exposure == 1 & mo_unemp_until_emp >= quantile(mo_unemp_until_emp[exposure == 1], c(0.5), na.rm=T) ~ 0,
                                                   TRUE ~ NA_real_),
                  ui_max_dur_less_med = case_when(moderator_max_weeks < med_max_duration ~ 1,
                                                  moderator_max_weeks >= med_max_duration ~ 0,
                                                  TRUE ~ NA_real_),
                  ui_max_gen_less_med = case_when(moderator_max_amount < med_max_generosity ~ 1,
                                                  moderator_max_amount >= med_max_generosity ~ 0,
                                                  TRUE ~ NA_real_),
                  ui_cat_dur = case_when(moderator_receipt == 0 ~ "None",
                                         moderator_receipt == 1 & ui_max_dur_less_med == 1 ~ "YesLessMedDur",
                                         moderator_receipt == 1 & ui_max_dur_less_med == 0 ~ "YesMoreMedDur",
                                         TRUE ~ NA_character_),
                  ui_cat_gen = case_when(moderator_receipt == 0 ~ "None",
                                         moderator_receipt == 1 & ui_max_gen_less_med == 1 ~ "YesLessMedGen",
                                         moderator_receipt == 1 & ui_max_gen_less_med == 0 ~ "YesMoreMedGen",
                                         TRUE ~ NA_character_)) -> temp
  #Creating empty results shell
  results <- data.frame(data = dataset, estimand = NA_character_,
                        #Matched Set Information for ATT
                        N_sets = NA_real_, N_sets_unmatched = NA_real_, N_sets_matched = NA_real_,
                        set_size_med = NA_real_, set_size_1qr = NA_real_, set_size_3qr = NA_real_,
                        #Estimate
                        est = NA_real_, se = NA_real_, lci = NA_real_, uci = NA_real_)
  
  #Creating the initial matched set
  match_object <- PanelMatch(lag = 2, time.id = "year", unit.id = "unique_id",
                             treatment = "exposure", refinement.method = "ps.match",
                             data = temp, match.missing = FALSE,
                             covs.formula = ~ Y_base + Y_base_year + age + male + black + other + 
                               hispanic + non_native + childhood_ses_poor + married + lesshs + 
                               college + srh_vgood_exc + family_income + family_wealth_no_home_equity + 
                               married_lag + lesshs_lag + college_lag + srh_vgood_exc_lag + family_income_lag + 
                               family_wealth_no_home_equity_lag + base_occ_farmforfish + base_occ_managerial +
                               base_occ_military + base_occ_opfablab + base_occ_precision + base_occ_profspec + 
                               base_occ_services + mo_unemp_pre_int_min + gsp_per_cap + unemp_rate_quart,
                             size.match = 10, qoi = "att", outcome.var = "outcome", exact.match.variables = c("state"),
                             lead = 0, forbid.treatment.reversal = FALSE,
                             use.diagonal.variance.matrix = TRUE)
  
  #Extracting the matched set object
  mset <- match_object$att
  set_details <- print(mset)
  
  #Populating ATT-specific results 'estimand', set' and 'set_size' variables
  att <- results
  att[1, c(2:8)] <- c("ATT",
                      (NROW(set_details)), 
                      (NROW(set_details[set_details$matched.set.size == 0,])),
                      (NROW(set_details[set_details$matched.set.size != 0,])),
                      (summary(set_details[set_details$matched.set.size >= 0,]$matched.set.size)[[3]]),
                      (summary(set_details[set_details$matched.set.size >= 0,]$matched.set.size)[[2]]),
                      (summary(set_details[set_details$matched.set.size >= 0,]$matched.set.size)[[5]]))
  
  #Estimating the ATT
  twfe_att <- PanelEstimate(sets = match_object, data = temp, se.method = "bootstrap", confidence.level = .95)
  
  #Populating ATT-specific results 'att' variables
  att[1,]$est <- summary(twfe_att)$summary[[1]]
  att[1,]$se <- summary(twfe_att)$summary[[2]]
  att[1,]$lci <- summary(twfe_att)$summary[[3]]
  att[1,]$uci <- summary(twfe_att)$summary[[4]]
  
  #Estimating CATT by Baseline EQ
  twfe_catt_eq <- PanelEstimate(sets = match_object, data = temp, moderator = "baseline_eq_quartile", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_eq_q1_sets <- print(twfe_catt_eq[[1]]$matched.sets)
  catt_eq_q2_sets <- print(twfe_catt_eq[[2]]$matched.sets)
  catt_eq_q3_sets <- print(twfe_catt_eq[[3]]$matched.sets)
  catt_eq_q4_sets <- print(twfe_catt_eq[[4]]$matched.sets)
  
  catt_eq_q1 <- results; catt_eq_q1$estimand <- "CATT_EQ_Q1"
  catt_eq_q2 <- results; catt_eq_q2$estimand <- "CATT_EQ_Q2"
  catt_eq_q3 <- results; catt_eq_q3$estimand <- "CATT_EQ_Q3"
  catt_eq_q4 <- results; catt_eq_q4$estimand <- "CATT_EQ_Q4"
  
  catt_eq_q1[1,3:8] <- c((NROW(catt_eq_q1_sets)),
                         (NROW(catt_eq_q1_sets$matched.set.size[catt_eq_q1_sets$matched.set.size == 0])),
                         (NROW(catt_eq_q1_sets$matched.set.size[catt_eq_q1_sets$matched.set.size > 0])),
                         (quantile(catt_eq_q1_sets$matched.set.size[catt_eq_q1_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_eq_q2[1,3:8] <- c((NROW(catt_eq_q2_sets)),
                         (NROW(catt_eq_q2_sets$matched.set.size[catt_eq_q2_sets$matched.set.size == 0])),
                         (NROW(catt_eq_q2_sets$matched.set.size[catt_eq_q2_sets$matched.set.size > 0])),
                         (quantile(catt_eq_q2_sets$matched.set.size[catt_eq_q2_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_eq_q3[1,3:8] <- c((NROW(catt_eq_q3_sets)),
                         (NROW(catt_eq_q3_sets$matched.set.size[catt_eq_q3_sets$matched.set.size == 0])),
                         (NROW(catt_eq_q3_sets$matched.set.size[catt_eq_q3_sets$matched.set.size > 0])),
                         (quantile(catt_eq_q3_sets$matched.set.size[catt_eq_q3_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_eq_q4[1,3:8] <- c((NROW(catt_eq_q4_sets)),
                         (NROW(catt_eq_q4_sets$matched.set.size[catt_eq_q4_sets$matched.set.size == 0])),
                         (NROW(catt_eq_q4_sets$matched.set.size[catt_eq_q4_sets$matched.set.size > 0])),
                         (quantile(catt_eq_q4_sets$matched.set.size[catt_eq_q4_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_eq_q1[1,9:12] <- summary(twfe_catt_eq[[1]])$summary[1:4]
  catt_eq_q2[1,9:12] <- summary(twfe_catt_eq[[2]])$summary[1:4]
  catt_eq_q3[1,9:12] <- summary(twfe_catt_eq[[3]])$summary[1:4]
  catt_eq_q4[1,9:12] <- summary(twfe_catt_eq[[4]])$summary[1:4]
  
  #Estimating CATT by Sex
  twfe_catt_sex <- PanelEstimate(sets = match_object, data = temp, moderator = "male", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_sex_f_sets <- print(twfe_catt_sex[["0"]]$matched.sets)
  catt_sex_m_sets <- print(twfe_catt_sex[["1"]]$matched.sets)
  
  catt_sex_f <- results; catt_sex_f$estimand <- "CATT_SEX_F"
  catt_sex_m <- results; catt_sex_m$estimand <- "CATT_SEX_M"
  
  catt_sex_f[1,3:8] <- c((NROW(catt_sex_f_sets)),
                         (NROW(catt_sex_f_sets$matched.set.size[catt_sex_f_sets$matched.set.size == 0])),
                         (NROW(catt_sex_f_sets$matched.set.size[catt_sex_f_sets$matched.set.size > 0])),
                         (quantile(catt_sex_f_sets$matched.set.size[catt_sex_f_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_sex_m[1,3:8] <- c((NROW(catt_sex_m_sets)),
                         (NROW(catt_sex_m_sets$matched.set.size[catt_sex_m_sets$matched.set.size == 0])),
                         (NROW(catt_sex_m_sets$matched.set.size[catt_sex_m_sets$matched.set.size > 0])),
                         (quantile(catt_sex_m_sets$matched.set.size[catt_sex_m_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_sex_f[1,9:12] <- summary(twfe_catt_sex[["0"]])$summary[1:4]
  catt_sex_m[1,9:12] <- summary(twfe_catt_sex[["1"]])$summary[1:4]
  
  #Estimating CATT by Race and Ethnicity
  twfe_catt_raceeth <- PanelEstimate(sets = match_object, data = temp, moderator = "race_eth", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_raceeth_nhwhite_sets <- print(twfe_catt_raceeth[["NHWhite"]]$matched.sets)
  catt_raceeth_nhblack_sets <- print(twfe_catt_raceeth[["NHBlack"]]$matched.sets)
  catt_raceeth_nhother_sets <- print(twfe_catt_raceeth[["NHOther"]]$matched.sets)
  catt_raceeth_hispanic_sets <- print(twfe_catt_raceeth[["Hispanic"]]$matched.sets)
  
  catt_raceeth_nhwhite <- results; catt_raceeth_nhwhite$estimand <- "CATT_RACEETH_NHWHITE"
  catt_raceeth_nhblack <- results; catt_raceeth_nhblack$estimand <- "CATT_RACEETH_NHBLACK"
  catt_raceeth_nhother <- results; catt_raceeth_nhother$estimand <- "CATT_RACEETH_NHOTHER"
  catt_raceeth_hispanic <- results; catt_raceeth_hispanic$estimand <- "CATT_RACEETH_HISPANIC"
  
  catt_raceeth_nhwhite[1,3:8] <- c((NROW(catt_raceeth_nhwhite_sets)),
                                   (NROW(catt_raceeth_nhwhite_sets$matched.set.size[catt_raceeth_nhwhite_sets$matched.set.size == 0])),
                                   (NROW(catt_raceeth_nhwhite_sets$matched.set.size[catt_raceeth_nhwhite_sets$matched.set.size > 0])),
                                   (quantile(catt_raceeth_nhwhite_sets$matched.set.size[catt_raceeth_nhwhite_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_raceeth_nhblack[1,3:8] <- c((NROW(catt_raceeth_nhblack_sets)),
                                   (NROW(catt_raceeth_nhblack_sets$matched.set.size[catt_raceeth_nhblack_sets$matched.set.size == 0])),
                                   (NROW(catt_raceeth_nhblack_sets$matched.set.size[catt_raceeth_nhblack_sets$matched.set.size > 0])),
                                   (quantile(catt_raceeth_nhblack_sets$matched.set.size[catt_raceeth_nhblack_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_raceeth_nhother[1,3:8] <- c((NROW(catt_raceeth_nhother_sets)),
                                   (NROW(catt_raceeth_nhother_sets$matched.set.size[catt_raceeth_nhother_sets$matched.set.size == 0])),
                                   (NROW(catt_raceeth_nhother_sets$matched.set.size[catt_raceeth_nhother_sets$matched.set.size > 0])),
                                   (quantile(catt_raceeth_nhother_sets$matched.set.size[catt_raceeth_nhother_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_raceeth_hispanic[1,3:8] <- c((NROW(catt_raceeth_hispanic_sets)),
                                    (NROW(catt_raceeth_hispanic_sets$matched.set.size[catt_raceeth_hispanic_sets$matched.set.size == 0])),
                                    (NROW(catt_raceeth_hispanic_sets$matched.set.size[catt_raceeth_hispanic_sets$matched.set.size > 0])),
                                    (quantile(catt_raceeth_hispanic_sets$matched.set.size[catt_raceeth_hispanic_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_raceeth_nhwhite[1,9:12] <- summary(twfe_catt_raceeth[["NHWhite"]])$summary[1:4]
  catt_raceeth_nhblack[1,9:12] <- summary(twfe_catt_raceeth[["NHBlack"]])$summary[1:4]
  catt_raceeth_nhother[1,9:12] <- summary(twfe_catt_raceeth[["NHOther"]])$summary[1:4]
  catt_raceeth_hispanic[1,9:12] <- summary(twfe_catt_raceeth[["Hispanic"]])$summary[1:4]
  
  #Estimating CATT by Length unemployed post-unemployment
  twfe_catt_length <- PanelEstimate(sets = match_object, data = temp, moderator = "lessmed_length_unemp", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_length_moremed_sets <- print(twfe_catt_length[["0"]]$matched.sets)
  catt_length_lessmed_sets <- print(twfe_catt_length[["1"]]$matched.sets)
  
  catt_length_moremed <- results; catt_length_moremed$estimand <- "CATT_LENGTH_MORESMED"
  catt_length_lessmed <- results; catt_length_lessmed$estimand <- "CATT_LENGTH_LESSMED"
  
  catt_length_moremed[1,3:8] <- c((NROW(catt_length_moremed_sets)),
                                  (NROW(catt_length_moremed_sets$matched.set.size[catt_length_moremed_sets$matched.set.size == 0])),
                                  (NROW(catt_length_moremed_sets$matched.set.size[catt_length_moremed_sets$matched.set.size > 0])),
                                  (quantile(catt_length_moremed_sets$matched.set.size[catt_length_moremed_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_length_lessmed[1,3:8] <- c((NROW(catt_length_lessmed_sets)), 
                                  (NROW(catt_length_lessmed_sets$matched.set.size[catt_length_lessmed_sets$matched.set.size == 0])),
                                  (NROW(catt_length_lessmed_sets$matched.set.size[catt_length_lessmed_sets$matched.set.size > 0])),
                                  (quantile(catt_length_lessmed_sets$matched.set.size[catt_length_lessmed_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_length_moremed[1,9:12] <- summary(twfe_catt_length[["0"]])$summary[1:4]
  catt_length_lessmed[1,9:12] <- summary(twfe_catt_length[["1"]])$summary[1:4]
  
  #Estimating CATT by UI Receipt
  twfe_catt_ui <- PanelEstimate(sets = match_object, data = temp, moderator = "moderator_receipt", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_ui_no_sets <- print(twfe_catt_ui[["0"]]$matched.sets)
  catt_ui_yes_sets <- print(twfe_catt_ui[["1"]]$matched.sets)
  
  catt_ui_no <- results; catt_ui_no$estimand <- "CATT_UI_NO"
  catt_ui_yes <- results; catt_ui_yes$estimand <- "CATT_UI_YES"
  
  catt_ui_no[1,3:8] <- c((NROW(catt_ui_no_sets)),
                         (NROW(catt_ui_no_sets$matched.set.size[catt_ui_no_sets$matched.set.size == 0])),
                         (NROW(catt_ui_no_sets$matched.set.size[catt_ui_no_sets$matched.set.size > 0])),
                         (quantile(catt_ui_no_sets$matched.set.size[catt_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_ui_yes[1,3:8] <- c((NROW(catt_ui_yes_sets)),
                          (NROW(catt_ui_yes_sets$matched.set.size[catt_ui_yes_sets$matched.set.size == 0])),
                          (NROW(catt_ui_yes_sets$matched.set.size[catt_ui_yes_sets$matched.set.size > 0])),
                          (quantile(catt_ui_yes_sets$matched.set.size[catt_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_ui_no[1,9:12] <- summary(twfe_catt_ui[["0"]])$summary[1:4]
  catt_ui_yes[1,9:12] <- summary(twfe_catt_ui[["1"]])$summary[1:4]
  
  #Estimating CATT by UI Duration
  twfe_catt_ui_dur <- PanelEstimate(sets = match_object, data = temp, moderator = "ui_cat_dur", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_ui_dur_none_sets <- print(twfe_catt_ui_dur[["None"]]$matched.sets)
  catt_ui_dur_yeslessmed_sets <- print(twfe_catt_ui_dur[["YesLessMedDur"]]$matched.sets)
  catt_ui_dur_yesmoremed_sets <- print(twfe_catt_ui_dur[["YesMoreMedDur"]]$matched.sets)
  
  catt_ui_dur_none <- results; catt_ui_dur_none$estimand <- "CATT_UI_DUR_NONE"
  catt_ui_dur_yeslessmed <- results; catt_ui_dur_yeslessmed$estimand <- "CATT_UI_DUR_LESSMED"
  catt_ui_dur_yesmoremed <- results; catt_ui_dur_yesmoremed$estimand <- "CATT_UI_DUR_MOREMED"
  
  catt_ui_dur_none[1,3:8] <- c((NROW(catt_ui_dur_none_sets)),
                               (NROW(catt_ui_dur_none_sets$matched.set.size[catt_ui_dur_none_sets$matched.set.size == 0])),
                               (NROW(catt_ui_dur_none_sets$matched.set.size[catt_ui_dur_none_sets$matched.set.size > 0])),
                               (quantile(catt_ui_dur_none_sets$matched.set.size[catt_ui_dur_none_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_ui_dur_yeslessmed[1,3:8] <- c((NROW(catt_ui_dur_yeslessmed_sets)),
                                     (NROW(catt_ui_dur_yeslessmed_sets$matched.set.size[catt_ui_dur_yeslessmed_sets$matched.set.size == 0])),
                                     (NROW(catt_ui_dur_yeslessmed_sets$matched.set.size[catt_ui_dur_yeslessmed_sets$matched.set.size > 0])),
                                     (quantile(catt_ui_dur_yeslessmed_sets$matched.set.size[catt_ui_dur_yeslessmed_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_ui_dur_yesmoremed[1,3:8] <- c((NROW(catt_ui_dur_yesmoremed_sets)),
                                     (NROW(catt_ui_dur_yesmoremed_sets$matched.set.size[catt_ui_dur_yesmoremed_sets$matched.set.size == 0])),
                                     (NROW(catt_ui_dur_yesmoremed_sets$matched.set.size[catt_ui_dur_yesmoremed_sets$matched.set.size > 0])),
                                     (quantile(catt_ui_dur_yesmoremed_sets$matched.set.size[catt_ui_dur_yesmoremed_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_ui_dur_none[1,9:12] <- summary(twfe_catt_ui_dur[["None"]])$summary[1:4]
  catt_ui_dur_yeslessmed[1,9:12] <- summary(twfe_catt_ui_dur[["YesLessMedDur"]])$summary[1:4]
  catt_ui_dur_yesmoremed[1,9:12] <- summary(twfe_catt_ui_dur[["YesMoreMedDur"]])$summary[1:4]
  
  #Estimating CATT by UI Generosity
  twfe_catt_ui_gen <- PanelEstimate(sets = match_object, data = temp, moderator = "ui_cat_gen", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_ui_gen_none_sets <- print(twfe_catt_ui_gen[["None"]]$matched.sets)
  catt_ui_gen_yeslessmed_sets <- print(twfe_catt_ui_gen[["YesLessMedGen"]]$matched.sets)
  catt_ui_gen_yesmoremed_sets <- print(twfe_catt_ui_gen[["YesMoreMedGen"]]$matched.sets)
  
  catt_ui_gen_none <- results; catt_ui_gen_none$estimand <- "CATT_UI_gen_NONE"
  catt_ui_gen_yeslessmed <- results; catt_ui_gen_yeslessmed$estimand <- "CATT_UI_GEN_LESSMED"
  catt_ui_gen_yesmoremed <- results; catt_ui_gen_yesmoremed$estimand <- "CATT_UI_GEN_MOREMED"
  
  catt_ui_gen_none[1,3:8] <- c((NROW(catt_ui_gen_none_sets)),
                               (NROW(catt_ui_gen_none_sets$matched.set.size[catt_ui_gen_none_sets$matched.set.size == 0])),
                               (NROW(catt_ui_gen_none_sets$matched.set.size[catt_ui_gen_none_sets$matched.set.size > 0])),
                               (quantile(catt_ui_gen_none_sets$matched.set.size[catt_ui_gen_none_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_ui_gen_yeslessmed[1,3:8] <- c((NROW(catt_ui_gen_yeslessmed_sets)),
                                     (NROW(catt_ui_gen_yeslessmed_sets$matched.set.size[catt_ui_gen_yeslessmed_sets$matched.set.size == 0])),
                                     (NROW(catt_ui_gen_yeslessmed_sets$matched.set.size[catt_ui_gen_yeslessmed_sets$matched.set.size > 0])),
                                     (quantile(catt_ui_gen_yeslessmed_sets$matched.set.size[catt_ui_gen_yeslessmed_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_ui_gen_yesmoremed[1,3:8] <- c((NROW(catt_ui_gen_yesmoremed_sets)),
                                     (NROW(catt_ui_gen_yesmoremed_sets$matched.set.size[catt_ui_gen_yesmoremed_sets$matched.set.size == 0])),
                                     (NROW(catt_ui_gen_yesmoremed_sets$matched.set.size[catt_ui_gen_yesmoremed_sets$matched.set.size > 0])),
                                     (quantile(catt_ui_gen_yesmoremed_sets$matched.set.size[catt_ui_gen_yesmoremed_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_ui_gen_none[1,9:12] <- summary(twfe_catt_ui_gen[["None"]])$summary[1:4]
  catt_ui_gen_yeslessmed[1,9:12] <- summary(twfe_catt_ui_gen[["YesLessMedGen"]])$summary[1:4]
  catt_ui_gen_yesmoremed[1,9:12] <- summary(twfe_catt_ui_gen[["YesMoreMedGen"]])$summary[1:4]
  
  #Creating moderator variably for UI Receipt by Baseline EQ (Q1 and Q4)
  temp %>% mutate(ui_eq = case_when(moderator_receipt == 0 & baseline_eq_quartile == 1 ~ "EQ1_NoUI",
                                    moderator_receipt == 1 & baseline_eq_quartile == 1 ~ "EQ1_YesUI",
                                    moderator_receipt == 0 & baseline_eq_quartile == 4 ~ "EQ4_NoUI",
                                    moderator_receipt == 1 & baseline_eq_quartile == 4 ~ "EQ4_YesUI",
                                    TRUE ~ NA_character_)) -> temp
  
  #Estimating CATT by UI Receipt and Baseline EQ (Q1 and Q4)
  twfe_catt_eq_ui <- PanelEstimate(sets = match_object, data = temp, moderator = "ui_eq", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_q1eq_ui_no_sets <- print(twfe_catt_eq_ui[["EQ1_NoUI"]]$matched.sets)
  catt_q1eq_ui_yes_sets <- print(twfe_catt_eq_ui[["EQ1_YesUI"]]$matched.sets)
  catt_q4eq_ui_no_sets <- print(twfe_catt_eq_ui[["EQ4_NoUI"]]$matched.sets)
  catt_q4eq_ui_yes_sets <- print(twfe_catt_eq_ui[["EQ4_YesUI"]]$matched.sets)
  
  catt_q1eq_ui_no <- results; catt_q1eq_ui_no$estimand <- "CATT_EQ_Q1_UI_NO"
  catt_q1eq_ui_yes <- results; catt_q1eq_ui_yes$estimand <- "CATT_EQ_Q1_UI_YES"
  catt_q4eq_ui_no <- results; catt_q4eq_ui_no$estimand <- "CATT_EQ_Q4_UI_NO"
  catt_q4eq_ui_yes <- results; catt_q4eq_ui_yes$estimand <- "CATT_EQ_Q4_UI_YES"
  
  catt_q1eq_ui_no[1,3:8] <- c((NROW(catt_q1eq_ui_no_sets)),
                              (NROW(catt_q1eq_ui_no_sets$matched.set.size[catt_q1eq_ui_no_sets$matched.set.size == 0])),
                              (NROW(catt_q1eq_ui_no_sets$matched.set.size[catt_q1eq_ui_no_sets$matched.set.size > 0])),
                              (quantile(catt_q1eq_ui_no_sets$matched.set.size[catt_q1eq_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_q1eq_ui_yes[1,3:8] <- c((NROW(catt_q1eq_ui_yes_sets)),
                               (NROW(catt_q1eq_ui_yes_sets$matched.set.size[catt_q1eq_ui_yes_sets$matched.set.size == 0])),
                               (NROW(catt_q1eq_ui_yes_sets$matched.set.size[catt_q1eq_ui_yes_sets$matched.set.size > 0])),
                               (quantile(catt_q1eq_ui_yes_sets$matched.set.size[catt_q1eq_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_q4eq_ui_no[1,3:8] <- c((NROW(catt_q4eq_ui_no_sets)),
                              (NROW(catt_q4eq_ui_no_sets$matched.set.size[catt_q4eq_ui_no_sets$matched.set.size == 0])),
                              (NROW(catt_q4eq_ui_no_sets$matched.set.size[catt_q4eq_ui_no_sets$matched.set.size > 0])),
                              (quantile(catt_q4eq_ui_no_sets$matched.set.size[catt_q4eq_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_q4eq_ui_yes[1,3:8] <- c((NROW(catt_q4eq_ui_yes_sets)),
                               (NROW(catt_q4eq_ui_yes_sets$matched.set.size[catt_q4eq_ui_yes_sets$matched.set.size == 0])),
                               (NROW(catt_q4eq_ui_yes_sets$matched.set.size[catt_q4eq_ui_yes_sets$matched.set.size > 0])),
                               (quantile(catt_q4eq_ui_yes_sets$matched.set.size[catt_q4eq_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_q1eq_ui_no[1,9:12] <- summary(twfe_catt_eq_ui[["EQ1_NoUI"]])$summary[1:4]
  catt_q1eq_ui_yes[1,9:12] <- summary(twfe_catt_eq_ui[["EQ1_YesUI"]])$summary[1:4]
  catt_q4eq_ui_no[1,9:12] <- summary(twfe_catt_eq_ui[["EQ4_NoUI"]])$summary[1:4]
  catt_q4eq_ui_yes[1,9:12] <- summary(twfe_catt_eq_ui[["EQ4_YesUI"]])$summary[1:4]
  
  #Creating moderator variably for UI Receipt by Sex
  temp %>% mutate(ui_sex = case_when(moderator_receipt == 0 & male == 1 ~ "Male_NoUI",
                                     moderator_receipt == 1 & male == 1 ~ "Male_YesUI",
                                     moderator_receipt == 0 & male == 0 ~ "Female_NoUI",
                                     moderator_receipt == 1 & male == 0 ~ "Female_YesUI",
                                     TRUE ~ NA_character_)) -> temp
  
  #Estimating CATT by UI Receipt and Sex
  twfe_catt_sex_ui <- PanelEstimate(sets = match_object, data = temp, moderator = "ui_sex", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_f_ui_no_sets <- print(twfe_catt_sex_ui[["Female_NoUI"]]$matched.sets)
  catt_f_ui_yes_sets <- print(twfe_catt_sex_ui[["Female_YesUI"]]$matched.sets)
  catt_m_ui_no_sets <- print(twfe_catt_sex_ui[["Male_NoUI"]]$matched.sets)
  catt_m_ui_yes_sets <- print(twfe_catt_sex_ui[["Male_YesUI"]]$matched.sets)
  
  catt_f_ui_no <- results; catt_f_ui_no$estimand <- "CATT_FEMALE_UI_NO"
  catt_f_ui_yes <- results; catt_f_ui_yes$estimand <- "CATT_FEMALE_UI_YES"
  catt_m_ui_no <- results; catt_m_ui_no$estimand <- "CATT_MALE_UI_NO"
  catt_m_ui_yes <- results; catt_m_ui_yes$estimand <- "CATT_MALE_UI_YES"
  
  catt_f_ui_no[1,3:8] <- c((NROW(catt_f_ui_no_sets)),
                           (NROW(catt_f_ui_no_sets$matched.set.size[catt_f_ui_no_sets$matched.set.size == 0])),
                           (NROW(catt_f_ui_no_sets$matched.set.size[catt_f_ui_no_sets$matched.set.size > 0])),
                           (quantile(catt_f_ui_no_sets$matched.set.size[catt_f_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_f_ui_yes[1,3:8] <- c((NROW(catt_f_ui_yes_sets)),
                            (NROW(catt_f_ui_yes_sets$matched.set.size[catt_f_ui_yes_sets$matched.set.size == 0])),
                            (NROW(catt_f_ui_yes_sets$matched.set.size[catt_f_ui_yes_sets$matched.set.size > 0])),
                            (quantile(catt_f_ui_yes_sets$matched.set.size[catt_f_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_m_ui_no[1,3:8] <- c((NROW(catt_m_ui_no_sets)),
                           (NROW(catt_m_ui_no_sets$matched.set.size[catt_m_ui_no_sets$matched.set.size == 0])),
                           (NROW(catt_m_ui_no_sets$matched.set.size[catt_m_ui_no_sets$matched.set.size > 0])),
                           (quantile(catt_m_ui_no_sets$matched.set.size[catt_m_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_m_ui_yes[1,3:8] <- c((NROW(catt_m_ui_yes_sets)),
                            (NROW(catt_m_ui_yes_sets$matched.set.size[catt_m_ui_yes_sets$matched.set.size == 0])),
                            (NROW(catt_m_ui_yes_sets$matched.set.size[catt_m_ui_yes_sets$matched.set.size > 0])),
                            (quantile(catt_m_ui_yes_sets$matched.set.size[catt_m_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_f_ui_no[1,9:12] <- summary(twfe_catt_sex_ui[["Female_NoUI"]])$summary[1:4]
  catt_f_ui_yes[1,9:12] <- summary(twfe_catt_sex_ui[["Female_YesUI"]])$summary[1:4]
  catt_m_ui_no[1,9:12] <- summary(twfe_catt_sex_ui[["Male_NoUI"]])$summary[1:4]
  catt_m_ui_yes[1,9:12] <- summary(twfe_catt_sex_ui[["Male_YesUI"]])$summary[1:4]
  
  #Creating moderator variably for UI Receipt by Race and Ethnicity
  temp %>% mutate(ui_raceeth = case_when(moderator_receipt == 0 & race_eth == "NHWhite" ~ "NHWhite_NoUI",
                                         moderator_receipt == 1 & race_eth == "NHWhite" ~ "NHWhite_YesUI",
                                         moderator_receipt == 0 & race_eth == "NHBlack" ~ "NHBlack_NoUI",
                                         moderator_receipt == 1 & race_eth == "NHBlack" ~ "NHBlack_YesUI",
                                         moderator_receipt == 0 & race_eth == "NHOther" ~ "NHOther_NoUI",
                                         moderator_receipt == 1 & race_eth == "NHOther" ~ "NHOther_YesUI",
                                         moderator_receipt == 0 & race_eth == "Hispanic" ~ "Hispanic_NoUI",
                                         moderator_receipt == 1 & race_eth == "Hispanic" ~ "Hispanic_YesUI",
                                         TRUE ~ NA_character_)) -> temp
  
  #Estimating CATT by UI Receipt and Race and Ethnicity
  twfe_catt_raceeth_ui <- PanelEstimate(sets = match_object, data = temp, moderator = "ui_raceeth", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_nhwhite_ui_no_sets <- print(twfe_catt_raceeth_ui[["NHWhite_NoUI"]]$matched.sets)
  catt_nhwhite_ui_yes_sets <- print(twfe_catt_raceeth_ui[["NHWhite_YesUI"]]$matched.sets)
  catt_nhblack_ui_no_sets <- print(twfe_catt_raceeth_ui[["NHBlack_NoUI"]]$matched.sets)
  catt_nhblack_ui_yes_sets <- print(twfe_catt_raceeth_ui[["NHBlack_YesUI"]]$matched.sets)
  catt_nhother_ui_no_sets <- print(twfe_catt_raceeth_ui[["NHOther_NoUI"]]$matched.sets)
  catt_nhother_ui_yes_sets <- print(twfe_catt_raceeth_ui[["NHOther_YesUI"]]$matched.sets)
  catt_hispanic_ui_no_sets <- print(twfe_catt_raceeth_ui[["Hispanic_NoUI"]]$matched.sets)
  catt_hispanic_ui_yes_sets <- print(twfe_catt_raceeth_ui[["Hispanic_YesUI"]]$matched.sets)
  
  catt_nhwhite_ui_no <- results; catt_nhwhite_ui_no$estimand <- "CATT_NHWHITE_UI_NO"
  catt_nhwhite_ui_yes <- results; catt_nhwhite_ui_yes$estimand <- "CATT_NHWHITE_UI_YES"
  catt_nhblack_ui_no <- results; catt_nhblack_ui_no$estimand <- "CATT_NHBLACK_UI_NO"
  catt_nhblack_ui_yes <- results; catt_nhblack_ui_yes$estimand <- "CATT_NHBLACK_UI_YES"
  catt_nhother_ui_no <- results; catt_nhother_ui_no$estimand <- "CATT_NHOTHER_UI_NO"
  catt_nhother_ui_yes <- results; catt_nhother_ui_yes$estimand <- "CATT_NHOTHER_UI_YES"
  catt_hispanic_ui_no <- results; catt_hispanic_ui_no$estimand <- "CATT_HISPANIC_UI_NO"
  catt_hispanic_ui_yes <- results; catt_hispanic_ui_yes$estimand <- "CATT_HISPANIC_UI_YES"
  
  catt_nhwhite_ui_no[1,3:8] <- c((NROW(catt_nhwhite_ui_no_sets)),
                                 (NROW(catt_nhwhite_ui_no_sets$matched.set.size[catt_nhwhite_ui_no_sets$matched.set.size == 0])),
                                 (NROW(catt_nhwhite_ui_no_sets$matched.set.size[catt_nhwhite_ui_no_sets$matched.set.size > 0])),
                                 (quantile(catt_nhwhite_ui_no_sets$matched.set.size[catt_nhwhite_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_nhwhite_ui_yes[1,3:8] <- c((NROW(catt_nhwhite_ui_yes_sets)),
                                  (NROW(catt_nhwhite_ui_yes_sets$matched.set.size[catt_nhwhite_ui_yes_sets$matched.set.size == 0])),
                                  (NROW(catt_nhwhite_ui_yes_sets$matched.set.size[catt_nhwhite_ui_yes_sets$matched.set.size > 0])),
                                  (quantile(catt_nhwhite_ui_yes_sets$matched.set.size[catt_nhwhite_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_nhblack_ui_no[1,3:8] <- c((NROW(catt_nhblack_ui_no_sets)),
                                 (NROW(catt_nhblack_ui_no_sets$matched.set.size[catt_nhblack_ui_no_sets$matched.set.size == 0])),
                                 (NROW(catt_nhblack_ui_no_sets$matched.set.size[catt_nhblack_ui_no_sets$matched.set.size > 0])),
                                 (quantile(catt_nhblack_ui_no_sets$matched.set.size[catt_nhblack_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_nhblack_ui_yes[1,3:8] <- c((NROW(catt_nhblack_ui_yes_sets)),
                                  (NROW(catt_nhblack_ui_yes_sets$matched.set.size[catt_nhblack_ui_yes_sets$matched.set.size == 0])),
                                  (NROW(catt_nhblack_ui_yes_sets$matched.set.size[catt_nhblack_ui_yes_sets$matched.set.size > 0])),
                                  (quantile(catt_nhblack_ui_yes_sets$matched.set.size[catt_nhblack_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_nhother_ui_no[1,3:8] <- c((NROW(catt_nhother_ui_no_sets)),
                                 (NROW(catt_nhother_ui_no_sets$matched.set.size[catt_nhother_ui_no_sets$matched.set.size == 0])),
                                 (NROW(catt_nhother_ui_no_sets$matched.set.size[catt_nhother_ui_no_sets$matched.set.size > 0])),
                                 (quantile(catt_nhother_ui_no_sets$matched.set.size[catt_nhother_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_nhother_ui_yes[1,3:8] <- c((NROW(catt_nhother_ui_yes_sets)),
                                  (NROW(catt_nhother_ui_yes_sets$matched.set.size[catt_nhother_ui_yes_sets$matched.set.size == 0])),
                                  (NROW(catt_nhother_ui_yes_sets$matched.set.size[catt_nhother_ui_yes_sets$matched.set.size > 0])),
                                  (quantile(catt_nhother_ui_yes_sets$matched.set.size[catt_nhother_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_hispanic_ui_no[1,3:8] <- c((NROW(catt_hispanic_ui_no_sets)),
                                  (NROW(catt_hispanic_ui_no_sets$matched.set.size[catt_hispanic_ui_no_sets$matched.set.size == 0])),
                                  (NROW(catt_hispanic_ui_no_sets$matched.set.size[catt_hispanic_ui_no_sets$matched.set.size > 0])),
                                  (quantile(catt_hispanic_ui_no_sets$matched.set.size[catt_hispanic_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_hispanic_ui_yes[1,3:8] <- c((NROW(catt_hispanic_ui_yes_sets)),
                                   (NROW(catt_hispanic_ui_yes_sets$matched.set.size[catt_hispanic_ui_yes_sets$matched.set.size == 0])),
                                   (NROW(catt_hispanic_ui_yes_sets$matched.set.size[catt_hispanic_ui_yes_sets$matched.set.size > 0])),
                                   (quantile(catt_hispanic_ui_yes_sets$matched.set.size[catt_hispanic_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_nhwhite_ui_no[1,9:12] <- summary(twfe_catt_raceeth_ui[["NHWhite_NoUI"]])$summary[1:4]
  catt_nhwhite_ui_yes[1,9:12] <- summary(twfe_catt_raceeth_ui[["NHWhite_YesUI"]])$summary[1:4]
  catt_nhblack_ui_no[1,9:12] <- summary(twfe_catt_raceeth_ui[["NHBlack_NoUI"]])$summary[1:4]
  catt_nhblack_ui_yes[1,9:12] <- summary(twfe_catt_raceeth_ui[["NHBlack_YesUI"]])$summary[1:4]
  catt_nhother_ui_no[1,9:12] <- summary(twfe_catt_raceeth_ui[["NHOther_NoUI"]])$summary[1:4]
  catt_nhother_ui_yes[1,9:12] <- summary(twfe_catt_raceeth_ui[["NHOther_YesUI"]])$summary[1:4]
  catt_hispanic_ui_no[1,9:12] <- summary(twfe_catt_raceeth_ui[["Hispanic_NoUI"]])$summary[1:4]
  catt_hispanic_ui_yes[1,9:12] <- summary(twfe_catt_raceeth_ui[["Hispanic_YesUI"]])$summary[1:4]
  
  #Creating moderator variably for UI Receipt by Length Unemployed Post-Unemployment
  temp %>% mutate(ui_length = case_when(moderator_receipt == 0 & lessmed_length_unemp == 1 ~ "LengthLessMed_NoUI",
                                        moderator_receipt == 1 & lessmed_length_unemp == 1 ~ "LengthLessMed_YesUI",
                                        moderator_receipt == 0 & lessmed_length_unemp == 0 ~ "LengthMoreMed_NoUI",
                                        moderator_receipt == 1 & lessmed_length_unemp == 0 ~ "LengthMoreMed_YesUI",
                                        TRUE ~ NA_character_)) -> temp
  
  #Estimating CATT by UI Receipt and Length Unemployed Post-Unemployment
  twfe_catt_length_ui <- PanelEstimate(sets = match_object, data = temp, moderator = "ui_length", se.method = "bootstrap", confidence.level = .95)
  
  #Populating CATT-specific results variables
  catt_lessmedlength_ui_no_sets <- print(twfe_catt_length_ui[["LengthLessMed_NoUI"]]$matched.sets)
  catt_lessmedlength_ui_yes_sets <- print(twfe_catt_length_ui[["LengthLessMed_YesUI"]]$matched.sets)
  catt_moremedlength_ui_no_sets <- print(twfe_catt_length_ui[["LengthMoreMed_NoUI"]]$matched.sets)
  catt_moremedlength_ui_yes_sets <- print(twfe_catt_length_ui[["LengthMoreMed_YesUI"]]$matched.sets)
  
  catt_lessmedlength_ui_no <- results; catt_lessmedlength_ui_no$estimand <- "CATT_LESSMEDLENGTH_UI_NO"
  catt_lessmedlength_ui_yes <- results; catt_lessmedlength_ui_yes$estimand <- "CATT_LESSMEDLENGTH_UI_YES"
  catt_moremedlength_ui_no <- results; catt_moremedlength_ui_no$estimand <- "CATT_MOREMEDLENGTH_UI_NO"
  catt_moremedlength_ui_yes <- results; catt_moremedlength_ui_yes$estimand <- "CATT_MOREMEDLENGTH_UI_YES"
  
  catt_lessmedlength_ui_no[1,3:8] <- c((NROW(catt_lessmedlength_ui_no_sets)),
                                       (NROW(catt_lessmedlength_ui_no_sets$matched.set.size[catt_lessmedlength_ui_no_sets$matched.set.size == 0])),
                                       (NROW(catt_lessmedlength_ui_no_sets$matched.set.size[catt_lessmedlength_ui_no_sets$matched.set.size > 0])),
                                       (quantile(catt_lessmedlength_ui_no_sets$matched.set.size[catt_lessmedlength_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_lessmedlength_ui_yes[1,3:8] <- c((NROW(catt_lessmedlength_ui_yes_sets)),
                                        (NROW(catt_lessmedlength_ui_yes_sets$matched.set.size[catt_lessmedlength_ui_yes_sets$matched.set.size == 0])),
                                        (NROW(catt_lessmedlength_ui_yes_sets$matched.set.size[catt_lessmedlength_ui_yes_sets$matched.set.size > 0])),
                                        (quantile(catt_lessmedlength_ui_yes_sets$matched.set.size[catt_lessmedlength_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_moremedlength_ui_no[1,3:8] <- c((NROW(catt_moremedlength_ui_no_sets)),
                                       (NROW(catt_moremedlength_ui_no_sets$matched.set.size[catt_moremedlength_ui_no_sets$matched.set.size == 0])),
                                       (NROW(catt_moremedlength_ui_no_sets$matched.set.size[catt_moremedlength_ui_no_sets$matched.set.size > 0])),
                                       (quantile(catt_moremedlength_ui_no_sets$matched.set.size[catt_moremedlength_ui_no_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  catt_moremedlength_ui_yes[1,3:8] <- c((NROW(catt_moremedlength_ui_yes_sets)),
                                        (NROW(catt_moremedlength_ui_yes_sets$matched.set.size[catt_moremedlength_ui_yes_sets$matched.set.size == 0])),
                                        (NROW(catt_moremedlength_ui_yes_sets$matched.set.size[catt_moremedlength_ui_yes_sets$matched.set.size > 0])),
                                        (quantile(catt_moremedlength_ui_yes_sets$matched.set.size[catt_moremedlength_ui_yes_sets$matched.set.size > 0], c(0.5, 0.25, 0.75))))
  
  catt_lessmedlength_ui_no[1,9:12] <- summary(twfe_catt_length_ui[["LengthLessMed_NoUI"]])$summary[1:4]
  catt_lessmedlength_ui_yes[1,9:12] <- summary(twfe_catt_length_ui[["LengthLessMed_YesUI"]])$summary[1:4]
  catt_moremedlength_ui_no[1,9:12] <- summary(twfe_catt_length_ui[["LengthMoreMed_NoUI"]])$summary[1:4]
  catt_moremedlength_ui_yes[1,9:12] <- summary(twfe_catt_length_ui[["LengthMoreMed_YesUI"]])$summary[1:4]
  
  #Appending all results
  results_temp <- rbind(att, 
                        catt_eq_q1, catt_eq_q2, catt_eq_q3, catt_eq_q4,
                        catt_sex_f, catt_sex_m,
                        catt_raceeth_nhwhite, catt_raceeth_nhblack, catt_raceeth_nhother, catt_raceeth_hispanic,
                        catt_length_lessmed, catt_length_moremed,
                        catt_ui_no, catt_ui_yes,
                        catt_ui_dur_none, catt_ui_dur_yeslessmed, catt_ui_dur_yesmoremed,
                        catt_ui_gen_none, catt_ui_gen_yeslessmed, catt_ui_gen_yesmoremed,
                        catt_q1eq_ui_no, catt_q1eq_ui_yes,
                        catt_q4eq_ui_no, catt_q4eq_ui_yes,
                        catt_m_ui_no, catt_m_ui_yes,
                        catt_f_ui_no, catt_f_ui_yes,
                        catt_nhwhite_ui_no, catt_nhwhite_ui_yes,
                        catt_nhblack_ui_no, catt_nhblack_ui_yes,
                        catt_nhother_ui_no, catt_nhother_ui_yes,
                        catt_hispanic_ui_no, catt_hispanic_ui_yes,
                        catt_lessmedlength_ui_no, catt_lessmedlength_ui_yes,
                        catt_moremedlength_ui_no, catt_moremedlength_ui_yes)
  
  assign(paste0(dataset, "_results"), results_temp)
  
  #Saving results
  qsave(results_temp, paste0(directory, dataset, "_results_220823.qs"))
  
  #Clearing objects for next dataset
  rm(dataset, temp, results, 
     match_object, mset, set_details, 
     att, twfe_att, 
     twfe_catt_eq, catt_eq_q1_sets, catt_eq_q2_sets, catt_eq_q3_sets, catt_eq_q4_sets, 
     catt_eq_q1, catt_eq_q2, catt_eq_q3, catt_eq_q4, 
     twfe_catt_sex, catt_sex_f_sets, catt_sex_m_sets, catt_sex_f, catt_sex_m, 
     twfe_catt_raceeth, catt_raceeth_nhwhite_sets, catt_raceeth_nhblack_sets, catt_raceeth_nhother_sets, catt_raceeth_hispanic_sets, 
     catt_raceeth_nhwhite, catt_raceeth_nhblack, catt_raceeth_nhother, catt_raceeth_hispanic, 
     twfe_catt_length, catt_length_moremed_sets, catt_length_lessmed_sets, catt_length_moremed, catt_length_lessmed, 
     twfe_catt_ui, catt_ui_no_sets, catt_ui_yes_sets, catt_ui_no, catt_ui_yes, 
     twfe_catt_ui_dur, catt_ui_dur_none_sets, catt_ui_dur_yeslessmed_sets, catt_ui_dur_yesmoremed_sets, 
     catt_ui_dur_none, catt_ui_dur_yeslessmed, catt_ui_dur_yesmoremed, 
     twfe_catt_ui_gen, catt_ui_gen_none_sets, catt_ui_gen_yeslessmed_sets, catt_ui_gen_yesmoremed_sets, 
     catt_ui_gen_none, catt_ui_gen_yeslessmed, catt_ui_gen_yesmoremed,
     twfe_catt_eq_ui, catt_q1eq_ui_no_sets, catt_q1eq_ui_yes_sets, catt_q4eq_ui_no_sets, catt_q4eq_ui_yes_sets,
     catt_q1eq_ui_no, catt_q1eq_ui_yes, catt_q4eq_ui_no, catt_q4eq_ui_yes,
     twfe_catt_sex_ui, catt_f_ui_no_sets, catt_f_ui_yes_sets, catt_m_ui_no_sets, catt_m_ui_yes_sets,
     catt_f_ui_no, catt_f_ui_yes, catt_m_ui_no, catt_m_ui_yes,
     twfe_catt_raceeth_ui, 
     catt_nhwhite_ui_no_sets, catt_nhwhite_ui_yes_sets, catt_nhblack_ui_no_sets, catt_nhblack_ui_yes_sets, 
     catt_nhother_ui_no_sets, catt_nhother_ui_yes_sets, catt_hispanic_ui_no_sets, catt_hispanic_ui_yes_sets,
     catt_nhwhite_ui_no, catt_nhwhite_ui_yes, catt_nhblack_ui_no, catt_nhblack_ui_yes, 
     catt_nhother_ui_no, catt_nhother_ui_yes, catt_hispanic_ui_no, catt_hispanic_ui_yes,
     twfe_catt_length_ui, catt_lessmedlength_ui_no_sets, catt_lessmedlength_ui_yes_sets, catt_moremedlength_ui_no_sets, catt_moremedlength_ui_yes_sets,
     catt_lessmedlength_ui_no, catt_lessmedlength_ui_yes, catt_moremedlength_ui_no, catt_moremedlength_ui_yes,
     results_temp)  
}

#### Step 4 - Creating a final summary dataset with CC results and MI results applying Rubin's Rules ####
#Creating an initial results shell with CC results included
final_results <- cc_data_results; cc_data_results$data <- "CC"

#Stitching all MI results together
for (mi in 1:35) {
  if (mi == 1) {
    mi_results <- mi1_data_results
  } else {
    eval(parse(text = paste0("mi_results <- rbind(mi_results, mi", mi, "_data_results)")))
  }
  rm(mi)
}

#Creating a reduced dataset with pooled mean estimates
mi_results %>% 
  mutate(N_sets = as.numeric(N_sets),
         N_sets_unmatched = as.numeric(N_sets_unmatched),
         N_sets_matched = as.numeric(N_sets_matched),
         set_size_med = as.numeric(set_size_med),
         set_size_1qr = as.numeric(set_size_1qr),
         set_size_3qr = as.numeric(set_size_3qr)) %>%
  group_by(estimand) %>%
  mutate(N_sets = mean(N_sets, na.rm=T),
         N_sets_unmatched = mean(N_sets_unmatched, na.rm=T),
         N_sets_matched = mean(N_sets_matched, na.rm=T),
         set_size_med = mean(set_size_med, na.rm=T),
         set_size_1qr = mean(set_size_1qr, na.rm=T),
         set_size_3qr = mean(set_size_3qr, na.rm=T),
         est_mu = mean(est, na.rm=T),
         se_within = (sum(se^2, na.rm=T))/35,
         se_between = (sum(((est - est_mu)^2), na.rm=T))/34,
         se_pooled = sqrt((se_within + se_between + (se_between/35))),
         lci_pooled = est_mu - 1.96*se_pooled,
         uci_pooled = est_mu + 1.96*se_pooled) %>%
  ungroup() %>% 
  mutate(data = "MI") -> mi_results

mi_results$est <- NULL; mi_results$se <- NULL; mi_results$lci <- NULL; mi_results$uci <- NULL
mi_results <- distinct(mi_results)
final_results_mi <- mi_results[,c("data", "estimand", 
                                  "N_sets", "N_sets_unmatched", "N_sets_matched", 
                                  "set_size_med", "set_size_1qr", "set_size_3qr", 
                                  "est_mu", "se_pooled", "lci_pooled", "uci_pooled")]
names(final_results_mi) <- names(final_results)
final_results <- rbind(final_results, final_results_mi)
rm(final_results_mi)

#Saving results
qsave(final_results, paste0(directory, "final_results_220823.qs"))

#### Step 5 - Creating summary characteristics for CC and MI datasets ####
#Creating complete case summary
cc_data %>% 
  select(unique_id, year, prepost, outcome, Y_base, moderator_receipt, age, male, black, other, hispanic, mo_unemp_until_emp) %>%
  mutate(age = (age*11.40922)+39.7,
         change_eq = outcome - Y_base) -> cc_data #Un-standardizing age ((Mu: 39.7, SD: 11.40922))
cc_CE <- cc_data[cc_data$prepost == "CE" & !is.na(cc_data$prepost), ]
cc_CE_inds <- cc_CE %>% select(unique_id, black, other, hispanic, male) %>% distinct()
cc_RU <- cc_data[cc_data$prepost == "RU" & !is.na(cc_data$prepost), ]
cc_RU_inds <- cc_RU %>% select(unique_id, black, other, hispanic, male) %>% distinct()

cc_CE_descript <- data.frame(data = "CC",
                             strata = "CE", 
                             n_ind = NROW(unique(cc_CE$unique_id)), n_obs = NROW(cc_CE),
                             age_mu = mean(cc_CE$age, na.rm = T), 
                             age_sd = sd(cc_CE$age, na.rm=T), 
                             age_na = NROW(cc_CE[is.na(cc_CE$age),]),
                             gender_male = NROW(cc_CE_inds[cc_CE_inds$male == 1 & !is.na(cc_CE_inds$male),]), 
                             gender_male_perc = mean(cc_CE_inds$male, na.rm=T),
                             gender_male_na = NROW(cc_CE_inds[is.na(cc_CE_inds$male),]),
                             race_nhwhite_n = NROW(cc_CE_inds[cc_CE_inds$black == 0 & cc_CE_inds$other == 0 & cc_CE_inds$hispanic == 0, ]), 
                             race_nhwhite_perc = (1-mean(cc_CE_inds$black, na.rm=T) - mean(cc_CE_inds$other, na.rm=T) - mean(cc_CE_inds$hispanic, na.rm=T))*100,
                             race_nhblack_n = NROW(cc_CE_inds[cc_CE_inds$black == 1,]),
                             race_nhblack_perc = mean(cc_CE_inds$black, na.rm=T),
                             race_nhother_n = NROW(cc_CE_inds[cc_CE_inds$other == 1,]),
                             race_nhother_perc = mean(cc_CE_inds$other, na.rm=T),
                             race_hispanic_n = NROW(cc_CE_inds[cc_CE_inds$hispanic == 1,]),
                             race_hispanic_perc = mean(cc_CE_inds$hispanic, na.rm=T),
                             base_eq_mu = mean(cc_CE$Y_base, na.rm = T), 
                             base_eq_sd = sd(cc_CE$Y_base, na.rm=T), 
                             base_eq_na = NROW(cc_CE[is.na(cc_CE$Y_base),]),
                             eq_change_mu = mean(cc_CE$change_eq, na.rm = T), 
                             eq_change_sd = sd(cc_CE$change_eq, na.rm=T), 
                             eq_change_na = NROW(cc_CE[is.na(cc_CE$change_eq),]),
                             moderator_n = NA_real_, moderator_perc = NA_real_, moderator_na = NA_real_,
                             mo_unemp_until_emp_mu = NA_real_, mo_unemp_until_emp_sd = NA_real_, mo_unemp_until_emp_na = NA_real_)
cc_RU_descript <- data.frame(data = "CC",
                             strata = "RU", 
                             n_ind = NROW(unique(cc_RU$unique_id)), n_obs = NROW(cc_RU),
                             age_mu = mean(cc_RU$age, na.rm = T), 
                             age_sd = sd(cc_RU$age, na.rm=T), 
                             age_na = NROW(cc_RU[is.na(cc_RU$age),]),
                             gender_male = NROW(cc_RU_inds[cc_RU_inds$male == 1 & !is.na(cc_RU_inds$male),]), 
                             gender_male_perc = mean(cc_RU_inds$male, na.rm=T),
                             gender_male_na = NROW(cc_RU_inds[is.na(cc_RU_inds$male),]),
                             race_nhwhite_n = NROW(cc_RU_inds[cc_RU_inds$black == 0 & cc_RU_inds$other == 0 & cc_RU_inds$hispanic == 0, ]), 
                             race_nhwhite_perc = (1-mean(cc_RU_inds$black, na.rm=T) - mean(cc_RU_inds$other, na.rm=T) - mean(cc_RU_inds$hispanic, na.rm=T))*100,
                             race_nhblack_n = NROW(cc_RU_inds[cc_RU_inds$black == 1,]),
                             race_nhblack_perc = mean(cc_RU_inds$black, na.rm=T),
                             race_nhother_n = NROW(cc_RU_inds[cc_RU_inds$other == 1,]),
                             race_nhother_perc = mean(cc_RU_inds$other, na.rm=T),
                             race_hispanic_n = NROW(cc_RU_inds[cc_RU_inds$hispanic == 1,]),
                             race_hispanic_perc = mean(cc_RU_inds$hispanic, na.rm=T),
                             base_eq_mu = mean(cc_RU$Y_base, na.rm = T), 
                             base_eq_sd = sd(cc_RU$Y_base, na.rm=T), 
                             base_eq_na = NROW(cc_RU[is.na(cc_RU$Y_base),]),
                             eq_change_mu = mean(cc_RU$change_eq, na.rm = T), 
                             eq_change_sd = sd(cc_RU$change_eq, na.rm=T), 
                             eq_change_na = NROW(cc_RU[is.na(cc_RU$change_eq),]),
                             moderator_n = NROW(cc_RU[cc_RU$moderator_receipt == 1 & !is.na(cc_RU$moderator_receipt),]), 
                             moderator_perc = mean(cc_RU$moderator_receipt, na.rm = T), 
                             moderator_na = NROW(cc_RU[is.na(cc_RU$moderator_receipt),]),
                             mo_unemp_until_emp_mu = mean(cc_RU$mo_unemp_until_emp, na.rm = T), 
                             mo_unemp_until_emp_sd = sd(cc_RU$mo_unemp_until_emp, na.rm=T), 
                             mo_unemp_until_emp_na = NROW(cc_RU[is.na(cc_RU$mo_unemp_until_emp),]))
cc_descriptives <- rbind(cc_CE_descript, cc_RU_descript)
rm(cc_CE, cc_CE_inds, cc_CE_descript, cc_RU, cc_RU_inds, cc_RU_descript)

descriptives <- cc_descriptives

#Creating MI summary
for (mi in 1:35) {
  print(paste0("MI: ", mi))
  temp <- eval(parse(text = paste0("mi", mi, "_data")))
  
  temp %>% 
    select(unique_id, year, prepost, outcome, Y_base, moderator_receipt, age, male, black, other, hispanic, mo_unemp_until_emp) %>%
    mutate(age = (age*11.40922)+39.7,
           change_eq = outcome - Y_base) -> temp #Un-standardizing age ((Mu: 39.7, SD: 11.40922))
  
  tempCE <- temp[temp$prepost == "CE" & !is.na(temp$prepost), ]
  tempCE_inds <- tempCE %>% select(unique_id, black, other, hispanic, male) %>% distinct()
  tempRU <- temp[temp$prepost == "RU" & !is.na(temp$prepost), ]
  tempRU_inds <- tempRU %>% select(unique_id, black, other, hispanic, male) %>% distinct()
  
  tempCE_descript <- data.frame(data = paste0("MI", mi),
                                strata = "CE", 
                                n_ind = NROW(unique(tempCE$unique_id)), n_obs = NROW(tempCE),
                                age_mu = mean(tempCE$age, na.rm = T), 
                                age_sd = sd(tempCE$age, na.rm=T), 
                                age_na = NROW(tempCE[is.na(tempCE$age),]),
                                gender_male = NROW(tempCE_inds[tempCE_inds$male == 1 & !is.na(tempCE_inds$male),]), 
                                gender_male_perc = mean(tempCE_inds$male, na.rm=T),
                                gender_male_na = NROW(tempCE_inds[is.na(tempCE_inds$male),]),
                                race_nhwhite_n = NROW(tempCE_inds[tempCE_inds$black == 0 & tempCE_inds$other == 0 & tempCE_inds$hispanic == 0, ]), 
                                race_nhwhite_perc = (1-mean(tempCE_inds$black, na.rm=T) - mean(tempCE_inds$other, na.rm=T) - mean(tempCE_inds$hispanic, na.rm=T))*100,
                                race_nhblack_n = NROW(tempCE_inds[tempCE_inds$black == 1,]),
                                race_nhblack_perc = mean(tempCE_inds$black, na.rm=T),
                                race_nhother_n = NROW(tempCE_inds[tempCE_inds$other == 1,]),
                                race_nhother_perc = mean(tempCE_inds$other, na.rm=T),
                                race_hispanic_n = NROW(tempCE_inds[tempCE_inds$hispanic == 1,]),
                                race_hispanic_perc = mean(tempCE_inds$hispanic, na.rm=T),
                                base_eq_mu = mean(tempCE$Y_base, na.rm = T), 
                                base_eq_sd = sd(tempCE$Y_base, na.rm=T), 
                                base_eq_na = NROW(tempCE[is.na(tempCE$Y_base),]),
                                eq_change_mu = mean(tempCE$change_eq, na.rm = T), 
                                eq_change_sd = sd(tempCE$change_eq, na.rm=T), 
                                eq_change_na = NROW(tempCE[is.na(tempCE$change_eq),]),
                                moderator_n = NA_real_, moderator_perc = NA_real_, moderator_na = NA_real_,
                                mo_unemp_until_emp_mu = NA_real_, mo_unemp_until_emp_sd = NA_real_, mo_unemp_until_emp_na = NA_real_)
  tempRU_descript <- data.frame(data = paste0("MI", mi),
                                strata = "RU", 
                                n_ind = NROW(unique(tempRU$unique_id)), n_obs = NROW(tempRU),
                                age_mu = mean(tempRU$age, na.rm = T), 
                                age_sd = sd(tempRU$age, na.rm=T), 
                                age_na = NROW(tempRU[is.na(tempRU$age),]),
                                gender_male = NROW(tempRU_inds[tempRU_inds$male == 1 & !is.na(tempRU_inds$male),]), 
                                gender_male_perc = mean(tempRU_inds$male, na.rm=T),
                                gender_male_na = NROW(tempRU_inds[is.na(tempRU_inds$male),]),
                                race_nhwhite_n = NROW(tempRU_inds[tempRU_inds$black == 0 & tempRU_inds$other == 0 & tempRU_inds$hispanic == 0, ]), 
                                race_nhwhite_perc = (1-mean(tempRU_inds$black, na.rm=T) - mean(tempRU_inds$other, na.rm=T) - mean(tempRU_inds$hispanic, na.rm=T))*100,
                                race_nhblack_n = NROW(tempRU_inds[tempRU_inds$black == 1,]),
                                race_nhblack_perc = mean(tempRU_inds$black, na.rm=T),
                                race_nhother_n = NROW(tempRU_inds[tempRU_inds$other == 1,]),
                                race_nhother_perc = mean(tempRU_inds$other, na.rm=T),
                                race_hispanic_n = NROW(tempRU_inds[tempRU_inds$hispanic == 1,]),
                                race_hispanic_perc = mean(tempRU_inds$hispanic, na.rm=T),
                                base_eq_mu = mean(tempRU$Y_base, na.rm = T), 
                                base_eq_sd = sd(tempRU$Y_base, na.rm=T), 
                                base_eq_na = NROW(tempRU[is.na(tempRU$Y_base),]),
                                eq_change_mu = mean(tempRU$change_eq, na.rm = T), 
                                eq_change_sd = sd(tempRU$change_eq, na.rm=T), 
                                eq_change_na = NROW(tempRU[is.na(tempRU$change_eq),]),
                                moderator_n = NROW(tempRU[tempRU$moderator_receipt == 1 & !is.na(tempRU$moderator_receipt),]), 
                                moderator_perc = mean(tempRU$moderator_receipt, na.rm = T), 
                                moderator_na = NROW(tempRU[is.na(tempRU$moderator_receipt),]),
                                mo_unemp_until_emp_mu = mean(tempRU$mo_unemp_until_emp, na.rm = T), 
                                mo_unemp_until_emp_sd = sd(tempRU$mo_unemp_until_emp, na.rm=T), 
                                mo_unemp_until_emp_na = NROW(tempRU[is.na(tempRU$mo_unemp_until_emp),]))
  if (mi == 1) {
    mi_descriptives <- rbind(tempCE_descript, tempRU_descript)
  } else {
    temp_descriptives <- rbind(tempCE_descript, tempRU_descript)
    mi_descriptives <- rbind(mi_descriptives, temp_descriptives)
  }
  rm(tempCE, tempCE_inds, tempCE_descript, tempRU, tempRU_inds, tempRU_descript, temp, mi, temp_descriptives)
}

mi_descriptives[,c(2:30)] %>%
  group_by(strata) %>%
  mutate(n_ind = mean(n_ind, na.rm=T), n_obs = mean(n_obs, na.rm=T),
         age_mu = mean(age_mu, na.rm=T), age_sd = mean(age_sd, na.rm=T), age_na = mean(age_na, na.rm=T),
         gender_male = mean(gender_male, na.rm=T), gender_male_perc = mean(gender_male_perc, na.rm=T), gender_male_na = mean(gender_male_na, na.rm=T),
         race_nhwhite_n = mean(race_nhwhite_n, na.rm=T), race_nhwhite_perc = mean(race_nhwhite_perc, na.rm=T), 
         race_nhblack_n = mean(race_nhblack_n, na.rm=T), race_nhblack_perc = mean(race_nhblack_perc, na.rm=T), 
         race_nhother_n = mean(race_nhother_n, na.rm=T), race_nhother_perc = mean(race_nhother_perc, na.rm=T), 
         race_hispanic_n = mean(race_hispanic_n, na.rm=T), race_hispanic_perc = mean(race_hispanic_perc, na.rm=T), 
         base_eq_mu = mean(base_eq_mu, na.rm=T), base_eq_sd = mean(base_eq_sd, na.rm=T), base_eq_na = mean(base_eq_na, na.rm=T),
         eq_change_mu = mean(eq_change_mu, na.rm=T), eq_change_sd = mean(eq_change_sd, na.rm=T), eq_change_na = mean(eq_change_na, na.rm=T),
         mo_unemp_until_emp_mu = mean(mo_unemp_until_emp_mu, na.rm=T), mo_unemp_until_emp_sd = mean(mo_unemp_until_emp_sd, na.rm=T), mo_unemp_until_emp_na = mean(mo_unemp_until_emp_na, na.rm=T),
         moderator_n = mean(moderator_n, na.rm=T), moderator_perc = mean(moderator_perc, na.rm=T), moderator_na = mean(moderator_na, na.rm=T),
         mo_unemp_until_emp_mu = mean(mo_unemp_until_emp_mu, na.rm=T), mo_unemp_until_emp_sd = mean(mo_unemp_until_emp_sd, na.rm=T), mo_unemp_until_emp_na = mean(mo_unemp_until_emp_na, na.rm=T)) %>%
  distinct() -> mi_descriptives

mi_descriptives$data <- "MI"
mi_descriptives <- mi_descriptives[,names(cc_descriptives)]

descriptives <- rbind(cc_descriptives, mi_descriptives)
rm(cc_descriptives, mi_descriptives)

#Comparing lengths unemployed until re-employment between those receiving and not receiving UI
by(mi1_data$mo_unemp_until_emp[mi1_data$prepost == "RU" & !is.na(mi1_data$prepost)], 
   INDICES = mi1_data$moderator_receipt[mi1_data$prepost == "RU" & !is.na(mi1_data$prepost)], FUN = summary)
by(mi1_data$mo_unemp_until_emp[mi1_data$prepost == "RU" & !is.na(mi1_data$prepost)], 
   INDICES = mi1_data$moderator_receipt[mi1_data$prepost == "RU" & !is.na(mi1_data$prepost)], FUN = function(x) {sd(x, na.rm=T)})

#Comparing occupational sectors pre- and post among MI dataset 1 
mi1_data %>% 
  arrange(unique_id, year) %>%
  group_by(unique_id) %>%
  mutate(next_occ_farmforfish = case_when(is.null(lead(base_occ_farmforfish,1)) | is.na(lead(base_occ_farmforfish,1)) ~ NA_real_,
                                          !is.na(lead(base_occ_farmforfish,1)) ~ lead(base_occ_farmforfish,1),
                                          TRUE ~ NA_real_),
         next_occ_managerial = case_when(is.null(lead(base_occ_managerial,1)) | is.na(lead(base_occ_managerial,1)) ~ NA_real_,
                                         !is.na(lead(base_occ_managerial,1)) ~ lead(base_occ_managerial,1),
                                         TRUE ~ NA_real_),
         next_occ_military = case_when(is.null(lead(base_occ_military,1)) | is.na(lead(base_occ_military,1)) ~ NA_real_,
                                       !is.na(lead(base_occ_military,1)) ~ lead(base_occ_military,1),
                                       TRUE ~ NA_real_),
         next_occ_opfablab = case_when(is.null(lead(base_occ_opfablab,1)) | is.na(lead(base_occ_opfablab,1)) ~ NA_real_,
                                       !is.na(lead(base_occ_opfablab,1)) ~ lead(base_occ_opfablab,1),
                                       TRUE ~ NA_real_),
         next_occ_precision = case_when(is.null(lead(base_occ_precision,1)) | is.na(lead(base_occ_precision,1)) ~ NA_real_,
                                        !is.na(lead(base_occ_precision,1)) ~ lead(base_occ_precision,1),
                                        TRUE ~ NA_real_),
         next_occ_profspec = case_when(is.null(lead(base_occ_profspec,1)) | is.na(lead(base_occ_profspec,1)) ~ NA_real_,
                                       !is.na(lead(base_occ_profspec,1)) ~ lead(base_occ_profspec,1),
                                       TRUE ~ NA_real_),
         next_occ_services = case_when(is.null(lead(base_occ_services,1)) | is.na(lead(base_occ_services,1)) ~ NA_real_,
                                       !is.na(lead(base_occ_services,1)) ~ lead(base_occ_services,1),
                                       TRUE ~ NA_real_)) %>% 
  ungroup() -> mi1_data

summary(mi1_data[mi1_data$moderator_receipt==0 & mi1_data$prepost == "CE" & !is.na(mi1_data$prepost) & !is.na(mi1_data$moderator_receipt), 
                 c("base_occ_farmforfish", "base_occ_managerial", "base_occ_military", "base_occ_opfablab", "base_occ_precision", "base_occ_profspec", "base_occ_services")])
summary(mi1_data[mi1_data$moderator_receipt==0 & mi1_data$prepost == "CE" & !is.na(mi1_data$prepost) & !is.na(mi1_data$moderator_receipt), 
                 c("next_occ_farmforfish", "next_occ_managerial", "next_occ_military", "next_occ_opfablab", "next_occ_precision", "next_occ_profspec", "next_occ_services")])
summary(mi1_data[mi1_data$moderator_receipt==0 & mi1_data$prepost == "RU" & !is.na(mi1_data$prepost) & !is.na(mi1_data$moderator_receipt), 
                 c("base_occ_farmforfish", "base_occ_managerial", "base_occ_military", "base_occ_opfablab", "base_occ_precision", "base_occ_profspec", "base_occ_services")])
summary(mi1_data[mi1_data$moderator_receipt==0 & mi1_data$prepost == "RU" & !is.na(mi1_data$prepost) & !is.na(mi1_data$moderator_receipt), 
                 c("next_occ_farmforfish", "next_occ_managerial", "next_occ_military", "next_occ_opfablab", "next_occ_precision", "next_occ_profspec", "next_occ_services")])
summary(mi1_data[mi1_data$moderator_receipt==1 & mi1_data$prepost == "RU" & !is.na(mi1_data$prepost) & !is.na(mi1_data$moderator_receipt), 
                 c("base_occ_farmforfish", "base_occ_managerial", "base_occ_military", "base_occ_opfablab", "base_occ_precision", "base_occ_profspec", "base_occ_services")])
summary(mi1_data[mi1_data$moderator_receipt==1 & mi1_data$prepost == "RU" & !is.na(mi1_data$prepost) & !is.na(mi1_data$moderator_receipt), 
                 c("next_occ_farmforfish", "next_occ_managerial", "next_occ_military", "next_occ_opfablab", "next_occ_precision", "next_occ_profspec", "next_occ_services")])

mi1_data %>% 
  mutate(base_occupations = case_when(base_occ_farmforfish == 1 ~ "farmforfish",
                                      base_occ_managerial == 1 ~ "managerial",
                                      base_occ_military == 1 ~ "military",
                                      base_occ_opfablab == 1 ~ "opfablab",
                                      base_occ_precision == 1 ~ "precision",
                                      base_occ_profspec == 1 ~ "profspec",
                                      base_occ_services == 1 ~ "services",
                                      is.na(base_occ_managerial) ~ NA_character_,
                                      TRUE ~ "techadmin"),
         next_occupations = case_when(next_occ_farmforfish == 1 ~ "farmforfish",
                                      next_occ_managerial == 1 ~ "managerial",
                                      next_occ_military == 1 ~ "military",
                                      next_occ_opfablab == 1 ~ "opfablab",
                                      next_occ_precision == 1 ~ "precision",
                                      next_occ_profspec == 1 ~ "profspec",
                                      next_occ_services == 1 ~ "services",
                                      is.na(next_occ_managerial) ~ NA_character_,
                                      TRUE ~ "techadmin")) %>%
  group_by(base_occupations) %>%
  mutate(mu_occupation_eq = mean(Y_base, na.rm=T)) %>%
  ungroup() -> mi1_data

mi1_data %>% 
  select(base_occupations, mu_occupation_eq) %>% 
  distinct() %>%
  filter(!is.na(base_occupations)) -> mi1_occupation_data
mi1_occupation_data_cc <- mi1_occupation_data; mi1_occupation_data_cc$group <- "CE"
mi1_occupation_data_ru_noui <- mi1_occupation_data; mi1_occupation_data_ru_noui$group <- "RU (No UI)"
mi1_occupation_data_ru_ui <- mi1_occupation_data; mi1_occupation_data_ru_ui$group <- "RU (UI)"
mi1_occupation_data <- rbind(mi1_occupation_data_cc, mi1_occupation_data_ru_noui, mi1_occupation_data_ru_ui)

mi1_occupation_data$prop_base_occ <- c(0.3044, 0.1179, 0.10067, 0.14677, 0.08695, 0.00975, 0.22893, 0.00463, 
                                       0.28435, 0.13972, 0.09931, 0.18476, 0.0892, 0.01415, 0.17985, 0.0087,
                                       0.3278, 0.1711, 0.0999, 0.1237, 0.1010, 0.0247, 0.1526, 0)
mi1_occupation_data$prop_next_occ <- c(0.30783, 0.11475, 0.10032, 0.14473, 0.08632, 0.01033, 0.23009, 0.00564,
                                       0.28835, 0.137186, 0.10988, 0.18641, 0.08574, 0.01207, 0.17371, 0.0067,
                                       0.3499, 0.1648, 0.1174, 0.1309, 0.0767, 0.0181, 0.140, 0.0022)

mi1_occupation_data %>% mutate(prepost_occ_prop_shift = (prop_next_occ - prop_base_occ)*100,
                               base_occ_eq_component = mu_occupation_eq*prop_base_occ,
                               next_occ_eq_component = mu_occupation_eq*prop_next_occ,
                               prepost_occ_eq_component_shift = next_occ_eq_component - base_occ_eq_component) -> mi1_occupation_data 
rm(mi1_occupation_data_cc, mi1_occupation_data_ru_noui, mi1_occupation_data_ru_ui)

ggplot(mi1_occupation_data, aes(x = prepost_occ_prop_shift, y = base_occupations, col = group)) +
  geom_point(size = 3) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  labs(x = "Pre-Post Sector Percentage Change (%)", y = "Occupational Sector", col = "Pre-Post Group") +
  scale_color_manual(values = c("grey", "purple", "purple4")) +
  theme_classic() + theme(legend.position = "bottom") + xlim(-2.5,2.5)

by(mi1_occupation_data$prepost_occ_eq_component_shift, INDICES = mi1_occupation_data$group, FUN = sum)

#### Step 6 - Creating summary figures ####
fig_data <- final_results[,c(1,2,9,11,12)]

#PSID Conference Poster Figure - MI Findings
poster_fig_data <- fig_data[fig_data$data == "MI" & fig_data$estimand %in% c("ATT", "CATT_UI_NO", "CATT_UI_YES", "CATT_UI_DUR_LESSMED", "CATT_UI_DUR_MOREMED", "CATT_UI_GEN_LESSMED", "CATT_UI_GEN_MOREMED"), ]
poster_fig_data$estimate <- c("Overall", "No UI", "UI (Any)", 
                              "UI (< Median Duration)", "UI (>= Median Duration)", 
                              "UI (< Median Generosity)", "UI (>= Median Generosity)")
poster_fig_data$estimate <- factor(poster_fig_data$estimate, levels = c("Overall", "No UI", "UI (Any)", 
                                                                        "UI (< Median Duration)", "UI (>= Median Duration)", 
                                                                        "UI (< Median Generosity)", "UI (>= Median Generosity)"))

ggplot(poster_fig_data, aes(x = est, y = fct_rev(estimate), col = estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey", linewidth = 1) +
  geom_point(size = 2) + geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 1) +
  labs(x = "ATT of Unemployment on EQ (95% CI)", y = "") +
  xlim(-2, 0.05) +
  theme_classic() +
  scale_color_manual(values = c("purple3", "purple3", "purple3", "purple3", "purple3", "purple3", "purple3")) +
  theme(legend.position = "none")

#PSID CC Finding Slide - Scarring Effects
cc_fig_scarring_data <- fig_data[fig_data$data == "CC" & fig_data$estimand %in% unique(fig_data$estimand)[c(1:13)], ]
cc_fig_scarring_data$strata <- c("Overall", 
                                 "EQ: Q1", "EQ: Q2", "EQ: Q3", "EQ: Q4",
                                 "Gender: Female", "Gender: Male", 
                                 "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                                 "Length Unemp. < Median", "Length Unemp. >= Median")
cc_fig_scarring_data$strata <- factor(cc_fig_scarring_data$strata, levels = c("Overall", 
                                                                              "EQ: Q1", "EQ: Q2", "EQ: Q3", "EQ: Q4",
                                                                              "Gender: Female", "Gender: Male", 
                                                                              "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                                                                              "Length Unemp. < Median", "Length Unemp. >= Median"))

ggplot(cc_fig_scarring_data, aes(x = est, y = fct_rev(strata))) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey", linewidth = 1) +
  geom_point(size = 2) + geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 1) +
  labs(x = "ATT of Unemployment on EQ (95% CI)", y = "") +
  theme_classic()

#PSID MI Finding Slide - Scarring Effects
mi_fig_scarring_data <- fig_data[fig_data$data == "MI" & fig_data$estimand %in% unique(fig_data$estimand)[c(1:13)], ]
mi_fig_scarring_data$strata <- c("Overall", 
                                 "EQ: Q1", "EQ: Q2", "EQ: Q3", "EQ: Q4",
                                 "Gender: Female", "Gender: Male", 
                                 "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                                 "Length Unemp. < Median", "Length Unemp. >= Median")
mi_fig_scarring_data$strata <- factor(mi_fig_scarring_data$strata, levels = c("Overall", 
                                                                              "EQ: Q1", "EQ: Q2", "EQ: Q3", "EQ: Q4",
                                                                              "Gender: Female", "Gender: Male", 
                                                                              "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                                                                              "Length Unemp. < Median", "Length Unemp. >= Median"))

ggplot(mi_fig_scarring_data, aes(x = est, y = fct_rev(strata))) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey", linewidth = 1) +
  geom_point(size = 2) + geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 1) +
  labs(x = "ATT of Unemployment on EQ (95% CI)", y = "") +
  theme_classic()

#PSID CC Finding Slide - UI Receipt
cc_fig_ui_data <- fig_data[fig_data$data == "CC" & fig_data$estimand %in% unique(fig_data$estimand)[c(14,15,22:41)], ]
cc_fig_ui_data$ui <- rep(c("No", "Yes"), times = 11)
cc_fig_ui_data$strata <- rep(c("Overall", 
                               "EQ: Q1", "EQ: Q4",
                               "Gender: Female", "Gender: Male", 
                               "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                               "Length Unemp. < Median", "Length Unemp. >= Median"), each = 2)
cc_fig_ui_data$strata <- factor(cc_fig_ui_data$strata, levels = c("Overall", 
                                                                  "EQ: Q1", "EQ: Q4",
                                                                  "Gender: Female", "Gender: Male", 
                                                                  "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                                                                  "Length Unemp. < Median", "Length Unemp. >= Median"))

ggplot(cc_fig_ui_data, aes(x = est, y = fct_rev(strata), col = ui)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey", linewidth = 1) +
  geom_point(size = 2, position = position_dodge(width=0.5)) + geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 1, position = position_dodge(width=0.5)) +
  labs(x = "ATT of Unemployment on EQ (95% CI)", y = "", col = "UI Receipt") +
  theme_classic() +
  scale_color_manual(values = rep(c("grey", "purple4"), times = 11))

#PSID MI Finding Slide - UI Receipt
mi_fig_ui_data <- fig_data[fig_data$data == "MI" & fig_data$estimand %in% unique(fig_data$estimand)[c(14,15,22:41)], ]
mi_fig_ui_data$ui <- rep(c("No", "Yes"), times = 11)
mi_fig_ui_data$strata <- rep(c("Overall", 
                               "EQ: Q1", "EQ: Q4",
                               "Gender: Female", "Gender: Male", 
                               "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                               "Length Unemp. < Median", "Length Unemp. >= Median"), each = 2)
mi_fig_ui_data$strata <- factor(mi_fig_ui_data$strata, levels = c("Overall", 
                                                                  "EQ: Q1", "EQ: Q4",
                                                                  "Gender: Female", "Gender: Male", 
                                                                  "Race & Eth: NH White", "Race & Eth: NH Black", "Race & Eth: NH Other", "Race & Eth: Hispanic",
                                                                  "Length Unemp. < Median", "Length Unemp. >= Median"))

ggplot(mi_fig_ui_data, aes(x = est, y = fct_rev(strata), col = ui)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey", linewidth = 1) +
  geom_point(size = 2, position = position_dodge(width=0.5)) + geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 1, position = position_dodge(width=0.5)) +
  labs(x = "ATT of Unemployment on EQ (95% CI)", y = "", col = "UI Receipt") +
  theme_classic() +
  scale_color_manual(values = rep(c("grey", "purple4"), times = 11))
