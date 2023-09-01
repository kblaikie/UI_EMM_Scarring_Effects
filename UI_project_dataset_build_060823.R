#### Unemployment-UI-EQ Project Dataset Build Script 3 ####
#Author - Kieran Blaikie
#Date - 08 June 2023

#Overview - This script aims to take monthly data on employment and UI receipt and:
#             1) Identify unemployment status in the 12 month period post t
#             2) Identify UI receipt in the 6 month period post-unemployment
#             3) Restrict to putatively eligible participants (aged 18-64, known in labor force/unemployed)
#             4) Code covariates as they'd be included in final analyses
#             4) Construct N multiply imputed datasets, based on % incomplete cases
#             5) For each N MI dataset, create final exposure and outcome variables

#NOTE - Data is collected in such a way that: - unemployment data (A) is from t-1 (-2 years) to t-0.5 (-1 year) in record t - so we use lead(A,1) to get A from t to t+0.5
#                                             - UI receipt data (M) is from t-0.5 (-1 year) to t-0.25 (-6 months) in record t - so we use lead(M,1) to get M from t+0.5 to t+0.75
#                                             - EQ data (Y) relates to t (current) in record t - so we use this as 'baseline' Y and use lead(Y,1+) for post-exposure Y 
#     - We aim to assume a lag L = 2 for A & Y, and = 1 for all covariates (Z) and M - i.e. we assume A and Y are independent of A-3, Y-3, etc 
#     - Therefore, per row we want variables for (A-2, A-1, A), (Y-2, Y-1, Y), (M-1, M), (Z-1, Z), C
#     - To have Z-1, Y-2, Y-1 variables in t=2001 (the earliest record for UI data), we should compute these in the demographic dataset before merging with the UI data

#Changes - Compared to UI_project_dataset_build_032723.R, this script:
#     - Incorporates multiple imputation & creating analysis variables
#     - No longer creates as many 'eligibility' indicators
#     - Fixed error in Step 3 where the wrong variables were summed to determine 12-month unemployment

#Script Structure - This script has the following steps:
#     - 0) Merging datasets and subsetting to needed variables
#     - 1) Creating variables for month-specific unemployment status in interview year t (T0) and t+1 (T1)
#     - 2) Creating indicator for any unknown monthly unemployment status information in 12 months post-interview (incl. interview month)
#     - 3) Creating indicator for any unemployment in 12 months post-interview
#     - 4-6) Creating variables for: - (continuous) first-month unemployed post-interview if relevant
#                                  - (binary) whether already unemployed pre-interview
#                                  - (continuous) known (and minimum) number of consecutive months unemployed pre-interview
#                                  - (continuous) cumulative continuous months of unemployment pre-re-employment (if employed at prior interview)
#                                  - (categorical) indicator for unknown or censored months of unemployment pre-re-employment
#     - 7) Creating variables for month-specific UI receipt status in post-interview year t+1 (T1)
#     - 8,11) Creating indicators for: - any unknown monthly UI receipt information in 6 months post-unemployment (not incl. month becoming unemployed)
#                                      - any UI receipt in 6 months post-unemployment
#                                      - number of months of UI receipt in 6 months post-unemployment
#     - 9) Conducting simple imputation of UI receipt during the 12-months post-interview (or unemployment, whichever later)
#     - 10) Merging in State max UI duration and benefit information for 6 months pre-unemployment (or interview-month if not unemployed)
#     - 12) Merging in State unemployment rate and GSP per Capita for 3 months pre-unemployment (or interview month, if not unemployed)
#     - 13) Restricting to the necessary set of variables for multiple imputation, determining wave eligibility, and modelling
#     - 14) Creating multiply imputed datasets
#     - XX) Restricting to waves where respondents were not self-employed pre-MI (note EQ would be missing for all those self-employed as a pre-PCA step in 'dataset_build_092921.R')
#     - XX) Constructing N multiply imputed datasets, where N = the % incomplete cases across model variables
#     - XX) For each MI dataset, we create: - pre-exposure and post-exposure EQ variables over different year lags/leads 
#                                           - 'previously most recent', 'most recent', and 'post-exposure' outcome variables

#Loading needed libraries
library(tidyverse) #Data-wrangling & manipulation
library(readstata13) # For reading in EQ dataset made in Part 1
library(qs) #Saving/loading
library(zoo) #For carry-forward/backward imputation
library(mice) #Multiple imputation
library(miceadds) #Longitudinal multiple imputation

#Setting directory
directory <- "R:/Project/precarityR01/PSID/kieran_ui_project/Data/"

#### Step 0 - Merging datasets ####
#Loading required data
demog_data <- qread(paste0(directory, "psid_dataset_110321.qs")) #Demographic data
ui_ind_data <- read.csv(paste0(directory, "psid_dataset_100422.csv")) #Individual UI data
ui_state_data <- read.csv(paste0(directory, "sarah_state_ui_dataset_050123.csv")) #State UI data
ui_state_elig_data <- read.csv(paste0(directory, "ui_eligibility_formatted_030823.csv")) #State UI Eligibility data
state_ur_data <- read.csv(paste0(directory, "BLS_LAUS_9419.csv")) #State Unemployment Rate data
state_gsp_data <- read.csv(paste0(directory, "BEA_GSP_9719.csv")) #State GSP per Capita data
eq_data <- read.dta13(paste0(directory, "eq_040223.dta")) #Multidim. EQ data. Note this dataset was made running dataset_build_092921, restricting age to 18:64 instead of 30:60
cpi_data <- read.csv(paste0(directory, "cpi.csv")) #Fed Reserve CPI

#Selecting required demographic variables 
demog_data <- demog_data %>% 
  select(
    #Identifiers
    unique_id, year, 
    #Fixed
    gender, race, ethnicity, nativity, parents_poor, 
    #Varying
    state, region, marital_status, age, education, disabl_limits_work, SRH, k6, 
    employment_status, self_employed, occupation_broad, emp_tenure_months,
    disabl_limits_work_partner, employment_status_partner, family_income,
    family_wealth_no_home_equity)

#Selecting required EQ variables (keeping dimension vars for exposure group comparisons, linear EQ var)
eq_data <- eq_data[,c(1:12,16)]

#Selecting required UI individual-level variables
ui_ind_data <- ui_ind_data[,c(2,3,7:22,36:59)]

#Selecting required UI state-level variables
ui_state_data <- ui_state_data[,c(1,2,4,6,7)]
state_gsp_data <- state_gsp_data[,c(2:5,8)]
state_ur_data <- state_ur_data[,c(2,4,8,9)]

#Merging together
#NOTE - Can't merge in state-level variables until I've created a variable for
#       month and year unemployment spell began in year post-interview 
data <- merge(demog_data, eq_data, by = c("unique_id", "year", "age"), all = T)
data <- merge(data, cpi_data, by = c("year"), all = T)
data <- merge(data, ui_ind_data, by = c("unique_id", "year"), all = T)
data <- data[!is.na(data$unique_id), ] #Removing empty rows
rm(demog_data, eq_data, cpi_data, ui_ind_data)

#Constructing Unemployment and UI variables post-current interview
#NOTE - To work out unemployment and UI vars post-current-interview, we use info
#       from subsequent wave assuming:
#         1) it is 2 years after, and 
#         2) current interviews occurred in typical interview years
#           > This is needed to actually have UI data post-interview

#### Step 1 - Creating 'unemp_status_tX_month_lead' vars for unemployment status per month in years t, t+1 ####
#NOTE - This requires 1) using lead() on emp_status_t-2/-1 variables to get emp_status_t,t+1,
#                     2) flip 0 and 1 coding emp_status_t,t+1 vars to reflect unemp_status_t,t+1
data <- data %>% arrange(unique_id, year)
library(data.table)
data <- as.data.table(data)
data[, ':=' (unemp_status_t0_jan_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_jan)), NA_integer_, fifelse(lead(emp_status_t2_jan) == 1, 0, 1)),
             unemp_status_t0_feb_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_feb)), NA_integer_, fifelse(lead(emp_status_t2_feb) == 1, 0, 1)),
             unemp_status_t0_mar_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_mar)), NA_integer_, fifelse(lead(emp_status_t2_mar) == 1, 0, 1)),
             unemp_status_t0_apr_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_apr)), NA_integer_, fifelse(lead(emp_status_t2_apr) == 1, 0, 1)),
             unemp_status_t0_may_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_may)), NA_integer_, fifelse(lead(emp_status_t2_may) == 1, 0, 1)),
             unemp_status_t0_jun_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_jun)), NA_integer_, fifelse(lead(emp_status_t2_jun) == 1, 0, 1)),
             unemp_status_t0_jul_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_jul)), NA_integer_, fifelse(lead(emp_status_t2_jul) == 1, 0, 1)),
             unemp_status_t0_aug_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_aug)), NA_integer_, fifelse(lead(emp_status_t2_aug) == 1, 0, 1)),
             unemp_status_t0_sep_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_sep)), NA_integer_, fifelse(lead(emp_status_t2_sep) == 1, 0, 1)),
             unemp_status_t0_oct_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_oct)), NA_integer_, fifelse(lead(emp_status_t2_oct) == 1, 0, 1)),
             unemp_status_t0_nov_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_nov)), NA_integer_, fifelse(lead(emp_status_t2_nov) == 1, 0, 1)),
             unemp_status_t0_dec_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t2_dec)), NA_integer_, fifelse(lead(emp_status_t2_dec) == 1, 0, 1)),
             unemp_status_t1_jan_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_jan)), NA_integer_, fifelse(lead(emp_status_t1_jan) == 1, 0, 1)),
             unemp_status_t1_feb_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_feb)), NA_integer_, fifelse(lead(emp_status_t1_feb) == 1, 0, 1)),
             unemp_status_t1_mar_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_mar)), NA_integer_, fifelse(lead(emp_status_t1_mar) == 1, 0, 1)),
             unemp_status_t1_apr_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_apr)), NA_integer_, fifelse(lead(emp_status_t1_apr) == 1, 0, 1)),
             unemp_status_t1_may_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_may)), NA_integer_, fifelse(lead(emp_status_t1_may) == 1, 0, 1)),
             unemp_status_t1_jun_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_jun)), NA_integer_, fifelse(lead(emp_status_t1_jun) == 1, 0, 1)),
             unemp_status_t1_jul_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_jul)), NA_integer_, fifelse(lead(emp_status_t1_jul) == 1, 0, 1)),
             unemp_status_t1_aug_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_aug)), NA_integer_, fifelse(lead(emp_status_t1_aug) == 1, 0, 1)),
             unemp_status_t1_sep_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_sep)), NA_integer_, fifelse(lead(emp_status_t1_sep) == 1, 0, 1)),
             unemp_status_t1_oct_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_oct)), NA_integer_, fifelse(lead(emp_status_t1_oct) == 1, 0, 1)),
             unemp_status_t1_nov_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_nov)), NA_integer_, fifelse(lead(emp_status_t1_nov) == 1, 0, 1)),
             unemp_status_t1_dec_lead = fifelse(lead(year) != (year+2) | is.na(lead(emp_status_t1_dec)), NA_integer_, fifelse(lead(emp_status_t1_dec) == 1, 0, 1))), by = .(unique_id)]
data <- as.data.frame(data)
detach("package:data.table", unload = TRUE)

#### Step 2 - Creating indicator for missing relevant unemployment indicators in 12 months post-interview (incl. interview month) ####
#NOTE - 19.8% of observations from 2001 onwards lack all relevant monthly employment information
data$unemp_status_period_na <- 0
data$unemp_status_period_na[data$interview_month == 1 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_jan_lead) | is.na(data$unemp_status_t0_feb_lead) | is.na(data$unemp_status_t0_mar_lead) | 
                                 is.na(data$unemp_status_t0_apr_lead) | is.na(data$unemp_status_t0_may_lead) | is.na(data$unemp_status_t0_jun_lead) | 
                                 is.na(data$unemp_status_t0_jul_lead) | is.na(data$unemp_status_t0_aug_lead) | is.na(data$unemp_status_t0_sep_lead) | 
                                 is.na(data$unemp_status_t0_oct_lead) | is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 2 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_feb_lead) | is.na(data$unemp_status_t0_mar_lead) | is.na(data$unemp_status_t0_apr_lead) | 
                                 is.na(data$unemp_status_t0_may_lead) | is.na(data$unemp_status_t0_jun_lead) | is.na(data$unemp_status_t0_jul_lead) | 
                                 is.na(data$unemp_status_t0_aug_lead) | is.na(data$unemp_status_t0_sep_lead) | is.na(data$unemp_status_t0_oct_lead) | 
                                 is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 3 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t1_mar_lead) | is.na(data$unemp_status_t0_apr_lead) | is.na(data$unemp_status_t0_may_lead) | 
                                 is.na(data$unemp_status_t0_jun_lead) | is.na(data$unemp_status_t0_jul_lead) | is.na(data$unemp_status_t0_aug_lead) | 
                                 is.na(data$unemp_status_t0_sep_lead) | is.na(data$unemp_status_t0_oct_lead) | is.na(data$unemp_status_t0_nov_lead) | 
                                 is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead) | is.na(data$unemp_status_t1_feb_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 4 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_apr_lead) | is.na(data$unemp_status_t0_may_lead) | is.na(data$unemp_status_t0_jun_lead) | 
                                 is.na(data$unemp_status_t0_jul_lead) | is.na(data$unemp_status_t0_aug_lead) | is.na(data$unemp_status_t0_sep_lead) | 
                                 is.na(data$unemp_status_t0_oct_lead) | is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead) | 
                                 is.na(data$unemp_status_t1_jan_lead) | is.na(data$unemp_status_t1_feb_lead) | is.na(data$unemp_status_t1_mar_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 5 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_may_lead) | is.na(data$unemp_status_t0_jun_lead) | is.na(data$unemp_status_t0_jul_lead) | 
                                 is.na(data$unemp_status_t0_aug_lead) | is.na(data$unemp_status_t0_sep_lead) | is.na(data$unemp_status_t0_oct_lead) | 
                                 is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead) | 
                                 is.na(data$unemp_status_t1_feb_lead) | is.na(data$unemp_status_t1_mar_lead) | is.na(data$unemp_status_t1_apr_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 6 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_jun_lead) | is.na(data$unemp_status_t0_jul_lead) | is.na(data$unemp_status_t0_aug_lead) | 
                                 is.na(data$unemp_status_t0_sep_lead) | is.na(data$unemp_status_t0_oct_lead) | is.na(data$unemp_status_t0_nov_lead) | 
                                 is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead) | is.na(data$unemp_status_t1_feb_lead) | 
                                 is.na(data$unemp_status_t1_mar_lead) | is.na(data$unemp_status_t1_apr_lead) | is.na(data$unemp_status_t1_may_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 7 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_jul_lead) | is.na(data$unemp_status_t0_aug_lead) | is.na(data$unemp_status_t0_sep_lead) | 
                                 is.na(data$unemp_status_t0_oct_lead) | is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead) |
                                 is.na(data$unemp_status_t1_jan_lead) | is.na(data$unemp_status_t1_feb_lead) | is.na(data$unemp_status_t1_mar_lead) | 
                                 is.na(data$unemp_status_t1_apr_lead) | is.na(data$unemp_status_t1_may_lead) | is.na(data$unemp_status_t1_jun_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 8 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_aug_lead) | is.na(data$unemp_status_t0_sep_lead) | is.na(data$unemp_status_t0_oct_lead) | 
                                 is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead) | 
                                 is.na(data$unemp_status_t1_feb_lead) | is.na(data$unemp_status_t1_mar_lead) | is.na(data$unemp_status_t1_apr_lead) | 
                                 is.na(data$unemp_status_t1_may_lead) | is.na(data$unemp_status_t1_jun_lead) | is.na(data$unemp_status_t1_jul_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 9 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_sep_lead) | is.na(data$unemp_status_t0_oct_lead) | is.na(data$unemp_status_t0_nov_lead) | 
                                 is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead) | is.na(data$unemp_status_t1_feb_lead) | 
                                 is.na(data$unemp_status_t1_mar_lead) | is.na(data$unemp_status_t1_apr_lead) | is.na(data$unemp_status_t1_may_lead) | 
                                 is.na(data$unemp_status_t1_jun_lead) | is.na(data$unemp_status_t1_jul_lead) | is.na(data$unemp_status_t1_aug_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 10 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_oct_lead) | is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead) | 
                                 is.na(data$unemp_status_t1_jan_lead) | is.na(data$unemp_status_t1_feb_lead) | is.na(data$unemp_status_t1_mar_lead) | 
                                 is.na(data$unemp_status_t1_apr_lead) | is.na(data$unemp_status_t1_may_lead) | is.na(data$unemp_status_t1_jun_lead) | 
                                 is.na(data$unemp_status_t1_jul_lead) | is.na(data$unemp_status_t1_aug_lead) | is.na(data$unemp_status_t1_sep_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 11 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_nov_lead) | is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead) | 
                                 is.na(data$unemp_status_t1_feb_lead) | is.na(data$unemp_status_t1_mar_lead) | is.na(data$unemp_status_t1_apr_lead) | 
                                 is.na(data$unemp_status_t1_may_lead) | is.na(data$unemp_status_t1_jun_lead) | is.na(data$unemp_status_t1_jul_lead) | 
                                 is.na(data$unemp_status_t1_aug_lead) | is.na(data$unemp_status_t1_sep_lead) | is.na(data$unemp_status_t1_oct_lead))] <- 1
data$unemp_status_period_na[data$interview_month == 12 & !is.na(data$interview_month) & 
                              (is.na(data$unemp_status_t0_dec_lead) | is.na(data$unemp_status_t1_jan_lead) | is.na(data$unemp_status_t1_feb_lead) | 
                                 is.na(data$unemp_status_t1_mar_lead) | is.na(data$unemp_status_t1_apr_lead) | is.na(data$unemp_status_t1_may_lead) | 
                                 is.na(data$unemp_status_t1_jun_lead) | is.na(data$unemp_status_t1_jul_lead) | is.na(data$unemp_status_t1_aug_lead) | 
                                 is.na(data$unemp_status_t1_sep_lead) | is.na(data$unemp_status_t1_oct_lead) | is.na(data$unemp_status_t1_nov_lead))] <- 1
data$unemp_status_period_na <- ifelse(data$year < 2001 | data$interview_year != data$year | is.na(data$interview_month), 1, data$unemp_status_period_na)

#### Step 3 - Creating 'unemp_status_lead' for unemployment status in year post-interview ####
#NOTE - This requires 1) checking if the sum of unemp_status_tX_month_lead is >0 and 
#                     2) that all relevant employment information is known
data_pre2001 <- data[data$year < 2001, ]
data_int_jan <- data[data$interview_month == 1 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_feb <- data[data$interview_month == 2 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_mar <- data[data$interview_month == 3 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_apr <- data[data$interview_month == 4 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_may <- data[data$interview_month == 5 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_jun <- data[data$interview_month == 6 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_jul <- data[data$interview_month == 7 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_aug <- data[data$interview_month == 8 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_sep <- data[data$interview_month == 9 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_oct <- data[data$interview_month == 10 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_nov <- data[data$interview_month == 11 & !is.na(data$interview_month) & data$year >= 2001, ]
data_int_dec <- data[data$interview_month == 12 & !is.na(data$interview_month) & data$year >= 2001, ]

data_pre2001$unemp_status_lead <- NA_integer_
data_int_jan$unemp_status_lead <- ifelse(rowSums(data_int_jan[,c(75:86)], na.rm=T) > 0, 1,
                                         ifelse(data_int_jan$unemp_status_period_na == 1, NA_integer_, 0))
data_int_feb$unemp_status_lead <- ifelse(rowSums(data_int_feb[,c(76:87)], na.rm=T) > 0, 1,
                                         ifelse(data_int_feb$unemp_status_period_na == 1, NA_integer_, 0))
data_int_mar$unemp_status_lead <- ifelse(rowSums(data_int_mar[,c(77:88)], na.rm=T) > 0, 1,
                                         ifelse(data_int_mar$unemp_status_period_na == 1, NA_integer_, 0))
data_int_apr$unemp_status_lead <- ifelse(rowSums(data_int_apr[,c(78:89)], na.rm=T) > 0, 1,
                                         ifelse(data_int_apr$unemp_status_period_na == 1, NA_integer_, 0))
data_int_may$unemp_status_lead <- ifelse(rowSums(data_int_may[,c(79:90)], na.rm=T) > 0, 1,
                                         ifelse(data_int_may$unemp_status_period_na == 1, NA_integer_, 0))
data_int_jun$unemp_status_lead <- ifelse(rowSums(data_int_jun[,c(80:91)], na.rm=T) > 0, 1,
                                         ifelse(data_int_jun$unemp_status_period_na == 1, NA_integer_, 0))
data_int_jul$unemp_status_lead <- ifelse(rowSums(data_int_jul[,c(81:92)], na.rm=T) > 0, 1,
                                         ifelse(data_int_jul$unemp_status_period_na == 1, NA_integer_, 0))
data_int_aug$unemp_status_lead <- ifelse(rowSums(data_int_aug[,c(82:93)], na.rm=T) > 0, 1,
                                         ifelse(data_int_aug$unemp_status_period_na == 1, NA_integer_, 0))
data_int_sep$unemp_status_lead <- ifelse(rowSums(data_int_sep[,c(83:94)], na.rm=T) > 0, 1,
                                         ifelse(data_int_sep$unemp_status_period_na == 1, NA_integer_, 0))
data_int_oct$unemp_status_lead <- ifelse(rowSums(data_int_oct[,c(84:95)], na.rm=T) > 0, 1,
                                         ifelse(data_int_oct$unemp_status_period_na == 1, NA_integer_, 0))
data_int_nov$unemp_status_lead <- ifelse(rowSums(data_int_nov[,c(85:96)], na.rm=T) > 0, 1,
                                         ifelse(data_int_nov$unemp_status_period_na == 1, NA_integer_, 0))
data_int_dec$unemp_status_lead <- ifelse(rowSums(data_int_dec[,c(86:97)], na.rm=T) > 0, 1,
                                         ifelse(data_int_dec$unemp_status_period_na == 1, NA_integer_, 0))

#### Step 4 - Creating 1) 'first_mo_unemp_post_int' for first month of unemployment post-interview ####
#                      2) 'already_unemp_pre_int' for whether already unemployed in month prior to interview
#                      3) 'mo_unemp_pre_int' for known months unemployed pre-interview if unemployed prior to interview
#                      3) 'mo_unemp_pre_int_min' for minimum months unemployed pre-interview if unemployed prior to interview, based on known monthly unemployment data 
#                      4) 'mo_unemp_pre_int_na' indicating whether all relevant unemployment information was known
#NOTE - For 'first_mo_unemp_post_int', assign each person the first month where they reported being unemployed between interview (t0m1) and t0m12
#         - If missing monthly unemployment information, only assign first month X if t0m to tXm are not NA 
#         - If someone is not unemployed in the 12 months post-interview or is missing all monthly unemployment status information, reassign to NA
#     - For 'already_unemp_pre_int', assign 1 if unemployed in month prior to interview, 0 if not, and NA if missing this information
#     - For 'mo_unemp_pre_int', assign known consecutive months unemployed pre-interview (not incl. interview month) for respondents 'first unemployed' in the month post-interview, or NA if unknown 
#     - For 'mo_unemp_pre_int_min', assign minimum known consecutive months unemployed pre-interview
#     - For 'mo_unemp_pre_int_na', assign '0' if fully known months unemployed, or 1 if censored due to missing data
data_pre2001 %>% mutate(first_mo_unemp_post_int = NA_integer_, 
                        already_unemp_pre_int = NA_integer_,
                        mo_unemp_pre_int = NA_integer_,
                        mo_unemp_pre_int_min = NA_integer_) -> data_pre2001
data_int_jan %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1) ~ NA_integer_,
                                                            unemp_status_t0_jan_lead == 1 ~ 1,
                                                            unemp_status_t0_feb_lead == 1 ~ 2, unemp_status_t0_mar_lead == 1 ~ 3,
                                                            unemp_status_t0_apr_lead == 1 ~ 4, unemp_status_t0_may_lead == 1 ~ 5,
                                                            unemp_status_t0_jun_lead == 1 ~ 6, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            unemp_status_t0_aug_lead == 1 ~ 8, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            unemp_status_t0_oct_lead == 1 ~ 10, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            unemp_status_t0_dec_lead == 1 ~ 12,
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 1 ~ 1,
                                                            is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 1 ~ 2,
                                                            is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 1 ~ 3,
                                                            is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 1 ~ 4,
                                                            is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 1 ~ 5,
                                                            is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 0, emp_status_t1_dec == 0 ~ 1, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 1,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 2,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 3,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 4,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 5,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 6,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 7,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 8,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 9,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 10,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 11,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 12,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 13,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 14,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 15,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 16,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 17,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 18,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 19,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 20,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 21,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 22,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 23, 
                                                     TRUE ~ 24),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 1,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 2,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 3,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 4,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 5,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 6,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 7,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 8,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 9,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 10,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 11,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 12,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 13,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 14,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 15,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 16,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 17,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 18,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 19,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 20,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 21,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 22,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 23, 
                                                         TRUE ~ 24)) -> data_int_jan
data_int_feb %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_feb_lead == 1 ~ 2,
                                                            unemp_status_t0_mar_lead == 1 ~ 3, unemp_status_t0_apr_lead == 1 ~ 4,
                                                            unemp_status_t0_may_lead == 1 ~ 5, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            unemp_status_t0_jul_lead == 1 ~ 7, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            unemp_status_t0_sep_lead == 1 ~ 9, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            unemp_status_t0_nov_lead == 1 ~ 11, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            unemp_status_t1_jan_lead == 1 ~ 13, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 1 ~ 2,
                                                            is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 1 ~ 3,
                                                            is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 1 ~ 4,
                                                            is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 1 ~ 5,
                                                            is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 1 ~ 1, unemp_status_t0_jan_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 1,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 2,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 3,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 4,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 5,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 6,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 7,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 8,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 9,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 10,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 11,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 12,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 13,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 14,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 15,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 16,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 17,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 18,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 19,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 20,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 21,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 22,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 23,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 24, 
                                                     TRUE ~ 25),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 1,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 2,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 3,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 4,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 5,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 6,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 7,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 8,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 9,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 10,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 11,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 12,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 13,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 14,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 15,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 16,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 17,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 18,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 19,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 20,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 21,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 22,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 23,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 24, 
                                                         TRUE ~ 25)) -> data_int_feb
data_int_mar %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_mar_lead == 1 ~ 3,
                                                            unemp_status_t0_apr_lead == 1 ~ 4, unemp_status_t0_may_lead == 1 ~ 5,
                                                            unemp_status_t0_jun_lead == 1 ~ 6, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            unemp_status_t0_aug_lead == 1 ~ 8, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            unemp_status_t0_oct_lead == 1 ~ 10, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            unemp_status_t0_dec_lead == 1 ~ 12, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            unemp_status_t1_feb_lead == 1 ~ 14,
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 1 ~ 3,
                                                            is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 1 ~ 4,
                                                            is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 1 ~ 5,
                                                            is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 1 ~ 1, unemp_status_t0_feb_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 1,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 2,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 3,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 4,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 5,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 6,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 7,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 8,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 9,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 10,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 11,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 12,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 13,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 14,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 15,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 16,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 17,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 18,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 19,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 20,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 21,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 22,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 23,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 24,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 25, 
                                                     TRUE ~ 26),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 1,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 2,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 3,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 4,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 5,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 6,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 7,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 8,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 9,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 10,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 11,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 12,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 13,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 14,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 15,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 16,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 17,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 18,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 19,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 20,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 21,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 22,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 23,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 24,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 25, 
                                                         TRUE ~ 26)) -> data_int_mar
data_int_apr %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_apr_lead == 1 ~ 4,
                                                            unemp_status_t0_may_lead == 1 ~ 5, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            unemp_status_t0_jul_lead == 1 ~ 7, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            unemp_status_t0_sep_lead == 1 ~ 9, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            unemp_status_t0_nov_lead == 1 ~ 11, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            unemp_status_t1_jan_lead == 1 ~ 13, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            unemp_status_t1_mar_lead == 1 ~ 15,
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 1 ~ 4,
                                                            is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 1 ~ 5,
                                                            is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 1 ~ 1, unemp_status_t0_mar_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 2,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 3,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 4,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 5,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 6,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 7,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 8,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 9,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 10,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 11,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 12,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 13,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 14,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 15,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 16,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 17,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 18,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 19,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 20,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 21,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 22,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 23,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 24,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 25,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 26, 
                                                     TRUE ~ 27),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 2,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 3,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 4,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 5,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 6,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 7,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 8,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 9,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 10,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 11,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 12,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 13,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 14,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 15,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 16,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 17,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 18,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 19,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 20,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 21,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 22,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 23,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 24,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 25,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 26, 
                                                         TRUE ~ 27)) -> data_int_apr
data_int_may %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_may_lead == 1 ~ 5,
                                                            unemp_status_t0_jun_lead == 1 ~ 6, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            unemp_status_t0_aug_lead == 1 ~ 8, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            unemp_status_t0_oct_lead == 1 ~ 10, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            unemp_status_t0_dec_lead == 1 ~ 12, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            unemp_status_t1_feb_lead == 1 ~ 14, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            unemp_status_t1_apr_lead == 1 ~ 16, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 1 ~ 5,
                                                            is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 1 ~ 1, unemp_status_t0_apr_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 3,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 4,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 5,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 6,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 7,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 8,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 9,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 10,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 11,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 12,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 13,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 14,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 15,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 16,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 17,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 18,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 19,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 20,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 21,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 22,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 23,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 24,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 25,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 26,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 27, 
                                                     TRUE ~ 28),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 3,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 4,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 5,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 6,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 7,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 8,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 9,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 10,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 11,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 12,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 13,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 14,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 15,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 16,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 17,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 18,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 19,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 20,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 21,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 22,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 23,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 24,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 25,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 26,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 27, 
                                                         TRUE ~ 28)) -> data_int_may
data_int_jun %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_jun_lead == 1 ~ 6,
                                                            unemp_status_t0_jul_lead == 1 ~ 7, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            unemp_status_t0_sep_lead == 1 ~ 9, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            unemp_status_t0_nov_lead == 1 ~ 11, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            unemp_status_t1_jan_lead == 1 ~ 13, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            unemp_status_t1_mar_lead == 1 ~ 15, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            unemp_status_t1_may_lead == 1 ~ 17, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 1 ~ 6,
                                                            is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 1 ~ 17,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 1 ~ 1, unemp_status_t0_may_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 3,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 4,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 5,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 6,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 7,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 8,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 9,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 10,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 11,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 12,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 13,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 14,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 15,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 16,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 17,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 18,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 19,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 20,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 21,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 22,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 23,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 24,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 25,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 26,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 27,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 28, 
                                                     TRUE ~ 29),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 3,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 4,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 5,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 6,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 7,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 8,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 9,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 10,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 11,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 12,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 13,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 14,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 15,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 16,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 17,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 18,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 19,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 20,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 21,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 22,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 23,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 24,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 25,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 26,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 27,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 28, 
                                                         TRUE ~ 29)) -> data_int_jun
data_int_jul %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_jul_lead == 1 ~ 7,
                                                            unemp_status_t0_aug_lead == 1 ~ 8, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            unemp_status_t0_oct_lead == 1 ~ 10, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            unemp_status_t0_dec_lead == 1 ~ 12, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            unemp_status_t1_feb_lead == 1 ~ 14, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            unemp_status_t1_apr_lead == 1 ~ 16, unemp_status_t1_may_lead == 1 ~ 17,
                                                            unemp_status_t1_jun_lead == 1 ~ 18, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 7,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 1 ~ 17,
                                                            is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 1 ~ 1, unemp_status_t0_jun_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 3,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 4,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 5,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 6,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 7,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 8,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 9,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 10,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 11,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 12,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 13,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 14,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 15,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 16,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 17,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 18,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 19,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 20,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 21,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 22,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 23,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 24,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 25,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 26,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 27,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 28,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 29, 
                                                     TRUE ~ 30),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 3,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 4,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 5,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 6,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 7,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 8,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 9,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 10,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 11,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 12,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 13,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 14,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 15,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 16,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 17,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 18,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 19,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 20,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 21,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 22,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 23,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 24,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 25,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 26,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 27,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 28,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 29, 
                                                         TRUE ~ 30)) -> data_int_jul
data_int_aug %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_aug_lead == 1 ~ 8,
                                                            unemp_status_t0_sep_lead == 1 ~ 9, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            unemp_status_t0_nov_lead == 1 ~ 11, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            unemp_status_t1_jan_lead == 1 ~ 13, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            unemp_status_t1_mar_lead == 1 ~ 15, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            unemp_status_t1_may_lead == 1 ~ 17, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            unemp_status_t1_jul_lead == 1 ~ 19, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 8,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 1 ~ 17,
                                                            is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 1 ~ 19,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 1 ~ 1, unemp_status_t0_jul_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 3,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 4,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 5,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 6,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 7,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 8,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 9,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 10,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 11,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 12,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 13,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 14,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 15,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 16,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 17,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 18,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 19,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 20,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 21,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 22,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 23,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 24,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 25,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 26,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 27,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 28,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 29,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 30, 
                                                     TRUE ~ 31),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 3,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 4,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 5,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 6,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 7,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 8,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 9,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 10,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 11,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 12,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 13,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 14,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 15,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 16,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 17,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 18,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 19,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 20,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 21,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 22,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 23,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 24,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 25,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 26,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 27,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 28,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 29,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 30, 
                                                         TRUE ~ 31)) -> data_int_aug
data_int_sep %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_sep_lead == 1 ~ 9,
                                                            unemp_status_t0_oct_lead == 1 ~ 10, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            unemp_status_t0_dec_lead == 1 ~ 12, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            unemp_status_t1_feb_lead == 1 ~ 14, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            unemp_status_t1_apr_lead == 1 ~ 16, unemp_status_t1_may_lead == 1 ~ 17,
                                                            unemp_status_t1_jun_lead == 1 ~ 18, unemp_status_t1_jul_lead == 1 ~ 19,
                                                            unemp_status_t1_aug_lead == 1 ~ 20, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 9,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 1 ~ 17,
                                                            is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 1 ~ 19,
                                                            is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 1 ~ 20,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 1 ~ 1, unemp_status_t0_aug_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 3,
                                                     is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 4,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 5,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 6,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 7,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 8,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 9,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 10,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 11,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 12,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 13,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 14,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 15,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 16,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 17,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 18,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 19,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 20,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 21,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 22,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 23,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 24,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 25,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 26,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 27,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 28,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 29,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 30,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 31, 
                                                     TRUE ~ 32),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 3,
                                                         is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 4,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 5,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 6,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 7,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 8,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 9,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 10,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 11,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 12,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 13,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 14,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 15,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 16,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 17,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 18,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 19,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 20,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 21,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 22,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 23,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 24,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 25,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 26,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 27,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 28,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 29,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 30,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 31, 
                                                         TRUE ~ 32)) -> data_int_sep
data_int_oct %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_oct_lead == 1 ~ 10,
                                                            unemp_status_t0_nov_lead == 1 ~ 11, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            unemp_status_t1_jan_lead == 1 ~ 13, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            unemp_status_t1_mar_lead == 1 ~ 15, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            unemp_status_t1_may_lead == 1 ~ 17, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            unemp_status_t1_jul_lead == 1 ~ 19, unemp_status_t1_aug_lead == 1 ~ 20,
                                                            unemp_status_t1_sep_lead == 1 ~ 21, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 10,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 1 ~ 17,
                                                            is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 1 ~ 19,
                                                            is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 1 ~ 20,
                                                            is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 1 ~ 21,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 1 ~ 1, unemp_status_t0_sep_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 3,
                                                     is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 4,
                                                     is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 5,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 6,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 7,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 8,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 9,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 10,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 11,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 12,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 13,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 14,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 15,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 16,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 17,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 18,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 19,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 20,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 21,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 22,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 23,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 24,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 25,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 26,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 27,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 28,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 29,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 30,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 31,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 32, 
                                                     TRUE ~ 33),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 3,
                                                         is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 4,
                                                         is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 5,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 6,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 7,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 8,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 9,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 10,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 11,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 12,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 13,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 14,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 15,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 16,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 17,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 18,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 19,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 20,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 21,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 22,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 23,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 24,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 25,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 26,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 27,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 28,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 29,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 30,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 31,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 32, 
                                                         TRUE ~ 33)) -> data_int_oct
data_int_nov %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_nov_lead == 1 ~ 11,
                                                            unemp_status_t0_dec_lead == 1 ~ 12, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            unemp_status_t1_feb_lead == 1 ~ 14, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            unemp_status_t1_apr_lead == 1 ~ 16, unemp_status_t1_may_lead == 1 ~ 17,
                                                            unemp_status_t1_jun_lead == 1 ~ 18, unemp_status_t1_jul_lead == 1 ~ 19,
                                                            unemp_status_t1_aug_lead == 1 ~ 20, unemp_status_t1_sep_lead == 1 ~ 21,
                                                            unemp_status_t1_oct_lead == 1 ~ 22, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 11,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 1 ~ 17,
                                                            is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 1 ~ 19,
                                                            is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 1 ~ 20,
                                                            is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 1 ~ 21,
                                                            is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 1 ~ 22,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 1 ~ 1, unemp_status_t0_oct_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 3,
                                                     is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 4,
                                                     is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 5,
                                                     is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 6,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 7,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 8,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 9,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 10,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 11,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 12,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 13,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 14,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 15,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 16,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 17,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 18,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 19,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 20,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 21,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 22,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 23,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 24,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 25,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 26,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 27,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 28,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 29,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 30,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 31,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 32,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 33, 
                                                     TRUE ~ 34),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 3,
                                                         is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 4,
                                                         is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 5,
                                                         is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 6,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 7,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 8,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 9,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 10,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 11,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 12,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 13,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 14,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 15,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 16,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 17,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 18,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 19,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 20,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 21,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 22,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 23,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 24,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 25,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 26,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 27,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 28,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 29,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 30,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 31,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 32,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 33, 
                                                         TRUE ~ 34)) -> data_int_nov
data_int_dec %>% mutate(first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_period_na == 1 & unemp_status_lead == 1)  ~ NA_integer_,
                                                            unemp_status_t0_dec_lead == 1 ~ 12,
                                                            unemp_status_t1_jan_lead == 1 ~ 13, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            unemp_status_t1_mar_lead == 1 ~ 15, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            unemp_status_t1_may_lead == 1 ~ 17, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            unemp_status_t1_jul_lead == 1 ~ 19, unemp_status_t1_aug_lead == 1 ~ 20,
                                                            unemp_status_t1_sep_lead == 1 ~ 21, unemp_status_t1_oct_lead == 1 ~ 22,
                                                            unemp_status_t1_nov_lead == 1 ~ 23, 
                                                            TRUE ~ NA_integer_),
                        first_mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | (unemp_status_lead == 1 & unemp_status_period_na == 0) ~ first_mo_unemp_post_int,
                                                            is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 1 ~ 12,
                                                            is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 1 ~ 13,
                                                            is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 1 ~ 14,
                                                            is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 1 ~ 15,
                                                            is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 1 ~ 16,
                                                            is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 1 ~ 17,
                                                            is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 1 ~ 18,
                                                            is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 1 ~ 19,
                                                            is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 1 ~ 20,
                                                            is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 1 ~ 21,
                                                            is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 1 ~ 22,
                                                            is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 1 ~ 23,
                                                            TRUE ~ NA_integer_),
                        already_unemp_pre_int = case_when(is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 1 ~ 1, unemp_status_t0_nov_lead == 0 ~ 0, 
                                                          TRUE ~ NA_integer_),
                        mo_unemp_pre_int = case_when(is.na(already_unemp_pre_int) ~ NA_integer_,
                                                     already_unemp_pre_int == 0 ~ 0,
                                                     is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 1,
                                                     is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 2,
                                                     is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 3,
                                                     is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 4,
                                                     is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 5,
                                                     is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 6,
                                                     is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 7,
                                                     is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 8,
                                                     is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 9,
                                                     is.na(unemp_status_t0_jan_lead) ~ NA_integer_, unemp_status_t0_jan_lead == 0 ~ 10,
                                                     is.na(emp_status_t1_dec) ~ NA_integer_, emp_status_t1_dec == 1 ~ 11,
                                                     is.na(emp_status_t1_nov) ~ NA_integer_, emp_status_t1_nov == 1 ~ 12,
                                                     is.na(emp_status_t1_oct) ~ NA_integer_, emp_status_t1_oct == 1 ~ 13,
                                                     is.na(emp_status_t1_sep) ~ NA_integer_, emp_status_t1_sep == 1 ~ 14,
                                                     is.na(emp_status_t1_aug) ~ NA_integer_, emp_status_t1_aug == 1 ~ 15,
                                                     is.na(emp_status_t1_jul) ~ NA_integer_, emp_status_t1_jul == 1 ~ 16,
                                                     is.na(emp_status_t1_jun) ~ NA_integer_, emp_status_t1_jun == 1 ~ 17,
                                                     is.na(emp_status_t1_may) ~ NA_integer_, emp_status_t1_may == 1 ~ 18,
                                                     is.na(emp_status_t1_apr) ~ NA_integer_, emp_status_t1_apr == 1 ~ 19,
                                                     is.na(emp_status_t1_mar) ~ NA_integer_, emp_status_t1_mar == 1 ~ 20,
                                                     is.na(emp_status_t1_feb) ~ NA_integer_, emp_status_t1_feb == 1 ~ 21,
                                                     is.na(emp_status_t1_jan) ~ NA_integer_, emp_status_t1_jan == 1 ~ 22,
                                                     is.na(emp_status_t2_dec) ~ NA_integer_, emp_status_t2_dec == 1 ~ 23,
                                                     is.na(emp_status_t2_nov) ~ NA_integer_, emp_status_t2_nov == 1 ~ 24,
                                                     is.na(emp_status_t2_oct) ~ NA_integer_, emp_status_t2_oct == 1 ~ 25,
                                                     is.na(emp_status_t2_sep) ~ NA_integer_, emp_status_t2_sep == 1 ~ 26,
                                                     is.na(emp_status_t2_aug) ~ NA_integer_, emp_status_t2_aug == 1 ~ 27,
                                                     is.na(emp_status_t2_jul) ~ NA_integer_, emp_status_t2_jul == 1 ~ 28,
                                                     is.na(emp_status_t2_jun) ~ NA_integer_, emp_status_t2_jun == 1 ~ 29,
                                                     is.na(emp_status_t2_may) ~ NA_integer_, emp_status_t2_may == 1 ~ 30,
                                                     is.na(emp_status_t2_apr) ~ NA_integer_, emp_status_t2_apr == 1 ~ 31,
                                                     is.na(emp_status_t2_mar) ~ NA_integer_, emp_status_t2_mar == 1 ~ 32,
                                                     is.na(emp_status_t2_feb) ~ NA_integer_, emp_status_t2_feb == 1 ~ 33,
                                                     is.na(emp_status_t2_jan) ~ NA_integer_, emp_status_t2_jan == 1 ~ 34, 
                                                     TRUE ~ 35),
                        mo_unemp_pre_int_min = case_when(already_unemp_pre_int == 0 | is.na(already_unemp_pre_int) ~ 0,
                                                         is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 1,
                                                         is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 2,
                                                         is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 3,
                                                         is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 4,
                                                         is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 5,
                                                         is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 6,
                                                         is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 7,
                                                         is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 8,
                                                         is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 9,
                                                         is.na(unemp_status_t0_jan_lead) | unemp_status_t0_jan_lead == 0 ~ 10,
                                                         is.na(emp_status_t1_dec) | emp_status_t1_dec == 1 ~ 11,
                                                         is.na(emp_status_t1_nov) | emp_status_t1_nov == 1 ~ 12,
                                                         is.na(emp_status_t1_oct) | emp_status_t1_oct == 1 ~ 13,
                                                         is.na(emp_status_t1_sep) | emp_status_t1_sep == 1 ~ 14,
                                                         is.na(emp_status_t1_aug) | emp_status_t1_aug == 1 ~ 15,
                                                         is.na(emp_status_t1_jul) | emp_status_t1_jul == 1 ~ 16,
                                                         is.na(emp_status_t1_jun) | emp_status_t1_jun == 1 ~ 17,
                                                         is.na(emp_status_t1_may) | emp_status_t1_may == 1 ~ 18,
                                                         is.na(emp_status_t1_apr) | emp_status_t1_apr == 1 ~ 19,
                                                         is.na(emp_status_t1_mar) | emp_status_t1_mar == 1 ~ 20,
                                                         is.na(emp_status_t1_feb) | emp_status_t1_feb == 1 ~ 21,
                                                         is.na(emp_status_t1_jan) | emp_status_t1_jan == 1 ~ 22,
                                                         is.na(emp_status_t2_dec) | emp_status_t2_dec == 1 ~ 23,
                                                         is.na(emp_status_t2_nov) | emp_status_t2_nov == 1 ~ 24,
                                                         is.na(emp_status_t2_oct) | emp_status_t2_oct == 1 ~ 25,
                                                         is.na(emp_status_t2_sep) | emp_status_t2_sep == 1 ~ 26,
                                                         is.na(emp_status_t2_aug) | emp_status_t2_aug == 1 ~ 27,
                                                         is.na(emp_status_t2_jul) | emp_status_t2_jul == 1 ~ 28,
                                                         is.na(emp_status_t2_jun) | emp_status_t2_jun == 1 ~ 29,
                                                         is.na(emp_status_t2_may) | emp_status_t2_may == 1 ~ 30,
                                                         is.na(emp_status_t2_apr) | emp_status_t2_apr == 1 ~ 31,
                                                         is.na(emp_status_t2_mar) | emp_status_t2_mar == 1 ~ 32,
                                                         is.na(emp_status_t2_feb) | emp_status_t2_feb == 1 ~ 33,
                                                         is.na(emp_status_t2_jan) | emp_status_t2_jan == 1 ~ 34, 
                                                         TRUE ~ 35)) -> data_int_dec

#Re-merging interview month-specific datasets to create full dataset  
data <- rbind(data_pre2001,
              data_int_jan, data_int_feb, data_int_mar, data_int_apr, data_int_may, data_int_jun,
              data_int_jul, data_int_aug, data_int_sep, data_int_oct, data_int_nov, data_int_dec)
data %>% arrange(unique_id, year) -> data
rm(data_pre2001,
   data_int_jan, data_int_feb, data_int_mar, data_int_apr, data_int_may, data_int_jun,
   data_int_jul, data_int_aug, data_int_sep, data_int_oct, data_int_nov, data_int_dec)

#Creating 'mo_unemp_pre_int_na', which = 0 if all relevant monthly employment data is known, and =1 otherwise 
#NOTE - Retrospective monthly unemployment information for T-2 is from 2003 onwards (so is only known for the 'interview year' in 2001)
#       and monthly unemployment information for T-1 is only from 2003 onwards as well (so is not available for T-1 of the interview year in 2001)
data %>% mutate(mo_unemp_pre_int_na = case_when((!is.na(already_unemp_pre_int) & already_unemp_pre_int == 0) | year < 2001 ~ NA_integer_,
                                                is.na(mo_unemp_pre_int) | is.na(already_unemp_pre_int) ~ 1,
                                                interview_month == 1 & mo_unemp_pre_int_min == 24 ~ 1,
                                                interview_month == 2 & mo_unemp_pre_int_min == 25 ~ 1,
                                                interview_month == 3 & mo_unemp_pre_int_min == 26 ~ 1,
                                                interview_month == 4 & mo_unemp_pre_int_min == 27 ~ 1,
                                                interview_month == 5 & mo_unemp_pre_int_min == 28 ~ 1,
                                                interview_month == 6 & mo_unemp_pre_int_min == 29 ~ 1,
                                                interview_month == 7 & mo_unemp_pre_int_min == 30 ~ 1,
                                                interview_month == 8 & mo_unemp_pre_int_min == 31 ~ 1,
                                                interview_month == 9 & mo_unemp_pre_int_min == 32 ~ 1,
                                                interview_month == 10 & mo_unemp_pre_int_min == 33 ~ 1,
                                                interview_month == 11 & mo_unemp_pre_int_min == 34 ~ 1,
                                                interview_month == 12 & mo_unemp_pre_int_min == 35 ~ 1,
                                                mo_unemp_pre_int == mo_unemp_pre_int_min ~ 0,
                                                TRUE ~ 99)) -> data #Noone should be left uncoded, so TRUE ~ 99 indicates an issue

#### Step 5 - Creating T+2 and T+3 monthly unemployment status variables ####
#NOTE - As this information is collected retrospectively, we only have information from T0 and T+1 up to 2017, so will only have information on T+2 and T+3 up to 2015 (requiring 2 post-waves)
data <- data %>% arrange(unique_id, year)
library(data.table)
data <- as.data.table(data)
data[, ':=' (unemp_status_t2_jan_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_jan,2)), NA_integer_, fifelse(lead(emp_status_t2_jan,2) == 1, 0, 1)),
             unemp_status_t2_feb_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_feb,2)), NA_integer_, fifelse(lead(emp_status_t2_feb,2) == 1, 0, 1)),
             unemp_status_t2_mar_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_mar,2)), NA_integer_, fifelse(lead(emp_status_t2_mar,2) == 1, 0, 1)),
             unemp_status_t2_apr_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_apr,2)), NA_integer_, fifelse(lead(emp_status_t2_apr,2) == 1, 0, 1)),
             unemp_status_t2_may_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_may,2)), NA_integer_, fifelse(lead(emp_status_t2_may,2) == 1, 0, 1)),
             unemp_status_t2_jun_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_jun,2)), NA_integer_, fifelse(lead(emp_status_t2_jun,2) == 1, 0, 1)),
             unemp_status_t2_jul_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_jul,2)), NA_integer_, fifelse(lead(emp_status_t2_jul,2) == 1, 0, 1)),
             unemp_status_t2_aug_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_aug,2)), NA_integer_, fifelse(lead(emp_status_t2_aug,2) == 1, 0, 1)),
             unemp_status_t2_sep_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_sep,2)), NA_integer_, fifelse(lead(emp_status_t2_sep,2) == 1, 0, 1)),
             unemp_status_t2_oct_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_oct,2)), NA_integer_, fifelse(lead(emp_status_t2_oct,2) == 1, 0, 1)),
             unemp_status_t2_nov_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_nov,2)), NA_integer_, fifelse(lead(emp_status_t2_nov,2) == 1, 0, 1)),
             unemp_status_t2_dec_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t2_dec,2)), NA_integer_, fifelse(lead(emp_status_t2_dec,2) == 1, 0, 1)),
             unemp_status_t3_jan_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_jan,2)), NA_integer_, fifelse(lead(emp_status_t1_jan,2) == 1, 0, 1)),
             unemp_status_t3_feb_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_feb,2)), NA_integer_, fifelse(lead(emp_status_t1_feb,2) == 1, 0, 1)),
             unemp_status_t3_mar_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_mar,2)), NA_integer_, fifelse(lead(emp_status_t1_mar,2) == 1, 0, 1)),
             unemp_status_t3_apr_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_apr,2)), NA_integer_, fifelse(lead(emp_status_t1_apr,2) == 1, 0, 1)),
             unemp_status_t3_may_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_may,2)), NA_integer_, fifelse(lead(emp_status_t1_may,2) == 1, 0, 1)),
             unemp_status_t3_jun_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_jun,2)), NA_integer_, fifelse(lead(emp_status_t1_jun,2) == 1, 0, 1)),
             unemp_status_t3_jul_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_jul,2)), NA_integer_, fifelse(lead(emp_status_t1_jul,2) == 1, 0, 1)),
             unemp_status_t3_aug_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_aug,2)), NA_integer_, fifelse(lead(emp_status_t1_aug,2) == 1, 0, 1)),
             unemp_status_t3_sep_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_sep,2)), NA_integer_, fifelse(lead(emp_status_t1_sep,2) == 1, 0, 1)),
             unemp_status_t3_oct_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_oct,2)), NA_integer_, fifelse(lead(emp_status_t1_oct,2) == 1, 0, 1)),
             unemp_status_t3_nov_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_nov,2)), NA_integer_, fifelse(lead(emp_status_t1_nov,2) == 1, 0, 1)),
             unemp_status_t3_dec_lead = fifelse(lead(year,2) != (year+4) | is.na(lead(emp_status_t1_dec,2)), NA_integer_, fifelse(lead(emp_status_t1_dec,2) == 1, 0, 1))), by = .(unique_id)]
data <- as.data.frame(data)
detach("package:data.table", unload = TRUE)

#### Step 6 - Working out months unemployed post-int and total months unemployed post-unemployment before re-employment ####
#             Creating 1) 'mo_unemp_post_int' for known months unemployed from first month unemployed post-interview (incl. interview month) if unemployed at interview
#                      2) 'mo_unemp_post_int_min' for minimum months unemployed post-interview (incl. interview month) if unemployed at interview, based on known monthly unemployment data 
#                      3) 'mo_unemp_post_int_na' indicating whether all relevant unemployment information was known
#                      4) 'mo_post_unemp_until_emp' for continuous months unemployed before becoming re-employed, if unemployed in the 12 months post-interview
#                      5) 'mo_post_unemp_until_emp_na' for whether exact months were known, or if not why
#     - For 'mo_unemp_until_emp_na', assign = NA if employed for 12 months post-interview or if this is unknown, 
#                                           = "actual" if all relevant unemployment information is known
#                                           = "cens_left" if information is only left-censored (where already unemployed pre-interview)
#                                           = "cens_right" if information is only right-censored (i.e. missing relevant information post-interview), or 
#                                           = "cens_both" if left and right-censored (missing or censored both pre- and post-interview)

#Splitting dataset by first month unemployed post-interview
#NOTE - Years <2001, not unemployed in the 12 months post-interview, and without known first month unemployed post-interview are assigned NA 
data_pre2001_or_emp <- data[data$year < 2001 | is.na(data$first_mo_unemp_post_int) | (!is.na(data$unemp_status_lead) & data$unemp_status_lead == 0) | (!is.na(data$unemp_status_lead) & data$unemp_status_lead == 1 & is.na(data$first_mo_unemp_post_int)), ]
data_first_unemp_m1 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 1, ]
data_first_unemp_m2 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 2, ]
data_first_unemp_m3 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 3, ]
data_first_unemp_m4 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 4, ]
data_first_unemp_m5 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 5, ]
data_first_unemp_m6 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 6, ]
data_first_unemp_m7 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 7, ]
data_first_unemp_m8 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 8, ]
data_first_unemp_m9 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 9, ]
data_first_unemp_m10 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 10, ]
data_first_unemp_m11 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 11, ]
data_first_unemp_m12 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 12, ]
data_first_unemp_m13 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 13, ]
data_first_unemp_m14 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 14, ]
data_first_unemp_m15 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 15, ]
data_first_unemp_m16 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 16, ]
data_first_unemp_m17 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 17, ]
data_first_unemp_m18 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 18, ]
data_first_unemp_m19 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 19, ]
data_first_unemp_m20 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 20, ]
data_first_unemp_m21 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 21, ]
data_first_unemp_m22 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 22, ]
data_first_unemp_m23 <- data[!is.na(data$first_mo_unemp_post_int) & data$first_mo_unemp_post_int == 23, ]

#Creating 'mo_unemp_post_int', 'mo_unemp_post_int_min' variables
data_pre2001_or_emp %>% mutate(mo_unemp_post_int = NA_integer_, mo_unemp_post_int_min = NA_integer_) -> data_pre2001_or_emp
data_first_unemp_m1 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_feb_lead) ~ NA_integer_, unemp_status_t0_feb_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 12,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 16, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 17, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 18, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 19, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 20, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 21, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 22, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 23,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 24,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 28, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 29, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 30, 
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 31, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 32, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 33, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 34, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 35,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 36,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 39, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 40, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 41, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 42, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 43, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 44, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 45, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 46,
                                                             TRUE ~ 47),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_feb_lead) | unemp_status_t0_feb_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 12,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 16, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 17, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 18, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 19, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 20, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 21, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 22, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 23,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 24,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 28, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 29, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 30, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 31, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 32, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 33, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 34, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 35,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 36,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 39, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 40, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 41, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 42, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 43, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 44, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 45, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 46,
                                                                 TRUE ~ 47)) -> data_first_unemp_m1
data_first_unemp_m2 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_mar_lead) ~ NA_integer_, unemp_status_t0_mar_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 11,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 16, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 17, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 18, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 19, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 20, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 21, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 22,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 23,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 24, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 28, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 29, 
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 30, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 31, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 32, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 33, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 34,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 35,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 39, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 40, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 41, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 42, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 43, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 44, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 45,
                                                             TRUE ~ 46),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_mar_lead) | unemp_status_t0_mar_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 11,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 16, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 17, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 18, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 19, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 20, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 21, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 22,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 23,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 28, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 29, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 30, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 31, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 32, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 33, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 34,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 35,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 39, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 40, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 41, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 42, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 43, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 44, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 45,
                                                                 TRUE ~ 46)) -> data_first_unemp_m2
data_first_unemp_m3 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_apr_lead) ~ NA_integer_, unemp_status_t0_apr_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 10,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 16, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 17, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 18, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 19, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 20, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 21,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 22,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 23, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 24, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 28, 
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 29, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 30, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 31, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 32, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 33,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 34,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 35, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 39, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 40, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 41, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 42, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 43, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 44,
                                                             TRUE ~ 45),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_apr_lead) | unemp_status_t0_apr_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 10,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 16, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 17, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 18, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 19, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 20, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 21,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 22,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 23, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 28, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 29, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 30, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 31, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 32, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 33,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 34,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 35, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 39, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 40, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 41, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 42, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 43, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 44,
                                                                 TRUE ~ 45)) -> data_first_unemp_m3
data_first_unemp_m4 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_may_lead) ~ NA_integer_, unemp_status_t0_may_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 9,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 16, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 17, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 18, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 19, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 20,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 21,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 22, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 23, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 24, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 28, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 29, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 30, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 31, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 32,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 33,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 34, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 35, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 39, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 40, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 41, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 42, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 43,
                                                             TRUE ~ 44),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_may_lead) | unemp_status_t0_may_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 9,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 16, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 17, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 18, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 19, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 20,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 21,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 22, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 23, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 28, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 29, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 30, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 31, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 32,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 33,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 34, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 35, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 39, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 40, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 41, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 42, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 43,
                                                                 TRUE ~ 44)) -> data_first_unemp_m4
data_first_unemp_m5 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_jun_lead) ~ NA_integer_, unemp_status_t0_jun_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 8,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 16, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 17, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 18, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 19,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 20,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 21, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 22, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 23, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 24, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 26,
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 28, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 29, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 30, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 31,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 32,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 33, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 34, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 35, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 39, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 40, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 41, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 42,
                                                             TRUE ~ 43),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_jun_lead) | unemp_status_t0_jun_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 8,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 16, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 17, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 18, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 19,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 20,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 21, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 22, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 23, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 28, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 29, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 30, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 31,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 32,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 33, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 34, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 35, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 39, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 40, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 41, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 42,
                                                                 TRUE ~ 43)) -> data_first_unemp_m5
data_first_unemp_m6 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_jul_lead) ~ NA_integer_, unemp_status_t0_jul_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 7,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 16, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 17, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 18,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 19,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 20, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 21, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 22, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 23, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 24, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 25,
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 28, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 29, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 30,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 31,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 32, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 33, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 34, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 35, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 39, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 40, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 41,
                                                             TRUE ~ 42),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_jul_lead) | unemp_status_t0_jul_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 7,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 16, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 17, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 18,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 19,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 20, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 21, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 22, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 23, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 28, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 29, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 30,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 31,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 32, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 33, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 34, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 35, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 39, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 40, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 41,
                                                                 TRUE ~ 42)) -> data_first_unemp_m6
data_first_unemp_m7 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_aug_lead) ~ NA_integer_, unemp_status_t0_aug_lead == 0 ~ 1,
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 6,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 16, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 17,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 18,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 19, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 20, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 21, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 22, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 23, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 24,
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 28, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 29,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 30,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 31, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 32, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 33, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 34, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 35, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 39, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 40,
                                                             TRUE ~ 41),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_aug_lead) | unemp_status_t0_aug_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 6,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 16, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 17,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 18,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 19, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 20, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 21, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 22, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 23, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 28, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 29,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 30,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 31, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 32, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 33, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 34, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 35, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 39, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 40,
                                                                 TRUE ~ 41)) -> data_first_unemp_m7
data_first_unemp_m8 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_sep_lead) ~ NA_integer_, unemp_status_t0_sep_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 4, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 5,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 15, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 16,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 17,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 18, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 19, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 20, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 21, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 22, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 23,
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 24, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 27, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 28,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 29,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 30, 
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 31, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 32, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 33, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 34, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 35, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 38, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 39,
                                                             TRUE ~ 40),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_sep_lead) | unemp_status_t0_sep_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 4, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 5,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 15, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 16,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 17,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 18, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 19, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 20, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 21, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 22, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 23, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 27, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 28,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 29,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 30, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 31, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 32, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 33, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 34, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 35, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 38, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 39,
                                                                 TRUE ~ 40)) -> data_first_unemp_m8
data_first_unemp_m9 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                             is.na(unemp_status_t0_oct_lead) ~ NA_integer_, unemp_status_t0_oct_lead == 0 ~ 1, 
                                                             is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 2, 
                                                             is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 3, 
                                                             is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 4,
                                                             is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 5, 
                                                             is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 6, 
                                                             is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 7, 
                                                             is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 8, 
                                                             is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 9, 
                                                             is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 10, 
                                                             is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 11, 
                                                             is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 12, 
                                                             is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 13, 
                                                             is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 14, 
                                                             is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 15,
                                                             is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 16,
                                                             is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 17, 
                                                             is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 18, 
                                                             is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 19, 
                                                             is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 20, 
                                                             is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 21, 
                                                             is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 22,
                                                             is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 23, 
                                                             is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 24, 
                                                             is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 25, 
                                                             is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 26, 
                                                             is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 27,
                                                             is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 28,
                                                             is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 29,
                                                             is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 30, 
                                                             is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 31, 
                                                             is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 32, 
                                                             is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 33, 
                                                             is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 34, 
                                                             is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 35, 
                                                             is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 36, 
                                                             is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 37, 
                                                             is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 38,
                                                             TRUE ~ 39),
                               mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                 is.na(unemp_status_t0_oct_lead) | unemp_status_t0_oct_lead == 0 ~ 1, 
                                                                 is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 2, 
                                                                 is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 3, 
                                                                 is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 4,
                                                                 is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 5, 
                                                                 is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 6, 
                                                                 is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 7, 
                                                                 is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 8, 
                                                                 is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 9, 
                                                                 is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 10, 
                                                                 is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 11, 
                                                                 is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 12, 
                                                                 is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 13, 
                                                                 is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 14, 
                                                                 is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 15,
                                                                 is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 16,
                                                                 is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 17, 
                                                                 is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 18, 
                                                                 is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 19, 
                                                                 is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 20, 
                                                                 is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 21, 
                                                                 is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 22, 
                                                                 is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 23, 
                                                                 is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 24, 
                                                                 is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 25, 
                                                                 is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 26, 
                                                                 is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 27,
                                                                 is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 28,
                                                                 is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 29, 
                                                                 is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 30, 
                                                                 is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 31, 
                                                                 is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 32, 
                                                                 is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 33, 
                                                                 is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 34, 
                                                                 is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 35, 
                                                                 is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 36, 
                                                                 is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 37, 
                                                                 is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 38,
                                                                 TRUE ~ 39)) -> data_first_unemp_m9
data_first_unemp_m10 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t0_nov_lead) ~ NA_integer_, unemp_status_t0_nov_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 3,
                                                              is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 12,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 14,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 15,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 21,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 26,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 27,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 28,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 30, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 31, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 32, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 33, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 34, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 35, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 36, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 37,
                                                              TRUE ~ 38),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t0_nov_lead) | unemp_status_t0_nov_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 3,
                                                                  is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 14,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 15,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 26,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 27,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 30, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 31, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 32, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 33, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 34, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 35, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 36, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 37,
                                                                  TRUE ~ 38)) -> data_first_unemp_m10
data_first_unemp_m11 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t0_dec_lead) ~ NA_integer_, unemp_status_t0_dec_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 2,
                                                              is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 11,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 13,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 14,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 20,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 25,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 26,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 27,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 30, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 31, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 32, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 33, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 34, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 35, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 36,
                                                              TRUE ~ 37),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t0_dec_lead) | unemp_status_t0_dec_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 2,
                                                                  is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 13,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 14,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 25,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 26,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 30, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 31, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 32, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 33, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 34, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 35, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 36,
                                                                  TRUE ~ 37)) -> data_first_unemp_m11
data_first_unemp_m12 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_jan_lead) ~ NA_integer_, unemp_status_t1_jan_lead == 0 ~ 1,
                                                              is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 10,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 12,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 13,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 19,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 24,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 25,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 26,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 30, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 31, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 32, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 33, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 34, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 35,
                                                              TRUE ~ 36),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_jan_lead) | unemp_status_t1_jan_lead == 0 ~ 1,
                                                                  is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 12,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 13,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 24,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 25,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 30, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 31, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 32, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 33, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 34, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 35,
                                                                  TRUE ~ 36)) -> data_first_unemp_m12
data_first_unemp_m13 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_feb_lead) ~ NA_integer_, unemp_status_t1_feb_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 9,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 11,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 12,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 18,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 23,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 24,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 25,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 30, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 31, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 32, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 33, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 34,
                                                              TRUE ~ 35),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_feb_lead) | unemp_status_t1_feb_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 11,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 12,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 23,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 24,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 30, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 31, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 32, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 33, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 34,
                                                                  TRUE ~ 35)) -> data_first_unemp_m13
data_first_unemp_m14 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_mar_lead) ~ NA_integer_, unemp_status_t1_mar_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 8,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 10,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 11,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 17,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 22,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 23,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 24,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 30, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 31, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 32, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 33,
                                                              TRUE ~ 34),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_mar_lead) | unemp_status_t1_mar_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 10,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 11,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 22,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 23,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 30, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 31, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 32, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 33,
                                                                  TRUE ~ 34)) -> data_first_unemp_m14
data_first_unemp_m15 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_apr_lead) ~ NA_integer_, unemp_status_t1_apr_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 7,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 9,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 10,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 16,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 21,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 22,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 23,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 30, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 31, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 32,
                                                              TRUE ~ 33),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_apr_lead) | unemp_status_t1_apr_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 9,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 10,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 21,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 22,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 30, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 31, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 32,
                                                                  TRUE ~ 33)) -> data_first_unemp_m15
data_first_unemp_m16 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_may_lead) ~ NA_integer_, unemp_status_t1_may_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 6,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 8,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 9,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 15,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 20,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 21,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 22,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 30, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 31,
                                                              TRUE ~ 32),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_may_lead) | unemp_status_t1_may_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 8,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 9,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 20,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 21,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 30, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 31,
                                                                  TRUE ~ 32)) -> data_first_unemp_m16
data_first_unemp_m17 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_jun_lead) ~ NA_integer_, unemp_status_t1_jun_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 5,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 7,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 8,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 14,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 19,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 20,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 21,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 29, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 30,
                                                              TRUE ~ 31),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_jun_lead) | unemp_status_t1_jun_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 7,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 8,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 19,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 20,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 29, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 30,
                                                                  TRUE ~ 31)) -> data_first_unemp_m17
data_first_unemp_m18 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_jul_lead) ~ NA_integer_, unemp_status_t1_jul_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 4,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 6,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 7,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 13,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 18,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 19,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 20,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 28, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 29,
                                                              TRUE ~ 30),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_jul_lead) | unemp_status_t1_jul_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 6,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 7,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 18,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 19,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 28, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 29,
                                                                  TRUE ~ 30)) -> data_first_unemp_m18
data_first_unemp_m19 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_aug_lead) ~ NA_integer_, unemp_status_t1_aug_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 3,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 5,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 6,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 12,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 17,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 18,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 19,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 27, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 28,
                                                              TRUE ~ 29),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_aug_lead) | unemp_status_t1_aug_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 5,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 6,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 10,
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 17,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 18,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 27, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 28,
                                                                  TRUE ~ 29)) -> data_first_unemp_m19
data_first_unemp_m20 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_sep_lead) ~ NA_integer_, unemp_status_t1_sep_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 2,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 4,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 5,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 11,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 15, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 16,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 17,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 18,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 26, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 27,
                                                              TRUE ~ 28),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_sep_lead) | unemp_status_t1_sep_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 4,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 5,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 9,
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 16,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 17,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 26, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 27,
                                                                  TRUE ~ 28)) -> data_first_unemp_m20
data_first_unemp_m21 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_oct_lead) ~ NA_integer_, unemp_status_t1_oct_lead == 0 ~ 1,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 2, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 3,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 4,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 10,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 14, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 15,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 16,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 17,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 25, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 26,
                                                              TRUE ~ 27),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_oct_lead) | unemp_status_t1_oct_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 2, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 3,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 4,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 8,
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 14, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 15,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 16,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 19,
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 25, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 26,
                                                                  TRUE ~ 27)) -> data_first_unemp_m21
data_first_unemp_m22 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_nov_lead) ~ NA_integer_, unemp_status_t1_nov_lead == 0 ~ 1, 
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 2,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 3,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 8, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 9,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 13, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 14,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 15,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 16,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 24, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 25,
                                                              TRUE ~ 26),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_nov_lead) | unemp_status_t1_nov_lead == 0 ~ 1, 
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 2,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 3,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 6, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 7,
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 13, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 14,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 15,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 17, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 18,
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 24, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 25,
                                                                  TRUE ~ 26)) -> data_first_unemp_m22
data_first_unemp_m23 %>% mutate(mo_unemp_post_int = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                              is.na(unemp_status_t1_dec_lead) ~ NA_integer_, unemp_status_t1_dec_lead == 0 ~ 1,
                                                              is.na(unemp_status_t2_jan_lead) ~ NA_integer_, unemp_status_t2_jan_lead == 0 ~ 2,
                                                              is.na(unemp_status_t2_feb_lead) ~ NA_integer_, unemp_status_t2_feb_lead == 0 ~ 3, 
                                                              is.na(unemp_status_t2_mar_lead) ~ NA_integer_, unemp_status_t2_mar_lead == 0 ~ 4, 
                                                              is.na(unemp_status_t2_apr_lead) ~ NA_integer_, unemp_status_t2_apr_lead == 0 ~ 5, 
                                                              is.na(unemp_status_t2_may_lead) ~ NA_integer_, unemp_status_t2_may_lead == 0 ~ 6, 
                                                              is.na(unemp_status_t2_jun_lead) ~ NA_integer_, unemp_status_t2_jun_lead == 0 ~ 7, 
                                                              is.na(unemp_status_t2_jul_lead) ~ NA_integer_, unemp_status_t2_jul_lead == 0 ~ 8,
                                                              is.na(unemp_status_t2_aug_lead) ~ NA_integer_, unemp_status_t2_aug_lead == 0 ~ 9, 
                                                              is.na(unemp_status_t2_sep_lead) ~ NA_integer_, unemp_status_t2_sep_lead == 0 ~ 10, 
                                                              is.na(unemp_status_t2_oct_lead) ~ NA_integer_, unemp_status_t2_oct_lead == 0 ~ 11, 
                                                              is.na(unemp_status_t2_nov_lead) ~ NA_integer_, unemp_status_t2_nov_lead == 0 ~ 12, 
                                                              is.na(unemp_status_t2_dec_lead) ~ NA_integer_, unemp_status_t2_dec_lead == 0 ~ 13,
                                                              is.na(unemp_status_t3_jan_lead) ~ NA_integer_, unemp_status_t3_jan_lead == 0 ~ 14,
                                                              is.na(unemp_status_t3_feb_lead) ~ NA_integer_, unemp_status_t3_feb_lead == 0 ~ 15,
                                                              is.na(unemp_status_t3_mar_lead) ~ NA_integer_, unemp_status_t3_mar_lead == 0 ~ 16, 
                                                              is.na(unemp_status_t3_apr_lead) ~ NA_integer_, unemp_status_t3_apr_lead == 0 ~ 17, 
                                                              is.na(unemp_status_t3_may_lead) ~ NA_integer_, unemp_status_t3_may_lead == 0 ~ 18, 
                                                              is.na(unemp_status_t3_jun_lead) ~ NA_integer_, unemp_status_t3_jun_lead == 0 ~ 19, 
                                                              is.na(unemp_status_t3_jul_lead) ~ NA_integer_, unemp_status_t3_jul_lead == 0 ~ 20, 
                                                              is.na(unemp_status_t3_aug_lead) ~ NA_integer_, unemp_status_t3_aug_lead == 0 ~ 21, 
                                                              is.na(unemp_status_t3_sep_lead) ~ NA_integer_, unemp_status_t3_sep_lead == 0 ~ 22, 
                                                              is.na(unemp_status_t3_oct_lead) ~ NA_integer_, unemp_status_t3_oct_lead == 0 ~ 23, 
                                                              is.na(unemp_status_t3_nov_lead) ~ NA_integer_, unemp_status_t3_nov_lead == 0 ~ 24,
                                                              TRUE ~ 25),
                                mo_unemp_post_int_min = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 ~ NA_integer_,
                                                                  is.na(unemp_status_t1_dec_lead) | unemp_status_t1_dec_lead == 0 ~ 1,
                                                                  is.na(unemp_status_t2_jan_lead) | unemp_status_t2_jan_lead == 0 ~ 2,
                                                                  is.na(unemp_status_t2_feb_lead) | unemp_status_t2_feb_lead == 0 ~ 3, 
                                                                  is.na(unemp_status_t2_mar_lead) | unemp_status_t2_mar_lead == 0 ~ 4, 
                                                                  is.na(unemp_status_t2_apr_lead) | unemp_status_t2_apr_lead == 0 ~ 5, 
                                                                  is.na(unemp_status_t2_may_lead) | unemp_status_t2_may_lead == 0 ~ 6,
                                                                  is.na(unemp_status_t2_jun_lead) | unemp_status_t2_jun_lead == 0 ~ 7, 
                                                                  is.na(unemp_status_t2_jul_lead) | unemp_status_t2_jul_lead == 0 ~ 8, 
                                                                  is.na(unemp_status_t2_aug_lead) | unemp_status_t2_aug_lead == 0 ~ 9, 
                                                                  is.na(unemp_status_t2_sep_lead) | unemp_status_t2_sep_lead == 0 ~ 10, 
                                                                  is.na(unemp_status_t2_oct_lead) | unemp_status_t2_oct_lead == 0 ~ 11, 
                                                                  is.na(unemp_status_t2_nov_lead) | unemp_status_t2_nov_lead == 0 ~ 12, 
                                                                  is.na(unemp_status_t2_dec_lead) | unemp_status_t2_dec_lead == 0 ~ 13,
                                                                  is.na(unemp_status_t3_jan_lead) | unemp_status_t3_jan_lead == 0 ~ 14,
                                                                  is.na(unemp_status_t3_feb_lead) | unemp_status_t3_feb_lead == 0 ~ 15, 
                                                                  is.na(unemp_status_t3_mar_lead) | unemp_status_t3_mar_lead == 0 ~ 16, 
                                                                  is.na(unemp_status_t3_apr_lead) | unemp_status_t3_apr_lead == 0 ~ 17,
                                                                  is.na(unemp_status_t3_may_lead) | unemp_status_t3_may_lead == 0 ~ 18, 
                                                                  is.na(unemp_status_t3_jun_lead) | unemp_status_t3_jun_lead == 0 ~ 19, 
                                                                  is.na(unemp_status_t3_jul_lead) | unemp_status_t3_jul_lead == 0 ~ 20, 
                                                                  is.na(unemp_status_t3_aug_lead) | unemp_status_t3_aug_lead == 0 ~ 21, 
                                                                  is.na(unemp_status_t3_sep_lead) | unemp_status_t3_sep_lead == 0 ~ 22, 
                                                                  is.na(unemp_status_t3_oct_lead) | unemp_status_t3_oct_lead == 0 ~ 23, 
                                                                  is.na(unemp_status_t3_nov_lead) | unemp_status_t3_nov_lead == 0 ~ 24,
                                                                  TRUE ~ 25)) -> data_first_unemp_m23

#Re-merging datasets
data <- rbind(data_pre2001_or_emp, 
              data_first_unemp_m1, data_first_unemp_m2, data_first_unemp_m3, data_first_unemp_m4,
              data_first_unemp_m5, data_first_unemp_m6, data_first_unemp_m7, data_first_unemp_m8,
              data_first_unemp_m9, data_first_unemp_m10, data_first_unemp_m11, data_first_unemp_m12,
              data_first_unemp_m13, data_first_unemp_m14, data_first_unemp_m15, data_first_unemp_m16,
              data_first_unemp_m17, data_first_unemp_m18, data_first_unemp_m19, data_first_unemp_m20,
              data_first_unemp_m21, data_first_unemp_m22, data_first_unemp_m23)
rm(data_pre2001_or_emp, 
   data_first_unemp_m1, data_first_unemp_m2, data_first_unemp_m3, data_first_unemp_m4,
   data_first_unemp_m5, data_first_unemp_m6, data_first_unemp_m7, data_first_unemp_m8,
   data_first_unemp_m9, data_first_unemp_m10, data_first_unemp_m11, data_first_unemp_m12,
   data_first_unemp_m13, data_first_unemp_m14, data_first_unemp_m15, data_first_unemp_m16,
   data_first_unemp_m17, data_first_unemp_m18, data_first_unemp_m19, data_first_unemp_m20,
   data_first_unemp_m21, data_first_unemp_m22, data_first_unemp_m23)

#Creating 'mo_unemp_post_int_na' variable
data %>% mutate(mo_unemp_post_int_na = case_when(is.na(first_mo_unemp_post_int) | year < 2001 ~ NA_integer_,
                                                 is.na(mo_unemp_post_int) ~ 1,
                                                 first_mo_unemp_post_int == 1 & mo_unemp_post_int_min == 47 ~ 1,
                                                 first_mo_unemp_post_int == 2 & mo_unemp_post_int_min == 46 ~ 1,
                                                 first_mo_unemp_post_int == 3 & mo_unemp_post_int_min == 45 ~ 1,
                                                 first_mo_unemp_post_int == 4 & mo_unemp_post_int_min == 44 ~ 1,
                                                 first_mo_unemp_post_int == 5 & mo_unemp_post_int_min == 43 ~ 1,
                                                 first_mo_unemp_post_int == 6 & mo_unemp_post_int_min == 42 ~ 1,
                                                 first_mo_unemp_post_int == 7 & mo_unemp_post_int_min == 41 ~ 1,
                                                 first_mo_unemp_post_int == 8 & mo_unemp_post_int_min == 40 ~ 1,
                                                 first_mo_unemp_post_int == 9 & mo_unemp_post_int_min == 39 ~ 1,
                                                 first_mo_unemp_post_int == 10 & mo_unemp_post_int_min == 38 ~ 1,
                                                 first_mo_unemp_post_int == 11 & mo_unemp_post_int_min == 37 ~ 1,
                                                 first_mo_unemp_post_int == 12 & mo_unemp_post_int_min == 36 ~ 1,
                                                 first_mo_unemp_post_int == 13 & mo_unemp_post_int_min == 35 ~ 1,
                                                 first_mo_unemp_post_int == 14 & mo_unemp_post_int_min == 34 ~ 1,
                                                 first_mo_unemp_post_int == 15 & mo_unemp_post_int_min == 33 ~ 1,
                                                 first_mo_unemp_post_int == 16 & mo_unemp_post_int_min == 32 ~ 1,
                                                 first_mo_unemp_post_int == 17 & mo_unemp_post_int_min == 31 ~ 1,
                                                 first_mo_unemp_post_int == 18 & mo_unemp_post_int_min == 30 ~ 1,
                                                 first_mo_unemp_post_int == 19 & mo_unemp_post_int_min == 29 ~ 1,
                                                 first_mo_unemp_post_int == 20 & mo_unemp_post_int_min == 28 ~ 1,
                                                 first_mo_unemp_post_int == 21 & mo_unemp_post_int_min == 27 ~ 1,
                                                 first_mo_unemp_post_int == 22 & mo_unemp_post_int_min == 26 ~ 1,
                                                 first_mo_unemp_post_int == 23 & mo_unemp_post_int_min == 25 ~ 1,
                                                 mo_unemp_post_int == mo_unemp_post_int_min ~ 0,
                                                 TRUE ~ 99)) -> data #Noone should be left uncoded, so TRUE ~ 99 indicates an issue

#Creating mo_unemp_until_emp' and 'mo_unemp_until_emp_na' variables
#NOTE - For 'mo_unemp_until_emp', this equals either 1) the number of consecutive months unemployed before re-employment, or
#                                                    2) the minimum consecutive months known unemployed before re-employment
#         - We need to use 'mo_unemp_until_emp_na' to determine which it is.
#         - For all possible scenarios based on 1) whether already unemployed (or known), 2) whether first month of unemployment was known, and 3) whether pre- and post-interview unemployment were known, as appropriate 
#     - For 'mo_unemp_until_emp_na', assign = NA if pre-2001, employed for 12 months post-int, or if this is NA
#                                           = 'actual' if all relevant unemployment information is known
#                                           = 'cens_left' if information is only left-censored (where unemployed pre-interview)
#                                           = 'cens_right' if information is only right-censored (where pre-interview info is known if already unemployed, or if not newly unemployed), or 
#                                           = 'cens_both' if left and right-censored (missing or censored both pre- and post-interview)
data %>% mutate(mo_unemp_until_emp = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | year < 2001 ~ NA_integer_,
                                               is.na(mo_unemp_pre_int_min) & is.na(mo_unemp_post_int_min) ~ 1, #To be counted as unemployed, they need to have been unemployed for at least 1 month
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int > interview_month) & !is.na(mo_unemp_post_int_min) ~ mo_unemp_post_int_min,
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int > interview_month) & is.na(mo_unemp_post_int_min) ~ 1,
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & !is.na(already_unemp_pre_int) & already_unemp_pre_int == 0 & !is.na(mo_unemp_post_int_min) ~ mo_unemp_post_int_min,
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & !is.na(already_unemp_pre_int) & already_unemp_pre_int == 0 & is.na(mo_unemp_post_int_min) ~ 1,
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) ~ (1 + mo_unemp_post_int_min),
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & is.na(mo_unemp_post_int_min) ~ (1 + mo_unemp_pre_int_min),
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) ~ (mo_unemp_pre_int_min + mo_unemp_post_int_min),
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & is.na(already_unemp_pre_int) & is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) ~ mo_unemp_post_int_min,
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & is.na(already_unemp_pre_int) & !is.na(mo_unemp_pre_int_min) & is.na(mo_unemp_post_int_min) ~ 1, 
                                               !is.na(first_mo_unemp_post_int) & (first_mo_unemp_post_int == interview_month) & is.na(already_unemp_pre_int) & !is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) ~ mo_unemp_post_int_min,
                                               is.na(first_mo_unemp_post_int) & !is.na(already_unemp_pre_int) & already_unemp_pre_int == 0 ~ 1, 
                                               is.na(first_mo_unemp_post_int) & !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 ~ (1 + mo_unemp_pre_int_min),
                                               is.na(first_mo_unemp_post_int) & is.na(already_unemp_pre_int) ~ 1, 
                                               TRUE ~ 999), #Noone should be left uncoded, so TRUE ~ 99 indicates an issue
                mo_unemp_until_emp_na = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | year < 2001 ~ NA_character_,
                                                  is.na(already_unemp_pre_int) & is.na(mo_unemp_post_int_min) ~ "cens_both",
                                                  is.na(already_unemp_pre_int) & !is.na(mo_unemp_post_int_min) & mo_unemp_post_int_na == 0 ~ "cens_left",
                                                  is.na(already_unemp_pre_int) & !is.na(mo_unemp_post_int_min) & mo_unemp_post_int_na == 1 ~ "cens_both",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 0 & is.na(mo_unemp_post_int_min) ~ "cens_right",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 0 & !is.na(mo_unemp_post_int_min) & mo_unemp_post_int_na == 0 ~ "actual",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 0 & !is.na(mo_unemp_post_int_min) & mo_unemp_post_int_na == 1 ~ "cens_right",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & is.na(mo_unemp_pre_int_min) & is.na(mo_unemp_post_int_min) ~ "cens_both",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) & mo_unemp_post_int_na == 0 ~ "cens_left",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) & mo_unemp_post_int_na == 1 ~ "cens_both",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & is.na(mo_unemp_post_int_min) & mo_unemp_pre_int_na == 0 ~ "cens_right",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & is.na(mo_unemp_post_int_min) & mo_unemp_pre_int_na == 1 ~ "cens_both",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) & mo_unemp_pre_int_na == 0 & mo_unemp_post_int_na == 0 ~ "actual",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) & mo_unemp_pre_int_na == 0 & mo_unemp_post_int_na == 1 ~ "cens_right",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) & mo_unemp_pre_int_na == 1 & mo_unemp_post_int_na == 0 ~ "cens_left",
                                                  !is.na(already_unemp_pre_int) & already_unemp_pre_int == 1 & !is.na(mo_unemp_pre_int_min) & !is.na(mo_unemp_post_int_min) & mo_unemp_pre_int_na == 1 & mo_unemp_post_int_na == 1 ~ "cens_both",
                                                  TRUE ~ "mistake")) -> data

#### Step 7 - Create ui_receipt_MONTH_lead variables ####
data <- data %>% arrange(unique_id, year)
library(data.table)
data <- as.data.table(data)
data[, ':=' (ui_receipt_jan_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_jan)), NA_integer_, lead(ui_receipt_jan)),
             ui_receipt_feb_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_feb)), NA_integer_, lead(ui_receipt_feb)),
             ui_receipt_mar_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_mar)), NA_integer_, lead(ui_receipt_mar)),
             ui_receipt_apr_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_apr)), NA_integer_, lead(ui_receipt_apr)),
             ui_receipt_may_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_may)), NA_integer_, lead(ui_receipt_may)),
             ui_receipt_jun_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_jun)), NA_integer_, lead(ui_receipt_jun)),
             ui_receipt_jul_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_jul)), NA_integer_, lead(ui_receipt_jul)),
             ui_receipt_aug_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_aug)), NA_integer_, lead(ui_receipt_aug)),
             ui_receipt_sep_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_sep)), NA_integer_, lead(ui_receipt_sep)),
             ui_receipt_oct_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_oct)), NA_integer_, lead(ui_receipt_oct)),
             ui_receipt_nov_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_nov)), NA_integer_, lead(ui_receipt_nov)),
             ui_receipt_dec_lead = fifelse(lead(year) != (year+2) | is.na(lead(ui_receipt_dec)), NA_integer_, lead(ui_receipt_dec))), by = .(unique_id)]
data <- as.data.frame(data)
detach("package:data.table", unload = TRUE)

#### Step 8 - Creating indicators for complete 6-month UI information and 6-month post-interview unempoyment UI receipt ####
#           - Creating 1) 'all_ui_receipt_info_available' indicator for complete UI information in 6 months post-unemployment (or interview month if unemployed in int month), and
#                      2) 'ui_receipt_6mo_post_unemp' indicator for UI received within 6 months post-unemployment
#NOTE - For 'all_ui_receipt_info_available', this must be =0 for all those first unemployed in months < 12 (Dec T0) or >18 (Jun T1)
#     - For 'ui_receipt_6mo_post_unemp', if 'all_ui_receipt_info_available' is 0 but at least 1 month's UI receipt indicator is 1 assign 1, otherwise NA
#         - Also assign NA if 1) not unemployed in the year post-interview, 2) unknown unemployment status in the year post-interview, or 3) interviewed in an atypical year

#Creating first unemployment month-specific datasets to work out 6 months post-unemployment UI information
data_checkUIavailability <- data[!is.na(data$first_mo_unemp_post_int), ]
data_checkUIavailability_int7 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 7, ]
data_checkUIavailability_int8 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 8, ]
data_checkUIavailability_int9 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 9, ]
data_checkUIavailability_int10 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 10, ]
data_checkUIavailability_int11 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 11, ]
data_checkUIavailability_int12 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 12, ]
data_checkUIavailability_int13 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 13, ]
data_checkUIavailability_int14 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 14, ]
data_checkUIavailability_int15 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 15, ]
data_checkUIavailability_int16 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 16, ]
data_checkUIavailability_int17 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 17, ]
data_checkUIavailability_int18 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 18, ]
data_checkUIavailability_int19 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 19, ]
data_checkUIavailability_int20 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 20, ]
data_checkUIavailability_int21 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 21, ]
data_checkUIavailability_int22 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 22, ]
data_checkUIavailability_int23 <- data_checkUIavailability[data_checkUIavailability$first_mo_unemp_post_int == 23, ]
rm(data_checkUIavailability)

#For each month-specific dataset, creating indicators for 1) complete UI information over the 6-month period and 2) UI receipt during this period
#NOTE - For ui_receipt_6mo_post_unemp, where some relevant UI month information is missing but what exists indicates UI receipt, assign 1, otherwise NA
data_checkUIavailability_int7$all_ui_receipt_info_available <- 0
data_checkUIavailability_int7$ui_receipt_6mo_post_unemp <- ifelse(!is.na(data_checkUIavailability_int7$ui_receipt_jan_lead) & data_checkUIavailability_int7$ui_receipt_jan_lead > 0, 1, NA_integer_)
data_checkUIavailability_int8$all_ui_receipt_info_available <- 0
data_checkUIavailability_int8$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int8[,c(135:136)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int9$all_ui_receipt_info_available <- 0
data_checkUIavailability_int9$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int9[,c(135:137)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int10$all_ui_receipt_info_available <- 0
data_checkUIavailability_int10$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int10[,c(135:138)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int11$all_ui_receipt_info_available <- 0
data_checkUIavailability_int11$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int11[,c(135:139)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int12$all_ui_receipt_info_available <- ifelse(is.na(data_checkUIavailability_int12$ui_receipt_jan_lead) | is.na(data_checkUIavailability_int12$ui_receipt_feb_lead) | 
                                                                         is.na(data_checkUIavailability_int12$ui_receipt_mar_lead) | is.na(data_checkUIavailability_int12$ui_receipt_apr_lead) | 
                                                                         is.na(data_checkUIavailability_int12$ui_receipt_may_lead) | is.na(data_checkUIavailability_int12$ui_receipt_jun_lead), 0, 1)
data_checkUIavailability_int12$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int12[,c(135:140)], na.rm=T) > 0, 1,
                                                                   ifelse(data_checkUIavailability_int12$all_ui_receipt_info_available == 0, NA_integer_, 0))
data_checkUIavailability_int13$all_ui_receipt_info_available <- ifelse(is.na(data_checkUIavailability_int13$ui_receipt_feb_lead) | is.na(data_checkUIavailability_int13$ui_receipt_mar_lead) | 
                                                                         is.na(data_checkUIavailability_int13$ui_receipt_apr_lead) | is.na(data_checkUIavailability_int13$ui_receipt_may_lead) | 
                                                                         is.na(data_checkUIavailability_int13$ui_receipt_jun_lead) | is.na(data_checkUIavailability_int13$ui_receipt_jul_lead), 0, 1)
data_checkUIavailability_int13$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int13[,c(136:141)], na.rm=T) > 0, 1,
                                                                   ifelse(data_checkUIavailability_int13$all_ui_receipt_info_available == 0, NA_integer_, 0))
data_checkUIavailability_int14$all_ui_receipt_info_available <- ifelse(is.na(data_checkUIavailability_int14$ui_receipt_mar_lead) | is.na(data_checkUIavailability_int14$ui_receipt_apr_lead) | 
                                                                         is.na(data_checkUIavailability_int14$ui_receipt_may_lead) | is.na(data_checkUIavailability_int14$ui_receipt_jun_lead) | 
                                                                         is.na(data_checkUIavailability_int14$ui_receipt_jul_lead) | is.na(data_checkUIavailability_int14$ui_receipt_aug_lead), 0, 1)
data_checkUIavailability_int14$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int14[,c(137:142)], na.rm=T) > 0, 1,
                                                                   ifelse(data_checkUIavailability_int14$all_ui_receipt_info_available == 0, NA_integer_, 0))
data_checkUIavailability_int15$all_ui_receipt_info_available <- ifelse(is.na(data_checkUIavailability_int15$ui_receipt_apr_lead) | is.na(data_checkUIavailability_int15$ui_receipt_may_lead) | 
                                                                         is.na(data_checkUIavailability_int15$ui_receipt_jun_lead) | is.na(data_checkUIavailability_int15$ui_receipt_jul_lead) | 
                                                                         is.na(data_checkUIavailability_int15$ui_receipt_aug_lead) | is.na(data_checkUIavailability_int15$ui_receipt_sep_lead), 0, 1)
data_checkUIavailability_int15$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int15[,c(138:143)], na.rm=T) > 0, 1,
                                                                   ifelse(data_checkUIavailability_int15$all_ui_receipt_info_available == 0, NA_integer_, 0))
data_checkUIavailability_int16$all_ui_receipt_info_available <- ifelse(is.na(data_checkUIavailability_int16$ui_receipt_may_lead) | is.na(data_checkUIavailability_int16$ui_receipt_jun_lead) | 
                                                                         is.na(data_checkUIavailability_int16$ui_receipt_jul_lead) | is.na(data_checkUIavailability_int16$ui_receipt_aug_lead) | 
                                                                         is.na(data_checkUIavailability_int16$ui_receipt_sep_lead) | is.na(data_checkUIavailability_int16$ui_receipt_oct_lead), 0, 1)
data_checkUIavailability_int16$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int16[,c(139:144)], na.rm=T) > 0, 1,
                                                                   ifelse(data_checkUIavailability_int16$all_ui_receipt_info_available == 0, NA_integer_, 0))
data_checkUIavailability_int17$all_ui_receipt_info_available <- ifelse(is.na(data_checkUIavailability_int17$ui_receipt_jun_lead) | is.na(data_checkUIavailability_int17$ui_receipt_jul_lead) | 
                                                                         is.na(data_checkUIavailability_int17$ui_receipt_aug_lead) | is.na(data_checkUIavailability_int17$ui_receipt_sep_lead) | 
                                                                         is.na(data_checkUIavailability_int17$ui_receipt_oct_lead) | is.na(data_checkUIavailability_int17$ui_receipt_nov_lead), 0, 1)
data_checkUIavailability_int17$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int17[,c(140:145)], na.rm=T) > 0, 1,
                                                                   ifelse(data_checkUIavailability_int17$all_ui_receipt_info_available == 0, NA_integer_, 0))
data_checkUIavailability_int18$all_ui_receipt_info_available <- ifelse(is.na(data_checkUIavailability_int18$ui_receipt_jul_lead) | is.na(data_checkUIavailability_int18$ui_receipt_aug_lead) | 
                                                                         is.na(data_checkUIavailability_int18$ui_receipt_sep_lead) | is.na(data_checkUIavailability_int18$ui_receipt_oct_lead) | 
                                                                         is.na(data_checkUIavailability_int18$ui_receipt_nov_lead) | is.na(data_checkUIavailability_int18$ui_receipt_dec_lead), 0, 1)
data_checkUIavailability_int18$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int18[,c(141:146)], na.rm=T) > 0, 1,
                                                                   ifelse(data_checkUIavailability_int18$all_ui_receipt_info_available == 0, NA_integer_, 0))
data_checkUIavailability_int19$all_ui_receipt_info_available <- 0
data_checkUIavailability_int19$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int19[,c(142:146)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int20$all_ui_receipt_info_available <- 0
data_checkUIavailability_int20$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int20[,c(143:146)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int21$all_ui_receipt_info_available <- 0
data_checkUIavailability_int21$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int21[,c(144:146)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int22$all_ui_receipt_info_available <- 0
data_checkUIavailability_int22$ui_receipt_6mo_post_unemp <- ifelse(rowSums(data_checkUIavailability_int22[,c(145:146)], na.rm=T) > 0, 1, NA_integer_)
data_checkUIavailability_int23$all_ui_receipt_info_available <- 0
data_checkUIavailability_int23$ui_receipt_6mo_post_unemp <- ifelse(!is.na(data_checkUIavailability_int23$ui_receipt_dec_lead) & data_checkUIavailability_int23$ui_receipt_dec_lead > 0, 1, NA_integer_)

#Re-merging month-specific datasets
data_ui_available <- rbind(data_checkUIavailability_int7, data_checkUIavailability_int8, data_checkUIavailability_int9,
                           data_checkUIavailability_int10, data_checkUIavailability_int11, data_checkUIavailability_int12, 
                           data_checkUIavailability_int13, data_checkUIavailability_int14, data_checkUIavailability_int15,
                           data_checkUIavailability_int16, data_checkUIavailability_int17, data_checkUIavailability_int18,
                           data_checkUIavailability_int19, data_checkUIavailability_int20, data_checkUIavailability_int21,
                           data_checkUIavailability_int22, data_checkUIavailability_int23)

#Adding these two indicators information back into the larger dataset
data_ui_available <- data_ui_available[, c("unique_id", "year", "all_ui_receipt_info_available", "ui_receipt_6mo_post_unemp")]
data <- merge(data, data_ui_available, by = c("unique_id", "year"), all = T)

#Assigning 0 (=No) for all those unemployed in the year post-interview whose month of first unemployment 
#precludes complete UI information, given UI information was only collected every second year
data$all_ui_receipt_info_available[!is.na(data$first_mo_unemp_post_int) & (data$first_mo_unemp_post_int < 7 | data$first_mo_unemp_post_int > 23)] <- 0 
rm(data_checkUIavailability_int7, data_checkUIavailability_int8, data_checkUIavailability_int9,
   data_checkUIavailability_int10, data_checkUIavailability_int11, data_checkUIavailability_int12, 
   data_checkUIavailability_int13, data_checkUIavailability_int14, data_checkUIavailability_int15,
   data_checkUIavailability_int16, data_checkUIavailability_int17, data_checkUIavailability_int18,
   data_checkUIavailability_int19, data_checkUIavailability_int20, data_checkUIavailability_int21,
   data_checkUIavailability_int22, data_checkUIavailability_int23, data_ui_available)

#### Step 9 - Simple imputation approach for UI receipt in 6 months post-interview (or unemployment, whichever is latest) ####
#NOTES - Most states determine eligibilty based on hours worked or wages earned during a 'base period' which is typically 
#        the first 4 of the last 5 complete calendar quarters. While we could work out the relevant base period for most 
#        individuals, we don't have monthly hours worked or income earned information, so cannot work out the relevant quantity
#        over the base period. Instead, we use the annual hours worked or annual labor income over the prior year, given 
#        this will include a significant portion of the base period per individual, or at least be highly correlated with 
#        base period income for that person.
#      - Related, the vast majority of individuals apply for UI within the first month of becoming unemployed. For this reason, 
#        if respondents were unemployed for at least the preceding 7 months (1 month wait + 6 months of possible receipt) before
#        their interview, we assume they would not receive any UI during the 12 months post-interview.
#      - We only impute *not* receiving UI based on these criteria. We do not impute receipt of UI.
#      - We assume where states are missing minimum base period weeks worked or average hours worked, this is not an eligibility requirement

#Creating indicator for 'base_period_year' (=year + 1, to match the observation's year reporting prior-year income in year -1)
#E.g. if the UI dataset obs. for TN in year 2002 was $3000, base_year_period would = 2003, 
#     as the main dataset observation containing 2002 annual income would be from the year 2003
ui_state_elig_data$base_period_year <- ui_state_elig_data$year + 1

#Creating average annual taxable_wage_base variable (~28% of state-year obs changed this value for Q1/2 vs. Q3/4, mostly in 2001-02, 2010-11, and 2017)
ui_state_elig_data %>% 
  select(state, base_period_year, base_period_wage_min, base_period_hours_min, base_period_weeks_worked_min) %>%
  filter(base_period_year >= 2001) -> ui_state_elig_data

#Merging in State information on base period year minimum eligible taxable income
#NOTE - This data only extends to 2017, which is fine for our purposes where we look at 2001-2017 data
data <- merge(data, ui_state_elig_data, by.x = c("year", "state"), by.y = c("base_period_year", "state"), all.x = T)
rm(ui_state_elig_data)

#Merging in total_labor_income from demog_data
demog_data <- qread(paste0(directory, "psid_dataset_110321.qs")) #Demographic data
demog_data <- demog_data[,c("unique_id", "year", "total_labor_income")]
data <- merge(data, demog_data, by = c("unique_id", "year"), all.x = T)
rm(demog_data)

#Re-coding ui_receipt_6mo_post_unemp as '0' where: 
# - 1) Unemployed pre-interview for more than 7 months, or
# - 2) Reported prior-year labor income (which we use as a proxy for base period income) is less than the average minimum taxable wage base for that year 
data %>% group_by(state, year) %>%
  mutate(ui_receipt_6mo_post_unemp = case_when(is.na(unemp_status_lead) | unemp_status_lead == 0 | !is.na(ui_receipt_6mo_post_unemp) ~ ui_receipt_6mo_post_unemp, #Keeping as already assigned
                                               unemp_status_lead == 1 & mo_unemp_pre_int_min > 7 ~ 0, #Max. likely period of UI receipt 
                                               unemp_status_lead == 1 & !is.na(total_labor_income) & total_labor_income < base_period_wage_min ~ 0,
                                               unemp_status_lead == 1 & !is.na(annual_hours_worked) & !is.na(base_period_hours_min) & annual_hours_worked < base_period_hours_min ~ 0,
                                               unemp_status_lead == 1 & !is.na(unemp_duration_months) & !is.na(base_period_weeks_worked_min) & (52 - (unemp_duration_months)*(4+1/3)) < base_period_weeks_worked_min ~ 0,
                                               TRUE ~ ui_receipt_6mo_post_unemp)) %>%
  ungroup() -> data

#### Step 10 - Creating 'state_ui_date' variable and merging in 'max_week_duration'  and 'max_weekly_ben_dep' ####
#             for state-specific max UI duration and weekly benefits w/ dependents in the year half 1 month prior to their interview
#Carry forward nearest maximum weekly benefit amount for state observations where maximum weekly benefit amount information is missing
#Excluding 2016 record with state=="ok" as a duplicate of the 2016 record with state=="OK"
ui_state_data <- ui_state_data[ui_state_data$state != "ok", ] 

#Creating separate month and year vars based on 'date' variable in UI state dataset
ui_state_data$ui_year <- as.numeric(substr(ui_state_data$date, start = 5, stop = 8))
ui_state_data$ui_month <- as.numeric(substr(ui_state_data$date, start = 1, stop = 1))

ui_state_data %>% arrange(state, ui_year) %>% group_by(state) %>%
  mutate(max_weekly_ben_dep = case_when(is.na(max_weekly_ben_dep) & !is.na(lag(max_weekly_ben_dep, 1)) ~ lag(max_weekly_ben_dep, 1),
                                        is.na(max_weekly_ben_dep) & !is.na(lag(max_weekly_ben_dep, 2)) ~ lag(max_weekly_ben_dep, 2),
                                        is.na(max_weekly_ben_dep) & !is.na(lag(max_weekly_ben_dep, 3)) ~ lag(max_weekly_ben_dep, 3),
                                        TRUE ~ max_weekly_ben_dep)) %>% ungroup() -> ui_state_data 

#Creating matching variable across datasets for UI year half
#NOTE - Where respondents are known unemployed but their first month unemployed is unknown, 
#       we assume the relevant state UI data year is the interview year, unless interviewed in January, in which case we assume the prior yaer
#     - Where respondents are known not unemployed in the year post-interview, we assign:
#         - The year they were interviewed if interviewed in February-December, or the prior year if interviewed in January
data %>% mutate(ui_year = case_when(year < 2001 | interview_year != year ~ year,
                                    is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month == 1 ~ (year-1),
                                    is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month > 1 ~ year,
                                    is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month == 1 ~ (year-1),
                                    is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month > 1 ~ year,
                                    first_mo_unemp_post_int == 1 ~ (year-1),
                                    first_mo_unemp_post_int > 1 & first_mo_unemp_post_int <= 13 ~ year,
                                    TRUE ~ (year+1)), 
                ui_month = case_when(year < 2001 ~ 1,
                                     is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month <=7 ~ 1,
                                     is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month >7 ~ 7,
                                     is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month <= 7 ~ 1,
                                     is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month > 7 ~ 7,
                                     first_mo_unemp_post_int == 1 ~ 7,
                                     first_mo_unemp_post_int >1 & first_mo_unemp_post_int <= 7 ~ 1,
                                     first_mo_unemp_post_int >7 & first_mo_unemp_post_int <=13 ~ 7,
                                     first_mo_unemp_post_int >13 & first_mo_unemp_post_int <=19 ~ 1,
                                     TRUE ~ 7)) -> data

#Merging in UI 'max_week_duration' variable
data <- merge(data, ui_state_data[,c("state", "max_weekly_ben_dep", "max_week_duration", "ui_year", "ui_month")], 
              by = c("state", "ui_year", "ui_month"), all.x = T)
rm(ui_state_data)

#### Step 11 - Creating 'ui_annual_amount_lead', 'ui_months_lead', and 'ui_monthly_amount_lead' to reflect average $ received per month ####
#For speed in this step, I convert to a data.table
#NOTE - 1 respondent reports a UI amount without reporting receiving UI in any month. This person has ui_monthly_amount_lead set to NA
#     - 166 respondents report UI without complete information on months received. I set months to the max known for these ids, or NA if max = 0
data %>% arrange(unique_id, year) -> data
library(data.table)
data <- as.data.table(data)
data[, ui_annual_amount_lead := fifelse((lead(year,1) != (year+2)) | is.na(lead(ui_annual_amount,1)), NA_integer_,
                                        lead(ui_annual_amount, 1)), by = .(unique_id)]
data[, ui_months_lead := fifelse(is.na(ui_receipt_jan_lead) | is.na(ui_receipt_feb_lead) | is.na(ui_receipt_mar_lead) | 
                                   is.na(ui_receipt_apr_lead) | is.na(ui_receipt_may_lead) | is.na(ui_receipt_jun_lead) | 
                                   is.na(ui_receipt_jul_lead) | is.na(ui_receipt_aug_lead) | is.na(ui_receipt_sep_lead) | 
                                   is.na(ui_receipt_oct_lead) | is.na(ui_receipt_nov_lead) | is.na(ui_receipt_dec_lead), NA_integer_,
                                 (ui_receipt_jan_lead + ui_receipt_feb_lead + ui_receipt_mar_lead +
                                    ui_receipt_apr_lead + ui_receipt_may_lead + ui_receipt_jun_lead +
                                    ui_receipt_jul_lead + ui_receipt_aug_lead + ui_receipt_sep_lead +
                                    ui_receipt_oct_lead + ui_receipt_nov_lead + ui_receipt_dec_lead)), by = .(unique_id)]
data[, ui_monthly_amount_lead := fifelse(is.na(ui_annual_amount_lead) | is.na(ui_months_lead), NA_real_,
                                         (ui_annual_amount_lead / ui_months_lead)), by = .(unique_id)]
data <- as.data.frame(data)
detach("package:data.table", unload = TRUE)

#Correcting UI months to = max known UI months where 1) annual UI is >0 and 2) there is incomplete UI information
data$ui_months_lead[data$ui_annual_amount_lead > 0 & 
                      !is.na(data$ui_annual_amount_lead) & is.na(data$ui_months_lead)] <- rowSums(data[data$ui_annual_amount_lead > 0 & 
                                                                                                         !is.na(data$ui_annual_amount_lead) & is.na(data$ui_months_lead), c(103:114)], na.rm = T)

#Setting to NA those who report 0 months of UI receipt
data$ui_monthly_amount_lead <- ifelse(data$ui_annual_amount_lead > 0 & !is.na(data$ui_annual_amount_lead) & 
                                        data$ui_months_lead == 0 & !is.na(data$ui_months_lead), NA_integer_,
                                      (data$ui_annual_amount_lead / data$ui_months_lead))

#### Step 12 - Merging in State GSP and Unemployment Rate Data ####
#NOTE - Per respondent observation, we take the State GSP and unemployment rate from the 3-month period pre-unemployment (or interview month, if employed)
#Creating 'state' variable with 2-letter abbreviation based on FIPS code for state_gsp_data and state_ur_data for merging
state_gsp_data %>% mutate(state = case_when(fips == 01 ~ "AL", fips == 02 ~ "AK", fips == 04 ~ "AZ", fips == 05 ~ "AR", 
                                            fips == 06 ~ "CA", fips == 08 ~ "CO", fips == 09 ~ "CT", fips == 10 ~ "DE",
                                            fips == 11 ~ "DC", fips == 12 ~ "FL", fips == 13 ~ "GA", fips == 15 ~ "HI",
                                            fips == 16 ~ "ID", fips == 17 ~ "IL", fips == 18 ~ "IN", fips == 19 ~ "IA",
                                            fips == 20 ~ "KS", fips == 21 ~ "KY", fips == 22 ~ "LA", fips == 23 ~ "ME",
                                            fips == 24 ~ "MD", fips == 25 ~ "MA", fips == 26 ~ "MI", fips == 27 ~ "MN",
                                            fips == 28 ~ "MS", fips == 29 ~ "MO", fips == 30 ~ "MT", fips == 31 ~ "NE",
                                            fips == 32 ~ "NV", fips == 33 ~ "NH", fips == 34 ~ "NJ", fips == 35 ~ "NM",
                                            fips == 36 ~ "NY", fips == 37 ~ "NC", fips == 38 ~ "ND", fips == 39 ~ "OH",
                                            fips == 40 ~ "OK", fips == 41 ~ "OR", fips == 42 ~ "PA", fips == 44 ~ "RI",
                                            fips == 45 ~ "SC", fips == 46 ~ "SD", fips == 47 ~ "TN", fips == 48 ~ "TX",
                                            fips == 49 ~ "UT", fips == 50 ~ "VT", fips == 51 ~ "VA", fips == 53 ~ "WA",
                                            fips == 54 ~ "WV", fips == 55 ~ "WI", fips == 56 ~ "WY", TRUE ~ "US Territory or Foreign Country")) -> state_gsp_data

state_ur_data %>% mutate(state = case_when(fips == 01 ~ "AL", fips == 02 ~ "AK", fips == 04 ~ "AZ", fips == 05 ~ "AR", 
                                           fips == 06 ~ "CA", fips == 08 ~ "CO", fips == 09 ~ "CT", fips == 10 ~ "DE",
                                           fips == 11 ~ "DC", fips == 12 ~ "FL", fips == 13 ~ "GA", fips == 15 ~ "HI",
                                           fips == 16 ~ "ID", fips == 17 ~ "IL", fips == 18 ~ "IN", fips == 19 ~ "IA",
                                           fips == 20 ~ "KS", fips == 21 ~ "KY", fips == 22 ~ "LA", fips == 23 ~ "ME",
                                           fips == 24 ~ "MD", fips == 25 ~ "MA", fips == 26 ~ "MI", fips == 27 ~ "MN",
                                           fips == 28 ~ "MS", fips == 29 ~ "MO", fips == 30 ~ "MT", fips == 31 ~ "NE",
                                           fips == 32 ~ "NV", fips == 33 ~ "NH", fips == 34 ~ "NJ", fips == 35 ~ "NM",
                                           fips == 36 ~ "NY", fips == 37 ~ "NC", fips == 38 ~ "ND", fips == 39 ~ "OH",
                                           fips == 40 ~ "OK", fips == 41 ~ "OR", fips == 42 ~ "PA", fips == 44 ~ "RI",
                                           fips == 45 ~ "SC", fips == 46 ~ "SD", fips == 47 ~ "TN", fips == 48 ~ "TX",
                                           fips == 49 ~ "UT", fips == 50 ~ "VT", fips == 51 ~ "VA", fips == 53 ~ "WA",
                                           fips == 54 ~ "WV", fips == 55 ~ "WI", fips == 56 ~ "WY", TRUE ~ "US Territory or Foreign Country")) -> state_ur_data

#Renaming 'year' as 'gsp_ur_year' in state_gsp and state_ur datasets
names(state_ur_data)[c(2,3)] <- c("gsp_ur_year", "gsp_ur_quarter")
names(state_gsp_data)[c(3,4)] <- c("gsp_ur_year", "gsp_ur_quarter")

#Creating matching variable across datasets for GSP/UR year quarter                                      
#NOTE - Where respondents are known unemployed but their first month unemployed is unknown, 
#       we assume the relevant state gsp_ur_year is the interview year, unless interviewed in January, in which case we assume the prior year
#     - Where respondents are known not unemployed in the year post-interview, we assign:
#         - The year they were interviewed if interviewed in February-December, or the prior year if interviewed in January
data %>% mutate(gsp_ur_year = case_when(year < 2001 | interview_year != year ~ year,
                                        is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month == 1 ~ (year-1),
                                        is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month > 1 ~ year,
                                        is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month == 1 ~ (year-1),
                                        is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month > 1 ~ year,
                                        first_mo_unemp_post_int == 1 ~ (year-1),
                                        first_mo_unemp_post_int > 1 & first_mo_unemp_post_int <= 13 ~ year,
                                        TRUE ~ (year+1)), 
                gsp_ur_quarter = case_when(year < 2001 ~ 1,
                                           is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month == 1 ~ 4,
                                           is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month <= 4 ~ 1,
                                           is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month <= 7 ~ 2,
                                           is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month <= 10 ~ 3,
                                           is.na(first_mo_unemp_post_int) & (unemp_status_lead == 0 | is.na(unemp_status_lead)) & interview_month <= 12 ~ 4,
                                           is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month == 1 ~ 4,
                                           is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month <= 4 ~ 1,
                                           is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month <= 7 ~ 2,
                                           is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month <= 10 ~ 3,
                                           is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month <= 12 ~ 4,
                                           first_mo_unemp_post_int == 1 ~ 4,
                                           first_mo_unemp_post_int %in% c(2:4,14:16) ~ 1,
                                           first_mo_unemp_post_int %in% c(5:7,17:19) ~ 2,
                                           first_mo_unemp_post_int %in% c(8:10,20:22) ~ 3,
                                           first_mo_unemp_post_int %in% c(11:13,23) ~ 4,
                                           TRUE ~ 99)) -> data #Noone should be left uncoded, so TRUE ~ 99 indicates an issue

#Merging in State Real GSP per Capita and State Unempoyment Rate variables
data <- merge(data, state_gsp_data[,c(6,3:5)], by = c("state", "gsp_ur_year", "gsp_ur_quarter"), all.x = T)
data <- merge(data, state_ur_data[,c(5,2:4)], by = c("state", "gsp_ur_year", "gsp_ur_quarter"), all.x = T)
data <- distinct(data)
rm(state_gsp_data, state_ur_data)

#### Step 13 - Restricting to the necessary set of variables for multiple imputation, determining wave eligibility, and modelling #### 
#Selecting relevant variables
id_vars <- names(data)[c(6,7,40,41)]
demog_vars <- names(data)[c(1,14,8:13,15:19,26,27)]
employment_vars <- names(data)[c(20:23,30:37,104:108,134,135,137,152,161)]
state_and_ui_vars <- names(data)[c(157,158,162,163)]
relevant_vars <- c(id_vars, demog_vars, employment_vars, state_and_ui_vars)
rm(id_vars, demog_vars, employment_vars, state_and_ui_vars)

#Restricting to the relevant set of variables
data <- data[, relevant_vars]
rm(relevant_vars)

#### Step 14 - Pre-MI Eligibility Exclusions ####
#Respondents eligibility is in-part determined by:
# - Year (2001-2017)
# - Age (18-64 years old)
# - Employment status (always either working, temp. layoff, or unemployed)
# - Self-employment status (respondents cannot be self-employed during follow-up)
# - Known gender, race, ethnicity, nativity, and childhood SES
#NOTE - Here we restrict to the maximum eligible set ob observations per respondent before MI

#Filling in race, ethnicity, gender, etc. via carry-forward/backward beforehand
data %>% 
  arrange(unique_id, year) %>% 
  group_by(unique_id) %>%
  mutate(gender = na.locf(gender, na.rm = F),
         gender = na.locf(gender, na.rm = F, fromLast = T),
         race = na.locf(race, na.rm = F),
         race = na.locf(race, na.rm = F, fromLast = T),
         ethnicity = na.locf(ethnicity, na.rm = F),
         ethnicity = na.locf(ethnicity, na.rm = F, fromLast = T),
         parents_poor = na.locf(parents_poor, na.rm = F),
         parents_poor = na.locf(parents_poor, na.rm = F, fromLast = T),
         nativity = na.locf(nativity, na.rm = F),
         nativity = na.locf(nativity, na.rm = F, fromLast = T)) %>% ungroup() -> data
data <- as.data.frame(data)

#Restrictions
data <- data[data$year %in% c(2001:2017), ] #Year
data <- data[!is.na(data$age) & data$age %in% c(18:64), ] #Age
data <- data[!is.na(data$employment_status) & data$employment_status %in% c("Working", "Temp. layoff", "Unemployed"), ] #Employment status
data <- data[data$employment_status == "Unemployed" | (!is.na(data$self_employed) & data$self_employed %in% c("Someone else", "Both someone else and self")), ] #Self-employed
data <- data[!is.na(data$gender) & !is.na(data$race) & !is.na(data$ethnicity) & !is.na(data$nativity) & !is.na(data$parents_poor), ] #Known key baseline variables
data <- as.data.frame(data)

#### Step 15 - Creating multiply-imputed datasets ####
#Calculating what percentages of cases are incomplete, to inform how many multiple imputed datasets are needed
#NOTE - As some variables are by-definition missing for employed or unemployed persons, we only check
#       whether observations are complete across the relevant range of variables per employment strata
#     - We no longer match using K-6, as it single-handedly accounted for a large proportion of our missing data
na_data <- data
na_data_unemp <- na_data[na_data$unemp_status_lead == 1 & !is.na(na_data$unemp_status_lead), ]
na_data_emp_or_na <- na_data[na_data$unemp_status_lead == 0 | is.na(na_data$unemp_status_lead), ]
unemp_complete <- NROW(na.omit(na_data_unemp[,c(1:16,18:20,32:40,42:45)]))
emp_or_na_complete <- NROW(na.omit(na_data_emp_or_na[,c(1:16,18:31,42:45)]))
percentage_na <- (100 - ((unemp_complete + emp_or_na_complete)/NROW(na_data))*100) #Indicates M > 31 is needed
rm(na_data, na_data_unemp, na_data_emp_or_na, emp_or_na_complete, unemp_complete)

#Recoding variables as will be used in analyses and to avoid MICE coding issues
#NOTE - For character/factor variables, we remove special characters which can cause errors with MICE
#     - For continuous variables where scaling can be an issue, we:
#         1) Truncate variables with large outlines (family_income and z_total_labor_income) at 2nd and 98th percentile
#         2) Recenter age so that 2001 = 0
#         2) Make a note of variable means and SDs to recode these back to their original format post-imputation
#             - Specifically only where the variable is used in analyses by itself (e.g. not EQ indicator variables)
#                - Age (Mu: 39.7, SD: 11.40922)
#                - Family income (Mu: 73755, SD: 54131.86)
#                - Months unemployed pre-interview minimum (Mu: 0.9022, SD: 4.220523)
#         3) Standardize variables
data %>% mutate(region = case_when(region == "US Territory or Foreign Country" ~ "South", TRUE ~ region),
                marital_status = str_replace_all(marital_status, "[ -,/]", ""),
                ethnicity = str_replace_all(ethnicity, "[ -,/]", ""),
                nativity = str_replace_all(nativity, "[ -,/]", ""),
                parents_poor = str_replace_all(parents_poor, "[ -,/]", ""),
                SRH = str_replace_all(SRH, "[ -,/]", ""),
                education = str_replace_all(education, "[ -,/]", ""),
                occupation_broad = str_replace_all(occupation_broad, "[ -,/]", ""),
                year = year - 2001,
                ethnicity = ifelse(ethnicity == "Non-Hispanic", "NonHispanic", ethnicity),
                education = ifelse(education == "<HS", "LessHS", education)) -> data 
data <- data %>% mutate(across(where(is.character), as.factor)) #Converts character variables to factors to avoid possible MI errors

data$family_income <- ifelse(data$family_income < quantile(data$family_income, c(0.02), na.rm=T),
                             quantile(data$family_income, c(0.02), na.rm=T), data$family_income)
data$family_income <- ifelse(data$family_income > quantile(data$family_income, c(0.98), na.rm=T),
                             quantile(data$family_income, c(0.98), na.rm=T), data$family_income)
data$z_total_labor_income <- ifelse(data$z_total_labor_income < quantile(data$z_total_labor_income, c(0.02), na.rm=T),
                                    quantile(data$z_total_labor_income, c(0.02), na.rm=T), data$z_total_labor_income)
data$z_total_labor_income <- ifelse(data$z_total_labor_income > quantile(data$z_total_labor_income, c(0.98), na.rm=T),
                                    quantile(data$z_total_labor_income, c(0.98), na.rm=T), data$z_total_labor_income)

data %>% mutate(age = (age - mean(age, na.rm=T))/sd(age, na.rm=T),
                family_income = (family_income - mean(family_income, na.rm=T))/sd(family_income, na.rm=T),
                annual_hours_worked = (annual_hours_worked - mean(annual_hours_worked, na.rm=T))/sd(annual_hours_worked, na.rm=T),
                unemp_duration_months = (unemp_duration_months - mean(unemp_duration_months, na.rm=T))/sd(unemp_duration_months, na.rm=T),
                mo_unemp_pre_int_min = (mo_unemp_pre_int_min - mean(mo_unemp_pre_int_min, na.rm=T))/sd(mo_unemp_pre_int_min, na.rm=T),
                gsp_per_cap = (gsp_per_cap - mean(gsp_per_cap, na.rm=T))/sd(gsp_per_cap, na.rm=T),
                unemp_rate_quart = (unemp_rate_quart - mean(unemp_rate_quart, na.rm=T))/sd(unemp_rate_quart, na.rm=T)) -> data

#Saving incomplete case dataset prior to multiple imputation
incomplete <- data
qsave(data, paste0(directory, "ui_preformatted_060823.qs"))

#Stratifying by employment status
#NOTE - The main reason we are conducting MI is to impute EQ. As we know unemployed persons
#       cannot have an EQ, we conduct our MI separately by current employment status.
#     - Likewise, to avoid complete conditional missingness of predictors by year, we do not:
#         1) Use as predictors variables which were only available over a given time-span attempt to
data_employed <- data[data$employment_status %in% c("Working", "Temp. layoff") & !is.na(data$employment_status), ]
data_employed$occupation_broad <- droplevels(data_employed$occupation_broad) #Could result in MICE non-convergence
data_unemployed <- data[data$employment_status == "Unemployed" & !is.na(data$employment_status), ]

# Constructing predictor matrix for those who were known to be currently employed
pred_matrix_employed <- make.predictorMatrix(data)
pred_matrix_employed[, c(1,3:5,17,19:21,23,33:35,37:39,41:45)] <- 0 #Excluding vars which are uninformative, irrelevant, or result in non-convergence (e.g. due to multi-collinearity) 
pred_matrix_employed[c(31:32), 1] <- -2  #Clustering EQ by ID

pred_matrix_unemployed <- make.predictorMatrix(data)
pred_matrix_unemployed[, c(1,3:5,17,19:35,37:39,41:45)] <- 0 #Excluding vars which are uninformative, irrelevant, or result in non-convergence (e.g. due to multi-collinearity) 

#Setting multiple imputation method per variable
mi_type_employed <- mice(data, maxit = 0)$method
mi_type_employed[c(31:32)] <- "2l.pmm"
mi_type_employed[c(1,3:5,17,19:21,23,33:35,37:39,41:45)] <- "" #Setting not to impute these variables

mi_type_unemployed <- mice(data, maxit = 0)$method
mi_type_unemployed[c(1,3:5,17,19:35,37:39,41:45)] <- "" #Setting not to impute these variables (including EQ)

#Performing imputation by employment status then saving. Will remerge employment strata to create complete datasets after
mi_datasets_employed <- mice(data_employed, method = mi_type_employed, pred = pred_matrix_employed, seed = 2023, maxit = 25, m = 35, printFlag = TRUE)
saveRDS(mi_datasets_employed, file = paste0(directory, "MI/ui_mids_employed_060823.rds"))

mi_datasets_unemployed <- mice(data_unemployed, method = mi_type_unemployed, pred = pred_matrix_unemployed, seed = 2023, maxit = 25, m = 35, printFlag = TRUE)
saveRDS(mi_datasets_unemployed, file = paste0(directory, "MI/ui_mids_unemployed_060823.rds"))

#### Step 16 - Re-merging employment-stratified MI datasets ####
#NOTE - After merging and before saving, we recode 'ui_receipt_6mo_post_unemp' as 0 if 'unemp_status_lead' is 0 (i.e. respondents are not unemployed in the 12 months post-interview)
for (mi in 1:35) {
  print(paste0("Current MI: ", mi))
  
  #Fixing UI receipt variable
  emp <- complete(mi_datasets_employed, mi)
  
  #Merging datasets
  unemp <- complete(mi_datasets_unemployed, mi)
  temp <- rbind(emp, unemp)
  temp %>% arrange(unique_id, year) -> temp
  
  #Saving datasets (will be over-written once post-formatting is done in Steps 17+)
  assign(paste0("mi_", mi), temp)
  qsave(temp, paste0(directory, "MI/ui_project_mi_", mi, "_060823.qs"))
  rm(emp, unemp, temp, mi)
}

#### Step 17 - Reformatting MI datasets ####
#NOTE - Here we: 1) impute directly known vars, 
#                2) change vars to minimum sensible coding (e.g. those imputed as 'not unemployed in 12 months' cannot receive UI in 6 months post-unemployment)
#                3) Excluding those with unknown state or in a given year residing in a foreign country or US territory
for (mi in 1:35) {
  print(paste0("Current MI: ", mi))
  eval(parse(text = paste0("temp <- mi_", mi)))
  
  #Resetting 12-month unemployment status, UI receipt, and months unemployed pre-interview to sensible minimums
  temp %>% mutate(unemp_status_lead = ifelse(is.na(unemp_status_lead) & employment_status == "Unemployed", 1, unemp_status_lead),
                  ui_receipt_6mo_post_unemp = ifelse(unemp_status_lead == 0, 0, ui_receipt_6mo_post_unemp),
                  mo_unemp_pre_int_min = ifelse(mo_unemp_pre_int_min < 0, 0, mo_unemp_pre_int_min)) -> temp
  
  #All those with missing Z-score GSP per capita are in states, GSP-relevant years, and GSP-relevant quarters with no non-missing observations
  #Therefore, we just remove the Z-score GSP per capita variable and re-merge in actual state-level GSP per capita information
  #Creating needed 'gsp_ur_year' and 'gsp_ur_quarter' variables for merging
  temp %>% mutate(gsp_ur_year = case_when(interview_year != (year+2001) ~ (year+2001),
                                          unemp_status_lead == 0 & interview_month == 1 ~ (year-1+2001),
                                          unemp_status_lead == 0 & interview_month > 1 ~ (year+2001),
                                          is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month == 1 ~ (year-1+2001),
                                          is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month > 1 ~ (year+2001),
                                          unemp_status_lead == 1 & first_mo_unemp_post_int == 1 ~ (year-1+2001),
                                          unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(2:13) ~ (year+2001),
                                          TRUE ~ (year+1+2001)),
                  gsp_ur_quarter = case_when(gsp_ur_year <= 2004 ~ 1, #We do not have monthly Real GSP per Capita information by state pre-2005 
                                             (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month == 1 ~ 4,
                                             (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month <= 4 ~ 1,
                                             (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month <= 7 ~ 2,
                                             (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month <= 10 ~ 3,
                                             (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month ~ 4,
                                             unemp_status_lead == 1 & first_mo_unemp_post_int == 1 ~ 4,
                                             unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(2:4,14:16) ~ 1,
                                             unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(5:7,17:19) ~ 2,
                                             unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(8:10,20:22) ~ 3,
                                             unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(11:13,23) ~ 4,
                                             TRUE ~ 99)) -> temp
  temp$gsp_per_cap <- NULL
  
  #Loading, formatting, merging gsp_per_cap with temp, then removing the state GSP dataset
  state_gsp_data <- read.csv(paste0(directory, "BEA_GSP_9719.csv")) #State GSP per Capita data
  state_gsp_data <- state_gsp_data[,c(2:5,8)]
  state_gsp_data %>% mutate(state = case_when(fips == 01 ~ "AL", fips == 02 ~ "AK", fips == 04 ~ "AZ", fips == 05 ~ "AR", 
                                              fips == 06 ~ "CA", fips == 08 ~ "CO", fips == 09 ~ "CT", fips == 10 ~ "DE",
                                              fips == 11 ~ "DC", fips == 12 ~ "FL", fips == 13 ~ "GA", fips == 15 ~ "HI",
                                              fips == 16 ~ "ID", fips == 17 ~ "IL", fips == 18 ~ "IN", fips == 19 ~ "IA",
                                              fips == 20 ~ "KS", fips == 21 ~ "KY", fips == 22 ~ "LA", fips == 23 ~ "ME",
                                              fips == 24 ~ "MD", fips == 25 ~ "MA", fips == 26 ~ "MI", fips == 27 ~ "MN",
                                              fips == 28 ~ "MS", fips == 29 ~ "MO", fips == 30 ~ "MT", fips == 31 ~ "NE",
                                              fips == 32 ~ "NV", fips == 33 ~ "NH", fips == 34 ~ "NJ", fips == 35 ~ "NM",
                                              fips == 36 ~ "NY", fips == 37 ~ "NC", fips == 38 ~ "ND", fips == 39 ~ "OH",
                                              fips == 40 ~ "OK", fips == 41 ~ "OR", fips == 42 ~ "PA", fips == 44 ~ "RI",
                                              fips == 45 ~ "SC", fips == 46 ~ "SD", fips == 47 ~ "TN", fips == 48 ~ "TX",
                                              fips == 49 ~ "UT", fips == 50 ~ "VT", fips == 51 ~ "VA", fips == 53 ~ "WA",
                                              fips == 54 ~ "WV", fips == 55 ~ "WI", fips == 56 ~ "WY", TRUE ~ "US Territory or Foreign Country")) -> state_gsp_data
  names(state_gsp_data)[c(3,4)] <- c("gsp_ur_year", "gsp_ur_quarter")
  temp <- merge(temp, state_gsp_data[,c(6,3:5)], by = c("state", "gsp_ur_year", "gsp_ur_quarter"), all.x = T)
  temp <- distinct(temp)
  rm(state_gsp_data)
  
  #Excluding empty rows or those residing in a given year in a US Territory or foreign country
  temp <- temp[!is.na(temp$unique_id), ]
  temp <- temp[temp$state != "US Territory or Foreign Country", ]
  temp <- temp[!is.na(temp$gsp_per_cap), ]
  
  #Renaming 'temp' as 'mi_X'
  assign(paste0("mi_", mi), temp)
  rm(temp, mi)
}

#### Step 18 - Creating wave-specific eligibility variable (eligible) and overall waves eligibile per participant indicator (include) based on known 1-lag and 2-lag waves 2 and 4 years prior ####
#NOTE - Waves are eligible if, for the same respondent, there are 2 prior waves of data which are from 2 and 4 years prior to the eligible wave
#     - The total set of eligible waves per respondent ('include') is the largest set of consecutive waves where 'eligible' = 1 as well as the 2 prior waves
## Creating 'eligible' indicator
for (mi in 1:35) {
  print(paste0("Current MI: ", mi))
  eval(parse(text = paste0("temp <- mi_", mi)))
  
  temp %>% arrange(unique_id, year) -> temp
  library(data.table)
  temp <- as.data.table(temp)
  temp[, eligible := fcase(is.na(lag(year,1)) | is.na(lag(year,2)), 0L,
                           lag(year,1) != (year-2) | lag(year,2) != (year-4), 0L,
                           default = 1L), by = .(unique_id)]
  temp <- as.data.frame(temp)
  detach("package:data.table", unload = TRUE) #data.table and tidyverse often conflict, so removing data.table when not in use
  
  #Renaming 'temp' as 'mi_X'
  assign(paste0("mi_", mi), temp)
  rm(temp, mi)
}

## Creating 'include' indicator
for (mi in 1:35) {
  #Creating temporary copy of mi_X
  eval(parse(text = paste0("temp <- mi_", mi)))
  
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
    print(paste0("MI: ", mi, ". Person: ", i, ". ", round((i/NROW(individuals)*100), 2), "%"))
    max_consecutive_count <- 0
    current_consecutive_count <- 0
    eligible_start_year <- NA
    eligible_end_year <- NA
    for (j in 1:nrow(temp[temp$unique_id == ind, ])) {
      if (temp[temp$unique_id == ind, ]$eligible[j] == 1) {
        current_consecutive_count <- current_consecutive_count + 1
        if (current_consecutive_count > max_consecutive_count) {
          max_consecutive_count <- current_consecutive_count
          eligible_start_year <- (temp[temp$unique_id == ind, ]$year[j - (current_consecutive_count - 1)]) - 4
          eligible_end_year <- temp[temp$unique_id == ind, ]$year[j]
        }
      } else {
        current_consecutive_count <- 0
      }
    }
    temp$max_period[temp$unique_id == ind] <- max_consecutive_count
    temp[temp$unique_id == ind, ]$include <- ifelse(temp[temp$unique_id == ind, ]$year >= eligible_start_year & temp[temp$unique_id == ind, ]$year <= eligible_end_year, 1, 0)
  }
  
  #Saving updated 'mi_X' after determining waves to keep for analyses
  assign(paste0("mi_", mi), temp)
  qsave(temp, paste0(directory, "MI/ui_project_mi_", mi, "_060823.qs"))
  rm(temp, mi, i, j)
}

#### Step 19 - Creating 'EQ_n1', 'EQ_0' and 'EQ_1' for second most recent, current/most recent, and subsequent Y in eligible years ####
#            - Creating 'EQ_n1_year', 'EQ_0_year', 'EQ_1_year' for years before/after current T
for (mi in 1:35) {
  print(paste0("Current MI: ", mi))
  eval(parse(text = paste0("temp <- mi_", mi)))
  
  temp %>% arrange(unique_id, year) -> temp
  library(data.table)
  temp <- as.data.table(temp)
  temp[, EQ_0 := ifelse(!is.na(eq), eq,               
                        ifelse(!is.na(lag(eq,1)), lag(eq,1), ifelse(!is.na(lag(eq,2)), lag(eq,2),
                                                                    ifelse(!is.na(lag(eq,3)), lag(eq,3), ifelse(!is.na(lag(eq,4)), lag(eq,4), ifelse(!is.na(lag(eq,5)), lag(eq,5),
                                                                                                                                                     ifelse(!is.na(lag(eq,6)), lag(eq,6), ifelse(!is.na(lag(eq,7)), lag(eq,7), ifelse(!is.na(lag(eq,8)), lag(eq,8), NA_real_))))))))), by = .(unique_id)]
  temp[, EQ_1 := ifelse(!is.na(lead(eq,1)), lead(eq,1), ifelse(!is.na(lead(eq,2)), lead(eq,2),
                                                               ifelse(!is.na(lead(eq,3)), lead(eq,3), ifelse(!is.na(lead(eq,4)), lead(eq,4), 
                                                                                                             ifelse(!is.na(lead(eq,5)), lead(eq,5), ifelse(!is.na(lead(eq,6)), lead(eq,6), 
                                                                                                                                                           ifelse(!is.na(lead(eq,7)), lead(eq,7), ifelse(!is.na(lead(eq,8)), lead(eq,8), NA_real_)))))))), by = .(unique_id)]
  temp[, EQ_0_year := ifelse(!is.na(eq), (year-year),               
                             ifelse(!is.na(lag(eq,1)), (lag(year,1) - year), 
                                    ifelse(!is.na(lag(eq,2)), (lag(year,2) - year), ifelse(!is.na(lag(eq,3)), (lag(year,3) - year), 
                                                                                           ifelse(!is.na(lag(eq,4)), (lag(year,4) - year), ifelse(!is.na(lag(eq,5)), (lag(year,5) - year),
                                                                                                                                                  ifelse(!is.na(lag(eq,6)), (lag(year,6) - year), ifelse(!is.na(lag(eq,7)), (lag(year,7) - year), 
                                                                                                                                                                                                         ifelse(!is.na(lag(eq,8)), (lag(year,8) - year), NA_real_))))))))), by = .(unique_id)]
  temp[, EQ_1_year := ifelse(!is.na(lead(eq,1)), (lead(year,1) - year), 
                             ifelse(!is.na(lead(eq,2)), (lead(year,2) - year),
                                    ifelse(!is.na(lead(eq,3)), (lead(year,3) - year), ifelse(!is.na(lead(eq,4)), (lead(year,4) - year), 
                                                                                             ifelse(!is.na(lead(eq,5)), (lead(year,5) - year), ifelse(!is.na(lead(eq,6)), (lead(year,6) - year), 
                                                                                                                                                      ifelse(!is.na(lead(eq,7)), (lead(year,7) - year), ifelse(!is.na(lead(eq,8)), (lead(year,8) - year), NA_real_)))))))), by = .(unique_id)]
  temp[, EQ_0_lags := ifelse(!is.na(eq), 0,        
                             ifelse(!is.na(lag(eq,1)), 1, ifelse(!is.na(lag(eq,2)), 2,
                                                                 ifelse(!is.na(lag(eq,3)), 3, ifelse(!is.na(lag(eq,4)), 4, ifelse(!is.na(lag(eq,5)), 5,
                                                                                                                                  ifelse(!is.na(lag(eq,6)), 6, ifelse(!is.na(lag(eq,7)), 7, ifelse(!is.na(lag(eq,8)), 8, NA_real_))))))))), by = .(unique_id)] #NOTE - EQ_0_lags only provides whole numbers 
  temp[, ':=' (EQ_n1 = fcase(EQ_0_lags %in% c(0) & !is.na(lag(eq, 1)), lag(eq,1),
                             EQ_0_lags %in% c(0:1) & !is.na(lag(eq, 2)), lag(eq, 2),
                             EQ_0_lags %in% c(0:2) & !is.na(lag(eq, 3)), lag(eq, 3),
                             EQ_0_lags %in% c(0:3) & !is.na(lag(eq, 4)), lag(eq, 4),
                             EQ_0_lags %in% c(0:4) & !is.na(lag(eq, 5)), lag(eq, 5),
                             EQ_0_lags %in% c(0:5) & !is.na(lag(eq, 6)), lag(eq, 6),
                             EQ_0_lags %in% c(0:6) & !is.na(lag(eq, 7)), lag(eq, 7),
                             EQ_0_lags %in% c(0:7) & !is.na(lag(eq, 8)), lag(eq, 8),
                             EQ_0_lags %in% c(1:8) & !is.na(lag(eq, 9)), lag(eq, 9),
                             EQ_0_lags %in% c(2:8) & !is.na(lag(eq, 10)), lag(eq, 10),
                             EQ_0_lags %in% c(3:8) & !is.na(lag(eq, 11)), lag(eq, 11),
                             EQ_0_lags %in% c(4:8) & !is.na(lag(eq, 12)), lag(eq, 12),
                             EQ_0_lags %in% c(5:8) & !is.na(lag(eq, 13)), lag(eq, 13),
                             EQ_0_lags %in% c(6:8) & !is.na(lag(eq, 14)), lag(eq, 14),
                             EQ_0_lags %in% c(7:8) & !is.na(lag(eq, 15)), lag(eq, 15),
                             EQ_0_lags %in% c(8) & !is.na(lag(eq, 16)), lag(eq, 16), default = NA_real_),
               EQ_n1_year = fcase(EQ_0_lags %in% c(0) & !is.na(lag(eq, 1)), (lag(year,1) - year),
                                  EQ_0_lags %in% c(0:1) & !is.na(lag(eq, 2)), (lag(year, 2) - year),
                                  EQ_0_lags %in% c(0:2) & !is.na(lag(eq, 3)), (lag(year, 3) - year),
                                  EQ_0_lags %in% c(0:3) & !is.na(lag(eq, 4)), (lag(year, 4) - year),
                                  EQ_0_lags %in% c(0:4) & !is.na(lag(eq, 5)), (lag(year, 5) - year),
                                  EQ_0_lags %in% c(0:5) & !is.na(lag(eq, 6)), (lag(year, 6) - year),
                                  EQ_0_lags %in% c(0:6) & !is.na(lag(eq, 7)), (lag(year, 7) - year),
                                  EQ_0_lags %in% c(0:7) & !is.na(lag(eq, 8)), (lag(year, 8) - year),
                                  EQ_0_lags %in% c(1:8) & !is.na(lag(eq, 9)), (lag(year, 9) - year),
                                  EQ_0_lags %in% c(2:8) & !is.na(lag(eq, 10)), (lag(year, 10) - year),
                                  EQ_0_lags %in% c(3:8) & !is.na(lag(eq, 11)), (lag(year, 11) - year),
                                  EQ_0_lags %in% c(4:8) & !is.na(lag(eq, 12)), (lag(year, 12) - year),
                                  EQ_0_lags %in% c(5:8) & !is.na(lag(eq, 13)), (lag(year, 13) - year),
                                  EQ_0_lags %in% c(6:8) & !is.na(lag(eq, 14)), (lag(year, 14) - year),
                                  EQ_0_lags %in% c(7:8) & !is.na(lag(eq, 15)), (lag(year, 15) - year),
                                  EQ_0_lags %in% c(8) & !is.na(lag(eq, 16)), (lag(year, 16) - year), default = NA_real_)), by = .(unique_id)]
  
  temp <- as.data.frame(temp)
  detach("package:data.table", unload = TRUE) #data.table and tidyverse often conflict, so removing data.table when not in use
  
  #Saving updated 'mi_X' after determining waves to keep for analyses
  assign(paste0("mi_", mi), temp)
  qsave(temp, paste0(directory, "MI/ui_project_mi_", mi, "_060823.qs"))
  rm(temp, mi)
}

#### Step 20 - Create 'Y_lead', 'Y_lead_year', 'EQ_0_lags', 'Y_base', and 'Y_base_year' variables ####
#NOTE - Ideally for all times t in the pre-post comparison, Y would come from T=t+1. 
#       If unemployed at T=t+1 there is no Y, though, so we draw from EQ_0 (nearest earlier) or EQ_1 (nearest later) as appropriate
#NOTE - Where A = {0,1} with 0=employed, 1=unemployed, we compare A(t0-t1) with our DiD = E[Yt2+ - Yt0|At1=1, At0=0] - E[Yt2 - Yt1|At1=0, At0=0]:
#     - Comparing A(t0-t1), for t0 (pre) obs. where A(t1) = 0, Y_lead = Y(t1)
#                           for t0 (pre) obs. where A(t1) = 1, Y_lead = nearest known Y(<t1)
#                           for t1 (post) obs. where A(t2) = 0, Y_lead = Y(t2)
#                           for t1 (post) obs. where A(t2) = 1, Y_lead = nearest known Y(>t2)
#     - For each i at t: 
#         - If (i,t-1) was a               'treated observation' with A{t-1,t,t+1} for (i,t) = {1,A,A}: Y_lead = EQ_0
#         - If (i,t) is a newly            'treated observation' with A{t-1,t,t+1} for (i,t) = {0,1,A}: Y_lead = EQ_1                 
#         - If (i,t+1) will be a newly     'treated observation' with A{t-1,t,t+1} for (i,t) = {0,0,1}: Y_lead = EQ_0
#         - if (i,t+1) will not be a newly 'treated observation' with A{t-1,t,t+1} for (i,t) = {0,0,0}: Y_lead = EQ_1
#     - Based on the above, I should code Y_lead for (i,t) as:
#         - Y_lead = fcase(lag(A,1) == 1, EQ_0, 
#                          A == 1, EQ_1,
#                          lead(A,1) == 1, EQ_0,
#                          lag(A,1) == 0 & A == 0 & lead(A,1) == 0, EQ_1,
#                          default = NA_real_)
#     - To have Y_base always precede Y_lead, if Y_lead = EQ_1, Y_base = EQ_0, and if Y_lead = EQ_0, Y_base = the latest Y(t) earlier than EQ_0
#         - 'EQ_0_lags' indicates how many lags prior EQ_0 comes from, so Y_base where Y_lead = EQ_0 must check Y(t) < Y(T<EQ_0_lags)
for (mi in 1:35) {
  print(paste0("Current MI: ", mi))
  eval(parse(text = paste0("temp <- mi_", mi)))
  
  temp %>% arrange(unique_id, year) -> temp
  
  library(data.table)
  temp <- as.data.table(temp)
  temp[, Y_lead := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_1,   #0,0,0 EQ_1
                         lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_0,   #0,0,1 EQ_0  
                         lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_1,                                    #0,1,A EQ_1 
                         lag(unemp_status_lead,1) == 1, EQ_0,                                                             #1,A,A EQ_0
                         lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_1, #0,0,- EQ_1
                         is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_0,                                  #-,0,0 EQ_0
                         is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_1,                                  #-,1,A EQ_1
                         default = NA_real_), by = .(unique_id)]
  temp[, Y_lead_year := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_1_year,   #0,0,0 EQ_1
                              lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_0_year,   #0,0,1 EQ_0  
                              lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_1_year,                                    #0,1,A EQ_1 
                              lag(unemp_status_lead,1) == 1, EQ_0_year,                                                             #1,A,A EQ_0
                              lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_1_year, #0,0,- EQ_1
                              is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_0_year,                                  #-,0,0 EQ_0
                              is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_1_year,                                  #-,1,A EQ_1
                              default = NA_real_), by = .(unique_id)]
  temp[, Y_base := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_0,   #0,0,0 EQ_1
                         lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_n1,  #0,0,1 EQ_0  
                         lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_0,                                    #0,1,A EQ_1 
                         lag(unemp_status_lead,1) == 1, EQ_n1,                                                            #1,A,A EQ_0
                         lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_0, #0,0,- EQ_1
                         is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_n1,                                 #-,0,0 EQ_0
                         is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_0,                                  #-,1,A EQ_1
                         default = NA_real_), by = .(unique_id)]
  temp[, Y_base_year := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_0_year,   #0,0,0 EQ_1
                              lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_n1_year,  #0,0,1 EQ_0  
                              lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_0_year,                                    #0,1,A EQ_1 
                              lag(unemp_status_lead,1) == 1, EQ_n1_year,                                                            #1,A,A EQ_0
                              lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_0_year, #0,0,- EQ_1
                              is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_n1_year,                                 #-,0,0 EQ_0
                              is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_0_year,                                  #-,1,A EQ_1
                              default = NA_real_), by = .(unique_id)]
  
  temp <- as.data.frame(temp)
  detach("package:data.table", unload = TRUE)
  
  #Saving updated 'mi_X' after determining waves to keep for analyses
  assign(paste0("mi_", mi), temp)
  qsave(temp, paste0(directory, "MI/ui_project_mi_", mi, "_060823.qs"))
  rm(temp, mi)
}

#### Step 21 - Create 1-lagged versions of time-varying confounding factors and 'pre-post' indicator ####
#             'prepost' indicator for A{t,t-1} 'CE'/'CU' = Consistently Employed/Unemployed, 'RU'/'RE' = Recently Unemployed/Employed
for (mi in 1:35) {
  print(paste0("Current MI: ", mi))
  eval(parse(text = paste0("temp <- mi_", mi)))
  
  temp %>% arrange(unique_id, year) -> temp
  
  library(data.table)
  temp <- as.data.table(temp)
  temp[, ':=' (SRH_lag = lag(SRH,1),
               marital_status_lag = lag(marital_status,1),
               education_lag = lag(education,1),
               disabl_limits_work_lag = lag(disabl_limits_work,1),
               occupation_broad_lag = lag(occupation_broad,1),
               family_income_lag = lag(family_income,1),
               family_wealth_no_home_equity_lag = lag(family_wealth_no_home_equity,1),
               prepost = fifelse(unemp_status_lead == 0 & lag(unemp_status_lead,1) == 0, "CE",
                                 fifelse(unemp_status_lead == 1 & lag(unemp_status_lead,1) == 0, "RU", 
                                         fifelse(unemp_status_lead == 0 & lag(unemp_status_lead,1) == 1, "RE",
                                                 fifelse(unemp_status_lead == 1 & lag(unemp_status_lead,1) == 1, "CU", "Unknown"))))), by = .(unique_id)]
  temp <- as.data.frame(temp)
  detach("package:data.table", unload = TRUE)
  
  #Saving updated 'mi_X' after determining waves to keep for analyses
  assign(paste0("mi_", mi), temp)
  qsave(temp, paste0(directory, "MI/ui_project_mi_", mi, "_060823.qs"))
  rm(temp, mi)
}

#### Step 22 - Repeating Steps 17-21 to create an equivalent 'complete case' dataset ####
## Part 1 - Reformatting the complete case dataset 'incomplete' 
#NOTE - Here we: 1) impute directly known vars, 
#                2) change vars to minimum sensible coding (e.g. those imputed as 'not unemployed in 12 months' cannot receive UI in 6 months post-unemployment)
#                3) Excluding those with unknown state or in a given year residing in a foreign country or US territory
#Creating a temporary dataset
temp <- incomplete

#Resetting 12-month unemployment status, UI receipt, and months unemployed pre-interview to sensible minimums
temp %>% mutate(unemp_status_lead = ifelse(is.na(unemp_status_lead) & employment_status == "Unemployed", 1, unemp_status_lead),
                ui_receipt_6mo_post_unemp = ifelse(unemp_status_lead == 0, 0, ui_receipt_6mo_post_unemp),
                mo_unemp_pre_int_min = ifelse(mo_unemp_pre_int_min < 0, 0, mo_unemp_pre_int_min)) -> temp

#All those with missing Z-score GSP per capita are in states, GSP-relevant years, and GSP-relevant quarters with no non-missing observations
#Therefore, we just remove the Z-score GSP per capita variable and re-merge in actual state-level GSP per capita information
#Creating needed 'gsp_ur_year' and 'gsp_ur_quarter' variables for merging
temp %>% mutate(gsp_ur_year = case_when(interview_year != (year+2001) ~ (year+2001),
                                        unemp_status_lead == 0 & interview_month == 1 ~ (year-1+2001),
                                        unemp_status_lead == 0 & interview_month > 1 ~ (year+2001),
                                        is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month == 1 ~ (year-1+2001),
                                        is.na(first_mo_unemp_post_int) & unemp_status_lead == 1 & interview_month > 1 ~ (year+2001),
                                        unemp_status_lead == 1 & first_mo_unemp_post_int == 1 ~ (year-1+2001),
                                        unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(2:13) ~ (year+2001),
                                        TRUE ~ (year+1+2001)),
                gsp_ur_quarter = case_when(gsp_ur_year <= 2004 ~ 1, #We do not have monthly Real GSP per Capita information by state pre-2005 
                                           (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month == 1 ~ 4,
                                           (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month <= 4 ~ 1,
                                           (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month <= 7 ~ 2,
                                           (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month <= 10 ~ 3,
                                           (unemp_status_lead == 0 | is.na(first_mo_unemp_post_int)) & interview_month ~ 4,
                                           unemp_status_lead == 1 & first_mo_unemp_post_int == 1 ~ 4,
                                           unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(2:4,14:16) ~ 1,
                                           unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(5:7,17:19) ~ 2,
                                           unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(8:10,20:22) ~ 3,
                                           unemp_status_lead == 1 & first_mo_unemp_post_int %in% c(11:13,23) ~ 4,
                                           TRUE ~ 99)) -> temp
temp$gsp_per_cap <- NULL

#Loading, formatting, merging gsp_per_cap with temp, then removing the state GSP dataset
state_gsp_data <- read.csv(paste0(directory, "BEA_GSP_9719.csv")) #State GSP per Capita data
state_gsp_data <- state_gsp_data[,c(2:5,8)]
state_gsp_data %>% mutate(state = case_when(fips == 01 ~ "AL", fips == 02 ~ "AK", fips == 04 ~ "AZ", fips == 05 ~ "AR", 
                                            fips == 06 ~ "CA", fips == 08 ~ "CO", fips == 09 ~ "CT", fips == 10 ~ "DE",
                                            fips == 11 ~ "DC", fips == 12 ~ "FL", fips == 13 ~ "GA", fips == 15 ~ "HI",
                                            fips == 16 ~ "ID", fips == 17 ~ "IL", fips == 18 ~ "IN", fips == 19 ~ "IA",
                                            fips == 20 ~ "KS", fips == 21 ~ "KY", fips == 22 ~ "LA", fips == 23 ~ "ME",
                                            fips == 24 ~ "MD", fips == 25 ~ "MA", fips == 26 ~ "MI", fips == 27 ~ "MN",
                                            fips == 28 ~ "MS", fips == 29 ~ "MO", fips == 30 ~ "MT", fips == 31 ~ "NE",
                                            fips == 32 ~ "NV", fips == 33 ~ "NH", fips == 34 ~ "NJ", fips == 35 ~ "NM",
                                            fips == 36 ~ "NY", fips == 37 ~ "NC", fips == 38 ~ "ND", fips == 39 ~ "OH",
                                            fips == 40 ~ "OK", fips == 41 ~ "OR", fips == 42 ~ "PA", fips == 44 ~ "RI",
                                            fips == 45 ~ "SC", fips == 46 ~ "SD", fips == 47 ~ "TN", fips == 48 ~ "TX",
                                            fips == 49 ~ "UT", fips == 50 ~ "VT", fips == 51 ~ "VA", fips == 53 ~ "WA",
                                            fips == 54 ~ "WV", fips == 55 ~ "WI", fips == 56 ~ "WY", TRUE ~ "US Territory or Foreign Country")) -> state_gsp_data
names(state_gsp_data)[c(3,4)] <- c("gsp_ur_year", "gsp_ur_quarter")
temp <- merge(temp, state_gsp_data[,c(6,3:5)], by = c("state", "gsp_ur_year", "gsp_ur_quarter"), all.x = T)
temp <- distinct(temp)
rm(state_gsp_data)

#Excluding empty rows or those residing in a given year in a US Territory or foreign country
temp <- temp[!is.na(temp$unique_id), ]
temp <- temp[temp$state != "US Territory or Foreign Country", ]
temp <- temp[!is.na(temp$gsp_per_cap), ]

#Renaming 'temp' as 'mi_X'
assign("incomplete", temp)
rm(temp)

## Part 2 - Creating wave-specific eligibility variable (eligible) and overall waves eligibile per participant indicator (include) based on known 1-lag and 2-lag waves 2 and 4 years prior
#NOTE - Waves are eligible if, for the same respondent, there are 2 prior waves of data which are from 2 and 4 years prior to the eligible wave
#     - The total set of eligible waves per respondent ('include') is the largest set of consecutive waves where 'eligible' = 1 as well as the 2 prior waves
## Creating 'eligible' indicator
temp <- incomplete

temp %>% arrange(unique_id, year) -> temp
library(data.table)
temp <- as.data.table(temp)
temp[, eligible := fcase(is.na(lag(year,1)) | is.na(lag(year,2)), 0L,
                         lag(year,1) != (year-2) | lag(year,2) != (year-4), 0L,
                         default = 1L), by = .(unique_id)]
temp <- as.data.frame(temp)
detach("package:data.table", unload = TRUE) #data.table and tidyverse often conflict, so removing data.table when not in use

#Renaming 'temp' as 'incomplete'
assign("incomplete", temp)
rm(temp)

## Creating 'include' indicator
temp <- incomplete

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
  print(paste0("Person: ", i, ". ", round((i/NROW(individuals)*100), 2), "%"))
  max_consecutive_count <- 0
  current_consecutive_count <- 0
  eligible_start_year <- NA
  eligible_end_year <- NA
  for (j in 1:nrow(temp[temp$unique_id == ind, ])) {
    if (temp[temp$unique_id == ind, ]$eligible[j] == 1) {
      current_consecutive_count <- current_consecutive_count + 1
      if (current_consecutive_count > max_consecutive_count) {
        max_consecutive_count <- current_consecutive_count
        eligible_start_year <- (temp[temp$unique_id == ind, ]$year[j - (current_consecutive_count - 1)]) - 4
        eligible_end_year <- temp[temp$unique_id == ind, ]$year[j]
      }
    } else {
      current_consecutive_count <- 0
    }
  }
  temp$max_period[temp$unique_id == ind] <- max_consecutive_count
  temp[temp$unique_id == ind, ]$include <- ifelse(temp[temp$unique_id == ind, ]$year >= eligible_start_year & temp[temp$unique_id == ind, ]$year <= eligible_end_year, 1, 0)
}

#Saving updated 'incomplete' after determining waves to keep for analyses
assign("incomplete", temp)
rm(temp, i, j, eligible_start_year, eligible_end_year, max_consecutive_count, current_consecutive_count, ind, individuals)

## Part 3 - Creating 'EQ_n1', 'EQ_0' and 'EQ_1' for second most recent, current/most recent, and subsequent Y in eligible years
#         - Creating 'EQ_n1_year', 'EQ_0_year', 'EQ_1_year' for years before/after current T
temp <- incomplete

temp %>% arrange(unique_id, year) -> temp
library(data.table)
temp <- as.data.table(temp)
temp[, EQ_0 := ifelse(!is.na(eq), eq,               
                      ifelse(!is.na(lag(eq,1)), lag(eq,1), ifelse(!is.na(lag(eq,2)), lag(eq,2),
                                                                  ifelse(!is.na(lag(eq,3)), lag(eq,3), ifelse(!is.na(lag(eq,4)), lag(eq,4), ifelse(!is.na(lag(eq,5)), lag(eq,5),
                                                                                                                                                   ifelse(!is.na(lag(eq,6)), lag(eq,6), ifelse(!is.na(lag(eq,7)), lag(eq,7), ifelse(!is.na(lag(eq,8)), lag(eq,8), NA_real_))))))))), by = .(unique_id)]
temp[, EQ_1 := ifelse(!is.na(lead(eq,1)), lead(eq,1), ifelse(!is.na(lead(eq,2)), lead(eq,2),
                                                             ifelse(!is.na(lead(eq,3)), lead(eq,3), ifelse(!is.na(lead(eq,4)), lead(eq,4), 
                                                                                                           ifelse(!is.na(lead(eq,5)), lead(eq,5), ifelse(!is.na(lead(eq,6)), lead(eq,6), 
                                                                                                                                                         ifelse(!is.na(lead(eq,7)), lead(eq,7), ifelse(!is.na(lead(eq,8)), lead(eq,8), NA_real_)))))))), by = .(unique_id)]
temp[, EQ_0_year := ifelse(!is.na(eq), (year-year),               
                           ifelse(!is.na(lag(eq,1)), (lag(year,1) - year), 
                                  ifelse(!is.na(lag(eq,2)), (lag(year,2) - year), ifelse(!is.na(lag(eq,3)), (lag(year,3) - year), 
                                                                                         ifelse(!is.na(lag(eq,4)), (lag(year,4) - year), ifelse(!is.na(lag(eq,5)), (lag(year,5) - year),
                                                                                                                                                ifelse(!is.na(lag(eq,6)), (lag(year,6) - year), ifelse(!is.na(lag(eq,7)), (lag(year,7) - year), 
                                                                                                                                                                                                       ifelse(!is.na(lag(eq,8)), (lag(year,8) - year), NA_real_))))))))), by = .(unique_id)]
temp[, EQ_1_year := ifelse(!is.na(lead(eq,1)), (lead(year,1) - year), 
                           ifelse(!is.na(lead(eq,2)), (lead(year,2) - year),
                                  ifelse(!is.na(lead(eq,3)), (lead(year,3) - year), ifelse(!is.na(lead(eq,4)), (lead(year,4) - year), 
                                                                                           ifelse(!is.na(lead(eq,5)), (lead(year,5) - year), ifelse(!is.na(lead(eq,6)), (lead(year,6) - year), 
                                                                                                                                                    ifelse(!is.na(lead(eq,7)), (lead(year,7) - year), ifelse(!is.na(lead(eq,8)), (lead(year,8) - year), NA_real_)))))))), by = .(unique_id)]
temp[, EQ_0_lags := ifelse(!is.na(eq), 0,        
                           ifelse(!is.na(lag(eq,1)), 1, ifelse(!is.na(lag(eq,2)), 2,
                                                               ifelse(!is.na(lag(eq,3)), 3, ifelse(!is.na(lag(eq,4)), 4, ifelse(!is.na(lag(eq,5)), 5,
                                                                                                                                ifelse(!is.na(lag(eq,6)), 6, ifelse(!is.na(lag(eq,7)), 7, ifelse(!is.na(lag(eq,8)), 8, NA_real_))))))))), by = .(unique_id)] #NOTE - EQ_0_lags only provides whole numbers 
temp[, ':=' (EQ_n1 = fcase(EQ_0_lags %in% c(0) & !is.na(lag(eq, 1)), lag(eq,1),
                           EQ_0_lags %in% c(0:1) & !is.na(lag(eq, 2)), lag(eq, 2),
                           EQ_0_lags %in% c(0:2) & !is.na(lag(eq, 3)), lag(eq, 3),
                           EQ_0_lags %in% c(0:3) & !is.na(lag(eq, 4)), lag(eq, 4),
                           EQ_0_lags %in% c(0:4) & !is.na(lag(eq, 5)), lag(eq, 5),
                           EQ_0_lags %in% c(0:5) & !is.na(lag(eq, 6)), lag(eq, 6),
                           EQ_0_lags %in% c(0:6) & !is.na(lag(eq, 7)), lag(eq, 7),
                           EQ_0_lags %in% c(0:7) & !is.na(lag(eq, 8)), lag(eq, 8),
                           EQ_0_lags %in% c(1:8) & !is.na(lag(eq, 9)), lag(eq, 9),
                           EQ_0_lags %in% c(2:8) & !is.na(lag(eq, 10)), lag(eq, 10),
                           EQ_0_lags %in% c(3:8) & !is.na(lag(eq, 11)), lag(eq, 11),
                           EQ_0_lags %in% c(4:8) & !is.na(lag(eq, 12)), lag(eq, 12),
                           EQ_0_lags %in% c(5:8) & !is.na(lag(eq, 13)), lag(eq, 13),
                           EQ_0_lags %in% c(6:8) & !is.na(lag(eq, 14)), lag(eq, 14),
                           EQ_0_lags %in% c(7:8) & !is.na(lag(eq, 15)), lag(eq, 15),
                           EQ_0_lags %in% c(8) & !is.na(lag(eq, 16)), lag(eq, 16), default = NA_real_),
             EQ_n1_year = fcase(EQ_0_lags %in% c(0) & !is.na(lag(eq, 1)), (lag(year,1) - year),
                                EQ_0_lags %in% c(0:1) & !is.na(lag(eq, 2)), (lag(year, 2) - year),
                                EQ_0_lags %in% c(0:2) & !is.na(lag(eq, 3)), (lag(year, 3) - year),
                                EQ_0_lags %in% c(0:3) & !is.na(lag(eq, 4)), (lag(year, 4) - year),
                                EQ_0_lags %in% c(0:4) & !is.na(lag(eq, 5)), (lag(year, 5) - year),
                                EQ_0_lags %in% c(0:5) & !is.na(lag(eq, 6)), (lag(year, 6) - year),
                                EQ_0_lags %in% c(0:6) & !is.na(lag(eq, 7)), (lag(year, 7) - year),
                                EQ_0_lags %in% c(0:7) & !is.na(lag(eq, 8)), (lag(year, 8) - year),
                                EQ_0_lags %in% c(1:8) & !is.na(lag(eq, 9)), (lag(year, 9) - year),
                                EQ_0_lags %in% c(2:8) & !is.na(lag(eq, 10)), (lag(year, 10) - year),
                                EQ_0_lags %in% c(3:8) & !is.na(lag(eq, 11)), (lag(year, 11) - year),
                                EQ_0_lags %in% c(4:8) & !is.na(lag(eq, 12)), (lag(year, 12) - year),
                                EQ_0_lags %in% c(5:8) & !is.na(lag(eq, 13)), (lag(year, 13) - year),
                                EQ_0_lags %in% c(6:8) & !is.na(lag(eq, 14)), (lag(year, 14) - year),
                                EQ_0_lags %in% c(7:8) & !is.na(lag(eq, 15)), (lag(year, 15) - year),
                                EQ_0_lags %in% c(8) & !is.na(lag(eq, 16)), (lag(year, 16) - year), default = NA_real_)), by = .(unique_id)]

temp <- as.data.frame(temp)
detach("package:data.table", unload = TRUE) #data.table and tidyverse often conflict, so removing data.table when not in use

#Replacing updated 'incomplete' after determining waves to keep for analyses
assign("incomplete", temp)
rm(temp)

## Step 4 - Create 'Y_lead', 'Y_lead_year', 'EQ_0_lags', 'Y_base', and 'Y_base_year' variables
#NOTE - Ideally for all times t in the pre-post comparison, Y would come from T=t+1. 
#       If unemployed at T=t+1 there is no Y, though, so we draw from EQ_0 (nearest earlier) or EQ_1 (nearest later) as appropriate
#NOTE - Where A = {0,1} with 0=employed, 1=unemployed, we compare A(t0-t1) with our DiD = E[Yt2+ - Yt0|At1=1, At0=0] - E[Yt2 - Yt1|At1=0, At0=0]:
#     - Comparing A(t0-t1), for t0 (pre) obs. where A(t1) = 0, Y_lead = Y(t1)
#                           for t0 (pre) obs. where A(t1) = 1, Y_lead = nearest known Y(<t1)
#                           for t1 (post) obs. where A(t2) = 0, Y_lead = Y(t2)
#                           for t1 (post) obs. where A(t2) = 1, Y_lead = nearest known Y(>t2)
#     - For each i at t: 
#         - If (i,t-1) was a               'treated observation' with A{t-1,t,t+1} for (i,t) = {1,A,A}: Y_lead = EQ_0
#         - If (i,t) is a newly            'treated observation' with A{t-1,t,t+1} for (i,t) = {0,1,A}: Y_lead = EQ_1                 
#         - If (i,t+1) will be a newly     'treated observation' with A{t-1,t,t+1} for (i,t) = {0,0,1}: Y_lead = EQ_0
#         - if (i,t+1) will not be a newly 'treated observation' with A{t-1,t,t+1} for (i,t) = {0,0,0}: Y_lead = EQ_1
#     - Based on the above, I should code Y_lead for (i,t) as:
#         - Y_lead = fcase(lag(A,1) == 1, EQ_0, 
#                          A == 1, EQ_1,
#                          lead(A,1) == 1, EQ_0,
#                          lag(A,1) == 0 & A == 0 & lead(A,1) == 0, EQ_1,
#                          default = NA_real_)
#     - To have Y_base always precede Y_lead, if Y_lead = EQ_1, Y_base = EQ_0, and if Y_lead = EQ_0, Y_base = the latest Y(t) earlier than EQ_0
#         - 'EQ_0_lags' indicates how many lags prior EQ_0 comes from, so Y_base where Y_lead = EQ_0 must check Y(t) < Y(T<EQ_0_lags)
temp <- incomplete

temp %>% arrange(unique_id, year) -> temp

library(data.table)
temp <- as.data.table(temp)
temp[, Y_lead := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_1,   #0,0,0 EQ_1
                       lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_0,   #0,0,1 EQ_0  
                       lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_1,                                    #0,1,A EQ_1 
                       lag(unemp_status_lead,1) == 1, EQ_0,                                                             #1,A,A EQ_0
                       lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_1, #0,0,- EQ_1
                       is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_0,                                  #-,0,0 EQ_0
                       is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_1,                                  #-,1,A EQ_1
                       default = NA_real_), by = .(unique_id)]
temp[, Y_lead_year := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_1_year,   #0,0,0 EQ_1
                            lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_0_year,   #0,0,1 EQ_0  
                            lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_1_year,                                    #0,1,A EQ_1 
                            lag(unemp_status_lead,1) == 1, EQ_0_year,                                                             #1,A,A EQ_0
                            lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_1_year, #0,0,- EQ_1
                            is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_0_year,                                  #-,0,0 EQ_0
                            is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_1_year,                                  #-,1,A EQ_1
                            default = NA_real_), by = .(unique_id)]
temp[, Y_base := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_0,   #0,0,0 EQ_1
                       lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_n1,  #0,0,1 EQ_0  
                       lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_0,                                    #0,1,A EQ_1 
                       lag(unemp_status_lead,1) == 1, EQ_n1,                                                            #1,A,A EQ_0
                       lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_0, #0,0,- EQ_1
                       is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_n1,                                 #-,0,0 EQ_0
                       is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_0,                                  #-,1,A EQ_1
                       default = NA_real_), by = .(unique_id)]
temp[, Y_base_year := fcase(lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 0, EQ_0_year,   #0,0,0 EQ_1
                            lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & lead(unemp_status_lead,1) == 1, EQ_n1_year,  #0,0,1 EQ_0  
                            lag(unemp_status_lead,1) == 0 & unemp_status_lead == 1, EQ_0_year,                                    #0,1,A EQ_1 
                            lag(unemp_status_lead,1) == 1, EQ_n1_year,                                                            #1,A,A EQ_0
                            lag(unemp_status_lead,1) == 0 & unemp_status_lead == 0 & is.na(lead(unemp_status_lead,1)), EQ_0_year, #0,0,- EQ_1
                            is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 0, EQ_n1_year,                                 #-,0,0 EQ_0
                            is.na(lag(unemp_status_lead,1)) & unemp_status_lead == 1, EQ_0_year,                                  #-,1,A EQ_1
                            default = NA_real_), by = .(unique_id)]

temp <- as.data.frame(temp)
detach("package:data.table", unload = TRUE)

#Replacing updated 'incomplete' after determining waves to keep for analyses
assign("incomplete", temp)
rm(temp)

## Part 5 - Create 1-lagged versions of time-varying confounding factors and 'pre-post' indicator
#         - 'prepost' indicator for A{t,t-1} 'CE'/'CU' = Consistently Employed/Unemployed, 'RU'/'RE' = Recently Unemployed/Employed
temp <- incomplete

temp %>% arrange(unique_id, year) -> temp

library(data.table)
temp <- as.data.table(temp)
temp[, ':=' (SRH_lag = lag(SRH,1),
             marital_status_lag = lag(marital_status,1),
             education_lag = lag(education,1),
             disabl_limits_work_lag = lag(disabl_limits_work,1),
             occupation_broad_lag = lag(occupation_broad,1),
             family_income_lag = lag(family_income,1),
             family_wealth_no_home_equity_lag = lag(family_wealth_no_home_equity,1),
             prepost = fifelse(unemp_status_lead == 0 & lag(unemp_status_lead,1) == 0, "CE",
                               fifelse(unemp_status_lead == 1 & lag(unemp_status_lead,1) == 0, "RU", 
                                       fifelse(unemp_status_lead == 0 & lag(unemp_status_lead,1) == 1, "RE",
                                               fifelse(unemp_status_lead == 1 & lag(unemp_status_lead,1) == 1, "CU", "Unknown"))))), by = .(unique_id)]
temp <- as.data.frame(temp)
detach("package:data.table", unload = TRUE)

#Saving updated 'mi_X' after determining waves to keep for analyses
assign("incomplete", temp)
qsave(incomplete, paste0(directory, "ui_formatted_complete_case_060823.qs"))
rm(temp)
