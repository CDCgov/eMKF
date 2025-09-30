# ------------------------------------------------------------------------------
# This file illustrates the use of the enhanced Modified Kalman Filter (eMKF) macro using data from the National Survey of Family 
# Growth (NSFG) in the following study to evaluate the utility of the eMKF in producing estimates of the risk of pregnancy loss 
# for subgroups of US women with small sample sizes to examine recent trends:
#
# Forrest SE, Rossen LM, Ahrens KA. Trends in Risk of Pregnancy Loss Among US Women by Metropolitan Status, 2000-2018. 
# Paediatric and Perinatal Epidemiology. 2025. DOI:10.1111/ppe.70066.
#
# It includes pregnancies reported in the 2006-2010, 2011-2013, 2013-2015, 2015-2017, and 2017-2019 NSFG survey periods.
# This file creates the eMKF input dataset to be read into Example-eMKF-NSFG-Data.sas to run the eMKF. The eMKF input dataset:
#
# Excludes:
# - Induced abortions and ongoing pregnancies
# - Pregnancies conceived before 2000 and after 2018
# - Pregnancies with outcomes occurring at ages younger than 15 or older than 44 years
#
# Tabulates pregnancy-level data by:
# - Maternal age group (15-19, 20-24, 25-29, 30-34, and 35-44 years)
# - Metropolitan (metro) status (metropolitan and nonmetropolitan)
# - Conception year interval (2000-2002, 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016, and 2017-2018)
#
# After running the eMKF, the eMKF output from Example-eMKF-NSFG-Data.sas is read back into this file to complete the analysis
#
# Technical guidance for using the enhanced MKF macro is available from:
#
# Talih M, Rossen LM, Patel P, Earp M, Parker JD. Technical Guidance for Using the Modified Kalman Filter 
# in Small Domain Estimation at the National Center for Health Statistics. National Center for Health Statistics. 
# Vital Health Stat 2(209). 2024. DOI: 10.15620/cdc:157496.
#
# The data is assumed to be in user directory: ..\eMKF\MKFdata
# ------------------------------------------------------------------------------

library(tidyverse)
library(readr)
library(survey)
library(rms)
library(foreach) 
library(ggplot2) 
library(patchwork) 
set.seed(4298)

# ------------------------------------------------------------------------------
# Create cohort
# ------------------------------------------------------------------------------

# Read in CSV file: combined 2006-2010, 2011-2013, 2013-2015, 2015-2017, and 2017-2019 NSFG data
nsfg <- read.csv("...\eMKF\MKFdata\nsfg_2006-2019.csv")

# Recode variables for consistency across survey periods ----------------------------

# Recode DATECON (year pregnancy began) for 2006-2010, 2011-2013, and 2013-2015 NSFG files
nsfg <- nsfg %>%
  mutate(year_con = case_when(
    survey %in% c("2006-2010", "2011-2013", "2013-2015") ~ case_when(
      DATECON >= 1153 & DATECON <= 1164 ~ 1996,
      DATECON >= 1165 & DATECON <= 1176 ~ 1997,
      DATECON >= 1177 & DATECON <= 1188 ~ 1998,
      DATECON >= 1189 & DATECON <= 1200 ~ 1999,
      DATECON >= 1201 & DATECON <= 1212 ~ 2000,
      DATECON >= 1213 & DATECON <= 1224 ~ 2001,
      DATECON >= 1225 & DATECON <= 1236 ~ 2002,
      DATECON >= 1237 & DATECON <= 1248 ~ 2003,
      DATECON >= 1249 & DATECON <= 1260 ~ 2004,
      DATECON >= 1261 & DATECON <= 1272 ~ 2005,
      DATECON >= 1273 & DATECON <= 1284 ~ 2006,
      DATECON >= 1285 & DATECON <= 1296 ~ 2007,
      DATECON >= 1297 & DATECON <= 1308 ~ 2008,
      DATECON >= 1309 & DATECON <= 1320 ~ 2009,
      DATECON >= 1321 & DATECON <= 1332 ~ 2010,
      DATECON >= 1333 & DATECON <= 1344 ~ 2011,
      DATECON >= 1345 & DATECON <= 1356 ~ 2012,
      DATECON >= 1357 & DATECON <= 1368 ~ 2013,
      DATECON >= 1369 & DATECON <= 1380 ~ 2014,
      DATECON >= 1381 & DATECON <= 1392 ~ 2015,
      DATECON >= 1393 & DATECON <= 1404 ~ 2016,
      DATECON >= 1405 & DATECON <= 1416 ~ 2017,
      DATECON >= 1417 & DATECON <= 1428 ~ 2018,
      DATECON >= 1429 & DATECON <= 1440 ~ 2019,
      TRUE ~ NA_real_),
    survey %in% c("2015-2017", "2017-2019") ~ DATECON,
    TRUE ~ NA_real_))

# Recode AGECON (age at conception) for 2006-2010, 2011-2013, and 2013-2015 NSFG files
nsfg <- nsfg %>%
  mutate(age_con = case_when(
    survey %in% c("2006-2010", "2011-2013", "2013-2015") ~ case_when(
      AGECON >= 1450 & AGECON <= 1599 ~ 15,
      AGECON >= 1600 & AGECON <= 1699 ~ 16,
      AGECON >= 1700 & AGECON <= 1799 ~ 17,
      AGECON >= 1800 & AGECON <= 1899 ~ 18,
      AGECON >= 1900 & AGECON <= 1999 ~ 19,
      AGECON >= 2000 & AGECON <= 2099 ~ 20,
      AGECON >= 2100 & AGECON <= 2199 ~ 21,
      AGECON >= 2200 & AGECON <= 2299 ~ 22,
      AGECON >= 2300 & AGECON <= 2399 ~ 23,
      AGECON >= 2400 & AGECON <= 2499 ~ 24,
      AGECON >= 2500 & AGECON <= 2599 ~ 25,
      AGECON >= 2600 & AGECON <= 2699 ~ 26,
      AGECON >= 2700 & AGECON <= 2799 ~ 27,
      AGECON >= 2800 & AGECON <= 2899 ~ 28,
      AGECON >= 2900 & AGECON <= 2999 ~ 29,
      AGECON >= 3000 & AGECON <= 3099 ~ 30,
      AGECON >= 3100 & AGECON <= 3199 ~ 31,
      AGECON >= 3200 & AGECON <= 3299 ~ 32,
      AGECON >= 3300 & AGECON <= 3399 ~ 33,
      AGECON >= 3400 & AGECON <= 3499 ~ 34,
      AGECON >= 3500 & AGECON <= 3599 ~ 35,
      AGECON >= 3600 & AGECON <= 3699 ~ 36,
      AGECON >= 3700 & AGECON <= 3799 ~ 37,
      AGECON >= 3800 & AGECON <= 3899 ~ 38,
      AGECON >= 3900 & AGECON <= 3999 ~ 39,
      AGECON >= 4000 & AGECON <= 4099 ~ 40,
      AGECON >= 4100 & AGECON <= 4199 ~ 41,
      AGECON >= 4200 & AGECON <= 4299 ~ 42,
      AGECON >= 4300 & AGECON <= 4399 ~ 43,
      AGECON >= 4400 & AGECON <= 4550 ~ 44,
      TRUE ~ NA_real_),
    survey %in% c("2015-2017", "2017-2019") ~ AGECON,
    TRUE ~ NA_real_))

# Subset using inclusion/exclusion criteria ------------------------------------

# Cohort for main analysis
nsfg_cohort <- nsfg %>%
  filter(
    (OUTCOME != 2 & OUTCOME != 6), # Pregnancy outcome is not an induced abortion or a current pregnancy
    (survey %in% c("2006-2010", "2011-2013", "2013-2015") & AGEPREG >= 1450 & AGEPREG <= 4499) | (survey %in% c("2015-2017", "2017-2019") & AGEPREG >= 15 & AGEPREG <= 44), # Age at pregnancy outcome is between 15 and 44
    (survey %in% c("2006-2010", "2011-2013", "2013-2015") & DATECON >= 1194 & DATECON <= 1388) | (survey %in% c("2015-2017", "2017-2019") & DATECON >= 2000 & DATECON <= 2018)) # Year of conception is between 2000 and 2018


# Cohorts for sensitivity analyses
# Restricted cohort to pregnancies that were intended or mistimed at conception 
# nsfg_cohort <- nsfg_cohort %>%
#   filter(WANTRESP != 5) # Pregnancy was not unwanted at conception

# Expanded cohort to include all completed pregnancies
# nsfg_cohort <- nsfg %>%
#   filter(OUTCOME != 6, # Pregnancy outcome is not a current pregnancy
#         (survey %in% c("2006-2010", "2011-2013", "2013-2015") & AGEPREG >= 1450 & AGEPREG <= 4499) | (survey %in% c("2015-2017", "2017-2019") & AGEPREG >= 15 & AGEPREG <= 44),
#         (survey %in% c("2006-2010", "2011-2013", "2013-2015") & DATECON >= 1194 & DATECON <= 1388) | (survey %in% c("2015-2017", "2017-2019") & DATECON >= 2000 & DATECON <= 2018))

# Recode and label variables ---------------------------------------------------

nsfg_cohort <- nsfg_cohort %>%
  mutate(
    # Pregnancy loss
    loss = ifelse(OUTCOME %in% c(3, 4, 5), 1, 0), # Stillbirth, miscarriage, and ectopic pregnancy
    
    # Conception year interval
    year_int = cut(year_con, 
                   breaks = c(2000, 2002, seq(2004, 2018, by = 2), 2020), 
                   include.lowest = TRUE),
    
    year_int = case_when(
      as.character(year_int) == "[2000,2002]" ~ "(1999,2002]", # Recode first conception year interval for formatting consistency across all intervals
      TRUE ~ as.character(year_int)),

    # Midpoint year of each conception year interval
    year = as.numeric(substr(year_int, 2, 5)) + 1,
    
    # Metro status at interview
    metro_int = ifelse(METRO <= 2, 1, 0),
    
    # Maternal age group at conception
    agegrp_con = case_when(
      age_con >= 15 & age_con <= 19 ~ 1,
      age_con >= 20 & age_con <= 24 ~ 2,
      age_con >= 25 & age_con <= 29 ~ 3,
      age_con >= 30 & age_con <= 34 ~ 4,
      age_con >= 35 & age_con <= 44 ~ 5,
      TRUE ~ NA_real_),
    
    # Hispanic origin and race
    raceth = HISPRACE2,
    
    # Educational attainment at interview
    edu_int = case_when(
      HIEDUC %in% c(5, 6, 7, 8) ~ 1,
      HIEDUC == 9 ~ 2,
      HIEDUC %in% c(10, 11) ~ 3,
      HIEDUC == 12 ~ 4,
      HIEDUC %in% c(13, 14, 15) ~ 5,
      TRUE ~ NA_real_),
    
    # Household income as a percentage of poverty level at interview
    poverty_int = case_when(
      POVERTY <= 99 ~ 1,
      POVERTY >= 100 & POVERTY <= 199 ~ 2,
      POVERTY >= 200 & POVERTY <= 299 ~ 3,
      POVERTY >= 300 & POVERTY <= 399 ~ 4,
      POVERTY > 399 ~ 5,
      TRUE ~ NA_real_),
    
    # Marital status at pregnancy end
    mar_preg = case_when(
      FMAROUT5 == 1 ~ 1,
      FMAROUT5 %in% c(2, 3, 4) ~ 2,
      FMAROUT5 == 5 ~ 3,
      TRUE ~ NA_real_),
    
    # Intendedness of pregnancy at conception
    intend_con = case_when(
      WANTRESP %in% c(1, 2, 4, 6) ~ 1,
      WANTRESP == 3 ~ 2,
      WANTRESP == 5 ~ 3,
      TRUE ~ NA_real_),
    
    # Gravidity
    firstpreg = ifelse(PREGORDR == 1, 1, 0),
    
    # Parity
    parity = case_when(
      PARITY == 0 ~ 0,
      PARITY == 1 ~ 1,
      PARITY == 2 ~ 2,
      PARITY > 2 ~ 3,
      TRUE ~ NA_real_))

# Subset to needed variables
nsfg_cohort <- nsfg_cohort %>%
  select(CASEID, stratvar, panelvar, weightvar, survey,
         loss, OUTCOME,
         year_con, year_int, year, DATECON,
         agegrp_con, AGECON,
         metro_int, METRO,
         raceth, HISPRACE2,
         edu_int, HIEDUC,
         poverty_int, POVERTY,
         mar_preg, FMAROUT5,
         intend_con, WANTRESP,
         firstpreg, PREGORDR,
         parity, PARITY)

# Assign labels to categorical variables
nsfg_cohort_labelled <- nsfg_cohort %>%
  mutate(
    metro_int = factor(metro_int, levels = c(0, 1), labels = c("Nonmetropolitan", "Metropolitan")),
    agegrp_con = factor(agegrp_con, levels = 1:5, labels = c("15-19", "20-24", "25-29", "30-34", "35-44")),
    raceth = factor(raceth, levels = 1:4, labels = c("Hispanic or Latina", "Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Other")),
    edu_int = factor(edu_int, levels = 1:5, labels = c("No high school diploma or GED", "High school diploma or GED", "Some college, no bachelor's degree", "Bachelor's degree", "Master's degree or higher")),
    poverty_int = factor(poverty_int, levels = 1:5, labels = c("Less than 100%", "100% - 199%", "200% - 299%", "300% - 399%", "400% or more")),
    mar_preg = factor(mar_preg, levels = 1:3, labels = c("Married", "Widowed, divorced, separated", "Never married")),
    intend_con = factor(intend_con, levels = 1:3, labels = c("Intended", "Mistimed", "Unwanted")),
    firstpreg = factor(firstpreg, levels = c(0, 1), labels = c("Not first pregnancy", "First pregnancy")),
    parity = factor(parity, levels = 0:3, labels = c("No children", "One child", "Two children", "Three or more children")))

# str(nsfg_cohort)
# str(nsfg_cohort_labelled)

# ------------------------------------------------------------------------------
# Create survey design object
# ------------------------------------------------------------------------------

nsfg_design <- svydesign(id = ~panelvar, 
                 weights = ~weightvar, 
                 strata = ~stratvar, 
                 nest = T, 
                 survey.lonely.psu = "adjust", 
                 data = nsfg_cohort_labelled)

# summary(nsfg_design)

# ------------------------------------------------------------------------------
# Custom function to calculate adjusted confidence intervals
#  -----------------------------------------------------------------------------

# Calculates the 95% confidence interval (CI) of a proportion using replicate weights in a complex survey design
# Modifies the Lumley Korn-Graubard method in svyciprop by incorporating adjustments for effective sample size and degrees of freedom (DF)

# Parameters:
# - formula: A formula specifying the outcome variable and predictor variables
# - design: A survey design object (created with svydesign)
# - method: The method to use for CI calculation
# - level: The confidence level for the interval (default is 0.95 for a 95% CI)
# - df: Degrees of freedom for the design

# Returns: A svyciprop object containing the estimated proportion and CIs

svyciprop.adj = function (formula, design, method = c("logit", "likelihood", "asin", "beta", "mean"), level = 0.95, df = degf(design), ...) {
  
  method = match.arg(method)
  
  # If the beta method is selected, perform specific calculations
  if (method == "beta") {
    
    m = eval(bquote(svymean(~as.numeric(.(formula[[2]])), design, ...))) 
    rval = coef(m)[1]
    n.eff = coef(m) * (1 - coef(m)) / vcov(m)
    attr(rval, "var") = vcov(m)
    alpha = 1 - level
    df = nrow(design) - 1
  
    if (df > 0) { 
      rat.squ = (qt(alpha / 2, nrow(design) - 1) / qt(alpha / 2, df))^2 
    } else {
      rat.squ = 0
    }
    
    if (rval > 0) { 
      n.eff = min(nrow(design), n.eff * rat.squ)
    } else {
      n.eff = nrow(design)
    }

    ci = c(qbeta(alpha / 2, n.eff * rval, n.eff * (1 - rval) + 1), 
           qbeta(1 - alpha / 2, n.eff * rval + 1, n.eff * (1 - rval)))
    
  # For other methods, call the standard svyciprop function
  } else {
    ci = svyciprop(formula, design, method, level, df, ...)
  }
  
  halfalpha = (1 - level) / 2
  names(ci) = paste(round(c(halfalpha, (1 - halfalpha)) * 100, 1), "%", sep = "")
  names(rval) = deparse(formula[[2]])
  attr(rval, "ci") = ci
  class(rval) = "svyciprop"
  
  rval # Return estimated proportion and CIs
}

# ------------------------------------------------------------------------------
# Calculate direct estimates (inputs for eMKF)
# ------------------------------------------------------------------------------

# Calculate the proportion of pregnancies ending in loss for each subgroup
# Subgroups are defined by conception year interval, metro status, and maternal age group
emkf_input <- svyby(~loss, ~year + metro_int + agegrp_con, nsfg_design, na.rm = T, deff = T, svymean)

# Calculate 95% CIs using the svyciprop.adj function for each subgroup
ci <- svyby(~I(loss == 1), ~year + metro_int + agegrp_con, nsfg_design, na.rm = T, method = "beta", svyciprop.adj, vartype = "ci")

# Add columns to emkf_input dataframe
emkf_input$ci_lower <- ci[,ncol(ci) - 1]
emkf_input$ci_upper <- ci[,ncol(ci)]
emkf_input$neff <- ((emkf_input$loss * (1 - emkf_input$loss)) / (emkf_input$se^2)) # Effective sample size = (p(1-p)/se^2); increases as the standard error decreases, and vice versa

# Convert survey design object back to a dataframe
nsfg_df <- as.data.frame(nsfg_design$variables) %>%
  group_by(agegrp_con, metro_int, year) %>%
  summarise(n = n())

# Merge dataframes by subgroup variables
emkf_input <- merge(emkf_input, nsfg_df, by = c("agegrp_con", "metro_int", "year"))

# Save eMKF input dataset
# write.csv(emkf_input,"...\eMKF\MKFdata\nsfg_2006-2019_pregnancy-loss.csv", row.names = F, na = "")

# ------------------------------------------------------------------------------
# Calculate model-based estimates (outputs from eMKF)
# ------------------------------------------------------------------------------

# Run the Example-eMKF-NSFG-Data.sas file (using the eMKF input dataset created above)
# After running, return to this file and read in the eMKF output dataset to complete the analysis

# Read in eMKF output dataset
emkf_output <- read.csv("...\eMKF\MKFdata\bmac_pred_nsfg_2006-2019_pregnancy-loss.csv")

# Make subgroup variables factors
emkf_output <- emkf_output %>%
  mutate(
    metro_int = as.factor(metro_int),
    agegrp_con = as.factor(agegrp_con))
  
# Combine eMKF input and output dataframes
emkf_combined <- merge(emkf_input, emkf_output, by = c("agegrp_con", "metro_int", "year")) 

# Remove duplicate estimates for loss, se, and neff in combined dataframe
emkf_combined <- emkf_combined %>%
  select(-loss.y, -se.y, -neff.y) %>% # Remove variables originally from emkf_output (model-based estimates)
  rename(loss = loss.x, # Rename variables originally from emkf_input (direct estimates)
         se = se.x,
         neff = neff.x)

# ------------------------------------------------------------------------------
# Estimate trends using weighted least squares log-binomial regression
# ------------------------------------------------------------------------------

# Create a weight variable to give more weight to observations with higher precision and less weight to observations with lower precision
emkf_combined <- emkf_combined %>%
  mutate(weight = 1 / predVar_Bayes_BMA_CUBIC) # Weight = inverse of the predicted variance from Bayesian model averaging

# Fit a weighted least squares log-binomial regression model
wls <- emkf_combined %>%
  group_by(agegrp_con, metro_int) %>%
  do(model = Glm(pred_Bayes_BMA_CUBIC ~ year * metro_int, data = . , weights = weight, family = binomial())) 

# Print model output for each subgroup by maternal age and metro status
foreach(x=1:dim(wls)[1]) %do% c(paste0(wls$agegrp_con[x], as.character(wls$metro_int[x])), print(wls$model[x]))

# ------------------------------------------------------------------------------
# Calculate improvements in precision
# ------------------------------------------------------------------------------

# Create CI variables using model-based estimates (eMKF output)
emkf_combined <- emkf_combined %>%
  mutate(year = as.numeric(year),
         pred_ci_lower = ifelse((pred_Bayes_BMA_CUBIC - 1.96 * predSE_Bayes_BMA_CUBIC) < 0, 0, pred_Bayes_BMA_CUBIC - 1.96 * predSE_Bayes_BMA_CUBIC), # Bounded lower CI = (p_hat - 1.96 * se(p_hat)), where p_hat = pred_Bayes_BMA_CUBIC and lower bound not < 1
         pred_ci_upper = ifelse((pred_Bayes_BMA_CUBIC + 1.96 * predSE_Bayes_BMA_CUBIC) > 1, 1, pred_Bayes_BMA_CUBIC + 1.96 * predSE_Bayes_BMA_CUBIC)) # Bounded upper CI = (p_hat + 1.96 * se(p_hat)), where p_hat = pred_Bayes_BMA_CUBIC and upper bound not > 1

# Calculate relative 95% CI widths and precision improvements
emkf_combined <- emkf_combined %>%
  mutate(rel_ci_model = 100 *(pred_ci_upper-pred_ci_lower) / pred_Bayes_BMA_CUBIC, # Relative width of the 95% CI for model-based estimates, expressed as a percentage
         rel_ci_direct = 100 * (ci_upper-ci_lower) / loss, # Relative width of the 95% CI for direct estimates, expressed as a percentage
         rel_ci_imp = 1 - (rel_ci_model / rel_ci_direct)) # Relative improvement in precision of the 95% CI for the model-based estimate compared to the 95% CI for the direct estimate, calculated as one minus the ratio of the relative widths

# Summarize median values for the relative 95% CI widths and precision improvements
emkf_combined %>%
  group_by(metro_int) %>% # By metro status only
  summarise(rel_ci_model = median(rel_ci_model),
            rel_ci_direct = median(rel_ci_direct),
            rel_ci_imp = median(rel_ci_imp))

emkf_combined %>%
  group_by(metro_int, agegrp_con) %>% # By metro status and maternal age group
  summarise(rel_ci_model = median(rel_ci_model),
            rel_ci_direct = median(rel_ci_direct),
            rel_ci_imp = median(rel_ci_imp))

# ------------------------------------------------------------------------------
# Calculate age-adjusted estimates for all-ages (standardized to the overall distribution of age groups)
# ------------------------------------------------------------------------------

# Calculate weighted frequencies for each subgroup using survey-weighted population size
weighted_freqs <- as.data.frame(svytable(~year + metro_int + agegrp_con, design = nsfg_design)) %>%
  rename(weighted_freq = Freq)

emkf_combined_adjusted <- merge(emkf_combined, weighted_freqs, by = c("year", "metro_int", "agegrp_con")) %>%
  group_by(metro_int, year) %>% 
  mutate(weighted_freq_sum = sum(weighted_freq), # Sum of weighted Ns for each subgroup by metro status and conception year interval
         freq_sum = sum(n)) %>% # Sum of unweighted Ns for each subgroup by metro status and conception year interval
  ungroup()

# Calculate age-specific weighted estimates for each subgroup
emkf_combined_adjusted <- emkf_combined_adjusted %>%
  mutate(weight = weighted_freq / weighted_freq_sum, # 0 < weight < 1; subgroup weights sum to 1
         weighted_est = pred_Bayes_BMA_CUBIC * weight)

# Calculate age-adjusted all-ages estimate using the weighted average of each age-specific weighted estimate for each subgroup
emkf_combined_adjusted <- emkf_combined_adjusted %>%
  group_by(metro_int, year) %>%
  mutate(
    age_adjusted_est = sum(weight * pred_Bayes_BMA_CUBIC)) %>% # Age-adjusted all-ages estimate
  ungroup()

# Calculate the variability of each age-adjusted all-ages estimate
emkf_combined_adjusted <- emkf_combined_adjusted %>%
  group_by(metro_int, year) %>%
  mutate(age_adjusted_var = sum(weight^2 * predVar_Bayes_BMA_CUBIC),
         age_adjusted_se = sqrt(age_adjusted_var)) %>%
  ungroup()

# Collapse dataframe for all ages
emkf_combined_adjusted %>%
  distinct(metro_int, year, age_adjusted_est, age_adjusted_var, age_adjusted_se, .keep_all = FALSE)

# Estimate trends for all-ages using log-binomial regression models ------------

emkf_combined_adjusted <- emkf_combined_adjusted %>%
  mutate(weight = 1/age_adjusted_var)

wls_adjusted <- emkf_combined_adjusted %>%
  group_by(metro_int) %>%
  do(model = Glm(pred_Bayes_BMA_CUBIC ~ year * metro_int, data = . , weights = weight, family = binomial())) 

foreach(x=1:dim(wls_adjusted)[1]) %do% c(as.character(wls_adjusted$metro_int[x]), print(wls_adjusted$model[x]))

# ------------------------------------------------------------------------------
# Plot direct and model-based estimates
#-------------------------------------------------------------------------------

color_palette <- c("Metropolitan" = "#0033A1", "Nonmetropolitan" = "#d06f1a") 

theme_elements <- theme(
  axis.line.x = element_line(linewidth = 0.3),
  axis.line.y = element_line(linewidth = 0.3),
  axis.text.x = element_text(angle = 90, hjust = 1),
  text = element_text(family = "Calibri", size = 12),
  strip.background = element_blank(),
  axis.line = element_line(),
  legend.position = "none",
  axis.title.x = element_text(margin = margin(t = 5)),
  axis.title.y = element_text(margin = margin(r = 5)))

# Plot for direct estimates only
plot_direct <- ggplot(emkf_combined) +
  geom_ribbon(aes(x = year, ymin = ci_lower, ymax = ci_upper, fill = as.factor(metro_int)), alpha = 0.2) +
  geom_smooth(aes(x = year, y = loss, group = as.factor(metro_int), color = as.factor(metro_int), linetype = as.factor(metro_int)), alpha = 1, linewidth = 0.75, method = "glm", se = FALSE) +
  geom_point(aes(x = year, y = loss, group = as.factor(metro_int), shape = as.factor(metro_int), color = as.factor(metro_int)), alpha = 1, size = 1.75) +
  theme_classic() + 
  theme_elements +
  scale_y_continuous("Risk of pregnancy loss", labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(name = "Conception year", breaks = seq(2000, 2018, by = 2), limits = c(2000, 2018)) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~agegrp_con, scales = "free_x", ncol = 5) +
  labs(title = "Direct estimates")

# Plot for model-based estimates only
plot_model <- ggplot(emkf_combined) +
  geom_ribbon(aes(x = year, ymin = pred_ci_lower, ymax = pred_ci_upper, fill = as.factor(metro_int)), alpha = 0.2) +
  geom_smooth(aes(x = year, y = pred_Bayes_BMA_CUBIC, group = as.factor(metro_int), color = as.factor(metro_int), linetype = as.factor(metro_int)), alpha = 1, linewidth = 0.75, method = "glm", se = FALSE) +
  geom_point(aes(x = year, y = pred_Bayes_BMA_CUBIC, group = as.factor(metro_int), shape = as.factor(metro_int), color = as.factor(metro_int)), alpha = 1, size = 1.75) +
  theme_classic() + 
  theme_elements +
  scale_y_continuous("Risk of pregnancy loss", labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(name = "Conception year", breaks = seq(2000, 2018, by = 2), limits = c(2000, 2018)) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~agegrp_con, scales = "free_x", ncol = 5) +
  labs(title = "Model-based estimates")

# Combined plot for both direct and model-based estimates
plot_combined <- (plot_direct + plot_model) + 
  plot_layout(nrow = 2)

# Save plot
# ggsave("plot.tiff", plot = final_plot, width = 9, height = 6, dpi = 600, compression = "lzw")