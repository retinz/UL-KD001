# DESCRIPTION 
# ===================== #

# Preprocessing questionnaire data for UL-KD001, starting with already
# pseudo-anonymised data.

# TASKS 
# ===================== #

# - Add the calculation of the no. of participants who finished pretest
# - Adjust the loading of the RDS from this script to taskswitch_prep

# - Check if the above has been done
# - Adjust data loading (using here())

# NOTES 
# ===================== #

# - Only the KD April 2024 cohort has ketone values from before Pretest. There are
# only 2 missing values for 2 different individuals --> no need to filter based 
# on missingness. But the missing values are imputed. 

# - CD have some additional columns. Some of the additional columns were assessed 
# but eventually removed for convenience. They also have only 15 (intervention / posttest) ketone
# measurements by default. CD Pretest is missing ketone measurements.

# - Exclusions in terms of age, BMI, sex, finishing, measuring ketones, previous diet, 
# and ketone missingness are performed at the cohort level. Exclusions in terms of ketone values
# themselves are performed at the study / final level. Ketone imputation is performed 
# at the cohort level.

# - In no cohort is there a case of both the 15th and 16th measurement missing 
# at the same time (same participant) --> no need for a special check of the 15th measurement
# with regard to posttest_date to see if imputation could be performed there. 

# - KD April 2024 has slightly different rules for filtering based on ketone 
# missingness and for imputation than the other cohorts. This is due to 
# a specific date for Posttest in the case of this cohort. 

# LIBRARIES #
# ===================== #

library(here)
library(zoo)
library(readxl)
library(writexl)
library(ggplot2)
library(tidyverse)
source(here('UL-KD001_analysis', 'analysis_helpers.R'))

# PATHS #
# ===================== #

# Data path
data_path <- here('UL-KD001_data', 'LBI_raw', 'raw_pseudoanonym')

# Plots path
plot_directory <- here('UL-KD001_analysis', 'plots')
if (!dir.exists(plot_directory)) {
  dir.create(plot_directory, recursive = TRUE)
}

# Preprocessing directory
prep_directory <- here('UL-KD001_preprocessing')

# LOAD PSEUDO-ANONYMISED DATA ----------------------------------

# KD Apr 2024 #
# --------------- #

kd_apr_24_pre <- read.csv(file.path(data_path, 'kd_apr_24_pre.csv'))
kd_apr_24_post <- read.csv(file.path(data_path,'kd_apr_24_post.csv')) %>% 
  relocate(participant_id, EndDate, Finished, Q75, Q74.1)

# KD Sep 2024 #
# --------------- #

kd_sep_24_pre <- read.csv(file.path(data_path,'kd_sep_24_pre.csv')) %>%
  relocate(participant_id, EndDate, Finished, Q75, Q74.1)
kd_sep_24_post <- read.csv(file.path(data_path,'kd_sep_24_post.csv')) %>%
  relocate(participant_id, Finished, EndDate, Q75, Q74.1)

# KD Jan 2025 #
# --------------- #

kd_jan_25_pre <- read.csv(file.path(data_path,'kd_jan_25_pre.csv')) %>%
  relocate(participant_id, EndDate, Finished, Q75, Q74.1)
kd_jan_25_post <- read.csv(file.path(data_path,'kd_jan_25_post.csv')) %>%
  relocate(participant_id, Finished, EndDate, Q75, Q74.1)

# CD 2025 #
# -------------- #

cd_25_pre <- read.csv(file.path(data_path,'cd_25_pre.csv')) %>%
  relocate(participant_id, EndDate, Finished, Q75, Q74.1)
cd_25_post <- read.csv(file.path(data_path,'cd_25_post.csv')) %>%
  relocate(participant_id, Finished, EndDate, Q74.1)
cd_25_post_ketones <- read.csv(file.path(data_path,'cd_25_post_ketones.csv')) %>%
  relocate(participant_id, Finished, Q2)

# CODEBOOKS ----------------------------

# Function to create a codebook from the first two rows of a dataframe
make_BOOK <- function(df) {
  df %>%
    slice_head(n = 2) %>%   
    mutate(row = row_number()) %>%
    pivot_longer(-row, names_to = 'var', values_to = 'val') %>%
    pivot_wider(names_from = row, values_from = val)
}

data_objs <- c(
  'kd_apr_24_pre',  'kd_apr_24_post',
  'kd_sep_24_pre',  'kd_sep_24_post',
  'kd_jan_25_pre',  'kd_jan_25_post',
  'cd_25_pre',      'cd_25_post', 
  'cd_25_post_ketones'
)

books <- data_objs %>%
  set_names() %>%
  purrr::map(\(nm) make_BOOK(get(nm))) %>%
  set_names(paste0(data_objs, '_BOOK'))

# Compare columns #
# ------------------- #

# Function for comparing pairs
diff_vars <- function(a, b) {
  list(
    only_in_a = setdiff(a$var, b$var),
    only_in_b = setdiff(b$var, a$var)
  )
}

# Define pairs
pairs <- tribble(
  ~a,               ~b,
  'kd_apr_24_pre_BOOK',  'kd_sep_24_pre_BOOK',
  'kd_apr_24_post_BOOK', 'kd_sep_24_post_BOOK',
  'kd_jan_25_pre_BOOK',  'kd_sep_24_pre_BOOK',
  'kd_jan_25_post_BOOK', 'kd_sep_24_post_BOOK',
  'kd_jan_25_pre_BOOK',  'cd_25_pre_BOOK',
  'kd_jan_25_post_BOOK', 'cd_25_post_BOOK'
)

# Compare desired pairs
var_diffs <- pairs %>%
  mutate(diff = map2(a, b, ~ diff_vars(books[[.x]], books[[.y]])))

# Save the codebooks to Excel
imap(books,
     ~ write_xlsx(.x, file.path(prep_directory, paste0(.y, '.xlsx'))))

# PREPROCESSING KD APRIL 2024 ---------------------------------
# - Pretest cleaning ---------------------------------

kd_apr_24_pre_CLEAN <- kd_apr_24_pre %>%
  
  # Remove unnecessary rows (contain metadata)
  slice(-c(1, 2)) %>%
  
  # Convert key columns
  mutate(
    weight = as.numeric(Q73.1),  # Ensure weight is numeric
    height = as.numeric(Q74.1) / 100,  # Convert height from cm to meters
    
    # Correct weight values --> weight in KG
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ),
    
    # Clean up sex column
    sex = case_when(
      Q75 == 'Vrouw (XX)' ~ 'female',
      Q75 == '' ~ NA_character_,  # Explicitly setting NA for missing values
      .default = NA_character_
    )
  ) %>%
  
  # Remove unnecessary columns
  select(-Q72, -Q73, -Q74, -RecordedDate, -DistributionChannel, -UserLanguage, 
         -Duration..in.seconds., -Status, -Progress, -ResponseId, -Q75) %>%
  
  # Dynamically rename ketone recording columns (Q76_#_1 → ketones_#)
  rename_with(~ str_replace(.x, '^Q76_(\\d+)_1$', 'ketones_\\1'), starts_with('Q76_')) %>%
  
  # Rename other necessary columns
  rename(
    age = Q72.1, measured_ketones = Q75.1,
    current_diet = Q77, current_diet_text = Q77_6_TEXT, 
    start_date = StartDate, end_date = EndDate, 
    finished = Finished
  ) %>%
  
  # Add computed columns
  mutate(
    BMI = weight / (height^2),
    
    # Convert 'measured_ketones' from text to logical
    measured_ketones = case_when(
      measured_ketones == 'Ik heb mijn ketowaarden niet gemeten' ~ FALSE,
      measured_ketones == 'Ik heb mijn ketowaarden gemeten' ~ TRUE,
      .default = NA  # Ensure it remains logical
    ),
    
    # Make age numeric
    age = as.numeric(age),
    
    # Convert 'finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  
  # Relocation
  relocate(participant_id, sex, weight, height, BMI) %>%
  # Final adjustments
  rename(pretest_date = end_date,
         pretest_time = end_time) %>%
  mutate(pretest_date = as.Date(pretest_date)) %>%  # Convert to short date format
  select(-start_date)

# Investigating duplicates #
# --------------------------------- #

duplicates_apr24_pre <- kd_apr_24_pre_CLEAN %>%
  filter(finished, measured_ketones) %>% # Select only cases where column == TRUE
  filter(BMI >= 18.5 & BMI < 30,
         sex == 'female',
         age >= 35 & age <= 45
  ) %>%
  # Get only duplicates again after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(pretest_date, pretest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# ------------------------------- #

kd_apr_24_pre_FRESH <- kd_apr_24_pre_CLEAN %>% 
  setdiff(duplicates_apr24_pre)  # Exclude duplicates from pretest

# - Posttest cleaning ----------------------------------

kd_apr_24_post_CLEAN <- kd_apr_24_post %>%
  select(-Status, -Progress, -Duration..in.seconds., -RecordedDate, 
         -ResponseId, -DistributionChannel, -UserLanguage, -Q72, -Q73, -Q74) %>%
  relocate(participant_id) %>%
  rename_with(~ str_replace(.x, '^Q76_(\\d+)_1$', 'ketones_\\1'), starts_with('Q76_')) %>%  # Dynamic renaming
  rename(
    weight = Q73.1, height = Q74.1, measured_ketones = Q75, 
    start_date = StartDate, end_date = EndDate, finished = Finished
  ) %>%
  # Remove unnecessary rows
  slice(-c(1, 2)) %>%
  mutate(
    weight = as.numeric(weight),  # Ensure weight is numeric
    height = as.numeric(height) / 100,  # Convert height from cm to meters
    
    # Correct weight values
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ), 
    
    # Adding BMI
    BMI = weight / (height^2),
    
    # Convert 'measured_ketones' from text to logical
    measured_ketones = case_when(
      measured_ketones == 'Ik heb mijn ketowaarden niet gemeten' ~ FALSE,
      measured_ketones == 'Ik heb mijn ketowaarden gemeten' ~ TRUE,
      .default = NA  # Ensure it remains logical
    ),
    
    # Convert 'Finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  # Final adjustments
  rename(posttest_date = end_date,
         posttest_time = end_time) %>%
  mutate(posttest_date = as.Date(posttest_date)) %>%
  select(-start_date)

# Investigating duplicates #
# --------------------------------- #

duplicates_apr24_post <- kd_apr_24_post_CLEAN %>%
  filter(finished, measured_ketones) %>% # Select only cases where column == TRUE
  # Get only duplicates after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(posttest_date, posttest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# -------------------------------- #

kd_apr_24_post_FRESH <- kd_apr_24_post_CLEAN %>% 
  setdiff(duplicates_apr24_post)  # Exclude duplicates from posttest

# - Combine datasets ------------------------------

# COMBINED DATA FRAME: Merging data sets
kd_apr_24_combined <- right_join(
  kd_apr_24_pre_FRESH, kd_apr_24_post_FRESH, 
  by = 'participant_id', 
  suffix = c('_pre', '_post')  # Automatically renames overlapping columns
) %>%
  # The following columns don't overlap --> need to add the suffix in this way
  rename_with(
    ~ gsub('ketones_(9|10|11|12|13|14|15|16)$', 'ketones_\\1_post', .),
    matches('^ketones_(9|10|11|12|13|14|15|16)$')
  ) %>%
  mutate(across(matches('^ketones_\\d+_(pre|post)$'),
                ~ {
                  ketone <- parse_number(str_replace_all(.x, ',', '.'))
                  round(
                    ifelse(ketone >= 10,
                           ketone / 10^(nchar(trunc(ketone)) - 1),
                           ketone),
                    1
                  )
                })) %>%
  mutate(current_diet = str_replace_all(current_diet, ' ', '_')) %>%
  
  # Adding metadata
  mutate(age = as.numeric(age), 
         group = 'KD',
         cohort = 'apr_24')
glimpse(kd_apr_24_combined)

# TIBBLE DESCRITION:
# - weight ~ in kg
# - height ~ in metres
# - ketones ~ in mmol/L

# FILTERING KD APRIL 2024 & IMPUTATION ---------------------------
  
# FILTERED DATASET (basic exclusions made) #
# ============================================== #

# All empty IDs from previous steps get filtered out due to one of the below filters; 
# no specific filter needed. Test IDs filtered out in an anonymiser script (not
# included here).

kd_apr_24_CLEAN <- kd_apr_24_combined %>%
  filter(finished_pre, finished_post,
         measured_ketones_pre, measured_ketones_post) %>% # Select only cases where column == TRUE
  filter(BMI_pre >= 18.5 & BMI_pre < 30,
         sex == 'female',
         age >= 35 & age <= 45,
  ) %>%
  # Exclude Paleo and Low carb diets
  filter(!current_diet %in% c('Paleo_dieet', 'Low_carb_dieet')) %>%
  arrange(participant_id)

# - Ketone missingness checks -------------------------------

# Check the nature of ketone missingness at Pretest
missing_pre_ketones_apr24 <- kd_apr_24_CLEAN %>%
  dplyr::select(participant_id, matches('^ketones_\\d+_pre$')) %>%
  pivot_longer(cols = matches('^ketones_\\d+_pre$'), 
               names_to = 'ketone_measurement', 
               values_to = 'ketone_level') %>%
  filter(is.na(ketone_level)) %>%
  group_by(participant_id) %>%
  summarise(missing_count = n(), .groups = 'drop') %>%
  arrange(desc(missing_count))
missing_pre_ketones_apr24

# Check the nature of ketone missingness at Posttest
missing_post_ketones_apr24 <- kd_apr_24_CLEAN %>%
  dplyr::select(participant_id, matches('^ketones_\\d+_post$')) %>%
  pivot_longer(cols = matches('^ketones_\\d+_post$'), 
               names_to = 'ketone_measurement', 
               values_to = 'ketone_level') %>%
  filter(is.na(ketone_level)) %>%
  group_by(participant_id) %>%
  summarise(missing_count = n(), .groups = 'drop') %>%
  arrange(desc(missing_count))
missing_post_ketones_apr24

# - Imputing ketones PRETEST ----------------------------------

# _int_ ~ interpolated

# 1. Identify ketone columns
ketone_apr24_cols_pre <- grep('^ketones_\\d+_pre$', names(kd_apr_24_CLEAN), value = TRUE)

ketones_apr24_imputed_pre <- kd_apr_24_CLEAN %>%
  
  # Make sure columns are numeric
  mutate(across(all_of(ketone_apr24_cols_pre), as.numeric)) %>%
  
  # 2. Go row-by-row
  rowwise() %>%
  
  # 3. For each row, create a single list-column with the interpolated vector
  mutate(
    ketones_pre_int = list({
      row_vals <- c_across(all_of(ketone_apr24_cols_pre))
      
      # If they're all NA, return them as-is instead of NULL
      if (all(is.na(row_vals))) {
        
        row_vals
        
      } else {
        
        # Interpolate missing values across that row of ketones
        ## Rule = 2 ~ interpolate both leading and trailing NAs
        ## Rule = 1 ~ don't interpolate either leading or trailing NAs
        as.numeric(na.approx(row_vals, na.rm = FALSE, rule = 2))
      }
    })
  ) %>%
  ungroup() %>%
  
  # 4. 'Unnest' that list of interpolated values back into columns
  unnest_wider(ketones_pre_int, names_sep = '_') %>%
  rename(ketones_pre = ketones_pre_int_8)

# - Imputing ketones POST --------------------------------------

# POST ~ intervention ketones
# _int_ ~ interpolated

# Filter:
# - Measurement 16 can be missing only if posttest_date <= 2024-05-06 [that's one day after 
# measurement 16]
# - No 2 adjacent measurements can be missing
# - Max. 3 / 8 measurements can be missing (from the last 2 weeks; the first 2 weeks not important)
# -- Max. 2 / 7 measurements can be missing if measurement 16 is missing (and if that's allowed)

# Imputation:
# - Linear interpolation in rows over max. 1 cell
# - Trailing NAs interpolated only when posttest_date > 2024-05-04 & posttest_date < 2024-05-07 (5 or 6 May)
# - Leading NAs interpolated

# Identify ALL posttest columns and sort them numerically
ketone_apr24_cols_post <- grep('^ketones_\\d+_post$', names(kd_apr_24_CLEAN), value = TRUE)
ketone_apr24_cols_post <- ketone_apr24_cols_post[order(as.numeric(gsub('.*_(\\d+)_post', '\\1', ketone_apr24_cols_post)))]

# Identify ONLY columns 9..16 for filtering and sort them numerically
filter_cols_9_16 <- grep('^ketones_(9|1[0-6])_post$', names(kd_apr_24_CLEAN), value = TRUE)
filter_cols_9_16 <- filter_cols_9_16[order(as.numeric(gsub('.*_(\\d+)_post', '\\1', filter_cols_9_16)))]

ketones_apr24_imputed_post <- kd_apr_24_CLEAN %>%
  # Drop baseline ketones columns
  select(-all_of(ketone_apr24_cols_pre)) %>%
  
  # Ensure all posttest columns are numeric
  mutate(across(all_of(ketone_apr24_cols_post), as.numeric)) %>%
  
  # Process each row individually
  rowwise() %>%
  
  # * Filter Step * - only for columns 9..16
  filter({
    row_vals_9_16 <- c_across(all_of(filter_cols_9_16))
    
    # Helper function: Are any two consecutive NA?
    is_consecutive_na <- function(x) {
      for (i in seq_along(x)[-1]) {
        if (is.na(x[i]) && is.na(x[i - 1])) {
          return(TRUE)
        }
      }
      return(FALSE)
    }
    
    # Count how many NAs
    n_missing <- sum(is.na(row_vals_9_16))
    
    # Check if measurement 16 is missing (position 8 in the vector)
    vals <- as.vector(row_vals_9_16)
    measurement_16_missing <- length(vals) >= 8 && is.na(vals[8])
    
    # Allow missing 16th measurement only if posttest_date <= 2024-05-06
    allow_missing_16 <- FALSE
    if (measurement_16_missing) {
      # Parse the posttest_date as Date if it's character
      if (is.character(posttest_date)) {
        date_to_check <- as.Date(posttest_date)
      } else {
        date_to_check <- posttest_date
      }
      
      # Check if the date condition is met
      allow_missing_16 <- date_to_check <= as.Date('2024-05-06')
      
      # If missing 16 is allowed due to date, don't count it as missing
      if (allow_missing_16) {
        n_missing <- n_missing - 1
      }
    }
    
    # Keep row only if:
    # - No consecutive NA in 9..16
    # - At most 3 total missing in 9..16
    # - Measurement 16 can be NA only if posttest_date <= 2024-05-06
    
    # Case when measurement 16 is missing
    if (measurement_16_missing) {
    !is_consecutive_na(row_vals_9_16) && 
      (n_missing <= 2) && 
      (allow_missing_16)} 
    # Case when measurement 16 is not missing
    else {
        !is_consecutive_na(row_vals_9_16) && 
          (n_missing <= 3)}
  }) %>%
  
  # Capture the date-based condition for use in interpolation decision
  mutate(
    # Store whether date condition is met for each row - between 2024-05-04 and 2024-05-07 (exclusive)
    date_in_range = if_else(
      is.character(posttest_date), 
      as.Date(posttest_date) > as.Date('2024-05-04') & as.Date(posttest_date) < as.Date('2024-05-07'),
      posttest_date > as.Date('2024-05-04') & posttest_date < as.Date('2024-05-07')
    )
  ) %>%
  
  # * Interpolation Step * - with conditional trailing NA handling
  mutate(
    ketones_post_int = list({
      row_vals <- c_across(all_of(ketone_apr24_cols_post))
      
      # If everything is NA across the entire row, just return them
      if (all(is.na(row_vals))) {
        row_vals
      } else {
        # Different interpolation based on date condition
        if (date_in_range) {
          # For dates between May 4 and May 7 (exclusive): interpolate both leading and trailing NAs
          as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 2)))
        } else {
          # For all other dates: only interpolate leading NAs, not trailing
          as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 1)))
        }
      }
    })
  ) %>%
  ungroup() %>%
  
  # Unnest the interpolated values into new columns
  unnest_wider(ketones_post_int, names_sep = '_') %>%
  
  # Remove the temporary flag
  select(-date_in_range)

# Checking the effects of posttest ketones imputation and filtering based on missingness
participants_out_post_ketones <- setdiff(kd_apr_24_CLEAN$participant_id, ketones_apr24_imputed_post$participant_id)
participants_out_post_ketones

# - Implementing imputed ketone values -------------------------------------

glimpse(ketones_apr24_imputed_pre)
glimpse(ketones_apr24_imputed_post)

# Tibble with the final (imputed) ketones
kd_apr_24_FINAL <- kd_apr_24_CLEAN %>%
  # Remove old ketone columns
  select(-matches('^ketones_\\d+_(pre|post)$')) %>%
  # Add imputed pretest ketones
  right_join(ketones_apr24_imputed_pre %>% select(participant_id, ketones_pre,
                                                  matches('^ketones_pre_int_\\d+$')),
             by = 'participant_id') %>%
  # Add imputed posttest ketones; imputed dataset excluded participants 
  # based on ketone missingness --> right join
  right_join(ketones_apr24_imputed_post %>% select(participant_id, matches('^ketones_post_int_\\d+$')), 
             by = 'participant_id') %>%
  # Drop unnecessary columns (original height and weight)
  select(-c(Q73.1, Q74.1))

# PREPROCESSING KD SEPTEMBER 2024 ----------------------------
# - Pretest cleaning ---------------------------------

kd_sep_24_pre_CLEAN <- kd_sep_24_pre %>%

  # Remove unnecessary rows (contain metadata)
  slice(-c(1, 2)) %>%
  
  # Convert key columns
  mutate(
    weight = as.numeric(Q73.1), # Ensure weight is numeric
    height = as.numeric(Q74.1) / 100, # Convert height from cm to meters
    
    # Correct weight values - weight in KG
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ),
    
    # Clean up sex column
    sex = case_when(
      Q75 == 'Vrouw (XX)' ~ 'female',
      Q75 == '' ~ NA_character_,  # Explicitly setting NA for missing values
      .default = NA_character_
    )
  ) %>%
  
  # Remove unnecessary columns
  select(-Q72, -Q73, -Q74, -RecordedDate, -DistributionChannel, -UserLanguage, 
         -Duration..in.seconds., -Status, -Progress, -ResponseId, -Q75) %>%
  
  # Rename other necessary columns
  rename(
    age = Q72.1, measured_ketones = Q75.1, ketones_pre = Q76_1_1,
    current_diet = Q77, current_diet_text = Q77_6_TEXT, 
    start_date = StartDate, end_date = EndDate, 
    finished = Finished
  ) %>%
  
  # Add computed columns
  mutate(
    BMI = weight / (height^2),
    
    # Convert 'measured_ketones' from text to logical
    measured_ketones = case_when(
      measured_ketones == 'Ik heb mijn ketowaarden gemeten / ik ga dit nu doen' ~ TRUE,
      measured_ketones == 'Ik heb mijn ketowaarden niet gemeten' ~ FALSE,
      .default = NA  # Ensure it remains logical
    ),
    
    # Make age numeric
    age = as.numeric(age),
    
    # Convert 'Finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  
  # Relocation
  relocate(participant_id, sex, weight, height, BMI, end_date) %>%
  # Final adjustments
  rename(pretest_date = end_date,
         pretest_time = end_time) %>%
  mutate(pretest_date = as.Date(pretest_date)) %>%  # Convert to short date format
  select(-start_date)

# Investigating duplicates #
# -------------------------------- #

duplicates_sep24_pre <- kd_sep_24_pre_CLEAN %>%
  filter(finished, measured_ketones) %>% # Select only cases where column == TRUE
  filter(BMI >= 18.5 & BMI < 30,
         sex == 'female',
         age >= 35 & age <= 45
  ) %>%
  # Get only duplicates again after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(pretest_date, pretest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# ------------------------------- #

kd_sep_24_pre_FRESH <- kd_sep_24_pre_CLEAN %>% 
  setdiff(duplicates_sep24_pre)  # Exclude duplicates from pretest

# - Posttest cleaning ----------------------------------

kd_sep_24_post_CLEAN <- kd_sep_24_post %>%
  
  # Remove unnecessary columns
  select(-Status, -Progress, -Duration..in.seconds., -RecordedDate, 
         -ResponseId, -DistributionChannel, -UserLanguage, -Q72, -Q73, -Q74) %>%
  relocate(participant_id) %>%
  rename_with(~ str_replace(.x, '^Q76_(\\d+)_1$', 'ketones_\\1'), starts_with('Q76_')) %>%  # Dynamic renaming
  rename(
    weight = Q73.1, height = Q74.1, measured_ketones = Q75, 
    start_date = StartDate, end_date = EndDate, finished = Finished
  ) %>%
  # Remove unnecessary rows
  slice(-c(1, 2)) %>%
  mutate(
    weight = as.numeric(weight),  # Ensure weight is numeric
    height = as.numeric(height) / 100,  # Convert height from cm to meters
    
    # Correct weight values
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ), 
    
    # Adding BMI
    BMI = weight / (height^2),
    
    # Convert 'measured_ketones' from text to logical
    measured_ketones = case_when(
      measured_ketones == 'Ik heb mijn ketowaarden gemeten' ~ TRUE,
      measured_ketones == 'Ik heb mijn ketowaarden niet gemeten' ~ FALSE,
      .default = NA  # Ensure it remains logical
    ),
    
    # Convert 'Finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  # Final adjustments
  rename(posttest_date = end_date,
         posttest_time = end_time) %>%
  mutate(posttest_date = as.Date(posttest_date)) %>%
  select(-start_date)

# Investigating duplicates #
# --------------------------------- #

duplicates_sep24_post <- kd_sep_24_post_CLEAN %>%
  filter(finished, measured_ketones) %>% # Select only cases where column == TRUE
  # Get only duplicates after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(posttest_date, posttest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# -------------------------------- #

kd_sep_24_post_FRESH <- kd_sep_24_post_CLEAN %>% 
  setdiff(duplicates_sep24_post)  # Exclude duplicates from posttest

# - Combine datasets ------------------------------

# COMBINED DATA FRAME: Merging data sets
kd_sep_24_combined <- kd_sep_24_pre_FRESH %>%
  right_join(kd_sep_24_post_FRESH, 
  by = 'participant_id', 
  suffix = c('_pre', '_post')  # Automatically renames overlapping columns
) %>%
  rename_with(
    ~ sub('^ketones_([1-9]|1[0-6])$', 'ketones_\\1_post', .),
    matches('^ketones_([1-9]|1[0-6])$')
  ) %>%
  mutate(across(matches('^ketones_\\d+_post$'),
                ~ {
                  ketone <- parse_number(str_replace_all(.x, ',', '.'))
                  round(
                    ifelse(ketone >= 10,
                           ketone / 10^(nchar(trunc(ketone)) - 1),
                           ketone),
                    1
                  )
                }),
         # Convert pretest ketones to numeric and round
         ketones_pre = parse_number(str_replace_all(ketones_pre, ',', '.')),
         ketones_pre = round(ifelse(ketones_pre >= 10, 
                                    ketones_pre / (10^(nchar(as.character(round(ketones_pre)))) - 1),
                                    ketones_pre), 1)
         ) %>%
  mutate(current_diet = str_replace_all(current_diet, ' ', '_')) %>%
  dplyr::select(-c(Q73.1, Q74.1)) %>%
  mutate(age = as.numeric(age), group = 'KD', cohort = 'sep_24')
glimpse(kd_sep_24_combined)

# TIBBLE DESCRITION:
# - weight ~ in kg
# - height ~ in metres
# - ketones ~ in mmol/L

# FILTERING KD SEPTEMBER 2024 & IMPUTATION ---------------------------

# FILTERED DATASET (basic exclusions made) #
# ============================================== #

kd_sep_24_CLEAN <- kd_sep_24_combined %>%
  filter(finished_pre, finished_post,
         measured_ketones_pre, measured_ketones_post) %>% # Select only cases where column == TRUE
  filter(BMI_pre >= 18.5 & BMI_pre < 30,
         sex == 'female',
         age >= 35 & age <= 45,
  ) %>%
  # Exclude Paleo and Low carb diets
  filter(!current_diet %in% c('Paleo_dieet', 'Low_carb_dieet')) %>%
  arrange(participant_id)
glimpse(kd_sep_24_CLEAN)

# - Ketone missingness checks -------------------------------

# Check ketone missingness at Pretest
missing_pre_ketones_sep24 <- kd_sep_24_CLEAN %>%
  filter(is.na(as.numeric(ketones_pre))) %>%
  summarise(missing_count = n())
missing_pre_ketones_sep24

# Check missingness during intervention
missing_post_ketones_sep24 <- kd_sep_24_CLEAN %>%
  dplyr::select(participant_id, matches('^ketones_\\d+_post$')) %>%
  pivot_longer(cols = matches('^ketones_\\d+_post$'), 
               names_to = 'ketone_measurement', 
               values_to = 'ketone_level') %>%
  filter(is.na(ketone_level)) %>%
  group_by(participant_id) %>%
  summarise(missing_count = n(), .groups = 'drop') %>%
  arrange(desc(missing_count))
missing_post_ketones_sep24

# - Imputing ketones POST ----------------------------------

# POST ~ intervention ketones
# _int ~ interpolated values

# Filter:
# - Last measurement (16) can be missing but is not interpolated; exception to no interpolation being 
# posttest_date > 2024-10-11
# - No 2 adjacent measurements can be missing
# - Max. 3 / 8 measurements can be missing (from the last 2 weeks; the first 2 not important)
# - Max. 2 / 7 measurements can be missing if measurement 16 is missing (allowed by default)

# Imputation:
# - Linear interpolation in rows over max. 1 cell
# - Trailing NAs not interpolated unless posttest_date > 2024-10-11
# - Leading NAs interpolated

# Identify ALL posttest columns and sort them numerically
ketone_sep24_cols_post <- grep('^ketones_\\d+_post$', names(kd_sep_24_CLEAN), value = TRUE)
ketone_sep24_cols_post <- ketone_sep24_cols_post[order(as.numeric(gsub('.*_(\\d+)_post', '\\1', ketone_sep24_cols_post)))]

# Identify ONLY columns 9..16 for filtering and sort them numerically
filter_cols_9_16 <- grep('^ketones_(9|1[0-6])_post$', names(kd_sep_24_CLEAN), value = TRUE)
filter_cols_9_16 <- filter_cols_9_16[order(as.numeric(gsub('.*_(\\d+)_post', '\\1', filter_cols_9_16)))]

ketones_sep24_imputed_post <- kd_sep_24_CLEAN %>%

  # Ensure all posttest columns are numeric
  mutate(across(all_of(ketone_sep24_cols_post), as.numeric)) %>%
  
  # Process each row individually
  rowwise() %>%
  
  # * Filter Step * - only for columns 9..16
  filter({
    row_vals_9_16 <- c_across(all_of(filter_cols_9_16))
    
    # Helper function: Are any two consecutive NA?
    is_consecutive_na <- function(x) {
      for (i in seq_along(x)[-1]) {
        if (is.na(x[i]) && is.na(x[i - 1])) {
          return(TRUE)
        }
      }
      return(FALSE)
    }
    
    # Count how many NAs
    n_missing <- sum(is.na(row_vals_9_16))
    
    # Check if measurement 16 is missing (position 8 in the vector)
    vals <- as.vector(row_vals_9_16)
    measurement_16_missing <- length(vals) >= 8 && is.na(vals[8])
    
    # Keep row only if:
    # - No consecutive NA in 9..16
    # - At most 3 total missing in 9..16
    # - At most 2 total missing in 9 .. 15 if measurement 16 is missing
    
    !is_consecutive_na(row_vals_9_16) && 
      (n_missing <= ifelse(measurement_16_missing, 2, 3))
  }) %>%
  
  # * Interpolation Step * - with conditional trailing NA handling
  mutate(
    ketones_post_int = list({
      row_vals <- c_across(all_of(ketone_sep24_cols_post))
      
      # If everything is NA across the entire row, just return them
      if (all(is.na(row_vals))) {
        row_vals
      } else if (posttest_date > as.Date('2024-10-11')) {
        # Interpolate leading and also trailing NAs
        as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 2)))
      }
        else {
        # Interpolate leading NAs, not trailing
        as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 1)))
      }
    })
  ) %>%
  
  ungroup() %>%
  # Unnest the interpolated values into new columns
  unnest_wider(ketones_post_int, names_sep = '_')

# - Implementing imputed ketone values -------------------------------------

glimpse(ketones_sep24_imputed_post)

# Tibble with the final (imputed) ketones
kd_sep_24_FINAL <- kd_sep_24_CLEAN %>%
  # Remove old ketone columns
  select(-matches('^ketones_\\d+_post$')) %>%
  # Add imputed posttest ketones; imputed dataset excluded participants 
  # based on ketone missingness --> right join
  right_join(ketones_sep24_imputed_post %>% select(participant_id, matches('^ketones_post_int_\\d+$')), 
             by = 'participant_id')
glimpse(kd_sep_24_FINAL)

# PREPROCESSING KD JANUARY 2025 -------------------------------
# - Pretest cleaning ---------------------------------

kd_jan_25_pre_CLEAN <- kd_jan_25_pre %>%
  
  # Remove unnecessary rows (contain metadata)
  slice(-c(1, 2)) %>%
  
  # Convert key columns
  mutate(
    weight = as.numeric(Q73.1),  # Ensure weight is numeric
    height = as.numeric(Q74.1) / 100,  # Convert height from cm to meters
    
    # Correct weight values - weight in KG
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ),
    
    # Clean up sex column
    sex = case_when(
      Q75 == 'Vrouw (XX)' ~ 'female',
      Q75 == '' ~ NA_character_,  # Explicitly setting NA for missing values
      .default = NA_character_
    )
  ) %>%
  
  # Remove unnecessary columns
  select(-Q72, -Q73, -Q74, -RecordedDate, -DistributionChannel, -UserLanguage, 
         -Duration..in.seconds., -Status, -Progress, -ResponseId, -Q75) %>%
  
  # Rename other necessary columns
  rename(
    age = Q72.1, measured_ketones = Q75.1, ketones_pre = Q76_1_1,
    current_diet = Q77, current_diet_text = Q77_6_TEXT, 
    start_date = StartDate, end_date = EndDate, 
    finished = Finished
  ) %>%
  
  # Add computed columns
  mutate(
    BMI = weight / (height^2),
    
    # Convert 'measured_ketones' from text to logical
    measured_ketones = case_when(
      measured_ketones == 'Ik heb mijn ketowaarden gemeten / ik ga dit nu doen' ~ TRUE,
      measured_ketones == 'Ik heb mijn ketowaarden niet gemeten' ~ FALSE,
      .default = NA  # Ensure it remains logical
    ),
    
    # Make age numeric
    age = as.numeric(age),
    
    # Convert 'Finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  
  # Relocation
  relocate(participant_id, sex, weight, height, BMI, end_date) %>%
  # Final adjustments
  rename(pretest_date = end_date,
         pretest_time = end_time) %>%
  mutate(pretest_date = as.Date(pretest_date)) %>%  # Convert to short date format
  select(-start_date)

# Investigating participant churn #
# ----------------------------------- #

kd_jan_25_pre_CLEAN_view <- kd_jan_25_pre_CLEAN %>%
  filter(finished, measured_ketones) %>%
  filter(BMI >= 18.5 & BMI < 30) %>%
  filter(sex == 'female') %>%
  filter(age >= 35 & age <= 45) %>%
  filter(!participant_id == '') %>%
  filter(!current_diet %in% c('Paleo dieet', 'Low carb dieet'))

drop_out_jan_25 <- setdiff(kd_jan_25_pre_CLEAN_view$participant_id, 
                                 kd_jan_25_post$participant_id)
length(drop_out_jan_25)

# NOTE: 
# - 5 participants met the criteria but two didn't finish the posttest --> 3 eligible participants
# at posttest

# Investigating duplicates #
# -------------------------------- #

duplicates_jan25_pre <- kd_jan_25_pre_CLEAN %>%
  filter(finished, measured_ketones) %>% # Select only cases where column == TRUE
  filter(BMI >= 18.5 & BMI < 30,
         sex == 'female',
         age >= 35 & age <= 45
  ) %>%
  # Get only duplicates again after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(pretest_date, pretest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# ------------------------------- #

kd_jan_25_pre_FRESH <- kd_jan_25_pre_CLEAN %>% 
  setdiff(duplicates_jan25_pre)  # Exclude duplicates from pretest

# - Posttest cleaning ----------------------------------

kd_jan_25_post_CLEAN <- kd_jan_25_post %>%
  # Remove unnecessary columns
  select(-Status, -Progress, -Duration..in.seconds., -RecordedDate, 
         -ResponseId, -DistributionChannel, -UserLanguage, -Q72, -Q73, -Q74) %>%
  relocate(participant_id) %>%
  select(-c(starts_with('Q77_'))) %>% # Duplicated ketone values
  rename_with(~ str_replace(.x, '^Q76_(\\d+)_1$', 'ketones_\\1'), starts_with('Q76_')) %>%  # Dynamic renaming
  rename(
    weight = Q73.1, height = Q74.1, measured_ketones = Q75, 
    start_date = StartDate, end_date = EndDate, finished = Finished
  ) %>%
  # Remove unnecessary rows
  slice(-c(1, 2)) %>%
  mutate(
    weight = as.numeric(weight),  # Ensure weight is numeric
    height = as.numeric(height) / 100,  # Convert height from cm to meters
    
    # Correct weight values
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ), 
    
    # Adding BMI
    BMI = weight / (height^2),
    
    # Convert 'measured_ketones' from text to logical
    measured_ketones = case_when(
      measured_ketones == 'Ik heb mijn ketowaarden gemeten' ~ TRUE,
      measured_ketones == 'Ik heb mijn ketowaarden niet gemeten' ~ FALSE,
      .default = NA  # Ensure it remains logical
    ),
    # Convert 'Finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
  ) %>%
  # Final adjustments
  rename(posttest_date = end_date,
         posttest_time = end_time) %>%
  mutate(posttest_date = as.Date(posttest_date)) %>%
  select(-start_date)

# Investigating duplicates #
# --------------------------------- #

duplicates_jan25_post <- kd_jan_25_post_CLEAN %>%
  filter(finished, measured_ketones) %>% # Select only cases where column == TRUE
  # Get only duplicates after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(posttest_date, posttest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# -------------------------------- #

kd_jan_25_post_FRESH <- kd_jan_25_post_CLEAN %>% 
  setdiff(duplicates_jan25_post)  # Exclude duplicates from posttest

# - Combine datasets ------------------------------

# COMBINED DATA FRAME: Merging data sets
kd_jan_25_combined <- kd_jan_25_pre_FRESH %>%
  right_join(kd_jan_25_post_FRESH, 
            by = 'participant_id', 
            suffix = c('_pre', '_post')  # Automatically renames overlapping columns
  ) %>%
  rename_with(
    ~ sub('^ketones_([1-9]|1[0-6])$', 'ketones_\\1_post', .),
    matches('^ketones_([1-9]|1[0-6])$')
  ) %>%
  mutate(across(matches('^ketones_\\d+_post$'),
                ~ {
                  ketone <- parse_number(str_replace_all(.x, ',', '.'))
                  round(
                    ifelse(ketone >= 10,
                           ketone / 10^(nchar(trunc(ketone)) - 1),
                           ketone),
                    1
                  )
                }),
         # Convert pretest ketones to numeric and round
         ketones_pre = parse_number(str_replace_all(ketones_pre, ',', '.')),
         ketones_pre = round(ifelse(ketones_pre >= 10, 
                                    ketones_pre / (10^(nchar(as.character(round(ketones_pre)))) - 1),
                                    ketones_pre), 1)
         ) %>%
  mutate(current_diet = str_replace_all(current_diet, ' ', '_')) %>%
  mutate(age = as.numeric(age),
           group = 'KD',
           cohort = 'jan_25') %>%
  dplyr::select(-c(Q73.1, Q74.1))
glimpse(kd_jan_25_combined)

# TIBBLE DESCRITION:
# - weight ~ in kg
# - height ~ in metres
# - ketones ~ in mmol/L

# FILTERING KD JANUARY 2025 & IMPUTATION ---------------------------

# FILTERED DATASET (basic exclusions made) #
# ============================================== #

kd_jan_25_CLEAN <- kd_jan_25_combined %>%
  filter(finished_pre, finished_post,
         measured_ketones_pre, measured_ketones_post) %>% # Select only cases where column == TRUE
  filter(BMI_pre >= 18.5 & BMI_pre < 30,
         sex == 'female',
         age >= 35 & age <= 45,
  ) %>%
  # Exclude Paleo and Low carb diets
  filter(!current_diet %in% c('Paleo_dieet', 'Low_carb_dieet')) %>%
  arrange(participant_id)
glimpse(kd_jan_25_CLEAN)

# - Ketone missingness checks -------------------------------

# Check ketone missingness at Pretest
missing_pre_ketones_jan25 <- kd_jan_25_CLEAN %>%
  filter(is.na(as.numeric(ketones_pre))) %>%
  summarise(missing_count = n())
missing_pre_ketones_jan25

# Check missingness during intervention
missing_post_ketones_jan25 <- kd_jan_25_CLEAN %>%
  dplyr::select(participant_id, matches('^ketones_\\d+_post$')) %>%
  pivot_longer(cols = matches('^ketones_\\d+_post$'), 
               names_to = 'ketone_measurement', 
               values_to = 'ketone_level') %>%
  filter(is.na(as.numeric(ketone_level))) %>%
  group_by(participant_id) %>%
  summarise(missing_count = n(), .groups = 'drop') %>%
  arrange(desc(missing_count))
missing_post_ketones_jan25

# - Imputing ketones POST ----------------------------------

# POST ~ intervention ketones

# Filter:
# - Last measurement (16) can be missing but is not interpolated; exception to no interpolation being 
# posttest_date > 2025-02-14
# - No 2 adjacent measurements can be missing
# - Max. 3 / 8 measurements can be missing (from the last 2 weeks; the first 2 not important)
# - Max. 2 / 7 measurements can be missing if measurement 16 is missing (allowed by default)

# Imputation:
# - Linear interpolation in rows over max. 1 cell
# - Trailing NAs not interpolated unless posttest_date > 2025-02-14
# - Leading NAs interpolated

# Identify ALL posttest columns and sort them numerically
ketone_jan25_cols_post <- grep('^ketones_\\d+_post$', names(kd_jan_25_CLEAN), value = TRUE)
ketone_jan25_cols_post <- ketone_jan25_cols_post[order(as.numeric(gsub('.*_(\\d+)_post', '\\1', ketone_jan25_cols_post)))]

# Identify ONLY columns 9..16 for filtering and sort them numerically
filter_cols_9_16 <- grep('^ketones_(9|1[0-6])_post$', names(kd_jan_25_CLEAN), value = TRUE)
filter_cols_9_16 <- filter_cols_9_16[order(as.numeric(gsub('.*_(\\d+)_post', '\\1', filter_cols_9_16)))]

ketones_jan25_imputed_post <- kd_jan_25_CLEAN %>%
  
  # Ensure all posttest columns are numeric
  mutate(across(all_of(ketone_jan25_cols_post), as.numeric)) %>%
  
  # Process each row individually
  rowwise() %>%
  
  # * Filter Step * - only for columns 9..16
  filter({
    row_vals_9_16 <- c_across(all_of(filter_cols_9_16))
    
    # Helper function: Are any two consecutive NA?
    is_consecutive_na <- function(x) {
      for (i in seq_along(x)[-1]) {
        if (is.na(x[i]) && is.na(x[i - 1])) {
          return(TRUE)
        }
      }
      return(FALSE)
    }
    
    # Count how many NAs
    n_missing <- sum(is.na(row_vals_9_16))
    
    # Check if measurement 16 is missing (position 8 in the vector)
    vals <- as.vector(row_vals_9_16)
    measurement_16_missing <- length(vals) >= 8 && is.na(vals[8])
    
    # Keep row only if:
    # - No consecutive NA in 9..16
    # - At most 3 total missing in 9..16
    # - At most 2 total missing in 9 .. 15 if measurement 16 is missing
    
    !is_consecutive_na(row_vals_9_16) && 
      (n_missing <= ifelse(measurement_16_missing, 2, 3))
  }) %>%
  
  # * Interpolation Step * - with conditional trailing NA handling
  mutate(
    ketones_post_int = list({
      row_vals <- c_across(all_of(ketone_jan25_cols_post))
      
      # If everything is NA across the entire row, just return them
      if (all(is.na(row_vals))) {
        row_vals
      } else if (posttest_date > as.Date('2025-02-14')) {
        # Interpolate leading and also trailing NAs
        as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 2)))
      }
      else {
        # Interpolate leading NAs, not trailing
        as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 1)))
      }
    })
  ) %>%
  
  ungroup() %>%
  # Unnest the interpolated values into new columns
  unnest_wider(ketones_post_int, names_sep = '_')
  
# - Implementing imputed ketone values -------------------------------------

glimpse(ketones_jan25_imputed_post)

# Tibble with the final (imputed) ketones
kd_jan_25_FINAL <- kd_jan_25_CLEAN %>%
  # Remove old ketone columns
  select(-matches('^ketones_\\d+_post$')) %>%
  # Add imputed posttest ketones; imputed dataset excluded participants 
  # based on ketone missingness --> right join
  right_join(ketones_jan25_imputed_post %>% select(participant_id, matches('^ketones_post_int_\\d+$')), 
             by = 'participant_id')
glimpse(kd_jan_25_FINAL)

# PREPROCESSING CD 2025 -------------------------------
# - Pretest cleaning ---------------------------------

cd_25_pre_CLEAN <- cd_25_pre %>%

  # Remove unnecessary rows (contain metadata)
  slice(-c(1, 2)) %>%
  
  # Convert key columns
  mutate(
    weight = as.numeric(Q73.1),  # Ensure weight is numeric
    height = as.numeric(Q74.1) / 100,  # Convert height from cm to meters
    
    # Correct weight values - weight in KG
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ),
    
    # Clean up sex column
    sex = case_when(
      Q75 == 'Vrouw (XX)' ~ 'female',
      Q75 == '' ~ NA_character_,  # Explicitly setting NA for missing values
      .default = NA_character_ 
    )
  ) %>%
  
  # Remove unnecessary columns
  select(-Q72, -Q73, -Q74, -RecordedDate, -DistributionChannel, -UserLanguage, 
         -Duration..in.seconds., -Status, -Progress, -ResponseId, -Q75) %>%
  
  # Rename other necessary columns
  rename(
    age = Q72.1, current_diet = Q77, current_diet_text = Q77_4_TEXT, 
    start_date = StartDate, end_date = EndDate, 
    finished = Finished
  ) %>%
  
  # Add computed columns
  mutate(
    BMI = weight / (height^2),
    
    # Convert 'Finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Make age numeric
    age = as.numeric(age),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  
  # Relocation
  relocate(participant_id, sex, weight, height, BMI, end_date) %>%
  # Final adjustments
  rename(pretest_date = end_date,
         pretest_time = end_time) %>%
  mutate(pretest_date = as.Date(pretest_date), # Convert to short date format
         ketones_pre = NA_real_) %>%  # Add missing values for consistency
  select(-start_date) %>%
  # Get rid of empty participant IDs
  mutate(participant_id = na_if(str_remove_all(participant_id, '\\s+'), '')) %>%
  filter(!is.na(participant_id)) # Remove empty participant IDs
glimpse(cd_25_pre_CLEAN)

# Investigating duplicates #
# -------------------------------- #

duplicates_cd25_pre <- cd_25_pre_CLEAN %>%
  filter(finished) %>% # Select only cases where column == TRUE
  filter(BMI >= 18.5 & BMI < 30,
         sex == 'female',
         age >= 35 & age <= 45
  ) %>%
  # Get only duplicates again after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(pretest_date, pretest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# ------------------------------- #

cd_25_pre_FRESH <- cd_25_pre_CLEAN %>% 
  setdiff(duplicates_cd25_pre)  # Exclude duplicates from pretest

# - Posttest cleaning ----------------------------------

cd_25_post_CLEAN <- cd_25_post %>%
  
  # Remove unnecessary columns
  select(-Status, -Progress, -Duration..in.seconds., -RecordedDate, 
         -ResponseId, -DistributionChannel, -UserLanguage, -Q72, -Q73, -Q74) %>%
  relocate(participant_id) %>%
  
  # Dynamic renaming
  rename_with(~ str_replace(.x, '^Q76_(\\d+)_1$', 'ketones_\\1'), starts_with('Q76_')) %>%  
  rename(
    weight = Q73.1, height = Q74.1, start_date = StartDate, end_date = EndDate, 
    finished = Finished
  ) %>%
  
  # Remove unnecessary rows
  slice(-c(1, 2)) %>%
  
  # Convert key columns
  mutate(
    weight = as.numeric(weight),  # Ensure weight is numeric
    height = as.numeric(height) / 100,  # Convert height from cm to meters
    
    # Correct weight values
    weight = case_when(
      !is.na(weight) & abs(weight) > 200 & nchar(as.character(round(weight))) > 2 ~ 
        weight / (10^(nchar(as.character(round(weight))) - 2)),  
      .default = weight  # Keep the original weight if conditions aren't met
    ), 
    
    # Adding BMI
    BMI = weight / (height^2),
    
    # Convert 'Finished' column from text to logical
    finished = case_when(
      finished == 'True' ~ TRUE,
      finished == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  
  # Additional adjustments
  rename(posttest_date = end_date,
         posttest_time = end_time) %>%
  mutate(posttest_date = as.Date(posttest_date)) %>%
  select(-start_date) %>%
  # Get rid of empty participant IDs
  mutate(participant_id = na_if(str_remove_all(participant_id, '\\s+'), '')) %>%
  filter(!is.na(participant_id)) %>% # Remove empty participant IDs
  rename(foodlog_starches = Q80_1.1,
         foodlog_bread = Q80_8,
         foodlog_sweetpotatoe = Q80_13,
         foodlog_dairy = Q80_2,
         foodlog_fats = Q80_3,
         foodlog_sweeteners = Q80_4,
         foodlog_legumes = Q80_6,
         foodlog_meat = Q80_7,
         foodlog_berries = Q80_10,
         foodlog_fruit = Q80_9,
         foodlog_greens = Q80_11,
         foodlog_vegetables = Q80_12,
         foodlog_complex_starches = Q80_14
         ) %>%
  # Substitute any NAs with 0
  mutate(across(starts_with('foodlog_'),
                ~ coalesce(parse_number(.x), 0)))

# Investigating duplicates #
# --------------------------------- #

duplicates_cd25_post <- cd_25_post_CLEAN %>%
  filter(finished) %>% # Select only cases where column == TRUE
  # Get only duplicates after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(posttest_date, posttest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# -------------------------------- #

cd_25_post_FRESH <- cd_25_post_CLEAN %>% 
  setdiff(duplicates_cd25_post)  # Exclude duplicates from posttest

# - Posttest ketone values cleaning --------------------------

# This is an extra step for CD 2025 because ketones were recorded 
# separately from the usual questions in this cohort

cd_25_post_ketones_CLEAN <- cd_25_post_ketones %>%
  
  # Remove unnecessary columns
  select(-Status, -Progress, -Duration..in.seconds., -RecordedDate, 
         -ResponseId, -DistributionChannel, -UserLanguage) %>%
  
  # Relocate participant_id
  relocate(participant_id) %>%
  
  # Dynamic renaming
  rename_with(~ str_replace(.x, '^Ketone.values_(\\d+)_1$', 'ketones_\\1'), 
              starts_with('Ketone.values_')) %>%  
  rename(
    measured_ketones_post = Q2, start_date = StartDate, end_date = EndDate, 
    finished_ketones = Finished
  ) %>%
  
  # Remove unnecessary rows
  slice(-c(1, 2)) %>%
  
  # Convert key columns
  mutate(
    
    measured_ketones_post = case_when(
      measured_ketones_post == 'Ik heb mijn ketowaarden gemeten' ~ TRUE,
      measured_ketones_post == 'Ik heb mijn ketowaarden niet gemeten' ~ FALSE,
      .default = NA  # Ensure it remains logical
    ),
    
    # Convert 'finished' column from text to logical
    finished_ketones = case_when(
      finished_ketones == 'True' ~ TRUE,
      finished_ketones == 'False' ~ FALSE,
      .default = NA
    ),
    
    # Create a time column
    end_time = hms::as_hms(ymd_hms(end_date))
    
  ) %>%
  
  # Final adjustments
  rename(posttest_date = end_date,
         posttest_time = end_time) %>%
  mutate(posttest_date = as.Date(posttest_date)) %>%
  select(-start_date) %>%
  # Get rid of empty participant IDs
  mutate(participant_id = na_if(str_remove_all(participant_id, '\\s+'), '')) %>%
  filter(!is.na(participant_id)) # Remove empty participant IDs
glimpse(cd_25_post_ketones_CLEAN)

# Investigating duplicates #
# --------------------------------- #

duplicates_cd25_post_ketones <- cd_25_post_ketones_CLEAN %>%
  filter(finished_ketones, measured_ketones_post) %>% # Select only cases where column == TRUE
  # Get only duplicates after the filtering above
  filter(duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)) %>%
  group_by(participant_id) %>%
  arrange(posttest_date, posttest_time, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Remove duplicates #
# -------------------------------- #

cd_25_post_ketones_FRESH <- cd_25_post_ketones_CLEAN %>% 
  setdiff(duplicates_cd25_post_ketones) %>% # Exclude duplicates from posttest
  # Filter and drop
  filter(finished_ketones, measured_ketones_post) %>% # Select only cases where column == TRUE
  select(-c(posttest_time, posttest_date)) # Drop time and date columns

# - Combine datasets ------------------------------

# COMBINED DATA FRAME: Merging data sets
cd_25_combined <- cd_25_post_FRESH %>%
  # Add ketone values to the Posttest
  full_join(cd_25_post_ketones_FRESH, by = 'participant_id') %>%
  # Add Pretest data
  left_join(cd_25_pre_FRESH, 
             by = 'participant_id', 
             suffix = c('_post', '_pre')  # Automatically renames overlapping columns
  ) %>%
  rename_with(
    ~ sub('^ketones_([1-9]|1[0-5])$', 'ketones_\\1_post', .),
    matches('^ketones_([1-9]|1[0-5])$')
  ) %>%
  # Turn decimal commas into dots
  mutate(across(
    matches('^ketones_\\d+_post$'),
    ~ {
      ketone <- parse_number(str_replace_all(.x, ',', '.'))
      round(
        ifelse(ketone >= 10,
               ketone / 10^(nchar(trunc(ketone)) - 1),
               ketone),
        1
      )
    }
  )) %>% 
  mutate(current_diet = str_replace_all(current_diet, ' ', '_')) %>%
  mutate(age = as.numeric(age),
         # Add metadata
         group = 'CD',
         cohort = 'cd_25',
         ) %>%
  dplyr::select(-c(Q73.1, Q74.1))
glimpse(cd_25_combined)

# TIBBLE DESCRITION:
# - weight ~ in kg
# - height ~ in metres
# - ketones ~ in mmol/L

# FILTERING CD 2025 & IMPUTATION ---------------------------

# FILTERED DATASET (basic exclusions made) #
# ============================================== #

cd_25_CLEAN <- cd_25_combined %>%
  filter(finished_pre, finished_post, finished_ketones, measured_ketones_post) %>%
  filter(BMI_pre >= 18.5 & BMI_pre < 30,
         sex == 'female',
         age >= 35 & age <= 45,
  ) %>%
  # Exclude Paleo and Low carb diets
  filter(!current_diet %in% c('Paleo_dieet', 'Low_carb_dieet')) %>%
  arrange(participant_id) %>%
  select(-c(finished_ketones)) %>%
  # Exclude based on food log
  group_by(participant_id) %>%
  mutate(normal_diet_score = sum(foodlog_starches, foodlog_bread, foodlog_sweetpotatoe, 
                foodlog_dairy, foodlog_sweeteners, 
                foodlog_legumes, foodlog_meat, foodlog_berries, 
                foodlog_fruit, foodlog_greens, foodlog_vegetables,
                foodlog_complex_starches),
         keto_score = sum(foodlog_fats, foodlog_dairy)) %>%
  ungroup() %>%
  filter(normal_diet_score > keto_score) %>%
  select(-c(normal_diet_score, keto_score), -c(starts_with('foodlog_')))
glimpse(cd_25_CLEAN)

# Check if anyone gets filtered out because of the diet filter
cd_25_combined %>%
  filter(finished_pre, finished_post, finished_ketones, measured_ketones_post) %>%
  filter(BMI_pre >= 18.5 & BMI_pre < 30,
         sex == 'female',
         age >= 35 & age <= 45,
  ) %>%
  # Exclude Paleo and Low carb diets
  filter(!current_diet %in% c('Paleo_dieet', 'Low_carb_dieet')) %>%
  arrange(participant_id) %>%
  select(-c(finished_ketones)) %>%
  # Exclude based on food log
  group_by(participant_id) %>%
  mutate(normal_diet_score = sum(foodlog_starches, foodlog_bread, foodlog_sweetpotatoe, 
                                 foodlog_dairy, foodlog_sweeteners, 
                                 foodlog_legumes, foodlog_meat, foodlog_berries, 
                                 foodlog_fruit, foodlog_greens, foodlog_vegetables,
                                 foodlog_complex_starches),
         keto_score = sum(foodlog_fats, foodlog_dairy)) %>%
  ungroup() %>%
  filter(!(normal_diet_score > keto_score)) %>%
  nrow()

# - Ketone missingness checks -------------------------------

# Check missingness during intervention
missing_post_ketones_cd25 <- cd_25_CLEAN %>%
  dplyr::select(participant_id, matches('^ketones_\\d+_post$')) %>%
  pivot_longer(cols = matches('^ketones_\\d+_post$'), 
               names_to = 'ketone_measurement',
               values_to = 'ketone_level') %>%
  filter(is.na(as.numeric(ketone_level))) %>%
  group_by(participant_id) %>%
  summarise(missing_count = n(), .groups = 'drop') %>%
  arrange(desc(missing_count))
missing_post_ketones_cd25

# - Imputing ketones POST ----------------------------------

# POST ~ intervention ketones

# Filter:
# - Last measurement (15) can be missing and but is not interpolated; exception to no interpolation being 
# posttest_date > 2025-04-02
# - No 2 adjacent measurements can be missing
# - Max. 6 / 15 measurements can be missing (from all 4 weeks)

# Imputation:
# - Linear interpolation across rows over max. 1 cell
# - Trailing NAs not interpolated unless posttest_date > 2025-04-02
# - Leading NAs interpolated

# Identify ALL posttest columns and sort them numerically
ketone_cd25_cols_post <- grep('^ketones_\\d+_post$', names(cd_25_CLEAN), value = TRUE)
ketone_cd25_cols_post <- ketone_cd25_cols_post[order(as.numeric(gsub('.*_(\\d+)_post', '\\1', ketone_cd25_cols_post)))]

ketones_cd25_imputed_post <- cd_25_CLEAN %>%
  
  # Ensure all posttest columns are numeric
  mutate(across(all_of(ketone_cd25_cols_post), as.numeric)) %>%
  
  # Process each row individually
  rowwise() %>%
  
  # * Filter Step * - all columns 1..15
  filter({
    row_vals <- c_across(all_of(ketone_cd25_cols_post))
    
    # Helper function: Are any two consecutive NA?
    is_consecutive_na <- function(x) {
      for (i in seq_along(x)[-1]) {
        if (is.na(x[i]) && is.na(x[i - 1])) {
          return(TRUE)
        }
      }
      return(FALSE)
    }
    
    # Count how many NAs
    n_missing <- sum(is.na(row_vals))
    
    # Keep row only if:
    # - No consecutive NA in the row
    # - At most 6 total missing in the row
    
    !is_consecutive_na(row_vals) && n_missing <= 6
  }) %>%
  
  # * Interpolation Step *
  mutate(
    ketones_post_int = list({
      row_vals <- c_across(all_of(ketone_cd25_cols_post))
      
      # If everything is NA across the entire row, just return them
      if (all(is.na(row_vals))) {
        row_vals
      } else if (posttest_date > as.Date('2025-04-02')) {
        # Interpolate leading and also trailing NAs
        as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 2)))
      }
      else {
        # Interpolate leading NAs, not trailing
        as.numeric(na.approx(row_vals, maxgap = 1, na.rm = FALSE, rule = c(2, 1)))
      }
    })
  ) %>%
  
  ungroup() %>%
  # Unnest the interpolated values into new columns
  unnest_wider(ketones_post_int, names_sep = '_')

# - Implementing imputed ketone values -------------------------------------

glimpse(ketones_cd25_imputed_post)

# Tibble with the final (imputed) ketones
cd_25_FINAL <- cd_25_CLEAN %>%
  # Remove old ketone columns
  select(-matches('^ketones_\\d+_post$')) %>%
  # Add imputed posttest ketones; imputed dataset excluded participants 
  # based on ketone missingness --> right join
  right_join(ketones_cd25_imputed_post %>% select(participant_id, matches('^ketones_post_int_\\d+$')), 
             by = 'participant_id') %>%
  mutate(
    ketones_post_int_16 = NA_real_  # Placeholder for the 16th measurement
    )
glimpse(cd_25_FINAL)

# - Compare final ketones to raw ketones ------------------------------

# Just comparing the final ketone values with the raw ones through 'view' side by side

cd_25_ketones_final <- cd_25_FINAL %>%
  select(participant_id, starts_with('ketones_post_int_')) %>%
  rename_with(~ str_replace(.x, '^ketones_post_int_(\\d+)$', 'ketones_\\1'), 
              starts_with('ketones_post_int_')) %>%
  mutate(identifier = 'final') %>%
  relocate(participant_id, identifier)

cd_25_ketones_raw <- cd_25_post_ketones %>%
  select(participant_id, starts_with('Ketone.values_')) %>%
  filter(participant_id %in% cd_25_ketones_final$participant_id) %>%
  rename_with(~ str_replace(.x, '^Ketone.values_(\\d+)_1$', 'ketones_\\1'), 
              starts_with('Ketone.values_')) %>%
  mutate(identifier = 'raw') %>%
  relocate(participant_id, identifier)

# COMBINE ALL DATASETS -------------------------------

dataset_start <- bind_rows(
  kd_apr_24_FINAL,
  kd_sep_24_FINAL,
  kd_jan_25_FINAL,
  cd_25_FINAL
) %>%
  # Arrange by cohort
  arrange(cohort) %>%
  # Q95 and Q94_1 are additional questions specific to CD
  select(-c(pretest_time, posttest_time, Q95, Q94_1))
glimpse(dataset_start)

# - Filter based on diet text -----------------------------

# Get unique current diet text
unique(dataset_start$current_diet_text)

dataset_diettext <- dataset_start %>%
  filter(!(current_diet_text == 'glutenvrij, en lactose en koolhydraat arm' | current_diet_text == 'Orthomoleculair paleo' | current_diet_text == 'Oersterk'))

# - Filter based on ketones -------------------------

# Filter Pretest:
# - Neither group in ketosis (< 0.5 mmol/L)
# - April 2024 KD cohort: additionally no measurements in ketosis 
# during baseline 

# Filter Posttest:
# - KD max. 1 measurement out of ketosis in Weeks 3 and 4 
# - KD in ketosis at Posttest
# - CD max. 1 measurement in ketosis in Weeks 3 and 4 
# - CD out of ketosis at Posttest

dataset_ketout <- dataset_diettext %>%
  rowwise() %>% # Operate one row at a time
  # Get the last actual value, not NA
  mutate(last_ketones_post_int =
           last(na.omit(c_across(starts_with('ketones_post_int_'))))) %>%
  ungroup() %>%
  filter(
    
    # Pretest conditions #
    # ------------------------------ #
    ketones_pre < 0.5 | is.na(ketones_pre),
    
    # Pretest: April 2024 KD cohort specific condition
    !(cohort == 'apr_24' & if_any(starts_with('ketones_pre_int_'), ~ coalesce(.x >= 0.5, FALSE))),
    
    # Posttest conditions #
    # ------------------------------ #
    
    # Acute KD ketosis at Posttest
    !(group == 'KD' & last_ketones_post_int < 0.5),
    # No CD acute ketosis at Posttest
    !(group == 'CD' & ketones_post_int_15 >= 0.5),
    
    # General filters on intervention ketones
    (group == 'KD' & # KD cohort
    rowSums( # How many fail?
      across(
        starts_with('ketones_post_int_') & # All post-int cols …
          !num_range('ketones_post_int_', 1:8), # … except 1–8
        ~ .x < 0.5 # TRUE when < 0.5
      ),
      na.rm = TRUE # NA ≡ no failure
    ) <= 1) 
    
    | # OR
      
    (group == 'CD' & # CD cohort
       rowSums( # How many fail?
         across(
           starts_with('ketones_post_int_')& # All post-int cols …
             !num_range('ketones_post_int_', 1:7), # … except 1–7
           ~ .x >= 0.5 # TRUE when >= 0.5
         ),
         na.rm = TRUE # NA ≡ no failure
       ) <= 1)
    ) %>%
  select(-last_ketones_post_int) # Remove the temporary column
glimpse(dataset_ketout)

# Check number of participants in each group after filtering
count_after_ketone_exclusions <- dataset_ketout %>%
  count(group) %>%
  ungroup()
count_after_ketone_exclusions

# FILTERING ASSESSMENT ---------------------------------------
# - Assessing individual filter effects --------------------------

# Count before ketone exclusions
count_before_ketone_exclusions <- dataset_diettext %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  ungroup()
count_before_ketone_exclusions

# Pretest ketones
ketout_check <- dataset_diettext %>%
  filter(
    ketones_pre < 0.5 | is.na(ketones_pre)) %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  ungroup()
ketout_check

# Pretest: April 2024 KD cohort specific condition
apr_24_ketones_check <- dataset_diettext %>%
  filter(!(cohort == 'apr_24' & if_any(starts_with('ketones_pre_int_'),
                                       ~ coalesce(.x >= 0.5, FALSE)))) %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  ungroup()
apr_24_ketones_check

# Posttest KD
KD_post_check <- dataset_diettext %>%
  filter(group == 'KD' & # KD cohort
           rowSums( # How many fail?
             across(
               starts_with('ketones_post_int_') & # All post-int cols …
                 !num_range('ketones_post_int_', 1:8), # … except 1–8
               ~ .x < 0.5 # TRUE when < 0.5
             ),
             na.rm = TRUE # NA ≡ no failure
           ) <= 1) %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  ungroup()
KD_post_check

# Posttest CD (0 values in ketosis allowed)
CD_post_check <- dataset_diettext %>%
  filter(group == 'CD' & if_all(starts_with('ketones_post_int_'),
                                ~ coalesce(.x < 0.5, TRUE))) %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  ungroup()
CD_post_check

# Posttest CD 2 (1 value in ketosis allowed)
CD_post_check_2 <- dataset_diettext %>%
  filter(group == 'CD' & 
           rowSums( # How many fail?
             across(
               starts_with('ketones_post_int_'), # All post-int cols
               ~ .x >= 0.5 # TRUE when >= 0.5
             ),
             na.rm = TRUE # NA ≡ no failure
           ) <= 2) %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  ungroup()
CD_post_check_2

# Posttest CD 3 (count only Weeks 3 and 4)
CD_post_check_3 <- dataset_diettext %>%
  filter(group == 'CD' & rowSums( # How many fail?
    across(
      starts_with('ketones_post_int_') & # All post-int cols …
        !num_range('ketones_post_int_', 1:7), ~ .x >= 0.5 # TRUE when >= 0.5
    ),
    na.rm = TRUE # NA ≡ no failure
  ) <= 1) %>%
  filter(ketones_post_int_15 <= 0.5) %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  ungroup()
CD_post_check_3

# - Assessing exclusions cohort-wise --------------------------

# KD APRIL 2024 #
# ============================== #

# Pretest before filtering
count_kd_apr_24_pre <- kd_apr_24_pre_FRESH %>%
  filter(finished) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Posttest before filtering
count_kd_apr_24_post <- kd_apr_24_combined %>%
  filter(finished_pre, finished_post) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Count after filtering
count_kd_apr_24_final <- dataset_ketout %>%
  filter(cohort == 'apr_24') %>%
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

cat('KD April 24: No. of participants at Pretest =', count_kd_apr_24_pre, '\n')
cat('KD April 24: No. of participants at Posttest =', count_kd_apr_24_post, '\n')
cat('KD April 24: Final no. of participants =', count_kd_apr_24_final, '\n')

# KD SEPTEMBER 2024 #
# =============================== #

# Pretest before filtering
count_kd_sep_24_pre <- kd_sep_24_pre_FRESH %>%
  filter(finished) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Posttest before filtering
count_kd_sep_24_post <- kd_sep_24_combined %>%
  filter(finished_pre, finished_post) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Count after filtering
count_kd_sep_24_final <- dataset_ketout %>%
  filter(cohort == 'sep_24') %>%
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

cat('KD September 24: No. of participants at Pretest =', count_kd_sep_24_pre, '\n')
cat('KD September 24: No. of participants at Posttest =', count_kd_sep_24_post, '\n')
cat('KD September 24: Final no. of participants =', count_kd_sep_24_final, '\n')

# KD JANUARY 2025 #
# =============================== #

# Pretest before filtering
count_kd_jan_25_pre <- kd_jan_25_pre_FRESH %>%
  filter(finished) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Posttest before filtering
count_kd_jan_25_post <- kd_jan_25_combined %>%
  filter(finished_pre, finished_post) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Count after filtering
count_kd_jan_25_final <- dataset_ketout %>%
  filter(cohort == 'jan_25') %>%
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

cat('KD January 25: No. of participants at Pretest =', count_kd_jan_25_pre, '\n')
cat('KD January 25: No. of participants at Posttest =', count_kd_jan_25_post, '\n')
cat('KD January 25: Final no. of participants =', count_kd_jan_25_final, '\n')

# CD 2025 #
# =============================== #

# Pretest before filtering
count_cd_25_pre <- cd_25_pre_FRESH %>%
  filter(finished) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Posttest before filtering
count_cd_25_post <- cd_25_combined %>%
  filter(finished_pre, finished_post) %>%
  filter(!(participant_id == '' | participant_id == ' ')) %>% # Remove empty participant IDs
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

# Count after filtering
count_cd_25_final <- dataset_ketout %>%
  filter(cohort == 'cd_25') %>%
  summarise(count_id = length(unique(participant_id))) %>%
  pull(count_id)

cat('CD 2025: No. of participants at Pretest =', count_cd_25_pre, '\n')
cat('CD 2025: No. of participants at Posttest =', count_cd_25_post, '\n')
cat('CD 2025: Final no. of participants =', count_cd_25_final, '\n')

# - Assessing total exclusions --------------------------

# General overview #
# ============================ #

# Total number of participants who finished Pretest
pretest_total_count <- tibble(
  count = c(sum(count_kd_apr_24_pre, count_kd_sep_24_pre,
                count_kd_jan_25_pre), count_cd_25_pre),
  group = c('KD', 'CD')
  )
pretest_total_count

# Total number of participants who finished Posttest
posttest_total_count <- tibble(
  count = c(sum(count_kd_apr_24_post, count_kd_sep_24_post,
                count_kd_jan_25_post), count_cd_25_post),
  group = c('KD', 'CD')
)
posttest_total_count

# Total number of participants after filtering
final_total_count <- tibble(
  count = c(sum(count_kd_apr_24_final, count_kd_sep_24_final,
                count_kd_jan_25_final), count_cd_25_final),
  group = c('KD', 'CD'))
final_total_count

# Participant loss
participant_loss <- tibble(
  count = c(sum(count_kd_apr_24_pre, count_kd_sep_24_pre,
                                             count_kd_jan_25_pre) -
              sum(count_kd_apr_24_final, count_kd_sep_24_final,
                  count_kd_jan_25_final),
            count_cd_25_pre - count_cd_25_final
            ),
  group = c('KD', 'CD')
)
participant_loss

# Detailed accounting #
# ============================ #

# Unfiltered datasets
unfiltered_datasets <- bind_rows(
  kd_apr_24_combined,
  kd_sep_24_combined,
  kd_jan_25_combined,
  cd_25_combined
) %>% filter(!(participant_id == '' | participant_id == ' '))

# Finished both Pretest and Posttest
finished <- unfiltered_datasets %>%
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       finished_pre & finished_post) |
      (cohort == 'cd_25' &
         finished_pre & finished_post)
  ) %>%
  count(group) %>%
  ungroup() %>%
  as_tibble()

# Measured ketones
measured_ketones <- unfiltered_datasets %>%
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       finished_pre & finished_post) |
      (cohort == 'cd_25' &
         finished_pre & finished_post)
  ) %>%
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       measured_ketones_pre & measured_ketones_post) |
      (cohort == 'cd_25' & measured_ketones_post)
  ) %>% 
  count(group) %>%
  ungroup() %>%
  as_tibble()

# BMI filter
BMI_filter <- unfiltered_datasets %>%
  # Finished, measured filter
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       finished_pre & finished_post) |
      (cohort == 'cd_25' &
         finished_pre & finished_post)
  ) %>%
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       measured_ketones_pre & measured_ketones_post) |
      (cohort == 'cd_25' & measured_ketones_post)
  ) %>% 
  # BMI filter
  filter(BMI_pre >= 18.5 & BMI_pre < 30) %>%
  count(group) %>%
  ungroup() %>%
  as_tibble()

# Sex and age filter
sexage_filter <- unfiltered_datasets %>%
  # Finished, measured filter
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       finished_pre & finished_post) |
      (cohort == 'cd_25' &
         finished_pre & finished_post)
  ) %>%
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       measured_ketones_pre & measured_ketones_post) |
      (cohort == 'cd_25' & measured_ketones_post)
  ) %>% 
  # BMI filter
  filter(BMI_pre >= 18.5 & BMI_pre < 30) %>%
  # Sex filter 
  filter(sex == 'female', age >= 35 & age <= 45) %>%
  count(group) %>%
  ungroup() %>%
  as_tibble()

# Diet filters
diet_filters <- unfiltered_datasets %>%
  # Finished, measured filter
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       finished_pre & finished_post) |
      (cohort == 'cd_25' &
         finished_pre & finished_post)
  ) %>%
  filter(
    (cohort %in% c('apr_24', 'sep_24', 'jan_25') &
       measured_ketones_pre & measured_ketones_post) |
      (cohort == 'cd_25' & measured_ketones_post)
  ) %>% 
  # BMI filter
  filter(BMI_pre >= 18.5 & BMI_pre < 30) %>%
  # Sex filter 
  filter(sex == 'female', age >= 35 & age <= 45) %>%
  # Diet exclusions
  filter(!current_diet %in% c('Paleo_dieet', 'Low_carb_dieet')) %>%
  filter(!(current_diet_text == 'glutenvrij, en lactose en koolhydraat arm' | current_diet_text == 'Orthomoleculair paleo' | current_diet_text == 'Oersterk')) %>%
  count(group) %>%
  ungroup() %>%
  as_tibble()

# Ketone missingness filter 
ketones_missing_out <- dataset_start %>%
  # Diet exclusions only - other filters already implemented (incl. ketone missingness)
  filter(!current_diet %in% c('Paleo_dieet', 'Low_carb_dieet')) %>%
  filter(!(current_diet_text == 'glutenvrij, en lactose en koolhydraat arm' | current_diet_text == 'Orthomoleculair paleo' | current_diet_text == 'Oersterk')) %>%
  count(group) %>%
  ungroup() %>%
  as_tibble()

# Ketone exclusion filter
ketones_out_filter <- dataset_ketout %>%
  count(group) %>%
  ungroup()

cat('No. of participants before filtering: \n')
finished
cat('No. of participants who measured ketones: \n')
measured_ketones
cat('No. of participants after BMI filter: \n')
BMI_filter
cat('No. of participants after sex and age filter: \n')
sexage_filter
cat('No. of participants after diet filter: \n')
diet_filters
cat('No. of participants after ketone missingness filter: \n')
ketones_missing_out
cat('No. of participants after ketone exclusion filter: \n')
ketones_out_filter

# - Descriptives CD, KD ketones ---------------------------

# Aggregate at the subject level
ketone_subjects <- dataset_diettext %>% 
  select(participant_id, group, starts_with("ketones_post_int_")) %>% 
  rowwise() %>% 
  mutate(
    mean_ketone = mean(c_across(starts_with("ketones_post_int_")), na.rm = TRUE),
    sd_ketone = sd(c_across(starts_with("ketones_post_int_")), na.rm = TRUE),
    n_ketone = sum(!is.na(c_across(starts_with("ketones_post_int_")))),
    se_ketone = sd_ketone / sqrt(n_ketone)
  ) %>% 
  ungroup() %>% 
  # Replace any NaN means (all-NA rows) with NA_real_
  mutate(mean_ketone = ifelse(is.nan(mean_ketone), NA_real_, mean_ketone),
         sd_ketone = ifelse(is.nan(sd_ketone), NA_real_, sd_ketone),
         se_ketone = ifelse(is.nan(se_ketone), NA_real_, se_ketone))
glimpse(ketone_subjects)

pal <- c(CD = '#600985', KD = '#004445')   # palette

# Plot subject-level means
ggplot(ketone_subjects, aes(x = group, y = mean_ketone, colour = group)) +
  geom_pointrange(aes(
    ymin = mean_ketone - se_ketone,
    ymax = mean_ketone + se_ketone
  ),
  position = position_jitter(width = 0.2, height = 0),
  alpha = 0.6,
  ) +
  scale_colour_manual(values = pal) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black") +
  labs(
    x = "Group",
    y = "Ketone mm/L",
    title = "Subject-average ketone values per group",
    colour = "Group"
  ) +
  theme_minimal()

# Ketone descriptives
ketone_descriptives <- ketone_subjects %>%
  group_by(group) %>% 
  summarise(
    n = sum(!is.na(mean_ketone)),
    ketone_mean = mean(mean_ketone, na.rm = TRUE),
    sd_ketone = sd(mean_ketone,  na.rm = TRUE),
    min_ketone = min(mean_ketone, na.rm = TRUE),
    max_ketone = max(mean_ketone, na.rm = TRUE),
    median_ketone = median(mean_ketone, na.rm = TRUE)
  )
ketone_descriptives

# Check the number of measurements in ketosis
ketone_descriptives_2 <- dataset_diettext %>% 
  rowwise() %>% 
  mutate(
    ketosis = sum(c_across(starts_with('ketones_post_int_') ) >= 0.5, na.rm = TRUE)
  ) %>% 
  ungroup() %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_ketosis = median(ketosis, na.rm = TRUE),
    IQR_ketosis = IQR(ketosis, na.rm = TRUE)
  ) %>%
  ungroup()
ketone_descriptives_2

# Check CD participants in ketosis
cd_ketosis_check <- dataset_diettext %>%
  filter(group == 'CD' & 
           rowSums(across(starts_with('ketones_post_int_'), ~ .x >= 0.5), na.rm = TRUE) > 0) %>%
  select(participant_id, starts_with('ketones_post_int_'))
view(cd_ketosis_check)

cd_ketosis_check_2 <- dataset_diettext %>% 
  rowwise() %>% 
  mutate(
    ketosis = sum(c_across(starts_with('ketones_post_int_') ) >= 0.5, na.rm = TRUE)
  ) %>% 
  ungroup() %>%
  select(participant_id, group, ketosis) %>%
  filter(group == 'CD' & ketosis > 0)
print(cd_ketosis_check_2, n = Inf)

inspect_cd_foodlog <- cd_25_post_CLEAN %>%
  select(participant_id, starts_with('foodlog_')) %>%
  filter(participant_id == 'R_295A0U8cU31POI9')

inspect_cd_ketosis <- dataset_diettext %>%
  select(participant_id, group, starts_with('ketones_post_int_')) %>%
  filter(participant_id == 'R_295A0U8cU31POI9')

# - Time series CD, KD ketones -----------------------------

# NOTE: These visualisations explore the sample before ketone-level filtering

keto_timeseries <- dataset_diettext %>%
  select(participant_id, group, starts_with('ketones_post_int_')) %>%
  rename_with(~ str_replace(.x, '^ketones_post_int_(\\d+)$', 'time_\\1'), 
              starts_with('ketones_post_int_')) %>%
  pivot_longer(
    cols = starts_with('time_'),
    names_to = 'timepoint',
    names_prefix = 'time_',
    values_to = 'ketone_level'
  ) %>%
  mutate(
    timepoint = as.integer(timepoint), # Convert timepoint to integer
    group = factor(group, levels = c('CD', 'KD'))  # Ensure group is a factor with specific levels
  )

# KD
ggplot(keto_timeseries %>% filter(group == 'KD'),
       aes(x = timepoint, y = ketone_level, colour = group)) +
  
  # Individual participant trajectories (faint)
  geom_line(aes(group = participant_id),
            alpha = 0.3, linewidth = 0.4) +
  
  # Group means
  stat_summary(fun = mean, geom = 'line', linewidth = 1.1) +
  stat_summary(fun = mean, geom = 'point', size = 2) +
  
  # ± SE error bars
  stat_summary(fun.data = ggplot2::mean_se, geom = 'errorbar',
               width = 0.15, linewidth = 0.6) +
  
  scale_colour_manual(values = pal) +
  labs(
    x = 'Time point',
    y = 'Ketone value (mmol/L)',
    colour = 'Group'
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'top')

# CD
ggplot(keto_timeseries %>% filter(group == 'CD'),
       aes(x = timepoint, y = ketone_level, colour = group)) +
  
  # Individual participant trajectories (faint)
  geom_line(aes(group = participant_id),
            alpha = 0.3, linewidth = 0.4) +
  
  # Group means
  stat_summary(fun = mean, geom = 'line', linewidth = 1.1) +
  stat_summary(fun = mean, geom = 'point', size = 2) +
  
  # ± SE error bars
  stat_summary(fun.data = ggplot2::mean_se, geom = 'errorbar',
               width = 0.15, linewidth = 0.6) +
  
  scale_colour_manual(values = pal) +
  labs(
    x = 'Time point',
    y = 'Ketone value (mmol/L)',
    colour = 'Group'
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'top')

# Individual regression lines
ggplot(keto_timeseries,
       aes(x = timepoint, y = ketone_level, colour = group)) +
  
  # one lm line per participant
  geom_smooth(aes(group = participant_id),
              method  = 'lm',
              se = FALSE,
              linewidth = 0.6,
              linetype = 'solid') +
  
  # keep your group colouring for the trajectories
  scale_colour_manual(values = alpha(pal, 0.6)) +
  labs(
    x = 'Time point',
    y = 'Ketone value (mmol/L)',
    colour = 'Group'
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'top')

# GLOBAL COLUMN CLEANING ------------------------------

dataset_clean <- dataset_ketout %>%
  
  # ASRS cleaning ---------------------------

rename_with(
  ~ str_replace(.x, '^Q([1-9]|1[0-8])_(pre|post)$', 'ASRS_\\1_\\2'),
  matches('^Q([1-9]|1[0-8])_(pre|post)$')
) %>%
  mutate(
    across(starts_with('ASRS_'), 
           ~ case_when(
             .x == 'Nooit' ~ 0, # Never
             .x == 'Zelden' ~ 1, # Rarely
             .x == 'Soms' ~ 2, # Sometimes
             .x == 'Vaak' ~ 3, # Often
             .x == 'Zeer vaak' ~ 4, # Very often
             .default = NA_real_  # Ensure it remains numeric
           )
    )
  ) %>%
  
  # BDI cleaning --------------------------

rename_with(
  ~ str_replace(.x, '^Q(4[1-9]|5\\d|6[01])_(pre|post)$', 'BDI_\\1_\\2'),
  matches('^Q(4[1-9]|5\\d|6[01])_(pre|post)$')
) %>%
  mutate(
    across(starts_with('BDI_'), ~ parse_number(.x)
    )) %>%
  
  # PSQI cleaning ---------------------------

rename(
  PSQI_1_pre = Q80_1_pre, PSQI_1_post = Q80_1_post,
  PSQI_2_pre = Q81_pre, PSQI_2_post = Q81_post,
  PSQI_3_pre = Q82_1_pre, PSQI_3_post = Q82_1_post,
  PSQI_4_pre = Q83_pre, PSQI_4_post = Q83_post
) %>%
  rename_with(
    ~ str_replace(.x, '^Q(2[3-9]|3[0-7])_(pre|post)$', 'PSQI_\\1_\\2'),
    matches('^Q(2[3-9]|3[0-7])_(pre|post)$')
  ) %>%
  # PSQI question 5 renaming
  rename(
    PSQI_5a_pre = PSQI_23_pre, PSQI_5a_post = PSQI_23_post,
    PSQI_5b_pre = PSQI_25_pre, PSQI_5b_post = PSQI_25_post,
    PSQI_5c_pre = PSQI_26_pre, PSQI_5c_post = PSQI_26_post,
    PSQI_5d_pre = PSQI_27_pre, PSQI_5d_post = PSQI_27_post,
    PSQI_5e_pre = PSQI_28_pre, PSQI_5e_post = PSQI_28_post,
    PSQI_5f_pre = PSQI_29_pre, PSQI_5f_post = PSQI_29_post,
    PSQI_5g_pre = PSQI_30_pre, PSQI_5g_post = PSQI_30_post,
    PSQI_5h_pre = PSQI_31_pre, PSQI_5h_post = PSQI_31_post,
    PSQI_5i_pre = PSQI_32_pre, PSQI_5i_post = PSQI_32_post
  ) %>%
  # PSQI 5j not usable for us (no generable scores) --> PSQI_33 drop
  select(-c(PSQI_33_pre, PSQI_33_post)) %>% 
  # PSQI remaining questions renaming (tongue-twister not intended); the columns
  # were reordered to match the order of the original questionnaire
  rename(
    PSQI_6_pre = PSQI_37_pre, PSQI_6_post = PSQI_37_post,
    PSQI_7_pre = PSQI_34_pre, PSQI_7_post = PSQI_34_post,
    PSQI_8_pre = PSQI_35_pre, PSQI_8_post = PSQI_35_post,
    PSQI_9_pre = PSQI_36_pre, PSQI_9_post = PSQI_36_post
  ) %>%
  mutate(
    across(c(PSQI_1_pre, PSQI_1_post, PSQI_3_pre, PSQI_3_post), 
         ~ lubridate::hm(.x)),
    across(c(PSQI_1_pre, PSQI_1_post),
           ~ {
             if_else(.x >= hours(7) & .x <= hours(12),
                     .x + hours(12),
                     .x)
           }
    ),
    across(c(PSQI_2_pre, PSQI_2_post, PSQI_4_pre, PSQI_4_post), 
         ~ as.numeric(.x)),
    across(c(PSQI_4_pre, PSQI_4_post), ~ if_else(.x > 15, .x / 10, .x)),
    across(c(PSQI_6_pre, PSQI_6_post), ~ case_when(
    .x == 'Zeer goed' ~ 0, # Very good
    .x == 'Redelijk goed' ~ 1, # Fairly good
    .x == 'Eerder slecht' ~ 2, # Fairly bad
    .x == 'Zeer slecht' ~ 3, # Very bad
    .default = NA_real_  
  )
  ),
  across(c(PSQI_2_pre, PSQI_2_post), ~ case_when(
    .x <= 15 ~ 0, 
    .x <= 30 ~ 1, 
    .x <= 60 ~ 2,
    TRUE ~ 3,
    is.na(.x) ~ NA_real_  
  )),
  # PSQI_4 intentionally skipped here because it is needed in two forms 
  # for the composite scores
  across(c(PSQI_5a_pre, PSQI_5a_post), ~ case_when(
    .x == 'Niet tijdens de afgelopen 2 weken' ~ 0, # Not during the past 2 weeks
    .x == 'Minder dan één maal per week' ~ 1, # Less than once a week
    .x == 'Één – of tweemaal per week' ~ 2, # 1-2 times a week
    .x == 'Drie – of meermaals per week' ~ 3, # 3 or more times a week
    .default = NA_real_
  )
  ),
  across(matches('^PSQI_[578]') & !matches('^PSQI_5a'), 
         # Though Q8 seems to have different answer options, in our Dutch version,
         # the answer options matched those of the other selected questions
         ~ case_when(
           .x %in% c('Niet tijdens de afgelopen 2 weken',
                     'Niet tijdens de 2 weken')  ~ 0, # Not during the past 2 weeks
           .x == 'Minder dan één maal per week' ~ 1, # Less than once a week
           .x == 'Één – of tweemaal per week' ~ 2, # 1-2 times a week
           .x == 'Drie – of meermaals per week' ~ 3, # 3 or more times a week
           .default = NA_real_
         )
  ),
  across(c(PSQI_9_pre, PSQI_9_post), 
         ~ case_when(
           .x == 'Geen enkel probleem' ~ 0, # No problem at all
           .x == 'Slechts een klein probleem' ~ 1, # Only a small problem
           .x == 'Eingszins een probleem' ~ 2, # Somewhat of a problem
           .x == 'Een heel groot probleem' ~ 3, # A very big problem
           .default = NA_real_  
         )
  )
  )


# COMPOSITE SCORES ----------------------------------

dataset_final <- dataset_clean %>%
  
  # ASRS composite scores ---------------------

# Using just the first 6 questions because they're the screener with a diagnostic 
# value
  mutate(
  ASRS_total_pre = rowSums(across(matches('^ASRS_([1-6])_pre$')), na.rm = TRUE),
  ASRS_total_post = rowSums(across(matches('ASRS_([1-6])_post$')), na.rm = TRUE)
  ) %>%
  
  # BDI composite scores ---------------------
  mutate(
  BDI_total_pre = rowSums(across(starts_with('BDI_') & ends_with('_pre')), na.rm = TRUE),
  BDI_total_post = rowSums(across(starts_with('BDI_') & ends_with('_post')), na.rm = TRUE)
  ) %>%
  
  # PSQI composite scores --------------------

  mutate(
  PSQI_comp_1_pre = PSQI_6_pre, PSQI_comp_1_post = PSQI_6_post,
  PSQI_comp_2_pre = {
    psqi_sum <- PSQI_2_pre + PSQI_5a_pre
    case_when(
      psqi_sum == 0 ~ 0,
      psqi_sum <= 2 ~ 1,
      psqi_sum <= 4 ~ 2,
      TRUE ~ 3,
      is.na(psqi_sum) ~ NA_real_
    )
    }, 
  PSQI_comp_2_post = {
    psqi_sum <- PSQI_2_post + PSQI_5a_post
    case_when(
      psqi_sum == 0 ~ 0,
      psqi_sum <= 2 ~ 1,
      psqi_sum <= 4 ~ 2,
      TRUE ~ 3,
      is.na(psqi_sum) ~ NA_real_ 
    )
    },
  PSQI_comp_3_pre = case_when(
    PSQI_4_pre >= 7 ~ 0, 
    PSQI_4_pre >= 6 ~ 1, 
    PSQI_4_pre >= 5 ~ 2, 
    TRUE ~ 3,
    is.na(PSQI_4_pre) ~ NA_real_
  ),
  PSQI_comp_3_post = case_when(
    PSQI_4_post >= 7 ~ 0, 
    PSQI_4_post >= 6 ~ 1,
    PSQI_4_post >= 5 ~ 2,  
    TRUE ~ 3, 
    is.na(PSQI_4_post) ~ NA_real_
  )) %>%
  rowwise() %>% # One row ⇒ one efficiency value
  mutate(
    PSQI_comp_4_pre = {
      # Time in bed as a positive PERIOD
      tib <- PSQI_3_pre - PSQI_1_pre
      tib <- if_else(tib < hours(0), tib + hours(24), tib)
      
      # Convert to numeric hours *after* the period is final
      tib_hrs <- time_length(as.duration(tib), 'hours')
      
      sleep_eff <- 100 * PSQI_4_pre / tib_hrs # Habitual sleep efficiency (%)
      ifelse(sleep_eff > 100, 100, sleep_eff) # Cap at 100%
    },
    
    PSQI_comp_4_post = {
      # Time in bed as a positive PERIOD
      tib <- PSQI_3_post - PSQI_1_post
      tib <- if_else(tib < hours(0), tib + hours(24), tib)
      
      # Convert to numeric hours *after* the period is final
      tib_hrs <- time_length(as.duration(tib), 'hours')
      
      sleep_eff <- 100 * PSQI_4_post / tib_hrs # Habitual sleep efficiency (%)
      ifelse(sleep_eff > 100, 100, sleep_eff) # Cap at 100%
    },
  ) %>% 
  ungroup() %>%
  mutate(across(c(PSQI_comp_4_pre, PSQI_comp_4_post), ~ case_when(
    .x >= 85 ~ 0,
    .x >= 75 ~ 1,
    .x >= 65 ~ 2,
    TRUE ~ 3,
    is.na(.x) ~ NA_real_
  ))
  ) %>%
  mutate(
    PSQI_comp_5_pre = { 
      psqi_sum <- rowSums(  
        across(matches('^PSQI_5[b-i]_pre$')),
        na.rm = TRUE
      )
      
      case_when(
        psqi_sum == 0 ~ 0,
        psqi_sum <= 9 ~ 1,
        psqi_sum <= 18 ~ 2,
        TRUE ~ 3,
        is.na(psqi_sum) ~ NA_real_  
      )
    },
    PSQI_comp_5_post = {
      psqi_sum <- rowSums(across(matches('^PSQI_5[b-i]_post$')), na.rm = TRUE)
      case_when(
        psqi_sum == 0 ~ 0,
        psqi_sum <= 9 ~ 1,
        psqi_sum <= 18 ~ 2,
        TRUE ~ 3,
        is.na(psqi_sum) ~ NA_real_ 
      )
    },
    PSQI_comp_6_pre = PSQI_7_pre, PSQI_comp_6_post = PSQI_7_post,
  ) %>%
  mutate(
    PSQI_comp_7_pre = {
      psqi_sum <- PSQI_8_pre + PSQI_9_pre
      case_when(
        psqi_sum == 0 ~ 0,
        psqi_sum <= 2 ~ 1,
        psqi_sum <= 4 ~ 2,
        TRUE ~ 3,
        is.na(psqi_sum) ~ NA_real_
      )
    },
    
    PSQI_comp_7_post = {
      psqi_sum <- PSQI_8_post + PSQI_9_post
      case_when(
        psqi_sum == 0 ~ 0,
        psqi_sum <= 2 ~ 1,
        psqi_sum <= 4 ~ 2,
        TRUE ~ 3,
        is.na(psqi_sum) ~ NA_real_
      )
    } 
  ) %>%
  
  # Total score
  mutate(PSQI_total_pre = rowSums(across(starts_with('PSQI_comp_') & ends_with('_pre')), 
                                  na.rm = TRUE),
         PSQI_total_post = rowSums(across(starts_with('PSQI_comp_') & ends_with('_post')), 
                                   na.rm = TRUE)
         ) %>%
  
  # Final cleaning ---------------------------

  select(-matches('^SC\\d+'))

  
# - Check questionnaire ranges -------------------

dataset_final_minmax <- dataset_final %>%
  summarise(across(
    c(ends_with('_total_post'), ends_with('_total_pre')),
    list(min = ~min(.x, na.rm = TRUE),
         max = ~max(.x, na.rm = TRUE)),
    .names = '{.col}_{.fn}'
  ))
print(dataset_final_minmax, width = Inf)

# - Clean time series CD, KD ketones -----------------------------

keto_timeseries_clean <- dataset_final %>%
  dplyr::select(participant_id, group, starts_with('ketones_post_int_')) %>%
  rename_with(~ str_replace(.x, '^ketones_post_int_(\\d+)$', 'time_\\1'), 
              starts_with('ketones_post_int_')) %>%
  pivot_longer(
    cols = starts_with('time_'),
    names_to = 'timepoint',
    names_prefix = 'time_',
    values_to = 'ketone_level'
  ) %>%
  mutate(
    timepoint = as.integer(timepoint), # Convert timepoint to integer
    group = factor(group, levels = c('CD', 'KD'))  # Ensure group is a factor with specific levels
  )

# KD
gg_clean_timeseries_kd <- ggplot(keto_timeseries_clean %>% filter(group == 'KD'),
       aes(x = timepoint, y = ketone_level, colour = group)) +
  
  # Individual participant trajectories (faint)
  geom_line(aes(group = participant_id),
            alpha = 0.3, linewidth = 0.4) +
  
  # Group means
  stat_summary(fun = mean, geom = 'line', linewidth = 1.1) +
  stat_summary(fun = mean, geom = 'point', size = 2) +
  
  # ± SE error bars
  stat_summary(fun.data = ggplot2::mean_se, geom = 'errorbar',
               width = 0.15, linewidth = 0.6) +
  
  scale_colour_manual(values = pal, labels = c('Ketogenic Diet')) +
  labs(
    x = 'Time point',
    y = 'Ketone value (mmol/L)',
    colour = 'Group'
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 16)) +
  theme_apa() +
  theme(legend.position = 'top')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'clean_timeseries_kd.pdf'),
  plot = gg_clean_timeseries_kd,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# CD
gg_clean_timeseries_cd <- ggplot(keto_timeseries_clean %>% filter(group == 'CD'),
       aes(x = timepoint, y = ketone_level, colour = group)) +
  
  # Individual participant trajectories (faint)
  geom_line(aes(group = participant_id),
            alpha = 0.3, linewidth = 0.4) +
  
  # Group means
  stat_summary(fun = mean, geom = 'line', linewidth = 1.1) +
  stat_summary(fun = mean, geom = 'point', size = 2) +
  
  # ± SE error bars
  stat_summary(fun.data = ggplot2::mean_se, geom = 'errorbar',
               width = 0.15, linewidth = 0.6) +
  
  scale_colour_manual(values = pal, labels = c('Clean Diet')) +
  labs(
    x = 'Time point',
    y = 'Ketone value (mmol/L)',
    colour = 'Group'
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 16)) +
  theme_apa() +
  theme(legend.position = 'top')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'clean_timeseries_cd.pdf'),
  plot = gg_clean_timeseries_cd,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Individual regression lines
ggplot(keto_timeseries,
       aes(x = timepoint, y = ketone_level, colour = group)) +
  
  # one lm line per participant
  geom_smooth(aes(group = participant_id),
              method  = 'lm',
              se = FALSE,
              linewidth = 0.6,
              linetype = 'solid') +
  
  # keep your group colouring for the trajectories
  scale_colour_manual(values = alpha(pal, 0.6)) +
  labs(
    x = 'Time point',
    y = 'Ketone value (mmol/L)',
    colour = 'Group'
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'black') +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'top')

# SAVE DATA ---------------------

# Save as RDS
saveRDS(dataset_final, 
        here('UL-KD001_data', 'LBI_clean', 'UL-KD001_questionnaires_pseudo.rds'))

# Save as CSV
write_csv(dataset_final, 
          here('UL-KD001_data', 'LBI_clean', 'UL-KD001_questionnaires_pseudo.csv'))