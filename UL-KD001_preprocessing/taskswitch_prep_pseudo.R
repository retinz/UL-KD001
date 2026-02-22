# DESCRIPTION 
# ===================== #

# Preprocessing task-switching data for UL-KD001.

# LIBRARIES #
# ===================== #

library(here)
library(ggplot2)
library(tidyverse)

# RAW DATA ------------------------

# Load raw data
data_raw <- read.csv(here('UL-KD001_data', 'LBI_raw', 'raw_pseudoanonym',
                          'cognitive_data_20250417.csv'), 
                     header = TRUE)

# CLEAN DATA --------------------------

# Specify which columns are relevant
relevant_columns <- c('participant_id', 'group_id', 'session', 'date_startdate',
                      'date_starttime', 'block_count', 'count_experimental_trial_sequence',
                      'cue_color', 'task', 'congruence', 'digit', 'response',
                      'correct_response', 'response_time', 'correct')

# - Trim dataset ------------------------

data_trim <- data_raw %>% 
  filter(practice == 'no') %>% 
  dplyr::select(all_of(relevant_columns)) %>%
  distinct() %>%
  mutate(
    error = ifelse(correct == 1, FALSE, TRUE)
  ) %>%
  dplyr::select(-c(correct, response, correct_response, 
            date_starttime))
glimpse(data_trim)

# - Filter participants -------------------------

# Load clean and ready questionnaires to extract meta-data
metadata_questionnaires <- readRDS(here('UL-KD001_data', 'LBI_clean', 
'UL-KD001_questionnaires_pseudo.rds')) %>%
  dplyr::select(participant_id, group, cohort)

data_select <- data_trim %>%
  right_join(metadata_questionnaires, by = 'participant_id') %>%
  relocate(group, cohort, .after = participant_id)
glimpse(data_select)

# Check dataset properties
check_trials <- data_select %>%
  group_by(participant_id, session) %>%
  summarise(
    max_trial_sequence = max(count_experimental_trial_sequence, na.rm = TRUE),
    trial_count = n()
  ) %>%
  ungroup() %>%
  arrange(participant_id, session)
print(check_trials, n = Inf)
unique(check_trials$max_trial_sequence)
unique(check_trials$trial_count)

# Keep only participants with both sessions
data_filter <- data_select %>%
  group_by(participant_id) %>%
  filter(n_distinct(session) == 2 & session %in% c(1, 2)) %>%
  ungroup()
glimpse(data_filter)

# Check dataset properties again
check_trials_2 <- data_filter %>%
  group_by(participant_id, session) %>%
  summarise(
    trial_count = n(),
    blocks = n_distinct(block_count)
  ) %>%
  ungroup() %>%
  arrange(participant_id, session)
print(check_trials_2, n = Inf)
  
# - Trial adjustments ---------------------------

data_chopped <- data_filter %>%
  group_by(participant_id, session, block_count) %>%
  arrange(count_experimental_trial_sequence, .by_group = TRUE) %>%
  mutate(
    # Add trials 
    trial = row_number(),
    # Create transition columns
    task_transition = ifelse(task == lag(task), 'repeat', 'switch'),
    digit_transition = ifelse(digit == lag(digit), 'repeat', 'switch'),
    cue_transition = ifelse(cue_color == lag(cue_color), 'repeat', 'switch'),
    # Identify post-error trials
    post_error = lag(error)
  ) %>%
  ungroup() %>%
  # Remove first trial from each block
  filter(!trial == 1,
         # Remove repeating digit trials
         !digit_transition == 'repeat',
         # Remove post-error trials
         post_error == FALSE) %>%
  # Drop unnecessary columns 
  dplyr::select(-group_id) %>%
  rename(trial_sequence = count_experimental_trial_sequence)
glimpse(data_chopped)

# - Check no. of participants ----------------------------

participants_taskswitch <- data_chopped %>% 
  group_by(group) %>%
  summarise(n_participants = n_distinct(participant_id))
participants_taskswitch

# - Check trial exclusions -------------------

# Digit repetitions
digit_repetitions <- data_filter %>%
  group_by(participant_id, session, block_count) %>%
  mutate(
    # Add trials 
    trial = row_number(),
    # Create transition columns
    task_transition = ifelse(task == lag(task), 'repeat', 'switch'),
    digit_transition = ifelse(digit == lag(digit), 'repeat', 'switch'),
    cue_transition = ifelse(cue_color == lag(cue_color), 'repeat', 'switch')
  ) %>%
  ungroup() %>%
  filter(!trial == 1) %>%
  summarise(digit_rep = mean(digit_transition == 'repeat', na.rm = TRUE)) %>%
  dplyr::pull(digit_rep)
digit_repetitions

# Post error trials
post_errors <- data_filter %>%
  group_by(participant_id, session, block_count) %>%
  arrange(count_experimental_trial_sequence) %>%
  mutate(
    # Identify post-error trials
    post_error = lag(error),
    # Add trials 
    trial = row_number(),
    # Create transition columns
    task_transition = ifelse(task == lag(task), 'repeat', 'switch'),
    digit_transition = ifelse(digit == lag(digit), 'repeat', 'switch'),
    cue_transition = ifelse(cue_color == lag(cue_color), 'repeat', 'switch')
  ) %>%
  ungroup() %>%
  # Remove first trial from each block
  filter(!trial == 1) %>%
  summarise(post_errors = mean(post_error == TRUE, na.rm = TRUE)) %>%
  dplyr::pull(post_errors)
post_errors

# SAVE DATA -----------------------

# Save as RDS
saveRDS(data_chopped, here('UL-KD001_data', 'LBI_clean', 'UL-KD001_taskswitch_pseudo.rds'))

# Save as CSV
write_csv(data_chopped, here('UL-KD001_data', 'LBI_clean', 'UL-KD001_taskswitch_pseudo.csv'))

