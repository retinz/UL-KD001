# DESCRIPTION #
# ===================== #

# The core analysis of UL-KD001 data. This script contains the most relevant analyses 
# for the study.

# TASKS #
# ===================== #

# - task-switching re-analysis:
# --- post-error marking (DONE) and exclusions (DONE)
# --- outlier exclusion at the trial level (DONE)
# --- same model specification (also log-link function) (DONE)
# --- task transition * cue transition interaction (DONE)
# --- possible task transition * congruence interaction (DONE)
# --- cross-validation (DONE)
# --- posthocs + plotting (DONE)
# --- RT and ACC baselines (DONE)
# --- p-value adjustments

# NOTES #
# ===================== #

# Task switching RT:
# - none of the fitted lme4 converges (properly)
# - glmmTMB solves the convergence issue
# - cue transition * task transition interaction rank deficient and automatically
# dropped by lme4 and glmmTMB (tested also with treatment contrasts - doesn't help)
# - taskswitch_1a7_rt the best RT model
# - taskswitch_2_er the best ACC model (cross-validated)

# LIBRARIES #
# ===================== #

# The order of libraries somewhat matters!

library(here)
library(ARTool)
library(WRS)
library(WRS2)
library(afex)
library(DHARMa)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(effectsize)
library(MatchIt)
library(cobalt)
library(correlation)
library(estimatr)
library(glmmTMB)
source(here('UL-KD001_analysis', 'analysis_helpers.R'))

# The WRS package doesn't come from CRAN but can be downloaded from: https://github.com/nicebread/WRS. 

# SESSION SETTINGS #
# ===================== #

# Set contrasts
options(contrasts = c('contr.sum', 'contr.poly'))

# Palette
# - Should be differentiable even by people with a colour-perception disorder
pal <- c(CD = '#a2bff4', KD = '#80c680')

# Dodge
dodge_tsmm <- 0.2
dodge_w <- 0.6
dodge_w_BMI <- 0.2

# Default directory for plots
plot_directory <- here('UL-KD001_analysis', 'plots')
if (!dir.exists(plot_directory)) {
  dir.create(plot_directory, recursive = TRUE)
}

# Default directory for tables
table_directory <- here('UL-KD001_analysis', 'tables')
if (!dir.exists(table_directory)) {
  dir.create(table_directory, recursive = TRUE)
}

# Data path
data_path <- here('UL-KD001_data', 'LBI_clean')

# Detect number of logical processors for parallel processing
n_cores = max(1L, as.integer(floor(parallel::detectCores())))

# LOAD DATA -------------------------

# QUESTIONNAIRES #
# ================================ #

# Load questionnaire RDS
questionnaires_data <- readRDS(file.path(data_path, 'UL-KD001_questionnaires_pseudo.rds'))

# Adjust columns
questionnaires_ready <- questionnaires_data %>%
  mutate(
    group = factor(group, levels = c('CD', 'KD')),
    cohort = factor(cohort, levels = c('apr_24', 'sep_24', 'jan_25', 'cd_25'))
  ) %>%
  mutate(ASRS_IA_pre = rowSums(across(matches('^ASRS_([1-4])_pre$')), na.rm = TRUE),
         ASRS_IA_post = rowSums(across(matches('^ASRS_([1-4])_post$')), na.rm = TRUE),
         ASRS_H_pre = rowSums(across(matches('^ASRS_([5-6])_pre$')), na.rm = TRUE),
         ASRS_H_post = rowSums(across(matches('^ASRS_([5-6])_post$')), na.rm = TRUE)
  ) %>%
  dplyr::select(participant_id, group, cohort, age, BMI_pre, BMI_post,
         ASRS_total_pre, ASRS_total_post, BDI_total_pre, BDI_total_post, PSQI_total_pre,
         PSQI_total_post, ASRS_IA_pre, ASRS_IA_post, ASRS_H_pre, ASRS_H_post)

# TASK-SWITCHING #
# ================================ #

# Load task-switching RDS 
taskswitch_data <- readRDS(file.path(data_path, 'UL-KD001_taskswitch_pseudo.rds'))

# Adjust columns
taskswitch_ready <- taskswitch_data %>%
  mutate(
    group = factor(group, levels = c('CD', 'KD')),
    cohort = factor(cohort, levels = c('apr_24', 'sep_24', 'jan_25', 'cd_25')),
    session = factor(session, levels = c(1, 2)),
    task = factor(task, levels = c('p', 'm'), labels = c('parity', 'magnitude')),
    congruence = factor(congruence, levels = c('c', 'i'), 
                        labels = c('congruent', 'incongruent')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch')),
    cue_transition = factor(cue_transition, levels = c('repeat' , 'switch')),
    trial_sequence_factor = glmmTMB::numFactor(trial_sequence)
    )

# QUESTIONNAIRES + KETONES #
# ================================ #

# Adjust columns
ketones_tibble <- questionnaires_data %>%
  mutate(
    group = factor(group, levels = c('CD', 'KD')),
    cohort = factor(cohort, levels = c('apr_24', 'sep_24', 'jan_25', 'cd_25')),
  ) %>%
  dplyr::select(participant_id, group, cohort, ketones_pre, starts_with('ketones_post_int_'))
glimpse(ketones_tibble)

# QUESTIONNAIRES PREPARATION ------------------

glimpse(questionnaires_ready)

# Reshape tibble
questionnaires_long <- questionnaires_ready %>%
  pivot_longer(cols = matches('_(pre|post)$'),
               names_to = c('measure', 'session'),
               names_pattern = '(.*)_(pre|post)$',
               values_to = 'value') %>%
  mutate(session = if_else(session == 'pre', 1, 2),
         measure = str_remove(measure, '_total$')) %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  ) %>%
  mutate(session = factor(session, levels = c(1, 2)))
glimpse(questionnaires_long)

# Get number of participants
participants_n_questionnaires <- questionnaires_long %>%
  group_by(group) %>%
  summarise(n_participants = n_distinct(participant_id)) %>%
  ungroup()
participants_n_questionnaires

# Get robust descriptives
questionnaires_trwin <- trim_winsorise(questionnaires_long, 
                                       dv_colnames = c('ASRS', 
                                                       'BDI',
                                                       'PSQI',
                                                       'BMI',
                                                       'ASRS_IA',
                                                       'ASRS_H'),
                                       within = 'session',
                                       between = 'group',
                                       id = 'participant_id',
                                       tr = 0.2)
print(questionnaires_trwin, width = Inf)

# DEMOGRAPHICS ---------------------------

glimpse(questionnaires_long)

# - Assess age -------------------

# Age distributions #
# ----------------------- #

ggplot(questionnaires_long, aes(x = age, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Age',
       y = 'Density',
       fill = 'Group',
       title = 'Age Distributions by Group and Session') +
  theme_minimal()

# Age descriptives #
# ---------------------- #

age_descriptives <- questionnaires_long %>%
  group_by(group) %>%
  summarise(mean_age = mean(age, na.rm = TRUE),
            sd_age = sd(age, na.rm = TRUE),
            median_age = median(age, na.rm = TRUE),
            IQR_age = IQR(age, na.rm = TRUE)) %>%
  ungroup()
age_descriptives

# Test difference #
# ---------------------- #

age_diff <- wilcox.test(
  x = questionnaires_long %>% filter(group == 'CD', 
                                     session == 1) %>% dplyr::pull(age),
  y = questionnaires_long %>% filter(group == 'KD',
                                     session == 1) %>% dplyr::pull(age),
                        exact = FALSE)
age_diff

# - Assess BMI ---------------

# BMI distributions #
# ----------------------- #

ggplot(questionnaires_long, aes(x = BMI, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'BMI',
       y = 'Density',
       fill = 'Group',
       title = 'BMI Distributions by Group and Session') +
  theme_minimal()

# BMI descriptives #
# ---------------------- #

BMI_descriptives <- questionnaires_long %>%
  group_by(group, session) %>%
  summarise(mean_BMI = mean(BMI, na.rm = TRUE),
            sd_BMI = sd(BMI, na.rm = TRUE),
            median_BMI = median(BMI, na.rm = TRUE),
            IQR_BMI = IQR(BMI, na.rm = TRUE)) %>%
  ungroup()
BMI_descriptives

# Testing differences #
# =========================== #

# Baseline #
# ---------------- #

# Yuen
BMI_baseline_test <- yuen(BMI ~ group, 
                          data = questionnaires_long %>% filter(session == 1),
                          tr = 0.2
                          )
BMI_baseline_test

# Regular ttest
BMI_baseline_rtest <- t.test(BMI ~ group, 
                          data = questionnaires_long %>% filter(session == 1))
BMI_baseline_rtest

# Interaction #
# ---------------- #

# 2x2 interaction (TRIMMED) 
BMI_tranova <- WRS2::bwtrim(BMI ~ group * session,
                            id = participant_id,
                            data = questionnaires_long,
                            tr = 0.2)
BMI_tranova

# 2x2 interaction (REGULAR)
BMI_anova <- aov_4(BMI ~ group * session + (session | participant_id),
                   data = questionnaires_long, 
                   anova_table = list(es = 'pes'))
BMI_anova

# Visualise interaction
afex_plot(
  BMI_anova,
  x = 'session',
  trace = 'group',
  error = 'within',
  mapping = c('fill', 'colour'),
  dodge = 0.30, 
  point_arg = list(size = 3),
  data_arg  = list(           
    size = 2,
    alpha = 0.6,
    position = position_jitterdodge(
      jitter.width = 0.2,     
      jitter.height = 0,
      dodge.width = 0.30  
    )
  )) +
  scale_fill_manual(values  = pal, name = NULL) +
  scale_colour_manual(values = pal, name = NULL) +
  theme_minimal()

# Visualise interaction 2
questionnaires_long_afex <- questionnaires_long %>%
  mutate(session_afex = factor(session, levels = c(1,2),
                               labels = c('X1', 'X2')))

BMI_afex <- afex_plot(
  BMI_anova,
  x = 'session',
  trace = 'group',
  error = 'within',
  mapping = c('fill', 'colour'),
  dodge = dodge_w_BMI,
  data_plot = FALSE,
  point_arg = list(size = 4)
) +
  # Participant trajectories
  geom_line(
    data = dplyr::filter(questionnaires_long_afex, group == 'KD'),
    aes(session_afex, BMI, group = participant_id, colour = group, 
        fill = NULL),
    position = position_nudge(x = +dodge_w_BMI / 4),
    alpha = 0.35, linewidth = 0.4
  ) +
  geom_line(
    data = dplyr::filter(questionnaires_long_afex, group == 'CD'),
    aes(session_afex, BMI, group = participant_id, colour = group,
        fill = NULL),
    position = position_nudge(x = -dodge_w_BMI / 4),
    alpha = 0.35, linewidth = 0.4
  ) +
  # Raw points: dodge only by KD/CD; no jitter
  geom_point(
    data = questionnaires_long_afex,
    aes(session_afex, BMI, colour = group, fill = group, group = group),
    position = position_jitterdodge(
      jitter.width = 0.06,
      jitter.height = 0,
      dodge.width = dodge_w_BMI,
      seed = 1
    ),
    size = 2, alpha = 0.4
  ) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  guides(colour = guide_legend(title = 'Group'),
         fill = guide_legend(title = 'Group')) +
  labs(x = 'Session') +
  scale_x_discrete(labels = c('Pretest','Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BMI_afex.pdf'),
  plot = BMI_afex,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Robust plot #
# ---------------- #

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_BMI_with <- ggplot(questionnaires_long,
                       aes(x = session, y = BMI, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = questionnaires_trwin %>% 
                  filter(dv == 'BMI'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(x = NULL, y = 'BMI', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0.1, 0.1))) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BMI_with.pdf'),
  plot = gg_BMI_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Post-hoc tests #
# ---------------- #

BMI_emm <- emmeans(BMI_anova, ~ group * session)
BMI_emm

# Difference in differences
BMI_did <- contrast(BMI_emm, interaction = c('pairwise', 'pairwise'))
BMI_did

# Declines in groups
BMI_sess <- emmeans(BMI_anova, ~ session | group) %>%
  contrast('pairwise')
BMI_sess

# Robust post-hoc tests #
# ----------------------------- #

# Create ordered list for analysis
BMI_post_tr_list <- build_split_list(questionnaires_long, 
                                      'BMI', c('group', 'session'))

# Contrasts on trimmed means
# - Gives warnings for some reason; reducing the no. of bootstrap samples
# doesn't help
BMI_post_tr <- bwmcp(2, 2, BMI_post_tr_list, tr=0.2,con=0, nboot=599)
BMI_post_tr

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
BMI_post_tr_int <- bwimcpES(2, 2, BMI_post_tr_list, tr=0.2, CI=TRUE, alpha=0.05)
BMI_post_tr_int

# Session changes KD
BMI_session_KD <- yuend(questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                          dplyr::pull(BMI),
                         questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                          dplyr::pull(BMI),
                         tr = 0.2)
BMI_session_KD

# Session change CD
BMI_session_CD <- yuend(questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                          dplyr::pull(BMI),
                         questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                          dplyr::pull(BMI),
                         tr = 0.2)
BMI_session_CD

# Effect size of session changes
BMI_session_eff <- bw.es.B(2, 2, BMI_post_tr_list, 
                            tr = 0.2, POOL = FALSE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
BMI_session_eff

# TASK SWITCHING PREPARATION -----------------------------

# Check current dataset
glimpse(taskswitch_ready)

# Filter data #
# ----------------------- #

# Select participants to be excluded due to too high an error rate
taskswitch_exclude <- taskswitch_ready %>%
  # Aggregate data
  group_by(group, session, participant_id) %>%
  summarise(error_rate = mean(error, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(error_rate >= 0.3)

# Clean the RT trials
taskswitch_crispy_rt <- taskswitch_ready %>%
  filter(
    # Exclude error trials
    error == FALSE, 
    # Filter out too fast and too slow trials
    response_time > 300 & response_time < 3000,
    # Filter out participants with too high an error rate
    !(participant_id %in% taskswitch_exclude$participant_id),
    # Exclude post error trials
    post_error == FALSE
  )

# Cleaning for ER analysis
taskswitch_crispy_er <- taskswitch_ready %>%
  filter(
    # Filter out too fast and too slow trials
    response_time > 300 & response_time < 3000,
    # Filter out participants with too high an error rate
    !(participant_id %in% taskswitch_exclude$participant_id),
    # Exclude post error trials
    post_error == FALSE
  ) %>%
  mutate(response_correct = ifelse(error == FALSE, 1, 0))

# Examine filtering results #
# -------------------------------- #

# Check number of participants after filtering
taskswitch_participants <- taskswitch_crispy_rt %>%
  group_by(group) %>%
  summarise(count = n_distinct(participant_id)) %>%
  ungroup()
taskswitch_participants

# Excluded RT trials because being too fast or too slow (RT analyses)
task_switch_speed <- taskswitch_ready %>%
  filter(error == FALSE,
         !(participant_id %in% taskswitch_exclude$participant_id)) %>%
  summarise(speed_off_prop = mean(response_time < 300 | response_time > 3000)) %>%
  dplyr::pull(speed_off_prop)
task_switch_speed

# Excluded RT trials because being too fast or too slow (ACC analyses)
task_switch_speed_acc <- taskswitch_ready %>%
  filter(!(participant_id %in% taskswitch_exclude$participant_id)) %>%
  summarise(speed_off_prop = mean(response_time < 300 | response_time > 3000)) %>%
  dplyr::pull(speed_off_prop)
task_switch_speed_acc

# - Variable trial filtering -----------------

# Assess subject-level RT distributions
inspect_subject_dist(
  taskswitch_crispy_rt,
  n = 15,
  pid_col = 'participant_id',
  col_dist = 'response_time',
  group_col = 'group',
  session_col = 'session',
  wrap = c('task_transition', 'congruence', 'cue_transition'),
  pal = pal,
  robust = 'MAD',
  mad_mult = 2.2,
  seed = NULL
)

# Detect and remove outliers (RT)
taskswitch_rt <- taskswitch_crispy_rt %>%
  group_by(participant_id, group, session, task_transition, congruence, cue_transition) %>%
  filter({
    m <- median(response_time, na.rm = TRUE)
    d <- mad(response_time, na.rm = TRUE)
    between(response_time, m - d * 2.2, m + d * 2.2)
  }) %>%
  ungroup()

# No outlier removal for ER analyses
taskswitch_er <- taskswitch_crispy_er

# Examine filtering results #
# -------------------------------- #

# Excluded RT trials because being too fast or too slow (RT analyses)
task_switch_outliers_rt <- taskswitch_crispy_rt %>%
  group_by(participant_id, group, session, task_transition, congruence, cue_transition) %>%
  mutate(outlier = {
    m <- median(response_time, na.rm = TRUE)
    d <- mad(response_time, na.rm = TRUE)
    !between(response_time, m - d * 2.2, m + d * 2.2)
  }) %>%
  ungroup() %>%
  summarise(outlier_prop = mean(outlier, na.rm = TRUE)) %>%
  dplyr::pull(outlier_prop)
task_switch_outliers_rt

# TASK SWITCHING RT ANALYSIS -----------
# - 4-way RE*+ ---------

# Model
# - Doesn't converge
taskswitch_1_rt <- lme4::glmer(
  response_time ~ group * session * task_transition * congruence + 
    task_transition * cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_1_rt)

# Different package
taskswitch_1a1_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence + 
    task_transition * cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_1a1_rt)

# Different package - no deficient ranks
taskswitch_1a2_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence
  + cue_transition +
  (task_transition * congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_1a2_rt)

# Different package - gamma family
# - Doesn't converge
taskswitch_1a3_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence
  + cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = Gamma(link = 'log'))
summary(taskswitch_1a3_rt)

# DIAGNOSTICS: taskswitch_1a2_rt #
# -------------------------------------- #

diagnose_model(taskswitch_1a2_rt)

# - 4-way RE*+ variance --------------

# Get predicted values
taskswitch_rt_pred <- taskswitch_rt %>%
  mutate(fitted_1a2_rt = predict(taskswitch_1a2_rt))
glimpse(taskswitch_rt_pred)

# Model
# - doesn't converge
taskswitch_1a4_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence
  + cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  dispformula = ~ poly(fitted_1a2_rt, 2) + session + group * task_transition,
  data = taskswitch_rt_pred, family = gaussian(link = 'log'))
summary(taskswitch_1a4_rt)

# Model
taskswitch_1a5_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence
  + cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  dispformula = ~ session * group + task_transition,
  data = taskswitch_rt_pred, family = gaussian(link = 'log'))
summary(taskswitch_1a5_rt)

# Model
taskswitch_1a6_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence
  + cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  dispformula = ~ session + group + task_transition,
  data = taskswitch_rt_pred, family = gaussian(link = 'log'),
  control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(taskswitch_1a6_rt)

# Diagnostics
diagnose_model(taskswitch_1a6_rt)

# Model
taskswitch_1a7_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence
  + cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  dispformula = ~ poly(fitted_1a2_rt, 2) + session + group + task_transition,
  data = taskswitch_rt_pred, family = gaussian(link = 'log'),
  control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(taskswitch_1a7_rt)

# Diagnostics
diagnose_model(taskswitch_1a7_rt)

# Cross-validation #
# ----------------------- #

# Cluster-based
cv_taskswitch_1a7_rt <- cv(taskswitch_1a7_rt, 
                         k = 5, 
                         clusterVariables = 'participant_id', 
                         ncores = n_cores,
                         reps = 10)
summary(cv_taskswitch_1a7_rt)
plot(cv_taskswitch_1a7_rt)

# - Autocorrelation check -------------

# Get model residuals
set_simcodes(taskswitch_1a7_rt, val = 'fix')
taskswitch_1a7_rt_res <- simulateResiduals(taskswitch_1a7_rt, plot = FALSE)

# Temporal autocorrelation per subject
taskswitch_1a7_rt_autocorr <- taskswitch_rt_pred %>%
  dplyr::select(participant_id, group, session, trial_sequence) %>%
  mutate(residual = taskswitch_1a7_rt_res$scaledResiduals) %>%
  group_by(participant_id, group, session) %>%
  summarise(
    dw = lmtest::dwtest(residual ~ 1, order.by = trial_sequence)$statistic,
    p_value = round(lmtest::dwtest(residual ~ 1, order.by = trial_sequence)$p.value, 3)
  ) %>%
  ungroup()

# Autocorrelation descriptives
taskswitch_1a7_rt_autocorr_desc <- taskswitch_1a7_rt_autocorr %>%
  group_by(group, session) %>%
  summarise(autocorr_prop = mean(dw < 1.5 | dw > 2.5),
            mean_dw = mean(dw, na.rm = TRUE),
            sd_dw = sd(dw, na.rm = TRUE)) %>%
  ungroup()
taskswitch_1a7_rt_autocorr_desc

# Plot
max_lag <- 40
acf_mm_1a7_rt <- taskswitch_rt_pred %>%
  mutate(residual = taskswitch_1a7_rt_res$scaledResiduals) %>%     
  group_by(participant_id, group, session) %>% 
  # Ensure strict temporal order before acf calculation
  arrange(trial_sequence, .by_group = TRUE) %>%
  nest() %>%                                                     
  mutate(acf_tbl = map(data, ~{
    ac <- acf(.x$residual, lag.max = max_lag, plot = FALSE)
    tibble(lag = as.numeric(ac$lag[-1]),                        
           acf = as.numeric(ac$acf[-1]))
  })) %>%
  dplyr::select(-data) %>%                                       
  unnest(acf_tbl) %>%
  ungroup()

taskswitch_1a7_rt_summary <- acf_mm_1a7_rt %>%
  group_by(group, session, lag) %>%
  summarise(mean_acf = mean(acf, na.rm = TRUE),
            median_acf = median(acf, na.rm = TRUE),
            .groups = 'drop')

ggplot(acf_mm_1a7_rt, aes(lag, acf, group = participant_id)) +
  geom_line(alpha = 0.25, linewidth = 0.3) +              
  geom_line(data = taskswitch_1a7_rt_summary,                            
            aes(lag, mean_acf),
            inherit.aes = FALSE,
            colour = 'black', linewidth = 1) +
  facet_grid(group ~ session) +
  geom_hline(yintercept = c(-.05, .05), linetype = 'dotted') +  # ±0.05 guides
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Lag', y = 'ACF',
       title = 'Participant ACFs with grand-mean overlay') +
  theme_minimal()

# - 4-way RE*+ variance, autocorr --------------

# Model
# - doesn't converge
taskswitch_1a8_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence
  + cue_transition +
    (task_transition * congruence + cue_transition | participant_id) +
  ou(trial_sequence_factor + 0 | participant_id:session),
  dispformula = ~ poly(fitted_1a2_rt, 2) + session + group + task_transition,
  data = taskswitch_rt_pred, family = gaussian(link = 'log'),
  control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(taskswitch_1a8_rt)

# Diagnostics
diagnose_model(taskswitch_1a8_rt)

# - 4-way RE+ ---------

# Model
# - Doesn't converge
taskswitch_1b_rt <- lme4::glmer(
  response_time ~ group * session * task_transition * congruence + 
    + cue_transition +
    (task_transition + congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_1b_rt)

# Different package
taskswitch_1b1_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence + 
    + cue_transition +
    (task_transition + congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_1b1_rt)

# - 4-way RE+|| ---------

# Model
taskswitch_1c_rt <- lme4::glmer(
  response_time ~ group * session * task_transition * congruence + 
    + cue_transition +
    (task_transition + congruence + cue_transition || participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_1c_rt)

# Different package
taskswitch_1c1_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * task_transition * congruence + 
    + cue_transition +
    (task_transition + congruence + cue_transition || participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_1c1_rt)

# - 3-way RE*+ --------

# Model
# - Doesn't converge
taskswitch_2_rt <- lme4::glmer(
  response_time ~ group * session * (task_transition + congruence) + 
    task_transition * congruence + cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2_rt)

# Different package
taskswitch_2a1_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * (task_transition + congruence) + 
    task_transition * congruence + cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2a1_rt)

# - 3-way RE** --------

# Different package
# - Doesn't converge
taskswitch_2a2_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * (task_transition + congruence) + 
    task_transition * congruence + cue_transition +
    (task_transition * congruence * cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2a2_rt)

# - 3-way RE+ --------

# Model
taskswitch_2b_rt <- lme4::glmer(
  response_time ~ group * session * (task_transition + congruence) + 
    task_transition * congruence + cue_transition +
    (task_transition + congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2b_rt)

# Different package
taskswitch_2b1_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * (task_transition + congruence) + 
    task_transition * congruence + cue_transition +
    (task_transition + congruence + cue_transition | participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2b1_rt)

# - 3-way RE+|| --------

# Model
taskswitch_2c_rt <- lme4::glmer(
  response_time ~ group * session * (task_transition + congruence) + 
    task_transition * congruence + cue_transition +
    (task_transition + congruence + cue_transition || participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2c_rt)

# Different package
taskswitch_2c1_rt <- glmmTMB::glmmTMB(
  response_time ~ group * session * (task_transition + congruence) + 
    task_transition * congruence + cue_transition +
    (task_transition + congruence + cue_transition || participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2c1_rt)

# - 3-way (simpler) RE+|| --------

# Model
taskswitch_2d_rt <- lme4::glmer(
  response_time ~ group * session * (task_transition + congruence) + cue_transition +
    (task_transition + congruence + cue_transition || participant_id),
  data = taskswitch_rt, family = gaussian(link = 'log'))
summary(taskswitch_2d_rt)

# MIXED-EFFECTS FOLLOW-UPS: RT ------------------------------------

# Anova table
car::Anova(taskswitch_1a7_rt, type = 'III')

# - Switch cost RT 4-way RE*+ variance  -----------------------

switch_1a7_rt <- emmeans(taskswitch_1a7_rt, 
                         ~ group * session * task_transition,
                         type = 'response')
switch_1a7_rt

# Tibble for plotting
switch_tibble_1a7_rt <- switch_1a7_rt %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Post-hocs 3-way
switch_3way_1a7_rt <- emmeans(taskswitch_1a7_rt,
                              ~ session * task_transition * group) %>%
  contrast(interaction = c('pairwise', 'pairwise', 'pairwise'))
switch_3way_1a7_rt

# Post-hocs 2-way
switch_2way_1a7_rt <- emmeans(taskswitch_1a7_rt,
                              ~ task_transition * session | group) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
switch_2way_1a7_rt

# Joint test of 2-way contrasts
switch_2way_1a7_rt_joint <- test(switch_2way_1a7_rt, by = NULL, joint = TRUE)
switch_2way_1a7_rt_joint

# Session given trial type
switch_session_1a7_rt <- emmeans(taskswitch_1a7_rt, 
                                 ~ group * session * task_transition) %>%
  contrast('pairwise', by = c('group', 'task_transition'), 
                                  combine = TRUE)
switch_session_1a7_rt

# Trial type given session
switch_trial_1a7_rt <- emmeans(taskswitch_1a7_rt, 
                               ~ group * session * task_transition) %>% 
  contrast('pairwise', by = c('group', 'session'), 
                                combine = TRUE)
switch_trial_1a7_rt

# Baseline comparison
switch_baseline_rt <- emmeans(taskswitch_1a7_rt,
                              ~ group * task_transition | session) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
switch_baseline_rt

# Within-subject SE plot #
# ----------------------------- #

data_rt_switch_1a7_rt <- afex_plot(
  taskswitch_1a7_rt,
  x = 'session',
  trace = 'task_transition',      
  panel = 'group',                
  id = 'participant_id',
  dv = 'response_time',
  data = taskswitch_rt_pred,
  within_vars = c('task_transition', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_rt_switch_1a7_rt_adj <- data_rt_switch_1a7_rt$data %>%
  mutate(
    session = factor(
      session,
      levels = c(1,2),
      labels = c('Pretest' , 'Posttest')),
    task_transition = factor(
      task_transition,
      levels = c('repeat','switch'),
      labels = c('Repeat', 'Switch')
    ))

# Afex means for SE
afex_means_switch_1a7_rt <- data_rt_switch_1a7_rt$means %>%
  as_tibble() %>%
  # Select only what we need to join
  dplyr::select(session, task_transition, group, SE) %>%
  # Ensure factors match the emmeans tibble exactly for the join
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'), 
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Join SEs into the main plotting data
switch_plot_data_1a7_rt <- switch_tibble_1a7_rt %>%
  # Drop SE from the emmeans tibble
  dplyr::select(-SE) %>%
  left_join(afex_means_switch_1a7_rt, by = c('group', 'session', 'task_transition'))

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_switch_rt <- ggplot(switch_plot_data_1a7_rt,
                             aes(x = session, y = response,
                                 colour = group,
                                 linetype = task_transition,
                                 group = interaction(group, task_transition))) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(aes(ymin = response - SE,
                    ymax = response + SE),
                width = .1, position = pos_dodge) +
  geom_jitter(
    data = data_rt_switch_1a7_rt_adj,
    aes(x = session, y = y, colour = group,
        group = interaction(group, task_transition)),
    position = position_jitterdodge(
      jitter.width = 0.05, jitter.height = 0, dodge.width = dodge_tsmm
    ),
    alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  scale_colour_manual(values = pal, guide = 'none') +
  scale_linetype_manual(values = c(Repeat = 'solid', Switch = 'dashed'),
                        name = 'Trial Type') +
  labs(x = NULL, y = 'Reaction Time (ms)') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'mixed_switch_rt.pdf'),
  plot = gg_mixed_switch_rt,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Incongruence cost RT 4-way RE*+ variance  -----------------------

incongr_1a7_rt <- emmeans(taskswitch_1a7_rt, 
                          ~ group * session * congruence,
                          type = 'response')
incongr_1a7_rt

# Tibble for plotting
incongr_tibble_1a7_rt <- incongr_1a7_rt %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'),
                        labels = c('Congruent', 'Incongruent')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Within-subject SE plot
# ----------------------------- #

data_rt_incongr_1a7_rt <- afex_plot(
  taskswitch_1a7_rt,
  x = 'session',
  trace = 'congruence',      
  panel = 'group',                
  id = 'participant_id',
  dv = 'response_time',
  data = taskswitch_rt_pred,
  within_vars = c('congruence', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_rt_incongr_1a7_rt_adj <- data_rt_incongr_1a7_rt$data %>%
  mutate(
    session = factor(
      session,
      levels = c(1,2),
      labels = c('Pretest' , 'Posttest')),
    congruence = factor(
      congruence,
      levels = c('congruent','incongruent'),
      labels = c('Congruent', 'Incongruent')
    ))

# Afex means for SE
afex_means_incongr_1a7_rt <- data_rt_incongr_1a7_rt$means %>%
  as_tibble() %>%
  # Select only what we need to join
  dplyr::select(session, congruence, group, SE) %>%
  # Ensure factors match the emmeans tibble exactly for the join
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'), 
                        labels = c('Congruent', 'Incongruent')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Join SEs into the main plotting data
incongr_plot_data_1a7_rt <- incongr_tibble_1a7_rt %>%
  # Drop SE from the emmeans tibble
  dplyr::select(-SE) %>%
  left_join(afex_means_incongr_1a7_rt, by = c('group', 'session', 'congruence'))

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_incongr_rt <- ggplot(incongr_plot_data_1a7_rt,
                              aes(x = session, y = response,
                                  colour = group,
                                  linetype = congruence,
                                  group = interaction(group, congruence))) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(aes(ymin = response - SE,
                    ymax = response + SE),
                width = .1, position = pos_dodge) +
  geom_jitter(
    data = data_rt_incongr_1a7_rt_adj,
    aes(x = session, y = y, colour = group,
        group = interaction(group, congruence)),
    position = position_jitterdodge(
      jitter.width = 0.05, jitter.height = 0, dodge.width = dodge_tsmm
    ),
    alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  scale_colour_manual(values = pal, guide = 'none') +
  scale_linetype_manual(values = c(Congruent = 'solid', Incongruent = 'dashed'),
                        name = 'Trial Type') +
  labs(x = NULL, y = 'Reaction Time (ms)') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'mixed_incongr_rt.pdf'),
  plot = gg_mixed_incongr_rt,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# TASK SWITCHING ACC ANALYSIS ----------------
# - 4-way RE*+ --------------

# Model
# - doesn't converge
taskswitch_1_er <- glmmTMB::glmmTMB(
  response_correct ~ 
    group * session * task_transition * congruence + 
    cue_transition +
    (task_transition * congruence + cue_transition | participant_id),
  data = taskswitch_er, family = binomial(),
  control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(taskswitch_1_er)

# - 4-way RE*+|| --------------

# Model
taskswitch_2_er <- glmmTMB::glmmTMB(
  response_correct ~ 
    group * session * task_transition * congruence + 
    cue_transition +
    (task_transition * congruence + cue_transition || participant_id),
  data = taskswitch_er, family = binomial(),
  control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(taskswitch_2_er)

# Diagnostics
diagnose_model(taskswitch_2_er)

# Cross-validation #
# ----------------------- #

# Cluster-based
cv_taskswitch_2_er <- cv(taskswitch_2_er, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_2_er)
plot(cv_taskswitch_2_er)

# MIXED-EFFECTS FOLLOW-UPS: ACC ---------

# Anova table
car::Anova(taskswitch_2_er, type = 'III')

# - Switch cost ACC --------

switch_2_er <- emmeans(taskswitch_2_er, 
                       ~ group * session * task_transition,
                       type = 'response')
switch_2_er

# Tibble for plotting
switch_tibble_2_er <- switch_2_er %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Within-subject SE plot #
# ----------------------------- #

# Within-subject SE plot
data_acc_switch_2_er <- afex_plot(
  taskswitch_2_er,
  x = 'session',
  trace = 'task_transition',      
  panel = 'group',                
  id = 'participant_id',
  dv = 'response_correct',
  data = taskswitch_er,
  within_vars = c('task_transition', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_acc_switch_2_er_adj <- data_acc_switch_2_er$data %>%
  mutate(
    session = factor(
      session,
      levels = c(1,2),
      labels = c('Pretest' , 'Posttest')),
    task_transition = factor(
      task_transition,
      levels = c('repeat','switch'),
      labels = c('Repeat', 'Switch')
    ))

# Afex means for SE
afex_means_switch_2_er <- data_acc_switch_2_er$means %>%
  as_tibble() %>%
  # Select only what we need to join
  dplyr::select(session, task_transition, group, SE) %>%
  # Ensure factors match the emmeans tibble exactly for the join
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'), 
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Join SEs into the main plotting data
switch_plot_data_2_er <- switch_tibble_2_er %>%
  # Drop SE from the emmeans tibble
  dplyr::select(-SE) %>%
  left_join(afex_means_switch_2_er, by = c('group', 'session', 'task_transition'))

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_switch_acc <- ggplot(switch_plot_data_2_er,
                              aes(x = session, y = prob * 100,
                                  colour = group,
                                  linetype = task_transition,
                                  group = interaction(group, task_transition))) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(aes(ymin = prob * 100 - SE * 100,
                    ymax = prob * 100 + SE * 100),
                width = .1, position = pos_dodge) +
  geom_jitter(
    data = data_acc_switch_2_er_adj,
    aes(x = session, y = y * 100, colour = group,
        group = interaction(group, task_transition)),
    position = position_jitterdodge(
      jitter.width = 0.05, jitter.height = 0, dodge.width = dodge_tsmm
    ),
    alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  scale_colour_manual(values = pal, guide = 'none') +
  scale_linetype_manual(values = c(Repeat = 'solid', Switch = 'dashed'),
                        name = 'Trial Type') +
  labs(x = NULL, y = 'Mean Accuracy (%)') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'mixed_switch_acc_2_er.pdf'),
  plot = gg_mixed_switch_acc,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Incongruence cost ACC --------

incongruence_2_er <- emmeans(taskswitch_2_er, 
                             ~ group * session * congruence,
                             type = 'response')
incongruence_2_er

# Tibble for plotting
incongruence_tibble_2_er <- incongruence_2_er %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'),
                        labels = c('Congruent', 'Incongruent')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Within-subject SE plot #
# ----------------------------- #

# Within-subject SE plot
data_acc_incongr_2_er <- afex_plot(
  taskswitch_2_er,
  x = 'session',
  trace = 'congruence',      
  panel = 'group',                
  id = 'participant_id',
  dv = 'response_correct',
  data = taskswitch_er,
  within_vars = c('congruence', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_acc_incongr_2_er_adj <- data_acc_incongr_2_er$data %>%
  mutate(
    session = factor(
      session,
      levels = c(1,2),
      labels = c('Pretest' , 'Posttest')),
    congruence = factor(
      congruence,
      levels = c('congruent','incongruent'),
      labels = c('Congruent', 'Incongruent')
    ))

# Afex means for SE
afex_means_incongr_2_er <- data_acc_incongr_2_er$means %>%
  as_tibble() %>%
  # Select only what we need to join
  dplyr::select(session, congruence, group, SE) %>%
  # Ensure factors match the emmeans tibble exactly for the join
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'), 
                        labels = c('Congruent', 'Incongruent')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Join SEs into the main plotting data
incongr_plot_data_2_er <- incongruence_tibble_2_er %>%
  # Drop SE from the emmeans tibble
  dplyr::select(-SE) %>%
  left_join(afex_means_incongr_2_er, by = c('group', 'session', 'congruence'))

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_incongr_acc <- ggplot(incongr_plot_data_2_er,
                               aes(x = session, y = prob * 100,
                                   colour = group,
                                   linetype = congruence,
                                   group = interaction(group, congruence))) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(aes(ymin = prob * 100 - SE * 100,
                    ymax = prob * 100 + SE * 100),
                width = .1, position = pos_dodge) +
  geom_jitter(
    data = data_acc_incongr_2_er_adj,
    aes(x = session, y = y * 100, colour = group,
        group = interaction(group, congruence)),
    position = position_jitterdodge(
      jitter.width = 0.05, jitter.height = 0, dodge.width = dodge_tsmm
    ),
    alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  scale_colour_manual(values = pal, guide = 'none') +
  scale_linetype_manual(values = c(Congruent = 'solid', Incongruent = 'dashed'),
                        name = 'Trial Type') +
  labs(x = NULL, y = 'Mean Accuracy (%)') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'mixed_incongr_acc_2_er.pdf'),
  plot = gg_mixed_incongr_acc,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# TASK SWITCHING RT & ACC BASELINE --------------------------

# RT #
# ========= #

# Grid
taskswitch_emm_base_rt <- emmeans(taskswitch_1a7_rt, ~ group | session)
taskswitch_emm_base_rt

# Comparison
taskswitch_base_rt <- contrast(taskswitch_emm_base_rt, 'pairwise')
taskswitch_base_rt

# ACC #
# ========= #

# Grid
taskswitch_emm_base_er <- emmeans(taskswitch_2_er, ~ group | session)
taskswitch_emm_base_er

# Comparison
taskswitch_base_er <- contrast(taskswitch_emm_base_er, 'pairwise')
taskswitch_base_er

# ASRS ANALYSIS ---------------------

glimpse(questionnaires_long)

# - Visualise data -------------------

# RAW VALUES #
# ===================== #

# Density
ggplot(questionnaires_long, aes(x = ASRS, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'ASRS',
       y = 'Density',
       fill = 'Group',
       title = 'ASRS Distributions by Group and Session') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_ASRS_with <- ggplot(questionnaires_long,
       aes(x = session, y = ASRS, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                                    'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = questionnaires_trwin %>% 
                  filter(dv == 'ASRS'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(x = NULL, y = 'ASRS Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = c(0, 24)) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_with.pdf'),
  plot = gg_ASRS_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_ASRS_without <- ggplot(questionnaires_long,
                        aes(x = session,
                            y = ASRS,
                            fill = group,
                            colour = group)) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA,
    width = 0.5,
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      jitter.height = 0,
      dodge.width = dodge_w,
      seed = 1
    ),
    alpha = 0.5,
    size = 1.8,
    shape = 21,
    stroke = 0.4
  ) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w),
               show.legend = FALSE) + 
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'ASRS Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_without.pdf'),
  plot = gg_ASRS_without,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# ROBUSTIFICATION #
# ===================== #

# Plot
pd <- position_dodge(width = 0.6) 
gg_ASRS_trbar <- ggplot(questionnaires_trwin %>% filter(dv == 'ASRS'),
       aes(x = session,
           y = mean_tr,
           fill = group,
           colour = group,
           group = group)) +
  geom_col(position = pd, width = 0.6, colour = NA) + 
  geom_errorbar(aes(ymin = mean_tr - se_tr_ws,
                    ymax = mean_tr + se_tr_ws),
                position = pd, width = 0.15, linewidth = 0.6) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Session',
       y = 'Trimmed Mean ASRS Score',
       fill = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_trbar.pdf'),
  plot = gg_ASRS_trbar,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# CHANGE SCORES #
# ===================== #

ASRS_wide <- questionnaires_ready %>%
  dplyr::select(participant_id, group, cohort, age, 
                ASRS_total_pre, ASRS_total_post) %>%
  rename(ASRS_pre = ASRS_total_pre,
         ASRS_post = ASRS_total_post) %>%
  mutate(ASRS_change = ASRS_pre - ASRS_post)
glimpse(ASRS_wide)

# Density
ggplot(ASRS_wide, aes(x = ASRS_change, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'ASRS Change Scores',
       y = 'Density',
       fill = 'Group',
       title = 'ASRS Change Scores Distributions by Group') +
  theme_minimal()

# Boxplots
gg_ASRS_change <- ggplot(ASRS_wide, aes(x = group,
                      y = ASRS_change,
                      fill = group,
                      colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(
    aes(fill = group, colour = group),
    width = 0.05,         
    height = 0,          
    alpha = 0.6,
    size = 2      
  ) +
  stat_summary(fun = function(z) mean(z, trim = 0.2, na.rm = TRUE), geom = 'point',
            shape = 18, size = 3.5, colour = 'black', 
            show.legend = FALSE) +
  stat_summary(fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05) + 
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, guide = 'none') + 
  labs(x = 'Group',
       y = 'ASRS Change Score',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_change.pdf'),
  plot = gg_ASRS_change,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Trimmed mixed ANOVA -----------------------

ASRS_tranova <- WRS2::bwtrim(ASRS ~ group * session, 
                             id = participant_id, 
                             data = questionnaires_long,
                             tr = 0.2
)
ASRS_tranova

# Baseline #
# ------------------- #

# Yuen's t-test
ASRS_baseline_test <- yuen(ASRS ~ group, 
                           data = questionnaires_long %>% filter(session == 1),
                           tr = 0.2)
ASRS_baseline_test

# Bootstrapped Yuen
ASRS_baseline_test_bt <- yuenbt(ASRS ~ group, 
                                data = questionnaires_long %>% filter(session == 1),
                                tr = 0.2,
                                nboot = 100000,
                                side = TRUE)
ASRS_baseline_test_bt

# Quantile comparison
ASRS_quant_baseline <- WRS2::qcomhd(ASRS ~ group, 
                                    data = questionnaires_long %>% 
                                      filter(session == 1), 
                                    q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                                    nboot = 50000, alpha = 0.05, ADJ.CI = TRUE)
ASRS_quant_baseline

# Post-hocs #
# -------------------- #

# Create ordered list for analysis
ASRS_post_tr_list <- build_split_list(questionnaires_long, 
                                      'ASRS', c('group', 'session'))

# Contrasts on trimmed means
# - Increasing nboot here can lead to NAs 
# - Bootstrap-t method
ASRS_post_tr <- WRS::bwmcp(2, 2, ASRS_post_tr_list, tr=0.2,con=0, nboot=1000)
ASRS_post_tr

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
# - No bootstrap
ASRS_post_tr_int <- WRS::bwimcpES(2, 2, ASRS_post_tr_list, tr=0.2, 
                                  CI=TRUE, alpha=0.05)
ASRS_post_tr_int

# Interaction contrast on trimmed difference scores with bootstrap
# - Only for comparison with bwimcpES
# - Compute difference scores and then trim
# - Percentile bootstrap method
ASRS_post_tr_int_boot <- WRS::spmcpi(2, 2, ASRS_post_tr_list,
                                     est=tmean, alpha= 0.05, nboot=50000, 
                                     SEED=TRUE,SR=FALSE,pr=TRUE)
ASRS_post_tr_int_boot

# NOTE:
# - AKP: robust version of Cohen's D
# - EP: explanatory measure of effect size
# - QS (median): quantile shift of the median
# - QStr: quantile shift of the trimmed mean
# - WMW: probability p = P(Xi1 < Xi2)
# - KMS: Kulinskaya and Staudte measure (not robust though)

# Session changes KD
ASRS_session_KD <- yuend(questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS),
                         questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS),
                         tr = 0.2)
ASRS_session_KD

# Session change CD
ASRS_session_CD <- yuend(questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS),
                         questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS),
                         tr = 0.2)
ASRS_session_CD

# Effect size of session changes
ASRS_session_eff <- bw.es.B(2, 2, ASRS_post_tr_list, 
                            tr = 0.2, POOL = FALSE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
ASRS_session_eff

# Pooled effect size (session)
ASRS_session_eff_pool <- bw.es.B(2, 2, ASRS_post_tr_list, 
                                 tr = 0.2, POOL = TRUE, OPT = FALSE, 
                                 CI = TRUE, SEED = TRUE, REL.MAG = NULL)
ASRS_session_eff_pool

# - Prepare data for balancing -------------------

# NOTE: the following tibble contains only the questionnaire data 
# and therefore all the participants who completed the questionnaires

explore_tibble_quest <- questionnaires_long %>%
  pivot_wider(
    id_cols = c(participant_id, group, cohort, age),
    names_from = session,
    values_from = -c(participant_id, group, cohort, age, session),
    names_glue = '{.value}_{session}'
  )
glimpse(explore_tibble_quest)

# Centre covariates
explore_tibble_quest_centred <- explore_tibble_quest %>%
  mutate(across(ends_with('_1'), ~ as.numeric(scale(.x, scale = FALSE))))
glimpse(explore_tibble_quest_centred)

# - ASRS robust ANCOVA ---------------------

# NOTE: After cardinality 'matching', the groups remain independent

# Balance groups
balanced_quest <- matchit(
  group ~ ASRS_1,
  data = explore_tibble_quest,
  method = 'cardinality',
  tols = .05,
  std.tols = TRUE)
summary(balanced_quest)

# Quick balancing assessment
bal.plot(balanced_quest, 
         var.name = 'ASRS_1', 
         which = 'both',
         type = 'density',
         colors = pal)

# Extract the matched data
balanced_ASRS <- match_data(balanced_quest) %>%
  dplyr::select(participant_id, group, ASRS_1, ASRS_2)
glimpse(balanced_ASRS)

# TWO-STEP ROBUST ANCOVA #
# ======================= #

# Robust regression
ASRS_reg <- MASS::rlm(ASRS_2 ~ ASRS_1, 
                      data = balanced_ASRS,
                      method = 'MM',
                      psi = MASS::psi.bisquare)
summary(ASRS_reg)

# Extract and embed weighted residuals
balanced_ASRS <- balanced_ASRS %>%
  mutate(ASRS_2_res = ASRS_reg$residuals,
         ASRS_2_res_log = log(ASRS_2_res + 10))

# Plot the residuals
plot(ASRS_reg$fitted.values, balanced_ASRS$ASRS_2_res)

# Plot baseline by group #
# ------------------------------ #

# Boxplots
gg_balanced_ASRS_baseline <- ggplot(balanced_ASRS,
                                    aes(x = group,
                                        y = ASRS_1,
                                        fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(
    aes(fill = group, colour = group),
    width = 0.2,
    height = 0,
    alpha = 0.6
  ) +
  stat_summary(fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w),
               show.legend = FALSE) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') + 
  labs(x = 'Group',
       y = 'ASRS Pretest',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'balanced_ASRS_baseline.pdf'),
  plot = gg_balanced_ASRS_baseline,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Plot the residuals by group #
# ---------------------------------- #

# Density
ggplot(balanced_ASRS, 
       aes(x = ASRS_2_res, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'ASRS Residuals',
       y = 'Density',
       fill = 'Group',
       title = 'ASRS Residuals Distributions by Group') +
  theme_minimal()

# Boxplots
gg_res_ASRS_posttest <- ggplot(balanced_ASRS,
                               aes(x = group,
                                   y = ASRS_2_res,
                                   fill = group,
                                   colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(
    aes(fill = group, colour = group),
    width = 0.2,
    height = 0,
    alpha = 0.6
  ) +
  stat_summary(fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w)) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') + 
  labs(x = 'Group',
       y = 'ASRS Posttest Residuals',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'res_ASRS_posttest.pdf'),
  plot = gg_res_ASRS_posttest,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Check baseline balance #
# -------------------------- #

cor.test(balanced_ASRS$ASRS_1, as.numeric(balanced_ASRS$group),
         exact = FALSE, method = 'pearson')

# Test the residuals #
# -------------------------- #

# Yuen's bootstrapped t-test
ASRS_respost_yuenbt <- yuenbt(ASRS_2_res ~ group,
                              data = balanced_ASRS,
                              tr = 0.2,
                              nboot = 100000,
                              side = TRUE)
ASRS_respost_yuenbt

# Effect size 
ASRS_respost_yuenbt_eff <- yuen.effect.ci(ASRS_2_res ~ group,
                                          data = balanced_ASRS,
                                          tr = 0.2,
                                          nboot = 100000)
ASRS_respost_yuenbt_eff

# ASRS INATTENTION ---------------
# - Visualise data -------------------

# RAW VALUES #
# ===================== #

# Density
ggplot(questionnaires_long, aes(x = ASRS_IA, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'ASRS INATTENTION SCALE',
       y = 'Density',
       fill = 'Group',
       title = 'ASRS INATTENTION SCALE Distributions by Group and Session') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_ASRS_IA_with <- ggplot(questionnaires_long,
                       aes(x = session, y = ASRS_IA, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = questionnaires_trwin %>% 
                  filter(dv == 'ASRS_IA'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(x = NULL, y = 'ASRS INATTENTION Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = c(0, 24)) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_IA_with.pdf'),
  plot = gg_ASRS_IA_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_ASRS_IA_without <- ggplot(questionnaires_long,
                          aes(x = session,
                              y = ASRS_IA,
                              fill = group,
                              colour = group)) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA,
    width = 0.5,
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      jitter.height = 0,
      dodge.width = dodge_w,
      seed = 1
    ),
    alpha = 0.5,
    size = 1.8,
    shape = 21,
    stroke = 0.4
  ) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w),
               show.legend = FALSE) + 
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'ASRS INATTENTION Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_IA_without.pdf'),
  plot = gg_ASRS_IA_without,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)



# - Trimmed mixed ANOVA -----------------------

ASRS_IA_tranova <- WRS2::bwtrim(ASRS_IA ~ group * session, 
                             id = participant_id, 
                             data = questionnaires_long,
                             tr = 0.2
                             )
ASRS_IA_tranova

# Baseline #
# ------------------- #

# Yuen's t-test
ASRS_IA_baseline_test <- yuen(ASRS_IA ~ group, 
                           data = questionnaires_long %>% filter(session == 1),
                           tr = 0.2)
ASRS_IA_baseline_test

# Posttest #
# ------------------- #

# Yuen's t-test
ASRS_IA_post_test <- yuen(ASRS_IA ~ group, 
                              data = questionnaires_long %>% filter(session == 2),
                              tr = 0.2)
ASRS_IA_post_test

# Post-hocs #
# -------------------- #

# Create ordered list for analysis
ASRS_IA_post_tr_list <- build_split_list(questionnaires_long, 
                                      'ASRS_IA', c('group', 'session'))

# Contrasts on trimmed means
# - Increasing nboot here can lead to NAs 
# - Bootstrap-t method
ASRS_IA_post_tr <- WRS::bwmcp(2, 2, ASRS_IA_post_tr_list, tr=0.2,con=0, nboot=1000)
ASRS_IA_post_tr

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
# - No bootstrap
ASRS_IA_post_tr_int <- WRS::bwimcpES(2, 2, ASRS_IA_post_tr_list, tr=0.2, 
                                  CI=TRUE, alpha=0.05)
ASRS_IA_post_tr_int

# NOTE:
# - AKP: robust version of Cohen's D
# - EP: explanatory measure of effect size
# - QS (median): quantile shift of the median
# - QStr: quantile shift of the trimmed mean
# - WMW: probability p = P(Xi1 < Xi2)
# - KMS: Kulinskaya and Staudte measure (not robust though)

# Session changes KD
ASRS_IA_session_KD <- yuend(questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS_IA),
                         questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS_IA),
                         tr = 0.2)
ASRS_IA_session_KD

# Session change CD
ASRS_IA_session_CD <- yuend(questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS_IA),
                         questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           dplyr::pull(ASRS_IA),
                         tr = 0.2)
ASRS_IA_session_CD

# Effect size of session changes
ASRS_IA_session_eff <- bw.es.B(2, 2, ASRS_IA_post_tr_list, 
                           tr = 0.2, POOL = FALSE, OPT = FALSE, 
                           CI = TRUE, SEED = TRUE, REL.MAG = NULL)
ASRS_IA_session_eff

# Pooled effect size (session)
ASRS_IA_session_eff_pool <- bw.es.B(2, 2, ASRS_IA_post_tr_list, 
                            tr = 0.2, POOL = TRUE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
ASRS_IA_session_eff_pool

# - ASRS robust ANCOVA ---------------------

# NOTE: After cardinality 'matching', the groups remain independent

# Balance groups
balanced_quest_IA <- matchit(
  group ~ ASRS_IA_1,
  data = explore_tibble_quest,
  method = 'cardinality',
  tols = .05,
  std.tols = TRUE)
summary(balanced_quest_IA)

# Quick balancing assessment
bal.plot(balanced_quest_IA, 
         var.name = 'ASRS_IA_1', 
         which = 'both',
         type = 'density',
         colors = pal)

# Extract the matched data
balanced_ASRS_IA <- match_data(balanced_quest_IA) %>%
  dplyr::select(participant_id, group, ASRS_IA_1, ASRS_IA_2)
glimpse(balanced_ASRS_IA)

# TWO-STEP ROBUST ANCOVA #
# ======================= #

# Robust regression
ASRS_IA_reg <- MASS::rlm(ASRS_IA_2 ~ ASRS_IA_1, 
                      data = balanced_ASRS_IA,
                      method = 'MM',
                      psi = MASS::psi.bisquare)
summary(ASRS_IA_reg)

# Extract and embed weighted residuals
balanced_ASRS_IA <- balanced_ASRS_IA %>%
  mutate(ASRS_IA_2_res = ASRS_IA_reg$residuals)

# Plot the residuals
plot(ASRS_IA_reg$fitted.values, balanced_ASRS_IA$ASRS_IA_2_res)

# Plot baseline by group #
# ------------------------------ #

# Boxplots
gg_balanced_ASRS_IA_baseline <- ggplot(balanced_ASRS_IA,
                                    aes(x = group,
                                        y = ASRS_IA_1,
                                        fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(
    aes(fill = group, colour = group),
    width = 0.2,
    height = 0,
    alpha = 0.6
  ) +
  stat_summary(fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w),
               show.legend = FALSE) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') + 
  labs(x = 'Group',
       y = 'ASRS Inattention Pretest',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'balanced_ASRS_IA_baseline.pdf'),
  plot = gg_balanced_ASRS_IA_baseline,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Plot the residuals by group #
# ---------------------------------- #

# Density
ggplot(balanced_ASRS_IA, 
       aes(x = ASRS_IA_2_res, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'ASRS Inattention Residuals',
       y = 'Density',
       fill = 'Group',
       title = 'ASRS Inattention Residuals Distributions by Group') +
  theme_minimal()

# Boxplots
gg_res_ASRS_IA_posttest <- ggplot(balanced_ASRS_IA,
                               aes(x = group,
                                   y = ASRS_IA_2_res,
                                   fill = group,
                                   colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(
    aes(fill = group, colour = group),
    width = 0.2,
    height = 0,
    alpha = 0.6
  ) +
  stat_summary(fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w)) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') + 
  labs(x = 'Group',
       y = 'ASRS Inattention Posttest Residuals',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'res_ASRS_IA_posttest.pdf'),
  plot = gg_res_ASRS_IA_posttest,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Check baseline balance #
# -------------------------- #

cor.test(balanced_ASRS_IA$ASRS_IA_1, as.numeric(balanced_ASRS_IA$group),
         exact = FALSE, method = 'pearson')

# Test the residuals #
# -------------------------- #

# Yuen's bootstrapped t-test
ASRS_IA_respost_yuenbt <- yuenbt(ASRS_IA_2_res ~ group,
                              data = balanced_ASRS_IA,
                              tr = 0.2,
                              nboot = 100000,
                              side = TRUE)
ASRS_IA_respost_yuenbt

# Effect size 
ASRS_IA_respost_yuenbt_eff <- yuen.effect.ci(ASRS_IA_2_res ~ group,
                                          data = balanced_ASRS_IA,
                                          tr = 0.2,
                                          nboot = 100000)
ASRS_IA_respost_yuenbt_eff

# ASRS HYPERACTIVITY ---------------
# - Visualise data -------------------

# RAW VALUES #
# ===================== #

# Density
ggplot(questionnaires_long, aes(x = ASRS_H, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'ASRS HYPERACTIVITY SCALE',
       y = 'Density',
       fill = 'Group',
       title = 'ASRS HYPERACTIVITY SCALE Distributions by Group and Session') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_ASRS_H_with <- ggplot(questionnaires_long,
                          aes(x = session, y = ASRS_H, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = questionnaires_trwin %>% 
                  filter(dv == 'ASRS_H'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(x = NULL, y = 'ASRS HYPERACTIVITY Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = c(0, 24)) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_H_with.pdf'),
  plot = gg_ASRS_H_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_ASRS_H_without <- ggplot(questionnaires_long,
                             aes(x = session,
                                 y = ASRS_H,
                                 fill = group,
                                 colour = group)) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA,
    width = 0.5,
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      jitter.height = 0,
      dodge.width = dodge_w,
      seed = 1
    ),
    alpha = 0.5,
    size = 1.8,
    shape = 21,
    stroke = 0.4
  ) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w),
               show.legend = FALSE) + 
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'ASRS HYPERACTIVITY Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ASRS_H_without.pdf'),
  plot = gg_ASRS_H_without,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Trimmed mixed ANOVA -----------------------

ASRS_H_tranova <- WRS2::bwtrim(ASRS_H ~ group * session, 
                                id = participant_id, 
                                data = questionnaires_long,
                                tr = 0.2
)
ASRS_H_tranova


# Post-hocs #
# -------------------- #

# Create ordered list for analysis
ASRS_H_post_tr_list <- build_split_list(questionnaires_long, 
                                         'ASRS_H', c('group', 'session'))

# Pooled effect size (session)
ASRS_H_session_eff_pool <- bw.es.B(2, 2, ASRS_H_post_tr_list, 
                                    tr = 0.2, POOL = TRUE, OPT = FALSE, 
                                    CI = TRUE, SEED = TRUE, REL.MAG = NULL)
ASRS_H_session_eff_pool

# PSQI ANALYSIS -----------------
# - Visualise data -------------------

# RAW VALUES #
# ===================== #

# Density
ggplot(questionnaires_long, aes(x = PSQI, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'PSQI',
       y = 'Density',
       fill = 'Group',
       title = 'PSQI Distributions by Group and Session') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_PSQI_with <- ggplot(questionnaires_long,
       aes(x = session, y = PSQI, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = questionnaires_trwin %>% 
                  filter(dv == 'PSQI'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = NULL, y = 'PSQI Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'PSQI_with.pdf'),
  plot = gg_PSQI_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Boxplots without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_PSQI_without <- ggplot(questionnaires_long,
       aes(x = session,
           y = PSQI,
           fill = group,
           colour = group)) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA,
    width = 0.5,
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      jitter.height = 0,
      dodge.width = dodge_w,
      seed = 1
    ),
    alpha = 0.5,
    size = 1.8,
    shape = 21,
    stroke = 0.4
  ) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w),
               show.legend = FALSE) + 
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'PSQI Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'PSQI_without.pdf'),
  plot = gg_PSQI_without,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# ROBUSTIFICATION #
# ===================== #

# Plot
pd <- position_dodge(width = 0.6)
gg_PSQI_trbar <- ggplot(questionnaires_trwin %>% filter(dv == 'PSQI'),
       aes(x = session,
           y = mean_tr,
           fill = group,
           colour = group,
           group = group)) +
  geom_col(position = pd, width = 0.6, colour = NA) + 
  geom_errorbar(aes(ymin = mean_tr - se_tr,
                    ymax = mean_tr + se_tr),
                position = pd, width = 0.15, linewidth = 0.6) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Session',
       y = 'Trimmed Mean PSQI Score',
       fill = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'PSQI_trbar.pdf'),
  plot = gg_PSQI_trbar,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# CHANGE SCORES #
# ===================== #

PSQI_wide <- questionnaires_ready %>%
  dplyr::select(participant_id, group, cohort, age, PSQI_total_pre, PSQI_total_post) %>%
  rename(PSQI_pre = PSQI_total_pre,
         PSQI_post = PSQI_total_post) %>%
  mutate(PSQI_change = PSQI_pre - PSQI_post)
glimpse(PSQI_wide)

# Density
ggplot(PSQI_wide, aes(x = PSQI_change, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'PSQI Change Scores',
       y = 'Density',
       fill = 'Group',
       title = 'PSQI Change Scores Distributions by Group') +
  theme_minimal()

# Boxplots
gg_PSQI_change <- ggplot(PSQI_wide, aes(x = group,
                      y = PSQI_change,
                      fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(
    aes(fill = group, colour = group),
    width = 0.05,         
    height = 0,          
    alpha = 0.6,
    size = 2      
  ) +
  stat_summary(fun = function(z) mean(z, trim = 0.2, na.rm = TRUE), geom = 'point',
               shape = 18, size = 3.5, colour = 'black', 
               show.legend = FALSE) +
  stat_summary(fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') + 
  labs(x = 'Group',
       y = 'PSQI Change Score',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() +
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'PSQI_change.pdf'),
  plot = gg_PSQI_change,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Trimmed mixed ANOVA -----------------------

PSQI_tranova <- WRS2::bwtrim(PSQI ~ group * session, 
                             id = participant_id, 
                             data = questionnaires_long,
                             tr = 0.2
)
PSQI_tranova

# Baseline #
# ------------------- #

PSQI_baseline_test <- yuen(PSQI ~ group, 
                           data = questionnaires_long %>% filter(session == 1),
                           tr = 0.2)
PSQI_baseline_test

# Post-hocs #
# -------------------- #

# Create ordered list for analysis
PSQI_post_tr_list <- build_split_list(questionnaires_long, 
                                      'PSQI', c('group', 'session'))

# Contrasts on trimmed means
PSQI_post_tr <- bwmcp(2, 2, PSQI_post_tr_list, tr=0.2,con=0, nboot=1000)
PSQI_post_tr

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
PSQI_post_tr_int <- bwimcpES(2, 2, PSQI_post_tr_list, tr=0.2, CI=TRUE, alpha=0.05)
PSQI_post_tr_int

# NOTE:
# - AKP: robust version of Cohen's D
# - EP: explanatory measure of effect size
# - QS (median): quantile shift of the median
# - QStr: quantile shift of the trimmed mean
# - WMW: probability p = P(Xi1 < Xi2) (not robust when heteroscedasticity)
# - KMS: Kulinskaya and Staudte measure (not robust though)

# Main effect of Session #
# ---------------------------- #

# Effect size of the main effect of session
PSQI_session_eff <- bw.es.B(2, 2, PSQI_post_tr_list, 
                            tr = 0.2, POOL = TRUE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
PSQI_session_eff

# Main effect of Group #
# ----------------------------- #

# Effect size of the main effect of group (only separately for 
# each time point)
PSQI_group_eff <- bw.es.A(2, 2, PSQI_post_tr_list, tr = 0.2, 
                          pr = TRUE, fun = ES.summary.CI)
PSQI_group_eff

# Pooled
PSQI_group_eff_pool <- yuen.effect.ci(PSQI ~ group,
                                      data = questionnaires_long,
                                      tr = 0.2,
                                      nboot = 1000)
PSQI_group_eff_pool

# BDI ANALYSIS -------------------------------------------------
# - Visualise data ------------------

# RAW VALUES #
# ===================== #

# Density
ggplot(questionnaires_long, aes(x = BDI, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'BDI',
       y = 'Density',
       fill = 'Group',
       title = 'BDI Distributions by Group and Session') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_BDI_with <- ggplot(questionnaires_long,
       aes(x = session, y = BDI, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, 
             labeller = labeller(group = c('CD' = 'Clean Diet', 
                                           'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = questionnaires_trwin %>% 
                  filter(dv == 'BDI'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = NULL, y = 'BDI Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa() + 
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_with.pdf'),
  plot = gg_BDI_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Boxplots without trajectories
gg_BDI_without <- ggplot(questionnaires_long,
       aes(x = session,
           y = BDI,
           fill = group,
           colour = group)) +
  geom_boxplot(
    alpha = 0.6,
    outlier.shape = NA,
    width = 0.5,
    position = position_dodge(width = dodge_w)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      jitter.height = 0,
      dodge.width = dodge_w,
      seed = 1
    ),
    alpha = 0.5,
    size = 1.8,
    shape = 21, 
    stroke = 0.4
  ) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point',
               shape = 18,
               size = 3.5,
               colour = 'black',
               position = position_dodge(width = dodge_w),
               show.legend = FALSE) + 
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'BDI Score', fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_without.pdf'),
  plot = gg_BDI_without,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# ROBUSTIFICATION #
# ===================== #

# Plot
pd <- position_dodge(width = 0.6) 
gg_BDI_trbar <- ggplot(questionnaires_trwin %>% filter(dv == 'BDI'),
       aes(x = session,
           y = mean_tr,
           fill = group,
           colour = group,
           group = group)) +
  geom_col(position = pd, width = 0.6, colour = NA) + 
  geom_errorbar(aes(ymin = mean_tr - se_tr,
                    ymax = mean_tr + se_tr),
                position = pd, width = 0.15, linewidth = 0.6) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Session',
       y = 'Trimmed Mean BDI Score',
       fill = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_trbar.pdf'),
  plot = gg_BDI_trbar,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# CHANGE SCORES #
# ===================== #

BDI_wide <- questionnaires_ready %>%
  dplyr::select(participant_id, group, cohort, age, BDI_total_pre, BDI_total_post) %>%
  rename(BDI_pre = BDI_total_pre,
         BDI_post = BDI_total_post) %>%
  mutate(BDI_change = BDI_pre - BDI_post)
glimpse(BDI_wide)

# Density
ggplot(BDI_wide, aes(x = BDI_change, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'BDI Change Scores',
       y = 'Density',
       fill = 'Group',
       title = 'BDI Change Scores Distributions by Group') +
  theme_minimal()

# Boxplots
gg_BDI_change <- ggplot(BDI_wide, aes(x = group,
                     y = BDI_change,
                     fill = group,
                     colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(
    aes(fill = group, colour = group),
    width = 0.05,         
    height = 0,          
    alpha = 0.6,
    size = 2      
  ) +
  stat_summary(fun = function(z) mean(z, trim = 0.2, na.rm = TRUE), geom = 'point',
               shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  stat_summary(fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') + 
  labs(x = 'Group',
       y = 'BDI Change Score',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() +
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_change.pdf'),
  plot = gg_BDI_change,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Trimmed mixed ANOVA -----------------------

BDI_tranova <- WRS2::bwtrim(BDI ~ group * session, 
                             id = participant_id, 
                             data = questionnaires_long,
                             tr = 0.2
)
BDI_tranova

# Baseline #
# ------------------- #

# Yuen's t-test
BDI_baseline_test <- yuen(BDI ~ group, 
                           data = questionnaires_long %>% filter(session == 1),
                           tr = 0.2)
BDI_baseline_test

# Comparing quantiles
BDI_quant_baseline <- WRS2::qcomhd(BDI ~ group, 
                          data = questionnaires_long %>% 
                            filter(session == 1), 
                          q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                          nboot = 2000, alpha = 0.05, ADJ.CI = TRUE)
BDI_quant_baseline

# Compare variances
BDI_perm <- permg(questionnaires_long %>% filter(session == 1,
                                                  group == 'CD') %>% dplyr::pull(BDI),
                   questionnaires_long %>% filter(session == 1, 
                                                  group == 'KD') %>% dplyr::pull(BDI),
                   alpha = 0.05, est = var, nboot = 10000)
BDI_perm

# Kolmogorov-Smirnov 
BDI_ks_baseline <- ks.test(questionnaires_long %>% 
                             filter(session == 1,
                                    group == 'CD') %>% 
                             dplyr::pull(BDI),
                           questionnaires_long %>% 
                             filter(session == 1,
                                    group == 'KD') 
                           %>% dplyr::pull(BDI),
                           alternative = 'two.sided',
                           simulate.p.value = TRUE,
                           B = 10000)
BDI_ks_baseline

# Post-hocs #
# -------------------- #

# Create ordered list for analysis
BDI_post_tr_list <- build_split_list(questionnaires_long, 
                                      'BDI', c('group', 'session'))

# Contrasts on trimmed means
BDI_post_tr <- bwmcp(2, 2, BDI_post_tr_list, tr=0.2,con=0, nboot=1000)
BDI_post_tr

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
BDI_post_tr_int <- bwimcpES(2, 2, BDI_post_tr_list, tr=0.2, CI=TRUE, alpha=0.05)
BDI_post_tr_int

# NOTE:
# - AKP: robust version of Cohen's D
# - EP: explanatory measure of effect size
# - QS (median): quantile shift of the median
# - QStr: quantile shift of the trimmed mean
# - WMW: probability p = P(Xi1 < Xi2)
# - KMS: Kulinskaya and Staudte measure (not robust though)

# Main effect of Session #
# ---------------------------- #

# Effect size of the main effect of session
BDI_session_eff <- bw.es.B(2, 2, BDI_post_tr_list, 
                            tr = 0.2, POOL = TRUE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
BDI_session_eff

# MEDIATION --------------------------

# MEDIATION: Do ketone changes and / or BMI changes mediate the greater 
# improvement in ASRS in the KD group? As ketones varied only in the KD
# group, only this group was used for these analyses. 

# - Prepare data -------------

# Weighing ketones #
# ------------------------ #

# Define weights for ALL 8 time points (9–16)
ketone_weights <- sqrt(c(1:7, 7))  # Repeat weight for time 16 = time 15

# Normalize to sum to 1
ketone_weights_norm <- ketone_weights / sum(ketone_weights)

# Create weight lookup table
ketone_weights_norm_tbl <- tibble(
  measurement = 9:16,
  weights_norm = ketone_weights_norm
)

# Compute weighted average
ketones_weighed <- ketones_tibble %>%
  dplyr::select(participant_id, group, cohort, ketones_pre, 
                matches('^ketones_post_int_(9|1[0-6])$')) %>%
  pivot_longer(cols = starts_with('ketones_post_int'),
               names_to = 'measurement',
               values_to = 'ketone_value',
               names_pattern = '_(\\d+)$',
               names_transform = list(measurement = as.integer)) %>%
  left_join(ketone_weights_norm_tbl, by = 'measurement') %>%
  mutate(ketones_post_norm = ketone_value * weights_norm) %>%
  group_by(participant_id, group, cohort) %>%
  summarise(
    ketones_post = sum(ketones_post_norm, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  left_join(ketones_tibble %>% 
              dplyr::select(-starts_with('ketones_post_int_')), 
            by = c('participant_id', 'group', 'cohort'))
glimpse(ketones_weighed)

# Implement weighing #
# ------------------------ #

# Merge questionnaires with the ketones
qt_ketones <- questionnaires_long %>%
  pivot_wider(
    id_cols = c(participant_id, group, cohort, age),
    names_from = session,
    values_from = -c(participant_id, group, cohort, age, session),
    names_glue = '{.value}_{session}'
  ) %>%
  left_join(ketones_weighed, by = c('participant_id', 'group', 'cohort')) %>%
  # Ketone mechanism only relevant for KD (stable for CD)
  filter(group == 'KD') %>%
  mutate(
    ketones_change = ketones_post - ketones_pre,
    BMI_change = BMI_1 - BMI_2,
    ASRS_change = ASRS_1 - ASRS_2,
    BDI_change = BDI_1 - BDI_2,
    PSQI_change = PSQI_1 - PSQI_2,
    across(c(ends_with('_change'), ends_with('_1'), 'ketones_pre'), 
           ~ as.numeric(scale(.x, scale = FALSE))))
glimpse(qt_ketones)

qt_ketones_both <- questionnaires_long %>%
  pivot_wider(
    id_cols = c(participant_id, group, cohort, age),
    names_from = session,
    values_from = -c(participant_id, group, cohort, age, session),
    names_glue = '{.value}_{session}'
  ) %>%
  left_join(ketones_weighed, by = c('participant_id', 'group', 'cohort')) %>%
  mutate(
    ketones_change = ketones_post - ketones_pre,
    BMI_change = BMI_1 - BMI_2,
    ASRS_change = ASRS_1 - ASRS_2,
    BDI_change = BDI_1 - BDI_2,
    PSQI_change = PSQI_1 - PSQI_2,
    across(c(ends_with('_change'), ends_with('_1'), 'ketones_pre'), 
           ~ as.numeric(scale(.x, scale = FALSE))))
glimpse(qt_ketones_both)

# - Visualise data ----------------

# Visualise weighted ketones #
# ================================ #

weighted_ketones_plot <- ggplot(ketones_weighed %>% filter(group == 'KD'),
                                aes(x = group, y = ketones_post, colour = group)) +
  geom_jitter(size = 2.6, width = 0.2, alpha = 0.7) +
  scale_colour_manual(values = pal, guide = 'none') +
  geom_hline(yintercept = 0.5, color = "black", 
             linetype = "dashed", linewidth = 0.5) +
  labs(x = NULL, y = 'WEIGHTED Intervention Ketones (mmol/L)') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'weighted_ketones.pdf'),
  plot = weighted_ketones_plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Visualise ketones change #
# ================================ #

ketones_change_plot <- ggplot(qt_ketones %>% filter(group == 'KD'),
                              aes(x = group, y = ketones_change, colour = group)) +
  geom_jitter(size = 2.6, width = 0.2, alpha = 0.7) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = NULL, y = 'Ketones Change (mmol/L) - Centred') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ketones_change.pdf'),
  plot = ketones_change_plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Visualise BMI change #
# ================================ #

BMI_change_plot <- ggplot(qt_ketones_both,
                          aes(x = group, y = BMI_change, colour = group)) +
  geom_jitter(size = 2.6, width = 0.2, alpha = 0.7) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = NULL, y = 'BMI Change - Centred') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BMI_change.pdf'),
  plot = BMI_change_plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - ASRS ---------------------

# ROBUST REGRESSION APPROACH #
# ================================== #

# Model 1 #
# ------------------ #

# Define predictors matrix
ASRS_keto_X <- model.matrix(~ ASRS_1 + ketones_change + ketones_pre, 
                            data = qt_ketones)[ , -1]

# Define regression model
ASRS_keto <- WRS::tshdreg(ASRS_keto_X, qt_ketones$ASRS_2, xout=TRUE, 
                          iter = 10000, outfun=outpro, corfun=pbcor, WARN = FALSE)
ASRS_keto$coef
ASRS_keto$Strength.Assoc
ASRS_keto$Explanatory.Power

# Omnibus test
ASRS_keto_omni <- WRS::regtestMC(ASRS_keto_X, 
                                 qt_ketones$ASRS_2, 
                                 regfun = tshdreg, nboot = 599, alpha = 0.05, 
                                 plotit = TRUE, xout = TRUE, outfun = outpro) 
ASRS_keto_omni

# Model 1b #
# ------------------- #

# Model matrix
ASRS_keto_X1b <- model.matrix(~ ASRS_1, data = qt_ketones)[ , -1, drop = FALSE]

# Define regression model
ASRS_keto_1b <- WRS::tshdreg(ASRS_keto_X1b, qt_ketones$ASRS_2, xout=TRUE, 
                             iter = 10000, outfun=outpro, corfun=pbcor, WARN = FALSE)
ASRS_keto_1b$coef
ASRS_keto_1b$Strength.Assoc
ASRS_keto_1b$Explanatory.Power

# Omnibus test
ASRS_keto_1b_omni <- WRS::regtestMC(ASRS_keto_X1b, qt_ketones$ASRS_2, 
                                    regfun = tshdreg, nboot = 599, alpha = 0.05, 
                                    plotit = TRUE, xout = TRUE, outfun = outpro) 
ASRS_keto_1b_omni

# Model 2 #
# ------------------- #

# Define predictors matrix
ASRS_keto_X2 <- model.matrix(~ ASRS_1 + ketones_change + 
                               ketones_pre + BMI_change + BMI_1, 
                             data = qt_ketones)[ , -1]

# Define regression model
ASRS_keto_2 <- WRS::tshdreg(ASRS_keto_X2, qt_ketones$ASRS_2, xout=TRUE, 
                            iter = 10000, outfun=outpro, corfun=pbcor, WARN = FALSE)
ASRS_keto_2$coef
ASRS_keto_2$Strength.Assoc
ASRS_keto_2$Explanatory.Power

# Omnibus test
ASRS_keto_2_omni <- WRS::regtestMC(ASRS_keto_X2, qt_ketones$ASRS_2, 
                                   iter = 10000, regfun = tshdreg, nboot = 599, 
                                   alpha = 0.05, plotit = TRUE, xout = TRUE, 
                                   outfun = outpro, SEED = TRUE) 
ASRS_keto_2_omni

# Test coefficients
ASRS_keto_2_coeffs <- regciMC_c(ASRS_keto_X2, qt_ketones$ASRS_2,
                                regfun = WRS::tshdreg, nboot = 599, 
                                iter = 100, alpha = 0.05, plotit = FALSE, pr = FALSE,
                                null.val = NULL, method = 'hoch', 
                                xlab = 'Predictor 1', ylab = 'Predictor 2', 
                                xout = TRUE, outfun = WRS::outpro, SEED = TRUE)
ASRS_keto_2_coeffs

# Model 3 #
# ------------------- #

# Define predictors matrix
ASRS_keto_X3 <- model.matrix(~ ASRS_1 + BMI_change + BMI_1,
                             data = qt_ketones)[ , -1]

# Define regression model
ASRS_keto_3 <- WRS::tshdreg(ASRS_keto_X3, qt_ketones$ASRS_2, xout=TRUE, 
                            outfun=outpro, corfun=pbcor, WARN = FALSE,
                            iter = 10000)
ASRS_keto_3$coef
ASRS_keto_3$Strength.Assoc
ASRS_keto_3$Explanatory.Power

# Omnibus test
ASRS_keto_3_omni <- WRS::regtestMC(ASRS_keto_X3, qt_ketones$ASRS_2, 
                                   iter = 10000, regfun = tshdreg, nboot = 599, 
                                   alpha = 0.05, plotit = TRUE, xout = TRUE, 
                                   outfun = outpro, corfun = pbcor) 
ASRS_keto_3_omni

# Test coefficients
ASRS_keto_3_coeffs <- regciMC_c(ASRS_keto_X3, qt_ketones$ASRS_2,
                                regfun = WRS::tshdreg, nboot = 599, 
                                iter = 100, alpha = 0.05, plotit = FALSE, pr = FALSE,
                                null.val = NULL, method = 'hoch', 
                                xlab = 'Predictor 1', ylab = 'Predictor 2', 
                                xout = TRUE, outfun = WRS::outpro, SEED = TRUE)
ASRS_keto_3_coeffs

# PLOTS #
# ===================== #

# ASRS_2 ~ ketone changes #
# ------------------------------ #

# Control variables matrix
ASRS_keto_X2_ctrl_ASRS <- ASRS_keto_X2[, c('ASRS_1', 'ketones_pre', 
                                           'BMI_change', 'BMI_1'), drop = FALSE]

# Coefficients of the control variables for ASRS_2
fit_y_X2 <- WRS::tshdreg(ASRS_keto_X2_ctrl_ASRS, qt_ketones$ASRS_2,
                         xout = TRUE, iter = 10000,
                         outfun = outpro, corfun = pbcor, WARN = FALSE)

# Remove the impact of the control variables from ASRS_2
y_X2_res <- qt_ketones$ASRS_2 - as.numeric(cbind(1, ASRS_keto_X2_ctrl_ASRS) %*% fit_y_X2$coef)

# Coefficients of the control variables for the mediator (ketones_change)
fit_x_X2 <- WRS::tshdreg(ASRS_keto_X2_ctrl_ASRS, ASRS_keto_X2[, 'ketones_change'],
                         xout = TRUE, iter = 10000,
                         outfun = outpro, corfun = pbcor, WARN = FALSE)

# Remove the impact of the control variables from the mediator (ketones_change)
x_X2_res <- ASRS_keto_X2[, 'ketones_change'] - as.numeric(cbind(1, ASRS_keto_X2_ctrl_ASRS) %*% fit_x_X2$coef)

# Plot
ASRS_X2_slope <- ASRS_keto_2$coef[3]
ASRS_X2_intercept <- ASRS_keto_2$coef[1]

tibble(x_res = x_X2_res, y_res = y_X2_res) %>% 
  ggplot(aes(x = x_res, y = y_res)) +
  geom_point(alpha = 0.75, size = 2, colour = pal['KD']) +
  geom_abline(intercept = 0, slope = ASRS_X2_slope,
              colour = 'black', linewidth = 0.9) +
  labs(x = 'Ketones change (residuals)',
       y = 'ASRS-2 (residuals)',
       title = 'ASRS Posttest as a Function of Ketone Changes') +
  theme_minimal()

# ASRS_2 ~ BMI changes #
# ------------------------------ #

ASRS_keto_X2_ctrl_BMI <- ASRS_keto_X2[, c('ASRS_1', 'ketones_pre', 
                                          'ketones_change', 'BMI_1'), drop = FALSE]

# Residualise ASRS_2 on the controls
fit_y_X2_BMI <- WRS::tshdreg(ASRS_keto_X2_ctrl_BMI, qt_ketones$ASRS_2,
                             xout = TRUE, iter = 10000,
                             outfun = outpro, corfun = pbcor, WARN = FALSE)
y_X2_res_BMI <- qt_ketones$ASRS_2 - as.numeric(cbind(1, ASRS_keto_X2_ctrl_BMI) %*% fit_y_X2_BMI$coef)

# Residualise BMI_change on the controls
fit_x_X2_BMI <- WRS::tshdreg(ASRS_keto_X2_ctrl_BMI, ASRS_keto_X2[, 'BMI_change'],
                             xout = TRUE, iter = 10000,
                             outfun = outpro, corfun = pbcor, WARN = FALSE)
x_X2_res_BMI <- ASRS_keto_X2[, 'BMI_change'] - as.numeric(cbind(1, ASRS_keto_X2_ctrl_BMI) %*% fit_x_X2_BMI$coef)

# Pull slope & intercept for BMI_change from the full model
ASRS_X2_slope_BMI <- ASRS_keto_2$coef[5]
ASRS_X2_intercept_BMI <- ASRS_keto_2$coef[1]

# Plot
tibble(x_res = x_X2_res_BMI, y_res = y_X2_res_BMI) %>% 
  ggplot(aes(x = x_res, y = y_res)) +
  geom_point(alpha = 0.75, size = 2, colour = pal['KD']) +
  geom_abline(intercept = 0, slope = ASRS_X2_slope_BMI,
              colour = 'black', linewidth = 0.9) +
  labs(x = 'BMI change (residuals)',
       y = 'ASRS-2 (residuals)',
       title = 'ASRS Posttest as a Function of BMI Changes') +
  theme_minimal()

# P-VALUE ADJUSTMENT ------------------

lab_key <- c(
  Qa = 'group',
  Qb = 'session',
  Qab = 'group:session'
)

# Collect bwtrims
bwtrim_list <- list(
  ASRS = ASRS_tranova,
  PSQI = PSQI_tranova,
  BDI = BDI_tranova
)

# Collect glmmTMB models
glmm_list <- list(
  ACC = taskswitch_2_er,
  RT = taskswitch_1a7_rt
)

# Extract bwtrim results
bwtrim_dfs <- imap(bwtrim_list, ~ {
  df <- tidy_WRS2(.x,
                  lab_key,
                  p.adjust = FALSE)
  df[['response']] <- .y 
  df
}) %>%
  list_rbind()

# Extract glmmTMB results
glmm_dfs <- imap(glmm_list, ~ {
  df <- broom::tidy(car::Anova(.x, type = 'III'))
  df[['response']] <- .y 
  df
}) %>%
  list_rbind() %>%
  rename(df1 = df)

# Combine all extracted results together
combined_dfs <- bind_rows(bwtrim_dfs, glmm_dfs)

# GROUP ONLY #
# ====================== #

effects_combined_group <- combined_dfs %>%
  filter(term == 'group') %>%
  relocate(response) %>%
  # Adjust p-values
  mutate(p.adj = p.adjust(p.value, 'BH')) %>%
  mutate(across(c(statistic, p.value, p.adj), ~ round(.x, 3)))
effects_combined_group

# SESSION ONLY #
# ====================== #

effects_combined_session <- combined_dfs %>%
  filter(term == 'session') %>%
  relocate(response) %>%
  # Adjust p-values
  mutate(p.adj = p.adjust(p.value, 'BH')) %>%
  mutate(across(c(statistic, p.value, p.adj), ~ round(.x, 3)))
effects_combined_session

# INT ONLY #
# ======================== #

effects_combined_int <- combined_dfs %>%
  filter(term %in% c('group:session', 'group:session:task_transition', 'group:session:congruence')) %>%
  relocate(response) %>%
  # Adjust p-values across all models' interactions
  mutate(p.adj = p.adjust(p.value, 'BH')) %>%
  mutate(across(c(statistic, p.value, p.adj), ~ round(.x, 3)))
print(effects_combined_int, n = Inf)

# Sensitivity change scores ASRS #
# -------------------------------- #

# Recreate the interaction table with the case-specific ASRS adjustment
effects_combined_change <- combined_dfs %>%
  filter(term %in% c('group:session', 'group:session:task_transition', 'group:session:congruence')) %>%
  relocate(response) %>%
  # Case-specific adjustment 
  mutate(
    statistic = case_when(
      response == 'ASRS' & term == 'group:session' ~ 2.251634,
      .default = statistic
    ),
    df1 = case_when(
      response == 'ASRS' & term == 'group:session' ~ NA_real_,
      .default = df1
    ),
    df2 = case_when(
      response == 'ASRS' & term == 'group:session' ~ NA_real_,
      .default = df2
    ),
    p.value = case_when(
      response == 'ASRS' & term == 'group:session' ~ 0.03314199,
      .default = p.value
    )
  ) %>%
  # Adjust p-values
  mutate(p.adj = p.adjust(p.value, 'BH')) %>%
  mutate(across(c(statistic, p.value, p.adj), ~ round(.x, 3)))
print(effects_combined_change, n = Inf)