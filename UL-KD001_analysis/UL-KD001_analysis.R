# DESCRIPTION #
# ===================== #

# Analysing UL-KD001 data. 

# TASKS #
# ===================== #

# - Mediation

# NOTES #
# ===================== #

# GENERALLY:
# --------------------- #

# - The task-switching data were analysed in two different ways 
# (mixed-effects models and regular aggregation analysis) so as to 
# see if the two different approaches would converge on the same substantive 
# conclusions. 

# - All analyses were performed in a robust way by default (trimmed means + 
# Winsorised variance), excluding the mixed-effects analysis. Some of the robust 
# alternatives also utilised bootstrapping and / or comparisons of other quantiles 
# such as median. All MAIN ('planned') analyses were also double-checked using 
# the conventional method (means + regular variance). 

# - When suitable / desirable, other methods were used too, e.g. ranked-based 
# or variance-based methods. 

# - 'conventional' = 'regular' = 'classic'
# - 'er' = 'acc' ... used rather interchangeably in names
# - '(c)variance' = 'complex variance' ... variance factors + predicted values
# - '(e)variance' = 'extra variance' ... additional variance factors (no predicted values)
# - '(ce)variance' = combination of the two above

# - Bayesian factors for main effects in the case of models with interactions
# are calculated based on the Type III Sum of Squares principle. This is in 
# contrast to the approach in JASP. 

# TASK-SWITCHING RT | MIXED-MODELS:
# ----------------------------------- #

# - The overall number of trials per group relatively balanced;
# no convincing differences despite a significant group * session interaction 
# in the case of repeat trials regardless of task conditioning (parity / magnitude). 

# - Somewhat of a group * session interaction in the case of switch trials 
# when collapsing over task (parity / magnitude) --> probably still not a point of 
# concern; controlling for trial proportions not a significant predictor of RT 
# in the mixed-effects models. 

# - diagnose_model() sometimes throws an error because of DHARMa;
# it doesn't seem to affect the diagnostics, however. 

# - None of the lmers converged except lmer log 4c
# - glmmTMB log 5 doesn't converge
# - glmmTMB with Gamma (link = 'identity') didn't run at all --> removed
# - Modelling dispersion (variance) seems to improve the fit of the MMs
# - glmmTMB log 19c, 10b, Gamma log 7g best fitting models (cross-validated)
# --- 19c: only random intercepts for participants and sessions (nested)
# --- 10b: random slopes and no nested session intercepts
# --- 7g: random slopes and no nested session intercepts, gamma family with
# log-link on log RT (double-log); perhaps the best fitting model overall 
# (residual fit), but difficult to interpret --> used only as a sensitivity check
# --- In terms of fitting parameters (e.g. AIC, BIC), 10b is by far the best

# - Controlling for error rate (at session but also participant level) 
# did not affect the main coefficients almost at all (+ perhaps slightly worse 
# fit as indicated by residuals); (currently) statistically plausible to test error 
# (as a covariate) at session level only when sessions are nested in participants.

# - brms Bayesian model takes too long to compute --> removed

# TASK-SWITCHING ER | MIXED-MODELS:
# ----------------------------------- #

# - glmmTMB 2 doesn't converge because of the RT covariates; same problems when 
# adding the covariates to simpler models, e.g. glmmTMB 5# (already removed)
# - glmmTMB 5 best fitting model (cross-validated)

# TASK-SWITCHING | AGGREGATED:
# -------------------------------- #

# - Only cue-switch trials were included in the aggregation analysis. Otherwise,
# one would have to collapse over cue switch and cue repeat trials. That may lead
# to more noise in the data, especially as each participant had a slightly 
# different number of each trial type. However, it is not clear if the substantive 
# conclusion would change; most likely not as the mixed-effects (cue 
# transition modelled there) and aggregation analyses converged on 
# the same conclusion. Another reason for using only cue-switch trials is that 
# the difference between switch and repeat trials is a clean estimate 
# of the switch cost. 

# MEDIATION REGRESSIONS:
# ------------------------------- #

# - No advantage of using the TSTS estimator instead of the Theil-Sen estimator. 

# CHANGE-BASELINE REGRESSIONS:
# ------------------------------- #

# - Blomqvist correction not possible (coeffs above 1); may have to do 
# with unequal variances (assumption of this correction)
# - Spearman and Pearson coefficients highly comparable

# LIBRARIES #
# ===================== #

# The order of libraries somewhat matters!

library(BayesFactor)
library(scales)
library(ggbeeswarm)
library(ARTool)
library(WRS)
library(WRS2)
library(afex)
library(lmtest)
library(DHARMa)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(robustlmm)
library(broom.mixed)
library(emmeans)
library(broom)
library(performance)
library(cv)
library(clubSandwich)
library(effectsize)
library(MatchIt)
library(cobalt)
library(correlation)
library(estimatr)
source('analysis_helpers.R')

# SESSION SETTINGS #
# ===================== #

# Set contrasts
options(contrasts = c('contr.sum', 'contr.poly'))

# Palette
# - Should be differentiable even by people with a colour-perception disorder
pal <- c(CD = '#600985', KD = '#004445')

# Dodge
dodge_tsmm <- 0.2
dodge_w <- 0.6
dodge_w_BMI <- 0.2

# Default directory for plots
plot_directory <- './plots'
if (!dir.exists(plot_directory)) {
  dir.create(plot_directory, recursive = TRUE)
}

table_directory <- './tables'
if (!dir.exists(table_directory)) {
  dir.create(table_directory, recursive = TRUE)
}

# Detect number of logical processors for parallel processing
n_cores = max(1L, as.integer(floor(parallel::detectCores())))

# LOAD DATA -------------------------

# QUESTIONNAIRES #
# ================================ #

# Load questionnaire RDS
questionnaires_data <- readRDS('~/UL-KD001/UL-KD001_data/LBI_clean/UL-KD001_questionnaires_pseudo.rds')

# Adjust columns
questionnaires_ready <- questionnaires_data %>%
  mutate(
    group = factor(group, levels = c('CD', 'KD')),
    cohort = factor(cohort, levels = c('apr_24', 'sep_24', 'jan_25', 'cd_25')),
  ) %>%
  dplyr::select(participant_id, group, cohort, age, BMI_pre, BMI_post,
         ASRS_total_pre, ASRS_total_post, BDI_total_pre, BDI_total_post, PSQI_total_pre,
         PSQI_total_post)

# TASK-SWITCHING #
# ================================ #

# Load task-switching RDS 
taskswitch_data <- readRDS('~/UL-KD001/UL-KD001_data/LBI_clean/UL-KD001_taskswitch_pseudo.rds')

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
    cue_transition = factor(cue_transition, levels = c('repeat' , 'switch'))
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
                                                       'BMI'),
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

age_diff <- wilcox.test(x = questionnaires_long %>% filter(group == 'CD',
                                                           session == 1) %>% pull(age),
                        y = questionnaires_long %>% filter(group == 'KD',
                                                           session == 1) %>% pull(age),
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

BMI_2 <- afex_plot(
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
  filename = file.path(plot_directory, 'BMI_2.pdf'),
  plot = BMI_2,
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
  filename = file.path(plot_directory, 'gg_BMI_with.pdf'),
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
                           pull(BMI),
                         questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           pull(BMI),
                         tr = 0.2)
BMI_session_KD

# Session change CD
BMI_session_CD <- yuend(questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                           pull(BMI),
                         questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           pull(BMI),
                         tr = 0.2)
BMI_session_CD

# Effect size of session changes
BMI_session_eff <- bw.es.B(2, 2, BMI_post_tr_list, 
                            tr = 0.2, POOL = FALSE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
BMI_session_eff

# TASK-SWITCHING PREPARATION -----------------------------

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
    !(participant_id %in% taskswitch_exclude$participant_id)
  ) %>%
  mutate(log_rt = log(response_time))

# Cleaning for ER analysis
taskswitch_crispy_er <- taskswitch_ready %>%
  filter(
    # Filter out too fast and too slow trials
    response_time > 300 & response_time < 3000,
    # Filter out participants with too high an error rate
    !(participant_id %in% taskswitch_exclude$participant_id)
  ) %>%
  mutate(log_rt = log(response_time))

# Examine filtering results #
# -------------------------------- #

# Check number of participants after filtering
taskswitch_participants <- taskswitch_crispy_rt %>%
  group_by(group) %>%
  summarise(count = n_distinct(participant_id)) %>%
  ungroup()
taskswitch_participants

# Create a tibble crossing over task, task transitions, congruence and cue transitions
# to get a better understanding of the data
taskswitch_balance <- taskswitch_crispy_rt %>%
  group_by(task_transition, congruence, cue_transition, task) %>%
  summarise(count = n()) %>%
  ungroup()
print(taskswitch_balance, n = Inf)

# Excluded RT trials because being too fast or too slow (RT analyses)
task_switch_speed <- taskswitch_ready %>%
  filter(error == FALSE,
         !(participant_id %in% taskswitch_exclude$participant_id)) %>%
  summarise(speed_off_prop = mean(response_time < 300 | response_time > 3000)) %>%
  pull(speed_off_prop)
task_switch_speed

# Excluded RT trials because being too fast or too slow (ACC analyses)
task_switch_speed_acc <- taskswitch_ready %>%
  filter(!(participant_id %in% taskswitch_exclude$participant_id)) %>%
  summarise(speed_off_prop = mean(response_time < 300 | response_time > 3000)) %>%
  pull(speed_off_prop)
task_switch_speed_acc

# - Trial comparison 1 ---------------------------

# Descriptives
comparison_taskswitch <- taskswitch_crispy_rt %>%
  group_by(session, group, task_transition, congruence, cue_transition, task) %>%
  summarise(count_corrected = n() / n_distinct(participant_id)) %>%
  ungroup() %>%
  pivot_wider(values_from = count_corrected,
              names_from = group)
comparison_taskswitch

# Form subject-level tibble
comparison_taskswitch_subjects <- taskswitch_crispy_rt %>%
  group_by(participant_id, group, session, task_transition, 
           congruence, cue_transition, task) %>%
  summarise(n_trials = n()) %>%
  unite(
    col = 'trial_combination',
    task_transition, congruence, cue_transition, task,
    sep = '_',
    remove = FALSE
  ) %>%
  ungroup() %>%
  mutate(trial_combination = as.factor(trial_combination))
glimpse(comparison_taskswitch_subjects)

# Session 1 distributions
ggplot(comparison_taskswitch_subjects %>% filter(session == 1),
       aes(x = n_trials, fill = group)) +
  geom_density(alpha = 0.5, adjust = 1) +
  facet_wrap(
    ~ task_transition + congruence + cue_transition + task,
    nrow = 2, ncol = 6,
    labeller = label_both
  ) +
  scale_fill_manual(values = pal) +
  labs(
    x = 'Number of Trials (Session 1)',
    y = 'Density',
    fill = 'Group'
  ) +
  theme_minimal()

# Session 2 distributions
ggplot(comparison_taskswitch_subjects %>% filter(session == 2), 
       aes(x = n_trials, fill = group)) +
  geom_density(alpha = 0.5, adjust = 1) +
  facet_wrap(
    ~ task_transition + congruence + cue_transition + task,
    nrow = 2, ncol = 6,
    labeller = label_both
  ) +
  scale_fill_manual(values = pal) +
  labs(
    x = 'Number of Trials (Session 2)',
    y = 'Density',
    fill = 'Group'
  ) +
  theme_minimal()

# Test for differences (repeat task transitions)
taskswitch_comp_model_1 <- aov_4(n_trials ~ group * session * trial_combination +
                                   (session * trial_combination | participant_id),
                                 data = comparison_taskswitch_subjects %>% 
                                   filter(task_transition == 'repeat'),
                                 anova_table = list(es = 'pes'))
taskswitch_comp_model_1
summary(taskswitch_comp_model_1$Anova)

# Test for differences (switch task transitions)
taskswitch_comp_model_2 <- aov_4(n_trials ~ group * session * trial_combination + 
                                   (session * trial_combination | participant_id),
                                 data = comparison_taskswitch_subjects %>% 
                                   filter(task_transition == 'switch'),
                                 anova_table = list(es = 'pes'))
taskswitch_comp_model_2
summary(taskswitch_comp_model_2$Anova)

# Repeats plotting #
# ------------------------- #

comparison_taskswitch_subjects_ggtibble_rep <- comparison_taskswitch_subjects %>% 
  filter(task_transition == 'repeat') %>%
  group_by(participant_id, session, group) %>%
  summarise(trial_count = sum(n_trials)) %>%
  ungroup()

ggplot(
  comparison_taskswitch_subjects_ggtibble_rep,
  aes(x = session, y = trial_count, fill = group) 
) +
  stat_summary(
    fun = mean,
    geom = 'col',
    position = position_dodge(width = 0.8),
    alpha = 0.9,
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = 'errorbar',
    aes(colour = group),  
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.7
  ) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +  
  labs(
    x = 'Session',
    y = 'Mean Number of Repeat Trials (±SE)',
    fill = 'Group',
    colour = 'Group'
  ) +
  theme_minimal(base_size = 14)

# Switch plotting #
# ------------------------- #

comparison_taskswitch_subjects_ggtibble_switch <- comparison_taskswitch_subjects %>% 
  filter(task_transition == 'switch') %>%
  group_by(participant_id, session, group) %>%
  summarise(trial_count = sum(n_trials)) %>%
  ungroup()

ggplot(
  comparison_taskswitch_subjects_ggtibble_switch,
  aes(x = session, y = trial_count, fill = group) 
) +
  stat_summary(
    fun = mean,
    geom = 'col',
    position = position_dodge(width = 0.8),
    alpha = 0.9,
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = 'errorbar',
    aes(colour = group),  
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.7
  ) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +  
  labs(
    x = 'Session',
    y = 'Mean Number of Switch Trials (±SE)',
    fill = 'Group',
    colour = 'Group'
  ) +
  theme_minimal(base_size = 14)

# - Trial comparison 2 ------------------------

# NOTE: Same as comparison 1, just collapsing over task

# Descriptives
comparison_taskswitch_2 <- taskswitch_crispy_rt %>%
  group_by(session, group, task_transition, congruence, cue_transition) %>%
  summarise(count_corrected = n() / n_distinct(participant_id)) %>%
  ungroup() %>%
  pivot_wider(values_from = count_corrected,
              names_from = group)
comparison_taskswitch_2

# Form subject-level tibble
comparison_taskswitch_subjects_2 <- taskswitch_crispy_rt %>%
  group_by(participant_id, group, session, task_transition, 
           congruence, cue_transition) %>%
  summarise(n_trials = n()) %>%
  unite(
    col = 'trial_combination',
    task_transition, congruence, cue_transition,
    sep = '_',
    remove = FALSE
  ) %>%
  ungroup() %>%
  mutate(trial_combination = as.factor(trial_combination))
glimpse(comparison_taskswitch_subjects_2)

# Distributions
ggplot(comparison_taskswitch_subjects_2, aes(x = n_trials, fill = group)) +
  geom_density(alpha = 0.5, adjust = 1) +
  facet_wrap(
    ~ session + task_transition + congruence + cue_transition,
    nrow = 2, ncol = 6,
    labeller = label_both
  ) +
  scale_fill_manual(values = pal) +
  labs(
    x = 'Number of Trials',
    y = 'Density',
    fill = 'Group'
  ) +
  theme_minimal()

# Test for differences (repeat task transitions)
taskswitch_comp_2_model_1 <- aov_4(n_trials ~ group * session * trial_combination + 
                                     (session * trial_combination | participant_id),
                                 data = comparison_taskswitch_subjects_2 %>% 
                                   filter(task_transition == 'repeat'),
                                 anova_table = list(es = 'pes'))
taskswitch_comp_2_model_1
summary(taskswitch_comp_2_model_1$Anova)

# Test for differences (switch task transitions)
taskswitch_comp_2_model_2 <- aov_4(n_trials ~ group * session * trial_combination + 
                                     (session * trial_combination | participant_id),
                                 data = comparison_taskswitch_subjects_2 %>% 
                                   filter(task_transition == 'switch'),
                                 anova_table = list(es = 'pes'))
taskswitch_comp_2_model_2
summary(taskswitch_comp_2_model_2$Anova)

# Switch plotting #
# ------------------------- #

comparison_taskswitch_subjects_ggtibble_switch_2 <- comparison_taskswitch_subjects_2 %>% 
  filter(task_transition == 'switch') %>%
  group_by(participant_id, session, group, trial_combination) %>%
  summarise(trial_count = sum(n_trials)) %>%
  ungroup()

ggplot(
  comparison_taskswitch_subjects_ggtibble_switch_2,
  aes(x = trial_combination,
      y = trial_count, fill = group) 
) +
  facet_wrap(~ session) +
  stat_summary(
    fun = mean,
    geom = 'col',
    position = position_dodge(width = 0.8),
    alpha = 0.9,
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = 'errorbar',
    aes(colour = group),  
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.7
  ) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +  
  labs(
    x = 'Session',
    y = 'Mean Number of Switch Trials (±SE)',
    fill = 'Group',
    colour = 'Group'
  ) +
  # Cue transitions are always switch in the case 
  # switch trials (task transition) --> only congruence left
  scale_x_discrete(labels = c('Congruent', 'Incongruent')) + 
  theme_minimal()

# - ER tibbles ----------------------

glimpse(taskswitch_crispy_er)

# Tibble with trial proportions
trial_proportions_er <- taskswitch_crispy_er %>%
  group_by(group, participant_id, session) %>%
  summarise(
    prop_switch = mean(task_transition == 'switch'),
    prop_incong = mean(congruence == 'incongruent'),
    prop_cue_switch = mean(cue_transition == 'switch')
  ) %>%
  ungroup() %>%
  # Centre
  mutate(across(starts_with(('prop_')), ~ as.numeric(scale(.x, scale = FALSE))))
glimpse(trial_proportions_er)

# Tibble with errors for RT mixed-models (covariates)
error_proportions <- taskswitch_crispy_er %>%
  group_by(participant_id, group, session) %>%
  summarise(
    error_prop = mean(error == TRUE, na.rm = TRUE), .groups = 'drop') %>%
  group_by(participant_id) %>%
  mutate(
    bs_error_prop = mean(error_prop, na.rm = TRUE),
    ss_error_prop = error_prop - bs_error_prop) %>%
  ungroup() %>%
  dplyr::select(-error_prop) %>%
  mutate(bs_error_prop = as.numeric(scale(bs_error_prop, scale = FALSE)))
glimpse(error_proportions)

# - RT tibbles ----------------------------

glimpse(taskswitch_crispy_rt)

# Tibble with trial proportions
trial_proportions_rt <- taskswitch_crispy_rt %>%
  group_by(group, participant_id, session) %>%
  summarise(
    prop_switch = mean(task_transition == 'switch'),
    prop_incong = mean(congruence == 'incongruent'),
    prop_cue_switch = mean(cue_transition == 'switch')
  ) %>%
  ungroup() %>%
  # Centre
  mutate(across(starts_with(('prop_')), ~ as.numeric(scale(.x, scale = FALSE))))
glimpse(trial_proportions_rt)

# Check distributions #
# ------------------------------ #

inspect_subject_dist(
  taskswitch_crispy_rt %>% filter(session == 1),
  n = 10,
  pid_col = 'participant_id',
  col_dist = 'log_rt',
  group_col = 'group',
  session_col = 'session',
  wrap = c('session'),
  pal = pal,
  robust = 'MAD',
  mad_mult = 2,
  seed = NULL
)

# RT decomposition #
# --------------------------- #

# Trial- and session-level pieces
rt_decomposed_ssws <- taskswitch_crispy_rt %>%
  group_by(participant_id) %>%
  mutate(bs_log_rt_tmp = mean(log_rt, na.rm = TRUE)) %>%  # Person mean over trials
  group_by(participant_id, group, session) %>%
  mutate(
    m_is = mean(log_rt, na.rm = TRUE), # Session mean
    m_is = ifelse(is.nan(m_is), NA_real_, m_is),
    ws_log_rt = log_rt - m_is, # Trial deviation from session mean
    ss_log_rt = m_is - bs_log_rt_tmp # Session mean minus person mean
  ) %>%
  ungroup() %>%
  dplyr::select(-log_rt, -m_is, -bs_log_rt_tmp)
glimpse(rt_decomposed_ssws)

# Between-person (grand-mean centred; trial-weighted)
rt_decomposed_bs <- taskswitch_crispy_rt %>%
  group_by(participant_id) %>%
  summarise(
    bs_raw = mean(log_rt, na.rm = TRUE),
    n_trials = sum(!is.na(log_rt)),
    .groups = 'drop'
  ) %>%
  mutate(
    grand_log_rt = weighted.mean(bs_raw, w = n_trials, na.rm = TRUE),
    bs_log_rt = bs_raw - grand_log_rt # Centred
  ) %>%
  dplyr::select(participant_id, bs_log_rt)
glimpse(rt_decomposed_bs)

# Final table
rt_decomposed <- rt_decomposed_ssws %>%
  left_join(rt_decomposed_bs, by = 'participant_id')
glimpse(rt_decomposed)

# - Data for mixed-effects ER ------------------------

taskswitch_mixed_er <- taskswitch_crispy_er %>%
  left_join(trial_proportions_er, by = c('participant_id', 'session', 'group')) %>%
  left_join(rt_decomposed %>% dplyr::select(-c(response_time, error)),
            by = c('participant_id', 'session', 'group', 'cohort', 'trial',
                                  'date_startdate', 'block_count',
                                  'trial_sequence', 'cue_color', 'task', 'congruence',
                                  'digit', 'task_transition', 'digit_transition', 
                                  'cue_transition')) %>%
  mutate(response_correct = ifelse(error == FALSE, 1, 0),
         response_incorrect = ifelse(error == TRUE, 1, 0),
         trial_sequence_factor = glmmTMB::numFactor(trial_sequence))
glimpse(taskswitch_mixed_er)

# - Data for mixed-effects RT ----------------------

# Tibble for mixed effects
taskswitch_mixed_rt <- taskswitch_crispy_rt %>%
  left_join(error_proportions, by = c('participant_id', 'session', 'group')) %>%
  left_join(trial_proportions_rt, by = c('participant_id', 'session', 'group')) %>%
  mutate(trial_sequence_factor = glmmTMB::numFactor(trial_sequence))
glimpse(taskswitch_mixed_rt)

# - Assess distributions ------------------------

# NOTE: The distributions of binomial (ER / ACC) responses cannot be assessed directly 
# like the RT response --> model fit is examined instead

# Plot density of all trial combinations (regular RT)
ggplot(taskswitch_mixed_rt, aes(x = response_time, fill = group)) +
  geom_density(alpha = 0.6, adjust = 1) + 
  facet_wrap(~ session + task_transition + congruence,
             nrow = 2, ncol = 4) + 
  scale_fill_manual(values = pal) +
  labs(
    x = 'RT (ms)',
    y = 'Density',
    fill = 'Group'
  ) +
  theme_minimal()

# Plot density of all trial combinations (log RT)
ggplot(taskswitch_mixed_rt, aes(x = log_rt, fill = group)) +
  geom_density(alpha = 0.6, adjust = 1) + 
  facet_wrap(~ session + task_transition + congruence,
             nrow = 2, ncol = 4) + 
  scale_fill_manual(values = pal) +
  labs(
    x = 'RT (ms)',
    y = 'Density',
    fill = 'Group'
  ) +
  theme_minimal()

# - Assess variances ---------------------------

# Tibble of subject-level dispersions for each trial type
subject_dispersion <- taskswitch_mixed_rt %>%
  group_by(participant_id, group, session, task_transition, congruence) %>%
  summarise(dispersion = var(log_rt, na.rm = TRUE)) %>%
  ungroup()

# Visualise dispersion
ggplot(subject_dispersion, aes(x = group, y = dispersion, fill = group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .15, alpha = .4, size = 1) +
  facet_grid(session ~ task_transition + congruence) + 
  scale_fill_manual(values = pal) + 
  labs(title = 'Dispersion per Group, Session, Trial-Type',
       x = NULL,
       y = 'Variance (log RT)')+
  theme_minimal()

# Test for differences
dispersion_anova <- aov_4(dispersion ~ group + (session * task_transition * congruence | participant_id), data = subject_dispersion, anova_table = list(es = 'pes'))
dispersion_anova

# MIXED-EFFECTS MODELS RT -------------------

# Clean up a bit of workspace
rm(taskswitch_data, taskswitch_ready, questionnaires_data, questionnaires_ready)

# - lmer raw 1 (4-way + all task covariates + RE*) -----------------------

taskswitch_mm_1_rt <- lme4::lmer(
  response_time ~ group * session * task_transition * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (task_transition * congruence | participant_id),
  data = taskswitch_mixed_rt, REML = FALSE)
summary(taskswitch_mm_1_rt)

# Diagnostics #
# -------------------- #

taskswitch_mm_1_rt_diag <- diagnose_lmer(taskswitch_mm_1_rt)
print(taskswitch_mm_1_rt_diag, show_plots = TRUE)

# - lmer log 2 (4-way + all task covariates + RE*) ---------------------

# Same as model 1; RT just log-transformed
taskswitch_mm_2_rt <- lme4::lmer(
  log_rt ~ group * session * task_transition * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (task_transition * congruence | participant_id),
  data = taskswitch_mixed_rt, REML = FALSE)
summary(taskswitch_mm_2_rt)
isSingular(taskswitch_mm_2_rt)
  
# - lmer log 3 (4-way + all task covariates + RE1+) ---------------------

# Model 2 but simpler random-effects structure
taskswitch_mm_3_rt <- lme4::lmer(
  log_rt ~ group * session * task_transition * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (task_transition + congruence | participant_id),
  data = taskswitch_mixed_rt, REML = FALSE)
summary(taskswitch_mm_3_rt)
isSingular(taskswitch_mm_3_rt)

# Diagnostics #
# ===================== #

taskswitch_mm_3_rt_diag <- diagnose_lmer(taskswitch_mm_3_rt)
print(taskswitch_mm_3_rt_diag, show_plots = TRUE)

diagnose_model(taskswitch_mm_3_rt)

# - lmer log 4 (3-way + all task covariates + RE1+) ---------------------

taskswitch_mm_4_rt <- lme4::lmer(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (task_transition + congruence | participant_id),
  data = taskswitch_mixed_rt, REML = FALSE)
summary(taskswitch_mm_4_rt)
isSingular(taskswitch_mm_4_rt)

# Diagnostics #
# ===================== #

taskswitch_mm_4_rt_diag <- diagnose_lmer(taskswitch_mm_4_rt)
print(taskswitch_mm_4_rt_diag, show_plots = TRUE)

diagnose_model(taskswitch_mm_4_rt, lme4 = TRUE)

# - lmer log 4b (3-way + all task covariates + RE1+||) ---------------------

taskswitch_mm_4b_rt <- lme4::lmer(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (task_transition + congruence || participant_id),
  data = taskswitch_mixed_rt, REML = FALSE)
summary(taskswitch_mm_4_rt)
isSingular(taskswitch_mm_4_rt)

# - lmer log 4c (3-way + all task covariates + RE1/) ---------------------

taskswitch_mm_4c_rt <- lme4::lmer(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (1 | participant_id / session), 
  data = taskswitch_mixed_rt, REML = FALSE)
summary(taskswitch_mm_4c_rt)
isSingular(taskswitch_mm_4_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_4c_rt, lme4 = TRUE)

# - glmmTMB log 5 (4-way + variance + all task covariates + RE1/) -----------------------

taskswitch_mm_5_rt <- glmmTMB(
  log_rt ~ group * session * task_transition * congruence +
    cue_transition + prop_switch + prop_incong + prop_cue_switch +
    (1 | participant_id / session),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_5_rt)

# - glmmTMB Gamma raw 6 (4-way + variance + all task covariates + RE1/) ----------------------------

taskswitch_mm_6_rt <- glmmTMB(
  response_time ~ group * session * task_transition * congruence +
      cue_transition + prop_switch + prop_incong + prop_cue_switch +
      (1 | participant_id / session),
      dispformula = ~ session + group * task_transition,
    data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_6_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_6_rt)

# - glmmTMB Gamma raw 7 (3-way + variance + all task covariates + RE1/) ----------------------------

taskswitch_mm_7_rt <- glmmTMB(
  response_time ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (1 | participant_id / session),
    dispformula = ~ session + group * task_transition,
    data = taskswitch_mixed_rt,
    family = Gamma(link = 'log')
)
summary(taskswitch_mm_7_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7_rt)

# - glmmTMB Gamma log 7b (3-way + variance + all task covariates + RE1/) ----------------------------

taskswitch_mm_7b_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (1 | participant_id / session),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7b_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7b_rt)

# - glmmTMB Gamma log 7c (3-way + variance + 1 covariate + RE1+) ----------------------------

# Singular
taskswitch_mm_7c_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (task_transition + congruence + cue_transition | participant_id),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7c_rt)
check_singularity(taskswitch_mm_7c_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7c_rt)

# - glmmTMB Gamma log 7d (3-way + variance + 1 covariate + RE1+||) ----------------------------

taskswitch_mm_7d_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (task_transition + congruence + cue_transition || participant_id),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7d_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7d_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_7d_rt, taskswitch_mm_7c_rt)

# - glmmTMB Gamma log 7e (3-way + (e)variance + 1 covariate + RE1+||) ----------------------------

taskswitch_mm_7e_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (task_transition + congruence + cue_transition || participant_id),
  dispformula = ~ session + group * task_transition + cue_transition,
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7e_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7e_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_7d_rt, taskswitch_mm_7e_rt)

# - glmmTMB Gamma log 7f (3-way + (e)variance + 1 covariate + RE1) ----------------------------

taskswitch_mm_7f_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (1 | participant_id),
  dispformula = ~ session + group * task_transition + cue_transition,
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7f_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7f_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_7f_rt, taskswitch_mm_7e_rt)

# - glmmTMB Gamma log 7g (3-way + (ce)variance + 1 covariate + RE1+||) ---------------

taskswitch_mixed_rt <- taskswitch_mixed_rt %>%
  mutate(fitted_mgam7e = predict(taskswitch_mm_7e_rt))
glimpse(taskswitch_mixed_rt)

taskswitch_mm_7g_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (task_transition + congruence + cue_transition || participant_id),
  dispformula = ~ poly(fitted_mgam7e, 2) + 
    session + group * task_transition + cue_transition,
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7g_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7g_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_7g_rt, taskswitch_mm_7e_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_7g_rt <- cv(taskswitch_mm_7g_rt, 
                              k = 5, 
                              clusterVariables = 'participant_id', 
                              ncores = n_cores,
                              reps = 10)
summary(cv_taskswitch_mm_7g_rt)
plot(cv_taskswitch_mm_7g_rt)

# - glmmTMB Gamma log 7h (3-way + poly + 1 covariate + RE1+||) ---------------

taskswitch_mm_7h_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (task_transition + congruence + cue_transition || participant_id),
  dispformula = ~ poly(fitted_mgam7e, 2),
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7h_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7h_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_7g_rt, taskswitch_mm_7h_rt)

# - glmmTMB Gamma log 7j (3-way + (ce)variance + 1 covariate + autocorr + RE1+||) ---------------

# Not as good quantile fit as glmmTMB Gamma log 7g
taskswitch_mm_7j_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (task_transition + congruence + cue_transition || participant_id) + 
    ou(trial_sequence_factor + 0 | participant_id:session),
  dispformula = ~ poly(fitted_mgam7e, 2) + 
    session + group * task_transition + cue_transition,
  data = taskswitch_mixed_rt,
  family = Gamma(link = 'log')
)
summary(taskswitch_mm_7j_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_7j_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_7j_rt, taskswitch_mm_7g_rt)

# - glmmTMB log 8 (3-way + variance + all task covariates + RE1/) -----------------------

taskswitch_mm_8_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch +
    (1 | participant_id / session),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_8_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_8_rt)

# - glmmTMB log 9 (3-way + variance + 1 covariate + RE1) -----------------------

taskswitch_mm_9_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + cue_transition +
    (1 | participant_id),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_9_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_9_rt, taskswitch_mm_8_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_9_rt)

# Cross validation #
# ========================= #

cv_taskswitch_mm_9_rt <- cv(taskswitch_mm_9_rt, k = 10, ncores = n_cores)
summary(cv_taskswitch_mm_9_rt)

# - glmmTMB log 9b (3-way + variance + 1 covariate + RE1/) -----------------------

taskswitch_mm_9b_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + cue_transition +
    (1 | participant_id / session),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_9b_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_9b_rt, taskswitch_mm_8_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_9b_rt)

# Cross validation #
# ========================= #

cv_taskswitch_mm_9_rt <- cv(taskswitch_mm_9_rt, k = 10, ncores = n_cores)
summary(cv_taskswitch_mm_9_rt)

# - glmmTMB log 9c (3-way + variance + 1 covariate + RE+||) -----------------------

taskswitch_mm_9c_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + group * session * congruence + cue_transition +
    (task_transition + congruence + cue_transition || participant_id),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_9c_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_9c_rt, taskswitch_mm_9_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_9c_rt)

# Cross validation #
# ========================= #

cv_taskswitch_mm_9_rt <- cv(taskswitch_mm_9_rt, k = 10, ncores = n_cores)
summary(cv_taskswitch_mm_9_rt)

# - Autocorrelation check ----------------------------

# Get model residuals
set_simcodes(taskswitch_mm_9c_rt, val = 'fix')
taskswitch_mm_9c_rt_res <- simulateResiduals(taskswitch_mm_9c_rt, plot = FALSE)

# Temporal autocorrelation per subject
taskswitch_mm_9c_rt_autocorr <- taskswitch_mixed_rt %>%
  dplyr::select(participant_id, group, session, trial_sequence) %>%
  mutate(residual = taskswitch_mm_9c_rt_res$scaledResiduals) %>%
  group_by(participant_id, group, session) %>%
  summarise(
    dw = dwtest(residual ~ 1, order.by = trial_sequence)$statistic,
    p_value = round(dwtest(residual ~ 1, order.by = trial_sequence)$p.value, 3)
  ) %>%
  ungroup()

# Autocorrelation descriptives
taskswitch_mm_9c_rt_autocorr_desc <- taskswitch_mm_9c_rt_autocorr %>%
  group_by(group, session) %>%
  summarise(autocorr_prop = mean(dw < 1.5 | dw > 2.5),
            mean_dw = mean(dw, na.rm = TRUE),
            sd_dw = sd(dw, na.rm = TRUE)) %>%
  ungroup()
taskswitch_mm_9c_rt_autocorr_desc

# Plot
max_lag <- 40
acf_mm_9c_rt <- taskswitch_mixed_rt %>%
  mutate(residual = taskswitch_mm_9c_rt_res$scaledResiduals) %>%     
  group_by(participant_id, group, session) %>%                   
  nest() %>%                                                     
  mutate(acf_tbl = map(data, ~{
    ac <- acf(.x$residual, lag.max = max_lag, plot = FALSE)
    tibble(lag = as.numeric(ac$lag[-1]),                       
           acf = as.numeric(ac$acf[-1]))
  })) %>%
  dplyr::select(-data) %>%                                             
  unnest(acf_tbl)                                               

acf_mm_9c_rt_summary <- acf_mm_9c_rt %>%
  group_by(group, session, lag) %>%
  summarise(mean_acf = mean(acf, na.rm = TRUE),
            median_acf = median(acf, na.rm = TRUE),
            .groups = 'drop')

ggplot(acf_mm_9c_rt, aes(lag, acf, group = participant_id)) +
  geom_line(alpha = 0.25, linewidth = 0.3) +              
  geom_line(data = acf_mm_9c_rt_summary,                            
            aes(lag, mean_acf),
            inherit.aes = FALSE,
            colour = 'black', linewidth = 1) +
  facet_grid(group ~ session) +
  geom_hline(yintercept = c(-.05, .05), linetype = 'dotted') +  # ±0.05 guides
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Lag', y = 'ACF',
       title = 'Participant ACFs with grand-mean overlay') +
  theme_minimal()

# - glmmTMB log 9d (3-way + custom variance + 1 covariate + RE+||) -------------------

taskswitch_mm_9d_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (task_transition + congruence + cue_transition || participant_id),
  dispformula = ~ session * group + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_9d_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_9d_rt, taskswitch_mm_9c_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_9d_rt)

# - glmmTMB log 9e (3-way + (c)variance + 1 covariate + RE+||) -----------------------

taskswitch_mixed_rt <- taskswitch_mixed_rt %>%
  mutate(fitted_m9c = predict(taskswitch_mm_9c_rt))
glimpse(taskswitch_mixed_rt)

taskswitch_mm_9e_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + cue_transition +
    (task_transition + congruence + cue_transition || participant_id),
  dispformula = ~ poly(fitted_m9c, 2) + session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_9e_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_9e_rt, taskswitch_mm_9c_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_9e_rt)

# - glmmTMB log 10 (3-way + (c)variance + 1 covariate + autocorr + RE+||) -------------------

taskswitch_mm_10_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition +
    (task_transition + congruence + cue_transition || participant_id) + 
    ou(trial_sequence_factor + 0 | participant_id:session),
  dispformula = ~ poly(fitted_m9c, 2) + session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_10_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_9e_rt, taskswitch_mm_10_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_10_rt)

set_simcodes(taskswitch_mm_10_rt, val = 'fix')
taskswitch_mm_10_rt_res <- simulateResiduals(taskswitch_mm_10_rt)
taskswitch_mm_10_rt_res_agg <- recalculateResiduals(taskswitch_mm_10_rt_res, 
                                                    group = taskswitch_mixed_rt$cue_transition)
testDispersion(taskswitch_mm_10_rt_res_agg)
testCategorical(taskswitch_mm_10_rt_res, catPred = taskswitch_mixed_rt$cue_transition)

rm(taskswitch_mm_10_rt_res, taskswitch_mm_10_rt_res_agg)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_10_rt <- cv(taskswitch_mm_10_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_10_rt)
plot(cv_taskswitch_mm_10_rt)

# - glmmTMB log 10b (3-way + (ce)variance + 1 covariate + autocorr + RE+||) -------------------

# More complex dispersion formula not needed; this one is maximal 
taskswitch_mm_10b_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + cue_transition +
    (task_transition + congruence + cue_transition || participant_id) + 
    ou(trial_sequence_factor + 0 | participant_id:session),
  dispformula = ~ poly(fitted_m9c, 2) + 
    session + task_transition * group + cue_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_10b_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_10b_rt, taskswitch_mm_10_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_10b_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_10b_rt <- cv(taskswitch_mm_10b_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_10b_rt)
plot(cv_taskswitch_mm_10b_rt)

# - glmmTMB log 10c (3-way + (ce)variance + all task covariates + autocorr + RE+||) -------------------

# More complex dispersion formula not needed; this one is maximal 
taskswitch_mm_10c_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + cue_transition +
    prop_switch + prop_incong + prop_cue_switch +
    (task_transition + congruence + cue_transition || participant_id) + 
    ou(trial_sequence_factor + 0 | participant_id:session),
  dispformula = ~ poly(fitted_m9c, 2) + 
    session + task_transition * group + cue_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_10c_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_10b_rt, taskswitch_mm_10c_rt)

# - glmmTMB log 11 (3-way + 1 covariate + RE1/) -------------------

taskswitch_mm_11_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + cue_transition + 
    (1 | participant_id / session),
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_11_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_11_rt, taskswitch_mm_9_rt)

# Cross validation #
# ========================= #

cv_taskswitch_mm_11_rt <- cv(taskswitch_mm_11_rt, 
                             k = 10, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 5)
summary(cv_taskswitch_mm_11_rt)
plot(cv_taskswitch_mm_11_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_11_rt)

# - glmmTMB log 12 (3-way + variance + 1 covariate + error(2) + RE1/) -----------------------

taskswitch_mm_12_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + bs_error_prop + ss_error_prop + (1 | participant_id / session),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_12_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_12_rt, taskswitch_mm_9_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_12_rt)

# - glmmTMB log 12b (3-way + variance + 1 covariate + error(1) + RE1/) -----------------------

taskswitch_mm_12b_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + ss_error_prop + 
    (1 | participant_id / session),
  dispformula = ~ session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_12b_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_12b_rt, taskswitch_mm_9_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_12b_rt)

# - glmmTMB log 13 (3-way + custom variance + 1 covariate + RE1/) -----------------------

taskswitch_mm_13_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + group * session * congruence + cue_transition +
  (1 | participant_id / session),
  dispformula = ~ session * group + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_13_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_9_rt, taskswitch_mm_13_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_13_rt)

# Cross validation #
# ========================= #

cv_taskswitch_mm_13_rt <- cv(taskswitch_mm_13_rt, 
                             k = 10, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 3)
summary(cv_taskswitch_mm_13_rt)
plot(cv_taskswitch_mm_13_rt)

# - glmmTMB log 14 (3-way + poly + 1 covariate + RE1/) -----------------------

taskswitch_mixed_rt <- taskswitch_mixed_rt %>%
  mutate(fitted_m11 = predict(taskswitch_mm_11_rt))
glimpse(taskswitch_mixed_rt)

taskswitch_mm_14_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition +
    (1 | participant_id / session),
  dispformula = ~ poly(fitted_m11, 2),
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_14_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_14_rt, taskswitch_mm_11_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_14_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_14_rt <- cv(taskswitch_mm_14_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_14_rt)
plot(cv_taskswitch_mm_14_rt)

# Row-based
cv_taskswitch_mm_14_rt_ROW <- cv(taskswitch_mm_14_rt, 
                             k = 10,
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_14_rt_ROW)
plot(cv_taskswitch_mm_14_rt_ROW)

# - glmmTMB log 15 (3-way + linear + 1 covariate + RE1/) -----------------------

taskswitch_mixed_rt <- taskswitch_mixed_rt %>%
  mutate(fitted_m11 = predict(taskswitch_mm_11_rt))
glimpse(taskswitch_mixed_rt)

taskswitch_mm_15_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition +
    (1 | participant_id / session),
  dispformula = ~ fitted_m11,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_15_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_15_rt, taskswitch_mm_14_rt, taskswitch_mm_11_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_15_rt)

# Cross validation #
# ========================= #

cv_taskswitch_mm_15_rt <- cv(taskswitch_mm_15_rt, 
                             k = 10, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores)
summary(cv_taskswitch_mm_15_rt)
plot(cv_taskswitch_mm_15_rt)

# - glmmTMB log 16 (2-way + 1 covariate + RE1/) -----------------------

taskswitch_mm_16_rt <- glmmTMB(
  log_rt ~ session + task_transition * group + 
    congruence * group + cue_transition +
    (1 | participant_id / session),
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_16_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_14_rt, taskswitch_mm_16_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_16_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_16_rt <- cv(taskswitch_mm_16_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_16_rt)
plot(cv_taskswitch_mm_16_rt)

# Row-based
cv_taskswitch_mm_16_rt_ROW <- cv(taskswitch_mm_16_rt, 
                             k = 5,
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_16_rt_ROW)
plot(cv_taskswitch_mm_16_rt_ROW)

# - glmmTMB log 17 (3-way + custom variance + 1 covariate + RE1/) -----------------------

taskswitch_mm_17_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition +
    (1 | participant_id / session),
  dispformula = ~ poly(fitted_m11, 2) + session * group + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_17_rt)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_14_rt, taskswitch_mm_17_rt)
anova(taskswitch_mm_13_rt, taskswitch_mm_17_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_17_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_17_rt <- cv(taskswitch_mm_17_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_17_rt)
plot(cv_taskswitch_mm_17_rt)

# - glmmTMB log 18 (3-way + (c)variance + 1 covariate + error(1) + RE1/) -----------------------

taskswitch_mm_18_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + ss_error_prop + 
    (1 | participant_id / session),
  dispformula = ~ poly(fitted_m11, 2) + session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_18_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_18_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_18_rt <- cv(taskswitch_mm_18_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_18_rt)
plot(cv_taskswitch_mm_18_rt)

# - glmmTMB log 19 (3-way + (c)variance + 1 covariate + RE1/) -----------------------

taskswitch_mixed_rt <- taskswitch_mixed_rt %>%
  mutate(fitted_m13 = predict(taskswitch_mm_13_rt))
glimpse(taskswitch_mixed_rt)

taskswitch_mm_19_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (1 | participant_id / session),
  dispformula = ~ poly(fitted_m13, 2) + session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_19_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_19_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_19_rt <- cv(taskswitch_mm_19_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_19_rt)
plot(cv_taskswitch_mm_19_rt)

# - glmmTMB log 19b (3-way + (c)variance + 1 covariate + error(1) + RE1/) -----------------------

taskswitch_mm_19b_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + ss_error_prop + 
    (1 | participant_id / session),
  dispformula = ~ poly(fitted_m13, 2) + session + group * task_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_19b_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_19b_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_19_rt <- cv(taskswitch_mm_19_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_19_rt)
plot(cv_taskswitch_mm_19_rt)

# - glmmTMB log 19c (3-way + (ce)variance + 1 covariate + RE1/) -----------------------

taskswitch_mm_19c_rt <- glmmTMB(
  log_rt ~ group * session * task_transition + 
    group * session * congruence + 
    cue_transition + 
    (1 | participant_id / session),
  dispformula = ~ poly(fitted_m13, 2) + session + 
    task_transition * group + cue_transition,
  data = taskswitch_mixed_rt,
  family = gaussian()
)
summary(taskswitch_mm_19c_rt)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_19c_rt)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_19c_rt <- cv(taskswitch_mm_19c_rt, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_19c_rt)
plot(cv_taskswitch_mm_19c_rt)

# MIXED-EFFECTS FOLLOW-UPS: RT ------------------------------------

# Anova table
car::Anova(taskswitch_mm_19c_rt, type = 'III')
car::Anova(taskswitch_mm_10b_rt, type = 'III')
car::Anova(taskswitch_mm_7g_rt, type = 'III')

# - Switch cost RT glmmTMB log 10b -----------------------

switch_cost_emm_rt_10b <- emmeans(taskswitch_mm_10b_rt, 
                                  ~ group * session * task_transition)
switch_cost_emm_rt_10b

# Tibble for plotting
switch_cost_emm_rt_tibble_10b <- switch_cost_emm_rt_10b %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Post-hocs 3-way
switch_cost_3way_rt_10b <- emmeans(taskswitch_mm_10b_rt,
                                   ~ session * task_transition * group) %>%
  contrast(interaction = c('pairwise', 'pairwise', 'pairwise'))
switch_cost_3way_rt_10b

# Post-hocs 2-way
switch_cost_2way_rt_10b <- emmeans(taskswitch_mm_10b_rt,
                                   ~ task_transition * session | group) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
switch_cost_2way_rt_10b

# Joint test of 2-way contrasts
switch_cost_2way_rt_10b_joint <- test(switch_cost_2way_rt_10b, by = NULL, joint = TRUE)
switch_cost_2way_rt_10b_joint

# Session given trial type
switch_cost_session_rt_10b <- contrast(switch_cost_emm_rt_10b, 'pairwise', 
                                   by = c('group', 'task_transition'), 
                                   combine = TRUE)
switch_cost_session_rt_10b

# Trial type given session
switch_cost_trial_rt_10b <- contrast(switch_cost_emm_rt_10b, 'pairwise', 
                                 by = c('group', 'session'), 
                                 combine = TRUE)
switch_cost_trial_rt_10b

# Baseline comparison
switch_cost_baseline_rt <- emmeans(taskswitch_mm_10b_rt,
                                   ~ group * task_transition | session) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
switch_cost_baseline_rt

# Within-subject SE plot #
# ----------------------------- #

data_mxmdl_rt_switch_10b <- afex_plot(
  taskswitch_mm_10b_rt,
  x = 'session',
  trace = 'task_transition',      
  panel = 'group',                
  id = 'participant_id',
  dv = 'log_rt',
  data = taskswitch_mixed_rt,
  within_vars = c('task_transition', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_mxmdl_rt_switch_adj_10b <- data_mxmdl_rt_switch_10b$data %>%
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

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_switch_rt <- ggplot(switch_cost_emm_rt_tibble_10b,
       aes(x = session, y = emmean,
           colour = group,
           linetype = task_transition,
           group = interaction(group, task_transition))) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(aes(ymin = emmean - data_mxmdl_rt_switch_10b$means$SE,
                    ymax = emmean + data_mxmdl_rt_switch_10b$means$SE),
                width = .1, position = pos_dodge) +
  geom_jitter(
    data = data_mxmdl_rt_switch_adj_10b,
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
  labs(x = NULL, y = 'Reaction Time (log ms)') +
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

# - Incongruence cost RT glmmTMB log 10b ------------------------

incongruence_cost_emm_10b <- emmeans(taskswitch_mm_10b_rt, 
                                     ~ group * session * congruence)
incongruence_cost_emm_10b

# Tibble for plotting
incongruence_cost_emm_tibble_10b <- incongruence_cost_emm_10b %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'),
                        labels = c('Congruent', 'Incongruent')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Post-hocs 3-way
incongruence_cost_int_10b <- emmeans(taskswitch_mm_10b_rt, 
                                     ~ session * congruence * group) %>%
  contrast(interaction = c('pairwise', 'pairwise', 'pairwise'))
incongruence_cost_int_10b

# Post-hocs 2-way
incongruence_cost_2way_10b <- emmeans(taskswitch_mm_10b_rt, 
                                      ~ congruence * session | group) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
incongruence_cost_2way_10b

# Joint test of 2-way contrasts
switch_cost_2way_10b_joint <- test(incongruence_cost_2way_10b, 
                                   by = NULL, joint = TRUE)
switch_cost_2way_10b_joint

# Session given trial type
incongruence_cost_session_10b <- contrast(incongruence_cost_emm_10b, 'pairwise', 
                                      by = c('group', 'congruence'), 
                                      combine = TRUE)
incongruence_cost_session_10b

# Trial type given session
incongruence_cost_trial_10b <- contrast(incongruence_cost_emm_10b, 'pairwise', 
                                    by = c('group', 'session'), 
                                    combine = TRUE)
incongruence_cost_trial_10b

# Baseline comparison
incongruence_cost_baseline_rt <- emmeans(taskswitch_mm_10b_rt,
                                   ~ group * congruence | session) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
incongruence_cost_baseline_rt

# Within-subject SE plot #
# ----------------------------- #

data_mxmdl_rt_cngr_10b <- afex_plot(
  taskswitch_mm_10b_rt,
  x = 'session',
  trace = 'congruence',      
  panel = 'group',                
  id = 'participant_id',
  dv = 'log_rt',
  data = taskswitch_mixed_rt,
  within_vars = c('congruence', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_mxmdl_rt_cngr_adj_10b <- data_mxmdl_rt_cngr_10b$data %>%
  mutate(
    session = factor(
      session,
      levels = c(1,2),
      labels = c('Pretest', 'Posttest')),
    congruence = factor(
      congruence,
      levels = c('congruent','incongruent'),
      labels = c('Congruent', 'Incongruent')
    ))

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_incongr_rt <- ggplot(incongruence_cost_emm_tibble_10b,
       aes(x = session, y = emmean,
           colour = group,
           linetype = congruence,
           group = interaction(group, congruence))) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(aes(ymin = emmean - data_mxmdl_rt_cngr_10b$means$SE,
                    ymax = emmean + data_mxmdl_rt_cngr_10b$means$SE),
                width = .1, position = pos_dodge) +
  geom_jitter(
    data = data_mxmdl_rt_cngr_adj_10b,
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
  labs(x = NULL, y = 'Reaction Time (log ms)') +
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

# - Switch cost RT glmmTMB Gamma log 7g -----------------

switch_cost_emm_rt_7g <- emmeans(taskswitch_mm_7g_rt, 
                                 ~ group * session * task_transition,
                                 type = 'response')
switch_cost_emm_rt_7g

# Tibble for plotting
switch_cost_emm_rt_tibble_7g <- switch_cost_emm_rt_7g %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Post-hocs 3-way
switch_cost_3way_rt_7g <- emmeans(taskswitch_mm_7g_rt,
                                   ~ session * task_transition * group) %>%
  contrast(interaction = c('pairwise', 'pairwise', 'pairwise'), type = 'response')
switch_cost_3way_rt_7g

# Post-hocs 2-way
switch_cost_2way_rt_7g <- emmeans(taskswitch_mm_7g_rt,
                                   ~ task_transition * session | group) %>%
  contrast(interaction = c('pairwise', 'pairwise'),
           type = 'response')
switch_cost_2way_rt_7g

# Joint test of 2-way contrasts
switch_cost_2way_rt_7g_joint <- test(switch_cost_2way_rt_7g, by = NULL, joint = TRUE)
switch_cost_2way_rt_7g_joint

# Session given trial type
switch_cost_session_rt_7g <- contrast(switch_cost_emm_rt_7g, 'pairwise', 
                                       by = c('group', 'task_transition'), 
                                       combine = TRUE, type = 'response')
switch_cost_session_rt_7g

# Trial type given session
switch_cost_trial_rt_7g <- contrast(switch_cost_emm_rt_7g, 'pairwise', 
                                     by = c('group', 'session'), 
                                     combine = TRUE, type = 'response')
switch_cost_trial_rt_7g

# Baseline comparison
switch_cost_baseline_rt_7g <- emmeans(taskswitch_mm_7g_rt,
                                   ~ group * task_transition | session) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
switch_cost_baseline_rt_7g

# Within-subject SE plot #
# ----------------------------- #

data_mxmdl_rt_switch_7g <- afex_plot(
  taskswitch_mm_7g_rt,
  x = 'session',
  trace = 'task_transition',      
  panel = 'group',                
  id = 'participant_id',
  dv = 'log_rt',
  data = taskswitch_mixed_rt,
  within_vars = c('task_transition', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_mxmdl_rt_switch_adj_7g <- data_mxmdl_rt_switch_7g$data %>%
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

# Plot
dodge_w <- 0.2
pos_dodge <- position_dodge(width = dodge_w)

ggplot(switch_cost_emm_rt_tibble_7g,
       aes(x = session, y = response,
           colour = group,
           linetype = task_transition,
           group = interaction(group, task_transition))) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(aes(ymin = response - data_mxmdl_rt_switch_7g$means$SE,
                    ymax = response + data_mxmdl_rt_switch_7g$means$SE),
                width = .1, position = pos_dodge) +
  geom_jitter(
    data = data_mxmdl_rt_switch_adj_7g,
    aes(x = session, y = y, colour = group,
        group = interaction(group, task_transition)),
    position = position_jitterdodge(
      jitter.width = 0.05, jitter.height = 0, dodge.width = dodge_w
    ),
    alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal, name = 'Diet group') +
  scale_linetype_manual(values = c(Repeat = 'solid', Switch = 'dashed'),
                        name = 'Trial type') +
  labs(x = NULL, y = 'Estimated marginal mean RT (log ms)',
       title = 'Switch cost by diet, session and task transition') +
  theme_minimal(base_size = 12)

# MIXED-EFFECTS MODELS ER (ACC) --------------------------
# - glmmTMB 1 (4-way + all task covariates + RE1/) --------------------

taskswitch_mm_1_er <- glmmTMB(
  response_correct ~ 
    group * session * task_transition * congruence + 
    cue_transition + prop_switch + prop_incong + prop_cue_switch + 
    (1 | participant_id / session), 
  data = taskswitch_mixed_er, family = binomial())
summary(taskswitch_mm_1_er)

# - glmmTMB 2 (4-way + all task covariates + RT(3) + RE1/) --------------------

taskswitch_mm_2_er <- glmmTMB(
  response_correct ~ 
    group * session * task_transition * congruence + 
    cue_transition + prop_switch + prop_incong + 
    prop_cue_switch + bs_log_rt + ss_log_rt + ws_log_rt + 
    (1 | participant_id / session), 
  data = taskswitch_mixed_er, 
  family = binomial())
summary(taskswitch_mm_2_er)

# - glmmTMB 3 (3-way + group * switch * congruence + all task covariates + RE1/) -----------------------

taskswitch_mm_3_er <- glmmTMB(
  response_correct ~ group * session * task_transition + 
    group * session * congruence + 
    group * task_transition * congruence + 
    cue_transition + 
    prop_switch + prop_incong + prop_cue_switch + 
    (1 | participant_id / session), 
  data = taskswitch_mixed_er, 
  family = binomial())
summary(taskswitch_mm_3_er)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_1_er, taskswitch_mm_3_er)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_3_er)

# - glmmTMB 4 (3-way + group * switch * congruence + some task covariates + RE1/) -----------------------

taskswitch_mm_4_er <- glmmTMB(
  response_correct ~ group * session * task_transition + 
    group * session * congruence + 
    group * task_transition * congruence + 
    cue_transition +
    prop_switch + prop_cue_switch + (1 | participant_id / session), 
  data = taskswitch_mixed_er, family = binomial())
summary(taskswitch_mm_4_er)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_4_er)

# - glmmTMB 5 (3-way + group * switch * congruence + RE1/) -----------------------

taskswitch_mm_5_er <- glmmTMB(
  response_correct ~ 
    group * session * task_transition + 
    group * session * congruence + 
    group * task_transition * congruence + 
    cue_transition +
    (1 | participant_id / session), 
  data = taskswitch_mixed_er, 
  family = binomial())
summary(taskswitch_mm_5_er)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_5_er, taskswitch_mm_4_er, taskswitch_mm_3_er)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_5_er)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_5_er <- cv(taskswitch_mm_5_er, 
                             k = 5, 
                             clusterVariables = 'participant_id', 
                             ncores = n_cores,
                             reps = 10)
summary(cv_taskswitch_mm_5_er)
plot(cv_taskswitch_mm_5_er)

# - glmmTMB 5b (3-way + group * switch * congruence + RE1+) -----------------------

taskswitch_mm_5b_er <- glmmTMB(
  response_correct ~ group * session * task_transition + 
    group * session * congruence + 
    group * task_transition * congruence + 
    cue_transition +
    (task_transition + congruence + cue_transition || participant_id), 
  data = taskswitch_mixed_er, 
  family = binomial())
summary(taskswitch_mm_5b_er)

# Model comparisons #
# ========================== #

anova(taskswitch_mm_5_er, taskswitch_mm_5b_er)

# Diagnostics #
# ========================= #

diagnose_model(taskswitch_mm_5b_er)

# Cross validation #
# ========================= #

# Cluster-based
cv_taskswitch_mm_5b_er <- cv(taskswitch_mm_5b_er, 
                            k = 5, 
                            clusterVariables = 'participant_id', 
                            ncores = n_cores,
                            reps = 10)
summary(cv_taskswitch_mm_5b_er)
plot(cv_taskswitch_mm_5b_er)

# MIXED-EFFECTS FOLLOW-UPS: ER ------------------------------------

# Anova table
car::Anova(taskswitch_mm_5_er, type = 'III')

# - Switch cost ER -----------------------

# NOTE: task transition and congruence interacted, so better to analyse per 
# congruence type

switch_cost_emm_er <- emmeans(taskswitch_mm_5_er, 
                              ~ group * session * task_transition | congruence,
                              type = 'response')
switch_cost_emm_er

# Tibble for plotting
switch_cost_emm_er_tibble <- switch_cost_emm_er %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'), 
                        labels = c('Congruent', 'Incongruent')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Baseline comparison
switch_cost_baseline_er <- emmeans(taskswitch_mm_5_er,
                                   ~ group * task_transition | session) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
switch_cost_baseline_er

# Within-subject ER plot #
# ------------------------------ #

data_mxmdl_er_switch <- afex_plot(
  taskswitch_mm_5_er,
  x = 'session',
  trace = 'task_transition',
  panel = c('group', 'congruence'),
  id = 'participant_id',
  dv = 'response_correct',
  data = taskswitch_mixed_er,
  within_vars = c('congruence','task_transition', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_mxmdl_er_switch_adj <- data_mxmdl_er_switch$data %>%
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'),
                        labels = c('Congruent', 'Incongruent')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_switch_acc <- ggplot(
  switch_cost_emm_er_tibble,
  aes(
    x = session, y = prob * 100,
    colour = group,
    linetype = task_transition,
    group = interaction(group, task_transition)
  )
) +
  geom_jitter(
    data = data_mxmdl_er_switch_adj,
    aes(x = session, y = y * 100, colour = group,
        group = interaction(group, task_transition)),
    position = position_jitterdodge(
      jitter.width = 0.05, jitter.height = 0, dodge.width = dodge_tsmm
    ),
    alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(
    aes(
      ymin = prob * 100 - data_mxmdl_er_switch$means$SE * 100,
      ymax = prob * 100 + data_mxmdl_er_switch$means$SE * 100
    ),
    width = .1, position = pos_dodge
  ) +
  facet_wrap( ~ congruence + group, 
              labeller = labeller(group = c('CD' = 'Clean Diet', 
                                            'KD' = 'Ketogenic Diet'),
                                  congruence = c('Congruent', 
                                                 'Incongruent'))) +
  scale_colour_manual(values = pal, guide = 'none') +
  scale_linetype_manual(
    values = c(Repeat = 'solid', Switch = 'dashed'),
    name = 'Trial type'
  ) +
  labs(
    x = NULL, y = 'Mean Accuracy (%)'
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'mixed_switch_acc.pdf'),
  plot = gg_mixed_switch_acc,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Incongruence cost ER ------------------------

# NOTE: task transition and congruence interacted, so better to analyse per 
# task transition type

incongruence_cost_emm_er <- emmeans(taskswitch_mm_5_er, ~ group * session * congruence | task_transition, type = 'response')
incongruence_cost_emm_er

# Tibble for plotting
incongruence_cost_emm_er_tibble <- incongruence_cost_emm_er %>% 
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'),
                        labels = c('Congruent', 'Incongruent')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Baseline comparison
incongruence_cost_baseline_er <- emmeans(taskswitch_mm_5_er,
                                   ~ group * congruence | session) %>%
  contrast(interaction = c('pairwise', 'pairwise'))
incongruence_cost_baseline_er

# Within-subject ER plot #
# ------------------------------ #

data_mxmdl_er_incongruence <- afex_plot(
  taskswitch_mm_5_er,
  x = 'session',
  trace = 'congruence',
  panel = c('group', 'task_transition'),
  id = 'participant_id',
  dv = 'response_correct',
  data = taskswitch_mixed_er,
  within_vars = c('congruence','task_transition', 'session'),
  between_vars = 'group',
  error = 'within',
  dodge = .15,
  point_arg = list(size = 3),
  line_arg  = list(linewidth = .8),
  return = 'data'
)

# Adjust afex data output
data_mxmdl_er_incongruence_adj <- data_mxmdl_er_incongruence$data %>%
  as_tibble() %>%
  mutate(
    session = factor(session, levels = c(1, 2), labels = c('Pretest', 'Posttest')),
    congruence = factor(congruence, levels = c('congruent', 'incongruent'),
                        labels = c('Congruent', 'Incongruent')),
    task_transition = factor(task_transition, levels = c('repeat', 'switch'),
                             labels = c('Repeat', 'Switch')),
    group = factor(group, levels = c('CD', 'KD'))
  )

# Plot
pos_dodge <- position_dodge(width = dodge_tsmm)
gg_mixed_incongr_acc <- ggplot(
  incongruence_cost_emm_er_tibble,
  aes(
    x = session, y = prob * 100,
    colour = group,
    linetype = congruence,
    group = interaction(group, congruence)
  )
) +
  geom_jitter(
    data = data_mxmdl_er_incongruence_adj,
    aes(x = session, y = y * 100, colour = group,
        group = interaction(group, congruence)),
    position = position_jitterdodge(
      jitter.width = 0.05, jitter.height = 0, dodge.width = dodge_tsmm
    ),
    alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_line(position = pos_dodge, linewidth = .8) +
  geom_point(position = pos_dodge, size = 3) +
  geom_errorbar(
    aes(
      ymin = prob * 100 - data_mxmdl_er_incongruence$means$SE * 100,
      ymax = prob * 100 + data_mxmdl_er_incongruence$means$SE * 100
    ),
    width = .1, position = pos_dodge
  ) +
  facet_wrap(~ task_transition + group, 
              labeller = labeller(group = c('CD' = 'Clean Diet', 
                                            'KD' = 'Ketogenic Diet'),
                                  task_transition = c('Repeat', 
                                                 'Switch'))) +
  scale_colour_manual(values = pal, guide = 'none') +
  scale_linetype_manual(
    values = c(Congruent = 'solid', Incongruent = 'dashed'),
    name = 'Trial Type'
  ) +
  labs(
    x = NULL, y = 'Mean Accuracy (%)'
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'mixed_incongr_acc.pdf'),
  plot = gg_mixed_incongr_acc,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# MIXED-EFFECTS BASELINE --------------------------

# RT #
# ========= #

# Grid
taskswitch_emm_base_rt <- emmeans(taskswitch_mm_10b_rt, ~ group | session)
taskswitch_emm_base_rt

# Comparison
taskswitch_base_rt <- contrast(taskswitch_emm_base_rt, 'pairwise')
taskswitch_base_rt

# ACC #
# ========= #

# Grid
taskswitch_emm_base_er <- emmeans(taskswitch_mm_5_er, ~ group | session)
taskswitch_emm_base_er

# Comparison
taskswitch_base_er <- contrast(taskswitch_emm_base_er, 'pairwise')
taskswitch_base_er

# TASK-SWITCHING AGGREGATION RE-ANALYSIS -------------------------

glimpse(taskswitch_mixed_rt)
glimpse(taskswitch_mixed_er)

# - Aggregate data --------------

# Assess subject-level RT distributions
inspect_subject_dist(
  taskswitch_mixed_rt %>% filter(cue_transition == 'switch'),
  n = 10,
  pid_col = 'participant_id',
  col_dist = 'response_time',
  group_col = 'group',
  session_col = 'session',
  wrap = c('task_transition', 'congruence'),
  pal = pal,
  robust = 'MAD',
  mad_mult = 2.2,
  seed = NULL
)

# Aggregate task-switching RT data
ts_agg_rt <- taskswitch_mixed_rt %>%
  filter(cue_transition == 'switch') %>%
  group_by(participant_id, group, session, task_transition, congruence) %>%
  filter({
    m <- median(response_time, na.rm = TRUE)
    d <- mad(response_time, na.rm = TRUE)
    between(response_time, m - d * 2.2, m + d * 2.2)
  }) %>%
  summarise(subject_mean_rt = mean(response_time, na.rm = TRUE),
            .groups = 'drop') %>%
  pivot_wider(names_from = c('task_transition', 'congruence'),
              values_from = 'subject_mean_rt') %>%
  rowwise() %>%
  mutate(
    switch_cost = mean(c_across(c('switch_congruent','switch_incongruent')),
                       na.rm = TRUE) -
      mean(c_across(c('repeat_congruent','repeat_incongruent')), na.rm = TRUE),
    incongruence_cost = mean(c_across(c('repeat_incongruent','switch_incongruent')),
                             na.rm = TRUE) -
      mean(c_across(c('repeat_congruent','switch_congruent')), na.rm = TRUE),
    average_rt = mean(c_across(c('repeat_congruent', 'repeat_incongruent', 
                                    'switch_congruent', 'switch_incongruent')),
                           na.rm = TRUE)
  ) %>%
  ungroup()
glimpse(ts_agg_rt)

# Aggregate task-switching ER data (in terms of ACC)
ts_agg_er <- taskswitch_mixed_er %>%
  filter(cue_transition == 'switch') %>%
  group_by(participant_id, group, session, task_transition, congruence) %>%
  summarise(subject_mean_acc = mean(response_correct, na.rm = TRUE),
            .groups = 'drop') %>%
  pivot_wider(names_from = c('task_transition', 'congruence'),
              values_from = 'subject_mean_acc') %>%
  rowwise() %>%
  mutate(
    # Costs on ACC
    switch_cost = mean(c_across(c('repeat_congruent', 'repeat_incongruent')), 
                       na.rm = TRUE) - mean(c_across(c('switch_congruent',
                                                       'switch_incongruent')),
                                            na.rm = TRUE),
    incongruence_cost = mean(c_across(c('repeat_congruent', 'switch_congruent')),
                             na.rm = TRUE) - mean(c_across(c('repeat_incongruent',
                                                             'switch_incongruent')),
                                                  na.rm = TRUE),
    average_acc = mean(c_across(c('repeat_congruent', 'repeat_incongruent', 
                                  'switch_congruent', 'switch_incongruent')),
                       na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(c('switch_cost', 'incongruence_cost', 'average_acc'),
                ~ .x * 100)) # Turn values into percentage points
glimpse(ts_agg_er)

# Robust descriptives RT
taskswitch_rt_trwin <- trim_winsorise(ts_agg_rt, 
                                      dv_colnames = c(
                                        'average_rt', 
                                        'switch_cost', 
                                        'incongruence_cost'),
                                      within = 'session',
                                      between = 'group',
                                      id = 'participant_id',
                                      tr = 0.2)
print(taskswitch_rt_trwin, width = Inf)

# Robust descriptives ER (ACC)
taskswitch_er_trwin <- trim_winsorise(ts_agg_er, 
                                      dv_colnames = c(
                                        'average_acc',
                                        'switch_cost', 
                                        'incongruence_cost'),
                                      within = 'session',
                                      between = 'group',
                                      id = 'participant_id',
                                      tr = 0.2)
print(taskswitch_er_trwin, width = Inf)

# Check trial-level exclusions #
# ---------------------------------- #

outlying_trials <- taskswitch_mixed_rt %>%
  filter(cue_transition == 'switch') %>%
  group_by(participant_id, group, session, task_transition, congruence) %>%
  mutate(outlier = {
    m <- median(response_time, na.rm = TRUE)
    d <- mad(response_time, na.rm = TRUE)
    !between(response_time, m - d * 2.2, m + d * 2.2)
  }) %>%
  ungroup() %>%
  summarise(outliers_prop = mean(outlier, na.rm = TRUE)) %>%
  pull(outliers_prop)
outlying_trials

# - Visualise average RT -----------------------

# Histogram
ggplot(ts_agg_rt, aes(x = average_rt, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Average Reaction Time',
       y = 'Density',
       fill = 'Group',
       title = 'Average Reaction Time Distributions') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
ggplot(ts_agg_rt,
       aes(x = session, y = average_rt, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, labeller = labeller(group = c('CD' = 'Clean Diet', 
                                                    'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = taskswitch_rt_trwin %>% 
                  filter(dv == 'average_rt'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  guides(fill = 'none', colour = 'none') +
  labs(x = NULL, y = 'Average Reaction Time (ms)', fill = 'Group', 
       colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_avg_rt_agg <- ggplot(ts_agg_rt,
       aes(x = session,
           y = average_rt,
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
  labs(x = 'Session', y = 'Average Reaction Time (ms)', fill = 'Group', 
       colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'avg_rt_agg.pdf'),
  plot = gg_avg_rt_agg,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Visualise switch cost RT -----------------------

# Density
ggplot(ts_agg_rt, aes(x = switch_cost, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Switch Cost Reaction Time (ms)',
       y = 'Density',
       fill = 'Group',
       title = 'Switch Cost Reaction Time Distributions') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_switch_rt_agg_with <- ggplot(ts_agg_rt,
       aes(x = session, y = switch_cost, fill = group, colour = group)) +
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
  geom_errorbar(data = taskswitch_rt_trwin %>% 
                  filter(dv == 'switch_cost'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'Switch Cost Reaction Time (ms)', 
       fill = 'Group', colour = 'Group',) +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switch_cost_rt_agg_with.pdf'),
  plot = gg_switch_rt_agg_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_switch_rt_agg_without <- ggplot(ts_agg_rt,
       aes(x = session,
           y = switch_cost,
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
  labs(x = 'Session', y = 'Switch Cost Reaction Time (ms)', 
       fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switch_cost_rt_agg_without.pdf'),
  plot = gg_switch_rt_agg_without,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# ROBUSTIFICATION #
# ===================== #

# Bar plot
# - Error bars: between-subject SE of trimmed mean
pd <- position_dodge(width = 0.6) 
gg_switchrt_trbar <- ggplot(taskswitch_rt_trwin %>% filter(dv == 'switch_cost'),
       aes(x = session,
           y = mean_tr,
           fill = group,
           colour = group,
           group = group)) +
  geom_col(position = pd, width = 0.6, colour = NA, alpha = 1) + 
  geom_errorbar(aes(ymin = mean_tr - se_tr,
                    ymax = mean_tr + se_tr),
                position = pd, width = 0.15, linewidth = 0.6,
                alpha = 1) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Session',
       y = 'Trimmed Mean Switch Cost (ms)',
       fill = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switchrt_trbar.pdf'),
  plot = gg_switchrt_trbar,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# CHANGE SCORES #
# ===================== #

# Make the dataset wide
ts_agg_rt_wide <- ts_agg_rt %>%
  dplyr::select(participant_id, group, session, switch_cost, incongruence_cost) %>%
  pivot_wider(names_from = 'session',
              values_from = c('switch_cost', 'incongruence_cost')) %>%
  mutate(switch_cost_change = switch_cost_1 - switch_cost_2,
         incongruence_cost_change = incongruence_cost_1 - incongruence_cost_2)
glimpse(ts_agg_rt_wide)

# Density
ggplot(ts_agg_rt_wide, aes(x = switch_cost_change, 
                           colour = group, 
                           fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Switch Cost Change Scores (ms)',
       y = 'Density',
       fill = 'Group',
       title = 'Switch Cost Change Scores Distributions by Group') +
  theme_minimal()

# Boxplots
# - Error bars: between-subject SE of trimmed mean
gg_switchrt_change <- ggplot(ts_agg_rt_wide, aes(x = group,
                     y = switch_cost_change,
                     colour = group,
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
               shape = 18, size = 4, colour = 'black',
               show.legend = FALSE) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Group',
       y = 'Reaction Time Switch Cost Change Score (ms)',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() +
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switchrt_change.pdf'),
  plot = gg_switchrt_change,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Visualise incongruence cost RT -----------------------

# Density
ggplot(ts_agg_rt, aes(x = incongruence_cost, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Incongruence Cost Reaction Time (ms)',
       y = 'Density',
       fill = 'Group',
       title = 'Incongruence Cost Reaction Time Distributions') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_incongr_cost_rt_agg_with <- ggplot(ts_agg_rt,
       aes(x = session, y = incongruence_cost, fill = group, colour = group)) +
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
  geom_errorbar(data = taskswitch_rt_trwin %>% 
                  filter(dv == 'incongruence_cost'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'Incongruence Cost Reaction Time (ms)', 
       fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongr_cost_rt_agg_with.pdf'),
  plot = gg_incongr_cost_rt_agg_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_incongr_cost_rt_agg <- ggplot(ts_agg_rt,
       aes(x = session,
           y = incongruence_cost,
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
  labs(x = 'Session', y = 'Incongruence Cost Reaction Time (ms)', 
       fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongr_cost_rt_agg.pdf'),
  plot = gg_incongr_cost_rt_agg,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# ROBUSTIFICATION #
# ===================== #

# Bar plot
# - Error bars: between-subject SE of trimmed mean
gg_incongrrt_trbar <- ggplot(taskswitch_rt_trwin %>% filter(dv == 'incongruence_cost'),
                            aes(x = session,
                                y = mean_tr,
                                fill = group,
                                colour = group,
                                group = group)) +
  geom_col(position = pd, width = 0.6, colour = NA, alpha = 1) + 
  geom_errorbar(aes(ymin = mean_tr - se_tr,
                    ymax = mean_tr + se_tr),
                position = pd, width = 0.15, linewidth = 0.6,
                alpha = 1) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Session',
       y = 'Trimmed Mean Incongruence Cost (ms)',
       fill = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'gg_incongrrt_trbar.pdf'),
  plot = gg_incongrrt_trbar,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# CHANGE SCORES #
# ===================== #

# Density
ggplot(ts_agg_rt_wide, aes(x = incongruence_cost_change, 
                           colour = group, 
                           fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Incongruence Cost Change Scores (ms)',
       y = 'Density',
       fill = 'Group',
       title = 'Incongruence Cost Change Scores Distributions by Group') +
  theme_minimal()

# Boxplots
# - Error bars: between-subject SE of trimmed mean
gg_incongrrt_change <- ggplot(ts_agg_rt_wide, aes(x = group,
                                                 y = incongruence_cost_change,
                                                 colour = group,
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
               shape = 18, size = 4, colour = 'black',
               show.legend = FALSE) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Group',
       y = 'Reaction Time Switch Cost Change Score (ms)',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() +
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrrt_change.pdf'),
  plot = gg_incongrrt_change,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Visualise average ER (ACC) -----------------------

# Histogram
ggplot(ts_agg_er, aes(x = average_acc, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Average Accuracy',
       y = 'Density',
       fill = 'Group',
       title = 'Average Accuracy Distributions') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
ggplot(ts_agg_er,
       aes(x = session, y = average_acc, fill = group, colour = group)) +
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
  geom_errorbar(data = taskswitch_er_trwin %>% 
                  filter(dv == 'average_acc'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(y = 'Average Accuracy (%)', fill = 'Group', 
       colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_avg_acc_agg <- ggplot(ts_agg_er,
       aes(x = session,
           y = average_acc,
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
    shape = 21,  # lets colour = outline, fill = interior
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
  labs(x = 'Session', y = 'Average Accuracy (%)', 
       fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'avg_acc_agg.pdf'),
  plot = gg_avg_acc_agg,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Visualise switch cost ER (ACC) -----------------------

# Density
ggplot(ts_agg_er, aes(x = switch_cost, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Switch Cost Accuracy (p.p.)',
       y = 'Density',
       fill = 'Group',
       title = 'Switch Cost Accuracy Distributions') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_switch_cost_acc_agg_with <- ggplot(ts_agg_er,
       aes(x = session, y = switch_cost, fill = group, colour = group)) +
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
  geom_errorbar(data = taskswitch_er_trwin %>% 
                  filter(dv == 'switch_cost'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Session', y = 'Switch Cost Accuracy (p.p.)', fill = 'Group', 
       colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switch_cost_acc_agg_with.pdf'),
  plot = gg_switch_cost_acc_agg_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_switch_cost_acc_agg <- ggplot(ts_agg_er,
       aes(x = session,
           y = switch_cost,
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
  labs(x = 'Session', y = 'Switch Cost Accuracy (p.p.)', 
       fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switch_cost_acc_agg.pdf'),
  plot = gg_switch_cost_acc_agg,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# CHANGE SCORES #
# ===================== #

# Make the dataset wide
ts_agg_er_wide <- ts_agg_er %>%
  dplyr::select(participant_id, group, session, switch_cost, incongruence_cost) %>%
  pivot_wider(names_from = 'session',
              values_from = c('switch_cost', 'incongruence_cost')) %>%
  mutate(switch_cost_change = switch_cost_1 - switch_cost_2,
         incongruence_cost_change = incongruence_cost_1 - incongruence_cost_2)
glimpse(ts_agg_rt_wide)

# Density
ggplot(ts_agg_er_wide, aes(x = switch_cost_change, 
                           colour = group, 
                           fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Switch Cost Change Scores (Accuracy)',
       y = 'Density',
       fill = 'Group',
       title = 'Switch Cost Change Scores Distributions by Group') +
  theme_minimal()

# Boxplots
# - Error bars: between-subject SE of trimmed mean
gg_switcher_change <- ggplot(ts_agg_er_wide, aes(x = group,
                                                  y = switch_cost_change,
                                                  colour = group,
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
               shape = 18, size = 4, colour = 'black',
               show.legend = FALSE) +
  stat_summary(aes(group = group),
               fun.data = mean_se_tr,
               fun.args = list(tr = 0.2),
               geom = 'errorbar',
               colour = 'black',
               width = 0.05,
               position = position_dodge(width = dodge_w)) + 
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = 'Group',
       y = 'Switch Cost Change Score (Accuracy)',
       fill = 'Group') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa() +
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switcher_change.pdf'),
  plot = gg_switcher_change,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Visualise incongruence cost ER (ACC) -----------------------

# Density
ggplot(ts_agg_er, aes(x = incongruence_cost, colour = group, fill = group)) +
  geom_density(alpha = 0.7, adjust = 1) +
  facet_wrap(~ group + session) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(colour = 'none') + 
  labs(x = 'Incongruence Cost Accuracy (p.p.)',
       y = 'Density',
       fill = 'Group',
       title = 'Incongruence Cost Accuracy Distributions') +
  theme_minimal()

# With trajectories
# - Error bars: within-subject SE of trimmed mean
gg_incongr_cost_acc_agg_with <- ggplot(ts_agg_er,
       aes(x = session, y = incongruence_cost, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.5, size = 1.8, shape = 21, stroke = 0.4) +
  facet_wrap(~ group, labeller = labeller(group = c('CD' = 'Clean Diet', 
                                                    'KD' = 'Ketogenic Diet'))) +
  geom_line(aes(group = participant_id),
            alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = group),
               fun = function(z) mean(z, trim = 0.2, na.rm = TRUE),
               geom = 'point', shape = 18, size = 3.5, colour = 'black',
               show.legend = FALSE) +
  geom_errorbar(data = taskswitch_er_trwin %>% 
                  filter(dv == 'incongruence_cost'),
                aes(y = mean_tr, 
                    ymin = mean_tr - se_tr_ws, 
                    ymax = mean_tr + se_tr_ws,
                    group = group),
                colour = 'black',
                width = 0.05) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = NULL, y = 'Incongruence Cost Accuracy (p.p.)', 
       fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongr_cost_acc_agg_with.pdf'),
  plot = gg_incongr_cost_acc_agg_with,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Without trajectories
# - Error bars: between-subject SE of trimmed mean
gg_incongr_cost_acc_agg <- ggplot(ts_agg_er,
       aes(x = session,
           y = incongruence_cost,
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
    shape = 21,  # lets colour = outline, fill = interior
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
  labs(x = 'Session', y = 'Incongruence Cost Accuracy (p.p.)', 
       fill = 'Group', colour = 'Group') +
  scale_x_discrete(labels = c('Pretest', 'Posttest')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongr_cost_acc_agg.pdf'),
  plot = gg_incongr_cost_acc_agg,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# - Trimmed average RT ------------------

# Trimmed ANOVA
avg_rt_tranova <- WRS2::bwtrim(average_rt ~ group * session, 
                                  id = participant_id, 
                                  data = ts_agg_rt,
                                  tr = 0.2
)
avg_rt_tranova

# Baseline #
# ----------------- #

# Basic Yuen
avg_rt_yuen <- yuen(average_rt ~ group * session,
                    data = ts_agg_rt %>% 
                      filter(session == 1),
                    tr = 0.2)
avg_rt_yuen

# Bootstrapped
avg_rt_yuen_bt <- yuenbt(average_rt ~ group,
                    data = ts_agg_rt %>% 
                      filter(session == 1),
                    tr = 0.2,
                    nboot = 100000,
                    side = TRUE)
avg_rt_yuen_bt

# Quantiles
avg_rt_quant_baseline <- WRS2::qcomhd(average_rt ~ group, 
                                    data = ts_agg_rt %>% 
                                      filter(session == 1), 
                                    q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                                    nboot = 20000, alpha = 0.05, ADJ.CI = FALSE)
avg_rt_quant_baseline

# Effect sizes #
# -------------------- #

# Create ordered list for analysis
avg_rt_post_tr_list <- build_split_list(ts_agg_rt, 
                                      'average_rt', c('group', 'session'))

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
# - No bootstrap
avg_rt_post_tr_int <- WRS::bwimcpES(2, 2, avg_rt_post_tr_list, tr=0.2, 
                                  CI=TRUE, alpha=0.05)
avg_rt_post_tr_int

# Effect size of session changes
avg_rt_session_eff <- bw.es.B(2, 2, avg_rt_post_tr_list, 
                            tr = 0.2, POOL = TRUE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
avg_rt_session_eff

# - Bayesian average RT -------------------

# BASELINE COMPARISON #
# ======================= # 

# Check that there are no NAs
print(ts_agg_rt %>%
        summarise(across(everything(),
                         ~ sum(is.na(.x)))), width = Inf)

# Baseline
avgrt_baseline_bf <- ttestBF(x = ts_agg_rt %>% 
                             filter(session == 1, 
                                    group == 'CD') %>%
                             pull(average_rt),
                           y = ts_agg_rt %>% 
                             filter(session == 1, 
                                    group == 'KD') %>%
                             pull(average_rt),
                           rscale = 'medium')
avgrt_baseline_bf

# Baseline (r = 0.5)
avgrt_baseline_bf_smaller <- ttestBF(x = ts_agg_rt %>% 
                               filter(session == 1, 
                                      group == 'CD') %>%
                               pull(average_rt),
                             y = ts_agg_rt %>% 
                               filter(session == 1, 
                                      group == 'KD') %>%
                               pull(average_rt),
                             rscale = 0.5)
avgrt_baseline_bf_smaller

# Baseline (r = 0.3)
avgrt_baseline_bf_small <- ttestBF(x = ts_agg_rt %>% 
                                       filter(session == 1, 
                                              group == 'CD') %>%
                                       pull(average_rt),
                                     y = ts_agg_rt %>% 
                                       filter(session == 1, 
                                              group == 'KD') %>%
                                       pull(average_rt),
                                     rscale = 0.3)
avgrt_baseline_bf_small

# BAYESIAN ANOVA #
# ======================= # 

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
avgrt_bfanova_medium <- anovaBF(
  average_rt ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
avgrt_bfanova_medium

# Interaction BF
avgrt_bf_int_medium <- avgrt_bfanova_medium['group + session + group:session + participant_id'] /
  avgrt_bfanova_medium['group + session + participant_id']
avgrt_bf_int_medium

# Session BF
avgrt_bf_sess_medium <- avgrt_bfanova_medium['group + session + group:session + participant_id'] /
  avgrt_bfanova_medium['group + group:session + participant_id']
avgrt_bf_sess_medium

# Group BF
avgrt_bf_group_medium <- avgrt_bfanova_medium['group + session + group:session + participant_id'] /
  avgrt_bfanova_medium['session + group:session + participant_id']
avgrt_bf_group_medium

# Wide prior #
# ---------------- #

avgrt_bfanova_wide <- anovaBF(
  average_rt ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
avgrt_bfanova_wide

# Interaction BF
avgrt_bf_int_wide <- avgrt_bfanova_wide['group + session + group:session + participant_id'] /
  avgrt_bfanova_wide['group + session + participant_id']
avgrt_bf_int_wide

# Session BF
avgrt_bf_sess_wide <- avgrt_bfanova_wide['group + session + group:session + participant_id'] /
  avgrt_bfanova_wide['group + group:session + participant_id']
avgrt_bf_sess_wide

# Group BF
avgrt_bf_group_wide <- avgrt_bfanova_wide['group + session + group:session + participant_id'] /
  avgrt_bfanova_wide['session + group:session + participant_id']
avgrt_bf_group_wide

# Small prior #
# ---------------- #

avgrt_bfanova_small <- anovaBF(
  average_rt ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 100000
)
avgrt_bfanova_small

# Interaction BF
avgrt_bf_int_small <- avgrt_bfanova_small['group + session + group:session + participant_id'] /
  avgrt_bfanova_small['group + session + participant_id']
avgrt_bf_int_small

# Session BF
avgrt_bf_sess_small <- avgrt_bfanova_small['group + session + group:session + participant_id'] /
  avgrt_bfanova_small['group + group:session + participant_id']
avgrt_bf_sess_small

# Group BF
avgrt_bf_group_small <- avgrt_bfanova_small['group + session + group:session + participant_id'] /
  avgrt_bfanova_small['session + group:session + participant_id']
avgrt_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
avgrt_bfratio_int <- bf_ratio_plot(
  Small = avgrt_bf_int_small,
  Medium = avgrt_bf_int_medium,
  Wide = avgrt_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
avgrt_bfratio_int

# Session
avgrt_bfratio_sess <- bf_ratio_plot(
  Small = avgrt_bf_sess_small,
  Medium = avgrt_bf_sess_medium,
  Wide = avgrt_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
avgrt_bfratio_sess

# Group
avgrt_bfratio_group <- bf_ratio_plot(
  Small = avgrt_bf_group_small,
  Medium = avgrt_bf_group_medium,
  Wide = avgrt_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
avgrt_bfratio_group

# Medium prior INT ~ posterior #
# ---------------------------- #

avgrt_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = ts_agg_rt,
  dv = 'average_rt',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD'
)

avgrt_bfanova_medium_posterior_int$summary
avgrt_bfanova_medium_posterior_int$plot

# Medium prior SESSION ~ posterior #
# ---------------------------- #

avgrt_bfanova_medium_posterior_session <- bf_contrast_posterior(
  ts_agg_rt, 'average_rt', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen', binwidth = 2
)

avgrt_bfanova_medium_posterior_session$summary
avgrt_bfanova_medium_posterior_session$plot

# Medium prior GROUP ~ posterior #
# ---------------------------- #

avgrt_bfanova_medium_posterior_group <- bf_contrast_posterior(
  ts_agg_rt, 'average_rt', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'coral', binwidth = 5
)

avgrt_bfanova_medium_posterior_group$summary
avgrt_bfanova_medium_posterior_group$plot

# - Trimmed switch cost RT --------------------

# Trimmed ANOVA
switch_rt_tranova <- WRS2::bwtrim(switch_cost ~ group * session, 
                             id = participant_id, 
                             data = ts_agg_rt,
                             tr = 0.2)
switch_rt_tranova

# Baseline comparison #
# --------------------------- #

# Yuen's t-test
switch_rt_baseline_test <- yuen(switch_cost ~ group, 
                           data = ts_agg_rt %>% filter(session == 1),
                           tr = 0.2)
switch_rt_baseline_test

# Bootstrapped
switch_rt_baseline_bt <- yuenbt(switch_cost ~ group, 
                                data = ts_agg_rt %>% filter(session == 1),
                                tr = 0.2,
                                nboot = 10000,
                                side = TRUE)
switch_rt_baseline_bt

# Quantiles
switch_rt_baseline_quant <- qcomhd(switch_cost ~ group, 
                                   data = ts_agg_rt %>% filter(session == 1),
                                   q = c(0.1, 0.25, 0.5, 0.75, 0.9),
       nboot = 10000, alpha = 0.05, ADJ.CI = FALSE)
switch_rt_baseline_quant

# Post-hocs / exploration #
# ------------------------------ #

# Create ordered list for analysis
switch_rt_post_tr_list <- build_split_list(ts_agg_rt, 
                                      'switch_cost', c('group', 'session'))

# Contrasts on trimmed means
switch_rt_post_tr <- bwmcp(2, 2, switch_rt_post_tr_list, tr=0.2,con=0, 
                           nboot=10000)
switch_rt_post_tr

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
switch_rt_post_tr_int <- bwimcpES(2, 2, switch_rt_post_tr_list, 
                                  tr=0.2, CI=TRUE, alpha=0.05)
switch_rt_post_tr_int

# Bootstrapped interaction contrast (DiD)
switch_rt_post_tr_int_bt <- yuenbt(switch_cost_change ~ group,
                                   data = ts_agg_rt_wide, 
                                   tr = 0.2,
                                   nboot = 10000,
                                   side = TRUE)
switch_rt_post_tr_int_bt

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
switch_session_eff <- bw.es.B(2, 2, switch_rt_post_tr_list, 
                            tr = 0.2, POOL = TRUE, OPT = FALSE, 
                            CI = TRUE, SEED = TRUE, REL.MAG = NULL)
switch_session_eff

# - Regular switch cost RT ----------------

switch_rt_anova <- aov_4(switch_cost ~ group * session + (session | participant_id),
                    data = ts_agg_rt,
                    anova_table = list(es = 'pes'))
switch_rt_anova

# - Ranked switch cost RT -------------------

# Perform aligned-rank transform (ART)
switch_rt_ranked <- art(switch_cost ~ group * session + (1 | participant_id),
                   data = ts_agg_rt)

# Check the result of the transformation 
summary(switch_rt_ranked)
anova(switch_rt_ranked, response = 'aligned')

# Run ANOVA on the ART data
switch_rt_ranked_anova <- anova(switch_rt_ranked, response = 'art', 
                           type = 'III', test = 'F')
switch_rt_ranked_anova

# - Bayesian switch cost RT --------------------

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
switchrt_bfanova_medium <- anovaBF(
  switch_cost ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
switchrt_bfanova_medium

# Interaction BF
switchrt_bf_int_medium <- switchrt_bfanova_medium['group + session + group:session + participant_id'] /
  switchrt_bfanova_medium['group + session + participant_id']
switchrt_bf_int_medium

# Session BF
switchrt_bf_sess_medium <- switchrt_bfanova_medium['group + session + group:session + participant_id'] /
  switchrt_bfanova_medium['group + group:session + participant_id']
switchrt_bf_sess_medium

# Group BF
switchrt_bf_group_medium <- switchrt_bfanova_medium['group + session + group:session + participant_id'] /
  switchrt_bfanova_medium['session + group:session + participant_id']
switchrt_bf_group_medium

# Wide prior #
# ---------------- #

switchrt_bfanova_wide <- anovaBF(
  switch_cost ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
switchrt_bfanova_wide

# Interaction BF
switchrt_bf_int_wide <- switchrt_bfanova_wide['group + session + group:session + participant_id'] /
  switchrt_bfanova_wide['group + session + participant_id']
switchrt_bf_int_wide

# Session BF
switchrt_bf_sess_wide <- switchrt_bfanova_wide['group + session + group:session + participant_id'] /
  switchrt_bfanova_wide['group + group:session + participant_id']
switchrt_bf_sess_wide

# Group BF
switchrt_bf_group_wide <- switchrt_bfanova_wide['group + session + group:session + participant_id'] /
  switchrt_bfanova_wide['session + group:session + participant_id']
switchrt_bf_group_wide

# Small prior #
# ---------------- #

switchrt_bfanova_small <- anovaBF(
  switch_cost ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 1000000
)
switchrt_bfanova_small

# Interaction BF
switchrt_bf_int_small <- switchrt_bfanova_small['group + session + group:session + participant_id'] /
  switchrt_bfanova_small['group + session + participant_id']
switchrt_bf_int_small

# Session BF
switchrt_bf_sess_small <- switchrt_bfanova_small['group + session + group:session + participant_id'] /
  switchrt_bfanova_small['group + group:session + participant_id']
switchrt_bf_sess_small

# Group BF
switchrt_bf_group_small <- switchrt_bfanova_small['group + session + group:session + participant_id'] /
  switchrt_bfanova_small['session + group:session + participant_id']
switchrt_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
switchrt_bfratio_int <- bf_ratio_plot(
  Small = switchrt_bf_int_small,
  Medium = switchrt_bf_int_medium,
  Wide = switchrt_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
switchrt_bfratio_int

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switchrt_bfratio_int.pdf'),
  plot = switchrt_bfratio_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Session
switchrt_bfratio_sess <- bf_ratio_plot(
  Small = switchrt_bf_sess_small,
  Medium = switchrt_bf_sess_medium,
  Wide = switchrt_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
switchrt_bfratio_sess

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switchrt_bfratio_sess.pdf'),
  plot = switchrt_bfratio_sess$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Group
switchrt_bfratio_group <- bf_ratio_plot(
  Small = switchrt_bf_group_small,
  Medium = switchrt_bf_group_medium,
  Wide = switchrt_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
switchrt_bfratio_group

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switchrt_bfratio_group.pdf'),
  plot = switchrt_bfratio_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior INT ~ posterior #
# ---------------------------- #

switchrt_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = ts_agg_rt,
  dv = 'switch_cost',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD',
  binwidth = 4
)

switchrt_bfanova_medium_posterior_int$summary
switchrt_bfanova_medium_posterior_int$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switchrt_bfanova_medium_posterior_int.pdf'),
  plot = switchrt_bfanova_medium_posterior_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

switchrt_bfanova_medium_posterior_session <- bf_contrast_posterior(
  ts_agg_rt, 'switch_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen', binwidth = 2.5
)

switchrt_bfanova_medium_posterior_session$summary
switchrt_bfanova_medium_posterior_session$plot

# Save plot
ggsave(
  filename = file.path(plot_directory,
                       'switchrt_bfanova_medium_posterior_session.pdf'),
  plot = switchrt_bfanova_medium_posterior_session$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior GROUP ~ posterior #
# ---------------------------- #

switchrt_bfanova_medium_posterior_group <- bf_contrast_posterior(
  ts_agg_rt, 'switch_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'orangered', binwidth = 2.5
)

switchrt_bfanova_medium_posterior_group$summary
switchrt_bfanova_medium_posterior_group$plot

# Save plot
ggsave(
  filename = file.path(plot_directory,
                       'switchrt_bfanova_medium_posterior_group.pdf'),
  plot = switchrt_bfanova_medium_posterior_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# BASELINE COMPARISON #
# ======================== #

# Baseline
switchrt_baseline_bf <- ttestBF(x = ts_agg_rt %>% 
                                  filter(session == 1, 
                                         group == 'CD') %>%
                                  pull(switch_cost),
                                y = ts_agg_rt %>% 
                                  filter(session == 1, 
                                         group == 'KD') %>%
                                  pull(switch_cost),
                                rscale = 'medium')
switchrt_baseline_bf

# - Trimmed incongruence cost RT --------------------

# Trimmed ANOVA
incongr_rt_tranova <- WRS2::bwtrim(incongruence_cost ~ group * session, 
                                  id = participant_id, 
                                  data = ts_agg_rt,
                                  tr = 0.2
)
incongr_rt_tranova

# Yuen's t-test
incongr_rt_baseline_test <- yuen(incongruence_cost ~ group, 
                                data = ts_agg_rt %>% filter(session == 1),
                                tr = 0.2)
incongr_rt_baseline_test

# Post-hocs / exploration #
# ------------------------------ #

# Create ordered list for analysis
incongr_rt_post_tr_list <- build_split_list(ts_agg_rt, 
                                           'incongruence_cost',
                                           c('group', 'session'))

# Contrasts on trimmed means
incongr_rt_post_tr <- bwmcp(2, 2, incongr_rt_post_tr_list, tr=0.2, 
                            con=0, nboot=8000)
incongr_rt_post_tr

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
incongr_rt_post_tr_int <- bwimcpES(2, 2, incongr_rt_post_tr_list, 
                                   tr=0.2, CI=TRUE, alpha=0.05)
incongr_rt_post_tr_int

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
incongr_session_eff_rt <- bw.es.B(2, 2, incongr_rt_post_tr_list, 
                              tr = 0.2, POOL = TRUE, OPT = FALSE, 
                              CI = TRUE, SEED = TRUE, REL.MAG = NULL)
incongr_session_eff_rt

# - Regular incongruence cost RT ----------------

incongr_rt_anova <- aov_4(
  incongruence_cost ~ group * session + (session | participant_id),
                        data = ts_agg_rt,
                        anova_table = list(es = 'pes'))
incongr_rt_anova

# - Bayesian incongruence cost RT -------------------

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
incongrrt_bfanova_medium <- anovaBF(
  incongruence_cost ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
incongrrt_bfanova_medium

# Interaction BF
incongrrt_bf_int_medium <- incongrrt_bfanova_medium['group + session + group:session + participant_id'] /
  incongrrt_bfanova_medium['group + session + participant_id']
incongrrt_bf_int_medium

# Session BF
incongrrt_bf_sess_medium <- incongrrt_bfanova_medium['group + session + group:session + participant_id'] /
  incongrrt_bfanova_medium['group + group:session + participant_id']
incongrrt_bf_sess_medium

# Group BF
incongrrt_bf_group_medium <- incongrrt_bfanova_medium['group + session + group:session + participant_id'] /
  incongrrt_bfanova_medium['session + group:session + participant_id']
incongrrt_bf_group_medium

# Wide prior #
# ---------------- #

incongrrt_bfanova_wide <- anovaBF(
  incongruence_cost ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
incongrrt_bfanova_wide

# Interaction BF
incongrrt_bf_int_wide <- incongrrt_bfanova_wide['group + session + group:session + participant_id'] /
  incongrrt_bfanova_wide['group + session + participant_id']
incongrrt_bf_int_wide

# Session BF
incongrrt_bf_sess_wide <- incongrrt_bfanova_wide['group + session + group:session + participant_id'] /
  incongrrt_bfanova_wide['group + group:session + participant_id']
incongrrt_bf_sess_wide

# Group BF
incongrrt_bf_group_wide <- incongrrt_bfanova_wide['group + session + group:session + participant_id'] /
  incongrrt_bfanova_wide['session + group:session + participant_id']
incongrrt_bf_group_wide

# Small prior #
# ---------------- #

incongrrt_bfanova_small <- anovaBF(
  incongruence_cost ~ group * session + participant_id,
  data = ts_agg_rt %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 1000000 # 1,000,000 draws on purpose --> needed for error 
  # management
)
incongrrt_bfanova_small

# Interaction BF
incongrrt_bf_int_small <- incongrrt_bfanova_small['group + session + group:session + participant_id'] /
  incongrrt_bfanova_small['group + session + participant_id']
incongrrt_bf_int_small

# Session BF
incongrrt_bf_sess_small <- incongrrt_bfanova_small['group + session + group:session + participant_id'] /
  incongrrt_bfanova_small['group + group:session + participant_id']
incongrrt_bf_sess_small

# Group BF
incongrrt_bf_group_small <- incongrrt_bfanova_small['group + session + group:session + participant_id'] /
  incongrrt_bfanova_small['session + group:session + participant_id']
incongrrt_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
incongrrt_bfratio_int <- bf_ratio_plot(
  Small = incongrrt_bf_int_small,
  Medium = incongrrt_bf_int_medium,
  Wide = incongrrt_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
incongrrt_bfratio_int

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrrt_bfratio_int.pdf'),
  plot = incongrrt_bfratio_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Session
incongrrt_bfratio_sess <- bf_ratio_plot(
  Small = incongrrt_bf_sess_small,
  Medium = incongrrt_bf_sess_medium,
  Wide = incongrrt_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
incongrrt_bfratio_sess

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrrt_bfratio_sess.pdf'),
  plot = incongrrt_bfratio_sess$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Group
incongrrt_bfratio_group <- bf_ratio_plot(
  Small = incongrrt_bf_group_small,
  Medium = incongrrt_bf_group_medium,
  Wide = incongrrt_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
incongrrt_bfratio_group

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrrt_bfratio_group.pdf'),
  plot = incongrrt_bfratio_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior INT ~ posterior #
# ---------------------------- #

incongrrt_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = ts_agg_rt,
  dv = 'incongruence_cost',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD',
  binwidth = 4
)

incongrrt_bfanova_medium_posterior_int$summary
incongrrt_bfanova_medium_posterior_int$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrrt_bfanova_medium_posterior_int.pdf'),
  plot = incongrrt_bfanova_medium_posterior_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

incongrrt_bfanova_medium_posterior_session <- bf_contrast_posterior(
  ts_agg_rt, 'incongruence_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen', binwidth = 2
)

incongrrt_bfanova_medium_posterior_session$summary
incongrrt_bfanova_medium_posterior_session$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrrt_bfanova_medium_posterior_session.pdf'),
  plot = incongrrt_bfanova_medium_posterior_session$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior GROUP ~ posterior #
# ---------------------------- #

incongrrt_bfanova_medium_posterior_group <- bf_contrast_posterior(
  ts_agg_rt, 'incongruence_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'orangered', binwidth = 2
)

incongrrt_bfanova_medium_posterior_group$summary
incongrrt_bfanova_medium_posterior_group$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrrt_bfanova_medium_posterior_group.pdf'),
  plot = incongrrt_bfanova_medium_posterior_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# BASELINE COMPARISON #
# ======================== #

# Baseline
switchrt_baseline_bf <- ttestBF(x = ts_agg_rt %>% 
                                  filter(session == 1, 
                                         group == 'CD') %>%
                                  pull(switch_cost),
                                y = ts_agg_rt %>% 
                                  filter(session == 1, 
                                         group == 'KD') %>%
                                  pull(switch_cost),
                                rscale = 'medium')
switchrt_baseline_bf

# - Trimmed average ER (ACC) -----------------

# Trimmed ANOVA
avg_acc_tranova <- WRS2::bwtrim(average_acc ~ group * session, 
                               id = participant_id, 
                               data = ts_agg_er,
                               tr = 0.2
)
avg_acc_tranova

# Baseline #
# ----------------- #

# Basic Yuen
avg_acc_yuen <- yuen(average_acc ~ group * session,
                    data = ts_agg_er %>% 
                      filter(session == 1),
                    tr = 0.2)
avg_acc_yuen

# Bootstrapped
avg_acc_yuen_bt <- yuenbt(average_acc ~ group,
                         data = ts_agg_er %>% 
                           filter(session == 1),
                         tr = 0.2,
                         nboot = 100000,
                         side = TRUE)
avg_acc_yuen_bt

# Quantiles
avg_acc_quant_baseline <- WRS2::qcomhd(average_acc ~ group, 
                                      data = ts_agg_er %>% 
                                        filter(session == 1), 
                                      q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                                      nboot = 20000, alpha = 0.05, ADJ.CI = FALSE)
avg_acc_quant_baseline

# Effect sizes #
# -------------------- #

# Create ordered list for analysis
avg_acc_post_tr_list <- build_split_list(ts_agg_er, 
                                        'average_acc', c('group', 'session'))

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
# - No bootstrap
avg_acc_post_tr_int <- WRS::bwimcpES(2, 2, avg_acc_post_tr_list, tr=0.2, 
                                    CI=TRUE, alpha=0.05)
avg_acc_post_tr_int

# Effect size of session changes
avg_acc_session_eff <- bw.es.B(2, 2, avg_acc_post_tr_list, 
                              tr = 0.2, POOL = TRUE, OPT = FALSE, 
                              CI = TRUE, SEED = TRUE, REL.MAG = NULL)
avg_acc_session_eff

# - Bayesian average ER (ACC) -------------------

# BASELINE COMPARISON #
# ======================== #

# Check that there are no NAs
print(ts_agg_er %>%
        summarise(across(everything(),
                         ~ sum(is.na(.x)))), width = Inf)

# Baseline
avger_baseline_bf <- ttestBF(x = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'CD') %>%
                               pull(average_acc),
                             y = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'KD') %>%
                               pull(average_acc),
                             rscale = 'medium')
avger_baseline_bf

# Baseline (r = 0.5)
avger_baseline_bf_smaller <- ttestBF(x = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'CD') %>%
                               pull(average_acc),
                             y = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'KD') %>%
                               pull(average_acc),
                             rscale = 0.5)
avger_baseline_bf_smaller

# Baseline (r = 0.3)
avger_baseline_bf_small <- ttestBF(x = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'CD') %>%
                               pull(average_acc),
                             y = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'KD') %>%
                               pull(average_acc),
                             rscale = 0.3)
avger_baseline_bf_small

# BAYESIAN ANOVA #
# ======================= # 

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
avger_bfanova_medium <- anovaBF(
  average_acc ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
avger_bfanova_medium

# Interaction BF
avger_bf_int_medium <- avger_bfanova_medium['group + session + group:session + participant_id'] /
  avger_bfanova_medium['group + session + participant_id']
avger_bf_int_medium

# Session BF
avger_bf_sess_medium <- avger_bfanova_medium['group + session + group:session + participant_id'] /
  avger_bfanova_medium['group + group:session + participant_id']
avger_bf_sess_medium

# Group BF
avger_bf_group_medium <- avger_bfanova_medium['group + session + group:session + participant_id'] /
  avger_bfanova_medium['session + group:session + participant_id']
avger_bf_group_medium

# Wide prior #
# ---------------- #

avger_bfanova_wide <- anovaBF(
  average_acc ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
avger_bfanova_wide

# Interaction BF
avger_bf_int_wide <- avger_bfanova_wide['group + session + group:session + participant_id'] /
  avger_bfanova_wide['group + session + participant_id']
avger_bf_int_wide

# Session BF
avger_bf_sess_wide <- avger_bfanova_wide['group + session + group:session + participant_id'] /
  avger_bfanova_wide['group + group:session + participant_id']
avger_bf_sess_wide

# Group BF
avger_bf_group_wide <- avger_bfanova_wide['group + session + group:session + participant_id'] /
  avger_bfanova_wide['session + group:session + participant_id']
avger_bf_group_wide

# Small prior #
# ---------------- #

avger_bfanova_small <- anovaBF(
  average_acc ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 100000
)
avger_bfanova_small

# Interaction BF
avger_bf_int_small <- avger_bfanova_small['group + session + group:session + participant_id'] /
  avger_bfanova_small['group + session + participant_id']
avger_bf_int_small

# Session BF
avger_bf_sess_small <- avger_bfanova_small['group + session + group:session + participant_id'] /
  avger_bfanova_small['group + group:session + participant_id']
avger_bf_sess_small

# Group BF
avger_bf_group_small <- avger_bfanova_small['group + session + group:session + participant_id'] /
  avger_bfanova_small['session + group:session + participant_id']
avger_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
avger_bfratio_int <- bf_ratio_plot(
  Small = avger_bf_int_small,
  Medium = avger_bf_int_medium,
  Wide = avger_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
avger_bfratio_int

# Session
avger_bfratio_sess <- bf_ratio_plot(
  Small = avger_bf_sess_small,
  Medium = avger_bf_sess_medium,
  Wide = avger_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
avger_bfratio_sess

# Group
avger_bfratio_group <- bf_ratio_plot(
  Small = avger_bf_group_small,
  Medium = avger_bf_group_medium,
  Wide = avger_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
avger_bfratio_group

# Medium prior INT ~ posterior #
# ---------------------------- #

avger_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = ts_agg_er,
  dv = 'average_acc',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD'
)

avger_bfanova_medium_posterior_int$summary
avger_bfanova_medium_posterior_int$plot

# Medium prior SESSION ~ posterior #
# ---------------------------- #

avger_bfanova_medium_posterior_session <- bf_contrast_posterior(
  ts_agg_er, 'average_acc', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen'
)

avger_bfanova_medium_posterior_session$summary
avger_bfanova_medium_posterior_session$plot

# Medium prior GROUP ~ posterior #
# ---------------------------- #

avger_bfanova_medium_posterior_group <- bf_contrast_posterior(
  ts_agg_er, 'average_acc', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'coral'
)

avger_bfanova_medium_posterior_group$summary
avger_bfanova_medium_posterior_group$plot

# - Trimmed switch cost ER (ACC) --------------------

# Trimmed ANOVA
switch_acc_tranova <- WRS2::bwtrim(switch_cost ~ group * session, 
                                  id = participant_id, 
                                  data = ts_agg_er,
                                  tr = 0.2
)
switch_acc_tranova

# Yuen's t-test
switch_acc_baseline_test <- yuen(switch_cost ~ group, 
                                data = ts_agg_er %>% filter(session == 1),
                                tr = 0.2)
switch_acc_baseline_test

# Post-hocs / exploration #
# ------------------------------ #

# Create ordered list for analysis
switch_acc_post_tr_list <- build_split_list(ts_agg_er, 
                                           'switch_cost', c('group', 'session'))

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
switch_acc_post_tr_int <- bwimcpES(2, 2, switch_acc_post_tr_list, 
                                  tr=0.2, CI=TRUE, alpha=0.05)
switch_acc_post_tr_int

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
switch_session_eff_acc <- bw.es.B(2, 2, switch_acc_post_tr_list, 
                              tr = 0.2, POOL = TRUE, OPT = FALSE, 
                              CI = TRUE, SEED = TRUE, REL.MAG = NULL)
switch_session_eff_acc

# - Regular switch cost ER (ACC) ----------------

switch_acc_anova <- aov_4(switch_cost ~ group * session + (session | participant_id),
                        data = ts_agg_er,
                        anova_table = list(es = 'pes'))
switch_acc_anova

# - Bayesian switch cost ER (ACC) -------------------

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
switcher_bfanova_medium <- anovaBF(
  switch_cost ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
switcher_bfanova_medium

# Interaction BF
switcher_bf_int_medium <- switcher_bfanova_medium['group + session + group:session + participant_id'] /
  switcher_bfanova_medium['group + session + participant_id']
switcher_bf_int_medium

# Session BF
switcher_bf_sess_medium <- switcher_bfanova_medium['group + session + group:session + participant_id'] /
  switcher_bfanova_medium['group + group:session + participant_id']
switcher_bf_sess_medium

# Group BF
switcher_bf_group_medium <- switcher_bfanova_medium['group + session + group:session + participant_id'] /
  switcher_bfanova_medium['session + group:session + participant_id']
switcher_bf_group_medium

# Wide prior #
# ---------------- #

switcher_bfanova_wide <- anovaBF(
  switch_cost ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
switcher_bfanova_wide

# Interaction BF
switcher_bf_int_wide <- switcher_bfanova_wide['group + session + group:session + participant_id'] /
  switcher_bfanova_wide['group + session + participant_id']
switcher_bf_int_wide

# Session BF
switcher_bf_sess_wide <- switcher_bfanova_wide['group + session + group:session + participant_id'] /
  switcher_bfanova_wide['group + group:session + participant_id']
switcher_bf_sess_wide

# Group BF
switcher_bf_group_wide <- switcher_bfanova_wide['group + session + group:session + participant_id'] /
  switcher_bfanova_wide['session + group:session + participant_id']
switcher_bf_group_wide

# Small prior #
# ---------------- #

switcher_bfanova_small <- anovaBF(
  switch_cost ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 100000
)
switcher_bfanova_small

# Interaction BF
switcher_bf_int_small <- switcher_bfanova_small['group + session + group:session + participant_id'] /
  switcher_bfanova_small['group + session + participant_id']
switcher_bf_int_small

# Session BF
switcher_bf_sess_small <- switcher_bfanova_small['group + session + group:session + participant_id'] /
  switcher_bfanova_small['group + group:session + participant_id']
switcher_bf_sess_small

# Group BF
switcher_bf_group_small <- switcher_bfanova_small['group + session + group:session + participant_id'] /
  switcher_bfanova_small['session + group:session + participant_id']
switcher_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
switcher_bfratio_int <- bf_ratio_plot(
  Small = switcher_bf_int_small,
  Medium = switcher_bf_int_medium,
  Wide = switcher_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
switcher_bfratio_int

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switcher_bfratio_int.pdf'),
  plot = switcher_bfratio_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Session
switcher_bfratio_sess <- bf_ratio_plot(
  Small = switcher_bf_sess_small,
  Medium = switcher_bf_sess_medium,
  Wide = switcher_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
switcher_bfratio_sess

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switcher_bfratio_sess.pdf'),
  plot = switcher_bfratio_sess$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Group
switcher_bfratio_group <- bf_ratio_plot(
  Small = switcher_bf_group_small,
  Medium = switcher_bf_group_medium,
  Wide = switcher_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
switcher_bfratio_group

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switcher_bfratio_group.pdf'),
  plot = switcher_bfratio_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior INT ~ posterior #
# ---------------------------- #

switcher_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = ts_agg_er,
  dv = 'switch_cost',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD'
)

switcher_bfanova_medium_posterior_int$summary
switcher_bfanova_medium_posterior_int$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switcher_bfanova_medium_posterior_int.pdf'),
  plot = switcher_bfanova_medium_posterior_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

switcher_bfanova_medium_posterior_session <- bf_contrast_posterior(
  ts_agg_er, 'switch_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen'
)

switcher_bfanova_medium_posterior_session$summary
switcher_bfanova_medium_posterior_session$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switcher_bfanova_medium_posterior_session.pdf'),
  plot = switcher_bfanova_medium_posterior_session$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior GROUP ~ posterior #
# ---------------------------- #

switcher_bfanova_medium_posterior_group <- bf_contrast_posterior(
  ts_agg_er, 'switch_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'orangered'
)

switcher_bfanova_medium_posterior_group$summary
switcher_bfanova_medium_posterior_group$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switcher_bfanova_medium_posterior_group.pdf'),
  plot = switcher_bfanova_medium_posterior_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# BASELINE COMPARISON #
# ======================== #

# Baseline
switcher_baseline_bf <- ttestBF(x = ts_agg_er %>% 
                                   filter(session == 1, 
                                          group == 'CD') %>%
                                   pull(switch_cost),
                                 y = ts_agg_er %>% 
                                   filter(session == 1, 
                                          group == 'KD') %>%
                                   pull(switch_cost),
                                 rscale = 'medium')
switcher_baseline_bf

# - Trimmed incongruence cost ER (ACC) --------------------

# Trimmed ANOVA
incongr_acc_tranova <- WRS2::bwtrim(incongruence_cost ~ group * session, 
                                   id = participant_id, 
                                   data = ts_agg_er,
                                   tr = 0.2
)
incongr_acc_tranova

# Yuen's t-test
incongr_acc_baseline_test <- yuen(incongruence_cost ~ group, 
                                 data = ts_agg_er %>% filter(session == 1),
                                 tr = 0.2)
incongr_acc_baseline_test

# Post-hocs / exploration #
# ------------------------------ #

# Create ordered list for analysis
incongr_acc_post_tr_list <- build_split_list(ts_agg_er, 
                                             'incongruence_cost', c('group', 'session'))

# Interaction contrast on trimmed difference scores
# - Compute difference scores and then trim
incongr_acc_post_tr_int <- bwimcpES(2, 2, incongr_acc_post_tr_list, 
                                    tr=0.2, CI=TRUE, alpha=0.05)
incongr_acc_post_tr_int

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
incongr_session_eff_acc <- bw.es.B(2, 2, incongr_acc_post_tr_list, 
                                  tr = 0.2, POOL = TRUE, OPT = FALSE, 
                                  CI = TRUE, SEED = TRUE, REL.MAG = NULL)
incongr_session_eff_acc

# - Regular incongruence cost ER (ACC) ----------------

incongr_acc_anova <- aov_4(
  incongruence_cost ~ group * session + (session | participant_id),
                         data = ts_agg_er,
                         anova_table = list(es = 'pes'))
incongr_acc_anova

# - Bayesian incongruence cost ER (ACC) -------------------

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
incongrer_bfanova_medium <- anovaBF(
  incongruence_cost ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
incongrer_bfanova_medium

# Interaction BF
incongrer_bf_int_medium <- incongrer_bfanova_medium['group + session + group:session + participant_id'] /
  incongrer_bfanova_medium['group + session + participant_id']
incongrer_bf_int_medium

# Session BF
incongrer_bf_sess_medium <- incongrer_bfanova_medium['group + session + group:session + participant_id'] /
  incongrer_bfanova_medium['group + group:session + participant_id']
incongrer_bf_sess_medium

# Group BF
incongrer_bf_group_medium <- incongrer_bfanova_medium['group + session + group:session + participant_id'] /
  incongrer_bfanova_medium['session + group:session + participant_id']
incongrer_bf_group_medium

# Wide prior #
# ---------------- #

incongrer_bfanova_wide <- anovaBF(
  incongruence_cost ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
incongrer_bfanova_wide

# Interaction BF
incongrer_bf_int_wide <- incongrer_bfanova_wide['group + session + group:session + participant_id'] /
  incongrer_bfanova_wide['group + session + participant_id']
incongrer_bf_int_wide

# Session BF
incongrer_bf_sess_wide <- incongrer_bfanova_wide['group + session + group:session + participant_id'] /
  incongrer_bfanova_wide['group + group:session + participant_id']
incongrer_bf_sess_wide

# Group BF
incongrer_bf_group_wide <- incongrer_bfanova_wide['group + session + group:session + participant_id'] /
  incongrer_bfanova_wide['session + group:session + participant_id']
incongrer_bf_group_wide

# Small prior #
# ---------------- #

incongrer_bfanova_small <- anovaBF(
  incongruence_cost ~ group * session + participant_id,
  data = ts_agg_er %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 100000
)
incongrer_bfanova_small

# Interaction BF
incongrer_bf_int_small <- incongrer_bfanova_small['group + session + group:session + participant_id'] /
  incongrer_bfanova_small['group + session + participant_id']
incongrer_bf_int_small

# Session BF
incongrer_bf_sess_small <- incongrer_bfanova_small['group + session + group:session + participant_id'] /
  incongrer_bfanova_small['group + group:session + participant_id']
incongrer_bf_sess_small

# Group BF
incongrer_bf_group_small <- incongrer_bfanova_small['group + session + group:session + participant_id'] /
  incongrer_bfanova_small['session + group:session + participant_id']
incongrer_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
incongrer_bfratio_int <- bf_ratio_plot(
  Small = incongrer_bf_int_small,
  Medium = incongrer_bf_int_medium,
  Wide = incongrer_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
incongrer_bfratio_int

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrer_bfratio_int.pdf'),
  plot = incongrer_bfratio_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Session
incongrer_bfratio_sess <- bf_ratio_plot(
  Small = incongrer_bf_sess_small,
  Medium = incongrer_bf_sess_medium,
  Wide = incongrer_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = FALSE
)
incongrer_bfratio_sess

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrer_bfratio_sess.pdf'),
  plot = incongrer_bfratio_sess$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Group
incongrer_bfratio_group <- bf_ratio_plot(
  Small = incongrer_bf_group_small,
  Medium = incongrer_bf_group_medium,
  Wide = incongrer_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
incongrer_bfratio_group

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongrer_bfratio_group.pdf'),
  plot = incongrer_bfratio_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior INT ~ posterior #
# ---------------------------- #

incongrer_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = ts_agg_er,
  dv = 'incongruence_cost',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD'
)

incongrer_bfanova_medium_posterior_int$summary
incongrer_bfanova_medium_posterior_int$plot

# Save plot
ggsave(
  filename = file.path(plot_directory,
                       'incongrer_bfanova_medium_posterior_int.pdf'),
  plot = incongrer_bfanova_medium_posterior_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

incongrer_bfanova_medium_posterior_session <- bf_contrast_posterior(
  ts_agg_er, 'incongruence_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen'
)

incongrer_bfanova_medium_posterior_session$summary
incongrer_bfanova_medium_posterior_session$plot

# Save plot
ggsave(
  filename = file.path(plot_directory,
                       'incongrer_bfanova_medium_posterior_session.pdf'),
  plot = incongrer_bfanova_medium_posterior_session$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior GROUP ~ posterior #
# ---------------------------- #

incongrer_bfanova_medium_posterior_group <- bf_contrast_posterior(
  ts_agg_er, 'incongruence_cost', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'orangered'
)

incongrer_bfanova_medium_posterior_group$summary
incongrer_bfanova_medium_posterior_group$plot

# Save plot
ggsave(
  filename = file.path(plot_directory,
                       'incongrer_bfanova_medium_posterior_group.pdf'),
  plot = incongrer_bfanova_medium_posterior_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# BASELINE COMPARISON #
# ======================== #

# Baseline
incongrer_baseline_bf <- ttestBF(x = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'CD') %>%
                               pull(incongruence_cost),
                             y = ts_agg_er %>% 
                               filter(session == 1, 
                                      group == 'KD') %>%
                               pull(incongruence_cost),
                             rscale = 'medium')
incongrer_baseline_bf

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
  filename = file.path(plot_directory, 'gg_ASRS_with.pdf'),
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
  filename = file.path(plot_directory, 'gg_ASRS_without.pdf'),
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
  filename = file.path(plot_directory, 'gg_ASRS_trbar.pdf'),
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
  filename = file.path(plot_directory, 'gg_ASRS_change.pdf'),
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
                           pull(ASRS),
                         questionnaires_long %>% 
                           filter(group == 'KD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           pull(ASRS),
                         tr = 0.2)
ASRS_session_KD

# Session change CD
ASRS_session_CD <- yuend(questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 1) %>%
                           arrange(participant_id) %>%
                           pull(ASRS),
                         questionnaires_long %>% 
                           filter(group == 'CD', 
                                  session == 2) %>%
                           arrange(participant_id) %>%
                           pull(ASRS),
                         tr = 0.2)
ASRS_session_CD

# Effect size of session changes
ASRS_session_eff <- bw.es.B(2, 2, ASRS_post_tr_list, 
                           tr = 0.2, POOL = FALSE, OPT = FALSE, 
                           CI = TRUE, SEED = TRUE, REL.MAG = NULL)
ASRS_session_eff

# - Regular mixed ANOVA ---------------------

ASRS_anova <- aov_4(ASRS ~ group * session + (session | participant_id),
                    data = questionnaires_long,
                    anova_table = list(es = 'pes'))
ASRS_anova

# Baseline comparison #
# --------------------------- #

ASRS_baseline_regtest <- t.test(ASRS ~ group, 
                                data = questionnaires_long %>% 
                                  filter(session == 1),
                                var.equal = FALSE)
ASRS_baseline_regtest

# - Ranked mixed ANOVA -------------------

# Perform aligned-rank transform (ART)
ASRS_ranked <- art(ASRS ~ group * session + (1 | participant_id),
                  data = questionnaires_long)

# Check the result of the transformation 
summary(ASRS_ranked)
anova(ASRS_ranked, response = 'aligned')

# Run ANOVA on the ART data
ASRS_ranked_anova <- anova(ASRS_ranked, response = 'art', 
                           type = 'III', test = 'F')
ASRS_ranked_anova

# Post-hocs #
# ----------------- #

# Simple effects
ASRS_ranked_post_simple <- art.con(ASRS_ranked, ~ group * session)
ASRS_ranked_post_simple

# Group
ASRS_ranked_group <- art.con(ASRS_ranked, ~ group) %>%
  as_tibble() %>%
  mutate(r = t_to_r(t = t.ratio, df_error = df)) %>%
  unnest_wider(r, names_sep = '_') %>%
  rename(r = r_r, CI = r_CI, CI_low = r_CI_low, CI_high = r_CI_high)
ASRS_ranked_group

# Session
ASRS_ranked_session <- art.con(ASRS_ranked, ~ session) %>%
  as_tibble() %>%
  mutate(r = t_to_r(t = t.ratio, df_error = df)) %>%
  unnest_wider(r, names_sep = '_') %>%
  rename(r = r_r, CI = r_CI, CI_low = r_CI_low, CI_high = r_CI_high)
ASRS_ranked_session

# Interaction
ASRS_ranked_post_int <- art.con(ASRS_ranked, ~ group * session, 
                                interaction = TRUE)
ASRS_ranked_post_int <- as_tibble(ASRS_ranked_post_int) %>%
  mutate(r = t_to_r(t = t.ratio, df_error = df)) %>% 
  unnest_wider(r, names_sep = '_') %>%
  dplyr::select(group_pairwise, session_pairwise, estimate, SE, df, t.ratio, p.value,
         r = r_r, CI = r_CI, CI_low = r_CI_low, CI_high = r_CI_high)
print(ASRS_ranked_post_int, width = Inf)

# - Bayesian analysis ----------------------

# Check that there are no NAs
questionnaires_long %>%
  summarise(across(c('PSQI', 'ASRS', 'BDI',
                     'group', 'session', 'participant_id'),
                   ~ sum(is.na(.x))))

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
ASRS_bfanova_medium <- anovaBF(
  ASRS ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
ASRS_bfanova_medium

# Interaction BF
ASRS_bf_int_medium <- ASRS_bfanova_medium['group + session + group:session + participant_id'] /
  ASRS_bfanova_medium['group + session + participant_id']
ASRS_bf_int_medium

# Session BF
ASRS_bf_sess_medium <- ASRS_bfanova_medium['group + session + group:session + participant_id'] /
  ASRS_bfanova_medium['group + group:session + participant_id']
ASRS_bf_sess_medium

# Group BF
ASRS_bf_group_medium <- ASRS_bfanova_medium['group + session + group:session + participant_id'] /
  ASRS_bfanova_medium['session + group:session + participant_id']
ASRS_bf_group_medium

# Wide prior #
# ---------------- #

ASRS_bfanova_wide <- anovaBF(
  ASRS ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
ASRS_bfanova_wide

# Interaction BF
ASRS_bf_int_wide <- ASRS_bfanova_wide['group + session + group:session + participant_id'] /
  ASRS_bfanova_wide['group + session + participant_id']
ASRS_bf_int_wide

# Session BF
ASRS_bf_sess_wide <- ASRS_bfanova_wide['group + session + group:session + participant_id'] /
  ASRS_bfanova_wide['group + group:session + participant_id']
ASRS_bf_sess_wide

# Group BF
ASRS_bf_group_wide <- ASRS_bfanova_wide['group + session + group:session + participant_id'] /
  ASRS_bfanova_wide['session + group:session + participant_id']
ASRS_bf_group_wide

# Small prior #
# ---------------- #

ASRS_bfanova_small <- anovaBF(
  ASRS ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 100000
)
ASRS_bfanova_small

# Interaction BF
ASRS_bf_int_small <- ASRS_bfanova_small['group + session + group:session + participant_id'] /
  ASRS_bfanova_small['group + session + participant_id']
ASRS_bf_int_small

# Session BF
ASRS_bf_sess_small <- ASRS_bfanova_small['group + session + group:session + participant_id'] /
  ASRS_bfanova_small['group + group:session + participant_id']
ASRS_bf_sess_small

# Group BF
ASRS_bf_group_small <- ASRS_bfanova_small['group + session + group:session + participant_id'] /
  ASRS_bfanova_small['session + group:session + participant_id']
ASRS_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
ASRS_bfratio_int <- bf_ratio_plot(
  Small = ASRS_bf_int_small,
  Medium = ASRS_bf_int_medium,
  Wide = ASRS_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
ASRS_bfratio_int

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'ASRS_bfratio_int.pdf'),
  plot = ASRS_bfratio_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Session
ASRS_bfratio_sess <- bf_ratio_plot(
  Small = ASRS_bf_sess_small,
  Medium = ASRS_bf_sess_medium,
  Wide = ASRS_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = FALSE
)
ASRS_bfratio_sess

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'ASRS_bfratio_sess.pdf'),
  plot = ASRS_bfratio_sess$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Group
ASRS_bfratio_group <- bf_ratio_plot(
  Small = ASRS_bf_group_small,
  Medium = ASRS_bf_group_medium,
  Wide = ASRS_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
ASRS_bfratio_group

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'ASRS_bfratio_group.pdf'),
  plot = ASRS_bfratio_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior INT ~ posterior #
# ---------------------------- #

ASRS_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = questionnaires_long,
  dv = 'ASRS',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD'
)

ASRS_bfanova_medium_posterior_int$summary
ASRS_bfanova_medium_posterior_int$plot

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'ASRS_bfanova_medium_posterior_int.pdf'),
  plot = ASRS_bfanova_medium_posterior_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

ASRS_bfanova_medium_posterior_session <- bf_contrast_posterior(
  questionnaires_long, 'ASRS', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen'
)

ASRS_bfanova_medium_posterior_session$summary
ASRS_bfanova_medium_posterior_session$plot

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'ASRS_bfanova_medium_posterior_session.pdf'),
  plot = ASRS_bfanova_medium_posterior_session$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior GROUP ~ posterior #
# ---------------------------- #

ASRS_bfanova_medium_posterior_group <- bf_contrast_posterior(
  questionnaires_long, 'ASRS', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'orangered'
)

ASRS_bfanova_medium_posterior_group$summary
ASRS_bfanova_medium_posterior_group$plot

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'ASRS_bfanova_medium_posterior_group.pdf'),
  plot = ASRS_bfanova_medium_posterior_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# BASELINE COMPARISON #
# ======================== #

# Baseline (medium)
ASRS_baseline_bf <- ttestBF(x = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'CD') %>%
                              pull(ASRS),
                            y = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'KD') %>%
                              pull(ASRS),
                            rscale = 'medium')
ASRS_baseline_bf

# Baseline (r = 0.5)
ASRS_baseline_bf <- ttestBF(x = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'CD') %>%
                              pull(ASRS),
                            y = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'KD') %>%
                              pull(ASRS),
                            rscale = 0.5)
ASRS_baseline_bf

# Baseline (r = 0.3)
ASRS_baseline_bf <- ttestBF(x = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'CD') %>%
                              pull(ASRS),
                            y = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'KD') %>%
                              pull(ASRS),
                            rscale = 0.3)
ASRS_baseline_bf

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
  filename = file.path(plot_directory, 'gg_PSQI_trbar.pdf'),
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
  filename = file.path(plot_directory, 'gg_PSQI_change.pdf'),
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
PSQI_group_eff <- bw.es.A(2, 2, PSQI_post_tr_list, 
                            tr = 0.2, pr = TRUE, fun = ES.summary.CI)
PSQI_group_eff

# - Regular mixed ANOVA ---------------------

PSQI_anova <- aov_4(PSQI ~ group * session + (session | participant_id),
                    data = questionnaires_long,
                    anova_table = list(es = 'pes'))
PSQI_anova

# Post-hocs for reporting the DiD #
# ----------------------------------- #

# Grid
PSQI_emm <- emmeans(PSQI_anova, ~ group * session)
PSQI_emm

# Comparisons
PSQI_post_int <- contrast(PSQI_emm, interaction = c('pairwise', 'pairwise'))
PSQI_post_int

# - Ranked mixed ANOVA -------------------

# Perform aligned-rank transform (ART)
PSQI_ranked <- art(PSQI ~ group * session + (1 | participant_id),
                   data = questionnaires_long)

# Check the result of the transformation 
summary(PSQI_ranked)
anova(PSQI_ranked, response = 'aligned')

# Run ANOVA on the ART data
PSQI_ranked_anova <- anova(PSQI_ranked, response = 'art', 
                           type = 'III', test = 'F')
PSQI_ranked_anova

# Compare difference scores distributions #
# --------------------------------------------- #

# Permutation test - testing the equality of distributions
PSQI_perm <- permg(PSQI_wide %>% filter(group == 'CD') %>% pull(PSQI_change),
                   PSQI_wide %>% filter(group == 'KD') %>% pull(PSQI_change),
                   alpha = 0.05, est = var, nboot = 10000)
PSQI_perm

# Comparing quantiles
PSQI_quant <- WRS2::qcomhd(PSQI_change ~ group, PSQI_wide, 
                           q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                           nboot = 2000, alpha = 0.05, ADJ.CI = TRUE)
PSQI_quant

# - Bayesian analysis -------------------------

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
PSQI_bfanova_medium <- anovaBF(
  PSQI ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
PSQI_bfanova_medium

# Interaction BF
PSQI_bf_int_medium <- PSQI_bfanova_medium['group + session + group:session + participant_id'] /
  PSQI_bfanova_medium['group + session + participant_id']
PSQI_bf_int_medium

# Session BF
PSQI_bf_sess_medium <- PSQI_bfanova_medium['group + session + group:session + participant_id'] /
  PSQI_bfanova_medium['group + group:session + participant_id']
PSQI_bf_sess_medium

# Group BF
PSQI_bf_group_medium <- PSQI_bfanova_medium['group + session + group:session + participant_id'] /
  PSQI_bfanova_medium['session + group:session + participant_id']
PSQI_bf_group_medium

# Wide prior #
# ---------------- #

PSQI_bfanova_wide <- anovaBF(
  PSQI ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
PSQI_bfanova_wide

# Interaction BF
PSQI_bf_int_wide <- PSQI_bfanova_wide['group + session + group:session + participant_id'] /
  PSQI_bfanova_wide['group + session + participant_id']
PSQI_bf_int_wide

# Session BF
PSQI_bf_sess_wide <- PSQI_bfanova_wide['group + session + group:session + participant_id'] /
  PSQI_bfanova_wide['group + group:session + participant_id']
PSQI_bf_sess_wide

# Group BF
PSQI_bf_group_wide <- PSQI_bfanova_wide['group + session + group:session + participant_id'] /
  PSQI_bfanova_wide['session + group:session + participant_id']
PSQI_bf_group_wide

# Small prior #
# ---------------- #

PSQI_bfanova_small <- anovaBF(
  PSQI ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 100000
)
PSQI_bfanova_small

# Interaction BF
PSQI_bf_int_small <- PSQI_bfanova_small['group + session + group:session + participant_id'] /
  PSQI_bfanova_small['group + session + participant_id']
PSQI_bf_int_small

# Session BF
PSQI_bf_sess_small <- PSQI_bfanova_small['group + session + group:session + participant_id'] /
  PSQI_bfanova_small['group + group:session + participant_id']
PSQI_bf_sess_small

# Group BF
PSQI_bf_group_small <- PSQI_bfanova_small['group + session + group:session + participant_id'] /
  PSQI_bfanova_small['session + group:session + participant_id']
PSQI_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
PSQI_bfratio_int <- bf_ratio_plot(
  Small = PSQI_bf_int_small,
  Medium = PSQI_bf_int_medium,
  Wide = PSQI_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
PSQI_bfratio_int

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'PSQI_bfratio_int.pdf'),
  plot = PSQI_bfratio_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Session
PSQI_bfratio_sess <- bf_ratio_plot(
  Small = PSQI_bf_sess_small,
  Medium = PSQI_bf_sess_medium,
  Wide = PSQI_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = FALSE
)
PSQI_bfratio_sess

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'PSQI_bfratio_sess.pdf'),
  plot = PSQI_bfratio_sess$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Group
PSQI_bfratio_group <- bf_ratio_plot(
  Small = PSQI_bf_group_small,
  Medium = PSQI_bf_group_medium,
  Wide = PSQI_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
PSQI_bfratio_group

# Save plot with optimal publication dimensions
ggsave(
  filename = file.path(plot_directory, 'PSQI_bfratio_group.pdf'),
  plot = PSQI_bfratio_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior INT ~ posterior #
# ---------------------------- #

PSQI_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = questionnaires_long,
  dv = 'PSQI',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD'
)

PSQI_bfanova_medium_posterior_int$summary
PSQI_bfanova_medium_posterior_int$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'PSQI_bfanova_medium_posterior_int.pdf'),
  plot = PSQI_bfanova_medium_posterior_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

PSQI_bfanova_medium_posterior_session <- bf_contrast_posterior(
  questionnaires_long, 'PSQI', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen'
)

PSQI_bfanova_medium_posterior_session$summary
PSQI_bfanova_medium_posterior_session$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'PSQI_bfanova_medium_posterior_session.pdf'),
  plot = PSQI_bfanova_medium_posterior_session$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

PSQI_bfanova_medium_posterior_group <- bf_contrast_posterior(
  questionnaires_long, 'PSQI', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'orangered2'
)

PSQI_bfanova_medium_posterior_group$summary
PSQI_bfanova_medium_posterior_group$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'PSQI_bfanova_medium_posterior_group.pdf'),
  plot = PSQI_bfanova_medium_posterior_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# BASELINE COMPARISON #
# ======================== #

PSQI_baseline_bf <- ttestBF(x = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'CD') %>%
                              pull(PSQI),
                            y = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'KD') %>%
                              pull(PSQI),
                            rscale = 'medium')
PSQI_baseline_bf

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
  filename = file.path(plot_directory, 'gg_BDI_with.pdf'),
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
  filename = file.path(plot_directory, 'gg_BDI_without.pdf'),
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
  filename = file.path(plot_directory, 'gg_BDI_trbar.pdf'),
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
  filename = file.path(plot_directory, 'gg_BDI_change.pdf'),
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
                                                  group == 'CD') %>% pull(BDI),
                   questionnaires_long %>% filter(session == 1, 
                                                  group == 'KD') %>% pull(BDI),
                   alpha = 0.05, est = var, nboot = 10000)
BDI_perm

# Kolmogorov-Smirnov 
BDI_ks_baseline <- ks.test(questionnaires_long %>% 
                             filter(session == 1,
                                    group == 'CD') %>% 
                             pull(BDI),
                           questionnaires_long %>% 
                             filter(session == 1,
                                    group == 'KD') 
                           %>% pull(BDI),
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

# - Regular mixed ANOVA ---------------------

BDI_anova <- aov_4(BDI ~ group * session + (session | participant_id),
                    data = questionnaires_long,
                    anova_table = list(es = 'pes'))
BDI_anova

# Post-hocs for reporting the DiD #
# ----------------------------------- #

# Grid
BDI_emm <- emmeans(BDI_anova, ~ group * session)
BDI_emm

# Comparisons
BDI_post_int <- contrast(BDI_emm, interaction = c('pairwise', 'pairwise'))
BDI_post_int

# Baseline
BDI_baseline_test_reg <- emmeans(BDI_anova, ~ group | session) %>%
  contrast('pairwise')
BDI_baseline_test_reg

# - Ranked mixed ANOVA -------------------

# Perform aligned-rank transform (ART)
BDI_ranked <- art(BDI ~ group * session + (1 | participant_id),
                   data = questionnaires_long)

# Check the result of the transformation 
summary(BDI_ranked)
anova(BDI_ranked, response = 'aligned')

# Run ANOVA on the ART data
BDI_ranked_anova <- anova(BDI_ranked, response = 'art', 
                           type = 'III', test = 'F')
BDI_ranked_anova

# Compare difference scores distributions #
# --------------------------------------------- #

# Comparing quantiles
BDI_quant <- WRS2::qcomhd(BDI_change ~ group, BDI_wide, q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                           nboot = 2000, alpha = 0.05, ADJ.CI = TRUE)
BDI_quant

# - Bayesian analysis -------------------------------

# Medium prior #
# ---------------- #

# Mixed ANOVA model set (random = participant)
BDI_bfanova_medium <- anovaBF(
  BDI ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'medium',
  rscaleRandom = 'nuisance',
  iterations = 100000,
)
BDI_bfanova_medium

# Interaction BF
BDI_bf_int_medium <- BDI_bfanova_medium['group + session + group:session + participant_id'] /
  BDI_bfanova_medium['group + session + participant_id']
BDI_bf_int_medium

# Session BF
BDI_bf_sess_medium <- BDI_bfanova_medium['group + session + group:session + participant_id'] /
  BDI_bfanova_medium['group + group:session + participant_id']
BDI_bf_sess_medium

# Group BF
BDI_bf_group_medium <- BDI_bfanova_medium['group + session + group:session + participant_id'] /
  BDI_bfanova_medium['session + group:session + participant_id']
BDI_bf_group_medium

# Wide prior #
# ---------------- #

BDI_bfanova_wide <- anovaBF(
  BDI ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 'wide',
  rscaleRandom = 'nuisance',
  iterations = 100000
)
BDI_bfanova_wide

# Interaction BF
BDI_bf_int_wide <- BDI_bfanova_wide['group + session + group:session + participant_id'] /
  BDI_bfanova_wide['group + session + participant_id']
BDI_bf_int_wide

# Session BF
BDI_bf_sess_wide <- BDI_bfanova_wide['group + session + group:session + participant_id'] /
  BDI_bfanova_wide['group + group:session + participant_id']
BDI_bf_sess_wide

# Group BF
BDI_bf_group_wide <- BDI_bfanova_wide['group + session + group:session + participant_id'] /
  BDI_bfanova_wide['session + group:session + participant_id']
BDI_bf_group_wide

# Small prior #
# ---------------- #

BDI_bfanova_small <- anovaBF(
  BDI ~ group * session + participant_id,
  data = questionnaires_long %>% 
    mutate(participant_id = factor(participant_id)) %>%
    as.data.frame(),
  whichRandom = 'participant_id',       
  whichModels = 'all',               
  rscaleFixed = 0.25,
  rscaleRandom = 'nuisance',
  iterations = 100000
)
BDI_bfanova_small

# Interaction BF
BDI_bf_int_small <- BDI_bfanova_small['group + session + group:session + participant_id'] /
  BDI_bfanova_small['group + session + participant_id']
BDI_bf_int_small

# Session BF
BDI_bf_sess_small <- BDI_bfanova_small['group + session + group:session + participant_id'] /
  BDI_bfanova_small['group + group:session + participant_id']
BDI_bf_sess_small

# Group BF
BDI_bf_group_small <- BDI_bfanova_small['group + session + group:session + participant_id'] /
  BDI_bfanova_small['session + group:session + participant_id']
BDI_bf_group_small

# Sensitivity ~ BFs #
# ----------------- #

# Interaction
BDI_bfratio_int <- bf_ratio_plot(
  Small = BDI_bf_int_small,
  Medium = BDI_bf_int_medium,
  Wide = BDI_bf_int_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
BDI_bfratio_int$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_bfratio_int.pdf'),
  plot = BDI_bfratio_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Session
BDI_bfratio_sess <- bf_ratio_plot(
  Small = BDI_bf_sess_small,
  Medium = BDI_bf_sess_medium,
  Wide = BDI_bf_sess_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = FALSE
)
BDI_bfratio_sess

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_bfratio_sess.pdf'),
  plot = BDI_bfratio_sess$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Group
BDI_bfratio_group <- bf_ratio_plot(
  Small = BDI_bf_group_small,
  Medium = BDI_bf_group_medium,
  Wide = BDI_bf_group_wide,
  log10_scale = TRUE,
  annotate = TRUE,
  show_bands = TRUE
)
BDI_bfratio_group

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_bfratio_group.pdf'),
  plot = BDI_bfratio_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior INT ~ posterior #
# ---------------------------- #

BDI_bfanova_medium_posterior_int <- bf_contrast_posterior(
  data = questionnaires_long,
  dv = 'BDI',
  factor_a = 'group',
  factor_b = 'session',
  id_var = 'participant_id',
  a_levels = c('CD', 'KD'),
  b_levels = c('1', '2'),
  contrast = 'interaction',
  rscale_fixed = 'medium',
  iterations = 100000,
  fill = 'steelblue',
  ci_color = 'red4',
  mean_color = 'black',
  title = 'Posterior of the population-average DiD'
)

BDI_bfanova_medium_posterior_int$summary
BDI_bfanova_medium_posterior_int$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_bfanova_medium_posterior_int.pdf'),
  plot = BDI_bfanova_medium_posterior_int$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior SESSION ~ posterior #
# ---------------------------- #

BDI_bfanova_medium_posterior_session <- bf_contrast_posterior(
  questionnaires_long, 'BDI', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_b', main_effect_source = 'auto',
  fill = 'darkgreen'
)

BDI_bfanova_medium_posterior_session$summary
BDI_bfanova_medium_posterior_session$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_bfanova_medium_posterior_session.pdf'),
  plot = BDI_bfanova_medium_posterior_session$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Medium prior GROUP ~ posterior #
# ---------------------------- #

BDI_bfanova_medium_posterior_group <- bf_contrast_posterior(
  questionnaires_long, 'BDI', 'group', 'session', 'participant_id',
  a_levels = c('CD', 'KD'), b_levels = c('1', '2'),
  contrast = 'main_a', main_effect_source = 'auto',
  fill = 'orangered2'
)

BDI_bfanova_medium_posterior_group$summary
BDI_bfanova_medium_posterior_group$plot

# Save plot
ggsave(
  filename = file.path(plot_directory, 'BDI_bfanova_medium_posterior_group.pdf'),
  plot = BDI_bfanova_medium_posterior_group$plot,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# BASELINE COMPARISON #
# ======================== #

BDI_baseline_bf <- ttestBF(x = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'CD') %>%
                              pull(BDI),
                            y = questionnaires_long %>% 
                              filter(session == 1, 
                                     group == 'KD') %>%
                              pull(BDI),
                            rscale = 'medium')
BDI_baseline_bf

# ANCOVA ANALYSIS ------------------
# - Prepare data for sweep -----------------

# NOTE: the following tibble contains only the participants who completed 
# both the questionnaires and the task-switching

# Create tibble for correlations
explore_tibble <- ts_agg_rt %>%
  left_join(ts_agg_er, by = c('participant_id', 'group', 'session'), 
            suffix = c('_rt', '_acc')) %>%
  left_join(questionnaires_long, by = c('participant_id', 'group', 'session')) %>%
  dplyr::select(-starts_with(c('repeat_congruent', 'repeat_incongruent', 
                               'switch_congruent', 'switch_incongruent'))) %>%
  pivot_wider(
    id_cols = c(participant_id, group, cohort, age),
    names_from = session,
    values_from = -c(participant_id, group, cohort, age, session),
    names_glue = '{.value}_{session}'
  )
glimpse(explore_tibble)

# - Sweep for covariates -----------------------

# Run the screen #
# ---------------------- #

# Select response variables
response_vars_covsweep <- colnames(explore_tibble) %>%
  str_subset('_(2)$')

# Run the screen
partial_coeffs_all <- run_screen(explore_tibble,
                                 response_vars = response_vars_covsweep,
                                 must_have = c('age', 'BMI_1'),
                                 baseline = FALSE,
                                 method_cor = 'spearman',
                                 pre_vars_reg = '_(1)$')

# Extract
partial_coeffs_extract <- partial_coeffs_all %>%
  filter(!is.na(cor_ttl)) %>%
  filter(abs(cor_ttl) >= 0.15 | p_ttl < 0.1) %>%
  arrange(outcome, desc(abs(cor_ttl)))
print(partial_coeffs_extract, n = Inf, width = Inf)

# Screen 2 #
# ---------------------- #

partial_coeffs_extract_2 <- run_screen(explore_tibble,
                                       response_vars = response_vars_covsweep,
                                       must_have = c('age', 
                                                     'BMI_1',
                                                     'ASRS_1', 
                                                     'BDI_1',
                                                     'PSQI_1'),
                                       baseline = FALSE,
                                       method_cor = 'spearman',
                                       pre_vars_reg = '_(1)$'
                                ) %>%
  filter(!is.na(cor_ttl)) %>%
  filter(abs(cor_ttl) >= 0.15 | p_ttl < 0.1) %>%
  arrange(outcome, desc(abs(cor_ttl)))
print(partial_coeffs_extract_2, n = Inf, width = Inf)

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

# Contains only the intersection of the participants (see previous section)
explore_tibble_centred <- explore_tibble %>%
  mutate(across(ends_with('_1'), ~ as.numeric(scale(.x, scale = FALSE))))
glimpse(explore_tibble_centred)

# - Balancing -------------------

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

# - ASRS adjusted ---------------------

# REGULAR ANCOVA #
# ======================= #

# ANCOVA model
# - Adding switch_cost_rt_1 doesn't change anything
ASRS_ancova <- lm_robust(ASRS_2 ~ ASRS_1 + group, 
                         data = explore_tibble_quest_centred,
                         se_type = 'HC3')
summary(ASRS_ancova)
car::Anova(ASRS_ancova, type = 'III', test.statistic = 'F')

# ANCOVA slopes check
ASRS_ancova_sl <- lm_robust(ASRS_2 ~ ASRS_1 * group, 
                         data = explore_tibble_quest_centred,
                         se_type = 'HC3')
summary(ASRS_ancova_sl)

# Assess fit #
# ----------------------- #

# Get residuals
ASRS_ancova_res <- as.numeric(scale(residuals(lm(ASRS_2 ~ ASRS_1 + group, 
                                data = explore_tibble_quest_centred))))
# Get predicted values
ASRS_ancova_fit <- fitted(lm(ASRS_2 ~ ASRS_1 + group, 
                             data = explore_tibble_quest_centred))

# Final plot
plot(ASRS_ancova_fit, ASRS_ancova_res)

# Marginal means #
# ----------------------- #

# Check estimated means
ASRS_anc_emm <- emmeans(ASRS_ancova, ~ group, vcov. = vcov(ASRS_ancova))
ASRS_anc_emm

ASRS_anc_emm_rng <- range(explore_tibble_quest_centred$ASRS_1, na.rm = TRUE)
ASRS_anc_emm_pred <- emmip(
  ASRS_ancova,
  formula = group ~ ASRS_1, 
  at = list(ASRS_1 = seq(ASRS_anc_emm_rng[1], ASRS_anc_emm_rng[2], length.out = 100)),
  CIs = TRUE,
  PIs = FALSE, 
  type = 'response',
  vcov. = vcov(ASRS_ancova),
  plotit = FALSE
)

ggplot(ASRS_anc_emm_pred, aes(x = ASRS_1, y = yvar, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = 'ASRS Pretest (Centred)',
       y = 'Predicted ASRS Posttest',
       colour = 'Group', fill = 'Group',
       title = 'Model-implied ASRS_2 by ASRS_1 (ANCOVA with HC3)') +
  theme_minimal()

# Partial eta squared #
# ----------------------------- #

# Group
ASRS_petasq_group <- 
  ASRS_ancova$statistic['group1']^2 / 
  (ASRS_ancova$statistic['group1']^2 + ASRS_ancova$df.residual)
ASRS_petasq_group

# Pretest
ASRS_petasq_pretest <- 
  ASRS_ancova$statistic['ASRS_1']^2 / 
  (ASRS_ancova$statistic['ASRS_1']^2 + ASRS_ancova$df.residual)
ASRS_petasq_pretest

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

# Regular Yuen
ASRS_respost_yuen <- yuen(ASRS_2_res ~ group,
                              data = balanced_ASRS,
                              tr = 0.2)
ASRS_respost_yuen

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

# Quantiles
ASRS_respost_quant <- qcomhd(ASRS_2_res ~ group,
                             data = balanced_ASRS,
                             q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                             nboot = 50000, ADJ.CI = TRUE)
ASRS_respost_quant

# Regular t-test
ASRS_respost_ttest <- t.test(balanced_ASRS %>% filter(group == 'CD') %>%
                               pull(ASRS_2_res),
                             balanced_ASRS %>% filter(group == 'KD') %>%
                               pull(ASRS_2_res))
ASRS_respost_ttest

# - PSQI adjusted ------------------

# ANCOVA model
PSQI_ancova <- lm_robust(PSQI_2 ~ PSQI_1 + group, 
                         data = explore_tibble_quest_centred,
                         se_type = 'HC3')
summary(PSQI_ancova)
car::Anova(PSQI_ancova, type = 'III', test.statistic = 'F')

# ANCOVA slopes check
PSQI_ancova_sl <- lm_robust(PSQI_2 ~ PSQI_1 * group, 
                            data = explore_tibble_quest_centred,
                            se_type = 'HC3')
summary(PSQI_ancova_sl)

# Assess fit #
# ----------------------- #

# Get residuals
PSQI_ancova_res <- as.numeric(scale(residuals(lm(PSQI_2 ~ PSQI_1 + group, 
                                                 data = explore_tibble_quest_centred))))
# Get predicted values
PSQI_ancova_fit <- fitted(lm(PSQI_2 ~ PSQI_1 + group, 
                             data = explore_tibble_quest_centred))

# Final plot
plot(PSQI_ancova_fit, PSQI_ancova_res)

# Estimated marginal means #
# ------------------------------- #

# Check estimated means
PSQI_anc_emm <- emmeans(PSQI_ancova, ~ group, vcov. = vcov(PSQI_ancova))
PSQI_anc_emm

PSQI_anc_emm_rng <- range(explore_tibble_quest_centred$PSQI_1, na.rm = TRUE)
PSQI_anc_emm_pred <- emmip(
  PSQI_ancova,
  formula = group ~ PSQI_1, 
  at = list(PSQI_1 = seq(PSQI_anc_emm_rng[1], PSQI_anc_emm_rng[2], length.out = 100)),
  CIs = TRUE,
  PIs = FALSE, 
  type = 'response',
  vcov. = vcov(PSQI_ancova),
  plotit = FALSE
)

# Plot
ggplot(PSQI_anc_emm_pred, aes(x = PSQI_1, y = yvar, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = 'PSQI Pretest (Centred)',
       y = 'Predicted PSQI Posttest',
       colour = 'Group', fill = 'Group',
       title = 'Model-implied PSQI_2 by PSQI_1 (ANCOVA with HC3)') +
  theme_minimal()

# Partial eta squared #
# ----------------------------- #

# Group
PSQI_petasq_group <- 
  PSQI_ancova$statistic['group1']^2 / 
  (PSQI_ancova$statistic['group1']^2 + PSQI_ancova$df.residual)
PSQI_petasq_group

# Pretest
PSQI_petasq_pretest <- 
  PSQI_ancova$statistic['PSQI_1']^2 / 
  (PSQI_ancova$statistic['PSQI_1']^2 + PSQI_ancova$df.residual)
PSQI_petasq_pretest

# - BDI adjusted --------------------

# ANCOVA model
BDI_ancova <- lm_robust(BDI_2 ~ BDI_1 + group, 
                         data = explore_tibble_quest_centred,
                         se_type = 'HC3')
summary(BDI_ancova)
car::Anova(BDI_ancova, type = 'III', test.statistic = 'F')

# ANCOVA slopes check
BDI_ancova_sl <- lm_robust(BDI_2 ~ BDI_1 * group, 
                            data = explore_tibble_quest_centred,
                            se_type = 'HC3')
summary(BDI_ancova_sl)

# Assess fit #
# ----------------------- #

# Get residuals
BDI_ancova_res <- as.numeric(scale(residuals(lm(BDI_2 ~ BDI_1 + group, 
                                                 data = explore_tibble_quest_centred))))
# Get predicted values
BDI_ancova_fit <- fitted(lm(BDI_2 ~ BDI_1 + group, 
                             data = explore_tibble_quest_centred))

# Final plot
plot(BDI_ancova_fit, BDI_ancova_res)

# Estimated marginal means #
# ------------------------------- #

# Check estimated means
BDI_anc_emm <- emmeans(BDI_ancova, ~ group, vcov. = vcov(BDI_ancova))
BDI_anc_emm

BDI_anc_emm_rng <- range(explore_tibble_quest_centred$BDI_1, na.rm = TRUE)
BDI_anc_emm_pred <- emmip(
  BDI_ancova,
  formula = group ~ BDI_1, 
  at = list(BDI_1 = seq(BDI_anc_emm_rng[1], BDI_anc_emm_rng[2], length.out = 100)),
  CIs = TRUE,
  PIs = FALSE, 
  type = 'response',
  vcov. = vcov(BDI_ancova),
  plotit = FALSE
)

ggplot(BDI_anc_emm_pred, aes(x = BDI_1, y = yvar, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = 'BDI Pretest (Centred)',
       y = 'Predicted BDI Posttest',
       colour = 'Group', fill = 'Group',
       title = 'Model-implied BDI_2 by BDI_1 (ANCOVA with HC3)') +
  theme_minimal()

# Partial eta squared #
# ----------------------------- #

# Group
BDI_petasq_group <- 
  BDI_ancova$statistic['group1']^2 / 
  (BDI_ancova$statistic['group1']^2 + BDI_ancova$df.residual)
BDI_petasq_group

# Pretest
BDI_petasq_pretest <- 
  BDI_ancova$statistic['BDI_1']^2 / 
  (BDI_ancova$statistic['BDI_1']^2 + BDI_ancova$df.residual)
BDI_petasq_pretest

# - Switch cost RT adjusted ------------------------

# ANCOVA model
# - Adding switch_cost_acc_1 doesn't improve the model
switch_cost_rt_ancova <- lm_robust(
  switch_cost_rt_2 ~ switch_cost_rt_1 + group,
  data = explore_tibble_centred,
  se_type = 'HC3'
)
summary(switch_cost_rt_ancova)
car::Anova(switch_cost_rt_ancova, type = 'III', test.statistic = 'F')

# ANCOVA slopes check
switch_cost_rt_ancova_sl <- lm_robust(
  switch_cost_rt_2 ~ switch_cost_rt_1 * group,
  data = explore_tibble_centred,
  se_type = 'HC3')
summary(switch_cost_rt_ancova_sl)

# Assess fit #
# ----------------------- #

# Get residuals
switch_cost_rt_ancova_res <- as.numeric(
  scale(residuals(lm(switch_cost_rt_2 ~ switch_cost_rt_1 + group,
                     data = explore_tibble_centred))))
# Get predicted values
switch_cost_rt_ancova_fit <- fitted(
  lm(switch_cost_rt_2 ~ switch_cost_rt_1 + group,
     data = explore_tibble_centred))

# Final plot
plot(switch_cost_rt_ancova_fit, switch_cost_rt_ancova_res)

# Estimated marginal means #
# ------------------------------- #

# Check estimated means
switch_cost_rt_anc_emm <- emmeans(switch_cost_rt_ancova, ~ group, 
                                        vcov. = vcov(switch_cost_rt_ancova))
switch_cost_rt_anc_emm

switch_anc_rt_emm_rng <- range(explore_tibble_centred$switch_cost_rt_1, 
                            na.rm = TRUE)
switch_anc_rt_emm_pred <- emmip(
  switch_cost_rt_ancova,
  formula = group ~ switch_cost_rt_1, 
  at = list(switch_cost_rt_1 = seq(switch_anc_rt_emm_rng[1], 
                                   switch_anc_rt_emm_rng[2], 
                                   length.out = 100)),
  CIs = TRUE,
  PIs = FALSE, 
  type = 'response',
  vcov. = vcov(switch_cost_rt_ancova),
  plotit = FALSE
)

# ANCOVA plot
gg_switch_anc_rt <- ggplot(switch_anc_rt_emm_pred, aes(x = switch_cost_rt_1, 
                                y = yvar, 
                                colour = group, 
                                fill = group)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Pretest Switch Cost Reaction Time (ms) (Centred)',
       y = 'Posttest Switch Cost Reaction Time (ms)',
       colour = 'Group', fill = 'Group') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switch_anc_rt.pdf'),
  plot = gg_switch_anc_rt,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Partial eta squared #
# ----------------------------- #

# Group
switchrt_petasq_group <- 
  switch_cost_rt_ancova$statistic['group1']^2 / 
  (switch_cost_rt_ancova$statistic['group1']^2 + switch_cost_rt_ancova$df.residual)
switchrt_petasq_group

# Pretest
switchrt_petasq_pretest <- 
  switch_cost_rt_ancova$statistic['switch_cost_rt_1']^2 / 
  (switch_cost_rt_ancova$statistic['switch_cost_rt_1']^2 + switch_cost_rt_ancova$df.residual)
switchrt_petasq_pretest

# - Switch cost ER (ACC) adjusted ------------------------

# ANCOVA model
switch_cost_acc_ancova <- lm_robust(
  switch_cost_acc_2 ~ switch_cost_acc_1 + group,
  data = explore_tibble_centred,
  se_type = 'HC3'
)
summary(switch_cost_acc_ancova)
car::Anova(switch_cost_acc_ancova, type = 'III', test.statistic = 'F')

# ANCOVA slopes check
switch_cost_acc_ancova_sl <- lm_robust(
  switch_cost_acc_2 ~ switch_cost_acc_1 * group,
  data = explore_tibble_centred,
  se_type = 'HC3')
summary(switch_cost_acc_ancova_sl)

# Assess fit #
# ----------------------- #

# Get residuals
switch_cost_acc_ancova_res <- as.numeric(
  scale(residuals(lm(switch_cost_acc_2 ~ switch_cost_acc_1 + group,
                     data = explore_tibble_centred))))
# Get predicted values
switch_cost_acc_ancova_fit <- fitted(
  lm(switch_cost_acc_2 ~ switch_cost_acc_1 + group,
     data = explore_tibble_centred))

# Final plot
plot(switch_cost_acc_ancova_fit, switch_cost_acc_ancova_res)

# Estimated marginal means #
# ------------------------------- #

# Check estimated means
switch_cost_acc_anc_emm <- emmeans(switch_cost_acc_ancova, ~ group, 
                                  vcov. = vcov(switch_cost_acc_ancova))
switch_cost_acc_anc_emm

switch_anc_acc_emm_rng <- range(explore_tibble_centred$switch_cost_acc_1, 
                            na.rm = TRUE)
switch_anc_acc_emm_pred <- emmip(
  switch_cost_acc_ancova,
  formula = group ~ switch_cost_acc_1, 
  at = list(switch_cost_acc_1 = seq(switch_anc_acc_emm_rng[1], 
                                    switch_anc_acc_emm_rng[2], 
                                   length.out = 100)),
  CIs = TRUE,
  PIs = FALSE, 
  type = 'response',
  vcov. = vcov(switch_cost_acc_ancova),
  plotit = FALSE
)

# ANCOVA plot
gg_switch_anc_acc <- ggplot(switch_anc_acc_emm_pred, aes(x = switch_cost_acc_1, 
                                y = yvar, 
                                colour = group, 
                                fill = group)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  scale_fill_manual(values = pal, labels = c('Clean Diet', 'Ketogenic Diet')) +
  labs(x = 'Pretest Switch Cost Accuracy (p.p.) (Centred)',
       y = 'Posttest Switch Cost Accuracy (p.p.)',
       colour = 'Group', fill = 'Group') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'switch_anc_acc.pdf'),
  plot = gg_switch_anc_acc,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Partial eta squared #
# ----------------------------- #

# Group
switcher_petasq_group <- 
  switch_cost_acc_ancova$statistic['group1']^2 / 
  (switch_cost_acc_ancova$statistic['group1']^2 + switch_cost_acc_ancova$df.residual)
switcher_petasq_group

# Pretest
switcher_petasq_pretest <- 
  switch_cost_acc_ancova$statistic['switch_cost_acc_1']^2 / 
  (switch_cost_acc_ancova$statistic['switch_cost_acc_1']^2 + switch_cost_acc_ancova$df.residual)
switcher_petasq_pretest

# - Incongruence cost RT adjusted ------------------------

# ANCOVA model
# - Adding switch_cost_acc_1 doesn't help; incongruence_cost_acc_1 renders 
# average_rt_1 significant (R2 doesn't really change)
incongruence_cost_rt_ancova <- lm_robust(
  incongruence_cost_rt_2 ~ 
    incongruence_cost_rt_1 +
    average_rt_1 +
    group,
  data = explore_tibble_centred,
  se_type = 'HC3'
)
summary(incongruence_cost_rt_ancova)
car::Anova(incongruence_cost_rt_ancova, type = 'III', test.statistic = 'F')

# ANCOVA slopes check
incongruence_cost_rt_ancova_sl <- lm_robust(
  incongruence_cost_rt_2 ~ 
    incongruence_cost_rt_1 * average_rt_1 * group,
    data = explore_tibble_centred,
    se_type = 'HC3')
summary(incongruence_cost_rt_ancova_sl)

# Assess fit #
# ----------------------- #

# Get residuals
incongruence_cost_rt_ancova_res <- as.numeric(
  scale(residuals(lm(incongruence_cost_rt_2 ~ 
                         incongruence_cost_rt_1 +
                         average_rt_1 +
                         group,
                     data = explore_tibble_centred))))

# Get predicted values
incongruence_cost_rt_ancova_fit <- fitted(lm(incongruence_cost_rt_2 ~ 
                                                 incongruence_cost_rt_1 +
                                                 average_rt_1 +
                                                 group, 
                             data = explore_tibble_centred))

# Final plot
plot(incongruence_cost_rt_ancova_fit, incongruence_cost_rt_ancova_res)

# Estimated marginal means #
# ------------------------------- #

# Check estimated means
incongruence_cost_rt_anc_emm <- emmeans(incongruence_cost_rt_ancova, ~ group, 
                                        vcov. = vcov(incongruence_cost_rt_ancova))
incongruence_cost_rt_anc_emm

# Normal bar plot (2 covariates)
incgr_cost_rt_anc_tbl <- as_tibble(incongruence_cost_rt_anc_emm)
gg_incongr_anc_rt <- ggplot(incgr_cost_rt_anc_tbl, 
                            aes(x = group, y = emmean, fill = group)) +
  geom_col(alpha = 0.9, width = 0.4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.08, linewidth = 0.75) +
  scale_fill_manual(values = pal, guide = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  coord_cartesian(ylim = c(0, NA)) +   
  labs(x = 'Group',
       y = 'Posttest Incongruence Cost Reaction Time (ms)') + 
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  theme(legend.position = 'none') +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongr_anc_rt.pdf'),
  plot = gg_incongr_anc_rt,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Partial eta squared #
# ----------------------------- #

# Group
incongrrt_petasq_group <- 
  incongruence_cost_rt_ancova$statistic['group1']^2 / 
  (incongruence_cost_rt_ancova$statistic['group1']^2 + incongruence_cost_rt_ancova$df.residual)
incongrrt_petasq_group

# Pretest
incongrrt_petasq_pretest <- 
  incongruence_cost_rt_ancova$statistic['incongruence_cost_rt_1']^2 / 
  (incongruence_cost_rt_ancova$statistic['incongruence_cost_rt_1']^2 + incongruence_cost_rt_ancova$df.residual)
incongrrt_petasq_pretest

# Covariate
incongrrt_petasq_cov <- 
  incongruence_cost_rt_ancova$statistic['average_rt_1']^2 / 
  (incongruence_cost_rt_ancova$statistic['average_rt_1']^2 + incongruence_cost_rt_ancova$df.residual)
incongrrt_petasq_cov

# - Incongruence cost ER (ACC) adjusted ------------------------

# ANCOVA model
incongruence_cost_acc_ancova <- lm_robust(
  incongruence_cost_acc_2 ~ 
    incongruence_cost_acc_1 +
    switch_cost_acc_1 +
    average_acc_1 +
    group,
  data = explore_tibble_centred,
  se_type = 'HC3'
)
summary(incongruence_cost_acc_ancova)
car::Anova(incongruence_cost_acc_ancova, type = 'III', test.statistic = 'F')

# ANCOVA slopes check
incongruence_cost_acc_ancova_sl <- lm_robust(
  incongruence_cost_acc_2 ~ 
    incongruence_cost_acc_1 * switch_cost_acc_1 * average_acc_1 * group,
  data = explore_tibble_centred,
  se_type = 'HC3')
summary(incongruence_cost_acc_ancova_sl)

# Assess fit #
# ----------------------- #

# Get residuals
incongruence_cost_acc_ancova_res <- as.numeric(
  scale(residuals(lm(incongruence_cost_acc_2 ~ 
                         incongruence_cost_acc_1 +
                         switch_cost_acc_1 +
                         average_acc_1 +
                         group,
                     data = explore_tibble_centred))))

# Get predicted values
incongruence_cost_acc_ancova_fit <- fitted(
  lm(incongruence_cost_acc_2 ~ 
       incongruence_cost_acc_1 +
       switch_cost_acc_1 +
       average_acc_1 +
       group,
     data = explore_tibble_centred))

# Final plot
plot(incongruence_cost_acc_ancova_fit, incongruence_cost_acc_ancova_res)

# Estimated marginal means #
# ------------------------------- #

# Check estimated means
incongruence_cost_acc_anc_emm <- emmeans(incongruence_cost_acc_ancova,
                                         ~ group, 
                                        vcov. = vcov(incongruence_cost_acc_ancova))
incongruence_cost_acc_anc_emm

# Normal bar plot (2 covariates)
incgr_cost_acc_anc_tbl <- as_tibble(incongruence_cost_acc_anc_emm)
gg_incongr_anc_acc <- ggplot(incgr_cost_acc_anc_tbl, aes(x = group, y = emmean, fill = group)) +
  geom_col(alpha = 0.9, width = 0.4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.08, linewidth = 0.75) +
  scale_fill_manual(values = pal, guide = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  coord_cartesian(ylim = c(0, NA)) +   
  labs(x = 'Group',
       y = 'Posttest Incongruence Cost Accuracy (p.p)') +
  scale_x_discrete(labels = c('Clean Diet', 'Ketogenic Diet')) +
  theme(legend.position = 'none') +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'incongr_anc_acc.pdf'),
  plot = gg_incongr_anc_acc,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Partial eta squared #
# ----------------------------- #

# Group
incongrer_petasq_group <- 
  incongruence_cost_acc_ancova$statistic['group1']^2 / 
  (incongruence_cost_acc_ancova$statistic['group1']^2 + incongruence_cost_acc_ancova$df.residual)
incongrer_petasq_group

# Pretest
incongrer_petasq_pretest <- 
  incongruence_cost_acc_ancova$statistic['incongruence_cost_acc_1']^2 / 
  (incongruence_cost_acc_ancova$statistic['incongruence_cost_acc_1']^2 + incongruence_cost_acc_ancova$df.residual)
incongrer_petasq_pretest

# Covariate 1
incongrer_petasq_cov1 <- 
  incongruence_cost_acc_ancova$statistic['average_acc_1']^2 / 
  (incongruence_cost_acc_ancova$statistic['average_acc_1']^2 + incongruence_cost_acc_ancova$df.residual)
incongrer_petasq_cov1

# Covariate 2
incongrer_petasq_cov2 <- 
  incongruence_cost_acc_ancova$statistic['switch_cost_acc_1']^2 / 
  (incongruence_cost_acc_ancova$statistic['switch_cost_acc_1']^2 + incongruence_cost_acc_ancova$df.residual)
incongrer_petasq_cov2

# MEDIATION / MODERATION --------------------------

# MEDIATION: Do ketone changes and / or BMI changes mediate the greater 
# improvement in ASRS in the KD group? As ketones varied only in the KD
# group, only this group was used for these analyses. 

# MODERATION: Does the magnitude of the change from baseline depend on the 
# level of baseline?

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
qt_ketones <- explore_tibble %>%
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

qt_ketones_both <- explore_tibble %>%
  left_join(ketones_weighed, by = c('participant_id', 'group', 'cohort')) %>%
  mutate(
    ketones_change = ketones_post - ketones_pre,
    BMI_change = BMI_1 - BMI_2,
    ASRS_change = ASRS_1 - ASRS_2,
    BDI_change = BDI_1 - BDI_2,
    PSQI_change = PSQI_1 - PSQI_2,
    across(c(ends_with('_change'), ends_with('_1'), 'ketones_pre'), 
           ~ as.numeric(scale(.x, scale = FALSE))))
glimpse(qt_ketones)

# Prepare tibble for change analysis #
# -------------------------------------- #

# Adjust tibble with task-switching 
# and questionnaires for change analysis
change_data <- explore_tibble %>%
  mutate(
    ASRS_change = ASRS_1 - ASRS_2,
    BDI_change = BDI_1 - BDI_2,
    PSQI_change = PSQI_1 - PSQI_2,
    switch_cost_rt_change = switch_cost_rt_1 - switch_cost_rt_2,
    incongruence_cost_rt_change = incongruence_cost_rt_1 - incongruence_cost_rt_2,
    switch_cost_acc_change = switch_cost_acc_1 - switch_cost_acc_2,
    incongruence_cost_acc_change = incongruence_cost_acc_1 - incongruence_cost_acc_2)
glimpse(change_data)

# - Visualise data ----------------

# Visualise weighted ketones #
# ================================ #

weighted_ketones <- ggplot(ketones_weighed %>% filter(group == 'KD'),
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
  plot = weighted_ketones,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Visualise ketones change #
# ================================ #

ketones_change <- ggplot(qt_ketones %>% filter(group == 'KD'),
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
  plot = ketones_change,
  device = cairo_pdf,
  width = 6.5, 
  height = 4.5,
  units = 'in'
)

# Visualise BMI change #
# ================================ #

BMI_change <- ggplot(qt_ketones_both,
                         aes(x = group, y = BMI_change, colour = group)) +
  geom_jitter(size = 2.6, width = 0.2, alpha = 0.7) +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = NULL, y = 'BMI Change - Centred') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = limits) +
  theme_apa()

# Save plot
ggsave(
  filename = file.path(plot_directory, 'ketones_change.pdf'),
  plot = ketones_change,
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
                                   outfun = outpro) 
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
fit_y_X2 <- WRS::tshdreg(ASRS_keto_X2_ctrl_BMI, qt_ketones$ASRS_2,
                         xout = TRUE, iter = 10000,
                         outfun = outpro, corfun = pbcor, WARN = FALSE)
y_X2_res <- qt_ketones$ASRS_2 - as.numeric(cbind(1, ASRS_keto_X2_ctrl_BMI) %*% fit_y_X2$coef)

# Residualise BMI_change on the controls
fit_x_X2 <- WRS::tshdreg(ASRS_keto_X2_ctrl_BMI, ASRS_keto_X2[, 'BMI_change'],
                         xout = TRUE, iter = 10000,
                         outfun = outpro, corfun = pbcor, WARN = FALSE)
x_X2_res <- ASRS_keto_X2[, 'BMI_change'] - as.numeric(cbind(1, ASRS_keto_X2_ctrl_BMI) %*% fit_x_X2$coef)

# Pull slope & intercept for BMI_change from the full model
ASRS_X2_slope <- ASRS_keto_2$coef[5]
ASRS_X2_intercept <- ASRS_keto_2$coef[1]

# Plot
tibble(x_res = x_X2_res, y_res = y_X2_res) %>% 
  ggplot(aes(x = x_res, y = y_res)) +
  geom_point(alpha = 0.75, size = 2, colour = pal['KD']) +
  geom_abline(intercept = 0, slope = ASRS_X2_slope,
              colour = 'black', linewidth = 0.9) +
  labs(x = 'BMI change (residuals)',
       y = 'ASRS-2 (residuals)',
       title = 'ASRS Posttest as a Function of BMI Changes') +
  theme_minimal()

# - ASRS sensitivity ---------------

# ROBUST REGRESSION APPROACH (TSTS estimator) #
# ================================== #

# Model 1 #
# ------------------ #

# Define predictors matrix
ASRS_keto_X_tsts <- model.matrix(~ ASRS_1 + ketones_change + ketones_pre, 
                            data = qt_ketones)[ , -1]

# Define regression model
ASRS_keto_tsts <- WRS::tstsreg(ASRS_keto_X, qt_ketones$ASRS_2, xout=TRUE, 
                          iter = 50, sc = pbvar, outfun=outpro, plotit = FALSE)
ASRS_keto_tsts$coef

# Omnibus test
ASRS_keto_omni_tsts <- regtestMC_c(ASRS_keto_X, 
                                 qt_ketones$ASRS_2, 
                                 regfun = tstsreg, nboot = 599, alpha = 0.05, 
                                 plotit = TRUE, xout = TRUE, outfun = outpro) 
ASRS_keto_omni_tsts

# Model 2 #
# ------------------- #

# Define predictors matrix
ASRS_keto_X2_tsts <- model.matrix(~ ASRS_1 + ketones_change + 
                               ketones_pre + BMI_change + BMI_1, 
                             data = qt_ketones)[ , -1]

# Define regression model
ASRS_keto_2_tsts <- WRS::tstsreg(ASRS_keto_X2, qt_ketones$ASRS_2, xout=TRUE, 
                            iter = 50, sc = pbvar, outfun=outpro, plotit = FALSE)
ASRS_keto_2_tsts$coef

# Omnibus test
ASRS_keto_2_omni_tsts <- regtestMC_c(ASRS_keto_X2_tsts, qt_ketones$ASRS_2, 
                                   iter = 50, regfun = tstsreg, nboot = 599, 
                                   alpha = 0.05, plotit = TRUE, xout = TRUE, 
                                   outfun = outpro) 
ASRS_keto_2_omni_tsts

# - Baselines & change (spearman) -------------------

# Select response variables
response_vars_change <- colnames(change_data) %>%
  str_subset('_change')

# Run base-change prediction analysis
basechange_coeffs_all_spear <- run_screen(change_data,
                                 response_vars = response_vars_change,
                                 must_have = c('age', 'BMI_1',
                                               'average_rt_1',
                                               'average_acc_1'),
                                 method_cor = 'spearman',
                                 baseline = TRUE,
                                 pre_vars_reg = '_(1)$')
view(basechange_coeffs_all_spear)

# - Baselines & change (pearson) -------------------

# Select response variables
response_vars_change <- colnames(change_data) %>%
  str_subset('_change')

# Run base-change prediction analysis
basechange_coeffs_all_pear <- run_screen(change_data,
                                    response_vars = response_vars_change,
                                    must_have = c('age', 'BMI_1',
                                                  'average_rt_1',
                                                  'average_acc_1'),
                                    method_cor = 'pearson',
                                    baseline = TRUE,
                                    pre_vars_reg = '_(1)$')
view(basechange_coeffs_all_pear)

# - Baselines & change (pearson + no covs) -------------------

# Select response variables
response_vars_change <- colnames(change_data) %>%
  str_subset('_change')

# Run base-change prediction analysis
basechange_coeffs_all_pear_c0 <- run_screen(change_data,
                                         response_vars = response_vars_change,
                                         must_have = character(),
                                         method_cor = 'pearson',
                                         baseline = TRUE,
                                         pre_vars_reg = '_(1)$')
view(basechange_coeffs_all_pear_c0)

# - Double-check spearman coefficients -------------------------

# Prepare data #
# --------------------- #

# Rank data (both groups)
change_data_rnkd <- change_data %>%
  mutate(
    across(c(age, contains('change'), ends_with('_1')), 
                ~ rank(.x, na.last = 'keep', ties.method = 'average'))
    )

# Rank data (KD)
change_data_rnkd_KD <- change_data %>%
  filter(group == 'KD') %>%
  mutate(
    across(c(age, contains('change'), ends_with('_1')), 
           ~ rank(.x, na.last = 'keep', ties.method = 'average'))
    )

# Rank data (CD)
change_data_rnkd_CD <- change_data %>%
  filter(group == 'CD') %>%
  mutate(
    across(c(age, contains('change'), ends_with('_1')), 
           ~ rank(.x, na.last = 'keep', ties.method = 'average'))
  )

# ASRS #
# =========== #

# Linear model on ranks (default)
ASRS_lm <- lm(
  ASRS_change ~ ASRS_1 + BMI_1 + age + average_rt_1 + average_acc_1, 
  data = change_data_rnkd)
summary(ASRS_lm)
car::Anova(ASRS_lm, type = 'III', white.adjust = TRUE)

# Check for group here
ASRS_lm_2 <- lm(
  ASRS_change ~ group + ASRS_1 + BMI_1 + age + average_rt_1 + average_acc_1, 
  data = change_data_rnkd)
summary(ASRS_lm_2)

# Linear model on ranks (KD)
ASRS_lm_3a <- lm(
  ASRS_change ~ ASRS_1 + BMI_1 + age + average_rt_1 + average_acc_1, 
  data = change_data_rnkd_KD)
summary(ASRS_lm_3a)
car::Anova(ASRS_lm_3a, type = 'III', white.adjust = TRUE)

# Linear model on ranks (CD)
ASRS_lm_3b <- lm(
  ASRS_change ~ ASRS_1 + BMI_1 + age + average_rt_1 + average_acc_1, 
  data = change_data_rnkd_CD)
summary(ASRS_lm_3b)
car::Anova(ASRS_lm_3b, type = 'III', white.adjust = TRUE)

# Scale-free like Spearman #
# ---------------------------- #

ASRS_z <- c('BMI_1','age','average_rt_1','average_acc_1')

# KD #
# **************** #

# Residualise ranked Y and X on ranked controls (keep NA positions)
ASRS_change_fit_KD <- lm(ASRS_change ~ .,
                         data = change_data %>% 
                           filter(group == 'KD') %>%
                           dplyr::select(ASRS_change, all_of(ASRS_z)))

ASRS_base_fit_KD <- lm(ASRS_1 ~ .,
                       data = change_data %>% 
                         filter(group == 'KD') %>%
                         dplyr::select(ASRS_1, all_of(ASRS_z)))

ASRS_change_res_KD <- resid(ASRS_change_fit_KD)
ASRS_base_res_KD  <- resid(ASRS_base_fit_KD)

# Spearman-on-ranks == Pearson on (z-scored) residuals
rho_KD <- cor.test(rank(ASRS_change_res_KD),
              rank(ASRS_base_res_KD),
              method = 'pearson',
              exact = FALSE)
rho_KD

# CD #
# **************** #

# Residualise ranked Y and X on ranked controls (keep NA positions)
ASRS_change_fit_CD <- lm(ASRS_change ~ .,
                         data = change_data %>% 
                           filter(group == 'CD') %>%
                           dplyr::select(ASRS_change, all_of(ASRS_z)))

ASRS_base_fit_CD <- lm(ASRS_1 ~ .,
                       data = change_data %>% 
                         filter(group == 'CD') %>%
                         dplyr::select(ASRS_1, all_of(ASRS_z)))

ASRS_change_res_CD <- resid(ASRS_change_fit_CD)
ASRS_base_res_CD  <- resid(ASRS_base_fit_CD)

# Spearman-on-ranks == Pearson on (z-scored) residuals
rho_CD <- cor.test(rank(ASRS_change_res_CD),
              rank(ASRS_base_res_CD),
              method = 'pearson',
              exact = FALSE)
rho_CD

# - Blomqvist's coeffs without covariates ---------------------

# NOTE: Adjusting partial slopes for measurement error (Chiorelo et al., 2013)

# Specify ICCs
ICCs <- list(
  ASRS = 0.74,
  PSQI = 0.85,
  BDI = 0.85
)

# basechange_covariates not even used below
basechange_covariates <- c('age', 'BMI_1',
                           'average_rt_1',
                           'average_acc_1')

# Both groups combined #
# ========================== #

basechange_coeffs_adj <- blomqvist_k(change_data, 
                                     change_dv = c('BDI_change',
                                                   'PSQI_change',
                                                   'ASRS_change'), 
                                     change_reg = '_change$', 
                                     base_reg = '_1', 
                                     ICCs = ICCs)
basechange_coeffs_adj

# KD #
# ============= #

basechange_coeffs_adj_KD <- blomqvist_k(change_data %>% 
                                          filter(group == 'KD'), 
                                     change_dv = c('BDI_change',
                                                   'PSQI_change',
                                                   'ASRS_change'), 
                                     change_reg = '_change$', 
                                     base_reg = '_1', 
                                     ICCs = ICCs)
basechange_coeffs_adj_KD


# CD #
# ============= #

basechange_coeffs_adj_CD <- blomqvist_k(change_data %>% 
                                          filter(group == 'CD'), 
                                        change_dv = c('BDI_change',
                                                      'PSQI_change',
                                                      'ASRS_change'), 
                                        change_reg = '_change$', 
                                        base_reg = '_1', 
                                        ICCs = ICCs)
basechange_coeffs_adj_CD

# Correct coefficients #
# =========================== #

basechange_coeffs_blom <- basechange_coeffs_all_pear_c0 %>%
  rename(baseline = X) %>%
  left_join(basechange_coeffs_adj %>% 
              dplyr::select(baseline, k), 
            by = 'baseline') %>%
  rename(k_ttl = k) %>%
  left_join(basechange_coeffs_adj_KD %>% 
              dplyr::select(baseline, k), 
            by = 'baseline') %>%
  rename(k_KD = k) %>%
  left_join(basechange_coeffs_adj_CD %>% 
              dplyr::select(baseline, k), 
            by = 'baseline') %>%
  rename(k_CD = k) %>%
  mutate(cor_ttl_adj = (cor_ttl + k_ttl) / (1 - k_ttl),
         cor_KD_adj = (cor_KD + k_KD) / (1 - k_KD),
         cor_CD_adj = (cor_CD + k_CD) / (1 - k_CD)) %>%
  dplyr::select(outcome, baseline, cor_ttl, cor_ttl_adj,
         cor_KD, cor_KD_adj,
         cor_CD, cor_CD_adj,
         everything()) %>%
  mutate(across(where(is.numeric), 
                ~ round(.x, 2)))
view(basechange_coeffs_blom)

# P-VALUE ADJUSTMENT ------------------
# - Effects trimmed ---------------

lab_key <- c(
  Qa = 'group',
  Qb = 'session',
  Qab = 'group:session'
)

# Collect bwtrims
bwtrim_list <- list(
  switch_rt = switch_rt_tranova,
  incongr_rt = incongr_rt_tranova,
  switch_acc = switch_acc_tranova,
  incongr_acc = incongr_acc_tranova,
  ASRS = ASRS_tranova,
  PSQI = PSQI_tranova,
  BDI = BDI_tranova
)

# SESSION ONLY #
# ====================== #

# Apply extraction function
effects_trimmed <- imap(bwtrim_list, ~ {
  df <- tidy_WRS2(.x,
            lab_key,
            p.adjust = FALSE)
  df[['response']] <- .y 
  df
  }) %>%
  list_rbind() %>%
  relocate(response) %>%
  mutate(across(c(statistic, p.value), ~ round(.x, 3))) %>%
  # Select only main effect of session
  filter(term == 'session') %>%
  # Remove ASRS as it had an interaction
  filter(response != 'ASRS') %>%
  # Adjust p-values
  mutate(p.adj = p.adjust(p.value, 'BH'))
effects_trimmed

# SESSION + INT ONLY #
# ====================== #

# Trimmed ANOVA ASRS #
# ----------------------- #

# Apply extraction function
effects_trimmed_sessint <- imap(bwtrim_list, ~ {
  df <- tidy_WRS2(.x,
                  lab_key,
                  p.adjust = FALSE)
  df[['response']] <- .y 
  df
}) %>%
  list_rbind() %>%
  relocate(response) %>%
  # Select main effect of session and interaction
  filter(term == 'session' | (response == 'ASRS' & term == 'group:session')) %>%
  # Remove ASRS Session as there was an ASRS interaction
  filter(!(response == 'ASRS' & term == 'session')) %>%
  # Adjust p-values
  mutate(p.adj = p.adjust(p.value, 'BH')) %>%
  mutate(across(c(statistic, p.value, p.adj), ~ round(.x, 3)))
effects_trimmed_sessint

# Trimmed change scores ASRS #
# -------------------------------- #

# Apply extraction function
effects_trimmed_sesschange <- imap(bwtrim_list, ~ {
  df <- tidy_WRS2(.x,
                  lab_key,
                  p.adjust = FALSE)
  df[['response']] <- .y 
  df
}) %>%
  list_rbind() %>%
  relocate(response) %>%
  # Select only main effect of session
  filter(term == 'session' | (response == 'ASRS' & term == 'group:session')) %>%
  # Remove ASRS Session as there was an ASRS interaction
  filter(!(response == 'ASRS' & term == 'session')) %>%
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
effects_trimmed_sesschange

# - Effects conventional ---------------

# NOTE: 'conventional' = 'regular' = 'classic' 

# Collect bwtrims
conventional_list <- list(
  switch_rt = switch_rt_anova,
  incongr_rt = incongr_rt_anova,
  switch_acc = switch_acc_anova,
  incongr_acc = incongr_acc_anova,
  ASRS = ASRS_anova,
  PSQI = PSQI_anova,
  BDI = BDI_anova
)

# SESSION ONLY #
# ====================== #

# Apply extraction function
effects_conventional <- imap(conventional_list, 
                             ~ {
                               df <- suppressWarnings(tidy(.x$anova_table))
                               df[['response']] <- .y 
                               df
                             }) %>%
  list_rbind() %>%
  relocate(response) %>%
  mutate(across(c(statistic, pes, p.value), ~ round(.x, 3))) %>%
  # Select only main effect of session
  filter(term == 'session') %>%
  # Adjust p-values
  mutate(p.adj = p.adjust(p.value, 'BH'))
effects_conventional