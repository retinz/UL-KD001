# DESCRIPTION #
# ===================== #

# Create objects for reporting the results of UL-KD001. 

# LIBRARIES #
# ===================== #

library(gt)
library(here)
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
source(here('UL-KD001_analysis', 'analysis_helpers.R'))

# LOAD WORKSPACE #
# ===================== #

# Load a generated workspace of the analysis
load('UL-KD001_analysis.RData')

# Create directory if it doesn't exist
if (!dir.exists(here('UL-KD001_analysis', 'tables'))) {
  dir.create(here('UL-KD001_analysis', 'tables'))
}

# AGE ----------------

age_descriptives %>%
  gt() %>%
  fmt_number(
    columns = c(mean_age, sd_age, median_age, IQR_age),
    decimals = 1
  ) %>%
  tab_header(
    title = 'Age Descriptives by Group'
  ) %>%
  cols_label(
    group = 'Group',
    mean_age = 'Mean',
    sd_age = 'SD',
    median_age = 'Median',
    IQR_age = 'IQR'
  ) %>%
tab_options(
  # General look
  table.font.color = 'black',
  table.border.top.color = 'black',
  heading.border.bottom.color = 'black',
  heading.title.font.size = 20,
  table_body.border.bottom.color = 'black',
  column_labels.border.top.color = 'black',
  column_labels.border.bottom.color = 'black',
  column_labels.padding.horizontal = px(20),
  data_row.padding.horizontal = px(20),
  table.font.size = 20,
  data_row.padding = px(10),
  table.width = pct(70),
  table.align = "center",
  table_body.hlines.width = px(0),
  stub.border.width = px(0),
  column_labels.padding = px(10)
)


# BMI ----------------

BMI_descriptives %>%
  gt() %>%
  fmt_number(
    columns = c(group, session, mean_BMI, sd_BMI, median_BMI, IQR_BMI),
    decimals = 1
  ) %>%
  tab_header(
    title = 'BMI Descriptives by Group'
  ) %>%
  cols_label(
    group = 'Group',
    session = 'Session',
    mean_BMI = 'Mean',
    sd_BMI = 'SD',
    median_BMI = 'Median',
    IQR_BMI = 'IQR'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(70),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )


# ASRS Mediation --------------------

ASRS_keto_2_table <- ASRS_keto_2_coeffs$regci %>%
  as_tibble() %>%
  rename(Coefficient = Estimate, SE = S.E., P = `p-value`,
         `P Adj.` = p.adj) %>%
  mutate(`Confidence Interval` = sprintf("[%.2f, %.2f]", ci.low, ci.up)) %>%
  dplyr::select(-ci.low, -ci.up) %>%
  mutate(Term = c('Intercept', 'ASRS Pretest',
                  'Ketone Change', 'Ketone Pretest',
                  'BMI Change', 'BMI Pretest')) %>%
  relocate(Term, Coefficient, SE, `Confidence Interval`) %>%
  mutate(across(c(Coefficient, SE), ~ round(.x, 2)),
         across(c(P, `P Adj.`), ~ ifelse(.x < 0.001, '<0.001', round(.x, 3)))) %>%
  gt() %>%
  tab_header(
    title = 'ASRS Mediation Model'
  ) %>%
  cols_label(
    Term = 'Term',
    Coefficient = 'Coefficient',
    SE = 'SE',
    `Confidence Interval` = 'Confidence Interval',
    P = 'P',
    `P Adj.` = 'P Adj.'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = ASRS_keto_2_table, 
       filename = 'ASRS_keto_2_table.png', path = here('UL-KD001_analysis', 'tables'))

# Screener no. 2 ------------------------

screener_2 <- partial_coeffs_extract_2 %>%
  dplyr::select(-n_controls) %>%
  gt() %>%
  tab_header(
    title = 'Partial Correlation Screener'
  ) %>%
  cols_label(
    outcome = 'Outcome',
    X = 'Predictor',
    cor_ttl = 'Part. Correlation',
    p_ttl = 'Part. Correlation Sig.',
    cor_KD = 'Corr. KD',
    p_KD = 'Corr. KD Sig.',
    cor_CD = 'Corr. CD',
    p_CD = 'Corr. CD Sig.',
    cor_diff = 'Corr. Difference',
    p_het = 'Difference Sig.',
    p_adj = 'Adjusted Sig.'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = screener_2, 
       filename = 'screener_2.png', path = here('UL-KD001_analysis', 'tables'))

# Moderation Spearman --------------------

moderation_spear <- basechange_coeffs_all_spear %>%
  dplyr::select(-n_controls) %>%
  mutate(outcome = case_when(
    outcome == 'ASRS_change' ~ 'ASRS Change',
    outcome == 'BDI_change' ~ 'BDI Change',
    outcome == 'PSQI_change' ~ 'PSQI Change',
    outcome == 'switch_cost_rt_change' ~ 'Switch RT Change',
    outcome == 'incongruence_cost_rt_change' ~ 'Incongr. RT Change',
    outcome == 'switch_cost_acc_change' ~ 'Switch ACC Change',
    outcome == 'incongruence_cost_acc_change' ~ 'Incongr. ACC Change',
    ),
    X = case_when(
      X == 'ASRS_1' ~ 'ASRS Pretest',
      X == 'BDI_1' ~ 'BDI Pretest',
      X == 'PSQI_1' ~ 'PSQI Pretest',
      X == 'switch_cost_rt_1' ~ 'Switch RT Pretest',
      X == 'incongruence_cost_rt_1' ~ 'Incongr. RT Pretest',
      X == 'switch_cost_acc_1' ~ 'Switch ACC Pretest',
      X == 'incongruence_cost_acc_1' ~ 'Incongr. ACC Pretest',
    )) %>%
  gt() %>%
  cols_label(
    outcome = 'Outcome',
    X = 'Predictor',
    cor_ttl = 'Part. Correlation',
    p_ttl = 'Part. Correlation Sig.',
    cor_KD = 'Corr. KD',
    p_KD = 'Corr. KD Sig.',
    cor_CD = 'Corr. CD',
    p_CD = 'Corr. CD Sig.',
    cor_diff = 'Corr. Difference',
    p_het = 'Difference Sig.',
    p_adj = 'Adjusted Sig.'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = moderation_spear, 
       filename = 'moderation_spear.png', path = here('UL-KD001_analysis', 'tables'))

# Mixed-Effects RT ---------------

mixed_eff_rt <- broom.mixed::tidy(taskswitch_mm_10b_rt,
                                  effects = 'fixed',
                                  component = 'cond') %>%
  mutate(Term = c('Intercept', 'Group', 'Session', 'Task Transition',
                  'Congruence', 'Cue Transition', 'Group * Session',
                  'Group * Task Transition', 'Session * Task Transition',
                  'Group * Congruence', 'Session * Congruence', 
                  'Group * Session * Task Transition',
                  'Group * Session * Congruence')) %>%
  dplyr::select(-term, -effect, -component) %>%
  rename(Estimate = estimate, SE = std.error, P = p.value,
         Statistic = statistic) %>%
  mutate(across(c(Estimate, SE, Statistic), 
                ~ round(.x, 2)),
         P = ifelse(P < 0.001, '< 0.001',
                    round(P, 3))) %>%
  relocate(Term) %>% 
  gt() %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = mixed_eff_rt, 
       filename = 'mixed_eff_rt.png', path = here('UL-KD001_analysis', 'tables'))

# Mixed-Effects ACC ---------------

mixed_eff_acc <- broom.mixed::tidy(taskswitch_mm_5_er,
                                  effects = 'fixed',
                                  component = 'cond') %>%
  mutate(Term = c('Intercept', 'Group', 'Session', 'Task Transition',
                  'Congruence', 'Cue Transition', 'Group * Session',
                  'Group * Task Transition', 'Session * Task Transition',
                  'Group * Congruence', 'Session * Congruence', 
                  'Task Transition * Congruence',
                  'Group * Session * Task Transition',
                  'Group * Session * Congruence',
                  'Group * Task Transition * Congruence')) %>%
  dplyr::select(-term, -effect, -component) %>%
  rename(Estimate = estimate, SE = std.error, P = p.value,
         Statistic = statistic) %>%
  mutate(across(c(Estimate, SE, Statistic), 
                ~ round(.x, 2)),
         P = ifelse(P < 0.001, '< 0.001',
                    round(P, 3))) %>%
  relocate(Term) %>% 
  gt() %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = mixed_eff_acc, 
       filename = 'mixed_eff_acc.png', path = here('UL-KD001_analysis', 'tables'))

# Questionnaires descriptives --------------------

quest_desc <- questionnaires_trwin %>%
  dplyr::select(-c(se_tr, se_tr_ws)) %>%
  mutate(mean_tr = round(mean_tr, 2),
         sd_win = round(sd_win, 2)) %>%
  gt() %>%
  cols_label(
    dv = 'Variable',
    group = 'Group',
    session = 'Session',
    n = 'N',
    mean_tr = 'M',
    sd_win = 'SD'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = quest_desc, 
       filename = 'quest_desc.png', path = here('UL-KD001_analysis', 'tables'))

# Task-switching RT descriptives ----------------

taskswitch_rt_desc <- taskswitch_rt_trwin %>%
  dplyr::select(-c(se_tr, se_tr_ws)) %>%
  mutate(mean_tr = round(mean_tr, 2),
         sd_win = round(sd_win, 2)) %>%
  mutate(
    dv = case_when(
      dv == 'average_rt' ~ 'Average RT',
      dv == 'incongruence_cost' ~ 'Incongruence Cost',
      dv == 'switch_cost' ~ 'Switch Cost'
    )) %>%
  gt() %>%
  cols_label(
    dv = 'Variable',
    group = 'Group',
    session = 'Session',
    n = 'N',
    mean_tr = 'M',
    sd_win = 'SD'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = taskswitch_rt_desc, 
       filename = 'taskswitch_rt_desc.png', path = here('UL-KD001_analysis', 'tables'))

# Task-switching ACC descriptives ----------------

taskswitch_acc_desc <- taskswitch_er_trwin %>%
  dplyr::select(-c(se_tr, se_tr_ws)) %>%
  mutate(mean_tr = round(mean_tr, 2),
         sd_win = round(sd_win, 2)) %>%
  mutate(
    dv = case_when(
      dv == 'average_acc' ~ 'Average Accuracy',
      dv == 'incongruence_cost' ~ 'Incongruence Cost',
      dv == 'switch_cost' ~ 'Switch Cost'
    )) %>%
  gt() %>%
  cols_label(
    dv = 'Variable',
    group = 'Group',
    session = 'Session',
    n = 'N',
    mean_tr = 'M',
    sd_win = 'SD'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(20),
    data_row.padding.horizontal = px(20),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(90),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = taskswitch_acc_desc, 
       filename = 'taskswitch_acc_desc.png', path = here('UL-KD001_analysis', 'tables'))
# Conventional ANOVAs table -------------

conventional_anovas <- effects_conventional_all %>%
  mutate(across(c('MSE', 'statistic', 'pes'), 
                ~ round(.x, 2)),
         p.adj = round(p.adj, 3)) %>%
  arrange(response, term) %>%
  gt() %>%
  cols_label(
    response = 'Response',
    term = 'Effect',
    num.Df = 'DF Numer.',
    den.Df = 'DF Denom.',
    MSE = 'MSE',
    statistic = 'F',
    pes = 'η',
    p.value = 'P-Value',
    p.adj = 'Adjusted P-Value'
  ) %>%
  tab_options(
    # General look
    table.font.color = 'black',
    table.border.top.color = 'black',
    heading.border.bottom.color = 'black',
    heading.title.font.size = 20,
    table_body.border.bottom.color = 'black',
    column_labels.border.top.color = 'black',
    column_labels.border.bottom.color = 'black',
    column_labels.padding.horizontal = px(15),
    data_row.padding.horizontal = px(15),
    table.font.size = 20,
    data_row.padding = px(10),
    table.width = pct(80),
    table.align = "center",
    table_body.hlines.width = px(0),
    stub.border.width = px(0),
    column_labels.padding = px(10)
  )

gtsave(data = conventional_anovas, 
       filename = 'conventional_anovas.png', path = here('UL-KD001_analysis', 'tables'))