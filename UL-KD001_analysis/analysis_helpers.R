# DESCRIPTION 
# ====================== #

# Helper functions for standard behavioural sciences analyses. 

# LIBRARIES #
# ====================== #

library(ggnewscale)
library(qqplotr)
library(ggplot2)
library(tidyverse)

# CREATNG LISTS ----------------

#' Split a numeric column into a named list over arbitrary factors
#'
#' @param data         A data.frame
#' @param value_col    String name of the numeric column to split out
#' @param factor_cols  Character vector of one or more factor column names;
#'                     result keys will be in the order factor_cols[1]…factor_cols[n]
#' @param drop         Passed to split(); if TRUE, empty factor combos are dropped
#' @return A named list of vectors, in the order of factor_cols[1]…factor_cols[n]
#' @examples
#' 
#' # 2-factor use:
#' HF_list <- build_split_list(ecg_light, 'HRV_HF', c('condition', 'order'))
#'
#' # 3-factor example:
#' # assume ecg has columns 'condition', 'order', 'session'
#' HF3_list <- build_split_list(ecg_light, 'HRV_HF',
#'                              c('condition','order','session'))

build_split_list <- function(data, value_col, factor_cols, drop = TRUE) {
  
  # Sanity
  if (!value_col %in% names(data))
    stop('Column to split ‘', value_col, '’ not found.')
  if (!all(factor_cols %in% names(data)))
    stop('Some factor columns not in data: ',
         paste(setdiff(factor_cols, names(data)), collapse = ', '))
  
  # 1) Do the split
  grouping <- lapply(factor_cols, function(f) data[[f]])
  vecs <- split(data[[value_col]], grouping, drop = drop)
  
  # 2) Get each factor’s levels (always via factor() so non-factors are handled)
  levs <- lapply(factor_cols, function(f) levels(factor(data[[f]])))
  
  # 3) Build the full grid—but pass the factors *in reverse* so that
  #    expand.grid’s last argument (factor_cols[n]) varies fastest
  rev_cols <- rev(factor_cols)
  args <- setNames(lapply(rev_cols, function(f)
    levs[[match(f, factor_cols)]]),
    rev_cols)
  grid_rev <- do.call(expand.grid,
                      c(args, stringsAsFactors = FALSE))
  
  # 4) Reorder the grid back into the original factor_cols order
  grid <- grid_rev[, factor_cols, drop = FALSE]
  
  # 5) Build the same “f1.f2.….fn” keys that split() used internally
  keys <- do.call(paste, c(grid, sep = '.'))
  
  # 6) Re-index and return
  vecs[keys]
}

# TIDY ANOVA OUTPUT ------------------------
# - WRS --------------

# A function to make a generic WRS ANOVA output human-readable
tidy_WRS <- function(anova_vec,
                       lab_key,
                       p.adjust = FALSE,
                       boot = FALSE) {
  # Validation
  if (missing(lab_key) || !is.character(lab_key) || is.null(names(lab_key))) {
    stop('`lab_key` must be a named character vector (e.g. c(Qa = `group`, …)).')
  }
  
  # Helper: apply p-adjust if requested
  maybe_adjust <- function(tbl) {
    if (!identical(p.adjust, FALSE)) {
      tbl <- dplyr::mutate(tbl,
                           p.value = stats::p.adjust(p.value, method = p.adjust))
    }
    tbl
  }
  
  if (boot) {
    # Bootstrap version
    # - Expected names like 'p.value.A', 'p.value.AB', …
    out <- tibble::enframe(anova_vec, name = 'raw', value = 'p.value') %>%
      dplyr::mutate(
        p.value = purrr::map_dbl(p.value, as.numeric),
        
        # Strip 'p.value.'  →  'A', 'AB', …
        boot_id = stringr::str_remove(raw, '^p\\.value\\.'),
        
        # Convert to canonical 'Qa', 'Qab', … expected by lab_key
        term_id = paste0('Q', stringr::str_to_lower(boot_id))
      ) %>%
      dplyr::mutate(term = unname(lab_key[term_id])) %>%
      dplyr::select(term, p.value) %>%
      maybe_adjust() %>%
      dplyr::mutate(across(p.value, ~ round(.x, 3)))
    
  } else {
    # Regular version
    out <- tibble::enframe(anova_vec, name = 'raw', value = 'value') %>%
      dplyr::mutate(
        value = purrr::map_dbl(value, as.numeric),
        term_id = stringr::str_remove(raw, '\\.p.value$'),
        metric = dplyr::if_else(stringr::str_detect(raw, '\\.p.value$'),
                                 'p.value', 'statistic')
      ) %>%
      dplyr::select(-raw) %>%                                   
      tidyr::pivot_wider(id_cols = term_id,
                         names_from = metric,
                         values_from = value) %>%
      dplyr::mutate(term = unname(lab_key[term_id])) %>%
      dplyr::select(-term_id) %>%
      dplyr::relocate(term) %>%
      maybe_adjust() %>%
      dplyr::mutate(across(c(statistic, p.value), ~ round(.x, 3)))
  }
  
  out
}

# - WRS2 ---------------

# A function to create a tidy ANOVA tibble from WRS2 objects
tidy_WRS2 <- function(WRS2_object,
                      lab_key,
                      p.adjust = FALSE,
                      boot = FALSE) { # boot not implemented
  
  # ---- Validation ----
  
  if (missing(lab_key) || !is.character(lab_key) || is.null(names(lab_key))) {
    stop('`lab_key` must be a named character vector, e.g. c(Qa = \'group\', Qb = \'session\', Qab = \'group:session\').')
  }
  # Allow only known class for now
  if (!inherits(WRS2_object, 'bwtrim')) {
    stop('Currently implemented for objects of class \'bwtrim\'.')
  }
  
  # ---- Helpers ----
  
  maybe_adjust <- function(tbl) {
    if (!identical(p.adjust, FALSE)) {
      tbl <- dplyr::mutate(tbl, p.value = stats::p.adjust(p.value,
                                                          method = p.adjust))
    }
    tbl
  }
  get_el <- function(nm) {
    x <- WRS2_object[[nm]]
    if (is.null(x)) return(NA_real_)
    x
  }
  
  # ---- Build key table: Qa/Qb/Qab -> A/B/AB and term labels ----
  
  key_tbl <- tibble::tibble(
    stat_name = names(lab_key),
    term = unname(lab_key),
    letter = toupper(sub('^Q', '', names(lab_key)))  # 'Qa' -> 'A', 'Qab' -> 'AB'
  )
  
  # ---- Extract pieces ----
  
  out <- key_tbl %>%
    dplyr::mutate(
      statistic = purrr::map_dbl(stat_name, ~ as.numeric(get_el(.x))),
      df_vec = purrr::map(letter, ~ {
        v <- WRS2_object[[paste0(.x, '.df')]]
        if (is.null(v)) return(rep(NA_real_, 2))
        as.numeric(v)
      }),
      p.value = purrr::map_dbl(letter, ~ as.numeric(get_el(paste0(.x, '.p.value'))))
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      df1 = df_vec[1],
      df2 = df_vec[2]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(term, statistic, df1, df2, p.value) %>%
    maybe_adjust()
  
  out
}

# MIXED-MODEL DIAGNOSTICS ------------------

# For lme4::lmer using package performance #
# ------------------------------------------------ #

diagnose_lmer <- function(model) {
  # Helper that does “run a check + grab its plot”
  run_check <- function(fun, ...) {
    obj  <- fun(model, ...)
    plt  <- plot(obj) 
    list(result = obj, plot = plt)
  }
  
  # Name = run_check(<function>, ...extra args if needed)
  out <- list(
    collinearity = run_check(check_collinearity),
    outliers = run_check(check_outliers),
    homoscedasticity = run_check(check_heteroscedasticity),
    homoscedasticity_bin = run_check(binned_residuals),
    residuals_norm = run_check(check_normality),
    ran_eff_norm = run_check(check_normality, effects = 'random')
  )
  
  # Make it easy to write print / plot methods later
  class(out) <- c('diagnose_mixed_model', class(out)) 
  out
}

# Custom print and plot function
print.diagnose_mixed_model <- function(x, 
                                show_results = TRUE, 
                                show_plots = FALSE,
                                ...) {
  n <- length(x)
  cat('<diagnose_mixed_model>  —  containing', n, 'checks\n')
  cat('Names:', paste(names(x), collapse = ', '), '\n\n')
  
  if (show_results) {
    imap(x, ~ {
      cat('>>', .y, '<<\n')
      print(.x$result)
      cat('\n')
    })
  }
  
  if (show_plots) {
    walk(x, ~ print(.x$plot))
  }
  
  invisible(x)
}

# For any glm(m) model using package DHARMa #
# ------------------------------------------------ #

#' Diagnose a (G)LMM with DHARMa
#'
#' @param model A fitted glmmTMB / lme4 / glm / gam … object
#' @param sim_val Passed to `set_simcodes()`; use 'fix' for conditional sims
#' @param nsim Number of simulated residual draws
#' @param quiet Suppress console output?
#' @return (Invisibly) a named list with the DHARMa objects
#' @examples
#' fit  <- lme4::lmer(y ~ x + (1|id), data = dat)
#' diag <- diagnose_model(fit, lme4 = TRUE)


diagnose_model <- function(model,
                           sim_val = 'fix',
                           nsim = 250, # DHARMa default
                           quiet = FALSE,
                           lme4 = FALSE) {

  
  # Simulate residuals -------------
  
  sim_model <- model
  
  if (lme4 == FALSE) {
    glmmTMB::set_simcodes(sim_model, val = sim_val)
  }
  
  res <- DHARMa::simulateResiduals(sim_model, plot = FALSE, n = nsim)
  
  # Dispersion test -------------
  
  qt <- DHARMa::testQuantiles(res)
  if (!quiet) print(qt)
  
  # Dispersion test -------------
  
  disp <- DHARMa::testDispersion(res)
  if (!quiet) print(disp)
  
  # Heteroscedasticity -------------
  
  DHARMa::plotResiduals(res,
                        absoluteDeviation = TRUE,
                        quantreg = TRUE)
  
  # Distribution of residuals
  DHARMa::plotQQunif(res)
  
  # Uniformity (KS) test -------------
  
  unif <- DHARMa::testUniformity(res, plot = FALSE)
  if (!quiet) print(unif)
  
  invisible(list(
    simulate_residuals = res,
    quantile_test = qt,
    dispersion_test = disp,
    uniformity_test = unif
  ))
}


# VISUALISE DISTRIBUTIONS ------------------
# - Subject-level -------------------

inspect_subject_dist <- function(
    data,
    pid_col,
    col_dist,
    group_col,
    session_col,
    wrap = NULL,   # NULL, character vector, or formula
    n = 1,
    pal = pal,
    robust = FALSE,  # or 'MAD'
    mad_mult = 2,
    seed = NULL,
    binwidth = NULL, # allow control
    bins = 30
) {
  if (!is.null(seed)) set.seed(seed)
  
  pids <- sample(unique(data[[pid_col]]), size = n)
  
  make_facet <- function(df, wrap) {
    if (is.null(wrap) || length(wrap) == 0) return(NULL)
    if (inherits(wrap, 'formula')) return(wrap)
    if (is.character(wrap)) {
      missing <- setdiff(wrap, names(df))
      if (length(missing) > 0) {
        warning('wrap columns not in data: ', paste(missing, collapse = ', '))
        wrap <- setdiff(wrap, missing)
      }
      if (length(wrap) == 0) return(NULL)
      return(stats::as.formula(paste('~', paste(wrap, collapse = ' + '))))
    }
    warning('wrap must be NULL, a character vector, or a formula; ignoring.')
    NULL
  }
  
  plots <- lapply(pids, function(pid) {
    
    # Pick a session that actually exists for THIS pid (fix)
    df_pid <- data %>% dplyr::filter(.data[[pid_col]] == pid)
    sess <- sample(unique(df_pid[[session_col]]), size = 1)
    df <- df_pid %>% dplyr::filter(.data[[session_col]] == sess)
    
    if (nrow(df) == 0L) {
      return(
        ggplot() +
          annotate('text', x = 0, y = 0, label = paste('No data for', pid), size = 5) +
          theme_void()
      )
    }
    
    # Robust cut points
    if (identical(robust, 'MAD')) {
      med <- stats::median(df[[col_dist]], na.rm = TRUE)
      margin <- stats::mad(df[[col_dist]], na.rm = TRUE) * mad_mult
      cuts <- c(med - margin, med + margin)
    } else {
      cuts <- range(df[[col_dist]], na.rm = TRUE)
    }
    
    # Trimmed mean (falls back to plain mean if trimming removes all rows)
    mean_rt <- df %>%
      dplyr::filter(.data[[col_dist]] > cuts[1], .data[[col_dist]] < cuts[2]) %>%
      dplyr::summarise(m = mean(.data[[col_dist]], na.rm = TRUE)) %>%
      dplyr::pull(m)
    if (length(mean_rt) == 0L || is.na(mean_rt)) {
      mean_rt <- mean(df[[col_dist]], na.rm = TRUE)
    }
    
    # Pick a sensible binwidth if not provided
    if (is.null(binwidth)) {
      rng <- diff(range(df[[col_dist]], na.rm = TRUE))
      binwidth <- if (is.finite(rng) && rng > 0) rng / bins else 0.5
    }
    
    p <- ggplot(df, aes(x = .data[[col_dist]])) +
      # Put histogram on density scale to match the label (fix)
      geom_histogram(aes(y = after_stat(density), fill = .data[[group_col]]),
                     binwidth = binwidth, colour = 'white', alpha = 0.5) +
      # Show raw values as a jittered baseline strip on a compatible y (fix)
      geom_jitter(aes(x = .data[[col_dist]], y = 0),
                  colour = 'black', fill = 'black',
                  height = 0.002, width = 0, alpha = 0.65, size = 1.4) +
      scale_fill_manual(values = pal, name = group_col) +
      labs(
        title = paste('Participant', pid, '| session', sess),
        x = col_dist,
        y = 'Density',
        fill = group_col
      ) +
      theme_minimal()
    
    facet_spec <- make_facet(df, wrap)
    if (!is.null(facet_spec)) p <- p + facet_wrap(facet_spec)
    
    if (identical(robust, 'MAD')) {
      p <- p +
        geom_vline(xintercept = cuts, colour = 'red', linetype = 'dashed', linewidth = 0.7)
    }
    if (!is.na(mean_rt)) {
      p <- p + geom_vline(xintercept = mean_rt, colour = 'black', linetype = 'solid', linewidth = 1)
    }
    
    p
  })
  
  names(plots) <- pids
  plots
}


# ROBUSTIFICATION ------------------
# - Robust estimation -----------------

trim_winsorise <- function(data, dv_colnames, within, 
                           between = NULL, id, tr = 0.2) {
  
  # Input validation
  stopifnot(
    is.data.frame(data),
    is.character(dv_colnames), all(dv_colnames %in% names(data)),
    is.character(within), length(within) == 1, within %in% names(data),
    is.null(between) || (is.character(between) && all(between %in% names(data))),
    is.character(id), length(id) == 1, id %in% names(data),
    is.numeric(tr), tr >= 0, tr < 0.5
  )
  
  # Pivot to long format with all DVs in one column
  long_data <- data %>%
    pivot_longer(cols = all_of(dv_colnames), 
                 names_to = "dv", 
                 values_to = "value") %>%
    filter(!is.na(value))
  
  # Number of within-subject levels for Morey correction
  M <- n_distinct(long_data[[within]])
  
  # Calculate normalized values for within-subject SE
  normalized_data <- long_data %>%
    group_by(dv, .data[[id]], across(all_of(between))) %>%
    mutate(participant_mean = mean(value, na.rm = TRUE)) %>%
    group_by(dv, across(all_of(between))) %>%
    mutate(
      grand_mean = mean(value, na.rm = TRUE),
      normalized_value = value - participant_mean + grand_mean
    ) %>%
    ungroup()
  
  # Calculate all statistics in one go
  result <- normalized_data %>%
    group_by(dv, across(all_of(c(between, within)))) %>%
    summarise(
      n = n(),
      mean_tr = mean(value, trim = tr),
      sd_win = sqrt(WRS2::winvar(value, tr = tr)),
      se_tr = WRS2::trimse(value, tr = tr),
      se_tr_ws = if(M > 1) {
        WRS2::trimse(normalized_value, tr = tr) * sqrt(M / (M - 1))
      } else {
        NA_real_
      },
      .groups = 'drop'
    )
  
  return(result)
}

# - Robust central tendency (boxplots) --------------

# Computes between-subject SE of the trimmed mean
mean_se_tr <- function(z, tr = 0.2) {
  # Remove NAs
  z_complete <- z[!is.na(z)]
  if (length(z_complete) == 0) {
    return(data.frame(y = NA, ymin = NA, ymax = NA))
  }
  
  # Calculate trimmed mean and SE using WRS2 functions
  tr_mean <- mean(z_complete, trim = tr)
  se <- WRS2::trimse(z_complete, tr = tr)
  
  # Return a data frame with y, ymin, and ymax
  data.frame(
    y = tr_mean,
    ymin = tr_mean - se,
    ymax = tr_mean + se
  )
}

# BAYESIAN ANALYSIS ------------------
# - Sensitivity analysis ------------------

# General-purpose plotter for BayesFactor ratios (e.g., interaction vs additive)
bf_ratio_plot <- function(...,
                          labels = NULL,
                          log10_scale = TRUE,
                          thresholds = c(1/10, 1/3, 1, 3, 10),
                          annotate = TRUE,
                          show_bands = TRUE) {
  
  stopifnot(length(list(...)) > 0)
  objs <- list(...)
  
  # Names: prefer names(objs), else user-supplied labels, else generic
  nm <- names(objs)
  if (is.null(nm) || any(nm == '')) {
    if (!is.null(labels)) {
      stopifnot(length(labels) == length(objs))
      nm <- labels
    } else {
      nm <- paste0('model_', seq_along(objs))
    }
  }
  
  # Extract BF and relative error from each BFBayesFactor ratio
  get_row <- function(obj, name) {
    info <- extractBF(obj)[1, c('bf', 'error')]
    tibble(
      label = name,
      BF10 = as.numeric(info$bf),
      rel_error = as.numeric(info$error)
    )
  }
  df <- map2_dfr(objs, nm, get_row)
  
  # Preserve the input order in the x-axis
  df$label <- factor(df$label, levels = nm)
  
  # Create y values based on scale choice
  if (log10_scale) {
    df$y_plot <- log10(df$BF10)
    y_breaks <- log10(thresholds)
    y_labels <- thresholds
  } else {
    df$y_plot <- df$BF10
    y_breaks <- thresholds
    y_labels <- thresholds
  }
  
  # Calculate band positions (needed for y-axis limits even when show_bands = FALSE)
  if (log10_scale) {
    yvals <- df$y_plot
    outer_low <- floor(min(yvals, na.rm = TRUE)) - 1
    outer_high <- ceiling(max(yvals, na.rm = TRUE)) + 1.5
    breaks <- c(outer_low, log10(1/10), log10(1/3), 0,
                log10(3), log10(10), outer_high)
  } else {
    yvals <- df$y_plot
    outer_low <- max(min(yvals, na.rm = TRUE) / 10, .Machine$double.eps)
    outer_high <- max(yvals, na.rm = TRUE) * 30
    breaks <- c(outer_low, 1/10, 1/3, 1, 3, 10, outer_high)
  }
  
  # Plot - use y_plot for the y aesthetic
  p <- ggplot(df, aes(x = label, y = y_plot))
  
  if (show_bands) {
    band_labels <- c('strong H0', 'moderate H0', 'weak H0',
                     'weak H1', 'moderate H1', 'strong H1')
    bands <- tibble(
      ymin = head(breaks, -1),
      ymax = tail(breaks, -1),
      lab = band_labels
    )
    
    p <- p +
      annotate(
        'text',
        x = Inf,
        y = (bands$ymin + bands$ymax) / 2,
        label = bands$lab,
        hjust = 1.02, vjust = 0.5, size = 3
      ) +
      geom_hline(yintercept = if(log10_scale) log10(thresholds) else thresholds, 
                 linetype = 'dotted') +
      geom_hline(yintercept = if(log10_scale) 0 else 1, linetype = 'dashed')
  }
  
  p <- p +
    geom_point(size = 3) +
    geom_line(aes(group = 1), linewidth = 0.6) +
    labs(
      x = 'Prior',
      y = 'Bayes factor (BF10)',
      title = 'Evidence for target effect across priors'
    ) +
    theme_apa()
  
  # Improved scale_y_continuous with scientific formatting for large numbers
  if (log10_scale) {
    if (show_bands) {
      p <- p + 
        scale_y_continuous(
          breaks = y_breaks,
          labels = function(x) {
            values <- 10^x
            # Use scientific notation for very large or very small numbers
            ifelse(abs(x) >= 3, 
                   scales::scientific_format()(values),
                   round(values, 1))
          }
        )
    } else {
      p <- p + 
        scale_y_continuous(
          labels = function(x) {
            values <- 10^x
            # Use scientific notation for very large or very small numbers
            ifelse(abs(x) >= 3, 
                   scales::scientific_format()(values),
                   round(values, 1))
          }
        )
    }
  } else {
    if (show_bands) {
      p <- p + 
        scale_y_continuous(
          breaks = y_breaks,
          labels = function(x) {
            # Use scientific notation for very large or very small numbers
            ifelse(abs(x) >= 1000 | (x != 0 & abs(x) <= 0.001), 
                   scales::scientific_format()(x),
                   round(x, 1))
          }
        )
    } else {
      p <- p + 
        scale_y_continuous(
          labels = function(x) {
            # Use scientific notation for very large or very small numbers
            ifelse(abs(x) >= 1000 | (x != 0 & abs(x) <= 0.001), 
                   scales::scientific_format()(x),
                   round(x, 1))
          }
        )
    }
  }
  
  if (annotate) {
    # Format annotation text - use scientific notation for large numbers
    ann_text <- sprintf('%s (±%.1f%%)', 
                        ifelse(df$BF10 >= 1000 | (df$BF10 != 0 & df$BF10 <= 0.001),
                               scales::scientific_format()(df$BF10),
                               round(df$BF10, 2)),
                        round(100 * df$rel_error, 2))
    p <- p + geom_text(aes(label = ann_text), vjust = -0.8, size = 3.2)
  }
  
  list(data = df, plot = p)
}

# - Posterior draws --------------------

#' Posterior contrast (2x2) with BayesFactor::lmBF + histogram
#'
#' @param data A data frame.
#' @param dv Outcome column name (character).
#' @param factor_a Name of factor A (character).
#' @param factor_b Name of factor B (character).
#' @param id_var Participant id variable (character), coerced to factor.
#' @param a_levels Length-2 character vector: levels of factor A in order c(A1, A2).
#' @param b_levels Length-2 character vector: levels of factor B in order c(B1, B2).
#' @param contrast One of c('interaction', 'main_a', 'main_b').
#' @param main_effect_source One of c('auto', 'cells', 'explicit').
#'        'explicit' uses columns '<factor>-<level>' (e.g., 'group-KD').
#'        'cells' averages the four cell columns '<A>:<B>-<a>.&.<b>'.
#'        'auto' prefers 'explicit' if available, else falls back to 'cells'.
#' @param iterations Posterior draws.
#' @param rscale_fixed Fixed-effects prior scale.
#' @param rscale_random Random-effects prior scale.
#' @param posterior_exclude_regex Regex to drop columns (default drops participant effects).
#' @param ci Credible interval level.
#' @param binwidth Histogram binwidth.
#' @param fill Histogram fill color.
#' @param outline Histogram outline color.
#' @param ci_color Color for CI lines.
#' @param mean_color Color for mean line.
#' @param title Plot title.
#'
#' @return list(model, draws, contrast_draws, summary, plot)

bf_contrast_posterior <- function(data,
                                  dv,
                                  factor_a,
                                  factor_b,
                                  id_var,
                                  a_levels = NULL,
                                  b_levels = NULL,
                                  contrast = c('interaction', 'main_a',
                                               'main_b'),
                                  main_effect_source = c('auto', 'cells',
                                                         'explicit'),
                                  iterations = 100000,
                                  rscale_fixed = 'medium',
                                  rscale_random = 'nuisance',
                                  posterior_exclude_regex = '^participant_id($|-)',
                                  ci = 0.95,
                                  contrast_setting = c('contr.sum', 'contr.poly'),
                                  binwidth = 0.1,
                                  fill = 'steelblue',
                                  outline = 'black',
                                  ci_color = 'red4',
                                  mean_color = 'black',
                                  title = 'Posterior contrast') {
  contrast <- match.arg(contrast)
  main_effect_source <- match.arg(main_effect_source)
  # Set the kind of contrast to be performed in the lmBF
  options(contrasts = contrast_setting)
  
  df <- data %>%
    mutate(
      !!id_var := factor(.data[[id_var]]),
      !!factor_a := factor(.data[[factor_a]]),
      !!factor_b := factor(.data[[factor_b]])
    )
  
  a_lvls <- if (is.null(a_levels)) levels(df[[factor_a]])[1:2] else a_levels
  b_lvls <- if (is.null(b_levels)) levels(df[[factor_b]])[1:2] else b_levels
  stopifnot(length(a_lvls) == 2, length(b_lvls) == 2)
  
  form <- stats::as.formula(paste(dv, '~', 
                                  paste0(factor_a, ' * ',
                                         factor_b, ' + ',
                                         id_var)))
  bf_mod <- BayesFactor::lmBF(
    formula = form,
    data = df %>% as.data.frame(),
    whichRandom = id_var,
    rscaleFixed = rscale_fixed,
    rscaleRandom = rscale_random,
    iterations = iterations
  )
  
  draws <- BayesFactor::posterior(bf_mod, iterations = iterations)
  
  # Remove unnecessary columns --> smaller object in workspace
  if (!is.null(posterior_exclude_regex)) {
    keep <- !grepl(posterior_exclude_regex, colnames(draws))
    draws <- draws[, keep, drop = FALSE]
  }
  
  # BayesFactor::posterior adopts the following naming for columns
  cell_col <- function(a, b) paste0(factor_a, ':', factor_b, '-', a, '.&.', b)
  main_col <- function(fctr, lvl) paste0(fctr, '-', lvl)
  
  A1B1 <- cell_col(a_lvls[1], b_lvls[1])
  A1B2 <- cell_col(a_lvls[1], b_lvls[2])
  A2B1 <- cell_col(a_lvls[2], b_lvls[1])
  A2B2 <- cell_col(a_lvls[2], b_lvls[2])
  
  # Helpers to pick main-effect source
  have_explicit_a <- all(main_col(factor_a, a_lvls) %in% colnames(draws))
  have_explicit_b <- all(main_col(factor_b, b_lvls) %in% colnames(draws))
  
  use_explicit_a <- switch(main_effect_source,
                           'explicit' = TRUE,
                           'cells' = FALSE,
                           'auto' = have_explicit_a)
  
  use_explicit_b <- switch(main_effect_source,
                           'explicit' = TRUE,
                           'cells' = FALSE,
                           'auto' = have_explicit_b)
  
  # Build the requested contrast
  if (contrast == 'interaction') {
    needed <- c(A1B1, A1B2, A2B1, A2B2)
    miss <- setdiff(needed, colnames(draws))
    if (length(miss) > 0) stop('Missing cell columns: ', 
                               paste(miss, collapse = ', '))
    contrast_draws <- (draws[, A1B1] - draws[, A1B2]) - (draws[, A2B1] - draws[, A2B2])
    xlab <- 'Posterior Draws of the Interaction (Difference in Differences)'
    subtitle <- paste0('A: ', a_lvls[1], ' vs ', a_lvls[2], ' | B: ',
                       b_lvls[1], ' → ', b_lvls[2])
    
  } else if (contrast == 'main_a') {
    if (use_explicit_a) {
      cols <- main_col(factor_a, a_lvls)
      miss <- setdiff(cols, colnames(draws))
      if (length(miss) > 0) stop('Missing explicit main-effect columns for A: ', paste(miss, collapse = ', '))
      contrast_draws <- draws[, cols[1]] - draws[, cols[2]]
    } else {
      needed <- c(A1B1, A1B2, A2B1, A2B2)
      miss <- setdiff(needed, colnames(draws))
      if (length(miss) > 0) stop('Missing cell columns: ',
                                 paste(miss, collapse = ', '))
      contrast_draws <- ((draws[, A1B1] + draws[, A1B2]) / 2) - ((draws[, A2B1] + draws[, A2B2]) / 2)
    }
    xlab <- paste0('Posterior Draws of the Main Effect of Group')
    subtitle <- if (use_explicit_a) 'Using explicit main-effect columns' else 'Marginal over B (cell average)'
    
  } else { # main_b
    if (use_explicit_b) {
      cols <- main_col(factor_b, b_lvls)
      miss <- setdiff(cols, colnames(draws))
      if (length(miss) > 0) stop('Missing explicit main-effect columns for B: ', paste(miss, collapse = ', '))
      contrast_draws <- draws[, cols[2]] - draws[, cols[1]]  # B2 − B1
    } else {
      needed <- c(A1B1, A1B2, A2B1, A2B2)
      miss <- setdiff(needed, colnames(draws))
      if (length(miss) > 0) stop('Missing cell columns: ', paste(miss, collapse = ', '))
      contrast_draws <- ((draws[, A1B2] + draws[, A2B2]) / 2) - ((draws[, A1B1] + draws[, A2B1]) / 2)
    }
    xlab <- paste0('Posterior Draws of the Main Effect of Session')
    subtitle <- if (use_explicit_b) 'Using explicit main-effect columns' else 'Marginal over A (cell average)'
  }
  
  alpha <- (1 - ci) / 2
  contr_tbl <- tibble(contrast = as.numeric(contrast_draws))
  contr_sum <- contr_tbl %>%
    summarise(
      mean = mean(contrast),
      median = median(contrast),
      q_lower = quantile(contrast, alpha),
      q_upper = quantile(contrast, 1 - alpha)
    ) %>%
    ungroup()
  
  p <- contr_tbl %>%
    ggplot(aes(contrast)) +
    geom_histogram(alpha = 0.7, binwidth = binwidth, fill = fill, colour = outline) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 1) +
    geom_vline(xintercept = c(contr_sum$q_lower, contr_sum$q_upper),
               linetype = 'dotted', colour = ci_color, linewidth = 0.9) +
    geom_vline(xintercept = contr_sum$mean, linetype = 'solid', colour = mean_color, linewidth = 1) +
    labs(x = xlab, y = 'Density', title = title, subtitle = subtitle) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_apa()
  
  list(
    model = bf_mod,
    draws = draws,
    contrast_draws = contrast_draws,
    summary = contr_sum,
    plot = p
  )
}

# COVARIATE DETECTION ---------------------
# - Screener ----------------

# Get baseline variables
baseline_of <- function(y2) str_replace(y2, '_2$', '_1')

# Compute partial Spearman for a single pair (y2, x), 
# conditioning ONLY on 'controls'
pcorr_one <- function(data, y2, x, controls = character(),
                      method_cor) {
  cols <- unique(c(y2, x, controls))
  dd <- data %>% dplyr::select(all_of(cols)) %>% tidyr::drop_na()
  if(nrow(dd) < 5) {
    return(tibble(cor = NA_real_, p = NA_real_, n = nrow(dd)))
  }
  res <- correlation(dd,
                     method = method_cor,
                     partial = length(controls) > 0,
                     p_adjust = 'none') %>%
    as_tibble() %>%
    filter((Parameter1 == y2 & Parameter2 == x) | (Parameter1 == x & Parameter2 == y2)) %>%
    slice(1) %>%
    mutate(cor = .data[[names(.)[3]]], p = p, n = n_Obs, .keep = 'none')
  if(nrow(res) == 0) res <- tibble(cor = NA_real_, p = NA_real_, n = nrow(dd))
  res
}

# Fisher z test for KD vs CD difference in (partial) correlations
fisher_diff_p <- function(r1, n1, r2, n2, n_controls) {
  if (any(is.na(c(r1, n1, r2, n2)))) return(NA_real_)
  df1 <- n1 - n_controls - 3
  df2 <- n2 - n_controls - 3
  if (min(df1, df2) <= 1) return(NA_real_)
  z1 <- atanh(pmin(pmax(r1, -0.999999), 0.999999))
  z2 <- atanh(pmin(pmax(r2, -0.999999), 0.999999))
  se <- sqrt(1 / df1 + 1 / df2)
  2 * pnorm(abs((z1 - z2) / se), lower.tail = FALSE)
}

# Screen one outcome Y2 against all pre X1
screen_outcome <- function(y2,
                           data,
                           method_cor,
                           baseline = FALSE,
                           pre_vars = pre_vars,
                           pre_vars_reg = pre_vars_reg,
                           must_have = character()) {
  if (baseline == TRUE) {
    controls <- intersect(colnames(data), unique(c(must_have)))
    y2_trunk <- stringr::str_remove(y2, '_change$')
    Xcands <- stringr::str_subset(pre_vars, paste0('^', y2_trunk, pre_vars_reg))
  } else {
    y1 <- baseline_of(y2)
    controls <- intersect(colnames(data), unique(c(y1, must_have)))
    Xcands <- setdiff(pre_vars, controls)
  }
  
  if (length(Xcands) == 0) return(tibble())
  
  rows <- purrr::map(Xcands, function(x) {
    overall <- pcorr_one(data, y2, x, controls, 
                         method_cor = method_cor)
    kd <- pcorr_one(dplyr::filter(data, group == 'KD'), y2, x, controls, 
                    method_cor = method_cor)
    cd <- pcorr_one(dplyr::filter(data, group == 'CD'), y2, x, controls, 
                    method_cor = method_cor)
    
    p_het <- fisher_diff_p(kd$cor, kd$n, cd$cor, cd$n, length(controls))
    
    tibble(
      outcome = y2,
      X = x,
      n_controls = length(controls),
      cor_ttl = overall$cor,
      p_ttl = overall$p,
      cor_KD = kd$cor,
      p_KD = kd$p,
      cor_CD = cd$cor,
      p_CD = cd$p,
      cor_diff = abs(cd$cor - kd$cor),
      p_het = p_het
    )
  })
  
  purrr::list_rbind(rows) %>%
    mutate(p_adj = p.adjust(p_ttl, method = 'fdr'),
           across(c('cor_ttl', 'cor_KD', 'cor_CD', 'cor_diff'), ~ round(.x, 2)),
           across(starts_with('p_'), ~ round(.x, 3))) %>%
    arrange(desc(abs(cor_ttl)))
}

# Run screen
run_screen <- function(data,
                       response_vars,
                       must_have,
                       method_cor,
                       method_p = 'fdr',
                       baseline = FALSE,
                       pre_vars_reg = '_(1)$') {
  
  all_names <- colnames(data)
  pre_vars <- all_names %>% stringr::str_subset(pre_vars_reg)
  
  # Keep only outcomes that exist also as a baseline
  response_vars <- response_vars[baseline_of(response_vars) %in% all_names]
  
  res <- purrr::map(response_vars, ~ screen_outcome(.x,
                                                    data,
                                                    method_cor = method_cor,
                                                    baseline = baseline,
                                                    pre_vars = pre_vars,
                                                    pre_vars_reg = pre_vars_reg,
                                                    must_have = must_have))
  if (length(res) == 0) return(tibble())
  purrr::list_rbind(res)
}

# - Blomqvist's coefficient correction ---------------------

# Calculate correction coefficients
# - data: data.frame / tibble
# - change_dv: character vector
# - change_reg: character (regex)
# - base_reg: character (regex)
# - ICCs: named list
# - covariates: character vector

blomqvist_k <- function(data, change_dv, 
                        change_reg, base_reg, ICCs, covariates = character()) {
  
  # Keep only DVs that exist
  change_dv <- intersect(colnames(data), change_dv)
  covariates <- intersect(covariates, colnames(data))
  
  map_dfr(change_dv, function(x) {
    change_trunk <- str_remove(x, change_reg)
    ICC <- pluck(ICCs, change_trunk, .default = NA_real_)
    
    # Find the matching baseline (expect exactly one)
    baseline <- str_subset(colnames(data), paste0('^', change_trunk, base_reg))
    
    if (length(baseline) != 1L) {
      return(tibble(
        change_dv = x,
        trunk = change_trunk,
        baseline = if (length(baseline) == 0L) NA_character_ 
          else paste(baseline, collapse = ','),
        ICC = ICC,
        var_baseline = NA_real_,
        var_resid = NA_real_,
        k = NA_real_,
        note = 'baseline not found or not unique'
      ))
    }
    
    # Data just for this computation
    cols <- c(baseline, covariates)
    dd <- data %>%
      dplyr::select(all_of(cols)) %>%
      drop_na()
    
    if (nrow(dd) < 3L) {
      return(tibble(
        change_dv = x,
        trunk = change_trunk,
        baseline = baseline,
        ICC = ICC,
        var_baseline = NA_real_,
        var_resid = NA_real_,
        k = NA_real_,
        note = 'insufficient rows after drop_na'
      ))
    }
    
    var_baseline <- var(dd[[baseline]], na.rm = TRUE)
    
    if (length(covariates) == 0L) {
      var_resid <- var_baseline
    } else {
      fml <- as.formula(paste(baseline, '~', paste(covariates, collapse = ' + ')))
      fit <- lm(fml, data = dd)
      var_resid <- var(residuals(fit), na.rm = TRUE)
    }
    
    k <- if (is.na(ICC) || is.na(var_resid) || var_resid == 0) NA_real_ 
      else (1 - ICC) * var_baseline / var_resid
    
    tibble(
      change_dv = x,
      trunk = change_trunk,
      baseline = baseline,
      ICC = ICC,
      var_baseline = var_baseline,
      var_resid = var_resid,
      k = k,
      note = NA_character_
    )
  })
}

# WRS IMPROVEMENTS -----------------

# Correctly paralelised regtestMC
regtestMC_c <- function (x, y, regfun = tsreg, 
                         nboot = 600, alpha = 0.05, plotit = TRUE, 
                       grp = c(1:ncol(x)), 
                       nullvec = c(rep(0, length(grp))), xout = FALSE, 
                       outfun = outpro, SEED = TRUE, pr = TRUE, ...) 
{
  library(parallel)
  x <- as.matrix(x)
  p1 <- ncol(x) + 1
  p <- ncol(x)
  xy <- cbind(x, y)
  xy <- elimna(xy)
  x <- xy[, 1:p]
  y <- xy[, p1]
  if (xout) {
    if (pr) 
      print("Default for outfun is now outpro")
    m <- cbind(x, y)
    if (identical(outfun, outblp)) 
      flag = outblp(x, y, plotit = FALSE)$keep
    else flag <- outfun(x, plotit = FALSE, ...)$keep
    m <- m[flag, ]
    x <- m[, 1:p]
    y <- m[, p1]
  }
  x <- as.matrix(x)
  if (length(grp) != length(nullvec)) 
    stop("The arguments grp and nullvec must have the same length.")
  
  # Create cluster for parallel processing
  cl <- makeCluster(detectCores() - 1)
  on.exit(stopCluster(cl))
  
  # Export necessary variables and functions to the cluster
  clusterExport(cl, c("regbootMC", "x", "y", "regfun"), envir = environment())
  clusterEvalQ(cl, library(WRS))  # Adjust if needed for your specific package dependencies
  
  if (SEED) {
    clusterSetRNGStream(cl, 2)  # Set seed for reproducibility
  }
  
  data <- matrix(sample(length(y), size = length(y) * nboot, 
                        replace = TRUE), nrow = nboot)
  data_list <- split(data, row(data))  # Convert to list for parLapply
  
  # Use parLapply instead of mclapply
  bvec <- parLapply(cl, data_list, function(idx) regbootMC(idx, x, y, regfun))
  bvec <- do.call(cbind, bvec)  # Adjust based on actual output structure
  
  grp <- grp + 1
  est <- regfun(x, y)$coef
  estsub <- est[grp]
  bsub <- t(bvec[grp, ])
  
  if (length(grp) == 1) {
    m1 <- sum((bvec[grp, ] - est)^2)/(length(y) - 1)
    dis <- (bsub - estsub)^2/m1
  }
  if (length(grp) > 1) {
    mvec <- apply(bsub, 2, FUN = mean)
    m1 <- var(t(t(bsub) - mvec + estsub))
    dis <- mahalanobis(bsub, estsub, m1)
  }
  dis2 <- order(dis)
  dis <- sort(dis)
  critn <- floor((1 - alpha) * nboot)
  crit <- dis[critn]
  test <- mahalanobis(t(estsub), nullvec, m1)
  sig.level <- 1 - sum(test > dis)/nboot
  if (length(grp) == 2 && plotit) {
    plot(bsub, xlab = "Parameter 1", ylab = "Parameter 2")
    points(nullvec[1], nullvec[2], pch = 0)
    xx <- bsub[dis2[1:critn], ]
    xord <- order(xx[, 1])
    xx <- xx[xord, ]
    temp <- chull(xx)
    lines(xx[temp, ])
    lines(xx[c(temp[1], temp[length(temp)]), ])
  }
  list(test = test, crit = crit, p.value = sig.level, nullvec = nullvec, 
       est = estsub)
}

# Correctly paralelised regciMC for Windows
regciMC_c <- function(x, y, regfun = WRS::tsreg, nboot = 599, 
                      alpha = 0.05, plotit = FALSE,
                      pr = FALSE, null.val = NULL, 
                      method = 'hoch', xlab = 'Predictor 1',
                      ylab = 'Predictor 2', xout = FALSE, 
                      outfun = WRS::outpro, SEED = TRUE, ...) {
  
  # setup
  x <- as.matrix(x)
  p1 <- ncol(x) + 1
  p <- ncol(x)
  xy <- cbind(x, y)
  xy <- xy[stats::complete.cases(xy), , drop = FALSE]
  x <- xy[, 1:p, drop = FALSE]
  y <- xy[, p1]
  nrem <- length(y)
  
  # optional outlier filtering
  if (xout) {
    if (pr) print('Default for argument outfun is now outpro')
    m <- cbind(x, y)
    if (identical(outfun, WRS::outblp)) {
      flag <- WRS::outblp(x, y, plotit = FALSE)$keep
    } else {
      flag <- outfun(x, plotit = FALSE, ...)$keep
    }
    m <- m[flag, , drop = FALSE]
    x <- m[, 1:p, drop = FALSE]
    y <- m[, p1]
  }
  
  # point estimates
  estit <- regfun(x, y, ...)$coef
  if (is.null(null.val)) null.val <- rep(0, p1)
  
  if (identical(regfun, WRS::tsreg) && pr) {
    if (sum(duplicated(y)) > 0) {
      print('Duplicate values detected; tshdreg might have more power than tsreg')
    }
  }
  
  # bootstrap indices
  if (SEED) set.seed(2)
  idx_mat <- matrix(sample(length(y), size = length(y) * nboot, replace = TRUE),
                    nrow = length(y))
  data_list <- lapply(seq_len(nboot), function(i) idx_mat[, i])
  
  # cluster (Windows-friendly)
  n_workers <- max(1L, as.integer(floor(parallel::detectCores())))
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  # libs + RNG on workers
  parallel::clusterEvalQ(cl, { library(WRS) })
  if (SEED) parallel::clusterSetRNGStream(cl, 2)
  
  # IMPORTANT: avoid naming an argument 'x' in '...'
  bvec_list <- parallel::parLapply(
    cl, data_list,
    function(idx, Xmat, yvec, regfun, ...) {
      WRS::regbootMC(idx, Xmat, yvec, regfun, xout = FALSE, ...)
    },
    Xmat = x, yvec = y, regfun = regfun, ...
  )
  bvec <- do.call(cbind, bvec_list)  # p1 x nboot
  
  # assemble results
  p1 <- ncol(x) + 1
  regci <- matrix(0, p1, 6)
  vlabs <- c('Intercept', paste('Slope', seq_len(p1 - 1)))
  dimnames(regci) <- list(vlabs, c('ci.low', 'ci.up', 'Estimate', 'S.E.', 'p-value', 'p.adj'))
  
  ilow <- round((alpha / 2) * nboot); ihi <- nboot - ilow; ilow <- ilow + 1
  se <- rep(NA_real_, p1)
  sig.level <- rep(NA_real_, p1)
  
  for (i in seq_len(p1)) {
    temp <- (sum(bvec[i, ] < null.val[i]) + 0.5 * sum(bvec[i, ] == null.val[i])) / nboot
    sig.level[i] <- 2 * (min(temp, 1 - temp))
    bsort <- sort(bvec[i, ])
    regci[i, 1] <- bsort[ilow]
    regci[i, 2] <- bsort[ihi]
    se[i] <- sqrt(stats::var(bvec[i, ]))
  }
  
  if (p1 == 3 && plotit) {
    plot(bvec[2, ], bvec[3, ], xlab = xlab, ylab = ylab)
  }
  
  regci[, 3] <- estit
  regci[, 4] <- se
  regci[, 5] <- sig.level
  regci[, 6] <- NA_real_
  if (p1 >= 2) regci[2:p1, 6] <- p.adjust(sig.level[2:p1], method = method)
  
  list(regci = regci, n = nrem, n.keep = length(y))
}


# APA PLOTS ---------------------
# - APA theme --------

theme_apa <- function() {
theme_classic() +
  theme(
    axis.title.x = element_text(size = 12, family = 'sans'),
    axis.title.y = element_text(size = 12, family = 'sans'),
    axis.text = element_text(size = 11, family = 'sans', color = 'black'),
    axis.line = element_line(linewidth = 0.5, color = 'black'),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, family = 'sans', face = 'bold'),
    legend.position = 'right',
    legend.text = element_text(size = 11, family = 'sans'),
    legend.title = element_text(size = 12, family = 'sans')
  )
} 

# - Y-axis limits ---------------

limits <- function(x) {
  rng <- range(x, na.rm = TRUE)
  scales::expand_range(rng, mul = 0.15) 
}