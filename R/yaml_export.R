####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Functions for extracting model results and exporting to YAML format for use in
##              the manuscript (Quarto variables). This module provides:
##              1. Formatting helpers with en-dash negatives, thousand separators, and % suffixes
##              2. Uncertainty interval extraction with uniform precision
##              3. Per-tier rate extraction for tiered burden models
##              4. Full variable extraction from results, validation, and config objects
##              5. Manual YAML writer producing a quoted, section-commented _variables.yml
##
## Functions:   - format_yaml_number():   Format numbers with commas for YAML output
##              - format_yaml_percent():  Format percentages with % sign
##              - format_yaml_coefficient(): Format coefficients with en-dash for negatives
##              - format_yaml_pvalue():   Format p-values with <0.001 threshold
##              - extract_ui():           Extract uncertainty interval (median, low, upp)
##              - extract_tier_rate_info(): Extract per-tier rate model info
##              - extract_variables_from_results(): Main extraction from all result objects
##              - yaml_line():            Helper to produce a single YAML key-value line
##              - yaml_ui():              Helper to produce a YAML key with median/low/upp children
##              - write_variables_yaml(): Write formatted _variables.yml file
##
## Requires:    - Model results from run_analysis() with uncertainty quantification
##              - Validation results from validate_model()
##              - Configuration objects (conf, conf_validation)
##
## Outputs:     - _variables.yml file with all model results formatted as quoted strings
##
####################################################################################################

# =============================================================================
# Formatting helpers
# =============================================================================

#' Format number with commas for thousands, matching _variables.yml format
#'
#' @param value Numeric value to format
#' @param digits Number of decimal places (default: 0)
#' @return Character string with formatted number
format_yaml_number <- function(value, digits = 0) {
  if (is.null(value)) return("NA")
  if (length(value) > 1) value <- value[1]
  if (is.na(value)) return("NA")

  # Thousands separator:
  # - Numbers >=5 digits use single SPACE as thousands separator (e.g., "12 345")
  # - Numbers <5 digits use NO separator (e.g., "1234")
  big_mark <- if (abs(value) >= 10000) " " else ""
  formatted <- formatC(round(value, digits), format = "f", digits = digits, big.mark = big_mark)
  # Use en-dash for negative numbers (consistent with coefficient formatting)
  formatted <- gsub("^-", "\u2013", formatted)
  return(formatted)
}

#' Format percentage as character with % sign
#'
#' @param value Numeric value to convert to percentage
#' @param digits Number of decimal places (default: 1)
#' @param already_percent Whether value is already in percentage form (default: FALSE)
#' @return Character string with formatted percentage
format_yaml_percent <- function(value, digits = 1, already_percent = FALSE) {
  if (is.null(value)) return("NA")
  if (length(value) > 1) value <- value[1]
  if (is.na(value)) return("NA")

  percent_val <- if (already_percent) value else value * 100
  formatted <- formatC(round(percent_val, digits), format = "f", digits = digits)
  return(paste0(formatted, "%"))
}

#' Format coefficient with en-dash for negative values
#'
#' @param value Numeric coefficient value
#' @param digits Number of decimal places (default: 3)
#' @return Character string with formatted coefficient (en-dash U+2013 for negatives)
format_yaml_coefficient <- function(value, digits = 3) {
  if (is.null(value)) return("NA")
  if (length(value) > 1) value <- value[1]
  if (is.na(value)) return("NA")

  formatted <- formatC(abs(value), format = "f", digits = digits)
  if (value < 0) {
    formatted <- paste0("\u2013", formatted)
  }
  return(formatted)
}

#' Format p-value with <0.001 threshold
#'
#' @param value Numeric p-value
#' @param digits Number of decimal places (default: 3)
#' @return Character string: "<0.001" when very small, otherwise rounded
format_yaml_pvalue <- function(value, digits = 3) {
  if (is.null(value)) return("NA")
  if (length(value) > 1) value <- value[1]
  if (is.na(value)) return("NA")

  threshold <- 10^(-digits)
  if (value < threshold) {
    return(paste0("<", formatC(threshold, format = "f", digits = digits)))
  }
  return(formatC(round(value, digits), format = "f", digits = digits))
}

#' Format a number to N significant figures (IHME count rule)
#'
#' Rounds via signif() then formats at the natural decimal precision of the
#' rounded value. Adaptive: 103.2 → "103", 21.41 → "21.4", 4.91 → "4.91",
#' 0.62 → "0.62". Numbers >=10000 receive space thousands separator
#' (delegated to format_yaml_number).
#'
#' @param value Numeric value
#' @param sig_figs Number of significant figures (default 3 per IHME §Counts)
#' @return Character string
format_yaml_signif <- function(value, sig_figs = 3, pad = FALSE) {
  if (is.null(value)) return("NA")
  if (length(value) > 1) value <- value[1]
  if (is.na(value)) return("NA")

  rounded <- signif(value, sig_figs)
  # Determine dp from the rounded value's natural string form (no scientific)
  s  <- format(abs(rounded), scientific = FALSE)
  dp <- if (grepl("\\.", s)) nchar(sub("^.*\\.", "", s)) else 0
  # `pad`: keep a significant trailing zero that the string form drops, so a
  # column shows uniform significant figures (0.0196 -> "0.020", not "0.02").
  # Opt-in, because padding is wrong for values already at fixed precision.
  if (pad && rounded != 0) {
    dp <- max(dp, sig_figs - 1 - floor(log10(abs(rounded))), 0)
  }
  format_yaml_number(rounded, digits = dp)
}

# Event rates (events/year) are reported to 2 significant figures, the
# resolution of a 10 000-draw bootstrap.
RATE_SIG_FIGS <- 2

#' Extract uncertainty interval values (median, low, upp)
#'
#' Pass `digits` for fixed decimal places (IHME ratios/percentages/p-values),
#' or `sig_figs` for significant-figure formatting (IHME counts). The two
#' modes are mutually exclusive; whichever is non-NULL wins. Defaults to
#' fixed 0 dp when neither is provided.
#'
#' @param median_val Median value
#' @param low_val Lower bound
#' @param upp_val Upper bound
#' @param digits Number of decimal places (default: 0); ignored if sig_figs set
#' @param sig_figs Number of significant figures (default: NULL → use digits)
#' @return List with formatted median, low, and upp values
extract_ui <- function(median_val, low_val, upp_val, digits = 0, sig_figs = NULL,
                       pad = FALSE) {
  fmt <- if (is.null(sig_figs)) {
    function(v) format_yaml_number(v, digits)
  } else {
    function(v) format_yaml_signif(v, sig_figs, pad = pad)
  }
  list(median = fmt(median_val), low = fmt(low_val), upp = fmt(upp_val))
}

# =============================================================================
# Tier rate helper
# =============================================================================

#' Extract rate information from a per-tier rate list
#'
#' @param tier_rate A tier rate info list (e.g. results$burden$yearly_rate[[1]]$tier1_rate)
#' @return Named list with model_type, mean, trend coefficients, breakpoint_year, or NULL
extract_tier_rate_info <- function(tier_rate) {
  if (is.null(tier_rate)) return(NULL)

  info <- list(
    model_type = tier_rate$model_type %||% "Unknown",
    # Tier rate (events/year) — rate class; see RATE_SIG_FIGS
    mean       = format_yaml_signif(tier_rate$mean, RATE_SIG_FIGS, pad = TRUE)
  )

  # Trend coefficient if present
  # Trend coefficients are annual rate ratios near unity (e.g., 1.002);
  # 3 dp preserves an effect size that 2 dp would round to 1.00.
  if (!is.null(tier_rate$trend_coefficient)) {
    info$trend_coefficient <- format_yaml_coefficient(tier_rate$trend_coefficient, 3)
  }
  if (!is.null(tier_rate$trend_coefficient_lower)) {
    info$trend_ci_low <- format_yaml_coefficient(tier_rate$trend_coefficient_lower, 3)
  }
  if (!is.null(tier_rate$trend_coefficient_upper)) {
    info$trend_ci_upp <- format_yaml_coefficient(tier_rate$trend_coefficient_upper, 3)
  }

  # Regression coefficients (β_1 slope + SE; p-value).
  # Ratios: 2 dp; p-values: 2 dp (with a <0.01 floor).
  if (!is.null(tier_rate$beta_1_coefficient)) {
    info$beta_1 <- format_yaml_coefficient(tier_rate$beta_1_coefficient, 2)
  }
  if (!is.null(tier_rate$beta_1_se)) {
    info$beta_1_se <- format_yaml_number(tier_rate$beta_1_se, 2)
  }
  if (!is.null(tier_rate$beta_1_pvalue)) {
    info$beta_1_pvalue <- format_yaml_pvalue(tier_rate$beta_1_pvalue, 2)
  }

  # Breakpoint year (piecewise models)
  if (!is.null(tier_rate$breakpoint_year)) {
    info$breakpoint_year <- as.character(tier_rate$breakpoint_year)
  }

  # Breakpoint decline — level-shift coefficient (β_2) and its log-to-percent
  # transformation. Only emitted for piecewise models with a valid β_2.
  if (!is.null(tier_rate$beta_2_coefficient) &&
      !is.na(tier_rate$beta_2_coefficient) &&
      grepl("Piecewise", tier_rate$model_type %||% "", ignore.case = TRUE)) {
    b2  <- tier_rate$beta_2_coefficient
    se2 <- tier_rate$beta_2_se %||% 0
    info$beta_2        <- format_yaml_coefficient(b2, 2)
    info$beta_2_se     <- format_yaml_number(se2, 2)
    info$beta_2_pvalue <- format_yaml_pvalue(tier_rate$beta_2_pvalue, 2)
    info$breakpoint_decline <- list(
      median = format_yaml_percent((1 - exp(b2)) * 100, 1, already_percent = TRUE),
      low    = format_yaml_percent((1 - exp(b2 + 1.96 * se2)) * 100, 1, already_percent = TRUE),
      upp    = format_yaml_percent((1 - exp(b2 - 1.96 * se2)) * 100, 1, already_percent = TRUE)
    )
  }

  # Specific year rates — rate class (matches `mean` above)
  if (!is.null(tier_rate$rate_2026)) {
    info$rate_2026 <- format_yaml_signif(tier_rate$rate_2026, RATE_SIG_FIGS, pad = TRUE)
  }
  if (!is.null(tier_rate$rate_2100)) {
    info$rate_2100 <- format_yaml_signif(tier_rate$rate_2100, RATE_SIG_FIGS, pad = TRUE)
  }

  return(info)
}

# =============================================================================
# Sensitivity summary extraction
# =============================================================================

#' Extract formatted sensitivity summary for YAML export
#'
#' @param sens_summary List with cutoff, threshold, year sublists
#'   each containing min, max, mean, min_at, max_at (values in millions)
#' @return Nested list formatted for YAML output
extract_sensitivity_summary <- function(sens_summary) {
  format_sens <- function(s, format_at_fn = function(x) format_yaml_number(x, 0)) {
    if (is.null(s)) return(NULL)
    list(
      min    = format_yaml_number(s$min, 1),
      max    = format_yaml_number(s$max, 1),
      mean   = format_yaml_number(s$mean, 1),
      min_at = format_at_fn(s$min_at),
      max_at = format_at_fn(s$max_at)
    )
  }
  # Years should not have thousand separators (1800, not 1,800)
  format_year <- function(x) as.character(as.integer(x))
  list(
    cutoff    = format_sens(sens_summary$cutoff),
    threshold = format_sens(sens_summary$threshold),
    year      = format_sens(sens_summary$year, format_at_fn = format_year)
  )
}

#' Extract formatted HIV-exclusion sensitivity summary for YAML export
#'
#' Consumes the headline RDS produced by `src/run_hiv_sensitivity.R` and
#' formats deltas in millions for emission under `sensitivity.hiv_exclusion`.
#'
#' @param hiv_summary List produced by run_hiv_sensitivity.R
#' @return Nested list ready for YAML emission
extract_hiv_exclusion_summary <- function(hiv_summary) {
  if (is.null(hiv_summary)) return(NULL)

  cd <- hiv_summary$cumulative_deaths
  tr <- hiv_summary$tier_rates
  vs <- hiv_summary$validation

  out <- list(
    n_events_dropped = format_yaml_number(hiv_summary$meta$n_events_dropped, 0),
    cumulative_deaths = list(
      reference    = extract_ui(cd$reference$median  / 1e6,
                                cd$reference$lower_95 / 1e6,
                                cd$reference$upper_95 / 1e6, digits = 1),
      hiv_excluded = extract_ui(cd$hiv_excluded$median  / 1e6,
                                cd$hiv_excluded$lower_95 / 1e6,
                                cd$hiv_excluded$upper_95 / 1e6, digits = 1),
      abs_delta    = extract_ui(cd$abs_delta$median  / 1e6,
                                cd$abs_delta$lower_95 / 1e6,
                                cd$abs_delta$upper_95 / 1e6, digits = 1),
      pct_change   = format_yaml_number(cd$pct_change, 1)
    ),
    # Tier rates — rate class (matches main pipeline tier rate exports)
    tier_rates = list(
      reference = list(
        tier1 = format_yaml_signif(tr$reference$tier1_rate_mean, RATE_SIG_FIGS, pad = TRUE),
        tier2 = format_yaml_signif(tr$reference$tier2_rate_mean, RATE_SIG_FIGS, pad = TRUE)
      ),
      hiv_excluded = list(
        tier1 = format_yaml_signif(tr$hiv_excluded$tier1_rate_mean, RATE_SIG_FIGS, pad = TRUE),
        tier2 = format_yaml_signif(tr$hiv_excluded$tier2_rate_mean, RATE_SIG_FIGS, pad = TRUE)
      ),
      tier1_pct_change = format_yaml_number(tr$tier1_rate_pct_change, 1),
      tier2_pct_change = format_yaml_number(tr$tier2_rate_pct_change, 1)
    )
  )

  if (!is.null(vs)) {
    out$validation <- list(
      end_year   = as.character(as.integer(vs$end_year)),
      predicted  = extract_ui(vs$predicted$median   / 1e6,
                              vs$predicted$lower_95 / 1e6,
                              vs$predicted$upper_95 / 1e6, digits = 1),
      observed_reference     = format_yaml_number(vs$observed$reference     / 1e6, 1),
      observed_hiv_excluded  = format_yaml_number(vs$observed$hiv_excluded  / 1e6, 1),
      observed_abs_delta     = format_yaml_number(vs$observed$abs_delta     / 1e6, 1),
      observed_pct_delta     = format_yaml_number(vs$observed$pct_delta, 1),
      obs_in_pred_ui_reference    = ifelse(isTRUE(vs$obs_in_pred_ui$reference),    "true", "false"),
      obs_in_pred_ui_hiv_excluded = ifelse(isTRUE(vs$obs_in_pred_ui$hiv_excluded), "true", "false")
    )
  }
  out
}

# =============================================================================
# Main extraction function
# =============================================================================

#' Extract all variables from model results and validation results
#'
#' @param results Main analysis results object from run_analysis()
#' @param validation_results Validation results from validate_model()
#' @param conf Main analysis configuration
#' @param conf_validation Validation analysis configuration
#' @return Nested list with all extracted and formatted variables
extract_variables_from_results <- function(results, validation_results, conf, conf_validation,
                                           output_dir = NULL) {

  # ---------------------------------------------------------------------------
  # historical
  # ---------------------------------------------------------------------------
  period_start  <- as.character(conf$base_year)
  period_end    <- as.character(conf$time$current_year)
  forecast_start <- as.character(conf$time$current_year + 1)
  forecast_end   <- as.character(conf$time$current_year + conf$time$future_years)

  oos_obs_end   <- as.character(conf_validation$train_year_end)
  oos_pred_start <- as.character(conf_validation$test_year_start + 1)
  oos_pred_end   <- as.character(conf_validation$test_year_end)

  # Total database events and the PLE-qualifying count are set in config.yml
  # (they are not derivable from the event file).
  total_database_events <- conf$historical$total_database_events %||% NA
  events_ple_qualified  <- conf$historical$events_ple_qualified  %||% NA

  pre_1950_period_years  <- 1950 - conf$base_year
  post_1950_period_years <- conf$time$current_year - 1950

  # Crude historical event rates derived from era-tier counts. Distinct from the
  # MODELED forecast-horizon rates exported under probability.tier{1,2}_rate.*.
  # Per-century: 1 dp. Per-year: 2 dp.
  tier1_pre_1950_per_century  <- 100 * results$specs$moderate_pre_1950      / pre_1950_period_years
  tier1_post_1950_per_century <- 100 * results$specs$moderate_post_1950     / post_1950_period_years
  tier2_pre_1950_per_century  <- 100 * results$specs$catastrophic_pre_1950  / pre_1950_period_years
  tier2_post_1950_per_century <- 100 * results$specs$catastrophic_post_1950 / post_1950_period_years
  tier2_overall_per_year      <- results$specs$PoT /
    (pre_1950_period_years + post_1950_period_years)

  historical <- list(
    total_database_events       = as.character(total_database_events),
    events_ple_qualified        = as.character(events_ple_qualified),
    events_post_base_year       = as.character(results$specs$n),
    events_above_threshold      = as.character(results$specs$PoB),
    moderate_tier_events        = as.character(results$specs$PoB - results$specs$PoT),
    moderate_pre_1950           = as.character(results$specs$moderate_pre_1950),
    moderate_post_1950          = as.character(results$specs$moderate_post_1950),
    catastrophic_tier_events    = as.character(results$specs$PoT),
    catastrophic_pre_1950       = as.character(results$specs$catastrophic_pre_1950),
    catastrophic_post_1950      = as.character(results$specs$catastrophic_post_1950),
    pre_1950_period_years       = as.character(pre_1950_period_years),
    post_1950_period_years      = as.character(post_1950_period_years),
    tier1_pre_1950_per_century  = format_yaml_number(tier1_pre_1950_per_century,  1),
    tier1_post_1950_per_century = format_yaml_number(tier1_post_1950_per_century, 1),
    tier2_pre_1950_per_century  = format_yaml_number(tier2_pre_1950_per_century,  1),
    tier2_post_1950_per_century = format_yaml_number(tier2_post_1950_per_century, 1),
    tier2_overall_per_year      = format_yaml_number(tier2_overall_per_year,      2),
    data_min                    = format_yaml_number(conf$pop_min),
    period_start             = period_start,
    period_end               = period_end,
    forecast_start           = forecast_start,
    forecast_end             = forecast_end,
    oos_obs_end              = oos_obs_end,
    oos_pred_start           = oos_pred_start,
    oos_pred_end             = oos_pred_end
  )

  # ---------------------------------------------------------------------------
  # model
  # ---------------------------------------------------------------------------
  tail_threshold_millions <- format_yaml_number(results$specs$Threshold / 1e6, 0)

  median_pandemic_size <- extract_ui(
    results$results$boot / 1e6,
    results$results$boot_low / 1e6,
    results$results$boot_upp / 1e6,
    digits = 2
  )

  # Tier-specific conditional expected severity (in millions), 1 dp.
  bulk_conditional_mean <- extract_ui(
    results$results$bulk_mean_boot     / 1e6,
    results$results$bulk_mean_boot_low / 1e6,
    results$results$bulk_mean_boot_upp / 1e6,
    digits = 1
  )

  tail_conditional_mean <- extract_ui(
    results$results$tail_mean_boot     / 1e6,
    results$results$tail_mean_boot_low / 1e6,
    results$results$tail_mean_boot_upp / 1e6,
    digits = 1
  )

  model_section <- list(
    threshold              = format_yaml_number(conf$thresholds$lower_cutoff),
    tail_threshold         = tail_threshold_millions,
    forecast_years         = as.character(conf$time$future_years),
    median_pandemic_size   = median_pandemic_size,
    bulk_conditional_mean  = bulk_conditional_mean,
    tail_conditional_mean  = tail_conditional_mean,
    bulk_distribution      = as.character(results$params$bulk_distribution %||% "unknown"),
    constrained_survival   = as.character(isTRUE(results$params$constrained_survival))
  )

  # ---------------------------------------------------------------------------
  # probability
  # ---------------------------------------------------------------------------
  # Try yearly_rate first, then yearly_rate_detailed as fallback (tiered burden
  # stores the detailed rate info in yearly_rate_detailed while yearly_rate may be NA)
  annual_rate <- results$burden$yearly_rate
  rate_info <- NULL

  # Try extracting from yearly_rate
  if (is.list(annual_rate)) {
    if (length(annual_rate) == 1 && is.list(annual_rate[[1]])) {
      rate_info <- annual_rate[[1]]
    } else {
      rate_info <- annual_rate
    }
  }

  # Fallback: try yearly_rate_detailed (used by tiered burden)
  if (is.null(rate_info) || (is.list(rate_info) && isTRUE(rate_info$model_type == "failed"))) {
    yr_det <- results$burden$yearly_rate_detailed
    if (is.list(yr_det)) {
      if (length(yr_det) == 1 && is.list(yr_det[[1]])) {
        rate_info <- yr_det[[1]]
      } else {
        rate_info <- yr_det
      }
    }
  }

  # If rate_info has mean=NA and tiered rates both failed, treat as unavailable
  if (!is.null(rate_info) && is.list(rate_info) && is.na(rate_info$mean %||% NA) &&
      isTRUE(rate_info$tier1_rate$model_type == "failed")) {
    rate_info <- NULL
  }

  if (!is.null(rate_info)) {
    # List / tiered rate case
    rate_2025 <- rate_info$rate_2026 %||% rate_info$first_year_rate %||% rate_info$mean
    rate_2100 <- rate_info$rate_2100 %||% rate_info$last_year_rate %||% rate_info$mean

    mean_rate_ui <- extract_ui(
      rate_info$mean,
      rate_info$lower_95,
      rate_info$upper_95,
      digits = 3
    )

    beta_1_coeff   <- rate_info$beta_1_coefficient %||% 0
    beta_1_se      <- rate_info$beta_1_se %||% 0
    beta_1_pvalue  <- rate_info$beta_1_pvalue %||% 1
    beta_1_ci_low  <- beta_1_coeff - 1.96 * beta_1_se
    beta_1_ci_upp  <- beta_1_coeff + 1.96 * beta_1_se

    trend_coeff <- rate_info$trend_coefficient
    if (!is.null(trend_coeff)) {
      annual_decrease_pct <- (1 - trend_coeff) * 100
    } else {
      annual_decrease_pct <- 0
    }

    model_type <- rate_info$model_type %||% "Unknown"
  } else {
    # Constant scalar rate
    rate_2025 <- annual_rate
    rate_2100 <- annual_rate
    mean_rate_ui <- extract_ui(annual_rate, annual_rate, annual_rate, 3)
    beta_1_coeff   <- 0
    beta_1_ci_low  <- 0
    beta_1_ci_upp  <- 0
    beta_1_pvalue  <- 1
    annual_decrease_pct <- 0
    model_type <- "Constant Rate"
  }

  probability <- list(
    annual_percent = format_yaml_percent(rate_2025),
    at_2025        = format_yaml_number(rate_2025, 3),
    at_2100        = format_yaml_number(rate_2100, 3),
    mean_2026_2100 = mean_rate_ui,
    trend = list(
      model_type     = model_type,
      coefficient    = format_yaml_coefficient(beta_1_coeff, 3),
      ci_low         = format_yaml_coefficient(beta_1_ci_low, 3),
      ci_upp         = format_yaml_coefficient(beta_1_ci_upp, 3),
      p_value        = format_yaml_pvalue(beta_1_pvalue, 3),
      annual_decrease = format_yaml_percent(annual_decrease_pct, 2, already_percent = TRUE)
    )
  )

  # Per-tier rate blocks (only when tiered)
  if (!is.null(rate_info) && isTRUE(rate_info$tiered)) {
    probability$tier1_rate <- extract_tier_rate_info(rate_info$tier1_rate)
    probability$tier2_rate <- extract_tier_rate_info(rate_info$tier2_rate)

    # Per-draw scenario probability UIs derived from GLM rate draws +/- SFA multiplier
    tier_rate_scen <- results$scenarios$tier_rate_scenarios
    for (tname in c("tier1", "tier2")) {
      scen <- tier_rate_scen[[tname]]
      if (is.null(scen)) next
      tier_key <- paste0(tname, "_rate")
      for (scen_key in c("reference", "best_plausible", "worst_plausible")) {
        ui <- scen[[scen_key]]
        if (is.null(ui)) next
        # Tier rate scenarios — rate class (matches mean/rate_2026/rate_2100)
        probability[[tier_key]][[scen_key]] <- extract_ui(ui$median, ui$low, ui$upp, sig_figs = RATE_SIG_FIGS, pad = TRUE)
      }
    }
  }

  # Exceedance probabilities from saved draws. Exported as percentage strings
  # (e.g., "53.1%") at 1 dp.
  exceedance_thresholds <- c(25, 50, 100, 200, 500) * 1e6
  for (draws_file in c("cum_draws_reference.rds", "cum_draws_worst.rds", "cum_draws_best.rds")) {
    draws_path <- file.path(output_dir, draws_file)
    if (!is.null(output_dir) && file.exists(draws_path)) {
      draws <- readRDS(draws_path)
      scenario_key <- gsub("cum_draws_(.*)\\.rds", "\\1", draws_file)
      exceedance_probs <- list()
      for (t in exceedance_thresholds) {
        key <- paste0("exceed_", t / 1e6, "M")
        exceedance_probs[[key]] <- format_yaml_percent(mean(draws > t), digits = 1)
      }
      probability[[paste0("exceedance_", scenario_key)]] <- exceedance_probs
    }
  }

  # ---------------------------------------------------------------------------
  # reference
  # ---------------------------------------------------------------------------
  n_forecast <- conf$time$future_years

  # Counts in millions use 3 significant figures (103.2 -> "103", 21.4 -> "21.4").
  # Average annual deaths use a fixed 2 dp so that the three values of an
  # interval share a decimal place.
  ref_annual_deaths <- extract_ui(
    (results$burden$cum_deaths / n_forecast) / 1e6,
    (results$burden$cum_deaths_low / n_forecast) / 1e6,
    (results$burden$cum_deaths_up / n_forecast) / 1e6,
    digits = 2
  )

  ref_cumulative_deaths <- extract_ui(
    results$burden$cum_deaths / 1e6,
    results$burden$cum_deaths_low / 1e6,
    results$burden$cum_deaths_up / 1e6,
    sig_figs = 3
  )

  if (is.null(results$dalys)) {
    stop("DALY calculations are required but not found in results. Check calc_dalys setting.")
  }
  # Cumulative DALYs (billions): fixed 1 dp.
  total_dalys <- extract_ui(
    results$dalys$cum_dalys / 1e9,
    results$dalys$cum_dalys_low / 1e9,
    results$dalys$cum_dalys_up / 1e9,
    digits = 1
  )

  reference <- list(
    annual_deaths     = ref_annual_deaths,
    cumulative_deaths = ref_cumulative_deaths,
    total_dalys       = total_dalys
  )

  # ---------------------------------------------------------------------------
  # scenarios
  # ---------------------------------------------------------------------------
  if (is.null(results$scenarios)) {
    stop("Scenarios are required but not found. Check calc_scenarios and sfa_bounds in config.yml")
  }

  best_deaths <- extract_ui(
    results$scenarios$best_case$cum_deaths / 1e6,
    results$scenarios$best_case$cum_deaths_low / 1e6,
    results$scenarios$best_case$cum_deaths_up / 1e6,
    sig_figs = 3
  )
  worst_deaths <- extract_ui(
    results$scenarios$worst_case$cum_deaths / 1e6,
    results$scenarios$worst_case$cum_deaths_low / 1e6,
    results$scenarios$worst_case$cum_deaths_up / 1e6,
    sig_figs = 3
  )
  # Per-draw mean annual deaths: analysis.R:560-565 already summarizes
  # colMeans(yearly_deaths_adjusted_draws) into yearly_deaths{,_low,_up}.
  best_annual_deaths <- extract_ui(
    results$scenarios$best_case$yearly_deaths / 1e6,
    results$scenarios$best_case$yearly_deaths_low / 1e6,
    results$scenarios$best_case$yearly_deaths_up / 1e6,
    digits = 2
  )
  worst_annual_deaths <- extract_ui(
    results$scenarios$worst_case$yearly_deaths / 1e6,
    results$scenarios$worst_case$yearly_deaths_low / 1e6,
    results$scenarios$worst_case$yearly_deaths_up / 1e6,
    digits = 2
  )

  if (is.null(results$scenarios$deaths_gap)) {
    stop("Burden gap calculations required. Check calc_scenarios and sfa_bounds in config.yml")
  }

  burden_gap_deaths <- extract_ui(
    results$scenarios$deaths_gap$absolute$median / 1e6,
    results$scenarios$deaths_gap$absolute$lower_95 / 1e6,
    results$scenarios$deaths_gap$absolute$upper_95 / 1e6,
    sig_figs = 3
  )

  if (is.null(results$scenarios$dalys_gap)) {
    stop("DALY gap calculations required. Check calc_dalys and calc_scenarios in config.yml")
  }
  burden_gap_dalys <- extract_ui(
    results$scenarios$dalys_gap$absolute$median / 1e9,
    results$scenarios$dalys_gap$absolute$lower_95 / 1e9,
    results$scenarios$dalys_gap$absolute$upper_95 / 1e9,
    digits = 1
  )

  # Relative increase: use pre-computed if available, otherwise derive from absolute gap
  if (!is.null(results$scenarios$deaths_gap$relative)) {
    relative_increase <- extract_ui(
      results$scenarios$deaths_gap$relative$median,
      results$scenarios$deaths_gap$relative$lower_95,
      results$scenarios$deaths_gap$relative$upper_95,
      digits = 0
    )
  } else {
    # Compute from absolute gap / worst case * 100
    gap_abs <- results$scenarios$deaths_gap$absolute
    worst <- results$scenarios$worst_case
    rel_median <- gap_abs$median / worst$cum_deaths * 100
    rel_low <- gap_abs$lower_95 / worst$cum_deaths * 100
    rel_upp <- gap_abs$upper_95 / worst$cum_deaths * 100
    relative_increase <- extract_ui(rel_median, rel_low, rel_upp, digits = 0)
  }
  # Append % to relative increase values
  relative_increase$median <- paste0(relative_increase$median, "%")
  relative_increase$low    <- paste0(relative_increase$low, "%")
  relative_increase$upp    <- paste0(relative_increase$upp, "%")

  best_dalys <- extract_ui(
    results$dalys$cum_dalys_best / 1e9,
    results$dalys$cum_dalys_best_low / 1e9,
    results$dalys$cum_dalys_best_up / 1e9,
    digits = 1
  )
  worst_dalys <- extract_ui(
    results$dalys$cum_dalys_worst / 1e9,
    results$dalys$cum_dalys_worst_low / 1e9,
    results$dalys$cum_dalys_worst_up / 1e9,
    digits = 1
  )

  scenarios <- list(
    best_plausible = list(deaths = best_deaths, annual_deaths = best_annual_deaths, dalys = best_dalys),
    worst_plausible = list(deaths = worst_deaths, annual_deaths = worst_annual_deaths, dalys = worst_dalys),
    burden_gap = list(
      deaths            = burden_gap_deaths,
      dalys             = burden_gap_dalys,
      relative_increase = relative_increase
    )
  )

  # ---------------------------------------------------------------------------
  # daly_ratio
  # ---------------------------------------------------------------------------
  # DALY-to-death ratio: 1 dp, the precision its uncertainty interval supports.
  daly_ratio <- extract_ui(
    results$dalys$daly_ratio,
    results$dalys$daly_ratio_low,
    results$dalys$daly_ratio_up,
    digits = 1
  )

  # ---------------------------------------------------------------------------
  # validation
  # ---------------------------------------------------------------------------
  if (is.null(validation_results$predicted) || is.null(validation_results$predicted$result)) {
    stop("Validation predicted results are required but not found.")
  }
  pred_cum <- validation_results$predicted$result$cum_deaths
  pred_low <- validation_results$predicted$result$cum_deaths_low
  pred_upp <- validation_results$predicted$result$cum_deaths_up

  validation_predicted <- extract_ui(pred_cum / 1e6, pred_low / 1e6, pred_upp / 1e6, sig_figs = 3)

  if (is.null(validation_results$observed)) {
    stop("Validation observed results are required but not found.")
  }
  total_obs <- max(validation_results$observed$obs, na.rm = TRUE)
  validation_observed <- format_yaml_signif(total_obs / 1e6, sig_figs = 3)

  # Validation model parameters (shape, scale, training_n, threshold, event size)
  val_model <- validation_results$model
  # Shape parameter xi: 2 dp.
  val_shape <- extract_ui(
    val_model$params$xi,
    val_model$param_summary$lower_95[val_model$param_summary$parameter == "xi"],
    val_model$param_summary$upper_95[val_model$param_summary$parameter == "xi"],
    digits = 2
  )
  # GPD scale parameter beta: 3 significant figures.
  val_scale <- extract_ui(
    val_model$params$beta,
    val_model$param_summary$lower_95[val_model$param_summary$parameter == "beta"],
    val_model$param_summary$upper_95[val_model$param_summary$parameter == "beta"],
    sig_figs = 3
  )
  val_event_size <- extract_ui(
    val_model$results$boot / 1e6,
    val_model$results$boot_low / 1e6,
    val_model$results$boot_upp / 1e6,
    sig_figs = 3
  )

  # Per-decade forecast table
  val_table <- list()
  if (!is.null(validation_results$forecast)) {
    fc <- validation_results$forecast
    obs <- validation_results$observed
    for (i in seq_len(nrow(fc))) {
      yr <- paste0("y", fc$year[i])
      obs_val <- obs$obs[obs$year == fc$year[i]]
      val_table[[yr]] <- list(
        observed  = format_yaml_signif(obs_val / 1e6,            sig_figs = 3),
        predicted = format_yaml_signif(fc$cum_deaths[i] / 1e6,   sig_figs = 3),
        low       = format_yaml_signif(fc$cum_deaths_low[i] / 1e6, sig_figs = 3),
        upp       = format_yaml_signif(fc$cum_deaths_up[i] / 1e6,  sig_figs = 3)
      )
    }
  }

  validation <- list(
    period      = paste0(oos_pred_start, "\u2013", oos_pred_end),
    predicted   = validation_predicted,
    observed    = validation_observed,
    # Pre-1950 events above the 100 000 2025-equivalent floor, from the main
    # model's era counts (moderate + catastrophic). This is not
    # val_model$specs$n, the number of events in the 1600-1950 training window
    # that enters the GPD scale adjustment sigma = beta * (PoT/n)^xi, nor
    # val_model$specs$PoB, which is counted against the floor in the 1950
    # population frame of the validation profile.
    training_n  = as.character(results$specs$moderate_pre_1950 +
                                 results$specs$catastrophic_pre_1950),
    threshold   = format_yaml_number(val_model$specs$Threshold, 0),
    shape       = val_shape,
    scale       = val_scale,
    event_size  = val_event_size,
    table       = val_table
  )

  # ---------------------------------------------------------------------------
  # summary
  # ---------------------------------------------------------------------------
  summary_section <- list(
    main_deaths      = ref_cumulative_deaths,
    main_dalys       = total_dalys,
    scenario_range   = paste(best_deaths$median, "to", worst_deaths$median),
    # Already-formatted median (fixed 1 dp, in billions)
    burden_gap_dalys = burden_gap_dalys$median
  )

  # ---------------------------------------------------------------------------
  # periods
  # ---------------------------------------------------------------------------
  periods <- list(
    projection = paste(forecast_start, "and", forecast_end),
    historical = paste0(period_start, "\u2013", period_end)
  )

  # ---------------------------------------------------------------------------
  # population
  # ---------------------------------------------------------------------------
  population <- list(
    reference_year = period_end,
    reference_size = format_yaml_number(conf$pop_reference / 1e9, 1)
  )

  # ---------------------------------------------------------------------------
  # statistics
  # ---------------------------------------------------------------------------
  statistics <- list(
    tail_threshold      = tail_threshold_millions,
    poisson_window      = as.character(conf$window$default_size),
    bootstrap_draws     = format_yaml_number(conf$draws),
    gbd_posterior_draws = format_yaml_number(length(calc_daly_ratio_draws()))
  )

  # ---------------------------------------------------------------------------
  # parameters — from param_summary tibble
  # ---------------------------------------------------------------------------
  parameters <- NULL
  if (!is.null(results$param_summary) && nrow(results$param_summary) > 0) {
    ps <- results$param_summary

    make_param_ui <- function(param_name, component, digits = 0, sig_figs = NULL) {
      row <- ps[ps$parameter == param_name & ps$component == component, ]
      if (nrow(row) == 0) return(NULL)
      extract_ui(row$median[1], row$lower_95[1], row$upper_95[1],
                 digits = digits, sig_figs = sig_figs)
    }

    # Ratios: 2 dp; proportions: 3 dp on the proportion scale (1 dp on the
    # percent scale); large counts: 3 significant figures.
    # xi is a GPD shape parameter (ratio); meanlog is the lognormal
    # log-scale mean (ratio); tail/bulk_probability are mixture proportions;
    # tail_beta is the GPD scale parameter (count, in raw deaths).
    parameters <- list(
      tail_xi          = make_param_ui("xi",        "tail", digits = 2),
      tail_beta        = make_param_ui("beta",      "tail", sig_figs = 3),
      tail_probability = make_param_ui("tail_prob", "tail", digits = 3),
      bulk_mean        = make_param_ui("meanlog",   "bulk", digits = 2),
      bulk_sdlog       = make_param_ui("sdlog",     "bulk", digits = 2),
      bulk_probability = make_param_ui("bulk_prob", "bulk", digits = 3)
    )

    # Component conditional means (counts in raw deaths) — 3 sig figs.
    if (!is.null(results$results$bulk_mean_boot)) {
      parameters$bulk_conditional_mean <- extract_ui(
        results$results$bulk_mean_boot,
        results$results$bulk_mean_boot_low,
        results$results$bulk_mean_boot_upp,
        sig_figs = 3
      )
    }
    if (!is.null(results$results$tail_mean_boot)) {
      parameters$tail_conditional_mean <- extract_ui(
        results$results$tail_mean_boot,
        results$results$tail_mean_boot_low,
        results$results$tail_mean_boot_upp,
        sig_figs = 3
      )
    }
    # Overall conditional severity E[X|X>L] (count in raw deaths) — 3 sig figs.
    if (!is.null(results$results$bulk_mean_boot) && !is.null(results$results$tail_mean_boot)) {
      p_b <- ps$median[ps$parameter == "bulk_prob"]
      p_t <- ps$median[ps$parameter == "tail_prob"]
      overall_cond <- p_b * results$results$bulk_mean_boot + p_t * results$results$tail_mean_boot
      parameters$conditional_severity <- format_yaml_signif(overall_cond, 3)
    }
    # Remove any NULLs
    parameters <- parameters[!vapply(parameters, is.null, logical(1))]
    if (length(parameters) == 0) parameters <- NULL
  }

  # ---------------------------------------------------------------------------
  # software — static references used in prose
  # ---------------------------------------------------------------------------
  software <- list(
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )

  # ---------------------------------------------------------------------------
  # daly — static literature-reference constants
  # ---------------------------------------------------------------------------
  daly <- list(
    gbd_year = "2023"
  )

  # ---------------------------------------------------------------------------
  # sfa — from scenarios$sfa_bounds
  # ---------------------------------------------------------------------------
  sfa <- NULL
  if (!is.null(results$scenarios$sfa_bounds)) {
    extract_sfa_tier <- function(tier_key) {
      bounds <- results$scenarios$sfa_bounds[[tier_key]]$bounds
      if (is.null(bounds)) return(NULL)
      list(
        # Ratios: 2 dp.
        upper_multiplier    = format_yaml_number(bounds$upper_multiplier, 2),
        lower_multiplier    = format_yaml_number(bounds$lower_multiplier, 2),
        residual_skewness   = format_yaml_number(bounds$residual_skewness, 2)
      )
    }
    sfa <- list(
      tier1 = extract_sfa_tier("tier1"),
      tier2 = extract_sfa_tier("tier2")
    )
    # Remove NULLs
    sfa <- sfa[!vapply(sfa, is.null, logical(1))]
    if (length(sfa) == 0) sfa <- NULL
  }

  # ---------------------------------------------------------------------------
  # Assemble final list
  # ---------------------------------------------------------------------------
  list(
    historical  = historical,
    model       = model_section,
    probability = probability,
    reference   = reference,
    scenarios   = scenarios,
    daly_ratio  = daly_ratio,
    validation  = validation,
    summary     = summary_section,
    periods     = periods,
    population  = population,
    statistics  = statistics,
    parameters  = parameters,
    sfa         = sfa,
    software    = software,
    daly        = daly
  )
}

# =============================================================================
# YAML writer helpers
# =============================================================================

#' Produce a single YAML key: "value" line
#'
#' @param key YAML key name
#' @param value Character value (will be quoted)
#' @param indent Number of leading spaces (default: 2)
#' @return Character string with the formatted line
yaml_line <- function(key, value, indent = 2) {
  spaces <- paste0(rep(" ", indent), collapse = "")
  paste0(spaces, key, ': "', value, '"')
}

#' Produce YAML lines for a key with median/low/upp children
#'
#' @param key YAML key name
#' @param ui List with median, low, upp elements
#' @param indent Number of leading spaces for the parent key (default: 2)
#' @return Character vector with 4 lines (parent + 3 children)
yaml_ui <- function(key, ui, indent = 2) {
  spaces <- paste0(rep(" ", indent), collapse = "")
  child_spaces <- paste0(rep(" ", indent + 2), collapse = "")
  c(
    paste0(spaces, key, ":"),
    paste0(child_spaces, 'median: "', ui$median, '"'),
    paste0(child_spaces, 'low: "', ui$low, '"'),
    paste0(child_spaces, 'upp: "', ui$upp, '"')
  )
}

# =============================================================================
# Main YAML writer
# =============================================================================

#' Write formatted YAML file with proper structure, comments, and quoting
#'
#' @param yaml_data Nested list from extract_variables_from_results()
#' @param output_file Path to output YAML file (default: "output/_variables.yml")
write_variables_yaml <- function(yaml_data, output_file = "output/_variables.yml") {

  # Ensure output directory exists
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  lines <- character()

  # --- historical ---
  lines <- c(lines,
    "# Historical data",
    "historical:",
    yaml_line("total_database_events", yaml_data$historical$total_database_events),
    yaml_line("events_ple_qualified", yaml_data$historical$events_ple_qualified),
    yaml_line("events_post_base_year", yaml_data$historical$events_post_base_year),
    yaml_line("events_above_threshold", yaml_data$historical$events_above_threshold),
    yaml_line("moderate_tier_events", yaml_data$historical$moderate_tier_events),
    yaml_line("moderate_pre_1950", yaml_data$historical$moderate_pre_1950),
    yaml_line("moderate_post_1950", yaml_data$historical$moderate_post_1950),
    yaml_line("catastrophic_tier_events", yaml_data$historical$catastrophic_tier_events),
    yaml_line("catastrophic_pre_1950", yaml_data$historical$catastrophic_pre_1950),
    yaml_line("catastrophic_post_1950", yaml_data$historical$catastrophic_post_1950),
    yaml_line("pre_1950_period_years", yaml_data$historical$pre_1950_period_years),
    yaml_line("post_1950_period_years", yaml_data$historical$post_1950_period_years),
    yaml_line("tier1_pre_1950_per_century", yaml_data$historical$tier1_pre_1950_per_century),
    yaml_line("tier1_post_1950_per_century", yaml_data$historical$tier1_post_1950_per_century),
    yaml_line("tier2_pre_1950_per_century", yaml_data$historical$tier2_pre_1950_per_century),
    yaml_line("tier2_post_1950_per_century", yaml_data$historical$tier2_post_1950_per_century),
    yaml_line("tier2_overall_per_year", yaml_data$historical$tier2_overall_per_year),
    yaml_line("data_min", yaml_data$historical$data_min),
    yaml_line("period_start", yaml_data$historical$period_start),
    yaml_line("period_end", yaml_data$historical$period_end),
    yaml_line("forecast_start", yaml_data$historical$forecast_start),
    yaml_line("forecast_end", yaml_data$historical$forecast_end),
    yaml_line("oos_obs_end", yaml_data$historical$oos_obs_end),
    yaml_line("oos_pred_start", yaml_data$historical$oos_pred_start),
    yaml_line("oos_pred_end", yaml_data$historical$oos_pred_end)
  )

  # --- model ---
  lines <- c(lines, "",
    "# Model parameters",
    "model:",
    yaml_line("threshold", yaml_data$model$threshold),
    yaml_line("tail_threshold", yaml_data$model$tail_threshold),
    yaml_line("forecast_years", yaml_data$model$forecast_years),
    yaml_ui("median_pandemic_size", yaml_data$model$median_pandemic_size),
    yaml_ui("bulk_conditional_mean", yaml_data$model$bulk_conditional_mean),
    yaml_ui("tail_conditional_mean", yaml_data$model$tail_conditional_mean),
    yaml_line("bulk_distribution", yaml_data$model$bulk_distribution),
    yaml_line("constrained_survival", yaml_data$model$constrained_survival)
  )

  # --- probability ---
  lines <- c(lines, "",
    "# Probability estimates",
    "probability:",
    yaml_line("annual_percent", yaml_data$probability$annual_percent),
    yaml_line("at_2025", yaml_data$probability$at_2025),
    yaml_line("at_2100", yaml_data$probability$at_2100),
    yaml_ui("mean_2026_2100", yaml_data$probability$mean_2026_2100),
    "  trend:",
    yaml_line("model_type", yaml_data$probability$trend$model_type, indent = 4),
    yaml_line("coefficient", yaml_data$probability$trend$coefficient, indent = 4),
    yaml_line("ci_low", yaml_data$probability$trend$ci_low, indent = 4),
    yaml_line("ci_upp", yaml_data$probability$trend$ci_upp, indent = 4),
    yaml_line("p_value", yaml_data$probability$trend$p_value, indent = 4),
    yaml_line("annual_decrease", yaml_data$probability$trend$annual_decrease, indent = 4)
  )

  # Optional tier1_rate / tier2_rate blocks. Scalar fields print inline;
  # nested UI dicts (e.g. reference/best_plausible/worst_plausible) use yaml_ui.
  for (tier_key in c("tier1_rate", "tier2_rate")) {
    tier <- yaml_data$probability[[tier_key]]
    if (!is.null(tier)) {
      lines <- c(lines, paste0("  ", tier_key, ":"))
      for (field_name in names(tier)) {
        value <- tier[[field_name]]
        if (is.list(value) && all(c("median", "low", "upp") %in% names(value))) {
          lines <- c(lines, yaml_ui(field_name, value, indent = 4))
        } else {
          lines <- c(lines, yaml_line(field_name, value, indent = 4))
        }
      }
    }
  }

  # Exceedance probability blocks (cumulative 2026-2100)
  for (scenario_key in c("reference", "best", "worst")) {
    block_key <- paste0("exceedance_", scenario_key)
    block <- yaml_data$probability[[block_key]]
    if (!is.null(block)) {
      lines <- c(lines, paste0("  ", block_key, ":"))
      for (field_name in names(block)) {
        lines <- c(lines, yaml_line(field_name, block[[field_name]], indent = 4))
      }
    }
  }

  # --- reference ---
  lines <- c(lines, "",
    "# Reference scenario",
    "reference:",
    yaml_ui("annual_deaths", yaml_data$reference$annual_deaths),
    yaml_ui("cumulative_deaths", yaml_data$reference$cumulative_deaths),
    yaml_ui("total_dalys", yaml_data$reference$total_dalys)
  )

  # --- scenarios ---
  lines <- c(lines, "",
    "# Scenario analysis",
    "scenarios:",
    "  best_plausible:",
    "    deaths:",
    yaml_line("median", yaml_data$scenarios$best_plausible$deaths$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$best_plausible$deaths$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$best_plausible$deaths$upp, indent = 6),
    "    annual_deaths:",
    yaml_line("median", yaml_data$scenarios$best_plausible$annual_deaths$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$best_plausible$annual_deaths$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$best_plausible$annual_deaths$upp, indent = 6),
    "    dalys:",
    yaml_line("median", yaml_data$scenarios$best_plausible$dalys$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$best_plausible$dalys$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$best_plausible$dalys$upp, indent = 6),
    "  worst_plausible:",
    "    deaths:",
    yaml_line("median", yaml_data$scenarios$worst_plausible$deaths$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$worst_plausible$deaths$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$worst_plausible$deaths$upp, indent = 6),
    "    annual_deaths:",
    yaml_line("median", yaml_data$scenarios$worst_plausible$annual_deaths$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$worst_plausible$annual_deaths$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$worst_plausible$annual_deaths$upp, indent = 6),
    "    dalys:",
    yaml_line("median", yaml_data$scenarios$worst_plausible$dalys$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$worst_plausible$dalys$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$worst_plausible$dalys$upp, indent = 6),
    "  burden_gap:",
    "    deaths:",
    yaml_line("median", yaml_data$scenarios$burden_gap$deaths$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$burden_gap$deaths$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$burden_gap$deaths$upp, indent = 6),
    "    dalys:",
    yaml_line("median", yaml_data$scenarios$burden_gap$dalys$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$burden_gap$dalys$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$burden_gap$dalys$upp, indent = 6),
    "    relative_increase:",
    yaml_line("median", yaml_data$scenarios$burden_gap$relative_increase$median, indent = 6),
    yaml_line("low", yaml_data$scenarios$burden_gap$relative_increase$low, indent = 6),
    yaml_line("upp", yaml_data$scenarios$burden_gap$relative_increase$upp, indent = 6)
  )

  # --- daly_ratio ---
  lines <- c(lines, "",
    "# DALY conversion",
    yaml_ui("daly_ratio", yaml_data$daly_ratio, indent = 0)
  )

  # --- validation ---
  lines <- c(lines, "",
    "# Validation",
    "validation:",
    yaml_line("period", yaml_data$validation$period),
    yaml_ui("predicted", yaml_data$validation$predicted),
    yaml_line("observed", yaml_data$validation$observed),
    yaml_line("training_n", yaml_data$validation$training_n),
    yaml_line("threshold", yaml_data$validation$threshold),
    yaml_ui("shape", yaml_data$validation$shape),
    yaml_ui("scale", yaml_data$validation$scale),
    yaml_ui("event_size", yaml_data$validation$event_size)
  )

  # Per-decade validation table
  if (!is.null(yaml_data$validation$table) && length(yaml_data$validation$table) > 0) {
    lines <- c(lines, "  table:")
    for (yr in names(yaml_data$validation$table)) {
      row <- yaml_data$validation$table[[yr]]
      lines <- c(lines,
        paste0("    ", yr, ":"),
        yaml_line("observed", row$observed, indent = 6),
        yaml_line("predicted", row$predicted, indent = 6),
        yaml_line("low", row$low, indent = 6),
        yaml_line("upp", row$upp, indent = 6)
      )
    }
  }

  # --- summary ---
  lines <- c(lines, "",
    "# Summary",
    "summary:",
    yaml_ui("main_deaths", yaml_data$summary$main_deaths),
    yaml_ui("main_dalys", yaml_data$summary$main_dalys),
    yaml_line("scenario_range", yaml_data$summary$scenario_range),
    yaml_line("burden_gap_dalys", yaml_data$summary$burden_gap_dalys)
  )

  # --- periods ---
  lines <- c(lines, "",
    "# Time periods",
    "periods:",
    yaml_line("projection", yaml_data$periods$projection),
    yaml_line("historical", yaml_data$periods$historical)
  )

  # --- population ---
  lines <- c(lines, "",
    "# Population",
    "population:",
    yaml_line("reference_year", yaml_data$population$reference_year),
    yaml_line("reference_size", yaml_data$population$reference_size)
  )

  # --- statistics ---
  lines <- c(lines, "",
    "# Statistical methods",
    "statistics:",
    yaml_line("tail_threshold", yaml_data$statistics$tail_threshold),
    yaml_line("poisson_window", yaml_data$statistics$poisson_window),
    yaml_line("bootstrap_draws", yaml_data$statistics$bootstrap_draws),
    yaml_line("gbd_posterior_draws", yaml_data$statistics$gbd_posterior_draws)
  )

  # --- software ---
  lines <- c(lines, "",
    "# Software",
    "software:",
    yaml_line("r_version", yaml_data$software$r_version)
  )

  # --- daly (literature constants) ---
  lines <- c(lines, "",
    "# DALY / burden literature references",
    "daly:",
    yaml_line("gbd_year", yaml_data$daly$gbd_year)
  )

  # --- parameters (optional) ---
  if (!is.null(yaml_data$parameters)) {
    lines <- c(lines, "",
      "# Model parameter estimates",
      "parameters:"
    )
    for (param_name in names(yaml_data$parameters)) {
      val <- yaml_data$parameters[[param_name]]
      if (is.list(val) && !is.null(val$median)) {
        lines <- c(lines, yaml_ui(param_name, val))
      } else {
        lines <- c(lines, yaml_line(param_name, val))
      }
    }
  }

  # --- sfa (optional) ---
  if (!is.null(yaml_data$sfa)) {
    lines <- c(lines, "",
      "# SFA frontier bounds",
      "sfa:"
    )
    for (tier_name in names(yaml_data$sfa)) {
      tier <- yaml_data$sfa[[tier_name]]
      lines <- c(lines, paste0("  ", tier_name, ":"))
      for (field_name in names(tier)) {
        lines <- c(lines, yaml_line(field_name, tier[[field_name]], indent = 4))
      }
    }
  }

  # --- sensitivity (optional) ---
  if (!is.null(yaml_data$sensitivity)) {
    lines <- c(lines, "",
      "# Sensitivity analysis (cumulative deaths in millions, 2026-2100)",
      "sensitivity:"
    )
    for (dim_key in c("cutoff", "threshold", "year")) {
      dim <- yaml_data$sensitivity[[dim_key]]
      if (!is.null(dim)) {
        lines <- c(lines,
          paste0("  ", dim_key, ":"),
          yaml_line("min", dim$min, indent = 4),
          yaml_line("max", dim$max, indent = 4),
          yaml_line("mean", dim$mean, indent = 4),
          yaml_line("min_at", dim$min_at, indent = 4),
          yaml_line("max_at", dim$max_at, indent = 4)
        )
      }
    }

    # HIV-exclusion sensitivity
    hiv <- yaml_data$sensitivity$hiv_exclusion
    if (!is.null(hiv)) {
      lines <- c(lines,
        "  hiv_exclusion:",
        yaml_line("n_events_dropped", hiv$n_events_dropped, indent = 4),
        "    cumulative_deaths_M:",
        "      reference:",
        yaml_line("median", hiv$cumulative_deaths$reference$median, indent = 8),
        yaml_line("low",    hiv$cumulative_deaths$reference$low,    indent = 8),
        yaml_line("upp",    hiv$cumulative_deaths$reference$upp,    indent = 8),
        "      hiv_excluded:",
        yaml_line("median", hiv$cumulative_deaths$hiv_excluded$median, indent = 8),
        yaml_line("low",    hiv$cumulative_deaths$hiv_excluded$low,    indent = 8),
        yaml_line("upp",    hiv$cumulative_deaths$hiv_excluded$upp,    indent = 8),
        "      abs_delta:",
        yaml_line("median", hiv$cumulative_deaths$abs_delta$median, indent = 8),
        yaml_line("low",    hiv$cumulative_deaths$abs_delta$low,    indent = 8),
        yaml_line("upp",    hiv$cumulative_deaths$abs_delta$upp,    indent = 8),
        yaml_line("pct_change", hiv$cumulative_deaths$pct_change, indent = 6),
        "    tier_rates_per_year:",
        "      reference:",
        yaml_line("tier1", hiv$tier_rates$reference$tier1, indent = 8),
        yaml_line("tier2", hiv$tier_rates$reference$tier2, indent = 8),
        "      hiv_excluded:",
        yaml_line("tier1", hiv$tier_rates$hiv_excluded$tier1, indent = 8),
        yaml_line("tier2", hiv$tier_rates$hiv_excluded$tier2, indent = 8),
        yaml_line("tier1_pct_change", hiv$tier_rates$tier1_pct_change, indent = 6),
        yaml_line("tier2_pct_change", hiv$tier_rates$tier2_pct_change, indent = 6)
      )
      if (!is.null(hiv$validation)) {
        v <- hiv$validation
        lines <- c(lines,
          "    validation:",
          yaml_line("end_year", v$end_year, indent = 6),
          "      predicted_M:",
          yaml_line("median", v$predicted$median, indent = 8),
          yaml_line("low",    v$predicted$low,    indent = 8),
          yaml_line("upp",    v$predicted$upp,    indent = 8),
          yaml_line("observed_reference_M",     v$observed_reference,     indent = 6),
          yaml_line("observed_hiv_excluded_M",  v$observed_hiv_excluded,  indent = 6),
          yaml_line("observed_abs_delta_M",     v$observed_abs_delta,     indent = 6),
          yaml_line("observed_pct_delta",       v$observed_pct_delta,     indent = 6),
          yaml_line("obs_in_pred_ui_reference",    v$obs_in_pred_ui_reference,    indent = 6),
          yaml_line("obs_in_pred_ui_hiv_excluded", v$obs_in_pred_ui_hiv_excluded, indent = 6)
        )
      }
    }
  }

  # Final trailing newline
  lines <- c(lines, "")

  # Write
  writeLines(lines, output_file)

  cat("Successfully generated", output_file, "\n")
  cat("Variables extracted from model run and formatted for YAML output.\n")
}
