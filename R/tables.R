#' ===================================================================
#' TABLE FORMATTING HELPER FUNCTIONS
#' ===================================================================
#' These functions apply the number formatting used in the output tables:
#' - Counts: 3 significant figures
#' - Rates/ratios/percentages: 1 decimal place
#' - Numbers ≥5 digits: use spaces (12 345)
#' - Numbers <5 digits: no spaces/commas (1234)
#' - Uncertainty intervals: (95% UI XXX–XXX) with en-dash
#' - First UI in paragraph includes "95% UI"
#' - Subsequent UIs just (XXX–XXX)

#' Format numbers for the output tables
#'
#' @param value Numeric value(s) to format - can be single value or vector
#' @param digits Number of decimal places (default: 0)
#' @param sig_figs Number of significant figures for counts (default: 3)
#' @param is_count Whether this is a count (uses sig figs) vs rate (uses decimal places)
#' @return Formatted string(s)
format_table_number <- function(value, digits = 0, sig_figs = 3, is_count = FALSE) {
  # Handle NULL or all NA values
  if (is.null(value) || all(is.na(value))) {
    return("NA")
  }

  # Vectorized function for formatting individual values
  format_single <- function(val) {
    if (is.na(val)) {
      return("NA")
    }

    if (is_count) {
      # For counts: use significant figures
      formatted <- signif(val, sig_figs)
      # Convert to string and format with spaces if ≥5 digits
      str_val <- formatC(formatted, format = "f", digits = 0, big.mark = "")
    } else {
      # For rates/ratios: use decimal places
      # For very small rates, use more decimal places to avoid rounding to zero
      adjusted_digits <- digits
      if (val > 0 && val < 0.1 && digits == 1) {
        adjusted_digits <- 3 # Use 3 decimal places for rates < 0.1
      }
      str_val <- formatC(val, format = "f", digits = adjusted_digits, big.mark = "")
    }

    # Spacing rule: spaces for numbers ≥5 digits
    if (nchar(gsub("[^0-9]", "", str_val)) >= 5) {
      # Add spaces for thousands separator
      final_digits <- ifelse(is_count, 0, adjusted_digits)
      str_val <- formatC(as.numeric(str_val), format = "f", digits = final_digits, big.mark = " ")
    }

    return(str_val)
  }

  # Apply formatting to each value in the vector
  sapply(value, format_single, USE.NAMES = FALSE)
}

#' Format uncertainty intervals for the output tables
#'
#' @param mean_val Mean/median value(s) - can be single value or vector
#' @param lower_val Lower bound(s) of uncertainty interval
#' @param upper_val Upper bound(s) of uncertainty interval
#' @param digits Number of decimal places (default: 0)
#' @param sig_figs Number of significant figures for counts (default: 3)
#' @param is_count Whether this is a count vs rate/ratio
#' @param include_ui_label Whether to include "95% UI" (TRUE for first in paragraph)
#' @param scale_suffix Suffix to add (e.g., "M", "B")
#' @return Formatted uncertainty interval string(s)
format_table_uncertainty <- function(mean_val, lower_val, upper_val,
                                      digits = 0, sig_figs = 3, is_count = FALSE,
                                      include_ui_label = TRUE, scale_suffix = "") {
  # Check for any NA values across all vectors
  if (is.null(mean_val) || is.null(lower_val) || is.null(upper_val)) {
    return("NA")
  }

  # Handle single NA check - use any() for vectors
  has_na <- any(is.na(mean_val)) || any(is.na(lower_val)) || any(is.na(upper_val))
  if (has_na && length(mean_val) == 1) {
    return("NA")
  }

  # Vectorized function for formatting individual uncertainty intervals
  format_single_ui <- function(m_val, l_val, u_val, include_label = include_ui_label) {
    if (is.na(m_val) || is.na(l_val) || is.na(u_val)) {
      return("NA")
    }

    # Format individual values
    mean_str <- format_table_number(m_val, digits, sig_figs, is_count)
    lower_str <- format_table_number(l_val, digits, sig_figs, is_count)
    upper_str <- format_table_number(u_val, digits, sig_figs, is_count)

    # Add suffix if provided
    if (scale_suffix != "") {
      mean_str <- paste0(mean_str, scale_suffix)
      lower_str <- paste0(lower_str, scale_suffix)
      upper_str <- paste0(upper_str, scale_suffix)
    }

    # Create uncertainty interval with en-dash
    if (include_label) {
      ui_text <- paste0("(95% UI ", lower_str, "–", upper_str, ")")
    } else {
      ui_text <- paste0("(", lower_str, "–", upper_str, ")")
    }

    return(paste0(mean_str, " ", ui_text))
  }

  # Handle vectorized input
  if (length(mean_val) > 1 || length(lower_val) > 1 || length(upper_val) > 1) {
    # For vectors, only include UI label on first element
    include_labels <- c(include_ui_label, rep(FALSE, max(length(mean_val), length(lower_val), length(upper_val)) - 1))
    return(mapply(format_single_ui, mean_val, lower_val, upper_val, include_labels, USE.NAMES = FALSE))
  } else {
    # For single values
    return(format_single_ui(mean_val, lower_val, upper_val, include_ui_label))
  }
}

#' Format numbers with scale suffixes (M, B)
#'
#' @param value Numeric value
#' @param lower_val Lower uncertainty bound
#' @param upper_val Upper uncertainty bound
#' @param is_count Whether this is a count vs rate
#' @param include_ui_label Whether to include "95% UI" label
#' @return Formatted string with appropriate scale
format_table_scaled <- function(value, lower_val, upper_val, is_count = TRUE, include_ui_label = TRUE) {
  if (any(is.na(c(value, lower_val, upper_val)))) {
    return("NA")
  }

  # Determine scale based on magnitude
  if (value >= 1e9) {
    scale_factor <- 1e9
    suffix <- "B"
    digits <- 2 # 2 decimals for billions
  } else if (value >= 1e6) {
    scale_factor <- 1e6
    suffix <- "M"
    digits <- 1 # 1 decimal for millions
  } else {
    scale_factor <- 1
    suffix <- ""
    digits <- ifelse(is_count, 0, 1)
  }

  # Scale values
  scaled_value <- value / scale_factor
  scaled_lower <- lower_val / scale_factor
  scaled_upper <- upper_val / scale_factor

  return(format_table_uncertainty(
    scaled_value, scaled_lower, scaled_upper,
    digits = digits, is_count = FALSE,
    include_ui_label = include_ui_label, scale_suffix = suffix
  ))
}

# Internal utility: choose a 1e6/1e9 scale factor and suffix; returns list(value, low, upp, suffix)
.format_scale_triplet <- function(value, low = NULL, upp = NULL) {
  if (is.null(value)) {
    return(list(value = NA, low = NA, upp = NA, suffix = ""))
  }
  if (is.list(value)) value <- value[[1]]
  if (!is.null(low) && is.list(low)) low <- low[[1]]
  if (!is.null(upp) && is.list(upp)) upp <- upp[[1]]
  if (length(value) > 1) value <- value[1]
  if (!is.null(low) && length(low) > 1) low <- low[1]
  if (!is.null(upp) && length(upp) > 1) upp <- upp[1]
  # Determine scale
  if (is.na(value)) {
    sf <- 1
    suf <- ""
  } else if (value >= 1e9) {
    sf <- 1e9
    suf <- "B"
  } else if (value >= 1e6) {
    sf <- 1e6
    suf <- "M"
  } else {
    sf <- 1
    suf <- ""
  }
  list(
    value = value / sf,
    low = if (!is.null(low) && !is.na(low)) low / sf else low,
    upp = if (!is.null(upp) && !is.na(upp)) upp / sf else upp,
    suffix = suf
  )
}

# Internal simplified wrapper for uncertainty formatting with automatic scaling
.format_scaled_uncertainty <- function(value, low = NULL, upp = NULL, digits_count = 0, include_ui_label = FALSE, is_count = TRUE) {
  trip <- .format_scale_triplet(value, low, upp)
  if (is.null(low) || is.null(upp)) {
    # point value only
    val_str <- format_table_number(trip$value, digits = ifelse(is_count && trip$suffix == "", 0, digits_count), is_count = is_count)
    return(paste0(val_str, trip$suffix))
  }
  # Digits: 2 dp for billions, 1 dp for millions; otherwise digits_count (0 when is_count)
  digits_use <- ifelse(trip$suffix == "B", 2, ifelse(trip$suffix == "M", 1, ifelse(is_count, 0, digits_count)))
  format_table_uncertainty(
    trip$value, trip$low, trip$upp,
    digits = digits_use, is_count = FALSE,
    include_ui_label = include_ui_label, scale_suffix = trip$suffix
  )
}

#' Generate a summary table of pandemic model projections
#'
#' Creates a markdown table summarizing key pandemic model projections including
#' expected pandemic size, annual deaths, and cumulative burden
#'
#' @param results Results object from run_PoT containing model outputs
#' @param period_years Number of years in the projection period (default: 75, for 2025-2100)
#' @param digits Number of decimal places for scaled values (default: 1)
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_summary_table <- function(results, period_years = NULL, digits = 1, conf = NULL, output_dir = NULL) {
  if (is.null(period_years)) {
    period_years <- if (!is.null(conf) && !is.null(conf$tables$projection_years)) conf$tables$projection_years else getOption("pandemic.tables.projection_years", 75)
  }
  # Validate input
  if (is.null(results$results) || is.null(results$burden)) {
    stop("Results object must contain model results and burden calculations")
  }

  # Extract required values
  single_event_mean <- results$results$boot
  single_event_low <- results$results$boot_low
  single_event_upp <- results$results$boot_upp

  yearly_deaths_mean <- results$burden$yearly_deaths
  yearly_deaths_low <- results$burden$yearly_deaths_low
  yearly_deaths_up <- results$burden$yearly_deaths_up

  cum_deaths_mean <- results$burden$cum_deaths
  cum_deaths_low <- results$burden$cum_deaths_low
  cum_deaths_up <- results$burden$cum_deaths_up

  # DALY calculations - first check if they exist
  if (!is.null(results$dalys)) {
    cum_dalys_mean <- results$dalys$cum_dalys
    cum_dalys_low <- results$dalys$cum_dalys_low
    cum_dalys_up <- results$dalys$cum_dalys_up
    has_dalys <- TRUE
  } else {
    has_dalys <- FALSE
  }

  # Format function with custom digits
  format_with_scale <- function(value, low, upp, include_ui_label = TRUE) {
    # Ensure we're working with single values, not vectors or lists
    if (is.list(value)) value <- value[[1]]
    if (is.list(low)) low <- low[[1]]
    if (is.list(upp)) upp <- upp[[1]]
    if (length(value) > 1) value <- value[1]
    if (length(low) > 1) low <- low[1]
    if (length(upp) > 1) upp <- upp[1]

    # Custom scaling with user-specified digits
    if (any(is.na(c(value, low, upp)))) {
      return("NA")
    }

    # Determine scale based on magnitude
    if (value >= 1e9) {
      scale_factor <- 1e9
      suffix <- "B"
    } else if (value >= 1e6) {
      scale_factor <- 1e6
      suffix <- "M"
    } else {
      scale_factor <- 1
      suffix <- ""
    }

    # Scale values
    scaled_value <- value / scale_factor
    scaled_lower <- low / scale_factor
    scaled_upper <- upp / scale_factor

    return(format_table_uncertainty(
      scaled_value, scaled_lower, scaled_upper,
      digits = digits, is_count = FALSE,
      include_ui_label = include_ui_label, scale_suffix = suffix
    ))
  }

  # Format for raw numbers with suffix and custom digits
  format_number_suffix <- function(value, low, upp, include_ui_label = FALSE) {
    .format_scaled_uncertainty(value, low, upp, digits_count = digits, include_ui_label = include_ui_label, is_count = TRUE)
  }

  # Prepare table rows
  rows <- c()

  # Single pandemic size (first UI includes label)
  rows <- c(rows, paste0(
    "|**Expected Size of a Single Pandemic**|",
    format_with_scale(single_event_mean, single_event_low, single_event_upp, include_ui_label = TRUE),
    "|"
  ))

  # Annual expected deaths (subsequent UI without label)
  rows <- c(rows, paste0(
    "|**Annual Expected Deaths**|",
    format_number_suffix(yearly_deaths_mean, yearly_deaths_low, yearly_deaths_up, include_ui_label = FALSE),
    "|"
  ))

  # Calculate average annual deaths over forecast period (2026-2100)
  avg_annual_deaths_mean <- cum_deaths_mean / period_years
  avg_annual_deaths_low <- cum_deaths_low / period_years
  avg_annual_deaths_up <- cum_deaths_up / period_years

  # Average annual deaths over forecast period (subsequent UI without label)
  rows <- c(rows, paste0(
    "|**Average Annual Deaths (2026-2100)**|",
    format_number_suffix(avg_annual_deaths_mean, avg_annual_deaths_low, avg_annual_deaths_up, include_ui_label = FALSE),
    "|"
  ))

  # Cumulative deaths (subsequent UI without label)
  rows <- c(rows, paste0(
    "|**Cumulative Deaths (2026-2100)**|",
    format_with_scale(cum_deaths_mean, cum_deaths_low, cum_deaths_up, include_ui_label = FALSE),
    "|"
  ))

  # Cumulative DALYs (if available) (subsequent UI without label)
  if (has_dalys) {
    rows <- c(rows, paste0(
      "|**Cumulative DALYs (2026-2100)**|",
      format_with_scale(cum_dalys_mean, cum_dalys_low, cum_dalys_up, include_ui_label = FALSE),
      "|"
    ))
  }

  # Create header with time period
  end_year <- 2100
  header <- c(
    paste0("# Summary of Pandemic Model Projections (2026-", end_year, ")"),
    "",
    "|Metric|Estimate (95% UI)|",
    "|---|---|"
  )

  # Add footnotes
  footnotes <- c(
    "",
    "_Annual Expected Deaths: rate × severity at reference population (8.2B). Average Annual Deaths: cumulative 2025-equivalent deaths divided by forecast period._"
  )

  # Combine header, rows, and footnotes into a single string with line breaks
  md_table <- paste(c(header, rows, footnotes), collapse = "\n")

  # Print the table to console for easy copying
  cat(md_table, "\n")

  save_table_if_enabled(md_table, "summary_table.md", output_dir, conf)
}

#' Generate a scenario comparison table for pandemic model results
#'
#' Creates a markdown table comparing baseline, best case, and worst case scenarios
#' for pandemic projections over the specified time period
#'
#' @param results Results object from run_PoT containing scenario outputs
#' @param period_years Number of years in the projection period (default: 75, for 2025-2100)
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_scenario_table <- function(results, period_years = NULL, conf = NULL, output_dir = NULL) {
  if (is.null(period_years)) {
    period_years <- if (!is.null(conf) && !is.null(conf$tables$projection_years)) conf$tables$projection_years else getOption("pandemic.tables.projection_years", 75)
  }
  # Validate input
  if (is.null(results$burden) || is.null(results$scenarios)) {
    stop("Results object must contain burden calculations and scenario results")
  }

  # Extract scenario data
  baseline <- results$burden
  best_case <- results$scenarios$best_case
  worst_case <- results$scenarios$worst_case

  # Format function for values with suffix
  format_number_suffix <- function(value, low = NULL, upp = NULL, include_ui_label = FALSE) {
    .format_scaled_uncertainty(value, low, upp, digits_count = 0, include_ui_label = include_ui_label, is_count = TRUE)
  }

  # Create the rows - starting with annual deaths
  rows <- c()

  # Annual exceedance rate (point values)
  baseline_rate_formatted <- format_table_number(baseline$yearly_rate, digits = 3, is_count = FALSE)
  best_case_rate_formatted <- format_table_number(best_case$yearly_rate, digits = 3, is_count = FALSE)
  worst_case_rate_formatted <- format_table_number(worst_case$yearly_rate, digits = 3, is_count = FALSE)

  rows <- c(rows, paste0(
    "|**Annual Exceedance Rate**|",
    baseline_rate_formatted, "|",
    best_case_rate_formatted, "|",
    worst_case_rate_formatted, "|"
  ))

  # Annual expected deaths (first UI includes label)
  rows <- c(rows, paste0(
    "|**Annual Deaths**|",
    format_number_suffix(baseline$yearly_deaths, baseline$yearly_deaths_low, baseline$yearly_deaths_up, include_ui_label = TRUE), "|",
    format_number_suffix(best_case$yearly_deaths, best_case$yearly_deaths_low, best_case$yearly_deaths_up, include_ui_label = FALSE), "|",
    format_number_suffix(worst_case$yearly_deaths, worst_case$yearly_deaths_low, worst_case$yearly_deaths_up, include_ui_label = FALSE), "|"
  ))

  # Cumulative deaths (subsequent UIs without label)
  rows <- c(rows, paste0(
    "|**Cumulative Deaths**|",
    format_number_suffix(baseline$cum_deaths, baseline$cum_deaths_low, baseline$cum_deaths_up, include_ui_label = FALSE), "|",
    format_number_suffix(best_case$cum_deaths, best_case$cum_deaths_low, best_case$cum_deaths_up, include_ui_label = FALSE), "|",
    format_number_suffix(worst_case$cum_deaths, worst_case$cum_deaths_low, worst_case$cum_deaths_up, include_ui_label = FALSE), "|"
  ))

  # Add DALY data if available
  if (!is.null(results$dalys)) {
    rows <- c(rows, paste0(
      "|**Cumulative DALYs**|",
      format_number_suffix(results$dalys$cum_dalys, results$dalys$cum_dalys_low, results$dalys$cum_dalys_up, include_ui_label = FALSE), "|",
      format_number_suffix(results$dalys$cum_dalys_best, results$dalys$cum_dalys_best_low, results$dalys$cum_dalys_best_up, include_ui_label = FALSE), "|",
      format_number_suffix(results$dalys$cum_dalys_worst, results$dalys$cum_dalys_worst_low, results$dalys$cum_dalys_worst_up, include_ui_label = FALSE), "|"
    ))
  }

  # Create header
  end_year <- 2100
  header <- c(
    paste0("# Pandemic Scenario Comparison (2026-", end_year, ")"),
    "",
    "|Metric|Baseline Projection (95% UI)|Best Case Scenario (95% UI)|Worst Case Scenario (95% UI)|",
    "|---|---|---|---|"
  )

  # Combine header and rows into a single string with line breaks
  md_table <- paste(c(header, rows), collapse = "\n")

  # Print the table to console for easy copying
  cat(md_table, "\n")

  save_table_if_enabled(md_table, "scenario_table.md", output_dir, conf)
}


#' Generate a validation summary table for pandemic model results
#'
#' Creates a markdown table summarizing model validation results including
#' year-specific observed vs predicted comparisons, model performance metrics,
#' and overall validation statistics
#'
#' @param validation_results Results object from validate_model containing validation outputs
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_validation_table <- function(validation_results, output_dir = NULL, conf = NULL) {
  # Validate input
  if (is.null(validation_results$model) || is.null(validation_results$observed) || is.null(validation_results$predicted)) {
    stop("Validation results object must contain model, observed, and predicted components")
  }

  # Extract components
  model_results <- validation_results$model
  observed_data <- validation_results$observed
  predicted_data <- validation_results$predicted

  # Format function for values with suffix
  format_number_suffix <- function(value, low = NULL, upp = NULL, include_ui_label = FALSE) {
    .format_scaled_uncertainty(value, low, upp, digits_count = 0, include_ui_label = include_ui_label, is_count = TRUE)
  }

  # Create model information section
  model_info <- c(
    "## Model Training Information",
    "",
    "|Parameter|Value|",
    "|---|---|",
    "|Training Period|1600-1950|",
    "|Validation Period|1950-2025|",
    # Fitted sample size: all training-window events entering the fit. The count
    # of events above the 100 000 2025-equivalent floor is exported separately
    # (see R/yaml_export.R).
    sprintf("|Fitted Sample Size (1600-1950 universe)|%s events|", format_table_number(model_results$specs$n, digits = 0, is_count = TRUE)),
    sprintf("|2025-equivalent Threshold|%s deaths|", format_number_suffix(model_results$specs$Threshold)),
    sprintf("|Lower Cutoff|%s deaths|", format_number_suffix(model_results$specs$Cutoff)),
    sprintf(
      "|Shape Parameter (ξ)|%s|",
      if (!is.null(model_results$params$xi_ci_lower) && !is.null(model_results$params$xi_ci_upper)) {
        # Use pre-calculated confidence intervals from analysis
        format_table_uncertainty(
          model_results$params$xi,
          model_results$params$xi_ci_lower,
          model_results$params$xi_ci_upper,
          digits = 3, is_count = FALSE, include_ui_label = TRUE
        )
      } else {
        # Fallback to point estimate with SE
        paste0(
          format_table_number(model_results$params$xi, digits = 3, is_count = FALSE),
          " (SE: ", format_table_number(model_results$params$xi_se, digits = 3, is_count = FALSE), ")"
        )
      }
    ),
    sprintf(
      "|Scale Parameter (β)|%s|",
      if (!is.null(model_results$params$beta_ci_lower) && !is.null(model_results$params$beta_ci_upper)) {
        # Use pre-calculated confidence intervals from analysis
        format_table_uncertainty(
          model_results$params$beta,
          model_results$params$beta_ci_lower,
          model_results$params$beta_ci_upper,
          digits = 0, is_count = TRUE, include_ui_label = FALSE
        )
      } else {
        # Fallback to point estimate only
        format_number_suffix(model_results$params$beta)
      }
    ),
    ""
  )

  # Get predicted values for observation years
  # For validation, we need to interpolate or calculate predicted cumulative values
  obs_years <- observed_data$year
  predicted_cum <- predicted_data$result$cum_deaths
  predicted_cum_low <- predicted_data$result$cum_deaths_low
  predicted_cum_up <- predicted_data$result$cum_deaths_up

  # Create year-specific comparison table
  comparison_header <- c(
    "## Year-Specific Observed vs Predicted Comparison",
    "",
    "|Year|Observed Cumulative Deaths|Predicted Cumulative Deaths (95% UI)|",
    "|---|---|---|"
  )

  # Create comparison rows
  comparison_rows <- c()
  for (i in seq_len(nrow(observed_data))) {
    year <- observed_data$year[i]
    obs_val <- observed_data$obs[i]

    # Calculate predicted cumulative up to this year (simple approximation)
    # This assumes linear accumulation over the validation period
    years_elapsed <- year - 1950
    total_years <- 2025 - 1950
    pred_fraction <- years_elapsed / total_years

    pred_val <- predicted_cum * pred_fraction
    pred_low <- predicted_cum_low * pred_fraction
    pred_up <- predicted_cum_up * pred_fraction

    obs_formatted <- format_number_suffix(obs_val)
    pred_formatted <- format_number_suffix(pred_val, pred_low, pred_up, include_ui_label = (i == 1))

    comparison_rows <- c(comparison_rows, paste0("|", year, "|", obs_formatted, "|", pred_formatted, "|"))
  }

  # Create overall performance summary
  total_observed <- max(observed_data$obs)
  performance_summary <- c(
    "",
    "## Overall Validation Summary",
    "",
    "|Metric|Value|",
    "|---|---|",
    sprintf("|Total Observed Deaths (1950-2025)|%s|", format_number_suffix(total_observed)),
    sprintf("|Total Predicted Deaths (95%% UI)|%s|", format_number_suffix(predicted_cum, predicted_cum_low, predicted_cum_up, include_ui_label = FALSE)),
    sprintf(
      "|Annual Predicted Rate|%s events/year|",
      if (is.list(predicted_data$result$yearly_rate)) {
        # Extract uncertainty intervals if available
        rate_mean <- predicted_data$result$yearly_rate$mean
        rate_low <- if (!is.null(predicted_data$result$yearly_rate$lower_95)) {
          predicted_data$result$yearly_rate$lower_95
        } else {
          NULL
        }
        rate_upp <- if (!is.null(predicted_data$result$yearly_rate$upper_95)) {
          predicted_data$result$yearly_rate$upper_95
        } else {
          NULL
        }

        # Format with uncertainty if available
        if (!is.null(rate_low) && !is.null(rate_upp)) {
          format_table_uncertainty(
            rate_mean, rate_low, rate_upp,
            digits = 3, is_count = FALSE, include_ui_label = TRUE
          )
        } else {
          format_table_number(rate_mean, digits = 3, is_count = FALSE)
        }
      } else {
        format_table_number(predicted_data$result$yearly_rate, digits = 3, is_count = FALSE)
      }
    ),
    sprintf("|Expected Single Event Size|%s deaths|", format_number_suffix(model_results$results$boot, model_results$results$boot_low, model_results$results$boot_upp, include_ui_label = FALSE)),
    ""
  )

  header <- c("# Pandemic Model Validation Results", "")

  # Combine all sections
  full_table <- c(
    header,
    model_info,
    comparison_header,
    comparison_rows,
    performance_summary
  )

  # Combine into a single string with line breaks
  md_table <- paste(full_table, collapse = "\n")

  # Print the table to console for easy copying
  cat(md_table, "\n")

  save_table_if_enabled(md_table, "validation_table.md", output_dir, conf)

  return(md_table)
}

#' Generate a health burden gap table between scenarios
#'
#' Creates a markdown table summarizing the health burden gap between worst and best
#' plausible scenarios, including both deaths and DALYs with uncertainty intervals
#'
#' @param results Results object from run_PoT containing scenario gap calculations
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_gap_table <- function(results, output_dir = NULL, conf = NULL) {
  # Validate input
  if (is.null(results$scenarios)) {
    stop("Results object must contain scenario calculations")
  }

  # Check if gap calculations exist
  if (is.null(results$scenarios$deaths_gap)) {
    stop("Gap calculations not found in the scenario results")
  }

  # Format function for values with suffix and uncertainty
  format_gap_value <- function(median_val, low_val, upp_val, is_percentage = FALSE, include_ui_label = FALSE) {
    if (is_percentage) {
      # Format as percentage with 1 decimal place
      return(format_table_uncertainty(
        median_val, low_val, upp_val,
        digits = 1, is_count = FALSE,
        include_ui_label = include_ui_label, scale_suffix = "%"
      ))
    } else {
      # Scaled formatting for absolute values
      return(format_table_scaled(median_val, low_val, upp_val, is_count = TRUE, include_ui_label = include_ui_label))
    }
  }

  # Extract deaths gap data
  deaths_gap <- results$scenarios$deaths_gap
  deaths_absolute <- format_gap_value(
    deaths_gap$absolute$median,
    deaths_gap$absolute$lower_95,
    deaths_gap$absolute$upper_95,
    is_percentage = FALSE,
    include_ui_label = TRUE
  )
  deaths_relative <- format_gap_value(
    deaths_gap$relative$median,
    deaths_gap$relative$lower_95,
    deaths_gap$relative$upper_95,
    is_percentage = TRUE,
    include_ui_label = FALSE
  )

  # Start building rows
  rows <- c()

  # Deaths gap row
  rows <- c(rows, paste0(
    "|**Deaths**|",
    deaths_absolute, "|",
    deaths_relative, "|"
  ))

  # Check if DALY gap exists and add if available
  if (!is.null(results$scenarios$dalys_gap)) {
    dalys_gap <- results$scenarios$dalys_gap
    dalys_absolute <- format_gap_value(
      dalys_gap$absolute$median,
      dalys_gap$absolute$lower_95,
      dalys_gap$absolute$upper_95,
      is_percentage = FALSE,
      include_ui_label = FALSE
    )
    dalys_relative <- format_gap_value(
      dalys_gap$relative$median,
      dalys_gap$relative$lower_95,
      dalys_gap$relative$upper_95,
      is_percentage = TRUE,
      include_ui_label = FALSE
    )

    rows <- c(rows, paste0(
      "|**DALYs**|",
      dalys_absolute, "|",
      dalys_relative, "|"
    ))
  }

  # Create header
  header <- c(
    "# Health Burden Gap: Worst vs Best Plausible Scenarios (2026-2100)",
    "",
    "|Health Metric|Absolute Gap (95% UI)|Relative Gap (95% UI)|",
    "|---|---|---|"
  )

  # Add note about interpretation
  note <- c(
    "",
    "_Note: Gap represents the additional burden in the worst-case scenario compared to the best-case scenario._"
  )

  # Combine all parts
  md_table <- paste(c(header, rows, note), collapse = "\n")

  # Print the table to console
  cat(md_table, "\n")

  save_table_if_enabled(md_table, "gap_table.md", output_dir, conf)
}

#' Generate a parameter summary table for pandemic model results
#'
#' Creates a markdown table summarizing key model parameters including shape and scale
#' parameters for the GPD tail distribution, DALY-to-death ratio, and linear trend
#' coefficients with their uncertainty intervals and statistical significance.
#'
#' @param results Results object from run_analysis containing model parameters and estimates
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_parameter_summary <- function(results, output_dir = NULL, conf = NULL) {
  # Validate input
  if (is.null(results$params)) {
    stop("Results object must contain model parameters")
  }

  # Initialize rows container
  rows <- c()

  #--------------------------------------------------
  # GPD TAIL PARAMETERS
  #--------------------------------------------------
  if (!is.null(results$params$xi) && !is.null(results$params$beta)) {
    # Shape parameter (ξ) - prefer bootstrap intervals if available
    xi_estimate <- results$params$xi
    xi_formatted <- NULL
    
    # Check for bootstrap-derived intervals first
    if (!is.null(results$params$xi_boot_median) && !is.null(results$params$xi_boot_lower) && !is.null(results$params$xi_boot_upper)) {
      xi_formatted <- format_table_uncertainty(
        results$params$xi_boot_median, results$params$xi_boot_lower, results$params$xi_boot_upper,
        digits = 3, is_count = FALSE, include_ui_label = TRUE
      )
    } else if (!is.null(results$params$xi_se)) {
      # Fallback to SE-based intervals
      xi_se <- results$params$xi_se
      xi_lower <- xi_estimate - 1.96 * xi_se
      xi_upper <- xi_estimate + 1.96 * xi_se
      xi_formatted <- format_table_uncertainty(
        xi_estimate, xi_lower, xi_upper,
        digits = 3, is_count = FALSE, include_ui_label = TRUE
      )
    } else {
      xi_formatted <- format_table_number(xi_estimate, digits = 3, is_count = FALSE)
    }

    rows <- c(rows, paste0(
      "|**Shape Parameter (ξ)**|",
      xi_formatted, "|",
      "GPD tail distribution|"
    ))

    # Scale parameter (β) - prefer bootstrap intervals if available
    beta_estimate <- results$params$beta
    beta_formatted <- NULL
    
    # Check for bootstrap-derived intervals first
    if (!is.null(results$params$beta_boot_median) && !is.null(results$params$beta_boot_lower) && !is.null(results$params$beta_boot_upper)) {
      beta_formatted <- format_table_uncertainty(
        results$params$beta_boot_median, results$params$beta_boot_lower, results$params$beta_boot_upper,
        digits = 0, is_count = TRUE, include_ui_label = FALSE
      )
    } else if (!is.null(results$params$beta_se)) {
      # Fallback to SE-based intervals
      beta_se <- results$params$beta_se
      beta_lower <- max(0, beta_estimate - 1.96 * beta_se)  # Ensure positive
      beta_upper <- beta_estimate + 1.96 * beta_se
      beta_formatted <- format_table_uncertainty(
        beta_estimate, beta_lower, beta_upper,
        digits = 0, is_count = TRUE, include_ui_label = FALSE
      )
    } else {
      beta_formatted <- format_table_number(beta_estimate, digits = 0, is_count = TRUE)
    }

    rows <- c(rows, paste0(
      "|**Scale Parameter (β)**|",
      beta_formatted, "|",
      "GPD tail distribution|"
    ))
  }

  #--------------------------------------------------
  # BULK DISTRIBUTION PARAMETERS
  #--------------------------------------------------
  # Add bulk distribution parameters if available from bootstrap
  if (!is.null(results$param_summary)) {
    bulk_params <- results$param_summary[results$param_summary$component == "bulk", ]
    
    if (nrow(bulk_params) > 0) {
      # Get the distribution type from the first bulk parameter
      dist_type <- if (!is.null(results$params$bulk_distribution)) {
        results$params$bulk_distribution
      } else {
        "Bulk distribution"
      }
      
      for (i in seq_len(nrow(bulk_params))) {
        param_name <- bulk_params$parameter[i]
        
        # Format parameter name for display
        display_name <- switch(param_name,
          "meanlog" = "Log Mean (μ)",
          "sdlog" = "Log SD (σ)", 
          "shape" = "Shape (α)",
          "scale" = "Scale (β)",
          "rate" = "Rate (β)",
          "bulk_prob" = "Bulk Probability",
          paste0(tools::toTitleCase(param_name))
        )
        
        # Format the parameter estimate with uncertainty
        param_formatted <- format_table_uncertainty(
          bulk_params$median[i], bulk_params$lower_95[i], bulk_params$upper_95[i],
          digits = if (param_name == "bulk_prob") 3 else 2,
          is_count = FALSE,
          include_ui_label = FALSE
        )
        
        rows <- c(rows, paste0(
          "|**", display_name, "**|",
          param_formatted, "|",
          paste0(tools::toTitleCase(dist_type), " bulk distribution|")
        ))
      }
    }
  }

  #--------------------------------------------------
  # DALY-TO-DEATH RATIO
  #--------------------------------------------------
  if (!is.null(results$dalys) && !is.null(results$dalys$daly_ratio)) {
    daly_ratio_formatted <- format_table_uncertainty(
      results$dalys$daly_ratio,
      results$dalys$daly_ratio_low,
      results$dalys$daly_ratio_up,
      digits = 1, is_count = FALSE, include_ui_label = FALSE
    )

    rows <- c(rows, paste0(
      "|**DALY-to-Death Ratio**|",
      daly_ratio_formatted, "|",
      "Burden conversion factor|"
    ))
  }

  #--------------------------------------------------
  # CREATE TABLE
  #--------------------------------------------------
  # Create header
  header <- c(
    "# Model Parameter Summary",
    "",
    "|Parameter|Estimate (95% UI)|Description|",
    "|---|---|---|"
  )

  # Add note about methodology
  note <- c(
    "",
    "_Note: Uncertainty intervals calculated using bootstrap quantiles (2.5%-97.5%) when available, otherwise normal approximation (parameter ± 1.96 × SE)._"
  )

  # Combine all parts
  md_table <- paste(c(header, rows, note), collapse = "\n")

  # Print the table to console
  cat(md_table, "\n")

  save_table_if_enabled(md_table, "parameter_summary.md", output_dir, conf)
}

#' Generate tier x pathogen x era composition table
#'
#' Crosstab showing event counts by tier (moderate/catastrophic),
#' pathogen group, and era (pre/post breakpoint).
#'
#' @param data Data frame with deaths_scaled, type, and start_year columns
#' @param variable Variable name for severity (e.g. "deaths_scaled")
#' @param tier_boundary Threshold separating moderate from catastrophic
#' @param lower_cutoff Lower inclusion threshold
#' @param breakpoint_year Year dividing eras (default 1950)
#' @param output_dir Output directory
#' @param conf Configuration object
#' @return Data frame with composition counts
generate_tier_composition_table <- function(data, variable, tier_boundary,
                                            lower_cutoff,
                                            breakpoint_year = 1950,
                                            output_dir = NULL, conf = NULL) {
  comp_data <- data %>%
    dplyr::filter(.data[[variable]] >= lower_cutoff) %>%
    dplyr::mutate(
      tier = dplyr::case_when(
        .data[[variable]] > tier_boundary ~ "Catastrophic",
        .data[[variable]] > lower_cutoff ~ "Moderate",
        TRUE ~ NA_character_
      ),
      era = ifelse(start_year < breakpoint_year, "Pre-1950", "Post-1950"),
      era = factor(era, levels = c("Pre-1950", "Post-1950")),
      pathogen_group = ifelse(is.na(type) | type == "", "Unknown", type)
    ) %>%
    dplyr::filter(!is.na(tier))

  composition <- comp_data %>%
    dplyr::count(tier, pathogen_group, era) %>%
    tidyr::pivot_wider(names_from = era, values_from = n, values_fill = 0) %>%
    dplyr::arrange(tier, dplyr::desc(`Pre-1950` + `Post-1950`))

  totals <- comp_data %>%
    dplyr::count(tier, era) %>%
    tidyr::pivot_wider(names_from = era, values_from = n, values_fill = 0) %>%
    dplyr::mutate(pathogen_group = "**Total**")

  composition <- dplyr::bind_rows(composition, totals) %>%
    dplyr::arrange(tier, pathogen_group == "**Total**", pathogen_group)

  md_lines <- c(
    "| Tier | Pathogen Group | Pre-1950 | Post-1950 |",
    "|------|---------------|----------|-----------|"
  )
  for (i in seq_len(nrow(composition))) {
    row <- composition[i, ]
    md_lines <- c(md_lines,
      paste0("| ", row$tier, " | ", row$pathogen_group,
             " | ", row$`Pre-1950`, " | ", row$`Post-1950`, " |"))
  }
  md_table <- paste(md_lines, collapse = "\n")

  save_table_if_enabled(md_table, "tier_composition.md", output_dir, conf)
  cat("\n", md_table, "\n")

  invisible(composition)
}

