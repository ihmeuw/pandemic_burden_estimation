#' Get DALY to Death Ratio Draws for Uncertainty Propagation
#'
#' Retrieves pre-generated draws of the DALY-to-death ratio for uncertainty
#' quantification in burden estimation. These draws capture the uncertainty in
#' the relationship between deaths and disability-adjusted life years.
#'
#' @details
#' The function loads a set of pre-generated draws that represent the uncertainty
#' in the DALY-to-death ratio. Each draw is a possible value of the ratio based on
#' sampling from the joint distribution of DALYs and deaths. These draws are used
#' for Monte Carlo simulation in the burden calculation pipeline to properly
#' propagate uncertainty from multiple sources.
#'
#' @param path Path to the RDS file containing the DALY-to-death ratio draws
#' @return A vector of DALY-to-death ratio draws
#' @export
calc_daly_ratio_draws <- function(path = "data/daly_death_ratio_draws.rds") {
  ratio_draws <- readRDS(path)
  return(ratio_draws)
}

#' Summarize draw matrix rows with median and 95% uncertainty interval
#'
#' @param draw_matrix Numeric matrix with draws in columns and time in rows
#' @return Tibble containing per-row median, 2.5th percentile, and 97.5th percentile
summarize_draw_rows <- function(draw_matrix) {
  if (is.null(draw_matrix)) {
    return(tibble())
  }

  if (is.null(dim(draw_matrix))) {
    draw_matrix <- matrix(draw_matrix, nrow = length(draw_matrix), ncol = 1)
  }

  med <- apply(draw_matrix, 1, stats::median, na.rm = TRUE)
  low <- apply(draw_matrix, 1, function(x) stats::quantile(x, 0.025, na.rm = TRUE))
  upp <- apply(draw_matrix, 1, function(x) stats::quantile(x, 0.975, na.rm = TRUE))

  tibble(
    median = med,
    lower = low,
    upper = upp
  )
}

#' Calculate Tiered Burden with Independent Rate Models per Tier
#'
#' Decomposes burden into two tiers matching the mixture model components:
#' Tier 1 (bulk: lower_cutoff to tier_boundary) and Tier 2 (tail: above tier_boundary).
#' Each tier gets its own rate model, then draws are summed for total burden.
#'
#' @param data Data frame with pandemic data
#' @param variable Variable name (e.g., "deaths_scaled")
#' @param conf Configuration object
#' @param bulk_severity_draws Numeric vector of bootstrap draws for E[X | lower < X <= boundary]
#' @param tail_severity_draws Numeric vector of bootstrap draws for E[X | X > boundary]
#' @param shadow_mean Combined shadow mean point estimate (for reporting)
#' @param shadow_low Combined shadow mean lower bound (for reporting)
#' @param shadow_up Combined shadow mean upper bound (for reporting)
#' @param output_dir Output directory for plots
#' @param validation Logical, if TRUE uses validation-mode population data
#' @param window_sizes Window sizes for rate estimation
#' @param return_severity_draws Logical, return per-tier severity draw matrices
#' @return List with combined result, per-tier results, and draw matrices
calc_tiered_burden <- function(
    data, variable, conf,
    bulk_severity_draws,
    tail_severity_draws,
    shadow_mean, shadow_low, shadow_up,
    output_dir = NULL,
    validation = FALSE,
    window_sizes = NULL,
    return_severity_draws = FALSE) {

  tier_boundary <- conf$thresholds[[variable]]
  if (is.null(tier_boundary)) tier_boundary <- conf$thresholds$deaths_scaled
  lower_cutoff <- conf$thresholds$lower_cutoff
  base_year <- conf$base_year
  current_year <- conf$time$current_year
  future_yrs <- conf$time$future_years

  n_draws <- min(length(bulk_severity_draws), length(tail_severity_draws))
  bulk_severity_draws <- bulk_severity_draws[seq_len(n_draws)]
  tail_severity_draws <- tail_severity_draws[seq_len(n_draws)]

  cat("\n=== TIERED BURDEN CALCULATION ===\n")
  cat("Tier 1: events in (", lower_cutoff, ", ", tier_boundary, "]\n")
  cat("Tier 2: events > ", tier_boundary, "\n")
  cat("Bulk severity draws:", n_draws, "draws, median =",
      median(bulk_severity_draws), "\n")
  cat("Tail severity draws:", n_draws, "draws, median =",
      median(tail_severity_draws), "\n\n")

  # --- Tier 1: Override conf$time$model_type for this tier ---
  # The tier-specific model_type allows each tier to use a different temporal
  # regression (e.g., piecewise for Tier 1, nb_linear for Tier 2).
  # We mutate a local copy of conf so downstream functions (calc_burden ->
  # analyze_pandemic_trend -> compare_temporal_models) pick it up without
  # any signature changes.
  tier1_conf <- conf
  tier1_conf$time$model_type <- conf$time$tier1_model_type %||% conf$time$model_type

  cat("--- Tier 1: Bulk rate model (", tier1_conf$time$model_type, ") ---\n")
  tier1 <- tryCatch({
    calc_burden(
      data = data,
      var = variable,
      thresh = lower_cutoff,
      upper_thresh = tier_boundary,
      start = base_year,
      end = current_year,
      future_yrs = future_yrs,
      shadow_mean = median(bulk_severity_draws),
      shadow_low = quantile(bulk_severity_draws, 0.025),
      shadow_up = quantile(bulk_severity_draws, 0.975),
      use_time_trend = conf$time$use_time_trend,
      severity_draws = bulk_severity_draws,
      conf = tier1_conf,
      return_rate_draws = TRUE,
      return_yearly_deaths_draws = return_severity_draws,
      return_yearly_deaths_adjusted_draws = TRUE,
      return_forecast = TRUE,
      validation = validation,
      window_sizes = window_sizes
    )
  }, error = function(e) {
    warning("Tier 1 burden calculation failed: ", e$message)
    NULL
  })

  # --- Tier 2: Override conf$time$model_type for this tier ---
  # Same pattern as Tier 1 above: local conf copy with tier-specific model_type.
  tier2_conf <- conf
  tier2_conf$time$model_type <- conf$time$tier2_model_type %||% conf$time$model_type

  cat("--- Tier 2: Tail rate model (", tier2_conf$time$model_type, ") ---\n")
  tier2 <- tryCatch({
    calc_burden(
      data = data,
      var = variable,
      thresh = tier_boundary,
      upper_thresh = NULL,
      start = base_year,
      end = current_year,
      future_yrs = future_yrs,
      shadow_mean = median(tail_severity_draws),
      shadow_low = quantile(tail_severity_draws, 0.025),
      shadow_up = quantile(tail_severity_draws, 0.975),
      use_time_trend = conf$time$use_time_trend,
      severity_draws = tail_severity_draws,
      conf = tier2_conf,
      return_rate_draws = TRUE,
      return_yearly_deaths_draws = return_severity_draws,
      return_yearly_deaths_adjusted_draws = TRUE,
      return_forecast = TRUE,
      validation = validation,
      window_sizes = window_sizes
    )
  }, error = function(e) {
    warning("Tier 2 burden calculation failed: ", e$message)
    NULL
  })

  if (is.null(tier1) && is.null(tier2)) {
    stop("Both tier calculations failed")
  }

  cat("\n--- Combining tiers ---\n")

  tier1_draws <- if (!is.null(tier1)) tier1$yearly_deaths_adjusted_draws else
    matrix(0, nrow = future_yrs, ncol = n_draws)
  tier2_draws <- if (!is.null(tier2)) tier2$yearly_deaths_adjusted_draws else
    matrix(0, nrow = future_yrs, ncol = n_draws)

  n_combined <- min(ncol(tier1_draws), ncol(tier2_draws))
  tier1_draws <- tier1_draws[, seq_len(n_combined), drop = FALSE]
  tier2_draws <- tier2_draws[, seq_len(n_combined), drop = FALSE]

  total_adjusted_draws <- tier1_draws + tier2_draws

  cum_total_draws <- apply(total_adjusted_draws, 2, cumsum)
  cum_final <- cum_total_draws[nrow(cum_total_draws), ]

  cum_deaths_median <- stats::median(cum_final)
  cum_deaths_low <- as.numeric(stats::quantile(cum_final, 0.025))
  cum_deaths_up <- as.numeric(stats::quantile(cum_final, 0.975))

  # Save cumulative draw vector for exceedance probability analysis
  if (!is.null(output_dir)) {
    saveRDS(cum_final, file.path(output_dir, "cum_draws_reference.rds"))
    cat("Saved", length(cum_final), "cumulative burden draws to cum_draws_reference.rds\n")
    # Also save full year×draw cumulative trajectory for sensitivity comparisons
    forecast_years_traj <- seq_len(nrow(cum_total_draws)) + (current_year)
    rownames(cum_total_draws) <- as.character(forecast_years_traj)
    saveRDS(cum_total_draws, file.path(output_dir, "cum_traj_reference.rds"))
    cat("Saved cumulative trajectory matrix (",
        nrow(cum_total_draws), "yrs ×", ncol(cum_total_draws),
        "draws) to cum_traj_reference.rds\n")
  }

  yearly_summary <- summarize_draw_rows(total_adjusted_draws)
  cum_summary <- summarize_draw_rows(cum_total_draws)

  cat("Tier 1 cumulative (median):", stats::median(colSums(tier1_draws)) / 1e6, "M\n")
  cat("Tier 2 cumulative (median):", stats::median(colSums(tier2_draws)) / 1e6, "M\n")
  cat("Combined cumulative (median):", cum_deaths_median / 1e6, "M\n\n")

  # Extract per-tier rate info from calc_burden() results.
  # calc_burden() stores yearly_rate as a list-column in the result tibble,
  # so $yearly_rate returns the list directly (not wrapped in another list).
  extract_tier_rate <- function(tier_result) {
    if (is.null(tier_result)) return(list(mean = NA, model_type = "failed"))
    yr <- tier_result$result$yearly_rate
    if (is.list(yr)) return(yr)
    list(mean = NA, model_type = "failed")
  }
  tier1_rate_info <- extract_tier_rate(tier1)
  tier2_rate_info <- extract_tier_rate(tier2)

  yearly_rate_detailed <- list(
    mean = sum(tier1_rate_info$mean %||% 0, tier2_rate_info$mean %||% 0),
    tier1_rate = tier1_rate_info,
    tier2_rate = tier2_rate_info,
    tiered = TRUE
  )

  # Average across all forecast years
  avg_yearly_median <- stats::median(rowMeans(total_adjusted_draws))
  avg_yearly_low <- stats::quantile(colMeans(total_adjusted_draws), 0.025)
  avg_yearly_up <- stats::quantile(colMeans(total_adjusted_draws), 0.975)

  result <- tibble::tibble(
    yearly_rate = list(yearly_rate_detailed),
    yearly_deaths = avg_yearly_median,
    yearly_deaths_low = avg_yearly_low,
    yearly_deaths_up = avg_yearly_up,
    cum_deaths = cum_deaths_median,
    cum_deaths_low = cum_deaths_low,
    cum_deaths_up = cum_deaths_up,
    dynamic_rates = TRUE
  )

  # Build forecast data frame for combined cumulative plot
  forecast_years <- seq(current_year + 1, current_year + future_yrs)
  data_forecast <- tibble::tibble(
    year = forecast_years,
    cum_deaths = cum_summary$median,
    cum_deaths_low = cum_summary$lower,
    cum_deaths_up = cum_summary$upper,
    pred_deaths = yearly_summary$median,
    pred_deaths_low = yearly_summary$lower,
    pred_deaths_upp = yearly_summary$upper
  )

  # Per-tier cumulative summaries for stacked-lines plot
  cum_tier1 <- apply(tier1_draws, 2, cumsum)
  cum_tier2 <- apply(tier2_draws, 2, cumsum)
  tier1_cum_summary <- summarize_draw_rows(cum_tier1)
  tier2_cum_summary <- summarize_draw_rows(cum_tier2)

  tier_plot_df <- dplyr::bind_rows(
    tibble::tibble(year = forecast_years, median = tier1_cum_summary$median,
                    lower = tier1_cum_summary$lower, upper = tier1_cum_summary$upper,
                    tier = "Moderate PLE"),
    tibble::tibble(year = forecast_years, median = tier2_cum_summary$median,
                    lower = tier2_cum_summary$lower, upper = tier2_cum_summary$upper,
                    tier = "Catastrophic PLE"),
    tibble::tibble(year = forecast_years, median = cum_summary$median,
                    lower = cum_summary$lower, upper = cum_summary$upper,
                    tier = "Combined")
  )

  tier_colors <- c(
    "Moderate PLE" = "#00468B",
    "Catastrophic PLE" = "#ED0000",
    "Combined" = "#925E9F"
  )

  forecast_plot <- ggplot(tier_plot_df, aes(x = year, colour = tier, fill = tier)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.9) +
    scale_colour_manual(values = tier_colors) +
    scale_fill_manual(values = tier_colors) +
    scale_y_log10(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = c(1e5, 1e6, 1e7, 1e8, 5e8),
      limits = c(5e4, 5e8)
    ) +
    scale_x_continuous(breaks = seq(2030, 2100, by = 10)) +
    labs(x = "Year", y = "Cumulative deaths", colour = NULL, fill = NULL) +
    theme_nature() +
    theme(
      axis.title = element_text(family = "Arial", size = 9),
      axis.text = element_text(family = "Arial", size = 8),
      axis.ticks = element_line(linewidth = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = "Arial", size = 9),
      plot.background = element_rect(fill = "white", color = NA)
    )

  # Per-tier rate plot from rate_data
  rate_plot <- NULL
  if (!is.null(tier1$rate_data) && !is.null(tier2$rate_data)) {
    rate_plot <- plot_tiered_rates(tier1$rate_data, tier2$rate_data,
                                   conf = conf, output_dir = output_dir)
  }

  # SFA-based scenario bounds
  scenario_results <- NULL
  if (isTRUE(conf$sfa_bounds) && !is.null(tier1$rate_data) && !is.null(tier2$rate_data)) {
    window_size <- conf$window$default_size %||% 15

    # COLS-style SFA bounds: estimate frontier shift from GLM residuals.
    # The bounds are derived from the same GLM that produces the fitted trend,
    # so lower <= fitted <= upper by construction (no crossing possible).
    sfa_bounds_list <- list()
    for (tier_info in list(
      list(rd = tier1$rate_data, sev = bulk_severity_draws, label = "tier1"),
      list(rd = tier2$rate_data, sev = tail_severity_draws, label = "tier2")
    )) {
      hist_rd <- tier_info$rd %>% dplyr::filter(type == "historical")
      bounds <- tryCatch(
        estimate_sfa_bounds(hist_rd$observed_rate, hist_rd$fitted_rate, window_size),
        error = function(e) { warning("SFA bounds failed for ", tier_info$label, ": ", e$message); NULL }
      )
      if (!is.null(bounds)) {
        cat("  ", tier_info$label, "SFA shift: E[u] =", round(bounds$E_u, 4),
            " upper_mult =", round(bounds$upper_multiplier, 3),
            " lower_mult =", round(bounds$lower_multiplier, 3),
            " skewness =", round(bounds$residual_skewness, 3), "\n")
      }
      sfa_bounds_list[[tier_info$label]] <- list(bounds = bounds, severity = tier_info$sev)
    }

    # Build rate matrices and compute scenario burdens per tier.
    # The frontier rates are the GLM fitted rates shifted by the COLS multiplier,
    # so upper = fitted * exp(E[u]) and lower = fitted * exp(-E[u]).
    run_sfa_scenario <- function(direction) {
      tier_rate_sources <- list(
        tier1 = list(rd = tier1$rate_data),
        tier2 = list(rd = tier2$rate_data)
      )
      tier_draws_list <- list()
      for (tname in c("tier1", "tier2")) {
        bounds <- sfa_bounds_list[[tname]]$bounds
        sev <- sfa_bounds_list[[tname]]$severity
        if (is.null(bounds)) next
        rd <- tier_rate_sources[[tname]]$rd
        forecast_rd <- rd %>% dplyr::filter(year %in% forecast_years)
        fitted_at_forecast <- approx(rd$year, rd$fitted_rate,
                                      xout = forecast_years, rule = 2)$y
        frontier_rates <- sfa_bound_predict(fitted_at_forecast, bounds, direction)
        n_sev <- length(sev)
        rate_matrix <- matrix(rep(frontier_rates, n_sev), nrow = future_yrs, ncol = n_sev)
        tier_burden <- calc_burden(
          data = data, var = variable, thresh = lower_cutoff,
          start = base_year, end = current_year, future_yrs = future_yrs,
          shadow_mean = median(sev),
          shadow_low = quantile(sev, 0.025), shadow_up = quantile(sev, 0.975),
          use_time_trend = FALSE,
          fixed_rate = rate_matrix, severity_draws = sev, conf = conf,
          return_yearly_deaths_adjusted_draws = TRUE, validation = validation
        )
        tier_draws_list[[tname]] <- tier_burden$yearly_deaths_adjusted_draws
      }
      if (length(tier_draws_list) == 0) return(NULL)
      draws_list <- lapply(tier_draws_list, function(d) d)
      nc <- min(sapply(draws_list, ncol))
      total <- Reduce("+", lapply(draws_list, function(d) d[, seq_len(nc), drop = FALSE]))
      cum_total <- apply(total, 2, cumsum)
      cum_final <- cum_total[nrow(cum_total), ]
      list(
        yearly_deaths_adjusted_draws = total,
        cum_deaths = stats::median(cum_final),
        cum_deaths_low = as.numeric(stats::quantile(cum_final, 0.025)),
        cum_deaths_up = as.numeric(stats::quantile(cum_final, 0.975))
      )
    }

    worst_case <- tryCatch(run_sfa_scenario("upper"), error = function(e) {
      warning("SFA worst case failed: ", e$message); NULL })
    best_case <- tryCatch(run_sfa_scenario("lower"), error = function(e) {
      warning("SFA best case failed: ", e$message); NULL })

    # Per-draw tier-rate scenario summaries (probability UIs).
    # Re-use the GLM rate draw matrix `tier$rate_draws` (future_yrs x n_draws)
    # and apply the SFA point multiplier per direction so that 95% UIs
    # propagate GLM coefficient uncertainty scaled by exp(+/-E[u]).
    summarize_rate_per_draw <- function(mat) {
      per_draw <- colMeans(mat)
      list(
        median = stats::median(per_draw),
        low    = as.numeric(stats::quantile(per_draw, 0.025)),
        upp    = as.numeric(stats::quantile(per_draw, 0.975))
      )
    }
    build_tier_rate_scenarios <- function(tier_obj, bounds) {
      if (is.null(tier_obj) || is.null(tier_obj$rate_draws) || is.null(bounds)) return(NULL)
      m <- tier_obj$rate_draws
      list(
        reference       = summarize_rate_per_draw(m),
        best_plausible  = summarize_rate_per_draw(m * bounds$lower_multiplier),
        worst_plausible = summarize_rate_per_draw(m * bounds$upper_multiplier)
      )
    }
    tier_rate_scenarios <- list(
      tier1 = build_tier_rate_scenarios(tier1, sfa_bounds_list$tier1$bounds),
      tier2 = build_tier_rate_scenarios(tier2, sfa_bounds_list$tier2$bounds)
    )

    if (!is.null(worst_case) && !is.null(best_case)) {
      gap_draws <- colSums(worst_case$yearly_deaths_adjusted_draws) -
                   colSums(best_case$yearly_deaths_adjusted_draws)
      scenario_results <- list(
        worst_case = list(cum_deaths = worst_case$cum_deaths,
                          cum_deaths_low = worst_case$cum_deaths_low,
                          cum_deaths_up = worst_case$cum_deaths_up),
        best_case = list(cum_deaths = best_case$cum_deaths,
                         cum_deaths_low = best_case$cum_deaths_low,
                         cum_deaths_up = best_case$cum_deaths_up),
        deaths_gap = list(
          absolute = list(median = median(gap_draws),
                          lower_95 = as.numeric(quantile(gap_draws, 0.025)),
                          upper_95 = as.numeric(quantile(gap_draws, 0.975)))
        ),
        sfa_bounds = sfa_bounds_list,
        tier_rate_scenarios = tier_rate_scenarios
      )
      cat("SFA Worst case cumulative:", worst_case$cum_deaths / 1e6, "M\n")
      cat("SFA Best case cumulative:", best_case$cum_deaths / 1e6, "M\n")
      cat("SFA Gap (median):", median(gap_draws) / 1e6, "M\n")
      # Save scenario draw vectors for exceedance analysis
      if (!is.null(output_dir)) {
        worst_cum <- colSums(worst_case$yearly_deaths_adjusted_draws)
        best_cum <- colSums(best_case$yearly_deaths_adjusted_draws)
        saveRDS(worst_cum, file.path(output_dir, "cum_draws_worst.rds"))
        saveRDS(best_cum, file.path(output_dir, "cum_draws_best.rds"))
        cat("Saved SFA scenario draws to .rds\n")
      }
    }

    # Generate combined SFA frontier plot
    sfa_tier_list <- list()
    for (tier_info in list(
      list(rd = tier1$rate_data, label = "Moderate PLE", key = "tier1"),
      list(rd = tier2$rate_data, label = "Catastrophic PLE", key = "tier2")
    )) {
      b <- sfa_bounds_list[[tier_info$key]]$bounds
      if (!is.null(b)) {
        sfa_tier_list <- c(sfa_tier_list, list(list(
          rate_data = tier_info$rd, bounds = b,
          label = tier_info$label
        )))
      }
    }
    if (length(sfa_tier_list) > 0) {
      sfa_plot <- plot_sfa_frontiers(sfa_tier_list,
        forecast_years = forecast_years, current_year = current_year,
        conf = conf, output_dir = output_dir)
      tryCatch(print(sfa_plot), error = function(e) NULL)
    }
  }

  save_plot_if_enabled(forecast_plot, "tiered_cumulative_forecast", output_dir, conf)

  return_list <- list(
    result = result,
    plot = forecast_plot,
    rate_plot = rate_plot,
    forecast = data_forecast,
    tier1 = tier1,
    tier2 = tier2,
    yearly_deaths_adjusted_draws = total_adjusted_draws,
    scenarios = scenario_results
  )

  if (!is.null(tier1$rate_data)) return_list$tier1_rate_data <- tier1$rate_data
  if (!is.null(tier2$rate_data)) return_list$tier2_rate_data <- tier2$rate_data
  if (exists("worst_case") && !is.null(worst_case)) return_list$worst_case_draws <- worst_case$yearly_deaths_adjusted_draws
  if (exists("best_case") && !is.null(best_case)) return_list$best_case_draws <- best_case$yearly_deaths_adjusted_draws

  return(return_list)
}

#' Calculate Expected Annual and Total Burden with Uncertainty Intervals
#'
#' @param data A data frame containing the event data.
#' @param var The variable of interest for which the burden is calculated.
#' @param thresh The threshold value for the variable of interest.
#' @param start The start year of the time period for analysis.
#' @param end The end year of the time period for analysis.
#' @param future_yrs The number of future years for which to forecast the burden.
#' @param shadow_mean The mean shadow value for the burden calculation.
#' @param shadow_low The lower bound shadow value for the burden calculation.
#' @param shadow_up The upper bound shadow value for the burden calculation.
#' @param use_time_trend Logical, if TRUE, fit the temporal rate model; if FALSE, use `fixed_rate`.
#' @param fixed_rate Matrix of annual rates (forecast years x draws) used when `use_time_trend` is FALSE (scenario analysis).
#' @param severity_draws Vector of severity draws for uncertainty propagation.
#' @param conf Configuration parameters from config.yml.
#' @param return_rate_draws Logical, if TRUE, returns the rate draws matrix (default: FALSE).
#' @param return_yearly_deaths_draws Logical, if TRUE, returns the unadjusted yearly deaths draws (default: FALSE).
#' @param return_yearly_deaths_adjusted_draws Logical, if TRUE, returns the population-adjusted yearly deaths draws (default: TRUE).
#'
#' @return A list containing:
#'   \item{result}{Tibble with burden summary statistics}
#'   \item{rate_data}{Rate data if use_time_trend=TRUE}
#'   \item{yearly_deaths_adjusted_draws}{Matrix of population-adjusted yearly deaths draws if return_yearly_deaths_adjusted_draws=TRUE}
#'   \item{yearly_deaths_draws}{Matrix of unadjusted yearly deaths draws if return_yearly_deaths_draws=TRUE}
#'   \item{rate_draws}{Matrix of rate draws if return_rate_draws=TRUE}
#'   \item{forecast}{Yearly forecast summary if return_forecast=TRUE}
#' @export
calc_burden <- function(
    data, var, thresh, upper_thresh = NULL, start, end, future_yrs,
    shadow_mean, shadow_low, shadow_up, use_time_trend = FALSE,
    fixed_rate = NULL, severity_draws, conf,
    return_rate_draws = FALSE,
    return_yearly_deaths_draws = FALSE,
    return_yearly_deaths_adjusted_draws = TRUE,
    return_forecast = FALSE,
    validation = FALSE,
    window_sizes = NULL) {
  # Filter data for the time period
  data_period <- data %>%
    dplyr::filter(start_year >= start, end_year <= end)

  n_draws <- conf$draws

  # Get population forecast data
  pop_forecast <- get_pop_forecast() %>% dplyr::select(1, 2)
  pop_draws <- get_pop_draws()

  # Annual rates: either supplied (scenario analysis) or estimated from the temporal model
  if (!use_time_trend) {
    # Year-specific rates supplied as a matrix (n_years x n_draws)
    actual_n_draws <- length(severity_draws)
    if (ncol(fixed_rate) != actual_n_draws) {
      idx <- sample(seq_len(ncol(fixed_rate)), actual_n_draws, replace = TRUE)
      fixed_rate <- fixed_rate[, idx, drop = FALSE]
    }
    rate_draws <- fixed_rate
    lambda <- mean(rate_draws)

    yearly_deaths_draws <- t(t(rate_draws) * severity_draws)

    idx_pop <- sample(seq_len(ncol(pop_draws) - 1), actual_n_draws, replace = TRUE)
    pop_draws_filtered <- pop_draws |>
      dplyr::filter(year > end) |>
      dplyr::select(-year) |>
      as.matrix()
    pop_draws_filtered <- pop_draws_filtered[, idx_pop]

    yearly_deaths_adjusted_draws <- yearly_deaths_draws * pop_draws_filtered / conf$pop_reference

    yearly_deaths_summary <- apply(
      yearly_deaths_adjusted_draws, MARGIN = 1,
      function(x) { c(median(x), quantile(x, c(0.025, 0.975))) }
    ) |> t() |> as_tibble(.name_repair = ~ c("mean", "lower", "upper")) |>
      mutate(year = (end + 1):(end + future_yrs), .before = mean)

    cum_deaths_draws <- apply(yearly_deaths_adjusted_draws, 2, cumsum)
    cum_draws_vec <- cum_deaths_draws[nrow(cum_deaths_draws), ]
    cum_deaths_median <- stats::median(cum_draws_vec)
    cum_deaths_low <- as.numeric(stats::quantile(cum_draws_vec, 0.025))
    cum_deaths_up <- as.numeric(stats::quantile(cum_draws_vec, 0.975))

    cum_deaths_summary <- t(apply(cum_deaths_draws, 1, function(x) {
      c(median = stats::median(x), lower = as.numeric(stats::quantile(x, 0.025)),
        upper = as.numeric(stats::quantile(x, 0.975)))
    }))

    data_forecast <- tibble(
      year = (end + 1):(end + future_yrs), cum_years = 1:future_yrs, rate = lambda,
      pred_deaths = yearly_deaths_summary$mean,
      pred_deaths_low = yearly_deaths_summary$lower,
      pred_deaths_upp = yearly_deaths_summary$upper,
      cum_deaths = cum_deaths_summary[, "median"],
      cum_deaths_low = cum_deaths_summary[, "lower"],
      cum_deaths_up = cum_deaths_summary[, "upper"]
    )

    yearly_deaths <- lambda * shadow_mean
    yearly_deaths_low <- lambda * shadow_low
    yearly_deaths_up <- lambda * shadow_up

    # Return results
    result <- tibble(
      yearly_rate = lambda,
      yearly_deaths,
      yearly_deaths_low,
      yearly_deaths_up,
      cum_deaths = cum_deaths_median,
      cum_deaths_low = cum_deaths_low,
      cum_deaths_up = cum_deaths_up
    )
  } else {
    # Time-variable approach - use analyze_pandemic_trend function
    set.seed(42)
    trend_analysis <- analyze_pandemic_trend(
      data = data_period,
      variable = var,
      lower_cutoff = thresh,
      upper_threshold = upper_thresh,
      start_year = start,
      end_year = end,
      window_sizes = window_sizes,
      future_years = future_yrs,
      model_type = conf$time$model_type,
      n_draws = n_draws,
      generate_draws = TRUE,
      conf = conf
    )

    # Extract draws of prediction rates from time trend analysis
    rate_draws <- trend_analysis$draws$annual_rates

    # Get year indicator from time trend analysis
    years <- trend_analysis$draws$years

    # Filter only prediction period (years >= current year)
    current_year <- conf$time$current_year
    forecast_filter <- years > current_year
    rate_draws <- rate_draws[forecast_filter, , drop = FALSE]
    forecast_years <- years[forecast_filter]

    # Get actual number of draws after potential xi filtering
    actual_n_draws <- length(severity_draws)

    # Ensure rate_draws matches the number of severity draws
    if (ncol(rate_draws) != actual_n_draws) {
      # Resample columns to match severity_draws length
      idx_rate <- sample(seq_len(ncol(rate_draws)), actual_n_draws, replace = TRUE)
      rate_draws <- rate_draws[, idx_rate, drop = FALSE]
    }

    # Multiply rate and severity draws, propagating uncertainty
    # This creates a matrix where each column is a draw of yearly deaths
    yearly_deaths_draws <- t(t(rate_draws) * severity_draws)

    if (validation) {
      pop_forecast <- pop_forecast |>
        filter(year > end, year <= end + future_yrs) |>
        pull(val)
      pop_draws_filtered <- matrix(pop_forecast, nrow = length(pop_forecast), ncol = actual_n_draws)
    } else if (!validation) {
      # Resample population draws to match the actual number of severity draws
      idx_pop <- sample(seq_len(ncol(pop_draws) - 1), actual_n_draws, replace = TRUE)
      pop_draws_filtered <- pop_draws |>
        filter(year > current_year) |>
        dplyr::select(-year) |>
        as.matrix()
      pop_draws_filtered <- pop_draws_filtered[, idx_pop]
    }

    # Estimate the draw-level population adjusted future yearly deaths
    # This accounts for changing global population over the forecast period
    yearly_deaths_adjusted_draws <- yearly_deaths_draws * pop_draws_filtered / conf$pop_reference

    # Compute cumulative deaths at draw level (correct uncertainty propagation)
    # Each column is a draw; cumsum within each draw preserves correlation structure
    cum_deaths_draws <- apply(yearly_deaths_adjusted_draws, 2, cumsum)

    # Final totals (sum across all years)
    cum_draws_vec <- cum_deaths_draws[nrow(cum_deaths_draws), ]
    cum_deaths_median <- stats::median(cum_draws_vec)
    cum_deaths_low <- as.numeric(stats::quantile(cum_draws_vec, 0.025))
    cum_deaths_up <- as.numeric(stats::quantile(cum_draws_vec, 0.975))

    # Year-by-year cumulative summary (for plotting)
    cum_deaths_summary <- t(apply(cum_deaths_draws, 1, function(x) {
      c(
        median = stats::median(x),
        lower = as.numeric(stats::quantile(x, 0.025)),
        upper = as.numeric(stats::quantile(x, 0.975))
      )
    }))

    # Calculate summary statistics for forecast rates
    forecast_rates <- apply(
      rate_draws,
      MARGIN = 1,
      function(x) {
        mean_val <- mean(x)
        quant <- quantile(x, prob = c(0.025, 0.975))
        c(mean_val, quant)
      }
    ) |>
      t() |>
      as_tibble(.name_repair = ~ c("fitted_rate", "lower_95", "upper_95")) |>
      mutate(year = forecast_years, .before = fitted_rate)

    # Calculate summary statistics for yearly deaths
    yearly_deaths_summary <- apply(
      yearly_deaths_adjusted_draws,
      MARGIN = 1,
      function(x) {
        # Using median for central tendency as deaths may be skewed
        mean_val <- median(x)
        quant <- quantile(x, prob = c(0.025, 0.975))
        c(mean_val, quant)
      }
    ) |>
      t() |>
      as_tibble(.name_repair = ~ c("mean", "lower", "upper")) |>
      mutate(year = forecast_years, .before = mean)


    # Calculate year-by-year burden using time-variable rates
    forecast_burden <- forecast_rates %>%
      mutate(
        yearly_deaths = yearly_deaths_summary$mean,
        yearly_deaths_low = yearly_deaths_summary$lower,
        yearly_deaths_up = yearly_deaths_summary$upper
      ) %>%
      mutate(
        # Scale by population
        pred_deaths = yearly_deaths_summary$mean,
        pred_deaths_low = yearly_deaths_summary$lower,
        pred_deaths_upp = yearly_deaths_summary$upper
      )

    # Calculate cumulative values using percentiles-of-cumsums (correct approach)
    # This properly propagates uncertainty by computing cumsum per draw first,
    # then taking percentiles across draws at each time point
    data_forecast <- forecast_burden %>%
      mutate(
        cum_deaths = cum_deaths_summary[, "median"],
        cum_deaths_low = cum_deaths_summary[, "lower"],
        cum_deaths_up = cum_deaths_summary[, "upper"],
        cum_years = seq_len(n())
      )

    # Average rate for summary statistics
    mean_rate <- mean(forecast_rates$fitted_rate)
    mean_rate_low <- mean(forecast_rates$lower_95)
    mean_rate_up <- mean(forecast_rates$upper_95)

    # Extract specific year rates
    # First year of forecast (typically 2026)
    first_year_rate <- forecast_rates$fitted_rate[1]
    first_year <- forecast_rates$year[1]
    
    # Extract uncertainty bounds for first year
    first_year_rate_lower <- forecast_rates$lower_95[1]
    first_year_rate_upper <- forecast_rates$upper_95[1]

    # Last year of forecast (typically 2100)
    last_year_rate <- forecast_rates$fitted_rate[nrow(forecast_rates)]
    last_year <- forecast_rates$year[nrow(forecast_rates)]
    
    # Extract uncertainty bounds for last year
    last_year_rate_lower <- forecast_rates$lower_95[nrow(forecast_rates)]
    last_year_rate_upper <- forecast_rates$upper_95[nrow(forecast_rates)]

    # Extract model information from trend_analysis
    best_model_name <- trend_analysis$model_comparison$Model[1]

    # Extract statistical model coefficients and p-values
    best_model <- NULL
    beta_0_coefficient <- NA
    beta_0_se <- NA
    beta_0_pvalue <- NA
    beta_1_coefficient <- NA
    beta_1_se <- NA
    beta_1_pvalue <- NA
    beta_2_coefficient <- NA
    beta_2_se <- NA
    beta_2_pvalue <- NA
    beta_3_coefficient <- NA
    beta_3_se <- NA
    beta_3_pvalue <- NA
    breakpoint_year <- NA

    # Get the best model object directly from trend_analysis
    if (!is.null(trend_analysis$best_model)) {
      best_model <- trend_analysis$best_model
    }

    # Extract coefficient information if we have a model
    if (!is.null(best_model) && (inherits(best_model, "glm") || inherits(best_model, "negbin"))) {
      model_summary <- summary(best_model)

      # Extract coefficients
      coef_names <- rownames(model_summary$coefficients)

      # Look for intercept coefficient (beta_0)
      intercept_idx <- which(coef_names == "(Intercept)")
      if (length(intercept_idx) > 0) {
        beta_0_coefficient <- model_summary$coefficients[intercept_idx, "Estimate"]
        beta_0_se <- model_summary$coefficients[intercept_idx, "Std. Error"]
        beta_0_pvalue <- model_summary$coefficients[intercept_idx, "Pr(>|z|)"]
      }

      # Look for linear trend coefficient (beta_1)
      # This could be named "year_scaled", "year_centered", or similar
      trend_coef_idx <- which(grepl("year", coef_names, ignore.case = TRUE) & !grepl("squared|quad", coef_names, ignore.case = TRUE))

      if (length(trend_coef_idx) > 0) {
        # Take the first linear year coefficient
        idx <- trend_coef_idx[1]
        beta_1_coefficient <- model_summary$coefficients[idx, "Estimate"]
        beta_1_se <- model_summary$coefficients[idx, "Std. Error"]
        beta_1_pvalue <- model_summary$coefficients[idx, "Pr(>|z|)"]
      }

      # Look for piecewise breakpoint shift coefficient (beta_2)
      # This is named "post_breakpoint" in piecewise models
      breakpoint_coef_idx <- which(coef_names == "post_breakpoint")
      if (length(breakpoint_coef_idx) > 0) {
        beta_2_coefficient <- model_summary$coefficients[breakpoint_coef_idx, "Estimate"]
        beta_2_se <- model_summary$coefficients[breakpoint_coef_idx, "Std. Error"]
        beta_2_pvalue <- model_summary$coefficients[breakpoint_coef_idx, "Pr(>|z|)"]
        breakpoint_year <- if (!is.null(conf$time$breakpoint_year)) conf$time$breakpoint_year else 1950
      }

      # Look for piecewise slope change coefficient (beta_3)
      # This is named "year_scaled_x_post" in full piecewise models
      slope_change_idx <- which(coef_names == "year_scaled_x_post")
      if (length(slope_change_idx) > 0) {
        beta_3_coefficient <- model_summary$coefficients[slope_change_idx, "Estimate"]
        beta_3_se <- model_summary$coefficients[slope_change_idx, "Std. Error"]
        beta_3_pvalue <- model_summary$coefficients[slope_change_idx, "Pr(>|z|)"]
      }
    }

    # Annualised trend coefficient: exp(beta_1 / historical span) when beta_1 is available; the annualised last/first forecast-rate ratio only as fallback
    # This gives us a measure of how much the rate changes over the forecast period
    trend_coefficient <- NA
    trend_coefficient_lower <- NA
    trend_coefficient_upper <- NA
    
    # First try to calculate from β₁ coefficient (parameter uncertainty only)
    # This avoids using prediction intervals that include random variation
    if (!is.na(beta_1_coefficient) && !is.na(beta_1_se) && beta_1_se > 0) {
      # The model uses year_scaled which goes from 0 to 1 over the historical period
      # Try to get actual historical span from config or use reasonable default
      historical_years_span <- if (!is.null(conf$time$current_year) && !is.null(conf$base_year)) {
        conf$time$current_year - conf$base_year
      } else {
        425  # Default assumption: 1600-2025
      }
      
      # Calculate annual rate change from β₁ coefficient directly
      # β₁ represents change in log(rate) as year_scaled increases by 1 (full historical span)
      # Annual change in log space = β₁ / historical_years_span
      annual_log_change <- beta_1_coefficient / historical_years_span
      trend_coefficient <- exp(annual_log_change)
      
      # Calculate confidence bounds using β₁ standard error (not prediction intervals)
      beta_1_lower <- beta_1_coefficient - 1.96 * beta_1_se
      beta_1_upper <- beta_1_coefficient + 1.96 * beta_1_se
      
      annual_log_change_lower <- beta_1_lower / historical_years_span
      annual_log_change_upper <- beta_1_upper / historical_years_span
      
      trend_coefficient_lower <- exp(annual_log_change_lower)
      trend_coefficient_upper <- exp(annual_log_change_upper)
      
    } else if (first_year_rate > 0) {
      # Fallback: Use forecast rate ratios (includes prediction uncertainty)
      # This method includes random variation and gives wider intervals
      total_change_ratio <- last_year_rate / first_year_rate

      # Annualized trend coefficient (geometric mean of growth)
      years_span <- last_year - first_year
      if (years_span > 0) {
        trend_coefficient <- total_change_ratio^(1 / years_span)
        
        # Calculate uncertainty bounds using prediction intervals
        # Note: This includes random variation, making intervals wider
        if (first_year_rate_upper > 0 && first_year_rate_lower > 0) {
          total_change_ratio_lower <- last_year_rate_lower / first_year_rate_upper
          total_change_ratio_upper <- last_year_rate_upper / first_year_rate_lower
          trend_coefficient_lower <- total_change_ratio_lower^(1 / years_span)
          trend_coefficient_upper <- total_change_ratio_upper^(1 / years_span)
        }
      } else {
        trend_coefficient <- 1.0
        trend_coefficient_lower <- 1.0
        trend_coefficient_upper <- 1.0
      }
    }

    # For constant rate models, the coefficient should be exactly 1
    if (best_model_name == "Constant Rate") {
      trend_coefficient <- 1.0
      trend_coefficient_lower <- 1.0
      trend_coefficient_upper <- 1.0
      beta_1_coefficient <- 0.0 # No trend for constant rate
      beta_1_se <- 0.0
      beta_1_pvalue <- 1.0 # Not significant
    }

    # Calculate rate uncertainty bounds as mean of forecast period uncertainty
    # This provides uncertainty intervals for the mean annual exceedance rate
    mean_rate_lower_95 <- mean(forecast_rates$lower_95)
    mean_rate_upper_95 <- mean(forecast_rates$upper_95)

    # Create detailed yearly_rate information
    yearly_rate_detailed <- list(
      mean = mean_rate,
      lower_95 = mean_rate_lower_95,
      upper_95 = mean_rate_upper_95,
      rate_2026 = if (first_year == 2026) first_year_rate else NA,
      rate_2100 = if (last_year == 2100) last_year_rate else NA,
      first_year_rate = first_year_rate,
      first_year = first_year,
      last_year_rate = last_year_rate,
      last_year = last_year,
      trend_coefficient = trend_coefficient,
      trend_coefficient_lower = trend_coefficient_lower,
      trend_coefficient_upper = trend_coefficient_upper,
      model_type = best_model_name,
      beta_0_coefficient = beta_0_coefficient,
      beta_0_se = beta_0_se,
      beta_0_pvalue = beta_0_pvalue,
      beta_1_coefficient = beta_1_coefficient,
      beta_1_se = beta_1_se,
      beta_1_pvalue = beta_1_pvalue,
      beta_2_coefficient = beta_2_coefficient,
      beta_2_se = beta_2_se,
      beta_2_pvalue = beta_2_pvalue,
      beta_3_coefficient = beta_3_coefficient,
      beta_3_se = beta_3_se,
      beta_3_pvalue = beta_3_pvalue,
      breakpoint_year = breakpoint_year
    )

    # Return results
    result <- tibble(
      yearly_rate = yearly_rate_detailed,
      yearly_deaths = mean_rate * median(severity_draws),
      yearly_deaths_low = mean_rate * quantile(severity_draws, prob = 0.025),
      yearly_deaths_up = mean_rate * quantile(severity_draws, prob = 0.975),
      cum_deaths = cum_deaths_median,
      cum_deaths_low = cum_deaths_low,
      cum_deaths_up = cum_deaths_up,
      dynamic_rates = TRUE
    )
  }

  return_list <- list(result = result)

  if (use_time_trend) {
    return_list$rate_data <- trend_analysis$exceedance_rate_df
  }
  if (return_yearly_deaths_adjusted_draws && !is.null(yearly_deaths_adjusted_draws)) {
    return_list$yearly_deaths_adjusted_draws <- yearly_deaths_adjusted_draws
  }
  if (return_yearly_deaths_draws && !is.null(yearly_deaths_draws)) {
    return_list$yearly_deaths_draws <- yearly_deaths_draws
  }
  if (return_rate_draws && !is.null(rate_draws)) {
    return_list$rate_draws <- rate_draws
  }
  if (return_forecast) {
    return_list$forecast <- data_forecast
  }

  return(return_list)
}
