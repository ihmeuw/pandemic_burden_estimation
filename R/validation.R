#' Run Model Validation
#'
#' Performs out-of-sample validation by training on pre-1950 data and validating on post-1950
#'
#' @param data Complete dataset
#' @param variable Variable to analyze
#' @param conf Configuration object
#' @return List containing validation results and plots
validate_model <- function(data, variable, conf, output_dir = NULL) {
  train_data <- data %>%
    filter(
      start_year >= conf$train_year_start,
      start_year <= conf$train_year_end
    )

  # The validation profile disables burden calculation inside run_analysis();
  # the expected burden is computed separately below.
  train_results <- run_analysis(
    train_data,
    variable,
    conf = conf,
    year = conf$train_year_start,
    return_severity_draws = TRUE
  )

  expected_burden <- calc_tiered_burden(
    data = train_data,
    variable = variable,
    conf = conf,
    bulk_severity_draws = train_results$param_draws$bulk_mean,
    tail_severity_draws = train_results$param_draws$tail_mean,
    shadow_mean = train_results$results$boot,
    shadow_low = train_results$results$boot_low,
    shadow_up = train_results$results$boot_upp,
    output_dir = output_dir,
    validation = TRUE,
    window_sizes = conf$window_size,
    return_severity_draws = FALSE
  )

  # Build forecast and plot from the combined tier draws
  if (!is.null(expected_burden$yearly_deaths_adjusted_draws)) {
    adj_draws <- expected_burden$yearly_deaths_adjusted_draws
    cum_draws <- apply(adj_draws, 2, cumsum)
    cum_summary <- summarize_draw_rows(cum_draws)
    yearly_summary <- summarize_draw_rows(adj_draws)
    forecast_years <- seq(conf$train_year_end + 1,
                          conf$train_year_end + nrow(adj_draws))
    expected_burden$forecast <- tibble::tibble(
      year = forecast_years,
      cum_deaths = cum_summary$median,
      cum_deaths_low = cum_summary$lower,
      cum_deaths_up = cum_summary$upper,
      pred_deaths = yearly_summary$median,
      pred_deaths_low = yearly_summary$lower,
      pred_deaths_upp = yearly_summary$upper
    )
    expected_burden$plot <- ggplot(expected_burden$forecast, aes(x = year)) +
      geom_line(aes(y = cum_deaths, colour = "Predicted"), linewidth = 1) +
      geom_ribbon(aes(ymin = cum_deaths_low, ymax = cum_deaths_up,
                       fill = "Predicted"), alpha = 0.2) +
      scale_y_log10(
        labels = scales::label_number(scale_cut = scales::cut_short_scale()),
        breaks = c(1e5, 1e6, 1e7, 1e8), limits = c(5e4, 4e8)
      ) +
      labs(x = "Year", y = "Cumulative Deaths") +
      theme_nature()
  }

  forecast_plot <- expected_burden$plot

  # Process observed data
  observed_cumulative <- process_observed_data(
    data,
    variable,
    start = conf$test_year_start,
    end = conf$test_year_end,
    obs_years = conf$observation_years
  )

  validation_plot <- create_validation_plot(
    forecast_plot,
    observed_cumulative,
    conf,
    variable
  )

  # --- SFA Three-Scenario Validation (additional output) ---
  sfa_validation_plot <- NULL
  sfa_diagnostic_plot <- NULL
  sfa_bounds_result <- NULL

  if (!is.null(expected_burden$tier1_rate_data) &&
      !is.null(expected_burden$tier2_rate_data)) {

    cat("\n=== Running SFA Three-Scenario Validation ===\n")

    sfa_results <- run_sfa_validation(
      expected_burden = expected_burden,
      bulk_severity_draws = train_results$param_draws$bulk_mean,
      tail_severity_draws = train_results$param_draws$tail_mean,
      conf = conf
    )

    if (!is.null(sfa_results) &&
        !is.null(sfa_results$lower) &&
        !is.null(sfa_results$fitted) &&
        !is.null(sfa_results$upper)) {

      sfa_validation_plot <- create_sfa_validation_plot(
        scenario_forecasts = sfa_results,
        observed_data = observed_cumulative,
        conf = conf
      )

      sfa_tier_list <- list()
      for (tier_info in list(
        list(rd = expected_burden$tier1_rate_data, label = "Moderate PLE", key = "tier1"),
        list(rd = expected_burden$tier2_rate_data, label = "Catastrophic PLE", key = "tier2")
      )) {
        b <- sfa_results$sfa_bounds[[tier_info$key]]$bounds
        if (!is.null(b)) {
          sfa_tier_list <- c(sfa_tier_list, list(list(
            rate_data = tier_info$rd, bounds = b,
            label = tier_info$label
          )))
        }
      }
      if (length(sfa_tier_list) > 0) {
        forecast_years_diag <- seq(conf$test_year_start + 1, conf$test_year_end)
        sfa_diagnostic_plot <- plot_sfa_frontiers(
          sfa_tier_list,
          forecast_years = forecast_years_diag,
          current_year = conf$train_year_end,
          conf = conf,
          output_dir = NULL
        )
      }

      sfa_bounds_result <- sfa_results$sfa_bounds
    }
  }

  return(list(
    model = train_results,
    plot = validation_plot,
    observed = observed_cumulative,
    predicted = expected_burden,
    forecast = expected_burden$forecast |>
      filter(year %in% observed_cumulative$year) |>
      dplyr::select(year, cum_deaths:cum_deaths_up),
    sfa_validation_plot = sfa_validation_plot,
    sfa_diagnostic_plot = sfa_diagnostic_plot,
    sfa_bounds = sfa_bounds_result
  ))
}

#' Process Observed Data for Validation
#'
#' @param data Complete dataset
#' @param variable Variable name to analyze
#' @param start Beginning of test period
#' @param end End of test period
#' @param obs_years Years to extract for comparison
#' @return Dataframe with cumulative observed values
process_observed_data <- function(data, variable, start, end, obs_years) {
  # Filter to test period
  observed_data <- data %>%
    filter(start_year >= start)

  pop <- get_pop_forecast() |>
    filter(year <= 2025) |>
    dplyr::select(year, pop = val)

  # Process observed data to get cumulative deaths by year
  observed_expanded <- observed_data %>%
    dplyr::select(name, start_year, end_year, deaths) |>
    mutate(length = end_year - start_year + 1) |>
    # Expand each event to yearly rows
    mutate(year = map2(start_year, end_year, seq)) %>%
    unnest(year) %>%
    # Calculate deaths per year (assuming uniform distribution within event)
    mutate(deaths_per_year = deaths / length) %>%
    dplyr::select(year, name, deaths_per_year) %>%
    # Transform to wide format with events as columns
    pivot_wider(names_from = name, values_from = deaths_per_year, values_fill = 0) |>
    arrange(year)

  # Create full yearly dataset and fill in missing years
  observed_yearly <- tibble(year = start:end) %>%
    left_join(observed_expanded, by = "year") %>%
    replace(is.na(.), 0)

  # Calculate cumulative deaths by year
  observed_cumulative <- observed_yearly %>%
    pivot_longer(-year) %>%
    group_by(year) %>%
    summarize(deaths = sum(value)) %>%
    mutate(obs = cumsum(deaths)) %>%
    dplyr::select(year, obs) |>
    # Filter to specific observation years for comparison
    filter(year %in% obs_years)

  return(observed_cumulative)
}

#' Create Enhanced Validation Plot
#'
#' @param forecast_plot Base forecast plot from burden calculation
#' @param observed_data Observed cumulative data
#' @param conf Configuration object
#' @param variable Variable name used in analysis
#' @return ggplot object with an enhanced validation plot
create_validation_plot <- function(forecast_plot, observed_data, conf, variable) {
  # Palette shared with the other figures: blue for the prediction (as the fitted
  # trend in plot_tiered_rates / plot_sfa_frontiers), red for the observed series.
  plot_palette <- c(
    "Predicted" = "#0072B2",
    "Observed" = "#C0392B"
  )

  # Clamp ribbon lower bound in the base plot data to avoid cutoff on log scale
  y_floor <- 5e4
  validation_plot <- forecast_plot
  validation_plot$data <- validation_plot$data %>%
    dplyr::mutate(cum_deaths_low = pmax(cum_deaths_low, y_floor))

  validation_plot <- validation_plot +
    # Add observed data points
    geom_point(
      data = observed_data,
      aes(x = year, y = obs, colour = "Observed"),
      shape = 21, # Filled circle with border
      size = 3,
      stroke = 0.5,
      fill = plot_palette["Observed"],
      alpha = 0.7 # More transparent
    ) +
    # Add connecting line for observed points
    geom_line(
      data = observed_data,
      aes(x = year, y = obs, colour = "Observed"),
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.7 # More transparent
    )

  validation_plot <- validation_plot +
    # Colour scales
    scale_color_manual(
      values = plot_palette,
      name = "Data Source",
      labels = c(
        "Predicted" = "Model prediction",
        "Observed" = "Historical data"
      )
    ) +
    scale_fill_manual(
      values = alpha(plot_palette["Predicted"], 0.3), # More transparent fill for uncertainty interval
      name = "Uncertainty",
      labels = c("Predicted" = "95% uncertainty interval")
    ) +
    # Conditionally add labels based on config
    labs(
      title = if (!is.null(conf$plots$validation_plot_labels) && conf$plots$validation_plot_labels) "Model Validation Results" else NULL,
      subtitle = if (!is.null(conf$plots$validation_plot_labels) && conf$plots$validation_plot_labels) {
        sprintf(
          "Training: %d–%d, Validation: %d–%d",
          conf$train_year_start,
          conf$train_year_end,
          conf$test_year_start + 1,
          conf$test_year_end
        )
      } else {
        NULL
      },
      caption = if (!is.null(conf$plots$validation_plot_labels) && conf$plots$validation_plot_labels) {
        sprintf(
          "Projection based on %d bootstrap samples",
          conf$draws
        )
      } else {
        NULL
      },
      x = "Year",
      y = "Cumulative pandemic-related deaths (log scale)",
      fill = NULL,
      color = "Source"
    ) +
    scale_y_log10(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = scales::breaks_log(5), limits = c(5e4, NA)
    ) +
    scale_x_continuous(
      breaks = seq(1950, 2025, by = 25)
    ) +
    theme_nature()

  return(validation_plot)
}

#' Run SFA Three-Scenario Validation
#'
#' Computes three cumulative death trajectories using SFA-derived frontier rates
#' (lower/fitted/upper) from training-only data, each with 95% UI.
#'
#' @param expected_burden Output from calc_tiered_burden(validation=TRUE)
#' @param bulk_severity_draws Numeric vector of Tier 1 severity draws
#' @param tail_severity_draws Numeric vector of Tier 2 severity draws
#' @param conf Validation configuration object
#' @return List with scenario forecasts (3 tibbles) and sfa_bounds (per-tier)
run_sfa_validation <- function(expected_burden, bulk_severity_draws,
                                tail_severity_draws, conf) {
  tier1_rd <- expected_burden$tier1_rate_data
  tier2_rd <- expected_burden$tier2_rate_data
  if (is.null(tier1_rd) || is.null(tier2_rd)) {
    warning("SFA validation requires per-tier rate_data from calc_tiered_burden")
    return(NULL)
  }

  window_size <- conf$window$default_size %||% conf$window_size %||% 15
  validation_years <- seq(conf$test_year_start + 1, conf$test_year_end)

  pop_forecast <- get_pop_forecast() %>%
    dplyr::filter(year %in% validation_years)

  sfa_bounds_list <- list()
  for (tier_info in list(
    list(rd = tier1_rd, sev = bulk_severity_draws, label = "tier1"),
    list(rd = tier2_rd, sev = tail_severity_draws, label = "tier2")
  )) {
    hist_rd <- tier_info$rd %>% dplyr::filter(type == "historical")
    bounds <- tryCatch(
      estimate_sfa_bounds(hist_rd$observed_rate, hist_rd$fitted_rate, window_size),
      error = function(e) {
        warning("SFA bounds failed for ", tier_info$label, ": ", e$message)
        NULL
      }
    )
    if (!is.null(bounds)) {
      cat("  ", tier_info$label, "SFA shift: E[u] =", round(bounds$E_u, 4),
          " upper_mult =", round(bounds$upper_multiplier, 3),
          " lower_mult =", round(bounds$lower_multiplier, 3), "\n")
    }
    sfa_bounds_list[[tier_info$label]] <- list(
      bounds = bounds, severity = tier_info$sev, rate_data = tier_info$rd
    )
  }

  compute_scenario_draws <- function(direction) {
    tier_draws_list <- list()
    for (tname in c("tier1", "tier2")) {
      bounds <- sfa_bounds_list[[tname]]$bounds
      sev <- sfa_bounds_list[[tname]]$severity
      rd <- sfa_bounds_list[[tname]]$rate_data
      if (is.null(bounds)) next

      fitted_at_val <- approx(rd$year, rd$fitted_rate,
                               xout = validation_years, rule = 2)$y
      if (direction == "fitted") {
        scenario_rates <- fitted_at_val
      } else {
        scenario_rates <- sfa_bound_predict(fitted_at_val, bounds, direction)
      }

      n_sev <- length(sev)
      rate_matrix <- matrix(rep(scenario_rates, n_sev),
                             nrow = length(validation_years), ncol = n_sev)
      annual_deaths <- t(t(rate_matrix) * sev)

      for (i in seq_along(validation_years)) {
        yr <- validation_years[i]
        pop_i <- pop_forecast$val[pop_forecast$year == yr]
        if (length(pop_i) == 0) pop_i <- conf$pop_reference
        annual_deaths[i, ] <- annual_deaths[i, ] * (pop_i / conf$pop_reference)
      }
      tier_draws_list[[tname]] <- annual_deaths
    }
    if (length(tier_draws_list) == 0) return(NULL)
    nc <- min(sapply(tier_draws_list, ncol))
    total <- Reduce("+", lapply(tier_draws_list, function(d) d[, seq_len(nc), drop = FALSE]))
    cum_total <- apply(total, 2, cumsum)
    tibble::tibble(
      year = validation_years,
      cum_deaths = apply(cum_total, 1, stats::median),
      cum_deaths_low = apply(cum_total, 1, function(x) stats::quantile(x, 0.025)),
      cum_deaths_up = apply(cum_total, 1, function(x) stats::quantile(x, 0.975))
    )
  }

  cat("\n=== SFA VALIDATION SCENARIOS ===\n")
  lower_forecast <- compute_scenario_draws("lower")
  fitted_forecast <- compute_scenario_draws("fitted")
  upper_forecast <- compute_scenario_draws("upper")

  if (!is.null(lower_forecast)) {
    cat("Lower frontier (best case) final:", tail(lower_forecast$cum_deaths, 1) / 1e6, "M\n")
  }
  if (!is.null(fitted_forecast)) {
    cat("Fitted trend (reference) final:", tail(fitted_forecast$cum_deaths, 1) / 1e6, "M\n")
  }
  if (!is.null(upper_forecast)) {
    cat("Upper frontier (worst case) final:", tail(upper_forecast$cum_deaths, 1) / 1e6, "M\n")
  }

  list(
    lower = lower_forecast,
    fitted = fitted_forecast,
    upper = upper_forecast,
    sfa_bounds = sfa_bounds_list
  )
}

#' Create SFA Three-Scenario Validation Plot (Forest Plot)
#'
#' Renders a horizontal point-range comparison of three SFA scenario cumulative
#' deaths (at end of validation period) with 95% UI, alongside observed values.
#'
#' @param scenario_forecasts List with lower, fitted, upper tibbles
#' @param observed_data Observed cumulative deaths (all pandemics)
#' @param conf Validation configuration
#' @return ggplot object
create_sfa_validation_plot <- function(scenario_forecasts, observed_data, conf = NULL) {
  final_lower <- tail(scenario_forecasts$lower, 1)
  final_fitted <- tail(scenario_forecasts$fitted, 1)
  final_upper <- tail(scenario_forecasts$upper, 1)

  obs_final <- max(observed_data$obs)

  plot_df <- tibble::tibble(
    label = c("Plausible-worst",
              "Fitted trend",
              "Plausible-best"),
    median = c(final_upper$cum_deaths, final_fitted$cum_deaths, final_lower$cum_deaths),
    lower = c(final_upper$cum_deaths_low, final_fitted$cum_deaths_low, final_lower$cum_deaths_low),
    upper = c(final_upper$cum_deaths_up, final_fitted$cum_deaths_up, final_lower$cum_deaths_up),
    type = "Scenario"
  )
  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))

  # Palette harmonised with plot_sfa_frontiers.
  scenario_colors <- c(
    "Plausible-worst" = "#C0392B",
    "Fitted trend"    = "#0072B2",
    "Plausible-best"  = "#27AE60"
  )

  obs_lines <- tibble::tibble(
    label = "Observed (all pandemics)",
    value = obs_final,
    linetype = "All pandemics"
  )

  # Dark-neutral grey for observed reference lines — reads as "data truth"
  # without competing with scenario accent colours.
  obs_line_colors <- c("All pandemics" = "grey25")

  p <- ggplot(plot_df, aes(y = label)) +
    geom_vline(data = obs_lines,
               aes(xintercept = value, colour = linetype),
               linewidth = 0.7, alpha = 0.8) +
    geom_linerange(aes(xmin = lower, xmax = upper, colour = label),
                   linewidth = 1.5, alpha = 0.35) +
    geom_point(aes(x = median, colour = label),
               size = 4, shape = 18) +
    scale_colour_manual(
      values = c(scenario_colors, obs_line_colors),
      name = NULL,
      breaks = c(names(obs_line_colors), names(scenario_colors))
    ) +
    scale_x_log10(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = c(10e6, 20e6, 50e6, 100e6, 200e6)
    ) +
    labs(x = "Cumulative pandemic deaths (1951\u20132025)", y = NULL) +
    theme_nature() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(family = "Arial", size = 8),
      axis.title.x = element_text(family = "Arial", size = 9, margin = margin(t = 10)),
      axis.text.y = element_text(family = "Arial", size = 9, face = "bold"),
      axis.text.x = element_text(family = "Arial", size = 8),
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    )

  return(p)
}
