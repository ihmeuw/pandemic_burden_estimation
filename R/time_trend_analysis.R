#' Time trend analysis functions
#'
#' Contains functions for analyzing pandemic exceedance rate trends and comparing temporal models.

#' Analyze pandemic exceedance rates with time trend
#'
#' Complete time series analysis framework for pandemic exceedance rates with
#' multiple model types and uncertainty quantification through Monte Carlo simulation.
#'
#' @param data Event-level dataset (requires start_year plus the analysis variable)
#' @param lower_cutoff Threshold value to define exceedance
#' @param end_year Ending year for the analysis (default: current year)
#' @param future_years Number of years to forecast (default: 75)
#' @param window_sizes Window size in years for the sliding-window exceedance counts (default from config)
#' @param n_sims Number of simulations for prediction intervals (default: 1000)
#' @param model_type Type of temporal model to use ("linear", "quadratic",
#' "constant", "nb_linear", "nb_quadratic", "piecewise", "nb_piecewise", "piecewise_full", "nb_piecewise_full", or "auto"). The default config uses "piecewise" for Tier 1 and "nb_linear" for Tier 2; the validation profile overrides both tiers to "linear".
#' @param conf Configuration object containing time.current_year for forecasting reference
#' @param generate_draws Whether to generate prediction draws for uncertainty propagation (default: FALSE)
#' @param n_draws Number of draws to generate for prediction intervals (default: 1000)
#'
#' @return List with exceedance_rate_df, optimal_window, model_comparison, best_model, forecast_summary, waiting_time, window_counts, and draws
analyze_pandemic_trend <- function(
    data,
    variable = "deaths_scaled",
    lower_cutoff,
    upper_threshold = NULL,
    start_year = NULL,
    end_year = NULL,
    future_years = if (!is.null(conf) && !is.null(conf$time$future_years)) conf$time$future_years else 75,
    window_sizes = conf$window$default_size,
    n_sims = if (!is.null(conf) && !is.null(conf$draws)) conf$draws else 1000,
    model_type = "linear",
    conf = NULL,
    generate_draws = FALSE,
    n_draws = if (!is.null(conf) && !is.null(conf$draws)) conf$draws else 1000) {
  if (is.null(start_year)) {
    start_year <- min(data$start_year)
  }
  if (is.null(end_year)) {
    end_year <- max(data$start_year, as.numeric(format(Sys.Date(), "%Y")))
  }
  filtered_data <- data %>%
    dplyr::filter(start_year >= start_year) %>%
    dplyr::arrange(start_year)

  exceedance_results <- analyze_exceedance_rates(
    data = filtered_data,
    variable = variable,
    threshold = lower_cutoff,
    base_year = start_year,
    window_sizes = window_sizes,
    upper_threshold = upper_threshold,
    conf = conf
  )
  model_comparison <- compare_temporal_models(
    data = filtered_data,
    variable = variable,
    threshold = lower_cutoff,
    window_size = exceedance_results$optimal_window,
    base_year = start_year,
    upper_threshold = upper_threshold,
    conf = conf
  )
  if (model_type != "auto") {
    valid_models <- names(model_comparison$models)
    model_map <- c(
      "constant" = "Constant Rate",
      "linear" = "Linear Trend",
      "quadratic" = "Quadratic Trend",
      "nb_linear" = "NB Linear",
      "nb_quadratic" = "NB Quadratic",
      "piecewise" = "Piecewise",
      "nb_piecewise" = "NB Piecewise",
      "piecewise_full" = "Piecewise Full",
      "nb_piecewise_full" = "NB Piecewise Full"
    )
    requested_model <- model_map[model_type]
    if (!requested_model %in% valid_models) {
      warning(paste0(
        "Requested model '", model_type, "' (", requested_model, ") not available. ",
        "Using best model: ", model_comparison$comparison$Model[1]
      ))
    } else {
      user_model <- model_comparison$models[[requested_model]]
      model_comparison$best_model <- user_model
      custom_order <- c(requested_model, setdiff(valid_models, requested_model))
      model_idx <- match(custom_order, model_comparison$comparison$Model)
      model_comparison$comparison <- model_comparison$comparison[model_idx, ]
      cat("Using user-specified model:", requested_model, "\n")
    }
  } else {
    cat("Using auto-selected best model:", model_comparison$comparison$Model[1], "\n")
  }
  forecast_results <- forecast_exceedance_rates(
    model_results = model_comparison,
    forecast_years = future_years,
    window_size = exceedance_results$optimal_window,
    prediction_interval = TRUE,
    n_sims = n_sims
  )
  historical_detailed <- exceedance_results$optimal_results$window_counts %>%
    dplyr::select(
      year = window_end,
      window_start,
      observed_rate = rate,
      exceedances,
      fitted_rate = annual_rate,
      predicted,
      ci_lower_rate,
      ci_upper_rate
    ) %>%
    dplyr::mutate(
      type = "historical",
      year_idx = year - min(year) + 1
    )
  historical_rates <- historical_detailed %>%
    dplyr::select(year, observed_rate, fitted_rate, type, year_idx)
  forecast_rates <- forecast_results$forecast_data %>%
    dplyr::select(year = window_end, fitted_rate = annual_rate, lower_95, upper_95) %>%
    dplyr::mutate(
      type = "forecast",
      year_idx = max(historical_rates$year_idx) + dplyr::row_number(),
      observed_rate = NA_real_
    )
  exceedance_rate_df <- dplyr::bind_rows(
    historical_detailed %>%
      dplyr::select(
        year, type, year_idx,
        observed_rate,
        fitted_rate,
        exceedances,
        lower_95 = ci_lower_rate,
        upper_95 = ci_upper_rate
      ),
    forecast_rates
  ) %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(
      lower_95 = dplyr::if_else(is.na(lower_95), fitted_rate * 0.5, lower_95),
      upper_95 = dplyr::if_else(is.na(upper_95), fitted_rate * 1.5, upper_95)
    )
  draws_results <- NULL
  if (generate_draws) {
    cat("Generating prediction draws for uncertainty quantification...\n")
    best_model <- model_comparison$best_model
    window_counts <- model_comparison$window_counts
    window_size <- exceedance_results$optimal_window
    last_observed_year <- max(window_counts$window_end)
    current_year <- ifelse(is.null(conf$time$current_year),
      as.numeric(format(Sys.Date(), "%Y")),
      conf$time$current_year
    )
    forecast_end <- current_year + future_years
    years_seq <- min(window_counts$window_end):forecast_end
    years_df <- data.frame(window_end = years_seq)
    years_df$window_start <- years_df$window_end - window_size + 1
    years_df$year_centered <- years_df$window_end - min(window_counts$window_end)
    years_df$year_scaled <- years_df$year_centered / max(window_counts$year_centered)
    years_df$year_squared <- years_df$year_scaled^2
    years_df$offset <- window_size
    # Add breakpoint indicator for piecewise models
    breakpoint_year <- if (!is.null(conf$time$breakpoint_year)) conf$time$breakpoint_year else 1950
    years_df$post_breakpoint <- as.numeric(years_df$window_end >= breakpoint_year)
    years_df$year_scaled_x_post <- years_df$year_scaled * years_df$post_breakpoint
    set.seed(conf$seed + 3)
    param_draws <- MASS::mvrnorm(
      n = n_draws,
      mu = coef(best_model),
      Sigma = vcov(best_model)
    )
    X <- model.matrix(formula(best_model)[-2], years_df)
    eta_draws <- X %*% t(param_draws)
    if (inherits(best_model, "negbin")) {
      theta <- best_model$theta
      annual_rates <- matrix(0, nrow = length(years_seq), ncol = n_draws)
      for (i in seq_len(n_draws)) {
        lambda <- exp(eta_draws[, i])
        if (!is.null(conf$seed)) {
          set.seed(conf$seed + 4 + i)
        }
        counts <- rnbinom(length(lambda), size = theta, mu = lambda)
        annual_rates[, i] <- counts / window_size
      }
    } else {
      lambda_draws <- exp(eta_draws)
      annual_rates <- matrix(0, nrow = length(years_seq), ncol = n_draws)
      for (i in seq_len(n_draws)) {
        counts <- rpois(length(years_seq), lambda = lambda_draws[, i])
        annual_rates[, i] <- counts / window_size
      }
    }
    draws_summary <- data.frame(
      year = years_seq,
      mean = rowMeans(annual_rates),
      median = apply(annual_rates, 1, median),
      lower_95 = apply(annual_rates, 1, quantile, 0.025),
      upper_95 = apply(annual_rates, 1, quantile, 0.975),
      lower_80 = apply(annual_rates, 1, quantile, 0.1),
      upper_80 = apply(annual_rates, 1, quantile, 0.9)
    )
    draws_summary$type <- ifelse(
      draws_summary$year <= current_year,
      "historical",
      "forecast"
    )

    # Confidence interval for the mean rate (parameter uncertainty only,
    # no Poisson/NB observation noise). This produces a narrower band that
    # reflects model uncertainty and is continuous across the forecast boundary.
    mean_rate_draws <- exp(eta_draws) / window_size
    draws_summary$ci_lower_95 <- apply(mean_rate_draws, 1, quantile, 0.025)
    draws_summary$ci_upper_95 <- apply(mean_rate_draws, 1, quantile, 0.975)
    # Smooth model-implied central rate (plotting only — avoids Monte Carlo
    # wiggle from per-year count sampling in `annual_rates`)
    draws_summary$fitted_rate_smooth <- rowMeans(mean_rate_draws)

    draws_results <- list(
      annual_rates = annual_rates,
      years = years_seq,
      param_draws = param_draws,
      summary = draws_summary,
      observed_years = exceedance_rate_df %>%
        dplyr::filter(type == "historical") %>%
        dplyr::mutate(observed = TRUE) %>%
        dplyr::select(year, observed_rate, observed)
    )

    # Overwrite fitted_rate with the selected model's mean predictions so that
    # historical and forecast periods use the same model (e.g., piecewise) rather
    # than the initial linear model from analyze_exceedance_rates for history only.
    match_idx <- match(exceedance_rate_df$year, draws_summary$year)
    valid <- !is.na(match_idx)
    exceedance_rate_df$fitted_rate[valid] <- draws_summary$mean[match_idx[valid]]
    exceedance_rate_df$lower_95[valid] <- draws_summary$ci_lower_95[match_idx[valid]]
    exceedance_rate_df$upper_95[valid] <- draws_summary$ci_upper_95[match_idx[valid]]
    # Carry the smooth rate column separately so plots can use it without
    # affecting downstream consumers that read `fitted_rate`.
    exceedance_rate_df$fitted_rate_smooth <- NA_real_
    exceedance_rate_df$fitted_rate_smooth[valid] <- draws_summary$fitted_rate_smooth[match_idx[valid]]
  }
  return(list(
    exceedance_rate_df = exceedance_rate_df,
    optimal_window = exceedance_results$optimal_window,
    model_comparison = model_comparison$comparison,
    best_model = model_comparison$best_model,
    forecast_summary = forecast_results$summary,
    waiting_time = forecast_results$waiting_summary,
    window_counts = model_comparison$window_counts,
    draws = draws_results
  ))
}

#' Analyze exceedance rates over time
#'
#' @return List containing analysis results and diagnostics
analyze_exceedance_rates <- function(data,
                                     variable = "deaths_scaled",
                                     threshold,
                                     base_year = 1700,
                                     window_sizes = c(10, 20, 30, 50, 75, 100),
                                     upper_threshold = NULL,
                                     conf = NULL) {
  filtered_data <- data %>%
    dplyr::filter(start_year >= base_year) %>%
    dplyr::arrange(start_year)
  min_year <- min(filtered_data$start_year)
  max_year <- max(filtered_data$start_year)
  current_year <- ifelse(is.null(conf$time$current_year),
    as.numeric(format(Sys.Date(), "%Y")),
    conf$time$current_year
  )
  effective_max_year <- max(max_year, current_year)
  results <- list()
  results$data <- filtered_data
  results$parameters <- list(
    variable = variable,
    threshold = threshold,
    base_year = base_year,
    window_sizes = window_sizes
  )
  window_results <- list()
  for (window_size in window_sizes) {
    window_ends <- seq(min_year + window_size - 1, effective_max_year, by = 1)
    window_starts <- window_ends - window_size + 1
    window_counts <- tibble::tibble(
      window_start = window_starts,
      window_end = window_ends,
      exceedances = 0,
      events = 0,
      rate = 0
    )
    for (i in seq_len(nrow(window_counts))) {
      window_start <- window_counts$window_start[i]
      window_end <- window_counts$window_end[i]
      window_events <- filtered_data %>%
        dplyr::filter(start_year >= window_start, start_year <= window_end)
      total_events <- nrow(window_events)
      if (!is.null(upper_threshold)) {
        exceedances <- sum(window_events[[variable]] > threshold &
                           window_events[[variable]] <= upper_threshold)
      } else {
        exceedances <- sum(window_events[[variable]] > threshold)
      }
      window_counts$exceedances[i] <- exceedances
      window_counts$events[i] <- total_events
      window_counts$rate[i] <- exceedances / window_size
    }
    window_stats <- tibble::tibble(
      window_size = window_size,
      total_exceedances = sum(window_counts$exceedances),
      total_events = sum(window_counts$events),
      mean_exceedances = mean(window_counts$exceedances),
      sd_exceedances = sd(window_counts$exceedances),
      var_exceedances = var(window_counts$exceedances),
      dispersion = var(window_counts$exceedances) / mean(window_counts$exceedances),
      zero_proportion = mean(window_counts$exceedances == 0),
      mean_rate = mean(window_counts$rate)
    )
    window_counts$year_center <- window_counts$window_end - min(window_counts$window_end)
    model <- glm(exceedances ~ year_center, family = poisson(link = "log"), data = window_counts)
    window_counts$predicted <- predict(model, type = "response")
    window_counts$annual_rate <- window_counts$predicted / window_size
    pred_ci <- predict(model,
      newdata = window_counts,
      type = "link", se.fit = TRUE
    )
    ci_data <- tibble::tibble(
      fit = pred_ci$fit,
      se = pred_ci$se.fit
    ) %>%
      dplyr::mutate(
        ci_lower = exp(fit - 1.96 * se),
        ci_upper = exp(fit + 1.96 * se),
        predicted_rate = exp(fit) / window_size,
        ci_lower_rate = ci_lower / window_size,
        ci_upper_rate = ci_upper / window_size
      )
    window_counts <- window_counts %>%
      dplyr::bind_cols(
        ci_data %>%
          dplyr::select(ci_lower, ci_upper, predicted_rate, ci_lower_rate, ci_upper_rate)
      )
    window_results[[as.character(window_size)]] <- list(
      window_counts = window_counts,
      window_stats = window_stats,
      model = model
    )
  }
  window_comparisons <- dplyr::bind_rows(
    lapply(window_results, function(x) x$window_stats)
  )
  if (length(window_sizes) == 1) {
    optimal_window <- window_sizes
  } else {
    dispersion_diff <- abs(window_comparisons$dispersion - 1)
    optimal_idx <- which.min(dispersion_diff)
    optimal_window <- window_comparisons$window_size[optimal_idx]
  }
  results$optimal_window <- optimal_window
  results$optimal_results <- window_results[[as.character(optimal_window)]]
  results$window_comparison <- window_comparisons
  optimal_data <- results$optimal_results$window_counts
  main_plot <- ggplot(optimal_data, aes(x = window_end)) +
    geom_line(aes(y = predicted_rate), color = "blue", linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lower_rate, ymax = ci_upper_rate), alpha = 0.2, fill = "blue") +
    geom_point(aes(y = rate), alpha = 0.6) +
    labs(
      title = paste0("Pandemic Exceedance Rate Over Time (Window: ", optimal_window, " years)"),
      subtitle = paste0(
        "Events exceeding ", format(threshold, scientific = TRUE, digits = 2),
        " in variable '", variable, "' (rates based on previous ", optimal_window, " years)"
      ),
      x = "Year (end of backwards-looking window)",
      y = "Annual Exceedance Rate"
    ) +
    theme_nature() +
    theme(plot.title = element_text(face = "bold"))
  results$main_plot <- main_plot
  return(results)
}

#' Fit and compare multiple temporal models
#'
#' @return List of model comparison results
compare_temporal_models <- function(data,
                                    variable = "deaths_scaled",
                                    threshold,
                                    window_size = 50,
                                    base_year = 1700,
                                    upper_threshold = NULL,
                                    conf = NULL) {
  filtered_data <- data %>%
    dplyr::filter(start_year >= base_year) %>%
    dplyr::arrange(start_year)
  min_year <- min(filtered_data$start_year)
  max_year <- max(filtered_data$start_year)
  current_year <- ifelse(is.null(conf$time$current_year),
    as.numeric(format(Sys.Date(), "%Y")),
    conf$time$current_year
  )
  effective_max_year <- max(max_year, current_year)
  window_ends <- seq(min_year + window_size - 1, effective_max_year, by = 1)
  window_starts <- window_ends - window_size + 1
  window_counts <- tibble::tibble(
    window_start = window_starts,
    window_end = window_ends,
    exceedances = 0,
    events = 0,
    rate = 0
  )
  for (i in seq_len(nrow(window_counts))) {
    window_start <- window_counts$window_start[i]
    window_end <- window_counts$window_end[i]
    window_events <- filtered_data %>%
      dplyr::filter(start_year >= window_start, start_year <= window_end)
    total_events <- nrow(window_events)
    if (!is.null(upper_threshold)) {
      exceedances <- sum(window_events[[variable]] > threshold &
                         window_events[[variable]] <= upper_threshold)
    } else {
      exceedances <- sum(window_events[[variable]] > threshold)
    }
    window_counts$exceedances[i] <- exceedances
    window_counts$events[i] <- total_events
    window_counts$rate[i] <- exceedances / window_size
  }
  window_counts <- window_counts %>%
    dplyr::mutate(
      year_centered = window_end - min(window_end),
      year_scaled = year_centered / max(year_centered),
      year_squared = year_scaled^2
    )
  null_model <- glm(exceedances ~ 1,
    family = poisson(link = "log"),
    data = window_counts
  )
  linear_model <- glm(exceedances ~ year_scaled,
    family = poisson(link = "log"),
    data = window_counts
  )
  quadratic_model <- glm(exceedances ~ year_scaled + year_squared,
    family = poisson(link = "log"),
    data = window_counts
  )
  nb_model <- try(
    MASS::glm.nb(exceedances ~ year_scaled, data = window_counts),
    silent = TRUE
  )
  nb_quad_model <- try(
    MASS::glm.nb(exceedances ~ year_scaled + year_squared, data = window_counts),
    silent = TRUE
  )

  # Piecewise models with intercept shift at breakpoint year
  breakpoint_year <- if (!is.null(conf$time$breakpoint_year)) conf$time$breakpoint_year else 1950
  window_counts$post_breakpoint <- as.numeric(window_counts$window_end >= breakpoint_year)
  window_counts$year_scaled_x_post <- window_counts$year_scaled * window_counts$post_breakpoint

  # Piecewise with intercept shift only
  piecewise_model <- try(
    glm(exceedances ~ year_scaled + post_breakpoint,
        family = poisson(link = "log"),
        data = window_counts),
    silent = TRUE
  )
  nb_piecewise_model <- try(
    MASS::glm.nb(exceedances ~ year_scaled + post_breakpoint, data = window_counts),
    silent = TRUE
  )

  # Full piecewise with both intercept and slope change
  piecewise_full_model <- try(
    glm(exceedances ~ year_scaled + post_breakpoint + year_scaled_x_post,
        family = poisson(link = "log"),
        data = window_counts),
    silent = TRUE
  )
  nb_piecewise_full_model <- try(
    MASS::glm.nb(exceedances ~ year_scaled + post_breakpoint + year_scaled_x_post, data = window_counts),
    silent = TRUE
  )

  models <- list(
    "Constant Rate" = null_model,
    "Linear Trend" = linear_model,
    "Quadratic Trend" = quadratic_model
  )
  if (!inherits(nb_model, "try-error")) {
    models[["NB Linear"]] <- nb_model
  }
  if (!inherits(nb_quad_model, "try-error")) {
    models[["NB Quadratic"]] <- nb_quad_model
  }
  if (!inherits(piecewise_model, "try-error")) {
    models[["Piecewise"]] <- piecewise_model
  }
  if (!inherits(nb_piecewise_model, "try-error")) {
    models[["NB Piecewise"]] <- nb_piecewise_model
  }
  if (!inherits(piecewise_full_model, "try-error")) {
    models[["Piecewise Full"]] <- piecewise_full_model
  }
  if (!inherits(nb_piecewise_full_model, "try-error")) {
    models[["NB Piecewise Full"]] <- nb_piecewise_full_model
  }
  model_comparison <- tibble::tibble(
    Model = names(models),
    AIC = sapply(models, AIC),
    BIC = sapply(models, BIC),
    LogLik = sapply(models, function(m) as.numeric(logLik(m)))
  ) %>%
    dplyr::arrange(AIC)
  best_model_name <- model_comparison$Model[1]
  best_model <- models[[best_model_name]]
  prediction_data <- window_counts %>%
    dplyr::select(window_start, window_end, year_scaled, year_squared, post_breakpoint, year_scaled_x_post)
  for (model_name in names(models)) {
    model <- models[[model_name]]
    pred <- predict(model, newdata = prediction_data, type = "response")
    pred_rate <- pred / window_size
    prediction_data[[paste0(model_name, "_pred")]] <- pred_rate
  }
  prediction_data$actual_rate <- window_counts$rate
  pred_cols <- grep("_pred$", names(prediction_data), value = TRUE)
  models_data <- prediction_data %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(pred_cols),
      names_to = "model_name",
      values_to = "predicted_rate"
    ) %>%
    dplyr::mutate(
      model_name = gsub("_pred$", "", model_name),
      is_best = model_name == best_model_name
    )
  model_plot <- ggplot() +
    geom_point(
      data = prediction_data, aes(x = window_end, y = actual_rate),
      alpha = 0.4, color = "black"
    ) +
    geom_line(
      data = models_data %>%
        dplyr::left_join(prediction_data %>% dplyr::select(window_start, window_end), by = "window_start"),
      aes(
        x = window_end, y = predicted_rate, color = model_name,
        size = is_best, linetype = is_best
      )
    ) +
    scale_size_manual(values = c("TRUE" = 1.2, "FALSE" = 0.8), guide = "none") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_color_manual(
      values = c(
        "Constant Rate" = "red",
        "Linear Trend" = "blue",
        "Quadratic Trend" = "green",
        "NB Linear" = "purple",
        "NB Quadratic" = "orange",
        "Piecewise" = "#00468B",
        "NB Piecewise" = "#0099B4",
        "Piecewise Full" = "#42B540",
        "NB Piecewise Full" = "#925E9F"
      ),
      name = "Model Type"
    ) +
    labs(
      title = "Model Comparison for Pandemic Exceedance Rate",
      subtitle = paste0("Window Size: ", window_size, " years (looking back from each year)"),
      x = "Year (end of backwards-looking window)",
      y = "Annual Exceedance Rate"
    ) +
    theme_nature() +
    theme(legend.position = "bottom")
  return(list(
    models = models,
    best_model = best_model,
    comparison = model_comparison,
    prediction_data = prediction_data,
    model_plot = model_plot,
    window_counts = window_counts,
    conf = conf
  ))
}
