#' Run Complete Pandemic Analysis on Dataset
#'
#' Performs a complete extreme value analysis on pandemic data using a two-stage
#' mixture model approach. This function orchestrates the entire analytical pipeline,
#' including data transformation, threshold selection, model fitting, bootstrap
#' simulation, burden calculation, and DALY estimation.
#'
#' @param data Data frame containing pandemic data with at least the columns:
#'   - deaths: Number of deaths in the pandemic event
#'   - start_year: Year when the pandemic started
#'   - pop_affected: Population affected by the pandemic (optional)
#' @param variable String, name of variable to analyze (typically "deaths")
#' @param conf Configuration object from config package containing analysis parameters:
#'   - pop_reference: Reference population for standardization
#'   - pop_min: Minimum population threshold for inclusion
#'   - base_year: Default starting year for analysis
#'   - thresholds: Various cutoff parameters for death counts
#'   - bootstrap: Configuration for uncertainty quantification
#'   - bulk_distribution: Distribution type for bulk modeling (lognormal, weibull)
#'   - use_constrained_model: Enable constrained survival mixture model
#'   - plots: Plot generation settings and options
#'   - calc_burden: Enable burden calculation
#'   - calc_dalys: Enable DALY estimation
#'   - calc_scenarios: Enable best/worst case scenario projections
#'   - time: Parameters for time trend modeling
#' @param threshold Optional manual threshold override for peaks-over-threshold analysis
#' @param year Optional year filter to restrict analysis to events after this year
#' @param lower_cutoff Optional override for the lower cutoff threshold
#' @param return_severity_draws Logical, if TRUE, the return object includes the bootstrap severity draws (default: FALSE)
#' @param output_dir Output directory for plots and saved draws
#' @return A list containing multiple analysis components:
#'   - specs: Data frame with model specifications and parameters
#'   - params: List of fitted model parameters (distribution types, shape, scale, etc.)
#'   - results: Expected values and point estimates with bootstrap uncertainty intervals
#'   - burden: Population-scaled death burden estimates (yearly and cumulative)
#'   - dalys: Global DALY estimates with uncertainty bounds
#'   - scenarios: Results from best/worst case scenario analysis (if calculated)
#'   - severity_draws: Bootstrap severity draws (if return_severity_draws=TRUE)
run_analysis <- function(
    data,
    variable,
    conf = config::get(),
    threshold = NULL,
    year = NULL,
    lower_cutoff = NULL,
    return_severity_draws = FALSE,
    output_dir = NULL) {
  #--------------------------------------------------
  # 1. DATA PREPARATION AND FILTERING
  #--------------------------------------------------
  # Set master random seed for reproducibility
  if (!is.null(conf$seed)) {
    set.seed(conf$seed)
    cat("Setting master random seed:", conf$seed, "\n")
  }

  # Keep copy of original data
  original_data <- data

  base_year <- ifelse(is.null(year), conf$base_year, year)

  data <- data |>
    filter(
      start_year >= base_year,
      deaths > conf$pop_min
    )

  # Generate descriptive plots if enabled
  if (conf$plots$descriptive) {
    tryCatch(print(create_descriptive_plots(data, variable, output_dir = output_dir, conf = conf)),
             error = function(e) message("Note: Could not display plot: ", e$message))
  }

  # Generate time-scaled deaths plot if enabled
  if (conf$plots$time_scaled_deaths) {
    min_scaled <- conf$plots$time_scaled_deaths_min %||% 5e4
    tryCatch(print(plot_time_scaled_deaths(data, conf$thresholds$lower_cutoff,
      reference_year = conf$time$current_year,
      include_labels = conf$plots$time_scaled_deaths_labels,
      min_deaths_scaled = min_scaled,
      output_dir = output_dir, conf = conf
    )), error = function(e) message("Note: Could not display plot: ", e$message))
    # Also generate dual-panel version (raw + scaled)
    plot_time_scaled_deaths(data, conf$thresholds$lower_cutoff,
      reference_year = conf$time$current_year,
      include_labels = FALSE,
      show_raw_panel = TRUE,
      min_deaths_scaled = min_scaled,
      output_dir = output_dir, conf = conf
    )
  }

  # Generate descriptive plots and tables if enabled
  if (isTRUE(conf$plots$descriptive)) {
    plot_severity_by_era(data, variable, conf, output_dir = output_dir)
    plot_severity_over_time(data,
      breakpoint_year = conf$time$breakpoint_year %||% 1950,
      output_dir = output_dir, conf = conf
    )
    generate_tier_composition_table(
      data, variable,
      tier_boundary = conf$thresholds[[variable]] %||% conf$thresholds$deaths_scaled,
      lower_cutoff = conf$thresholds$lower_cutoff,
      breakpoint_year = conf$time$breakpoint_year %||% 1950,
      output_dir = output_dir, conf = conf
    )
  }

  #--------------------------------------------------
  # 2. DATA TRANSFORMATION AND THRESHOLD SELECTION
  #--------------------------------------------------
  # Apply dual transformation to standardize pandemic size relative to population
  data$dual <- dual_transform(data, variable, h = conf$pop_reference, l = conf$pop_min)

  u_best <- determine_threshold(data, variable, threshold, conf)

  # Determine lower cutoff
  variable_specific_cutoff_name <- paste0("lower_cutoff_", variable)

  if (!is.null(lower_cutoff)) {
    # User-specified cutoff takes precedence
    conf$thresholds$lower_cutoff <- lower_cutoff
  } else if (!is.null(conf$thresholds[[variable_specific_cutoff_name]])) {
    # Variable-specific cutoff from config
    conf$thresholds$lower_cutoff <- conf$thresholds[[variable_specific_cutoff_name]]
  }
  # Else use the default lower_cutoff already in conf

  #--------------------------------------------------
  # 3. INITIAL SUMMARY STATISTICS
  #--------------------------------------------------
  # Create initial summary of data characteristics and thresholds
  specs <- tibble(
    Variable = variable,
    Year = base_year,
    Threshold = inv_dual_transform(u_best, h = conf$pop_reference, l = conf$pop_min),
    Cutoff = conf$thresholds$lower_cutoff,
    n = nrow(data),
    PoT = sum(data$dual > u_best),
    PcoT = mean(data$dual > u_best),
    PoB = sum(data$dual > conf$thresholds$lower_cutoff)
  )

  # Pre/post-1950 era counts by tier, for events above the analytical floor
  # (lower cutoff); events at or above the tail threshold are catastrophic.
  era_data <- data %>% dplyr::filter(.data[[variable]] >= conf$thresholds$lower_cutoff)
  catastrophic_cutoff <- inv_dual_transform(u_best, h = conf$pop_reference, l = conf$pop_min)
  specs$moderate_pre_1950      <- sum(era_data$start_year <  1950 & era_data[[variable]] <  catastrophic_cutoff)
  specs$moderate_post_1950     <- sum(era_data$start_year >= 1950 & era_data[[variable]] <  catastrophic_cutoff)
  specs$catastrophic_pre_1950  <- sum(era_data$start_year <  1950 & era_data[[variable]] >= catastrophic_cutoff)
  specs$catastrophic_post_1950 <- sum(era_data$start_year >= 1950 & era_data[[variable]] >= catastrophic_cutoff)

  #--------------------------------------------------
  # 4. INITIALIZE ANALYSIS OBJECTS
  #--------------------------------------------------
  model_results <- NULL # Will store fitted model components
  boot_results <- NULL # Will store bootstrap results
  burden <- NULL # Will store burden calculation results
  params <- list() # Will store model parameters
  results <- list() # Will store key result metrics

  # Ensure scenario-related objects always exist (validation may disable burden/scenarios)
  scenario_results <- NULL
  worst_case <- NULL
  best_case <- NULL

  #--------------------------------------------------
  # 5. MIXTURE MODEL FITTING
  #--------------------------------------------------
  # Fit two-stage mixture model with flexible bulk and GPD tail
  model_results <- fit_mixture_model(
    data, variable, conf, u_best, conf$bulk_distribution,
    constrain_survival = conf$use_constrained_model,
    summary = specs
  )

  if (!is.null(model_results)) {
    #--------------------------------------------------
    # 6. EXTRACT MODEL PARAMETERS
    #--------------------------------------------------
    params <- list(
      model_type = ifelse(is.null(model_results$bulk), "gpd_only", "mixture"),
      xi = model_results$tail$xi,
      xi_se = model_results$tail$xi_se,
      beta = model_results$tail$beta,
      beta_se = model_results$tail$beta_se
    )

    # Add confidence intervals using confint() if model fit is available
    if (!is.null(model_results$tail$fit)) {
      tryCatch(
        {
          ci <- confint(model_results$tail$fit)
          if (!is.null(ci) && nrow(ci) >= 2) {
            # Extract confidence intervals for both parameters
            params$xi_ci_lower <- ci[1, 1] # Shape parameter lower bound
            params$xi_ci_upper <- ci[1, 2] # Shape parameter upper bound
            params$beta_ci_lower <- ci[2, 1] # Scale parameter lower bound
            params$beta_ci_upper <- ci[2, 2] # Scale parameter upper bound
          }
        },
        error = function(e) {
          # If confint fails, calculate using normal approximation with constraints
          params$xi_ci_lower <<- model_results$tail$xi - 1.96 * model_results$tail$xi_se
          params$xi_ci_upper <<- model_results$tail$xi + 1.96 * model_results$tail$xi_se
          # Ensure beta confidence intervals are positive
          params$beta_ci_lower <<- max(0, model_results$tail$beta - 1.96 * model_results$tail$beta_se)
          params$beta_ci_upper <<- model_results$tail$beta + 1.96 * model_results$tail$beta_se
        }
      )
    }

    # Add bulk parameters if it's a mixture model
    if (!is.null(model_results$bulk)) {
      params$bulk_distribution <- model_results$bulk$dist_type
      params <- c(params, model_results$bulk$parameters)
      params$constrained_survival <- conf$use_constrained_model
    }

    #--------------------------------------------------
    # 7. UNCERTAINTY QUANTIFICATION
    #--------------------------------------------------
    # Initialize results list with shadow mean (expected value)
    results$shadow <- model_results$mean
    
    # Add point estimates for individual component means
    if (!is.null(model_results$bulk)) {
      results$bulk_mean <- model_results$bulk$conditional_mean
    } else {
      results$bulk_mean <- 0  # No bulk component in GPD-only model
    }
    results$tail_mean <- model_results$tail$shadow

    # Run bootstrap simulation if enabled for uncertainty estimation
    if (conf$bootstrap$enabled) {
      boot_results <- bootstrap_mixture_mean(
        data = data,
        variable = variable,
        conf = conf,
        mixture_pe = model_results$mean,
        u_best = u_best,
        summary = specs,
        output_dir = output_dir
      )

      results$boot <- boot_results$estimate
      results$boot_low <- boot_results$ci_lower
      results$boot_upp <- boot_results$ci_upper
      
      # Add bootstrap results for individual component means
      results$bulk_mean_boot <- boot_results$bulk_mean_estimate
      results$bulk_mean_boot_low <- boot_results$bulk_mean_ci_lower
      results$bulk_mean_boot_upp <- boot_results$bulk_mean_ci_upper
      results$tail_mean_boot <- boot_results$tail_mean_estimate
      results$tail_mean_boot_low <- boot_results$tail_mean_ci_lower
      results$tail_mean_boot_upp <- boot_results$tail_mean_ci_upper
      
      # Process parameter draws and create summaries (store separately to avoid tibble conversion issues)
      bootstrap_param_draws <- NULL
      bootstrap_param_summary <- NULL
      if (!is.null(boot_results$param_draws)) {
        bootstrap_param_summary <- calculate_parameter_summaries(boot_results$param_draws)
        bootstrap_param_draws <- boot_results$param_draws
        
        # Update params with bootstrap-derived uncertainty intervals
        params <- merge_bootstrap_params(params, bootstrap_param_summary)
      }
    }

    # Generate the mixture model plot once bootstrap results are available
    if (conf$plots$tail_mixture) {
      truncated_bulk_tail_plot(
        data = data,
        variable = variable,
        mixture_fit = model_results,
        pop_ref = conf$pop_reference,
        log10 = conf$plots$log_scale,
        y_title = conf$plots$tail_mixture_y_title,
        param_draws = bootstrap_param_draws,
        show_uncertainty = conf$bootstrap$enabled && !is.null(bootstrap_param_draws),
        point_size = if (!is.null(conf$plots$point_size)) conf$plots$point_size else 3,
        output_dir = output_dir, conf = conf
      )
    }

    #--------------------------------------------------
    # 8. BURDEN CALCULATION
    #--------------------------------------------------
    # Initialize burden variables for death forecasts
    burden <- NULL
    burden_results <- tibble()

    # Calculate population-scaled burden if enabled in config
    if (isTRUE(conf$calc_burden)) {
      burden <- calc_tiered_burden(
        data = original_data,
        variable = variable,
        conf = conf,
        bulk_severity_draws = boot_results$param_draws$bulk_mean,
        tail_severity_draws = boot_results$param_draws$tail_mean,
        shadow_mean = results$boot,
        shadow_low = results$boot_low,
        shadow_up = results$boot_upp,
        output_dir = output_dir,
        validation = FALSE,
        window_sizes = conf$window$default_size,
        return_severity_draws = return_severity_draws
      )

      if (conf$plots$forecast && !is.null(burden$plot)) {
        tryCatch(print(burden$plot), error = function(e) message("Note: Could not display plot: ", e$message))
      }

      if (conf$plots$time_trend && !is.null(burden$rate_plot)) {
        tryCatch(print(burden$rate_plot), error = function(e) message("Note: Could not display plot: ", e$message))
      }

      results$yearly_rate <- burden$result$yearly_rate

      burden_results <- tibble(
        yearly_rate = burden$result$yearly_rate[[1]]$mean,
        yearly_rate_detailed = burden$result$yearly_rate,
        yearly_deaths = burden$result$yearly_deaths,
        yearly_deaths_low = burden$result$yearly_deaths_low,
        yearly_deaths_up = burden$result$yearly_deaths_up,
        cum_deaths = burden$result$cum_deaths,
        cum_deaths_low = burden$result$cum_deaths_low,
        cum_deaths_up = burden$result$cum_deaths_up
      )
    }

    #--------------------------------------------------
    # 9. SCENARIO ANALYSIS
    #--------------------------------------------------
    # Plausible-best and plausible-worst scenarios from the frontier bounds on the tier rates
    if (!is.null(burden$scenarios)) {
      scenario_results <- burden$scenarios

      # Compute yearly deaths summaries from scenario draws
      build_scenario_result <- function(sc, draws_mat) {
        if (is.null(draws_mat)) {
          return(list(result = tibble::tibble(
            yearly_rate = NA, yearly_deaths = NA,
            yearly_deaths_low = NA, yearly_deaths_up = NA,
            cum_deaths = sc$cum_deaths, cum_deaths_low = sc$cum_deaths_low,
            cum_deaths_up = sc$cum_deaths_up),
            yearly_deaths_adjusted_draws = NULL))
        }
        yearly_means <- colMeans(draws_mat)
        list(result = tibble::tibble(
          yearly_rate = NA,
          yearly_deaths = stats::median(yearly_means),
          yearly_deaths_low = as.numeric(stats::quantile(yearly_means, 0.025)),
          yearly_deaths_up = as.numeric(stats::quantile(yearly_means, 0.975)),
          cum_deaths = sc$cum_deaths,
          cum_deaths_low = sc$cum_deaths_low,
          cum_deaths_up = sc$cum_deaths_up),
          yearly_deaths_adjusted_draws = draws_mat)
      }

      worst_case <- build_scenario_result(
        scenario_results$worst_case, burden$worst_case_draws)
      best_case <- build_scenario_result(
        scenario_results$best_case, burden$best_case_draws)

      # Update scenario_results with enriched data for table generation
      scenario_results$worst_case <- worst_case$result
      scenario_results$best_case <- best_case$result
    }
  }

  #--------------------------------------------------
  # 10. DALY CALCULATIONS
  #--------------------------------------------------
  # Initialize DALY result containers
  daly_results <- NULL

  if (isTRUE(conf$calc_burden) && isTRUE(conf$calc_dalys) && isTRUE(conf$calc_scenarios) && !is.null(burden_results)) {
    # Calculate DALY to death ratio
    daly_ratio_draws <- calc_daly_ratio_draws()
    cumulative_draws <- as.vector(colSums(burden$yearly_deaths_adjusted_draws))
    cumulative_worst_draws <- as.vector(colSums(worst_case$yearly_deaths_adjusted_draws))
    cumulative_best_draws <- as.vector(colSums(best_case$yearly_deaths_adjusted_draws))

    # Resample daly-to-death ratio draws (with reproducible seed)
    if (!is.null(conf$seed)) {
      set.seed(conf$seed + 2) # Offset to avoid interfering with other random operations
    }
    daly_ratio_draws <- sample(daly_ratio_draws, conf$draws, replace = TRUE)

    # Multiply DALYs while propagating uncertainty
    cum_dalys_base_draws <- cumulative_draws * daly_ratio_draws
    cum_dalys_worst_draws <- cumulative_worst_draws * daly_ratio_draws
    cum_dalys_best_draws <- cumulative_best_draws * daly_ratio_draws

    # Summarise the DALY draws: median and 95% uncertainty interval
    daly_results <- tibble(
      cum_dalys = median(cum_dalys_base_draws),
      cum_dalys_low = quantile(cum_dalys_base_draws, prob = 0.025),
      cum_dalys_up = quantile(cum_dalys_base_draws, prob = 0.975),
      daly_ratio = median(daly_ratio_draws),
      daly_ratio_low = quantile(daly_ratio_draws, prob = 0.025),
      daly_ratio_up = quantile(daly_ratio_draws, prob = 0.975),
      cum_dalys_best = median(cum_dalys_best_draws),
      cum_dalys_best_low = quantile(cum_dalys_best_draws, prob = 0.025),
      cum_dalys_best_up = quantile(cum_dalys_best_draws, prob = 0.975),
      cum_dalys_worst = median(cum_dalys_worst_draws),
      cum_dalys_worst_low = quantile(cum_dalys_worst_draws, prob = 0.025),
      cum_dalys_worst_up = quantile(cum_dalys_worst_draws, prob = 0.975)
    )

    # Calculate DALY gap between scenarios at draw level
    if (!is.null(scenario_results) && !is.null(cum_dalys_worst_draws) && !is.null(cum_dalys_best_draws)) {
      # Calculate gap for each draw (worst - best)
      dalys_gap_draws <- cum_dalys_worst_draws - cum_dalys_best_draws

      # Calculate summary statistics for absolute gap
      dalys_gap_median <- median(dalys_gap_draws)
      dalys_gap_lower <- quantile(dalys_gap_draws, 0.025)
      dalys_gap_upper <- quantile(dalys_gap_draws, 0.975)

      # Calculate relative gap (as percentage of worst case)
      dalys_relative_gap_draws <- dalys_gap_draws / cum_dalys_worst_draws * 100

      # Add DALY gap calculations to scenario results
      scenario_results$dalys_gap <- list(
        absolute = list(
          median = dalys_gap_median,
          lower_95 = dalys_gap_lower,
          upper_95 = dalys_gap_upper
        ),
        relative = list(
          median = median(dalys_relative_gap_draws),
          lower_95 = quantile(dalys_relative_gap_draws, 0.025),
          upper_95 = quantile(dalys_relative_gap_draws, 0.975)
        )
      )
    }
  }

  #--------------------------------------------------
  # 11. RESULT COMPILATION AND RETURN
  #--------------------------------------------------
  # Create and return the structured output
  result_list <- list(
    specs = specs,
    params = as_tibble(params),
    results = as_tibble(results),
    burden = burden_results,
    dalys = daly_results,
    scenarios = scenario_results
  )

  # Add parameter draws and summaries if available
  if (!is.null(bootstrap_param_draws)) {
    result_list$param_draws <- bootstrap_param_draws
  }
  if (!is.null(bootstrap_param_summary)) {
    result_list$param_summary <- bootstrap_param_summary
  }

  # Add bootstrap severity draws if requested
  if (return_severity_draws && !is.null(boot_results) && !is.null(boot_results$severity_draws)) {
    result_list$severity_draws <- boot_results$severity_draws
  }

  return(result_list)
}

#' Calculate Parameter Summaries from Bootstrap Draws
#'
#' @param param_draws List containing parameter draws from bootstrap
#' @return Tibble with parameter summaries including median and 95% UI
calculate_parameter_summaries <- function(param_draws) {
  # Initialize results
  summaries <- tibble()
  
  # Tail parameters
  tail_params <- c("xi", "beta", "sigma", "tail_prob")
  for (param in tail_params) {
    if (param %in% names(param_draws) && length(param_draws[[param]]) > 0) {
      values <- param_draws[[param]][!is.na(param_draws[[param]])]
      if (length(values) > 0) {
        summaries <- bind_rows(summaries, tibble(
          parameter = param,
          component = "tail",
          median = median(values, na.rm = TRUE),
          lower_95 = quantile(values, 0.025, na.rm = TRUE),
          upper_95 = quantile(values, 0.975, na.rm = TRUE),
          n_valid = length(values)
        ))
      }
    }
  }
  
  # Bulk probability 
  if ("bulk_prob" %in% names(param_draws) && length(param_draws$bulk_prob) > 0) {
    values <- param_draws$bulk_prob[!is.na(param_draws$bulk_prob)]
    if (length(values) > 0) {
      summaries <- bind_rows(summaries, tibble(
        parameter = "bulk_prob",
        component = "bulk",
        median = median(values, na.rm = TRUE),
        lower_95 = quantile(values, 0.025, na.rm = TRUE),
        upper_95 = quantile(values, 0.975, na.rm = TRUE),
        n_valid = length(values)
      ))
    }
  }
  
  # Bulk distribution parameters
  if ("bulk_params" %in% names(param_draws) && length(param_draws$bulk_params) > 0) {
    for (bulk_param in names(param_draws$bulk_params)) {
      values <- param_draws$bulk_params[[bulk_param]]
      values <- values[!is.na(values) & values != 0]  # Remove NA and zero values
      if (length(values) > 0) {
        summaries <- bind_rows(summaries, tibble(
          parameter = bulk_param,
          component = "bulk",
          median = median(values, na.rm = TRUE),
          lower_95 = quantile(values, 0.025, na.rm = TRUE),
          upper_95 = quantile(values, 0.975, na.rm = TRUE),
          n_valid = length(values)
        ))
      }
    }
  }
  
  return(summaries)
}

#' Merge Bootstrap Parameter Summaries with Point Estimates
#'
#' @param params Original parameter list from point estimation
#' @param param_summary Parameter summaries from bootstrap draws
#' @return Updated parameter list with bootstrap-derived intervals
merge_bootstrap_params <- function(params, param_summary) {
  # Update tail parameters with bootstrap intervals
  tail_param_map <- list(
    "xi" = c("xi_boot_median", "xi_boot_lower", "xi_boot_upper"),
    "beta" = c("beta_boot_median", "beta_boot_lower", "beta_boot_upper"),
    "sigma" = c("sigma_boot_median", "sigma_boot_lower", "sigma_boot_upper"),
    "tail_prob" = c("tail_prob_boot_median", "tail_prob_boot_lower", "tail_prob_boot_upper")
  )
  
  for (param in names(tail_param_map)) {
    param_row <- param_summary[param_summary$parameter == param & param_summary$component == "tail", ]
    if (nrow(param_row) > 0) {
      param_names <- tail_param_map[[param]]
      params[[param_names[1]]] <- param_row$median
      params[[param_names[2]]] <- param_row$lower_95
      params[[param_names[3]]] <- param_row$upper_95
    }
  }
  
  # Add bulk parameter summaries
  bulk_params <- param_summary[param_summary$component == "bulk", ]
  for (i in seq_len(nrow(bulk_params))) {
    param_name <- bulk_params$parameter[i]
    base_name <- paste0(param_name, "_boot")
    params[[paste0(base_name, "_median")]] <- bulk_params$median[i]
    params[[paste0(base_name, "_lower")]] <- bulk_params$lower_95[i]
    params[[paste0(base_name, "_upper")]] <- bulk_params$upper_95[i]
  }
  
  return(params)
}
