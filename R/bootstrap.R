#' Bootstrap Mixture Model Mean
#'
#' @param data Data frame containing the input data
#' @param variable Name of the variable to analyze
#' @param conf Configuration list containing bootstrap parameters
#' @param mixture_pe Original point estimate from mixture model for comparison
#'
#' @return List with the bootstrap point estimate (median), 95% uncertainty interval, bulk/tail component means and their intervals, the severity draws and the full parameter draws
#' @export
bootstrap_mixture_mean <- function(data, variable, conf, mixture_pe, u_best, summary, output_dir = NULL) {
  # Input validation
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("Invalid input data: must be a non-empty data frame")
  }
  if (is.null(conf$draws) || conf$draws < 1) {
    stop("Invalid draws configuration")
  }

  boot_means <- numeric(conf$draws)
  boot_models  <- list()

  # Set seed for reproducible bootstrap sampling
  if (!is.null(conf$seed)) {
    set.seed(conf$seed + 1)  # Offset to avoid interfering with other random operations
  }

  # Initialize containers for parameter draws
  param_draws <- list(
    # Tail parameters
    xi = numeric(conf$draws),
    beta = numeric(conf$draws),
    sigma = numeric(conf$draws),
    tail_prob = numeric(conf$draws),
    
    # Bulk parameters (distribution-dependent)
    bulk_params = list(),
    bulk_prob = numeric(conf$draws),
    bulk_dist_type = character(conf$draws),
    
    # Component means
    bulk_mean = numeric(conf$draws),
    tail_mean = numeric(conf$draws),
    
    # Model identifiers
    model_type = character(conf$draws),
    constrained = logical(conf$draws)
  )

  # Perform bootstrap sampling
  for (i in seq_len(conf$draws)) {
    # Resample data with replacement
    boot_sample <- data[sample(nrow(data), replace = TRUE), ]

    # Fit mixture model for bootstrap sample with error handling
    boot_result <- tryCatch(
      {
        boot_model <- fit_mixture_model(
          boot_sample,
          variable,
          conf,
          u_best,
          dist_type = conf$bulk_distribution,
          constrain_survival = conf$use_constrained_model,
          summary = summary
        )
        # Validate that mean exists and is non-NULL
        if (is.null(boot_model$mean) || length(boot_model$mean) == 0) {
          NULL
        } else {
          boot_model
        }
      },
      error = function(e) {
        NULL
      }
    )

    # Handle failed bootstrap iterations
    if (is.null(boot_result)) {
      boot_means[i] <- NA
      boot_models[[i]] <- NULL
      next
    }

    boot_means[i] <- boot_result$mean
    boot_models[[i]] <- boot_result

    # Extract and store all parameters from this draw
    # Tail parameters
    param_draws$xi[i] <- boot_result$tail$xi
    param_draws$beta[i] <- boot_result$tail$beta
    param_draws$sigma[i] <- boot_result$tail$sigma
    param_draws$tail_prob[i] <- boot_result$tail$tail_prob

    # Model characteristics
    param_draws$model_type[i] <- ifelse(is.null(boot_result$bulk), "gpd_only", "mixture")
    param_draws$constrained[i] <- boot_result$constrained

    # Bulk parameters (if mixture model)
    if (!is.null(boot_result$bulk)) {
      param_draws$bulk_prob[i] <- boot_result$bulk$bulk_prob
      param_draws$bulk_dist_type[i] <- boot_result$bulk$dist_type

      # Store bulk distribution parameters
      bulk_params_i <- boot_result$bulk$parameters
      for (param_name in names(bulk_params_i)) {
        if (!param_name %in% names(param_draws$bulk_params)) {
          param_draws$bulk_params[[param_name]] <- numeric(conf$draws)
        }
        param_draws$bulk_params[[param_name]][i] <- bulk_params_i[[param_name]]
      }
    } else {
      # GPD-only model
      param_draws$bulk_prob[i] <- 0
      param_draws$bulk_dist_type[i] <- "none"
    }

    # Extract individual component means
    param_draws$tail_mean[i] <- boot_result$tail$shadow
    if (!is.null(boot_result$bulk)) {
      param_draws$bulk_mean[i] <- boot_result$bulk$conditional_mean
    } else {
      param_draws$bulk_mean[i] <- 0  # No bulk component in GPD-only model
    }
  }

  # Filter out failed bootstrap iterations (NA values from fitting failures)
  fitting_valid_idx <- !is.na(boot_means)
  n_fitting_failed <- sum(!fitting_valid_idx)
  if (n_fitting_failed > 0) {
    message(sprintf("Note: %d of %d bootstrap iterations failed due to insufficient exceedances",
                    n_fitting_failed, conf$draws))
  }

  boot_means <- boot_means[fitting_valid_idx]
  boot_models <- boot_models[fitting_valid_idx]

  # Filter parameter draws for fitting failures
  param_draws$xi <- param_draws$xi[fitting_valid_idx]
  param_draws$beta <- param_draws$beta[fitting_valid_idx]
  param_draws$sigma <- param_draws$sigma[fitting_valid_idx]
  param_draws$tail_prob <- param_draws$tail_prob[fitting_valid_idx]
  param_draws$bulk_prob <- param_draws$bulk_prob[fitting_valid_idx]
  param_draws$bulk_dist_type <- param_draws$bulk_dist_type[fitting_valid_idx]
  param_draws$bulk_mean <- param_draws$bulk_mean[fitting_valid_idx]
  param_draws$tail_mean <- param_draws$tail_mean[fitting_valid_idx]
  param_draws$model_type <- param_draws$model_type[fitting_valid_idx]
  param_draws$constrained <- param_draws$constrained[fitting_valid_idx]

  for (param_name in names(param_draws$bulk_params)) {
    param_draws$bulk_params[[param_name]] <- param_draws$bulk_params[[param_name]][fitting_valid_idx]
  }

  # Filter out invalid xi values (xi == -1)
  xi <- vapply(boot_models, function(m) m$tail$xi, numeric(1))
  xi_valid_idx <- xi != -1

  boot_means <- boot_means[xi_valid_idx]
  xi <- xi[xi_valid_idx]

  # Filter parameter draws to remove invalid models (xi == -1)
  param_draws$xi <- param_draws$xi[xi_valid_idx]
  param_draws$beta <- param_draws$beta[xi_valid_idx]
  param_draws$sigma <- param_draws$sigma[xi_valid_idx]
  param_draws$tail_prob <- param_draws$tail_prob[xi_valid_idx]
  param_draws$bulk_prob <- param_draws$bulk_prob[xi_valid_idx]
  param_draws$bulk_dist_type <- param_draws$bulk_dist_type[xi_valid_idx]
  param_draws$bulk_mean <- param_draws$bulk_mean[xi_valid_idx]
  param_draws$tail_mean <- param_draws$tail_mean[xi_valid_idx]
  param_draws$model_type <- param_draws$model_type[xi_valid_idx]
  param_draws$constrained <- param_draws$constrained[xi_valid_idx]

  # Filter bulk parameters
  for (param_name in names(param_draws$bulk_params)) {
    param_draws$bulk_params[[param_name]] <- param_draws$bulk_params[[param_name]][xi_valid_idx]
  }

  # Calculate statistics
  point_estimate <- stats::median(boot_means, na.rm = TRUE)
  ci <- stats::quantile(
    boot_means,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  # Calculate bulk and tail mean statistics
  bulk_mean_estimate <- stats::median(param_draws$bulk_mean, na.rm = TRUE)
  bulk_mean_ci <- stats::quantile(
    param_draws$bulk_mean,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  tail_mean_estimate <- stats::median(param_draws$tail_mean, na.rm = TRUE)
  tail_mean_ci <- stats::quantile(
    param_draws$tail_mean,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )

  # Generate diagnostic plot if requested
  if (conf$bootstrap$plot) {
    create_bootstrap_plot(
      boot_means = boot_means,
      point_estimate = point_estimate,
      mixture_pe = mixture_pe,
      ci = ci,
      variable = variable,
      output_dir = output_dir,
      conf = conf
    )
  }

  list(
    estimate = point_estimate,
    ci_lower = ci[1],
    ci_upper = ci[2],
    bulk_mean_estimate = bulk_mean_estimate,
    bulk_mean_ci_lower = bulk_mean_ci[1],
    bulk_mean_ci_upper = bulk_mean_ci[2],
    tail_mean_estimate = tail_mean_estimate,
    tail_mean_ci_lower = tail_mean_ci[1],
    tail_mean_ci_upper = tail_mean_ci[2],
    severity_draws = boot_means,
    param_draws = param_draws
  )
}

#' Create Diagnostic Plot for Bootstrap Results
#'
#' @param boot_means Vector of bootstrap means
#' @param point_estimate Calculated point estimate
#' @param mixture_pe Original point estimate
#' @param ci Numeric vector of the 2.5th and 97.5th bootstrap percentiles
#' @param variable Name of the variable being analyzed
#'
#' @noRd
create_bootstrap_plot <- function(boot_means, point_estimate, mixture_pe, ci, variable, output_dir = NULL, conf = NULL) {
  # Filter extreme values for better visualization
  plot_means <- boot_means[boot_means < (ci[2] * 1.1)]

  # Create data frame for ggplot
  plot_data <- data.frame(values = plot_means)
  
  # Create ggplot histogram
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = values)) +
    ggplot2::geom_histogram(bins = 50, fill = "lightgray", color = "black", alpha = 0.7) +
    ggplot2::geom_vline(aes(xintercept = point_estimate, color = "Bootstrapped Estimate"), size = 1) +
    ggplot2::geom_vline(aes(xintercept = mixture_pe, color = "Point Estimate"), size = 1) +
    ggplot2::geom_vline(aes(xintercept = ci[1], color = "95% CI"), linetype = "dashed", size = 1) +
    ggplot2::geom_vline(aes(xintercept = ci[2], color = "95% CI"), linetype = "dashed", size = 1) +
    ggplot2::scale_color_manual(
      name = "Reference Lines",
      values = c("Bootstrapped Estimate" = "red", "Point Estimate" = "purple", "95% CI" = "blue")
    ) +
    ggplot2::labs(
      title = paste0("Bootstrap Distribution of ", variable, " (Mixture Model)"),
      x = "Mixture Mean",
      y = "Frequency"
    ) +
    theme_nature() +
    ggplot2::theme(legend.position = "top")
  
  # Save plot if output directory is provided
  save_plot_if_enabled(p, "bootstrap_distribution", output_dir, conf)
  
  # Print plot for interactive use
  tryCatch(print(p), error = function(e) message("Note: Could not display plot: ", e$message))
  
  invisible(p)
}
