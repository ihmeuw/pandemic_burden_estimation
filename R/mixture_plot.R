#' Plot the Two-Component Severity Model with Bootstrap Parameter Uncertainty
#'
#' Draws the fitted composite severity model against the empirical exceedance
#' probabilities: a truncated lognormal for moderate events (lower cutoff to tail
#' threshold) and a generalised Pareto distribution (GPD) for events above the
#' threshold, each with a bootstrap uncertainty band.
#'
#' @details
#' The bulk curve is the fitted lognormal survival function; continuity with the
#' tail at the threshold is enforced inside fit_bulk_component(), so no rescaling
#' is applied here. The uncertainty bands are quantiles, across bootstrap draws,
#' of the model survival curves:
#' - Bulk: S_bulk,i(x) = 1 - F_lognormal(x; mu_i, sigma_i)
#' - Tail: S_tail,i(x) = S_bulk,i(threshold) * (1 - F_GPD(x; xi_i, beta_i, threshold))
#' so each draw is continuous at the threshold. The bands show parameter
#' uncertainty given the observed data.
#'
#' @param data Data frame containing the empirical data points
#' @param variable String, name of the variable to plot
#' @param mixture_fit Fitted mixture model containing bulk and tail components
#' @param pop_ref Numeric, population reference for the dual transformation
#' @param y_title String, y-axis label
#' @param log10 Logical, whether to use a logarithmic y-axis. Default is FALSE
#' @param param_draws List of bootstrap parameter draws for the uncertainty bands:
#'   - xi, beta: vectors of GPD parameter draws
#'   - bulk_params: list with the lognormal parameter draws (meanlog, sdlog)
#' @param show_uncertainty Logical, whether to display uncertainty bands. Default is TRUE
#' @param ci_level Numeric, level of the uncertainty bands (0-1). Default is 0.95
#' @param ci_alpha Numeric, transparency of the uncertainty bands (0-1). Default is 0.28
#' @param point_size Numeric, size of empirical data points. Default is 3
#' @param min_display Numeric, smallest event size drawn as a point
#' @param x_axis_start Numeric, lower limit of the x-axis
#' @param font_scale Numeric, multiplier applied to all text sizes
#' @param output_dir Output directory
#' @param conf Configuration object
#'
#' @return A ggplot object
#'
#' @references
#' Cirillo, P., & Taleb, N. N. (2020). Tail risk of contagious diseases. Nature Physics, 16(6), 606-613.
#'
#' @export
truncated_bulk_tail_plot <- function(data, variable, mixture_fit, pop_ref,
                                     y_title = "Historical probability of exceeding 2025-equivalent deaths",
                                     log10 = FALSE,
                                     param_draws = NULL, show_uncertainty = TRUE,
                                     ci_level = 0.95, ci_alpha = 0.28, point_size = 3,
                                     min_display = 5e4,
                                     x_axis_start = 5e4,
                                     font_scale = 1.3,
                                     output_dir = NULL, conf = NULL) {
  # Extract parameters from mixture fit
  lower_cutoff <- mixture_fit$lower_cutoff # Lower threshold
  u_best <- mixture_fit$tail_threshold # Upper threshold (GPD threshold)

  # Bulk (lognormal) parameters
  params <- mixture_fit$bulk$parameters

  # Maximum plot range
  plot_max <- max(data[[variable]]) * 1.5

  # Empirical data: sort, compute empirical survival
  all_data <- data %>%
    arrange(.data[[variable]]) %>%
    mutate(
      x = .data[[variable]],
      y = 1 - ecdf(.data$x)(.data$x)
    )

  # Avoid zero probability at minimum
  all_data$y[which.min(all_data$y)] <- 0.006

  # Extract tail parameters
  xi <- mixture_fit$tail$xi
  beta <- mixture_fit$tail$beta
  fit <- mixture_fit$tail$fit
  shadow <- mixture_fit$tail$shadow

  # Bulk component: fitted lognormal survival function up to the threshold
  bulk_x <- seq(lower_cutoff, u_best, length.out = 500)
  bulk_surv <- 1 - plnorm(bulk_x, params$meanlog, params$sdlog)

  # Tail component: GPD survival function from the threshold
  tail_x <- seq(u_best, plot_max, length.out = 500)
  prob_less_thresh <- fit$p.less.thresh
  tail_surv <- (1 - prob_less_thresh) * (1 - evir::pgpd(tail_x, xi, u_best, beta))

  bulk_data <- data.frame(x = bulk_x, y = bulk_surv)
  tail_data <- data.frame(x = tail_x, y = tail_surv)

  # =============================================================================
  # BOOTSTRAP PARAMETER UNCERTAINTY
  # =============================================================================
  bulk_uncertainty <- NULL
  tail_uncertainty <- NULL

  if (!is.null(param_draws) && show_uncertainty && length(param_draws$xi) > 1) {
    alpha <- 1 - ci_level
    lower_q <- alpha / 2
    upper_q <- 1 - alpha / 2

    # Bulk band: survival curve of each bootstrap lognormal draw
    if (!is.null(param_draws$bulk_params) && length(param_draws$bulk_params) > 0) {
      n_draws <- length(param_draws$xi)

      # Rows = x-values, columns = bootstrap draws
      bulk_surv_draws <- matrix(NA, nrow = length(bulk_x), ncol = n_draws)

      for (i in 1:n_draws) {
        meanlog_i <- param_draws$bulk_params$meanlog[i]
        sdlog_i <- param_draws$bulk_params$sdlog[i]

        if (!is.na(meanlog_i) && !is.na(sdlog_i) && sdlog_i > 0) {
          bulk_surv_draws[, i] <- 1 - plnorm(bulk_x, meanlog_i, sdlog_i)
        }
      }

      if (sum(!is.na(bulk_surv_draws[1, ])) > 2) { # Need at least 3 valid draws
        bulk_lower <- apply(bulk_surv_draws, 1, quantile, probs = lower_q, na.rm = TRUE)
        bulk_upper <- apply(bulk_surv_draws, 1, quantile, probs = upper_q, na.rm = TRUE)

        bulk_uncertainty <- data.frame(
          x = bulk_x,
          y_lower = bulk_lower,
          y_upper = bulk_upper
        )
      }
    }

    # Tail band: each GPD draw starts where its own bulk curve ends at the
    # threshold, so every draw is continuous there
    if (length(param_draws$xi) > 2 && length(param_draws$beta) > 2) {
      n_draws <- length(param_draws$xi)
      tail_surv_draws <- matrix(NA, nrow = length(tail_x), ncol = n_draws)

      # Step 1: S_bulk,i(threshold) for each draw
      bulk_at_threshold_draws <- numeric(n_draws)

      for (i in 1:n_draws) {
        meanlog_i <- param_draws$bulk_params$meanlog[i]
        sdlog_i <- param_draws$bulk_params$sdlog[i]

        if (!is.na(meanlog_i) && !is.na(sdlog_i) && sdlog_i > 0) {
          bulk_at_threshold_draws[i] <- 1 - plnorm(u_best, meanlog_i, sdlog_i)
        }
      }

      # Step 2: S_tail,i(x) = S_bulk,i(threshold) * S_GPD(x; xi_i, beta_i, threshold)
      for (i in 1:n_draws) {
        xi_i <- param_draws$xi[i]
        beta_i <- param_draws$beta[i]

        if (!is.na(xi_i) && !is.na(beta_i) && beta_i > 0 && !is.na(bulk_at_threshold_draws[i])) {
          gpd_surv_i <- 1 - evir::pgpd(tail_x, xi_i, u_best, beta_i)
          tail_surv_draws[, i] <- bulk_at_threshold_draws[i] * gpd_surv_i
        }
      }

      # Step 3: quantiles across bootstrap draws
      if (sum(!is.na(tail_surv_draws[1, ])) > 2) { # Need at least 3 valid draws
        tail_lower <- apply(tail_surv_draws, 1, quantile, probs = lower_q, na.rm = TRUE)
        tail_upper <- apply(tail_surv_draws, 1, quantile, probs = upper_q, na.rm = TRUE)

        tail_uncertainty <- data.frame(
          x = tail_x,
          y_lower = tail_lower,
          y_upper = tail_upper
        )
      }
    }
  }

  # Labels for the two thresholds
  lower_cutoff_label <- paste0(
    scales::comma(round(lower_cutoff)),
    " deaths"
  )

  threshold_label <- paste0(
    scales::comma(round(inv_dual_transform(u_best, h = pop_ref, l = 1e4))),
    " deaths"
  )

  # Colour palette
  plot_palette <- c(
    "Primary" = "#00468B", # Primary blue
    "Highlight" = "#ED0000", # Signal red
    "Secondary" = "#42B540", # Secondary green
    "Tertiary" = "#0099B4", # Tertiary blue
    "Quaternary" = "#925E9F" # Quaternary purple
  )

  bulk_model_label <- "Log-normal model"
  gpd_model_label <- "GPD model"

  # Construct the plot (all data points shown)
  p <- all_data %>%
    ggplot() +

    # Uncertainty bands for bulk component
    {
      if (!is.null(bulk_uncertainty) && show_uncertainty) {
        geom_ribbon(
          data = bulk_uncertainty,
          aes(x = .data$x, ymin = .data$y_lower, ymax = .data$y_upper),
          fill = plot_palette[["Quaternary"]], # Purple for log-normal
          alpha = ci_alpha
        )
      }
    } +

    # Uncertainty bands for tail component
    {
      if (!is.null(tail_uncertainty) && show_uncertainty) {
        geom_ribbon(
          data = tail_uncertainty,
          aes(x = .data$x, ymin = .data$y_lower, ymax = .data$y_upper),
          fill = plot_palette[["Secondary"]], # Green for GPD
          alpha = ci_alpha
        )
      }
    } +

    # Empirical points (events below min_display are not drawn; curves/ecdf unchanged)
    geom_point(
      data = all_data %>% dplyr::filter(.data$x >= min_display),
      aes(x = .data$x, y = y),
      size = point_size, alpha = 0.4
    ) +

    # Bulk distribution curve (truncated at threshold)
    geom_line(
      data = bulk_data,
      aes(x = .data$x, y = y, color = bulk_model_label),
      linewidth = 1
    ) +

    # GPD tail curve
    geom_line(
      data = tail_data,
      aes(x = .data$x, y = y, color = "GPD model"),
      linewidth = 1
    ) +

    # Vertical reference lines for the two thresholds
    geom_vline(
      aes(xintercept = lower_cutoff, color = "Lower Cutoff"),
      linetype = "dashed",
      linewidth = 0.25
    ) +
    geom_vline(
      aes(xintercept = u_best, color = "Tail Threshold"),
      linetype = "dashed",
      linewidth = 0.25
    ) +

    # Log10 scale on x-axis
    scale_x_log10(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = c(1e5, 1e6, 1e7, 1e8),
      limits = c(x_axis_start, NA),
      oob = scales::oob_keep
    ) +

    labs(
      x = "2025-equivalent deaths (log scale)",
      y = y_title,
      color = NULL,
      fill = NULL
    ) +

    # Threshold text annotations
    annotate(
      "text",
      x = lower_cutoff * 1.2, y = 1.1,
      label = lower_cutoff_label,
      color = plot_palette["Quaternary"],
      hjust = 0,
      size = 2.8 * font_scale,
      family = "Arial"
    ) +
    annotate(
      "text",
      x = u_best * 1.2, y = 1.1,
      label = threshold_label,
      color = plot_palette["Secondary"],
      hjust = 0,
      size = 2.8 * font_scale,
      family = "Arial"
    ) +

    # Repel labels for extreme points. Text is 2.99 mm on the 12-in canvas, about
    # 5 pt at a 183 mm print width; the seed fixes label placement across formats.
    ggrepel::geom_text_repel(
      data = subset(all_data, all_data$x > shadow),
      aes(x = .data$x, y = y, label = .data$name),
      size = 2.3 * font_scale,
      seed = 42,
      family = "Arial",
      box.padding = 0.5,
      max.overlaps = 40,
      color = "#666666",
      segment.size = 0.2,
      segment.color = "#999999",
      min.segment.length = 0
    )

  color_values <- c(
    "Lower Cutoff" = plot_palette["Primary"],
    "Tail Threshold" = plot_palette["Highlight"],
    "GPD Tail" = plot_palette["Highlight"]
  )
  color_values[bulk_model_label] <- plot_palette["Quaternary"]
  color_values[gpd_model_label] <- plot_palette["Secondary"]

  p <- p +
    theme_nature() +
    theme(
      legend.key.size = unit(1.2, "lines"),
      axis.title   = element_text(size = 11 * font_scale),
      axis.text    = element_text(size = 8  * font_scale),
      legend.title = element_text(size = 9  * font_scale, face = "bold"),
      legend.text  = element_text(size = 8  * font_scale),
      plot.title   = element_text(size = 10 * font_scale, face = "bold")
    ) +
    scale_color_manual(values = color_values) +
    {
      fill_values <- c()
      fill_values[bulk_model_label] <- plot_palette["Quaternary"]
      fill_values["GPD model"] <- plot_palette["Secondary"]
      scale_fill_manual(values = fill_values)
    }

  # Optionally use log10 scale on y-axis
  if (log10) {
    p <- p + scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
  }

  save_plot_if_enabled(p, "mixture_model_plot", output_dir, conf, width = 12, height = 8)

  tryCatch(print(p), error = function(e) message("Note: Could not display plot: ", e$message))

  invisible(p)
}
