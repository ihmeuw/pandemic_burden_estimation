#' Estimate SFA Bounds from GLM Residuals (COLS-style)
#'
#' Derives upper and lower frontier shifts from the GLM's fitted rates using
#' method-of-moments on composed-error residuals (Aigner, Lovell & Schmidt 1977;
#' COLS tradition). The bounds are structurally tied to the GLM mean, so
#' lower <= mean <= upper by construction.
#'
#' On the log-rate scale, residuals e_i = log(observed + eps) - log(fitted + eps)
#' decompose into symmetric noise v ~ N(0, sigma_v^2) and one-sided inefficiency
#' u ~ HalfNormal(0, sigma_u^2). The frontier shift E[u] = sigma_u * sqrt(2/pi)
#' quantifies how far the envelope sits from the mean.
#'
#' @param observed_rates Numeric vector of observed windowed exceedance rates
#' @param fitted_rates Numeric vector of GLM-fitted rates (same length)
#' @param window_size Window size for epsilon floor (default 15)
#' @return List with sigma_u, sigma_v, E_u (frontier shift on log scale),
#'         and multipliers (exp(E_u) for upper, exp(-E_u) for lower)
estimate_sfa_bounds <- function(observed_rates, fitted_rates, window_size = 15) {
  valid <- !is.na(observed_rates) & !is.na(fitted_rates)
  obs <- observed_rates[valid]
  fit <- fitted_rates[valid]

  epsilon <- 0.5 / window_size
  e <- log(obs + epsilon) - log(fit + epsilon)

  m2 <- var(e)
  m3 <- mean((e - mean(e))^3)

  # Method-of-moments: sigma_u from third central moment of residuals
  # For the half-normal composed error, m3 = -sigma_u^3 * sqrt(2/pi) * (1 - 4/pi)
  # so sigma_u^3 = |m3| / (sqrt(2/pi) * (1 - 4/pi))
  c3 <- sqrt(2 / pi) * (1 - 4 / pi)
  sigma_u_cubed <- abs(m3) / abs(c3)
  sigma_u <- sigma_u_cubed^(1/3)

  # sigma_v from remaining variance: m2 = sigma_v^2 + sigma_u^2 * (1 - 2/pi)
  sigma_v_sq <- m2 - sigma_u^2 * (1 - 2 / pi)
  sigma_v <- if (sigma_v_sq > 0) sqrt(sigma_v_sq) else 0

  E_u <- sigma_u * sqrt(2 / pi)

  list(
    sigma_u = sigma_u,
    sigma_v = sigma_v,
    E_u = E_u,
    upper_multiplier = exp(E_u),
    lower_multiplier = exp(-E_u),
    n_obs = sum(valid),
    residual_skewness = m3 / m2^(3/2)
  )
}

#' Predict SFA Bounds at Given Years
#'
#' Applies the COLS-derived frontier shift to the GLM fitted rates.
#' upper = fitted * exp(+E[u]), lower = fitted * exp(-E[u]).
#'
#' @param fitted_rates Numeric vector of GLM-fitted rates at prediction years
#' @param bounds Output from estimate_sfa_bounds()
#' @param direction "upper" or "lower"
#' @return Numeric vector of frontier rates
sfa_bound_predict <- function(fitted_rates, bounds, direction = "upper") {
  mult <- if (direction == "upper") bounds$upper_multiplier else bounds$lower_multiplier
  fitted_rates * mult
}

#' Plot SFA Frontiers for Multiple Tiers (Faceted)
#'
#' @param tiers List of tier specifications, each with: rate_data, bounds, label
#' @param forecast_years Vector of forecast years to extend frontiers
#' @param current_year Year dividing historical/forecast
#' @param conf Configuration object
#' @param output_dir Output directory for saving
#' @return ggplot object (faceted by tier)
plot_sfa_frontiers <- function(tiers, forecast_years = NULL, current_year = 2025,
                                conf = NULL, output_dir = NULL) {
  all_hist <- list()
  all_frontier <- list()

  for (tier in tiers) {
    rd <- tier$rate_data
    hist_data <- rd %>% dplyr::filter(type == "historical")
    all_years <- sort(unique(c(rd$year, forecast_years)))

    # Central trend and bounds both use the smooth model-implied rate directly
    # when available. That column (fitted_rate_smooth) is rowMeans(exp(eta_draws)
    # / window_size), i.e. the actual parametric shape of the selected model for
    # this tier (piecewise, nb_linear, etc.).
    central_col <- if ("fitted_rate_smooth" %in% names(rd)) "fitted_rate_smooth" else "fitted_rate"
    fitted_trend <- approx(rd$year, rd[[central_col]],
                            xout = all_years, rule = 2)$y

    upper_rates <- sfa_bound_predict(fitted_trend, tier$bounds, "upper")
    lower_rates <- sfa_bound_predict(fitted_trend, tier$bounds, "lower")

    all_hist[[tier$label]] <- hist_data %>%
      dplyr::mutate(tier = tier$label)

    all_frontier[[tier$label]] <- tibble::tibble(
      year = rep(all_years, 3),
      rate = c(upper_rates, lower_rates, fitted_trend),
      line = rep(c("Plausible-worst",
                    "Plausible-best",
                    "Fitted trend"), each = length(all_years)),
      tier = tier$label
    )
  }

  hist_df <- dplyr::bind_rows(all_hist)
  frontier_df <- dplyr::bind_rows(all_frontier)

  tier_levels <- sapply(tiers, `[[`, "label")
  hist_df$tier <- factor(hist_df$tier, levels = tier_levels)
  frontier_df$tier <- factor(frontier_df$tier, levels = tier_levels)

  # Preserve a meaningful legend order: worst on top, central in middle, best at bottom.
  line_levels <- c("Plausible-worst", "Fitted trend", "Plausible-best")
  frontier_df$line <- factor(frontier_df$line, levels = line_levels)

  # Toned red/green/blue palette: muted red and emerald green keep the
  # intuitive "worst / best" semantics.
  line_colors <- c(
    "Plausible-worst" = "#C0392B",  # muted red (firebrick-ish)
    "Fitted trend"    = "#0072B2",  # blue (matches plot_tiered_rates)
    "Plausible-best"  = "#27AE60"   # emerald green, toned
  )
  # All lines solid; visual hierarchy carried by weight + alpha rather than
  # linetype (dashed at <0.5 pt fragments into dots in print).
  line_types <- c(
    "Plausible-worst" = "solid",
    "Fitted trend"    = "solid",
    "Plausible-best"  = "solid"
  )

  year_breaks <- seq(
    floor(min(frontier_df$year) / 50) * 50,
    ceiling(max(frontier_df$year) / 50) * 50,
    by = 50
  )

  # Panel letters a, b, ... for caption references.
  panel_tags <- data.frame(
    tier  = factor(tier_levels, levels = tier_levels),
    label = letters[seq_along(tier_levels)]
  )

  # Separate fitted-trend layer so we can render bounds at lower weight + alpha.
  trend_df <- frontier_df %>% dplyr::filter(line == "Fitted trend")
  bounds_df <- frontier_df %>% dplyr::filter(line != "Fitted trend")

  p <- ggplot() +
    # Forecast region shading (both panels) — replaces orphan 2025 vline.
    annotate("rect",
             xmin = current_year, xmax = Inf,
             ymin = -Inf, ymax = Inf,
             fill = "grey92", alpha = 0.5) +
    geom_point(data = hist_df, aes(x = year, y = observed_rate),
               colour = "grey32", alpha = 0.35, size = 0.8) +
    geom_line(data = bounds_df, aes(x = year, y = rate, colour = line, linetype = line),
              linewidth = 0.6, alpha = 0.7) +
    geom_line(data = trend_df,  aes(x = year, y = rate, colour = line, linetype = line),
              linewidth = 0.9) +
    geom_text(
      data = panel_tags,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.5, vjust = 1.3, size = 8 / .pt, fontface = "bold", family = "Arial"
    ) +
    facet_wrap(~ tier, ncol = 1, scales = "free_y",
               axes = "all_x", axis.labels = "margins") +
    scale_colour_manual(values = line_colors, breaks = line_levels) +
    scale_linetype_manual(values = line_types, breaks = line_levels) +
    scale_x_continuous(breaks = year_breaks) +
    labs(x = "Year", y = "Exceedance rate (events per year)",
         colour = NULL, linetype = NULL) +
    theme_nature()

  save_plot_if_enabled(p, "sfa_frontiers", output_dir, conf, scale = 0.8)
  p
}
