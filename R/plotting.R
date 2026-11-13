#' @importFrom magrittr %>%
#' @importFrom dplyr filter
utils::globalVariables(c(".data", "%>%", "dual_transform", "determine_threshold"))

#' Nature-style theme for publication-ready figures
#'
#' Shared theme function replacing per-plot inline theme blocks.
#' All text is Arial (no Times New Roman). Designed for Nature submission.
#'
#' @param base_size Base font size (default 10)
#' @return ggplot2 theme object
theme_nature <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "Arial"),
      plot.title       = element_text(family = "Arial", size = 10, face = "bold",
                                       hjust = 0, margin = margin(b = 10)),
      plot.subtitle    = element_text(family = "Arial", size = 9, hjust = 0,
                                       margin = margin(b = 15)),
      plot.caption     = element_text(family = "Arial", size = 8, color = "gray30",
                                       hjust = 0, margin = margin(t = 10)),
      axis.title       = element_text(family = "Arial", size = 11),
      axis.title.x     = element_text(margin = margin(t = 10)),
      axis.title.y     = element_text(angle = 90, margin = margin(r = 10)),
      axis.text        = element_text(family = "Arial", size = 8),
      axis.line        = element_line(color = "black", linewidth = 0.5),
      axis.ticks       = element_line(linewidth = 0.5),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "bottom",
      legend.title     = element_text(family = "Arial", size = 9, face = "bold"),
      legend.text      = element_text(family = "Arial", size = 8),
      legend.key.size  = unit(1.2, "lines"),
      legend.key       = element_rect(fill = "white", color = NA),
      strip.text       = element_text(family = "Arial", size = 9, face = "bold",
                                       hjust = 0, margin = margin(b = 5)),
      panel.spacing    = unit(1.2, "lines"),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

#' Plot Histogram of Variable
#'
#' This function creates a histogram of the specified variable from the data.
#'
#' @param data Data frame containing the data.
#' @param variable String, name of the variable to plot.
#' @param bin Numeric, number of bins in the histogram. Default is 40.
#' @param log Logical, whether to use logarithmic scale for x-axis. Default is FALSE.
#' @return A ggplot object representing the histogram.
plot_histogram <- function(data, variable, bin = 40, log = FALSE) {
  p <- ggplot2::ggplot(data, aes(x = .data[[variable]])) +
    geom_histogram(
      bins = bin,
      fill = "skyblue",
      color = "black"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = paste(variable, "(in thousands)"), y = "Number of Events") +
    coord_cartesian(clip = "off") +
    theme_nature() +
    theme(axis.text.x = element_text(margin = margin(t = 10)))

  if (log == TRUE) {
    p <- p +
      scale_x_log10(labels = trans_format("log10", math_format())) +
      annotation_logticks(sides = "b", outside = TRUE)
  }

  return(p)
}

#' Create Zipf Plot
#'
#' This function creates a Zipf plot (log-log survival plot) of the specified variable.
#'
#' @param data Data frame containing the data.
#' @param variable String, name of the variable to plot.
#' @return A ggplot object representing the Zipf plot.
zipf_plot <- function(data, variable) {
  x <- sort(as.numeric(data[[variable]]))
  ypoints <- 1 - ppoints(data[[variable]])

  temp <- tibble(x, ypoints)

  temp %>%
    ggplot(aes(x = x, y = ypoints)) +
    geom_point(
      size = 1,
      alpha = 0.8,
      col = "purple"
    ) +
    scale_x_log10(labels = trans_format("log10", math_format())) +
    scale_y_log10() +
    labs(x = paste(variable, "(in thousands)"), y = "Empirical survival (log scale)") +
    annotation_logticks(sides = "bl", outside = TRUE) +
    coord_cartesian(clip = "off") +
    theme_nature() +
    theme(
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 10))
    )
}


#' Create Mean Excess Plot
#'
#' This function creates a mean excess plot of the specified variable.
#'
#' @param data Data frame containing the data.
#' @param variable String, name of the variable to plot.
#' @return A ggplot object representing the mean excess plot.
mean_excess_plot <- function(data, variable) {
  x <- as.numeric(data[[variable]])
  myrank <- function(x, na_last = TRUE) {
    ranks <- sort.list(sort.list(x, na.last = na_last))
    if (is.na(na_last)) {
      x <- x[!is.na(x)]
    }
    for (i in unique(x[duplicated(x)])) {
      which <- x == i & !is.na(x)
      ranks[which] <- max(ranks[which])
    }
    ranks
  }
  x <- sort(x)
  n_excess <- unique(floor(length(x) - myrank(x)))
  points <- unique(x)
  nl <- length(points)
  n_excess <- n_excess[-nl]
  points <- points[-nl]
  excess <- cumsum(rev(x))[n_excess] - n_excess * points
  y <- excess / n_excess
  omit <- 3
  xx <- points[1:(nl - omit)]
  yy <- y[1:(nl - omit)]

  ggplot(dplyr::tibble(xx, yy), aes(x = xx, y = yy)) +
    geom_point(
      size = 1,
      alpha = 0.8,
      col = "purple"
    ) +
    labs(x = "Threshold", y = "Mean Excess") +
    coord_cartesian(clip = "off") +
    theme_nature() +
    theme(
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 10))
    )
}

#' Create Combined Descriptive Plots
#'
#' Creates a combined view of descriptive plots including histogram (normal and log scale),
#' Zipf plot, and mean excess plot.
#'
#' @param data Data frame containing the analysis data
#' @param variable String, name of variable to analyze
#' @param bins Integer, number of bins for histograms (default: 40)
#' @return A ggplot object containing the combined plots
#' @export
create_descriptive_plots <- function(data, variable, bins = 40, output_dir = NULL, conf = NULL) {
  # Generate individual plots
  p1 <- plot_histogram(data, variable, bins, log = FALSE) +
    labs(title = "Regular Histogram")

  p2 <- plot_histogram(data, variable, bins, log = TRUE) +
    labs(title = "Log-scale Histogram")

  p3 <- zipf_plot(data, variable) +
    labs(title = "Zipf Plot")

  p4 <- mean_excess_plot(data, variable) +
    labs(title = "Mean Excess Plot")

  # Combine plots using patchwork
  combined_plots <- (p1 + p2) / (p3 + p4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste("Descriptive Analysis:", variable),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
    )

  # Save plot if output directory is provided
  save_plot_if_enabled(combined_plots, "descriptive_plots", output_dir, conf, width = 12, height = 10)
  
  return(combined_plots)
}

#' Plot Time-Scaled Deaths
#'
#' Creates a professional visualization of pandemic and epidemic deaths over time, with deaths scaled
#' to the reference population year. Shows events through time with text labels
#' for the largest events. Uses enhanced visual styling for improved readability and presentation.
#'
#' @param data Data frame containing the pandemic/epidemic data
#' @param lower_cutoff Numeric threshold for highlighting events with text labels
#' @param label_threshold Numeric threshold for applying text labels to events (defaults to 1e7)
#' @param min_year Numeric, earliest year to include in the visualization (defaults to 1600)
#' @param reference_year Numeric, population reference year for scaling (defaults to 2025)
#' @param include_labels Logical, whether to include title, subtitle, and caption (defaults to TRUE)
#' @return A ggplot object representing the enhanced time-scaled deaths plot
#' @export
plot_time_scaled_deaths <- function(data, lower_cutoff = 2e5, label_threshold = 1e7,
                                    min_year = 1600, reference_year = 2025,
                                    include_labels = TRUE,
                                    show_raw_panel = FALSE,
                                    min_deaths_scaled = 5e4,
                                    output_dir = NULL, conf = NULL) {
  # Filter data to events from min_year onwards
  data_time_full <- data |>
    dplyr::select(name, start_year, end_year, deaths, deaths_scaled, type) |>
    dplyr::filter(start_year >= min_year)

  # Compute y-axis range from full data (before min_deaths_scaled filter)
  # so the axis reflects where smaller events would have been
  y_range_scaled <- range(data_time_full$deaths_scaled, na.rm = TRUE)
  y_range_raw <- range(data_time_full$deaths, na.rm = TRUE)

  # Apply visibility filter: drop events smaller than min_deaths_scaled
  data_time <- data_time_full |>
    dplyr::filter(deaths_scaled >= min_deaths_scaled)

  # Add period categorization for visual encoding
  data_time <- data_time |>
    dplyr::mutate(
      # Create period groupings
      period = case_when(
        start_year < 1800 ~ "Pre-1800",
        start_year < 1900 ~ "19th Century",
        start_year < 1950 ~ "1900-1950",
        start_year < 2000 ~ "1950-2000",
        TRUE ~ "21st Century"
      ),
      # Order periods chronologically
      period = factor(period, levels = c("Pre-1800", "19th Century", "1900-1950", "1950-2000", "21st Century")),
      # Calculate deaths in millions for labels
      deaths_millions = deaths_scaled / 1e6,
      # Create formatted label with deaths
      event_label = if_else(
        deaths_scaled > label_threshold,
        paste0(name, "\n(", round(deaths_millions, 1), "M)"),
        name
      )
    )

  # Define color palette based on periods
  period_colors <- c(
    "Pre-1800" = "#925E9F", # Quaternary (purple) for oldest events
    "19th Century" = "#0099B4", # Tertiary (teal)
    "1900-1950" = "#42B540", # Secondary (green)
    "1950-2000" = "#00468B", # Primary (blue)
    "21st Century" = "#ED0000" # Highlight (red) for most recent events
  )

  # Define alpha and size based on recency and magnitude
  point_size_range <- c(2, 8)
  max_deaths <- max(data_time$deaths_scaled, na.rm = TRUE)

  # Create the enhanced plot
  p <- ggplot(data_time, aes(x = start_year, y = deaths_scaled)) +
    # Add minimal reference lines for major time periods (hair-line thin per guidelines)
    geom_vline(
      xintercept = c(1800, 1900, 1950, 2000), linetype = "dotted",
      color = "gray50", linewidth = 0.25
    ) +

    # Add threshold reference lines (hair-line thin per guidelines)
    geom_hline(yintercept = lower_cutoff,
      linetype = "dashed", color = "grey50", linewidth = 0.25
    ) +
    annotate("text",
      x = min(data_time$start_year), y = lower_cutoff * 0.8,
      label = "Analysis threshold", color = "grey50", hjust = 0, fontface = "italic", size = 3
    ) +
    geom_hline(yintercept = 3e6,
      linetype = "dashed", color = "grey50", linewidth = 0.25
    ) +
    annotate("text",
      x = min(data_time$start_year), y = 3e6 * 0.8,
      label = "Catastrophic PLE threshold", color = "grey50", hjust = 0, fontface = "italic", size = 3
    ) +


    # Add points with size proportional to deaths and color by period
    geom_point(aes(
      color = period, size = deaths_scaled, fill = period
    ), shape = 21, stroke = 1, alpha = 0.7) +

    # Add text labels for significant events
    ggrepel::geom_text_repel(
      data = data_time %>% filter(deaths_scaled > label_threshold),
      aes(label = event_label, color = period),
      box.padding = 0.8,
      point.padding = 0.5,
      min.segment.length = 0.2,
      max.overlaps = 20,
      segment.color = "gray40",
      segment.alpha = 0.7,
      fontface = "bold",
      size = 3,
      show.legend = FALSE
    ) +

    # Set scales
    scale_color_manual(values = period_colors, name = "Time Period") +
    scale_fill_manual(values = period_colors, guide = "none") +
    scale_size_continuous(range = point_size_range, guide = "none") +
    scale_y_log10(
      name = "2025-equivalent deaths (log scale)",
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = c(1e5, 1e6, 1e7, 1e8, 1e9),
      limits = y_range_scaled
    ) +
    scale_x_continuous(
      name = "Year",
      breaks = seq(min_year, 2025, by = 50), # Use evenly spaced, human-readable tick marks
      minor_breaks = NULL, # Remove minor breaks
      expand = expansion(mult = c(0.02, 0.02))
    ) +

    # Set labels conditionally
    labs(
      title = if (include_labels) "Historical Pandemic Mortality Impact (1600-Present)" else NULL,
      subtitle = if (include_labels) paste0("Death tolls scaled to ", reference_year, " world population") else NULL,
      x = "Year",
      caption = if (include_labels) {
        paste0(
          "Note: Events with < ", scales::label_number(scale_cut = scales::cut_short_scale())(label_threshold),
          " deaths may not be labeled. Size of points corresponds to magnitude of impact."
        )
      } else {
        NULL
      }
    ) +

    theme_nature() +
    guides(color = guide_legend(override.aes = list(
      shape = 21,
      fill = period_colors,
      size = 4,
      stroke = 1.2
    )))

  if (show_raw_panel) {
    # Use full-data ranges (including dropped events) so axis reflects
    # where smaller events would have been
    y_limits <- range(c(y_range_raw, y_range_scaled), na.rm = TRUE)
    catastrophic_thresh <- if (!is.null(conf$thresholds$deaths_scaled)) conf$thresholds$deaths_scaled else 3e6

    # Build raw deaths panel (Panel A) with same y scale as scaled panel
    p_raw <- ggplot(data_time, aes(x = start_year, y = deaths)) +
      geom_vline(
        xintercept = c(1800, 1900, 1950, 2000), linetype = "dotted",
        color = "gray50", linewidth = 0.25
      ) +
      geom_hline(yintercept = lower_cutoff,
        linetype = "dashed", color = "grey50", linewidth = 0.25
      ) +
      annotate("text",
        x = min(data_time$start_year), y = lower_cutoff * 0.8,
        label = "Analysis threshold", color = "grey50", hjust = 0, fontface = "italic", size = 3
      ) +
      geom_hline(yintercept = catastrophic_thresh,
        linetype = "dashed", color = "grey50", linewidth = 0.25
      ) +
      annotate("text",
        x = min(data_time$start_year), y = catastrophic_thresh * 0.8,
        label = "Catastrophic PLE threshold", color = "grey50", hjust = 0, fontface = "italic", size = 3
      ) +
      geom_point(aes(color = period, size = deaths_scaled, fill = period),
                 shape = 21, stroke = 1, alpha = 0.7) +
      scale_color_manual(values = period_colors, name = "Time Period") +
      scale_fill_manual(values = period_colors, guide = "none") +
      scale_size_continuous(range = point_size_range, guide = "none") +
      scale_y_log10(
        name = "Reported deaths (log scale)",
        labels = scales::label_number(scale_cut = scales::cut_short_scale()),
        breaks = c(1e5, 1e6, 1e7, 1e8, 1e9),
        limits = y_limits
      ) +
      scale_x_continuous(
        name = "Year",
        breaks = seq(min_year, 2025, by = 50),
        minor_breaks = NULL,
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      labs(title = NULL) +
      theme_nature() +
      # Match p_scaled's legend keys (filled circles) so patchwork dedupes to one guide
      guides(color = guide_legend(override.aes = list(
        shape = 21, fill = period_colors, size = 4, stroke = 1.2
      )))

    # Rename scaled panel as Panel B. Keep its legend (identical override applied in p)
    # so patchwork's guides = "collect" deduplicates the two into a single legend.
    p_scaled <- p + labs(title = NULL)

    # Combine with patchwork (legend collected and placed at the bottom)
    dual <- p_raw + p_scaled +
      patchwork::plot_layout(ncol = 2, guides = "collect") &
      theme(legend.position = "bottom")
    dual <- dual + patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(family = "Arial", size = 12, face = "bold")
      )
    )

    save_plot_if_enabled(dual, "time_scaled_deaths_dual", output_dir, conf,
                         width = 16, height = 8)
  }

  # Always save single-panel version
  save_plot_if_enabled(p, "time_scaled_deaths", output_dir, conf)

  return(p)
}

# Sensitivity plots

#' Format number as label with appropriate suffix (M for millions, k for thousands)
#' @param value Numeric value to format
#' @return Character string with formatted label
format_number_label <- function(value) {
  if (value >= 1e6) {
    # For millions, show one decimal place if not a whole number
    millions <- value / 1e6
    if (millions == round(millions)) {
      paste0(round(millions), "M")
    } else {
      paste0(round(millions, 1), "M")
    }
  } else if (value >= 1e3) {
    # For thousands, show whole numbers
    paste0(round(value / 1e3), "k")
  } else {
    as.character(round(value))
  }
}

#' Plot cumulative deaths sensitivity analysis with enhanced design
#' @param results Data frame with sensitivity analysis results
#' @param var Variable analyzed ("deaths" or "deaths_scaled")
#' @param analysis_type Type of analysis ("year", "threshold", or "cutoff")
#' @param config Configuration object containing analysis parameters
plot_cumulative_deaths <- function(results, var, analysis_type, config = NULL, output_dir = NULL, conf = NULL) {
  # Define the x-axis column based on analysis type
  x_col <- case_when(
    analysis_type == "year" ~ "year",
    analysis_type == "threshold" ~ "threshold",
    analysis_type == "cutoff" ~ "cutoff",
    TRUE ~ "value"
  )

  # Only proceed if we have the required columns
  req_cols <- c(x_col, "cum_deaths", "cum_deaths_low", "cum_deaths_up")
  if (!all(req_cols %in% colnames(results))) {
    warning(
      "Missing required columns: ",
      paste(setdiff(req_cols, colnames(results)), collapse = ", ")
    )
    return(NULL)
  }

  # Define professional color palette
  main_color <- "#377EB8" # Blue for main line
  ribbon_color <- "#C6DBEF" # Light blue for uncertainty ribbon
  accent_color <- "gray50" # Gray for vertical line indicators
  grid_color <- "gray95" # Light gray for grid lines

  # Display unit based on values
  max_value <- max(results$cum_deaths_up, na.rm = TRUE)
  if (max_value > 5000) {
    scale_factor <- "billions"
    results <- results %>%
      mutate(across(c(cum_deaths, cum_deaths_low, cum_deaths_up), ~ .x / 1000))
  } else {
    scale_factor <- "millions"
  }

  # Generate x-axis label
  x_label <- case_when(
    analysis_type == "year" ~ "Year",
    analysis_type == "threshold" ~ "Upper threshold (2025-equivalent deaths)",
    analysis_type == "cutoff" ~ "Lower cutoff (2025-equivalent deaths)",
    TRUE ~ stringr::str_to_title(analysis_type)
  )

  # Generate title with variable formatting
  title <- case_when(
    var == "deaths" ~ "Pandemic Cumulative Death Burden Sensitivity",
    var == "deaths_scaled" ~ "2025-equivalent Pandemic Death Burden Sensitivity",
    TRUE ~ paste0("Sensitivity Analysis for ", stringr::str_to_title(var))
  )

  # Generate subtitle with analysis type information
  subtitle <- case_when(
    analysis_type == "year" ~ "Impact of varying the historical time window on cumulative death estimates",
    analysis_type == "threshold" ~ "Impact of threshold selection on cumulative death estimates",
    analysis_type == "cutoff" ~ "Impact of lower cutoff selection on cumulative death estimates",
    TRUE ~ paste0("Sensitivity by ", stringr::str_to_title(analysis_type))
  )

  # Generate formatted caption
  caption <- paste0(
    "Note: Shaded area represents 95% confidence interval. ",
    if (analysis_type %in% c("threshold", "cutoff")) "Both axes use logarithmic scale. " else "Y-axis uses logarithmic scale. ",
    "Values are shown in ", scale_factor, "."
  )

  # Build plot with professional styling
  p <- ggplot(results, aes(x = .data[[x_col]])) +
    # Add ribbon for uncertainty bounds with semi-transparent fill
    geom_ribbon(
      aes(
        ymin = cum_deaths_low,
        ymax = cum_deaths_up
      ),
      fill = ribbon_color,
      alpha = 0.6
    ) +
    # Add line with purposeful styling
    geom_line(
      aes(y = cum_deaths),
      color = main_color,
      linewidth = 1
    ) +
    # Add points for emphasis
    geom_point(
      aes(y = cum_deaths),
      color = main_color,
      size = 2,
      alpha = 0.7
    ) +
    # Set appropriate scales
    scale_y_continuous(
      name = paste0("Cumulative deaths (", scale_factor, ", log scale)"),
      labels = scales::label_number(accuracy = 1),
      breaks = if (!is.null(config) && !is.null(config$plots$log_breaks_major)) config$plots$log_breaks_major else c(1, 10, 100, 1000),
      trans = "log1p",
      limits = c(0, max(results$cum_deaths_up, na.rm = TRUE) * 1.5)
    ) +
    # Apply conditional x-axis scaling
    {
      if (analysis_type %in% c("threshold", "cutoff")) {
        scale_x_log10(
          name = x_label,
          labels = scales::label_number(accuracy = 1, big.mark = " "),
          breaks = scales::breaks_log(n = 5)
        )
      } else {
        scale_x_continuous(
          name = x_label,
          breaks = scales::pretty_breaks(n = 8)
        )
      }
    } +
    # Apply professional styling
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    # Enhanced theme with attention to detail
    theme_nature() +
    theme(
      # Typography
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20), color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray30", hjust = 0, margin = margin(t = 15)),

      # Axis styling
      axis.title.x = element_text(margin = margin(t = 10), size = 12),
      axis.title.y = element_text(margin = margin(r = 10), size = 12),
      axis.text = element_text(size = 10),

      # Panel and border styling for distinct plotting area
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),

      # Grid styling for cleaner appearance
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = grid_color, linewidth = 0.3),

      # Overall spacing
      plot.margin = margin(20, 20, 20, 20)
    )


  # Add reference line for significant value if available in config
  if (analysis_type == "year" && !is.null(results$year) && !is.null(config)) {
    baseline_year <- if (!is.null(config$sens$year$cutoff)) config$sens$year$cutoff else config$sens$year$base_year
    if (!is.null(baseline_year)) {
      p <- p +
        geom_vline(xintercept = baseline_year, linetype = "dashed", color = accent_color, linewidth = 0.5) +
        annotate("text",
          x = baseline_year, y = min(results$cum_deaths_low, na.rm = TRUE) * 1.5,
          label = baseline_year, color = accent_color, hjust = -0.2, vjust = 0, size = 3, fontface = "italic"
        )
    }
  }
  if (analysis_type == "threshold" && !is.null(results$threshold) && !is.null(config)) {
    baseline_threshold <- config$sens$threshold$base_thresholds$deaths_scaled
    if (!is.null(baseline_threshold)) {
      threshold_label <- format_number_label(baseline_threshold)
      p <- p +
        geom_vline(xintercept = baseline_threshold, linetype = "dashed", color = accent_color, linewidth = 0.5) +
        annotate("text",
          x = baseline_threshold, y = min(results$cum_deaths_low, na.rm = TRUE) * 0.5,
          label = threshold_label, color = accent_color, hjust = -0.2, vjust = 0, size = 3, fontface = "italic"
        )
    }
  }
  if (analysis_type == "cutoff" && !is.null(results$cutoff) && !is.null(config)) {
    baseline_cutoff <- config$thresholds$lower_cutoff_deaths_scaled
    if (!is.null(baseline_cutoff)) {
      cutoff_label <- format_number_label(baseline_cutoff)
      p <- p +
        geom_vline(xintercept = baseline_cutoff, linetype = "dashed", color = accent_color, linewidth = 0.5) +
        annotate("text",
          x = baseline_cutoff, y = min(results$cum_deaths_low, na.rm = TRUE) * 1.5,
          label = cutoff_label, color = accent_color, hjust = -0.2, vjust = 0, size = 3, fontface = "italic"
        )
    }
  }

  # Save plot if output directory is provided
  save_plot_if_enabled(p, paste0("cumulative_deaths_", var, "_", analysis_type), output_dir, conf)
  
  # Also print for interactive use
  tryCatch(print(p), error = function(e) message("Note: Could not display plot: ", e$message))
  
  invisible(p)
}

#' Plot Severity Distribution by Era and Tier
#'
#' Creates a two-facet boxplot with jittered event points showing
#' population-scaled death distributions across historical eras,
#' split by severity tier. Demonstrates stationarity of severity
#' distribution over time.
#'
#' @param data Data frame with pandemic events
#' @param variable Column name for population-scaled deaths
#' @param conf Configuration object (uses thresholds)
#' @param output_dir Output directory for saving
plot_severity_by_era <- function(data, variable, conf, output_dir = NULL) {
  lower_cutoff <- conf$thresholds$lower_cutoff %||%
    conf$thresholds[[paste0("lower_cutoff_", variable)]] %||% 1e5
  tier_boundary <- conf$thresholds[[variable]] %||% 3e6
  base_year <- conf$base_year %||% 1600

  plot_data <- data %>%
    dplyr::filter(start_year >= base_year, .data[[variable]] >= lower_cutoff) %>%
    dplyr::mutate(
      era = dplyr::case_when(
        start_year < 1800 ~ "Pre-1800",
        start_year < 1900 ~ "19th Century",
        start_year < 1950 ~ "1900\u20131950",
        TRUE ~ "1950\u20132025"
      ),
      era = factor(era, levels = c("Pre-1800", "19th Century",
                                    "1900\u20131950", "1950\u20132025")),
      tier = ifelse(.data[[variable]] > tier_boundary,
                    "Catastrophic PLE (>3M 2025-equivalent deaths)",
                    "Moderate PLE (100K\u20133M 2025-equivalent deaths)"),
      tier = factor(tier, levels = c("Moderate PLE (100K\u20133M 2025-equivalent deaths)",
                                      "Catastrophic PLE (>3M 2025-equivalent deaths)"))
    )

  tier_colors <- c(
    "Moderate PLE (100K\u20133M 2025-equivalent deaths)" = "#00468B",
    "Catastrophic PLE (>3M 2025-equivalent deaths)" = "#ED0000"
  )

  counts <- plot_data %>%
    dplyr::count(era, tier) %>%
    dplyr::rename(n_events = n)

  y_mins <- plot_data %>%
    dplyr::group_by(tier) %>%
    dplyr::summarise(y_min = min(.data[[variable]]), .groups = "drop")
  counts <- dplyr::left_join(counts, y_mins, by = "tier")

  # Box plots are drawn only for groups with n >= 3; every event is still drawn as a point.
  box_data <- plot_data %>%
    dplyr::add_count(era, tier, name = "n_group") %>%
    dplyr::filter(n_group >= 3)

  set.seed(42)
  p <- ggplot(plot_data, aes(x = era, y = .data[[variable]])) +
    geom_boxplot(
      data = box_data,
      aes(group = era),
      fill = "white",
      colour = "black",
      linewidth = 0.4,
      width = 0.6,
      outlier.shape = NA
    ) +
    geom_jitter(
      aes(colour = tier),
      width = 0.18,
      size = 2.2,
      alpha = 0.7,
      shape = 16
    ) +
    geom_text(
      data = counts,
      aes(x = era, y = y_min * 0.45, label = paste0("n=", n_events)),
      size = 2.8,
      colour = "grey45",
      family = "Arial"
    ) +
    facet_wrap(~ tier, ncol = 1, scales = "free_y") +
    scale_y_log10(
      labels = function(x) {
        ifelse(x >= 1e6,
               paste0(x / 1e6, "M"),
               paste0(x / 1e3, "K"))
      }
    ) +
    scale_colour_manual(values = tier_colors, guide = "none") +
    labs(x = NULL, y = "2025-equivalent deaths") +
    theme_nature() +
    theme(
      strip.text = element_text(family = "Arial", size = 9, face = "bold",
                                 hjust = 0, margin = margin(b = 5)),
      axis.title.y = element_text(family = "Arial", size = 9, margin = margin(r = 10)),
      axis.text.x = element_text(family = "Arial", size = 8),
      axis.text.y = element_text(family = "Arial", size = 8),
      axis.ticks = element_line(linewidth = 0.5),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(1.2, "lines"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15)
    )

  save_plot_if_enabled(p, "severity_by_era", output_dir, conf)
  tryCatch(print(p), error = function(e) message("Note: Could not display plot: ", e$message))
  invisible(p)
}

#' Plot tiered exceedance rates (Moderate vs Catastrophic PLE)
#'
#' Two-panel faceted figure showing observed rates, fitted trend,
#' uncertainty ribbon, and forecast extrapolation per tier.
#'
#' @param tier1_rate_data Rate data frame for Moderate tier
#' @param tier2_rate_data Rate data frame for Catastrophic tier
#' @param conf Configuration object
#' @param output_dir Output directory for saving
#' @return ggplot object
plot_tiered_rates <- function(tier1_rate_data, tier2_rate_data, conf = NULL,
                              output_dir = NULL) {
  # Okabe-Ito palette (colour-blind safe).
  rate_colors <- c(
    "Moderate PLE"     = "#0072B2",  # Okabe-Ito blue
    "Catastrophic PLE" = "#D55E00"   # Okabe-Ito vermillion
  )

  # Make tier an explicit factor so the panel order (Catastrophic on top,
  # Moderate on bottom) is stable regardless of input row order.
  tier_levels <- c("Catastrophic PLE", "Moderate PLE")
  tier1_rd <- tier1_rate_data %>% dplyr::mutate(tier = factor("Moderate PLE",     levels = tier_levels))
  tier2_rd <- tier2_rate_data %>% dplyr::mutate(tier = factor("Catastrophic PLE", levels = tier_levels))
  combined_rd <- dplyr::bind_rows(tier1_rd, tier2_rd)

  # Use smooth model-implied rate for the central line when available;
  # fall back to `fitted_rate` otherwise. The 95% ribbon is unchanged.
  if (!"fitted_rate_smooth" %in% names(combined_rd)) {
    combined_rd$fitted_rate_smooth <- combined_rd$fitted_rate
  } else {
    combined_rd$fitted_rate_smooth <- dplyr::coalesce(
      combined_rd$fitted_rate_smooth, combined_rd$fitted_rate
    )
  }

  current_year <- conf$time$current_year %||% 2025
  breakpoint_year <- conf$time$breakpoint_year %||% 1950

  year_breaks <- seq(
    floor(min(c(tier1_rd$year, tier2_rd$year)) / 50) * 50,
    ceiling(max(c(tier1_rd$year, tier2_rd$year)) / 50) * 50,
    by = 50
  )

  # Panel letters a/b for caption references. Placed via geom_text
  # at -Inf/Inf so they sit in each panel's top-left corner after facetting.
  panel_tags <- data.frame(
    tier  = factor(c("Catastrophic PLE", "Moderate PLE"), levels = tier_levels),
    label = c("a", "b")
  )

  p <- ggplot(combined_rd, aes(x = year)) +
    # Forecast region: faint grey shading from current_year onward (both panels).
    annotate("rect",
             xmin = current_year, xmax = Inf,
             ymin = -Inf, ymax = Inf,
             fill = "grey92", alpha = 0.5) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = tier),
                alpha = 0.28, colour = NA) +
    geom_point(
      data = combined_rd %>% dplyr::filter(type == "historical"),
      aes(y = observed_rate), colour = "grey32", alpha = 0.35, size = 0.8
    ) +
    geom_line(aes(y = fitted_rate_smooth, colour = tier), linewidth = 0.7) +
    # Piecewise model breakpoint (Moderate panel only). Unified grey50 to match
    # the rest of the auxiliary chrome; label sits just inside the panel top.
    geom_vline(
      data = data.frame(tier = factor("Moderate PLE", levels = tier_levels)),
      aes(xintercept = breakpoint_year),
      linetype = "dotted", linewidth = 0.4, color = "grey50"
    ) +
    geom_text(
      data = data.frame(
        tier = factor("Moderate PLE", levels = tier_levels),
        x = breakpoint_year + 2,
        y = Inf,
        label = paste0("Breakpoint: ", breakpoint_year)
      ),
      aes(x = x, y = y, label = label),
      hjust = 0, vjust = 1.5, size = 8 / .pt, family = "Arial", color = "grey50"
    ) +
    geom_text(
      data = panel_tags,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.5, vjust = 1.3, size = 8 / .pt, fontface = "bold", family = "Arial"
    ) +
    facet_wrap(~ tier, ncol = 1, scales = "free_y",
               axes = "all_x", axis.labels = "margins") +
    scale_colour_manual(values = rate_colors) +
    scale_fill_manual(values = rate_colors) +
    scale_x_continuous(breaks = year_breaks) +
    labs(x = "Year", y = "Exceedance rate (events per year)") +
    theme_nature() +
    theme(legend.position = "none")

  save_plot_if_enabled(p, "tiered_rate_comparison", output_dir, conf, scale = 0.8)
  invisible(p)
}


#' Plot severity over time (supplementary figure)
#'
#' Scatter of 2025-equivalent deaths by start year, colored by pre/post-1950.
#' Supports the fixed-severity assumption.
#'
#' @param data Data frame with deaths_scaled and start_year columns
#' @param breakpoint_year Year to split eras (default 1950)
#' @param output_dir Output directory
#' @param conf Configuration object
#' @return ggplot object
plot_severity_over_time <- function(data, breakpoint_year = 1950,
                                    output_dir = NULL, conf = NULL) {
  colors_period <- c("Pre-1950" = "#00468B", "Post-1950" = "#ED0000")

  plot_data <- data %>%
    dplyr::filter(!is.na(deaths_scaled), start_year >= 1600) %>%
    dplyr::mutate(
      period = ifelse(start_year < breakpoint_year, "Pre-1950", "Post-1950"),
      period = factor(period, levels = c("Pre-1950", "Post-1950"))
    )

  p <- ggplot(plot_data, aes(x = start_year, y = deaths_scaled,
                              color = period, size = deaths_scaled)) +
    geom_point(alpha = 0.7, shape = 21, aes(fill = period), stroke = 0.5) +
    geom_vline(xintercept = breakpoint_year, linetype = "dashed",
               color = "gray40", linewidth = 0.5) +
    annotate("text", x = breakpoint_year, y = max(plot_data$deaths_scaled) * 0.9,
             label = as.character(breakpoint_year), hjust = -0.2,
             fontface = "bold", family = "Arial", size = 3) +
    scale_color_manual(values = colors_period, name = "Period") +
    scale_fill_manual(values = colors_period, guide = "none") +
    scale_y_log10(
      name = "2025-equivalent deaths (log scale)",
      labels = scales::label_number(scale = 1e-6, suffix = "M")
    ) +
    scale_size_continuous(range = c(2, 10), guide = "none") +
    scale_x_continuous(
      name = "Start year",
      breaks = seq(1600, 2025, by = 50)
    ) +
    theme_nature()

  save_plot_if_enabled(p, "severity_over_time", output_dir, conf)
  invisible(p)
}
