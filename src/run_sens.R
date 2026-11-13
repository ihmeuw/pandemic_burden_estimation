####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Orchestration script for pandemic model sensitivity analysis.
##              This script coordinates the entire sensitivity analysis workflow by:
##              1. Setting up the environment with sensitivity configuration
##              2. Running sensitivity analyses across multiple dimensions:
##                 a. Time period sensitivity - how results change with data inclusion years
##                 b. Threshold sensitivity - impact of varying the tail threshold
##                 c. Cutoff sensitivity - effect of different lower cutoff values
##              3. Generating visualizations for each sensitivity dimension
##
## Requires:    - Configuration file ("config.yml") with sensitivity parameters
##
## Outputs:     - Sensitivity analysis results for each dimension
##              - Visualizations showing parameter impact on:
##                  * Expected annual burden estimates
##                  * Model parameter stability
##                  * Cumulative death projections
##              - Plots
##
####################################################################################################

# Load project functions & packages explicitly
source("R/load_all.R")

# ===================================================================================
# Setup environment and configuration
# ===================================================================================

Sys.setenv(R_CONFIG_ACTIVE = "sensitivity")
conf <- load_config()

# Read and transform the pandemic data
data <- read_transform_data(conf)

# ===================================================================================
# Setup output management
# ===================================================================================
# Setup script-specific timestamped output directory
run_info <- setup_script_run("run_sens", conf)
output_dir <- run_info$output_dir

# Ensure required packages
for (pkg in c("future", "future.apply", "parallel", "viridis")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Package", pkg, "is not available."))
  }
}

# Set parallelization strategy from config
parallel_strategy <- ifelse(
  isTRUE(conf$parallel$enabled),
  conf$parallel$strategy,
  "none"
)
workers <- conf$parallel$workers # NULL means use detectCores() - 1

# Log parallelization settings
cat(
  "Parallelization:",
  ifelse(parallel_strategy == "none", "Disabled", "Enabled"),
  ifelse(parallel_strategy != "none", paste0("using '", parallel_strategy, "' strategy"), ""),
  "\n"
)

# Performance timing - start time
start_time <- Sys.time()
cat("Starting sensitivity analysis at:", format(start_time), "\n")

# ===================================================================================
# Time period sensitivity analysis
# ===================================================================================
cat("Running time period sensitivity analysis...\n")

# Generate year range for sensitivity testing
years <- get_year_range(data, conf)

# Configure analysis parameters
sensitivity_grid <- expand.grid(
  analysis_type = "year",
  plot_type = c("cumulative"),
  var = conf$variables,
  stringsAsFactors = FALSE
)

# Execute year-based sensitivity analysis
results_year <- sensitivity_grid |>
  pmap(
    function(var, analysis_type, plot_type) {
      result <- run_parallel_sensitivity(
        data = data,
        var = var,
        analysis_type = analysis_type,
        values = years,
        conf = conf,
        strategy = parallel_strategy,
        workers = workers,
        output_dir = output_dir
      )

      result
    }
  )

# ===================================================================================
# Threshold sensitivity analysis
# ===================================================================================
cat("Running threshold sensitivity analysis...\n")

# Extract threshold ranges for testing.
# deaths_scaled uses indices [13:55] → ~196K to ~11M, encompassing the canonical
# 200K–10M range cited in the main text. Manuscript Extended Data Figure 3
# (panel b) is built from this sweep.
thresholds <- list(
  deaths        = threshold_extract(data, "deaths",        conf$base_year, conf)[25:55],
  deaths_scaled = threshold_extract(data, "deaths_scaled", conf$base_year, conf)[13:55]
)

# Configure analysis parameters
sensitivity_grid <- expand.grid(
  analysis_type = "threshold",
  plot_type = c("cumulative"),
  var = conf$variables,
  stringsAsFactors = FALSE
)

# Execute threshold-based sensitivity analysis
results_threshold <- sensitivity_grid |>
  pmap(
    function(var, analysis_type, plot_type) {
      result <- run_parallel_sensitivity(
        data = data,
        var = var,
        analysis_type = analysis_type,
        values = thresholds[[var]],
        conf = conf,
        strategy = parallel_strategy,
        workers = workers,
        output_dir = output_dir
      )

      result
    }
  )

# ===================================================================================
# Lower cutoff sensitivity analysis
# ===================================================================================
cat("Running lower cutoff sensitivity analysis...\n")

# Extract cutoff ranges for testing
cutoffs <- list(
  deaths = threshold_extract(data, "deaths", conf$base_year, conf),
  deaths_scaled = threshold_extract(data, "deaths_scaled", conf$base_year, conf)
) |>
  map(~ .x[1:30])

# Configure analysis parameters
sensitivity_grid <- expand.grid(
  analysis_type = "cutoff",
  plot_type = c("cumulative"),
  var = conf$variables,
  stringsAsFactors = FALSE
)

# Execute cutoff-based sensitivity analysis
results_cutoff <- sensitivity_grid |>
  pmap(
    function(var, analysis_type, plot_type) {
      result <- run_parallel_sensitivity(
        data = data,
        var = var,
        analysis_type = analysis_type,
        values = cutoffs[[var]],
        conf = conf,
        strategy = parallel_strategy,
        workers = workers,
        output_dir = output_dir
      )

      result
    }
  )

# ===================================================================================
# Summarise sensitivity results for _variables.yml export
# ===================================================================================
cat("Computing sensitivity summary for YAML export...\n")

# Helper: extract summary from a results list for a given variable
summarise_sens <- function(results_list, grid, target_var, param_col) {
  idx <- which(grid$var == target_var)
  if (length(idx) == 0) return(NULL)
  tbl <- results_list[[idx[1]]]$results
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)

  list(
    min      = min(tbl$cum_deaths, na.rm = TRUE),
    max      = max(tbl$cum_deaths, na.rm = TRUE),
    mean     = mean(tbl$cum_deaths, na.rm = TRUE),
    min_at   = tbl[[param_col]][which.min(tbl$cum_deaths)],
    max_at   = tbl[[param_col]][which.max(tbl$cum_deaths)],
    n_values = nrow(tbl)
  )
}

# Build grids matching result ordering above
cutoff_grid <- expand.grid(
  analysis_type = "cutoff", plot_type = "cumulative",
  var = conf$variables, stringsAsFactors = FALSE
)
threshold_grid <- expand.grid(
  analysis_type = "threshold", plot_type = "cumulative",
  var = conf$variables, stringsAsFactors = FALSE
)
year_grid <- expand.grid(
  analysis_type = "year", plot_type = "cumulative",
  var = conf$variables, stringsAsFactors = FALSE
)

sensitivity_summary <- list(
  cutoff    = summarise_sens(results_cutoff, cutoff_grid, "deaths_scaled", "cutoff"),
  threshold = summarise_sens(results_threshold, threshold_grid, "deaths_scaled", "threshold"),
  year      = summarise_sens(results_year, year_grid, "deaths_scaled", "year")
)

# Save summary for yaml_export
saveRDS(sensitivity_summary, file.path(output_dir, "sensitivity_summary.rds"))
cat("Sensitivity summary saved to:", file.path(output_dir, "sensitivity_summary.rds"), "\n")

# Also save as CSV for human inspection
if (!is.null(output_dir)) {
  summary_df <- data.frame(
    dimension = c("cutoff", "threshold", "year"),
    min_millions = c(sensitivity_summary$cutoff$min,
                     sensitivity_summary$threshold$min,
                     sensitivity_summary$year$min),
    max_millions = c(sensitivity_summary$cutoff$max,
                     sensitivity_summary$threshold$max,
                     sensitivity_summary$year$max),
    mean_millions = c(sensitivity_summary$cutoff$mean,
                      sensitivity_summary$threshold$mean,
                      sensitivity_summary$year$mean),
    min_at = c(sensitivity_summary$cutoff$min_at,
               sensitivity_summary$threshold$min_at,
               sensitivity_summary$year$min_at),
    max_at = c(sensitivity_summary$cutoff$max_at,
               sensitivity_summary$threshold$max_at,
               sensitivity_summary$year$max_at)
  )
  save_table_if_enabled(summary_df, "sensitivity_summary.csv", output_dir, conf)
}

# ===================================================================================
# Multi-panel composite figure for manuscript ED Figure 3
# ===================================================================================
# Combines the three sensitivity ggplots (cutoff, threshold, year) into one
# vertically-stacked composite. Re-builds each panel via plot_cumulative_deaths()
# against the captured results tibbles, applies a tighter y-axis lower bound
# (20 M) to remove dead whitespace below the data, adds patchwork panel labels
# a/b/c, and saves alongside the standalone single-panel outputs.
cat("Generating multi-panel composite (sensitivity_combined)...\n")

# Capture index of deaths_scaled (the manuscript variable) in each results list
idx_cutoff    <- which(expand.grid(analysis_type = "cutoff",    plot_type = "cumulative",
                                   var = conf$variables, stringsAsFactors = FALSE)$var == "deaths_scaled")[1]
idx_threshold <- which(expand.grid(analysis_type = "threshold", plot_type = "cumulative",
                                   var = conf$variables, stringsAsFactors = FALSE)$var == "deaths_scaled")[1]
idx_year      <- which(expand.grid(analysis_type = "year",      plot_type = "cumulative",
                                   var = conf$variables, stringsAsFactors = FALSE)$var == "deaths_scaled")[1]

# Re-build each ggplot (writes the standalone single-panel file again, but
# returns the ggplot object so we can compose). To avoid an extra disk
# write of the standalone files, temporarily clear save_plots — restored
# after the three calls.
.saved_save_plots <- conf$output$save_plots
conf$output$save_plots <- FALSE

p_cutoff    <- plot_cumulative_deaths(results_cutoff[[idx_cutoff]]$results,
                                      "deaths_scaled", "cutoff",
                                      config = conf, output_dir = output_dir, conf = conf)
p_threshold <- plot_cumulative_deaths(results_threshold[[idx_threshold]]$results,
                                      "deaths_scaled", "threshold",
                                      config = conf, output_dir = output_dir, conf = conf)
p_year      <- plot_cumulative_deaths(results_year[[idx_year]]$results,
                                      "deaths_scaled", "year",
                                      config = conf, output_dir = output_dir, conf = conf)

conf$output$save_plots <- .saved_save_plots

# Trim y-axis lower bound to 20 M (data floor is ~25 M; preserves a small margin).
# coord_cartesian zooms without clipping the geom_ribbon data, avoiding warnings.
ylim_lower <- 20
p_cutoff_trim    <- p_cutoff    + coord_cartesian(ylim = c(ylim_lower, NA))
p_threshold_trim <- p_threshold + coord_cartesian(ylim = c(ylim_lower, NA))
p_year_trim      <- p_year      + coord_cartesian(ylim = c(ylim_lower, NA))

# Strip the per-panel title/subtitle/caption blocks (the merged figure carries
# one shared legend in the manuscript caption — internal titles add visual noise
# in a 3-panel vertical stack) and tighten margins.
strip_panel_chrome <- function(p) {
  p +
    labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme(plot.margin = margin(10, 10, 5, 10))
}
p_cutoff_clean    <- strip_panel_chrome(p_cutoff_trim)
p_threshold_clean <- strip_panel_chrome(p_threshold_trim)
p_year_clean      <- strip_panel_chrome(p_year_trim)

# Vertical stack via patchwork with bold lowercase panel tags
sensitivity_combined <- (p_cutoff_clean / p_threshold_clean / p_year_clean) +
  patchwork::plot_annotation(
    tag_levels = list(c("a", "b", "c")),
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

save_plot_if_enabled(sensitivity_combined, "sensitivity_combined", output_dir, conf,
                     width = 7, height = 12)

cat("Multi-panel composite saved as sensitivity_combined.{png,svg}\n")
