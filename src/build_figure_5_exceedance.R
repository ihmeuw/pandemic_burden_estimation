#!/usr/bin/env Rscript
# =============================================================================
# build_figure_5_exceedance.R
#
# Build main-text Figure 5: probability that cumulative pandemic deaths reach
# reference totals, 2026-2100, under best-plausible / reference /
# worst-plausible scenarios.
#
# - Panel a: empirical survival curve S(x) = P(X > x) on a log10 x-axis
#   from 25M to 1B cumulative deaths.
# - Panel b: probability matrix (grey-ramp tiles labelled with percentages) at five reference totals
#   25M / 50M / 100M / 200M / 500M.
#
# Inputs (sibling run_all directory):
#   <run_all>/run_pot/cum_draws_{best,reference,worst}.rds
#   <run_all>/run_yaml_export/_variables.yml  (validation)
#
# Outputs (under setup_script_run() session dir):
#   <output_dir>/plots/pandemics_figure5.{png,svg}
#   <output_dir>/tables/pandemics_figure5_{curve,threshold_matrix}.csv
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
  library(yaml)
})

source("R/load_all.R")
Sys.setenv(R_CONFIG_ACTIVE = "default")
conf <- load_config()

run_info <- setup_script_run("build_figure_5_exceedance", conf)
output_dir <- run_info$output_dir

# Locate the run_all session containing this script's prerequisites.
# Under run_all_scripts(): output_dir is .../run_all_<ts>/build_figure_5_exceedance,
#                         so the sibling run_pot is at dirname(output_dir)/run_pot.
# Standalone: pick the latest run_all_* directory.
parent_dir <- dirname(output_dir)
run_dir <- if (grepl("/run_all_[0-9]{4}-", parent_dir)) {
  parent_dir
} else {
  candidates <- list.dirs(file.path(getwd(), "output"), recursive = FALSE)
  candidates <- candidates[grepl("/run_all_[0-9]{4}-", candidates)]
  if (length(candidates) == 0) {
    stop("No run_all_<timestamp> directories found")
  }
  sort(candidates)[length(candidates)]
}

pot_dir         <- file.path(run_dir, "run_pot")
yaml_export_dir <- file.path(run_dir, "run_yaml_export")

stopifnot(
  "run_pot/ not found in selected run dir" =
    dir.exists(pot_dir),
  "cum_draws_reference.rds not found" =
    file.exists(file.path(pot_dir, "cum_draws_reference.rds"))
)

message("Using run dir: ", run_dir)

# ---- Load draw-level cumulative deaths --------------------------------------

read_draws <- function(name, file) {
  vec <- readRDS(file)
  if (!is.numeric(vec)) stop("Expected numeric draws, got ", class(vec)[1])
  vec <- as.numeric(vec)
  vec <- vec[is.finite(vec) & vec > 0]
  message(sprintf("  %-15s n=%d  median=%.2fM  range=[%.2fM, %.2fM]",
                  name, length(vec), median(vec) / 1e6,
                  min(vec) / 1e6, max(vec) / 1e6))
  vec
}

message("Loading draw-level cumulative deaths:")
draws <- list(
  best_plausible  = read_draws("best",      file.path(pot_dir, "cum_draws_best.rds")),
  reference       = read_draws("reference", file.path(pot_dir, "cum_draws_reference.rds")),
  worst_plausible = read_draws("worst",     file.path(pot_dir, "cum_draws_worst.rds"))
)

# ---- Sanity-check units: must be raw deaths in the 1e7-1e9 range -----------

ref_med <- median(draws$reference)
if (ref_med >= 50 && ref_med <= 500) {
  warning("Reference median ", round(ref_med, 1),
          " looks like millions, not raw deaths. Multiplying by 1e6.")
  draws <- lapply(draws, function(v) v * 1e6)
} else if (ref_med < 5e7 || ref_med > 5e8) {
  stop("Reference median is ", ref_med,
       " which is outside expected raw-deaths range [5e7, 5e8]. Inspect units.")
}

# ---- Compute curve data (Panel a) -------------------------------------------

x_grid_deaths <- exp(seq(log(25e6), log(1e9), length.out = 400))

curve_data <- bind_rows(lapply(names(draws), function(scen) {
  v <- draws[[scen]]
  tibble(
    scenario              = scen,
    threshold_deaths      = x_grid_deaths,
    threshold_millions    = x_grid_deaths / 1e6,
    exceedance_probability = vapply(x_grid_deaths,
                                    function(t) mean(v > t),
                                    numeric(1))
  )
}))

# ---- Compute threshold ladder (Panel b) -------------------------------------

fixed_thresholds <- c(25e6, 50e6, 100e6, 200e6, 500e6)
threshold_labels <- c("25 million", "50 million", "100 million",
                      "200 million", "500 million")

format_pct <- function(p) {
  ifelse(p >= 0.995, "≈100%",
         ifelse(p < 0.001, "<0.1%",
                sprintf("%.1f%%", p * 100)))
}

ladder_data <- bind_rows(lapply(names(draws), function(scen) {
  v <- draws[[scen]]
  tibble(
    scenario               = scen,
    threshold_deaths       = fixed_thresholds,
    threshold_label        = threshold_labels,
    exceedance_probability = vapply(fixed_thresholds,
                                    function(t) mean(v > t),
                                    numeric(1))
  )
})) %>%
  mutate(percent_label = format_pct(exceedance_probability))

# ---- Validate against _variables.yml ----------------------------------------

vars_file <- file.path(yaml_export_dir, "_variables.yml")
if (file.exists(vars_file)) {
  vars <- yaml::read_yaml(vars_file)
  expect <- list(
    best_plausible  = vars$probability$exceedance_best,
    reference       = vars$probability$exceedance_reference,
    worst_plausible = vars$probability$exceedance_worst
  )
  thr_keys <- c("exceed_25M", "exceed_50M", "exceed_100M",
                "exceed_200M", "exceed_500M")
  message("Validating against ", vars_file)
  max_diff <- 0
  for (scen in names(expect)) {
    for (i in seq_along(thr_keys)) {
      yaml_val <- as.numeric(expect[[scen]][[thr_keys[i]]])
      if (is.na(yaml_val)) next
      computed <- ladder_data$exceedance_probability[
        ladder_data$scenario == scen &
          ladder_data$threshold_deaths == fixed_thresholds[i]
      ]
      d <- abs(computed - yaml_val)
      max_diff <- max(max_diff, d)
      if (d > 0.005) {
        warning(sprintf(
          "Mismatch %s @ %s: computed=%.4f yaml=%.4f (diff=%.4f)",
          scen, thr_keys[i], computed, yaml_val, d))
      }
    }
  }
  message(sprintf("  max abs diff vs _variables.yml: %.5f", max_diff))
} else {
  warning("Skipping validation: ", vars_file, " not found")
}

# ---- Write source-data CSVs -------------------------------------------------

# setup_script_run() already created plots/ and tables/ subdirectories.
curve_csv  <- file.path(output_dir, "tables", "pandemics_figure5_curve.csv")
matrix_csv <- file.path(output_dir, "tables", "pandemics_figure5_threshold_matrix.csv")

write_csv(curve_data,  curve_csv)
write_csv(ladder_data, matrix_csv)
message("Wrote ", curve_csv)
message("Wrote ", matrix_csv)

# ---- Plotting ---------------------------------------------------------------

# Muted, colour-blind-safe palette
scenario_levels <- c("best_plausible", "reference", "worst_plausible")
scenario_pretty <- c(best_plausible  = "Best-plausible",
                     reference       = "Reference",
                     worst_plausible = "Worst-plausible")
scenario_colors <- c(best_plausible  = "#4A7BAA",   # muted blue
                     reference       = "#2F2F2F",   # charcoal
                     worst_plausible = "#B5482A")   # muted vermillion

curve_data <- curve_data %>%
  mutate(scenario = factor(scenario, levels = scenario_levels))
ladder_data <- ladder_data %>%
  mutate(
    scenario        = factor(scenario, levels = scenario_levels),
    threshold_label = factor(threshold_label, levels = threshold_labels)
  )

theme_fig5 <- function(base_size = 9) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "Arial"),
      axis.title       = element_text(size = base_size + 0.5),
      axis.text        = element_text(size = base_size - 0.5, color = "black"),
      axis.line.x      = element_line(color = "black", linewidth = 0.4),
      axis.line.y      = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = base_size, face = "plain",
                                      hjust = 0.5, margin = margin(b = 6)),
      plot.tag         = element_text(size = base_size + 2, face = "bold",
                                      hjust = 0, vjust = 1),
      plot.tag.position = c(0.005, 0.99),
      legend.position  = "none",
      plot.margin      = margin(8, 10, 6, 10)
    )
}

# --- Panel a: exceedance curve with direct labels ----------------------------

# Hand-tuned anchor positions: place each scenario's label on the steep
# descent of its own curve. Coordinates chosen so the three labels are well
# separated horizontally (each curve descends at a different x) and sit just
# above the curve.
direct_label_data <- tibble(
  scenario     = factor(scenario_levels, levels = scenario_levels),
  pretty_label = unname(scenario_pretty[scenario_levels]),
  label_x      = c(25,   130,  300),    # millions
  label_y      = c(0.77, 0.52, 0.37),   # probability
  label_hjust  = c(0,    0.5,  0.5)     # left-align Best at x=25, centre others
)

panel_a <- ggplot(curve_data,
                  aes(x = threshold_millions,
                      y = exceedance_probability,
                      colour = scenario)) +
  geom_vline(xintercept = c(100, 200),
             colour = "grey60", linewidth = 0.35, linetype = "22") +
  annotate("text", x = 100, y = 1.02, label = "100 million",
           colour = "grey35", size = 2.5, family = "Arial",
           hjust = 1.05, vjust = 0) +
  annotate("text", x = 200, y = 1.02, label = "200 million",
           colour = "grey35", size = 2.5, family = "Arial",
           hjust = -0.05, vjust = 0) +
  geom_line(linewidth = 1.0) +
  geom_text(data = direct_label_data,
            aes(x = label_x, y = label_y,
                label = pretty_label, colour = scenario,
                hjust = label_hjust),
            family = "Arial", size = 2.95, fontface = "bold",
            vjust = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = scenario_colors, guide = "none") +
  scale_x_log10(
    name   = "Cumulative pandemic deaths, 2026–2100 (millions)",
    breaks = c(25, 50, 100, 200, 500, 1000),
    labels = c("25", "50", "100", "200", "500", "1,000"),
    limits = c(25, 1000),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name   = "P(cumulative deaths ≥ x)",
    limits = c(0, 1.07),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  coord_cartesian(clip = "off") +
  labs(tag = "a") +
  theme_fig5()

# --- Panel b: threshold-risk probability matrix ------------------------------

# Ensure cells stack 25M at top, 500M at bottom; columns left→right.
# Wrap long scenario names onto two lines so column headers don't collide.
scenario_wrapped <- c(best_plausible  = "Best-\nplausible",
                      reference       = "Reference",
                      worst_plausible = "Worst-\nplausible")
matrix_data <- ladder_data %>%
  mutate(
    scenario_col  = factor(scenario_wrapped[as.character(scenario)],
                           levels = unname(scenario_wrapped[scenario_levels])),
    threshold_row = factor(threshold_label, levels = rev(threshold_labels))
  )

# Subtle grey ramp; never dark enough to require white text
matrix_fill <- function(p) {
  # map probability [0,1] to grey ramp from 97% to 78% lightness
  scales::seq_gradient_pal(low = "grey97", high = "grey78")(p)
}
matrix_data$fill_col <- matrix_fill(matrix_data$exceedance_probability)

panel_b <- ggplot(matrix_data,
                  aes(x = scenario_col, y = threshold_row)) +
  geom_tile(aes(fill = fill_col),
            colour = "grey80", linewidth = 0.35,
            width = 1, height = 1) +
  geom_text(aes(label = percent_label),
            family = "Arial", size = 2.85, colour = "black") +
  scale_fill_identity() +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Probability by threshold",
       tag = "b") +
  theme_fig5() +
  theme(
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    axis.text.x.top = element_text(face = "bold", size = 8,
                                   lineheight = 0.9,
                                   margin = margin(b = 4), colour = "black"),
    axis.text.y     = element_text(size = 8.2,
                                   margin = margin(r = 4), colour = "black"),
    panel.grid      = element_blank(),
    plot.title      = element_text(size = 8.5, hjust = 0.5,
                                   margin = margin(b = 8))
  )

# --- Combine -----------------------------------------------------------------

figure_5 <- panel_a + panel_b +
  plot_layout(widths = c(2.05, 1))   # ~67% / 33%

# --- Save --------------------------------------------------------------------

png_path <- file.path(output_dir, "plots", "pandemics_figure5.png")
svg_path <- file.path(output_dir, "plots", "pandemics_figure5.svg")

# Two-column width is ~180mm = 7.08in; trim height since no legend.
# `scale` > 1 makes text/lines render smaller relative to the plot area
# without touching any theme/font sizes (renders on a larger virtual canvas
# and rescales to the requested width/height). 1.4 = ~30% smaller elements.
fig_w     <- 7.2
fig_h     <- 3.2
fig_scale <- 1.4

ggsave(png_path, figure_5, width = fig_w, height = fig_h,
       scale = fig_scale, dpi = 450,
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
ggsave(svg_path, figure_5, width = fig_w, height = fig_h,
       scale = fig_scale,
       device = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else NULL)

message("Wrote ", png_path)
message("Wrote ", svg_path)

# ---- Quality-check log ------------------------------------------------------

stopifnot(
  "Reference curve must be monotonically non-increasing" =
    all(diff(curve_data$exceedance_probability[
      curve_data$scenario == "reference"]) <= 1e-12),
  "Best curve must be below reference at all thresholds" =
    all(
      curve_data$exceedance_probability[curve_data$scenario == "best_plausible"] <=
        curve_data$exceedance_probability[curve_data$scenario == "reference"] + 1e-9
    ),
  "Worst curve must be above reference at all thresholds" =
    all(
      curve_data$exceedance_probability[curve_data$scenario == "worst_plausible"] >=
        curve_data$exceedance_probability[curve_data$scenario == "reference"] - 1e-9
    )
)

message("Figure 5 build complete.")
