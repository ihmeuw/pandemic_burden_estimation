####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: HIV/AIDS-exclusion sensitivity analysis.
##              Re-fits the entire pipeline (severity model + tier-specific rate model with
##              post-1950 breakpoint + burden forecast + out-of-sample validation) on a PLE
##              set with HIV/AIDS removed, and compares projections against the reference run.
##
## Outputs:     <session_dir>/    (session managed by setup_script_run)
##                ├── reference/
##                │     ├── (run_PoT artifacts)
##                │     └── validation/              run_validation outputs under validation cfg
##                ├── hiv_excluded/
##                │     ├── (run_PoT artifacts)
##                │     └── validation/              run_validation outputs (HIV-excluded)
##                ├── plots/hiv_exclusion_comparison.{png,svg}     forecast-period comparison
##                ├── plots/hiv_exclusion_validation.{png,svg}     validation-period comparison
##                ├── tables/hiv_exclusion_summary.csv
##                └── hiv_exclusion_summary.rds      headline deltas for run_yaml_export
##
##              session_dir resolves to:
##                output/run_hiv_sensitivity_<ts>/                 (standalone)
##                output/run_all_<ts>/run_hiv_sensitivity/         (under run_all)
##
####################################################################################################

source("R/load_all.R")

# ---------------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------------

run_branch <- function(config_env, branch_dir, exclude_hiv) {
  message("\n============================================================")
  message("HIV sensitivity branch: ", config_env,
          " (exclude_hiv=", exclude_hiv, ") -> ", branch_dir)
  message("============================================================")

  # ---- run_PoT step ----
  Sys.setenv(R_CONFIG_ACTIVE = config_env)
  conf <- load_config()

  data <- read_transform_data(conf)
  message("[", config_env, "] events after filters: ", nrow(data))

  dir.create(branch_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(branch_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(branch_dir, "tables"), showWarnings = FALSE)

  results <- run_analysis(
    data,
    conf$variables,
    year = conf$base_year,
    return_severity_draws = TRUE,
    output_dir = branch_dir
  )

  saveRDS(results, file.path(branch_dir, "results.rds"))
  saveRDS(conf, file.path(branch_dir, "conf.rds"))
  saveRDS(data, file.path(branch_dir, "data.rds"))

  # ---- run_validation step ----
  # Validation always inherits from `validation` profile; HIV-exclusion is applied
  # via runtime override (R config inherit: only walks one level, so a dedicated
  # config block can't pick up validation's overrides via inherit chaining).
  validation_dir <- file.path(branch_dir, "validation")
  dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(validation_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(validation_dir, "tables"), showWarnings = FALSE)

  Sys.setenv(R_CONFIG_ACTIVE = "validation")
  conf_val <- load_config()
  conf_val$data_filters$exclude_hiv <- exclude_hiv

  data_val <- read_transform_data(conf_val)
  message("[validation, exclude_hiv=", exclude_hiv, "] events after filters: ",
          nrow(data_val))

  variable <- conf_val$variables[1]
  validation_results <- validate_model(
    data_val, variable, conf_val,
    output_dir = validation_dir
  )

  saveRDS(validation_results, file.path(validation_dir, "validation_results.rds"))
  saveRDS(conf_val,           file.path(validation_dir, "conf_validation.rds"))
  saveRDS(data_val,           file.path(validation_dir, "data.rds"))
  message("Saved validation_results.rds to ", validation_dir)

  list(
    conf = conf, data = data, results = results, dir = branch_dir,
    validation = list(
      conf = conf_val, data = data_val,
      results = validation_results, dir = validation_dir
    )
  )
}

summarize_traj <- function(traj_path) {
  if (!file.exists(traj_path)) {
    stop("Cumulative trajectory not found: ", traj_path,
         " — confirm cum_traj_reference.rds is being saved by R/burden.R")
  }
  m <- readRDS(traj_path)
  tibble::tibble(
    year   = as.integer(rownames(m)),
    median = apply(m, 1, stats::median),
    lower  = apply(m, 1, stats::quantile, probs = 0.025),
    upper  = apply(m, 1, stats::quantile, probs = 0.975)
  )
}

format_M <- function(x, digits = 1) {
  formatC(x / 1e6, format = "f", digits = digits)
}

extract_tier_rate_summary <- function(results) {
  yr_list <- results$burden$yearly_rate_detailed[[1]]
  list(
    tier1_rate_mean = yr_list$tier1_rate$mean,
    tier1_rate_low  = yr_list$tier1_rate$lower,
    tier1_rate_up   = yr_list$tier1_rate$upper,
    tier2_rate_mean = yr_list$tier2_rate$mean,
    tier2_rate_low  = yr_list$tier2_rate$lower,
    tier2_rate_up   = yr_list$tier2_rate$upper
  )
}

# ---------------------------------------------------------------------------------------------
# Setup output directory
# ---------------------------------------------------------------------------------------------

# setup_script_run() handles both modes: standalone creates
# output/run_hiv_sensitivity_<ts>/; under run_all_scripts() the harness
# override redirects to output/run_all_<ts>/run_hiv_sensitivity/. Either way
# plots/ and tables/ subdirs are pre-created.
Sys.setenv(R_CONFIG_ACTIVE = "default")
conf_session <- load_config()
run_info <- setup_script_run("run_hiv_sensitivity", conf_session)
session_dir <- run_info$output_dir
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

ref_dir  <- file.path(session_dir, "reference")
excl_dir <- file.path(session_dir, "hiv_excluded")

# ---------------------------------------------------------------------------------------------
# Run reference and HIV-excluded branches
# ---------------------------------------------------------------------------------------------

ref  <- run_branch("default",                 ref_dir,  exclude_hiv = FALSE)
excl <- run_branch("sensitivity_hiv_excluded", excl_dir, exclude_hiv = TRUE)

# ---------------------------------------------------------------------------------------------
# Build comparison figure
# ---------------------------------------------------------------------------------------------

ref_traj  <- summarize_traj(file.path(ref_dir,  "cum_traj_reference.rds"))
excl_traj <- summarize_traj(file.path(excl_dir, "cum_traj_reference.rds"))

ref_traj$scenario  <- "Reference"
excl_traj$scenario <- "HIV/AIDS-excluded"

# Floor for log scale: bootstrap sampling produces lower=0 in the first 1-2 forecast
# years, which is uninformative and breaks log scale. Floor at 0.01M = 10k deaths.
log_floor_M <- 0.01

traj_combined <- dplyr::bind_rows(ref_traj, excl_traj) |>
  dplyr::mutate(
    median_M = pmax(median / 1e6, log_floor_M),
    lower_M  = pmax(lower  / 1e6, log_floor_M),
    upper_M  = pmax(upper  / 1e6, log_floor_M),
    scenario = factor(scenario, levels = c("Reference", "HIV/AIDS-excluded"))
  )

scenario_colors <- c(
  "Reference"          = "#0072B2",   # blue
  "HIV/AIDS-excluded"  = "#D55E00"    # orange-red
)

p <- ggplot(traj_combined, aes(x = year, y = median_M,
                               colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower_M, ymax = upper_M), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = scenario_colors, name = NULL) +
  scale_fill_manual(values = scenario_colors,   name = NULL) +
  scale_y_log10(
    labels = scales::label_comma(),
    breaks = c(0.1, 1, 10, 100, 1000),
    minor_breaks = c(seq(0.1, 1, 0.1), seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100)),
    limits = c(log_floor_M, NA)
  ) +
  annotation_logticks(sides = "l", outside = FALSE,
                      short = unit(0.05, "cm"),
                      mid   = unit(0.1, "cm"),
                      long  = unit(0.15, "cm")) +
  labs(
    x = "Year",
    y = "Cumulative deaths (millions, log scale)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    plot.title = element_text(face = "bold")
  )

fig_basename <- file.path(session_dir, "plots", "hiv_exclusion_comparison")
ggsave(paste0(fig_basename, ".png"), p, width = 7.5, height = 4.5,
       dpi = 600, bg = "white")
ggsave(paste0(fig_basename, ".svg"), p, width = 7.5, height = 4.5,
       device = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else NULL,
       bg = "white")

message("Saved forecast comparison figure: ", fig_basename, ".{png,svg}")

# ---------------------------------------------------------------------------------------------
# Validation comparison figure (1950–2025): same predicted curve (training is pre-1950),
# two observed lines (with vs without HIV/AIDS).
# ---------------------------------------------------------------------------------------------

ref_val_forecast  <- ref$validation$results$predicted$forecast
ref_val_observed  <- ref$validation$results$observed
excl_val_observed <- excl$validation$results$observed

val_pred <- ref_val_forecast |>
  dplyr::transmute(
    year,
    median_M = pmax(cum_deaths     / 1e6, log_floor_M),
    lower_M  = pmax(cum_deaths_low / 1e6, log_floor_M),
    upper_M  = pmax(cum_deaths_up  / 1e6, log_floor_M)
  )

obs_combined <- dplyr::bind_rows(
  dplyr::mutate(ref_val_observed,  scenario = "Observed (incl. HIV/AIDS)",
                obs_M = pmax(obs / 1e6, log_floor_M)),
  dplyr::mutate(excl_val_observed, scenario = "Observed (excl. HIV/AIDS)",
                obs_M = pmax(obs / 1e6, log_floor_M))
) |>
  dplyr::mutate(scenario = factor(scenario,
    levels = c("Observed (incl. HIV/AIDS)", "Observed (excl. HIV/AIDS)")))

obs_colors <- c(
  "Observed (incl. HIV/AIDS)" = "#0072B2",
  "Observed (excl. HIV/AIDS)" = "#D55E00"
)

p_val <- ggplot() +
  geom_ribbon(data = val_pred,
              aes(x = year, ymin = lower_M, ymax = upper_M,
                  fill = "Predicted (95% UI)"), alpha = 0.20, colour = NA) +
  geom_line(data = val_pred,
            aes(x = year, y = median_M, colour = "Predicted (median)"),
            linewidth = 0.9) +
  geom_line(data = obs_combined,
            aes(x = year, y = obs_M, colour = scenario),
            linewidth = 0.9) +
  geom_point(data = obs_combined,
             aes(x = year, y = obs_M, colour = scenario),
             size = 2) +
  scale_colour_manual(values = c(
    "Predicted (median)" = "grey20",
    obs_colors
  ), name = NULL,
  breaks = c("Predicted (median)",
             "Observed (incl. HIV/AIDS)", "Observed (excl. HIV/AIDS)")) +
  scale_fill_manual(values = c("Predicted (95% UI)" = "grey60"), name = NULL) +
  scale_y_log10(
    labels = scales::label_comma(),
    breaks = c(0.1, 1, 10, 100, 1000),
    minor_breaks = c(seq(0.1, 1, 0.1), seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100)),
    limits = c(log_floor_M, NA)
  ) +
  annotation_logticks(sides = "l", outside = FALSE,
                      short = unit(0.05, "cm"),
                      mid   = unit(0.1, "cm"),
                      long  = unit(0.15, "cm")) +
  labs(
    x = "Year",
    y = "Cumulative deaths (millions, log scale)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm")
  )

val_basename <- file.path(session_dir, "plots", "hiv_exclusion_validation")
ggsave(paste0(val_basename, ".png"), p_val, width = 7.5, height = 4.5,
       dpi = 600, bg = "white")
ggsave(paste0(val_basename, ".svg"), p_val, width = 7.5, height = 4.5,
       device = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else NULL,
       bg = "white")

message("Saved validation comparison figure: ", val_basename, ".{png,svg}")

# ---------------------------------------------------------------------------------------------
# Headline deltas (final-year cumulative deaths)
# ---------------------------------------------------------------------------------------------

ref_cum_final  <- readRDS(file.path(ref_dir,  "cum_draws_reference.rds"))
excl_cum_final <- readRDS(file.path(excl_dir, "cum_draws_reference.rds"))

ref_med  <- stats::median(ref_cum_final);   ref_lo <- as.numeric(stats::quantile(ref_cum_final,  0.025)); ref_up <- as.numeric(stats::quantile(ref_cum_final,  0.975))
excl_med <- stats::median(excl_cum_final);  excl_lo <- as.numeric(stats::quantile(excl_cum_final, 0.025)); excl_up <- as.numeric(stats::quantile(excl_cum_final, 0.975))

abs_delta_med   <- ref_med - excl_med
abs_delta_low   <- ref_lo  - excl_lo
abs_delta_up    <- ref_up  - excl_up
pct_change_med  <- 100 * (excl_med - ref_med) / ref_med

ref_tier  <- extract_tier_rate_summary(ref$results)
excl_tier <- extract_tier_rate_summary(excl$results)

n_events_ref  <- nrow(ref$data)
n_events_excl <- nrow(excl$data)
n_events_dropped <- n_events_ref - n_events_excl

# ---------------------------------------------------------------------------------------------
# Validation diagnostics: predicted vs observed at end of test period (2025)
# ---------------------------------------------------------------------------------------------

extract_val_endpoint <- function(validation_results) {
  fc  <- validation_results$predicted$forecast
  obs <- validation_results$observed
  end_year <- max(obs$year)
  pred_row <- dplyr::filter(fc,  year == end_year)
  obs_row  <- dplyr::filter(obs, year == end_year)
  list(
    end_year      = end_year,
    pred_median   = pred_row$cum_deaths,
    pred_low      = pred_row$cum_deaths_low,
    pred_up       = pred_row$cum_deaths_up,
    observed      = obs_row$obs,
    obs_in_pred_ui = (obs_row$obs >= pred_row$cum_deaths_low &
                      obs_row$obs <= pred_row$cum_deaths_up)
  )
}

ref_val_endpoint  <- extract_val_endpoint(ref$validation$results)
excl_val_endpoint <- extract_val_endpoint(excl$validation$results)

summary_obj <- list(
  meta = list(
    timestamp        = timestamp,
    n_events_ref     = n_events_ref,
    n_events_excl    = n_events_excl,
    n_events_dropped = n_events_dropped,
    draws            = ref$conf$draws,
    forecast_window  = ref$conf$time$future_years,
    forecast_start   = ref$conf$time$current_year + 1,
    forecast_end     = ref$conf$time$current_year + ref$conf$time$future_years,
    validation       = list(
      train_start = ref$validation$conf$train_year_start,
      train_end   = ref$validation$conf$train_year_end,
      test_start  = ref$validation$conf$test_year_start,
      test_end    = ref$validation$conf$test_year_end
    )
  ),
  cumulative_deaths = list(
    reference   = list(median = ref_med,  lower_95 = ref_lo,  upper_95 = ref_up),
    hiv_excluded = list(median = excl_med, lower_95 = excl_lo, upper_95 = excl_up),
    abs_delta   = list(median = abs_delta_med, lower_95 = abs_delta_low, upper_95 = abs_delta_up),
    pct_change  = pct_change_med
  ),
  tier_rates = list(
    reference   = ref_tier,
    hiv_excluded = excl_tier,
    tier2_rate_pct_change = if (!is.na(ref_tier$tier2_rate_mean) && ref_tier$tier2_rate_mean > 0)
      100 * (excl_tier$tier2_rate_mean - ref_tier$tier2_rate_mean) / ref_tier$tier2_rate_mean
    else NA_real_,
    tier1_rate_pct_change = if (!is.na(ref_tier$tier1_rate_mean) && ref_tier$tier1_rate_mean > 0)
      100 * (excl_tier$tier1_rate_mean - ref_tier$tier1_rate_mean) / ref_tier$tier1_rate_mean
    else NA_real_
  ),
  validation = list(
    end_year     = ref_val_endpoint$end_year,
    predicted    = list(
      median   = ref_val_endpoint$pred_median,
      lower_95 = ref_val_endpoint$pred_low,
      upper_95 = ref_val_endpoint$pred_up
    ),
    observed = list(
      reference    = ref_val_endpoint$observed,
      hiv_excluded = excl_val_endpoint$observed,
      abs_delta    = ref_val_endpoint$observed - excl_val_endpoint$observed,
      pct_delta    = if (excl_val_endpoint$observed > 0)
        100 * (ref_val_endpoint$observed - excl_val_endpoint$observed) / ref_val_endpoint$observed
      else NA_real_
    ),
    obs_in_pred_ui = list(
      reference    = ref_val_endpoint$obs_in_pred_ui,
      hiv_excluded = excl_val_endpoint$obs_in_pred_ui
    )
  )
)

saveRDS(summary_obj, file.path(session_dir, "hiv_exclusion_summary.rds"))

summary_table <- tibble::tibble(
  metric = c(
    "n_events", "draws",
    "cum_deaths_M_median", "cum_deaths_M_low", "cum_deaths_M_up",
    "tier1_rate_per_year_mean", "tier2_rate_per_year_mean",
    paste0("validation_predicted_M_median_at_", ref_val_endpoint$end_year),
    paste0("validation_observed_M_at_",         ref_val_endpoint$end_year)
  ),
  reference = c(
    n_events_ref, ref$conf$draws,
    ref_med  / 1e6, ref_lo  / 1e6, ref_up  / 1e6,
    ref_tier$tier1_rate_mean, ref_tier$tier2_rate_mean,
    ref_val_endpoint$pred_median  / 1e6,
    ref_val_endpoint$observed     / 1e6
  ),
  hiv_excluded = c(
    n_events_excl, ref$conf$draws,
    excl_med / 1e6, excl_lo / 1e6, excl_up / 1e6,
    excl_tier$tier1_rate_mean, excl_tier$tier2_rate_mean,
    excl_val_endpoint$pred_median / 1e6,
    excl_val_endpoint$observed    / 1e6
  )
) |>
  dplyr::mutate(delta = hiv_excluded - reference)

write.csv(summary_table,
          file.path(session_dir, "tables", "hiv_exclusion_summary.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------------------------

cat("\n=========================  HIV-EXCLUSION SENSITIVITY  =========================\n")
cat(sprintf("Events:          reference = %d, HIV-excluded = %d (dropped %d)\n",
            n_events_ref, n_events_excl, n_events_dropped))
cat(sprintf("Forecast window: %d–%d  (%d-year horizon)\n",
            summary_obj$meta$forecast_start, summary_obj$meta$forecast_end,
            summary_obj$meta$forecast_window))
cat(sprintf("Validation:      train %d–%d, test %d–%d\n",
            summary_obj$meta$validation$train_start, summary_obj$meta$validation$train_end,
            summary_obj$meta$validation$test_start,  summary_obj$meta$validation$test_end))
cat(sprintf("Draws:           %d\n", summary_obj$meta$draws))
cat("-------------------------------------------------------------------------------\n")
cat(sprintf("Forecast cum deaths (M) reference     : %s   (95%% UI %s – %s)\n",
            format_M(ref_med),  format_M(ref_lo),  format_M(ref_up)))
cat(sprintf("                        HIV-excluded  : %s   (95%% UI %s – %s)\n",
            format_M(excl_med), format_M(excl_lo), format_M(excl_up)))
cat(sprintf("                        delta (median): %s M  (%.1f%% change)\n",
            format_M(abs_delta_med), pct_change_med))
cat("-------------------------------------------------------------------------------\n")
cat(sprintf("Tier 1 rate (per year)  reference     : %.4f\n", ref_tier$tier1_rate_mean))
cat(sprintf("                        HIV-excluded  : %.4f  (%.1f%% change)\n",
            excl_tier$tier1_rate_mean, summary_obj$tier_rates$tier1_rate_pct_change))
cat(sprintf("Tier 2 rate (per year)  reference     : %.4f\n", ref_tier$tier2_rate_mean))
cat(sprintf("                        HIV-excluded  : %.4f  (%.1f%% change)\n",
            excl_tier$tier2_rate_mean, summary_obj$tier_rates$tier2_rate_pct_change))
cat("-------------------------------------------------------------------------------\n")
cat(sprintf("Validation @ %d (M)     predicted     : %s   (95%% UI %s – %s)\n",
            ref_val_endpoint$end_year,
            format_M(ref_val_endpoint$pred_median),
            format_M(ref_val_endpoint$pred_low),
            format_M(ref_val_endpoint$pred_up)))
cat(sprintf("                        observed (incl. HIV) : %s   in 95%% UI: %s\n",
            format_M(ref_val_endpoint$observed),
            ifelse(ref_val_endpoint$obs_in_pred_ui, "YES", "NO")))
cat(sprintf("                        observed (excl. HIV) : %s   in 95%% UI: %s\n",
            format_M(excl_val_endpoint$observed),
            ifelse(excl_val_endpoint$obs_in_pred_ui, "YES", "NO")))
cat(sprintf("                        observed delta       : %s M  (%.1f%% of incl-HIV)\n",
            format_M(summary_obj$validation$observed$abs_delta),
            summary_obj$validation$observed$pct_delta))
cat("===============================================================================\n")

if (excl_med > ref_med) {
  warning("HIV-excluded median exceeds reference median — sanity check failed. ",
          "Inspect tier-rate refit (post-1950 catastrophic-tier breakpoint) before reporting.")
} else {
  message("Sanity check passed: HIV-excluded median is lower than reference median.")
}

message("\nSession directory: ", session_dir)
