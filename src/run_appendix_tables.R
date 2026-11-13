####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Generates appendix tables for the manuscript:
##              1. Temporal model comparison (Linear vs Piecewise vs Piecewise Full)
##                 for moderate (100K-3M) and catastrophic (>3M) tiers, with AIC/BIC
##                 and full coefficient estimates, standard errors, and p-values.
##              2. Breakpoint sensitivity for the moderate tier at 1940, 1950, 1960.
##
## Requires:    - Configuration file ("config.yml") with analysis parameters
##
## Outputs:     - tables/model_comparison_table.md (appendix table)
##              - tables/breakpoint_sensitivity_table.md (appendix table)
##
####################################################################################################

source("R/load_all.R")

Sys.setenv(R_CONFIG_ACTIVE = "default")
conf <- load_config()

data <- read_transform_data(conf)

run_info <- setup_script_run("run_appendix_tables", conf)
output_dir <- run_info$output_dir
tables_dir <- file.path(output_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# ===================================================================================
# Helpers
# ===================================================================================

format_p <- function(p) {
  if (is.na(p)) return("—")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

format_coef <- function(est, se, p) {
  if (is.na(est)) return("—")
  sprintf("%.3f (%.3f, %s)", est, se, format_p(p))
}

# Extract coefficient row for a given term from a fitted glm/glm.nb
get_term <- function(model, term) {
  s <- tryCatch(summary(model), error = function(e) NULL)
  if (is.null(s)) return(c(NA, NA, NA))
  co <- s$coefficients
  if (!(term %in% rownames(co))) return(c(NA, NA, NA))
  c(co[term, "Estimate"], co[term, "Std. Error"], co[term, 4])
}

# Build model comparison row for a single fitted model
model_row <- function(name, model, terms) {
  if (is.null(model) || inherits(model, "try-error")) {
    return(NULL)
  }
  aic <- AIC(model)
  bic <- BIC(model)
  ll <- as.numeric(logLik(model))
  term_vals <- lapply(terms, function(t) get_term(model, t))
  tibble::tibble(
    Model = name,
    AIC = round(aic, 1),
    BIC = round(bic, 1),
    LogLik = round(ll, 1),
    !!terms[1] := format_coef(term_vals[[1]][1], term_vals[[1]][2], term_vals[[1]][3]),
    !!terms[2] := format_coef(term_vals[[2]][1], term_vals[[2]][2], term_vals[[2]][3]),
    !!terms[3] := format_coef(term_vals[[3]][1], term_vals[[3]][2], term_vals[[3]][3])
  )
}

# Build tier-specific windowed Poisson GLM inputs (same logic as compare_temporal_models)
build_window_counts <- function(data, variable, lower, upper = NULL,
                                 base_year, window_size, current_year) {
  filtered <- data %>%
    dplyr::filter(start_year >= base_year) %>%
    dplyr::arrange(start_year)
  min_year <- min(filtered$start_year)
  max_year <- max(max(filtered$start_year), current_year)
  window_ends <- seq(min_year + window_size - 1, max_year, by = 1)
  window_starts <- window_ends - window_size + 1
  wc <- tibble::tibble(
    window_start = window_starts,
    window_end = window_ends,
    exceedances = 0
  )
  for (i in seq_len(nrow(wc))) {
    events <- filtered %>%
      dplyr::filter(start_year >= wc$window_start[i], start_year <= wc$window_end[i])
    if (!is.null(upper)) {
      wc$exceedances[i] <- sum(events[[variable]] > lower & events[[variable]] <= upper)
    } else {
      wc$exceedances[i] <- sum(events[[variable]] > lower)
    }
  }
  wc %>%
    dplyr::mutate(
      year_centered = window_end - min(window_end),
      year_scaled = year_centered / max(year_centered)
    )
}

# Fit Linear / Piecewise / Piecewise Full for a given breakpoint on window counts
fit_models_at_breakpoint <- function(wc, breakpoint_year) {
  wc$post_breakpoint <- as.numeric(wc$window_end >= breakpoint_year)
  wc$year_scaled_x_post <- wc$year_scaled * wc$post_breakpoint

  linear <- glm(exceedances ~ year_scaled,
                family = poisson(link = "log"), data = wc)
  piecewise <- glm(exceedances ~ year_scaled + post_breakpoint,
                   family = poisson(link = "log"), data = wc)
  piecewise_full <- glm(
    exceedances ~ year_scaled + post_breakpoint + year_scaled_x_post,
    family = poisson(link = "log"), data = wc
  )
  list(linear = linear, piecewise = piecewise, piecewise_full = piecewise_full)
}

# Tibble -> GitHub-flavored markdown table
tibble_to_md <- function(df) {
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

# ===================================================================================
# Config extraction
# ===================================================================================

variable <- conf$variables[1]
lower_cutoff <- conf$thresholds$lower_cutoff_deaths_scaled
tail_threshold <- conf$thresholds$deaths_scaled
base_year <- conf$base_year
current_year <- conf$time$current_year
window_size <- conf$window$default_size
breakpoint_year_default <- conf$time$breakpoint_year

cat("=== Appendix Tables Configuration ===\n")
cat("Variable:", variable, "\n")
cat("Moderate tier: (", lower_cutoff, ",", tail_threshold, "]\n")
cat("Catastrophic tier: (", tail_threshold, ", Inf)\n")
cat("Base year:", base_year, "| Current year:", current_year, "\n")
cat("Window size:", window_size, "| Default breakpoint:", breakpoint_year_default, "\n\n")

# ===================================================================================
# TABLE 1: Temporal model comparison for each tier
# ===================================================================================

cat("=== Table 1: Temporal Model Comparison ===\n\n")

terms <- c("year_scaled", "post_breakpoint", "year_scaled_x_post")
term_labels <- c("β (year_scaled)", "β (post-1950)", "β (year × post)")

build_tier_comparison <- function(tier_name, lower, upper) {
  wc <- build_window_counts(data, variable, lower, upper,
                             base_year, window_size, current_year)
  models <- fit_models_at_breakpoint(wc, breakpoint_year_default)
  rows <- dplyr::bind_rows(
    model_row("Linear", models$linear, terms),
    model_row("Piecewise-Level", models$piecewise, terms),
    model_row("Piecewise-Full", models$piecewise_full, terms)
  )
  # Rename coefficient columns with pretty labels
  names(rows)[5:7] <- term_labels
  rows$Tier <- tier_name
  rows[, c("Tier", setdiff(names(rows), "Tier"))]
}

tier1_comp <- build_tier_comparison(
  "Moderate (100K-3M)", lower_cutoff, tail_threshold
)
tier2_comp <- build_tier_comparison(
  "Catastrophic (>3M)", tail_threshold, NULL
)

model_comparison <- dplyr::bind_rows(tier1_comp, tier2_comp)
print(model_comparison)

md_model_comp <- paste0(
  "## Table S: Temporal Model Comparison\n\n",
  "Poisson GLM fits on sliding-window exceedance counts (window = ",
  window_size, " years, base year = ", base_year,
  "). Breakpoint fixed at ", breakpoint_year_default,
  ". Coefficients reported as estimate (SE, p-value).\n\n",
  tibble_to_md(model_comparison),
  "\n\n",
  "_Linear: rate ~ year. Piecewise-Level: rate ~ year + I(year >= ",
  breakpoint_year_default, "). ",
  "Piecewise-Full: rate ~ year + I(year >= ", breakpoint_year_default,
  ") + year:I(year >= ", breakpoint_year_default, "). ",
  "Lower AIC/BIC indicates better fit._\n"
)
writeLines(md_model_comp, file.path(tables_dir, "model_comparison_table.md"))
cat("\nSaved:", file.path(tables_dir, "model_comparison_table.md"), "\n\n")

# ===================================================================================
# TABLE 2: Breakpoint sensitivity for moderate tier
# ===================================================================================

cat("=== Table 2: Breakpoint Sensitivity (Moderate Tier) ===\n\n")

breakpoints <- c(1940, 1950, 1960)
wc_mod <- build_window_counts(data, variable, lower_cutoff, tail_threshold,
                              base_year, window_size, current_year)

# Linear baseline (no breakpoint) for ΔAIC reference
linear_baseline <- glm(exceedances ~ year_scaled,
                        family = poisson(link = "log"), data = wc_mod)
baseline_aic <- AIC(linear_baseline)
baseline_bic <- BIC(linear_baseline)

sens_rows <- list()
for (bp in breakpoints) {
  models <- fit_models_at_breakpoint(wc_mod, bp)
  for (mtype in c("piecewise", "piecewise_full")) {
    m <- models[[mtype]]
    year_coef <- get_term(m, "year_scaled")
    post_coef <- get_term(m, "post_breakpoint")
    interact_coef <- get_term(m, "year_scaled_x_post")
    sens_rows[[length(sens_rows) + 1]] <- tibble::tibble(
      Breakpoint = bp,
      Model = ifelse(mtype == "piecewise", "Piecewise-Level", "Piecewise-Full"),
      AIC = round(AIC(m), 1),
      BIC = round(BIC(m), 1),
      `ΔAIC vs Linear` = round(AIC(m) - baseline_aic, 1),
      `β (year)` = format_coef(year_coef[1], year_coef[2], year_coef[3]),
      `β (post)` = format_coef(post_coef[1], post_coef[2], post_coef[3]),
      `β (year × post)` = format_coef(interact_coef[1], interact_coef[2], interact_coef[3])
    )
  }
}
breakpoint_sens <- dplyr::bind_rows(sens_rows)
print(breakpoint_sens)

md_bp_sens <- paste0(
  "## Table S: Breakpoint Sensitivity — Moderate Tier\n\n",
  "Piecewise Poisson GLM fits for the moderate tier (",
  format(lower_cutoff, big.mark = ","), "–",
  format(tail_threshold, big.mark = ","),
  " 2025-equivalent deaths) using alternative breakpoint years. ",
  "Baseline linear model AIC = ", round(baseline_aic, 1),
  ", BIC = ", round(baseline_bic, 1), ".\n\n",
  tibble_to_md(breakpoint_sens),
  "\n\n",
  "_Negative ΔAIC indicates the piecewise specification improves fit over the linear baseline. ",
  "Results should be robust to breakpoint choice within the 1940–1960 window._\n"
)
writeLines(md_bp_sens, file.path(tables_dir, "breakpoint_sensitivity_table.md"))
cat("\nSaved:", file.path(tables_dir, "breakpoint_sensitivity_table.md"), "\n\n")

# ===================================================================================
# Save raw tibbles for reference (not consumed by any downstream script)
# ===================================================================================

saveRDS(list(
  model_comparison = model_comparison,
  breakpoint_sensitivity = breakpoint_sens,
  baseline_linear_aic = baseline_aic,
  baseline_linear_bic = baseline_bic
), file.path(output_dir, "appendix_tables.rds"))

cat("=== Appendix tables complete ===\n")
