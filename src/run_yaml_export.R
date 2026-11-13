####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Export pipeline results to _variables.yml for Quarto manuscript integration.
##              Loads saved .rds files from run_pot and run_validation, extracts all variables,
##              and writes a formatted YAML file.
##
## Requires:    - results.rds and conf.rds from run_PoT.R
##              - validation_results.rds and conf_validation.rds from run_validation.R
##
## Outputs:     - _variables.yml with all model results formatted for Quarto
##
####################################################################################################

# Load project functions & packages explicitly
source("R/load_all.R")

# ===================================================================================
# Setup environment and configuration
# ===================================================================================
Sys.setenv(R_CONFIG_ACTIVE = "default")
conf <- load_config()

# ===================================================================================
# Setup output management
# ===================================================================================
run_info <- setup_script_run("run_yaml_export", conf)
output_dir <- run_info$output_dir

# ===================================================================================
# Locate .rds files from upstream scripts
# ===================================================================================
# When run via run_all, sibling directories are under the same parent
parent_dir <- dirname(output_dir)

# Try sibling directories first (run_all mode), then fall back to latest output dirs
find_rds <- function(sibling_name, filename) {
  # Try sibling directory (run_all mode)
  sibling_path <- file.path(parent_dir, sibling_name, filename)
  if (file.exists(sibling_path)) return(sibling_path)

  # Fall back: search output/ for most recent directory containing the file
  output_base <- file.path(getwd(), "output")
  if (!dir.exists(output_base)) stop("Cannot find ", filename, " - no output directory exists")

  # Look in run_all_* directories first
  candidates <- list.files(output_base, pattern = "^run_all_", full.names = TRUE)
  for (candidate in rev(sort(candidates))) {
    path <- file.path(candidate, sibling_name, filename)
    if (file.exists(path)) return(path)
  }

  # Look in script-specific timestamped directories
  script_dirs <- list.files(output_base, pattern = paste0("^", sibling_name, "_"),
                            full.names = TRUE)
  for (d in rev(sort(script_dirs))) {
    path <- file.path(d, filename)
    if (file.exists(path)) return(path)
  }

  stop("Cannot find ", filename, " from ", sibling_name,
       ". Run run_PoT.R and run_validation.R first.")
}

results_path <- find_rds("run_pot", "results.rds")
conf_main_path <- find_rds("run_pot", "conf.rds")
validation_path <- find_rds("run_validation", "validation_results.rds")
conf_validation_path <- find_rds("run_validation", "conf_validation.rds")

message("Loading results from: ", results_path)
message("Loading validation from: ", validation_path)

# ===================================================================================
# Load results
# ===================================================================================
results <- readRDS(results_path)
conf_main <- readRDS(conf_main_path)
validation_results <- readRDS(validation_path)
conf_validation <- readRDS(conf_validation_path)

# Take conf_main$historical from config.yml when the saved configuration lacks it.
if (is.null(conf_main$historical) && !is.null(conf$historical)) {
  conf_main$historical <- conf$historical
} else if (!is.null(conf$historical)) {
  for (key in names(conf$historical)) {
    if (is.null(conf_main$historical[[key]])) {
      conf_main$historical[[key]] <- conf$historical[[key]]
    }
  }
}

# ===================================================================================
# Load sensitivity summary (optional — from run_sens)
# ===================================================================================
sens_summary <- tryCatch({
  sens_path <- find_rds("run_sens", "sensitivity_summary.rds")
  message("Loading sensitivity summary from: ", sens_path)
  readRDS(sens_path)
}, error = function(e) {
  message("Note: No sensitivity summary found (", e$message, "). ",
          "Sensitivity section will be omitted from _variables.yml.")
  NULL
})

# ===================================================================================
# Extract and write YAML
# ===================================================================================
yaml_data <- extract_variables_from_results(results, validation_results,
                                             conf_main, conf_validation,
                                             output_dir = dirname(results_path))

# Merge sensitivity summary into yaml_data
if (!is.null(sens_summary)) {
  yaml_data$sensitivity <- extract_sensitivity_summary(sens_summary)
}

# ===================================================================================
# Load HIV-exclusion sensitivity summary (optional — from run_hiv_sensitivity)
# ===================================================================================
hiv_summary <- tryCatch({
  output_base <- file.path(getwd(), "output")
  # Search both layouts: standalone (output/run_hiv_sensitivity_<ts>/) and
  # under-run_all (output/run_all_<ts>/run_hiv_sensitivity/). Pick the
  # newest hiv_exclusion_summary.rds by mtime so the most recent HIV
  # sensitivity run is picked up regardless of where it was produced.
  standalone <- Sys.glob(file.path(
    output_base, "run_hiv_sensitivity_*", "hiv_exclusion_summary.rds"
  ))
  under_run_all <- Sys.glob(file.path(
    output_base, "run_all_*", "run_hiv_sensitivity", "hiv_exclusion_summary.rds"
  ))
  candidates <- c(standalone, under_run_all)
  if (length(candidates) == 0) {
    stop("no hiv_exclusion_summary.rds found in standalone or run_all sessions")
  }
  hit <- candidates[which.max(file.info(candidates)$mtime)]
  message("Loading HIV-exclusion summary from: ", hit)
  readRDS(hit)
}, error = function(e) {
  message("Note: No HIV-exclusion summary found (", e$message, "). ",
          "sensitivity.hiv_exclusion block will be omitted from _variables.yml.")
  NULL
})

if (!is.null(hiv_summary)) {
  if (is.null(yaml_data$sensitivity)) yaml_data$sensitivity <- list()
  yaml_data$sensitivity$hiv_exclusion <- extract_hiv_exclusion_summary(hiv_summary)
}

output_file <- file.path(output_dir, "_variables.yml")
write_variables_yaml(yaml_data, output_file)

message("YAML export complete: ", output_file)
