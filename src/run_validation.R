####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Orchestration script for pandemic model validation.
##              This script coordinates the model validation workflow by:
##              1. Setting up the environment with validation configuration
##              2. Reading and transforming historical pandemic data
##              3. Splitting data into training and validation periods
##              4. Training the model on historical data and validating on recent data
##              5. Generating visualization comparing predicted vs observed pandemic patterns
##
## Requires:    - Configuration file ("config.yml") with validation parameters
##
## Outputs:     - Model validation results measuring prediction accuracy
##              - Visualization comparing predicted vs. observed pandemic patterns
##              - Model performance metrics
##              - Plots
##
####################################################################################################

# Load project functions & packages explicitly
source("R/load_all.R")

# ===================================================================================
# Setup environment and configuration
# ===================================================================================
Sys.setenv(R_CONFIG_ACTIVE = "validation")
conf <- load_config()

# ===================================================================================
# Data preparation
# ===================================================================================
data <- read_transform_data(conf)

# ===================================================================================
# Model validation execution
# ===================================================================================
# Select variable for validation (first configured variable)
variable <- conf$variables[1]

# ===================================================================================
# Setup output management
# ===================================================================================
# Setup script-specific timestamped output directory
run_info <- setup_script_run("run_validation", conf)
output_dir <- run_info$output_dir

# Run model validation process
validation_results <- validate_model(data, variable, conf, output_dir = output_dir)

# ===================================================================================
# Save results for downstream scripts (e.g., YAML export)
# ===================================================================================
# Save immediately after analysis, before plots that may fail
saveRDS(validation_results, file.path(output_dir, "validation_results.rds"))
saveRDS(conf, file.path(output_dir, "conf_validation.rds"))
message("Saved validation_results.rds and conf_validation.rds to ", output_dir)

# ===================================================================================
# Results generation and output
# ===================================================================================
# Generate validation summary table
generate_validation_table(validation_results, output_dir = output_dir, conf = conf)

# Display and save validation comparison plot
if (!is.null(validation_results$plot)) {
  save_plot_if_enabled(validation_results$plot, "validation_comparison_plot", output_dir, conf, scale = 0.8)
  tryCatch(print(validation_results$plot), error = function(e) message("Note: Could not display plot: ", e$message))
}

# Display and save SFA validation plots
if (!is.null(validation_results$sfa_validation_plot)) {
  save_plot_if_enabled(validation_results$sfa_validation_plot,
                       "sfa_validation_scenarios", output_dir, conf)
  tryCatch(print(validation_results$sfa_validation_plot), error = function(e) message("Note: Could not display SFA validation plot: ", e$message))
}
if (!is.null(validation_results$sfa_diagnostic_plot)) {
  save_plot_if_enabled(validation_results$sfa_diagnostic_plot,
                       "sfa_validation_diagnostic_rates", output_dir, conf)
  tryCatch(print(validation_results$sfa_diagnostic_plot), error = function(e) message("Note: Could not display SFA diagnostic plot: ", e$message))
}
