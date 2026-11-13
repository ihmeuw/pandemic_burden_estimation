####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Main orchestration script for Peaks-over-Threshold (PoT) mixture model analysis.
##              This script coordinates the entire pandemic estimation workflow by:
##              1. Setting up the environment and loading all required functions
##              2. Reading and transforming pandemic death data
##              3. Running the main analysis with the mixture model approach
##              4. Generating summary tables and visualizations of results
##
## Requires:    - Configuration file ("config.yml") with analysis parameters
##
## Outputs:     - Statistical summary of pandemic severity estimates
##              - Expected annual burden projections
##              - Scenario comparison tables
##              - Visualization plots
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
# Data preparation
# ===================================================================================
data <- read_transform_data(conf)

# ===================================================================================
# Setup output management
# ===================================================================================
# Setup script-specific timestamped output directory
run_info <- setup_script_run("run_pot", conf)
output_dir <- run_info$output_dir

# ===================================================================================
# Main analysis execution
# ===================================================================================
# Run analysis for all configured variables
results <- run_analysis(data, conf$variables,
  year = conf$base_year,
  return_severity_draws = TRUE,
  output_dir = output_dir
)

# ===================================================================================
# Save results for downstream scripts (e.g., YAML export)
# ===================================================================================
# Save immediately after analysis, before tables/plots that may fail
saveRDS(results, file.path(output_dir, "results.rds"))
saveRDS(conf, file.path(output_dir, "conf.rds"))
message("Saved results.rds and conf.rds to ", output_dir)

# ===================================================================================
# Results generation and output
# ===================================================================================
# Generate summary tables
generate_summary_table(results, digits = 3, output_dir = output_dir, conf = conf)
generate_parameter_summary(results, output_dir = output_dir, conf = conf)
generate_scenario_table(results, output_dir = output_dir, conf = conf)
generate_gap_table(results, output_dir = output_dir, conf = conf)
