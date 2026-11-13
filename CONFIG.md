# Configuration Parameter Reference

This document describes the parameters in `config.yml`. The file itself is commented line by line and is the authoritative reference; all values are set to those used for the paper.

## Data & Input Parameters

`file_path` (string) - Path to input pandemic dataset CSV file  
`pop_reference` (number) - Reference population size for rescaling historical deaths to 2025-equivalent counts  
`pop_min` (number) - Minimum population threshold value  
`base_year` (number) - Starting year for historical analysis  
`variables` (array) - Primary variables to analyze [deaths_scaled]  
`deaths_unit` (number) - Unit multiplier for death counts  
`historical.total_database_events` (number) - Number of events in the full database, quoted in the paper  
`historical.events_ple_qualified` (number) - Number of events qualifying as pandemic-level events, quoted in the paper  
`data_filters.included_only` (boolean) - Restrict to events flagged `included == 1`  
`data_filters.exclude_hiv` (boolean) - Drop HIV/AIDS events (case-insensitive `name` match on `HIV`/`AIDS`); default `false`. Used by the `sensitivity_hiv_excluded` config profile and the `src/run_hiv_sensitivity.R` driver to test the influence of HIV/AIDS classification on burden estimates.  

## Analysis Method Parameters

`method` (string) - Extreme value statistical method [mev, evir]  
`bulk_distribution` (string) - Distribution for bulk component [lognormal, weibull]  
`use_constrained_model` (boolean) - Enable constrained mixture model fitting  

## Threshold Parameters

`thresholds.deaths_scaled` (number) - 2025-equivalent death threshold for tail distribution  
`thresholds.lower_cutoff_deaths_scaled` (number) - Lower cutoff for mixture model (2025-equivalent deaths)  
`thresholds.lower_cutoff` (number) - Analytical floor (100 000 2025-equivalent deaths) read throughout the pipeline; `run_analysis()` overwrites it in memory with `lower_cutoff_<variable>` when that key exists  
`bulk_diagnostics.tail_threshold_fixed` (number) - Fixed tail threshold for diagnostics  
`bulk_diagnostics.lower_thresholds` (array) - Candidate lower cutoffs for testing  

## Bootstrap & Uncertainty Parameters

`draws` (number) - Number of bootstrap draws for uncertainty quantification  
`seed` (number) - Random seed for reproducibility  
`bootstrap.enabled` (boolean) - Enable bootstrap analysis  
`bootstrap.plot` (boolean) - Generate bootstrap diagnostic plots  

## Time Analysis Parameters

`time.future_years` (number) - Years to forecast into future  
`time.use_time_trend` (boolean) - Include time trend in forecasts  
`time.model_type` (string) - Time trend model type [constant, linear, quadratic, nb_linear, nb_quadratic, piecewise, nb_piecewise, auto]  
`time.tier1_model_type`, `time.tier2_model_type` (string) - Rate model for the moderate and catastrophic tiers; take precedence over `time.model_type`  
`time.breakpoint_year` (number) - Year of the intercept shift in piecewise models  
`time.current_year` (number) - Reference year for forecasts  
`window.default_size` (number) - Window size in years for the sliding-window exceedance counts  

## Burden Calculation Parameters

`calc_burden` (boolean) - Enable burden calculation  
`calc_dalys` (boolean) - Enable DALY estimation  
`calc_scenarios` (boolean) - Enable best/worst scenario calculations  
`sfa_bounds` (boolean) - Derive the plausible-best and plausible-worst scenario rates from frontier bounds on the tier rate models  

## Plotting Parameters

`plots.descriptive` (boolean) - Generate descriptive plots  
`plots.tail_mixture` (boolean) - Generate mixture model plots  
`plots.tail_mixture_y_title` (string) - Y-axis title for mixture plots  
`plots.log_scale` (boolean) - Use logarithmic scale  
`plots.forecast` (boolean) - Generate forecast plots  
`plots.time_trend` (boolean) - Plot time trends  
`plots.time_scaled_deaths` (boolean) - Generate time-scaled death plots  
`plots.time_scaled_deaths_labels` (boolean) - Include labels in time plots  
`plots.time_scaled_deaths_min` (number) - Smallest event drawn in the time-scaled deaths plot  
`plots.validation_plot_labels` (boolean) - Include labels in validation plots  
`plots.log_breaks_major` (array) - Major breaks for log-scaled axes  
`plots.point_size` (number) - Size of empirical data points  

## Output Management Parameters

`output.save_plots` (boolean) - Save plots to files  
`output.save_tables` (boolean) - Save tables to files  
`output.plot_formats` (array) - Plot output formats [svg, png, pdf, eps]  

## Parallel Processing Parameters

`parallel.enabled` (boolean) - Enable parallel processing  
`parallel.strategy` (string) - Parallelization strategy [future, none]  
`parallel.workers` (number) - Number of parallel workers  

## Tables Parameters

`tables.projection_years` (number) - Default projection period for summary tables  

## Validation-Specific Parameters

`train_year_start` (number) - Training period start year  
`train_year_end` (number) - Training period end year  
`test_year_start` (number) - Testing period start year  
`test_year_end` (number) - Testing period end year  
`window_size` (number) - Window size for validation  
`observation_years` (array) - Years for observation comparison  

## Sensitivity Analysis Parameters

`sens.year.start_idx` (number) - Start index for year range sensitivity  
`sens.year.end_idx` (number) - End index for year range sensitivity  
`sens.year.cutoff` (number) - Year cutoff value  
`sens.year.base_year` (number) - Base year for analysis  
`sens.threshold.base_thresholds.deaths_scaled` (number) - Baseline threshold marked in the threshold-sensitivity plot  

## Configuration Inheritance

Configuration sections can inherit from other sections using the `inherit` parameter:
- `inherit: default` - Inherit all parameters from the default section and override specific values
- Available in the validation, sensitivity and sensitivity_hiv_excluded configurations
