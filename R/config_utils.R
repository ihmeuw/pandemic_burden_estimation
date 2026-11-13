# Required Packages --------------------------------------------------------------
library(tidyverse)
library(evir)
library(gridExtra)
library(scales)
library(purrr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)
library(config)
library(data.table)

#' Load Configuration Settings
#'
#' @param config_file Path to config file
#' @param environment Configuration environment to use
#' @return List containing configuration settings
#' @export
load_config <- function(
    config_file = "config.yml",
    environment = Sys.getenv("R_CONFIG_ACTIVE", "default")) {
  conf <- config::get(file = config_file, config = environment)

  # Apply population-scaled threshold calculations for validation configuration
  if (environment == "validation") {
    # Load default config for reference values
    default_conf <- config::get(file = config_file, config = "default")

    conf <- calculate_population_scaled_thresholds(default_conf, conf)
  }

  return(conf)
}

#' Read and Transform Analysis Data
#'
#' Reads the data file and applies necessary transformations according to configuration.
#' Filters to events flagged `included == 1` and, when configured, drops HIV/AIDS events.
#'
#' @param conf List containing configuration settings
#' @return A tibble with transformed data
#' @export
read_transform_data <- function(conf) {
  data <- read_csv(conf$file_path, show_col_types = FALSE)

  # Apply data filters if columns exist and filters are configured
 if (!is.null(conf$data_filters)) {
    if ("included" %in% names(data) && isTRUE(conf$data_filters$included_only)) {
      data <- data %>% filter(included == 1)
    }

    # Exclude HIV/AIDS events (case-insensitive whole-word match) if filter is
    # enabled. Word boundaries prevent false-positives on substrings like
    # "Maids" or "archival" — relevant if source-citation text leaks into
    # the name column.
    if ("name" %in% names(data) && isTRUE(conf$data_filters$exclude_hiv)) {
      n_before <- nrow(data)
      data <- data %>% filter(!grepl("\\b(HIV|AIDS)\\b", name, ignore.case = TRUE))
      message("[exclude_hiv] dropped ", n_before - nrow(data),
              " HIV/AIDS event(s); ", nrow(data), " events remain")
    }
  }

  data %>%
    mutate(deaths = deaths * conf$deaths_unit) %>%
    mutate(
      length = end_year - start_year + 1,
      deaths_scaled = deaths * (conf$pop_reference / pop),
      deaths_scaled_year = (deaths * (conf$pop_reference / pop)) / length
    )
}

#' Get IHME Population Forecast Data
#'
#' Retrieves and processes IHME's population forecast global values.
#'
#' @return A tibble containing processed population forecast data
#' @export
get_pop_forecast <- function() {
  readr::read_csv("data/pop_forecast.csv", show_col_types = FALSE) %>%
    dplyr::select(
      year = .data$year_id,
      .data$sex,
      .data$val,
      .data$upper,
      .data$lower
    ) %>%
    dplyr::group_by(.data$year) %>%
    dplyr::summarize(
      val = sum(.data$val),
      upper = sum(.data$upper),
      lower = sum(.data$lower),
      .groups = "drop"
    )
}


#' Get IHME Population Forecast Draws
#'
#' Retrieves and processes IHME's population forecast global values with draws.
#' These draws represent uncertainty in future population projections and are
#' used for propagating uncertainty in burden calculations.
#'
#' @return A tibble containing processed population forecast draws, with years as rows
#'   and draws as columns (draw_1, draw_2, etc.)
#' @export
get_pop_draws <- function() {
  readr::read_csv("data/pop_draws.csv", show_col_types = FALSE) %>%
    dplyr::select(
      year = year_id,
      draw,
      value
    ) |>
    pivot_wider(names_from = draw, values_from = value, names_prefix = "draw_")
}

#' Calculate Population-Scaled Thresholds
#'
#' Automatically calculates threshold values for configurations that need to scale
#' thresholds based on population reference ratios. This is particularly useful
#' for validation configurations that should inherit default thresholds but scale
#' them according to their different population reference values.
#'
#' @param default_config Configuration object containing default threshold values
#' @param target_config Configuration object that needs scaled thresholds
#' @param threshold_names Character vector of threshold names to scale (default: c("deaths_scaled"))
#' @return Modified target_config with calculated scaled thresholds
#' @export
calculate_population_scaled_thresholds <- function(
    default_config,
    target_config,
    threshold_names = c("deaths_scaled")) {
  population_ratio <- target_config$pop_reference / default_config$pop_reference

  for (threshold_name in threshold_names) {
    if (!is.null(default_config$thresholds[[threshold_name]])) {
      target_config$thresholds[[threshold_name]] <-
        default_config$thresholds[[threshold_name]] * population_ratio
    }
  }

  return(target_config)
}

