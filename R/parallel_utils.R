#' Initialize parallel processing
#'
#' Sets up parallel backend for computation with automatic detection of optimal worker count.
#' @param workers Number of parallel workers to use (default: NULL, which uses available cores minus 1)
#' @param strategy Character, parallel strategy to use: 'future' or 'none'
#' @return Invisible logical indicating if parallel backend was successfully initialized
#' @export
setup_parallel <- function(workers = NULL, strategy = "future") {
  strategy <- match.arg(strategy, c("future", "none"))

  if (strategy == "none") {
    message("Using sequential processing (no parallelization)")
    return(invisible(FALSE))
  }

  if (is.null(workers)) {
    # Default to number of available cores minus 1 (leave one for system)
    workers <- max(1, parallel::detectCores() - 1)
  }

  message(paste0("Setting up parallel processing with ", workers, " workers using ", strategy, " strategy"))

  if (!requireNamespace("future", quietly = TRUE) ||
    !requireNamespace("future.apply", quietly = TRUE)) {
    warning("Required packages 'future' and 'future.apply' not available. Using sequential processing.")
    return(invisible(FALSE))
  }

  future::plan(future::multisession, workers = workers)
  invisible(TRUE)
}

#' Stop parallel processing
#'
#' Shuts down parallel backend and cleans up resources.
#' @return Invisible NULL
#' @export
stop_parallel <- function() {
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }

  return(invisible(NULL))
}

#' Run parallel sensitivity analysis
#'
#' Executes a sensitivity analysis across multiple parameter values in parallel
#' @param data Analysis dataset
#' @param var Variable to analyze ("deaths" or "deaths_scaled")
#' @param analysis_type Type of sensitivity analysis ("year", "threshold", "cutoff")
#' @param values Vector of values to analyze over (years, thresholds, or cutoff values)
#' @param conf Configuration list
#' @param strategy Parallel strategy to use: "future" or "none"
#' @param workers Number of workers for parallel processing (default: NULL uses available cores - 1)
#' @return A combined result object from the parallel runs
#' @export
run_parallel_sensitivity <- function(
    data,
    var,
    analysis_type = "year",
    values,
    conf,
    strategy = "future",
    workers = NULL,
    output_dir = NULL) {
  # ---------------------------------------------------------------------------
  # Internal helpers (non-exported) to reduce repeated code blocks
  # ---------------------------------------------------------------------------
  get_analysis_fn <- function(analysis_type, data, var, conf) {
    switch(analysis_type,
      year = function(x) run_analysis(data, var, year = x, conf = conf),
      threshold = function(x) run_analysis(data, var, threshold = x, conf = conf),
      cutoff = function(x) {
        # Only vary lower_cutoff; keep GPD threshold (u) and tier boundary fixed
        run_analysis(
          data, var,
          lower_cutoff = x,
          conf = conf
        )
      },
      stop("Unsupported analysis_type: ", analysis_type)
    )
  }

  process_result <- function(result, analysis_type) {
    param_name <- switch(analysis_type,
      year = "year",
      threshold = "threshold",
      cutoff = "cutoff",
      "value"
    )

    base_tibble <- tibble::tibble(
      !!param_name := switch(analysis_type,
        year = result$specs$Year,
        threshold = result$specs$Threshold,
        cutoff = result$specs$Cutoff,
        NA_real_
      ),
      shadow = ifelse(!is.null(result$results$shadow), result$results$shadow, NA_real_),
      boot = ifelse(!is.null(result$results$boot), result$results$boot, NA_real_),
      boot_low = ifelse(!is.null(result$results$boot_low), result$results$boot_low, NA_real_),
      boot_upp = ifelse(!is.null(result$results$boot_upp), result$results$boot_upp, NA_real_),
      Xi = ifelse(!is.null(result$params$xi), result$params$xi, NA_real_),
      Beta = ifelse(!is.null(result$params$beta), result$params$beta, NA_real_)
    )

    if (!is.null(result$burden)) {
      burden_fields <- tibble::tibble(
        yearly_deaths = result$burden$yearly_deaths,
        yearly_deaths_low = result$burden$yearly_deaths_low,
        yearly_deaths_up = result$burden$yearly_deaths_up,
        cum_deaths = result$burden$cum_deaths,
        cum_deaths_low = result$burden$cum_deaths_low,
        cum_deaths_up = result$burden$cum_deaths_up,
        exp_annual_burden = result$burden$yearly_deaths
      )
      base_tibble <- dplyr::bind_cols(base_tibble, burden_fields)
    }

    base_tibble
  }

  # Wrapper used by all strategies
  analysis_fn <- get_analysis_fn(analysis_type, data, var, conf)
  run_single <- function(value) {
    res <- analysis_fn(value)
    process_result(res, analysis_type)
  }
  is_parallel <- setup_parallel(workers, strategy)


  if (strategy == "future" && is_parallel) {
    raw_results <- future.apply::future_lapply(
      values,
      FUN = function(value) run_single(value),
      future.seed = TRUE,
      future.packages = c("dplyr", "tibble")
    )
  } else {
    # Fall back to sequential processing
    raw_results <- lapply(values, run_single)
  }

  results <- dplyr::bind_rows(raw_results) %>%
    # Convert deaths/burden columns from counts to millions (incl. the GPD scale parameter Beta, which is in death units)
    dplyr::mutate(dplyr::across(
      c(
        shadow, boot, boot_low, boot_upp,
        yearly_deaths, yearly_deaths_low, yearly_deaths_up,
        cum_deaths, cum_deaths_low, cum_deaths_up,
        exp_annual_burden, Beta
      ),
      ~ ifelse(!is.na(.x), .x / 1e6, .x)
    ))

  if (is_parallel) {
    stop_parallel()
  }

  plot <- plot_cumulative_deaths(results, var, analysis_type, config = conf, output_dir = output_dir, conf = conf)

  list(results = results, plot = plot)
}
