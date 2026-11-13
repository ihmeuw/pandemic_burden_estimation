#' Get Year Range from Data
#'
#' @param data The input dataset containing start_year column
#' @param conf Configuration list; the year range is sliced with conf$sens$year$start_idx:conf$sens$year$end_idx (config.yml sensitivity: 10 and 50)
#' @return Vector of ordered years within the specified range
#' @export
get_year_range <- function(data, conf) {
  data %>%
    distinct(start_year) %>%
    arrange(start_year) %>%
    slice(conf$sens$year$start_idx:conf$sens$year$end_idx) %>%
    pull(start_year)
}
