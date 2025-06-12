#' Summarize wind data by direction and intensity
#'
#' @param wind_data  A data frame containing at least the columns `wind_direction` and `wind_speed`
#' @param wind_percentile either 'mean' or a numeric between 0 and 1 specifying the percentile of wind speed to calculate
#'
#' @return A data frame summarizing wind by direction, with columns:
#'   \describe{
#'     \item{directions}{Wind direction (degrees, with 360 converted to 0)}
#'     \item{n}{Count of observations in each direction}
#'     \item{proportion}{Percentage of total observations in each direction}
#'     \item{speed}{Mean or percentile wind speed for each direction}
#'   }
#'
#' @examples
#'
#' wind_data <- data.frame(
#'   wind_direction = sample(c(90, 180, 270, 360), size = 25, replace = TRUE),
#'   wind_speed = runif(25, 0, 25)
#' )
#'
#' # Summarize using the 95th percentile for wind speed
#' summarize_wind_data(wind_data, 0.95)
#'
#' # Summarize using the mean wind speed
#' summarize_wind_data(wind_data, "mean")
#'
#' @export
summarize_wind_data <- function(wind_data, wind_percentile) {
  # ensure wind from the north is labelled as 0 not 360 and filter out calm conditions
  wind_data <- wind_data %>%
    dplyr::mutate(wind_direction = ifelse(.data$wind_direction == 360, 0, .data$wind_direction)) %>%
    dplyr::filter(.data$wind_speed != 0)

  # user determines how the wind will be summarized - by mean or percentile
  if (wind_percentile == "mean") {
    wind_data_summary <- wind_data %>%
      dplyr::mutate(directions = .data$wind_direction) %>%
      dplyr::group_by(.data$directions) %>%
      dplyr::summarize(
        n = dplyr::n(),
        proportion = 100 * .data$n / nrow(wind_data),
        speed = base::mean(.data$wind_speed, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    wind_data_summary <- wind_data %>%
      dplyr::mutate(directions = .data$wind_direction) %>%
      dplyr::group_by(.data$directions) %>%
      dplyr::summarize(
        n = dplyr::n(),
        proportion = 100 * .data$n / nrow(wind_data),
        speed = stats::quantile(.data$wind_speed, wind_percentile, na.rm = TRUE),
        .groups = "drop"
      )
  }

  return(wind_data_summary)
}
