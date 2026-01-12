#' Retrieve and clean wind data for WEMo
#'
#' This function is a wrapper around a workflow to find and download wind data
#' using the \link[worldmet:import_ghcn_stations]{import_ghcn_stations()} and
#' \link[worldmet:import_ghcn_hourly]{import_ghcn_hourly()} functions from the
#' **worldmet** package. It finds up to 5 nearby NOAA GHCN weather stations based
#' on a geographic point and downloads hourly wind data.
#'
#' @details This function utilizes data from **Global Historical Climatology
#'   Network (GHCN)** dataset. GHCN aggregates hourly meteorological
#'   observations from numerous fixed, land-based stations. maintained by NOAA,
#'   the U.S. Air Force, and many other meteorological agencies (Met Services)
#'   around the world.
#'
#' @param site_point A spatial point object (either `sf` or coercible to `sf`)
#'   indicating the target location.
#' @param years A vector of years (e.g., `2002:2024`) for which to download wind
#'   data.
#' @param which_station Either:
#'   - `"ask"`: interactively choose from 5 closest stations;
#'   - an integer (1-5): pick the nth closest station;
#'   - a string station code in the form "USAF-WBAN".
#'
#' @return A data frame with cleaned wind direction and wind speed, along with
#'   station code, timestamp, and date components (year, month, day).
#'
#' @examples
#' # Create a POINT geometry from coordinates and cast explicitly
#' site_point <- sf::st_sf(
#'   geometry = sf::st_sfc(sf::st_point(c(-76.67587, 34.71413))),
#'   crs = 4326
#'   )
#'
#' # 1. Prompt user to select a station interactively
#' \dontrun{
#' wind_data_ask <- get_wind_data(
#'   site_point,
#'   years = 2022:2023,
#'   which_station = "ask"
#'  )
#' }
#'
#' # 2. Automatically use the 2nd closest station
#' \dontrun{
#' wind_data_index <- get_wind_data(
#'   site_point,
#'   years = 2022:2023,
#'   which_station = 2
#'  )
#' }
#'
#' # 3. Manually specify a station code
#' \dontrun{
#' # this way you don't need a site_point
#' wind_data_manual <- get_wind_data(
#'   site_point = NULL,
#'   years = 2022:2023,
#'   which_station = "USW00093765"
#'  )
#' }
#'
#' @references For more information on the GHCN and to view an interactive
#'   map of stations, see
#'   https://www.ncei.noaa.gov/products/global-historical-climatology-network-hourly
#'
#' @export
get_wind_data <- function(site_point, years, which_station = 'ask') {
  if(missing(site_point)&(which_station %in% c("ask", 1:5))){
    stop('"site_point" must be supplied when station code isn\'t provided')
  }
  # If user asks to choose or uses index 1-5, begin by locating closest stations
  if (which_station %in% c("ask", 1:5)) {
    # Ensure site_point is an sf object
    if (!inherits(site_point, "sf")) {
      tryCatch({
        site_point <- sf::st_as_sf(site_point)
      }, error = function(e) {
        stop('"site_point" needs to be sf object or coerable to sf: ', e$message)
      })
    }

    # Extract latitude and longitude
    LAT <- sf::st_coordinates(site_point)[[1,2]]
    LON <- sf::st_coordinates(site_point)[[1,1]]

    station <- worldmet::import_ghcn_stations(lat = LAT, lng = LON, crs = sf::st_crs(site_point), n_max = 5, return = 'table')

    if(which_station == "ask"){
      inventory_raw <- worldmet::import_ghcn_inventory(database = 'hourly')

      inventory <- inventory_raw %>%
        dplyr::filter(.data$id %in% station$id) %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarize(
          start_year = min(.data$year, na.rm = TRUE),
          end_year = max(.data$year, na.rm = TRUE),
          # Use the helper to get the ranges or an empty string
          gap_string = summarize_years(setdiff(min(.data$year):max(.data$year), .data$year)),
          .groups = "drop"
        ) %>%
        # Join back to get names/dist if needed, then:
        dplyr::mutate(
          gap_label = ifelse(.data$gap_string == "",
                             "",
                             paste0(" [missing years: ", .data$gap_string, "]"))
        )

      inventory <- dplyr::right_join(station, inventory, by = dplyr::join_by(.data$id)) %>%
        dplyr::arrange(.data$distance)

      worldmet::import_ghcn_stations(lat = LAT, lng = LON, crs = sf::st_crs(site_point), n_max = 5, return = 'map')
    }

    if(which_station == "ask"){

      # Build the display options
      options <- paste0(
        inventory$id, " ", inventory$name,
        " (", round(inventory$distance, 1), " km), ",
        inventory$start_year, " to ", inventory$end_year,
        inventory$gap_label
      )

      # Prompt user
      cat("Choose a GHCN station or enter blank to cancel\n")
      selection <- utils::menu(options, title = "Available stations:")

      # Handle cancel
      if (!(selection %in% c(1:5))) stop("invalid station selected.")
    }else{
      # Use specified index (1-5) without prompt
      cat(
        "Using ", which_station, c("st", "nd", "rd", "th", "th")[which_station], " closest GHCN station",
        "\nstation: ", station$id[which_station], " ", station$name[which_station],
        ", (", station$lat[which_station], ", ", station$lng[which_station], "), ", round(station$distance[which_station], 1), " km from site_point\n",
        sep = ""
      )
      selection <- which_station
    }
    # Construct the NOAA station code (e.g., "723037-93765")
    station_code <- station$id[selection]
  } else {
    # If user directly passed a station code, use it as-is
    station_code <- as.character(which_station)
  }

  # Download NOAA met data for the specified years
  met_data <- worldmet::import_ghcn_hourly(station = station_code, year = years, append_codes = FALSE)

  # check if returned met data is blank
  if (is.null(met_data) || nrow(met_data) == 0) {
    message(
      "\nNo data found for station '", station_code, "' in year(s): ",
      paste(years, collapse = ", "),
      "\nPlease try a different station or a different year range."
    )
    return(NULL)
  }

  # Clean and restructure the data
  wind <- met_data %>%
    dplyr::mutate(id = .data$station_id,
           time = (.data$date), .before = 0) %>%
    dplyr::select("station_id", "time", wind_direction = "wd", wind_speed = "ws") %>%
    dplyr::mutate(
      year = lubridate::year(.data$time),
      month = lubridate::month(.data$time),
      day = lubridate::day(.data$time),
      .before = 3
    )


  return(wind)
}

#' Summarize a vector of years into ranges
#'
#' Internal helper to convert a numeric vector of years into a
#' human-readable string of ranges (e.g., "2000-2005, 2010").
#'
#' @param years A numeric vector of years.
#' @return A character string. Returns "" if the input is empty.
#' @noRd
summarize_years <- function(years) {
  if (length(years) == 0) return("") # Return blank for no gaps

  years <- sort(unique(years))
  groups <- cumsum(c(1, diff(years) != 1))

  out <- tapply(years, groups, function(x) {
    if (length(x) > 1) paste0(min(x), "-", max(x)) else as.character(x)
  })

  return(paste(out, collapse = ", "))
}
