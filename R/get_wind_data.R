#' Retrieve and clean wind data for WEMo
#'
#' This function is a wrapper around the [worldmet] package's [worldmet::getMeta] and [worldmet::importNOAA]
#' functions. It finds up to 5 nearby NOAA weather stations based on a geographic point, allows the user
#' to select a station (interactively or by index/code), and downloads wind direction and speed data
#' for specified years.
#'
#' @param site_point A spatial point object (either `sf` or coercible to `sf`) indicating the target location.
#' @param years A vector of years (e.g., 2002:2024) for which to download wind data.
#' @param which_station Either:
#'   - "ask": interactively choose from 5 closest stations;
#'   - an integer (1–5): pick the nth closest station;
#'   - a string station code in the form "USAF-WBAN".
#'
#' @return A data frame with cleaned wind direction and wind speed, along with station code, timestamp,
#'         and date components (year, month, day).
#'
#' @examples
#' require(sf)
#' require(dplyr)
#'
#' # Create a POINT geometry from coordinates and cast explicitly
#' site_point <- st_sf(
#'   geometry = st_sfc(st_point(c(-76.67587, 34.71413))),
#'   crs = 4326
#'   )
#'
#' # 1. Prompt user to select a station interactively
#' \dontrun{
#' wind_data_ask <- get_wind_data(site_point, years = 2022:2023, which_station = "ask")
#' }
#'
#' # 2. Automatically use the 2nd closest station
#' \dontrun{
#' wind_data_index <- get_wind_data(site_point, years = 2022:2023, which_station = 2)
#' }
#'
#' # 3. Manually specify a station code using getMeta
#' \dontrun{
#' # this way you don't need a site_point
#' wind_data_manual <- get_wind_data(site_point = NULL, years = 2022:2023, which_station = "723037-93765")
#' }
#'
#' @export
get_wind_data <- function(site_point, years, which_station = 'ask') {
  # If user asks to choose or uses index 1–5, begin by locating closest stations
  if (which_station %in% c("ask", 1:5)) {
    # Ensure site_point is an sf object
    if (!inherits(site_point, "sf")) {
      tryCatch({
        site_point <- sf::st_as_sf(site_point)
      }, error = function(e) {
        stop("`site_point` needs to be sf object or coerable to sf: ", e$message)
      })
    }

    # Reproject to WGS84 for compatibility with worldmet
    target_crs <- sf::st_crs("EPSG:4326")
    if (sf::st_crs(site_point) != target_crs) {
      site_point <- sf::st_transform(site_point, target_crs)
    }

    # Extract latitude and longitude
    LAT <- sf::st_coordinates(site_point)[[1,2]]
    LON <- sf::st_coordinates(site_point)[[1,1]]

    # Get metadata for the 5 closest NOAA stations
    station <- worldmet::getMeta(lat = LAT, lon = LON, n = 5)
    if(which_station == "ask"){
      # Display options for user to choose from
      options <- paste0(
        station$usaf, '-', station$wban, " ",station$station, " (", round(station$dist, 1), " km), ",
        station$begin, " to ", station$end
      )

      # Prompt user
      cat("Choose a met station\n")
      selection <- menu(options, title = "Available stations:")

      # Handle cancel
      if (!(selection %in% c(1:5))) stop("invalid station selected.")
    }else{
      # Use specified index (1–5) without prompt
      cat("Using Met station", which_station, paste(station$usaf[selection], station$wban[selection], sep = '-'), paste0(
        station$station[which_station], " (", round(station$dist[which_station], 1), " km), ",
        station$begin[which_station], " to ", station$end[which_station]
      ), "\n")
      selection <- which_station
    }
    # Construct the NOAA station code (e.g., "723037-93765")
    station_code <- paste(station$usaf[selection], station$wban[selection], sep = '-')
  } else {
    # If user directly passed a station code, use it as-is
    station_code <- as.character(which_station)
  }

  # Download NOAA met data for the specified years
  met_data <- worldmet::importNOAA(code = station_code, year = years)

  # Clean and restructure the data
  wind <- met_data %>%
    dplyr::mutate(id = code,
           time = (date), .before = 0) %>%
    dplyr::select(code, time, wind_direction = wd, wind_speed = ws) %>%
    dplyr::mutate(
      year = lubridate::year(time),
      month = lubridate::month(time),
      day = lubridate::day(time),
      .before = 3
    )

  # wind <- wind %>%
  #   mutate(wind_direction = round(wind_direction/(360/36), 0)*(360/36)) %>%
  #   filter()

  return(wind)
}

