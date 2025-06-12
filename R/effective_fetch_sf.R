#' Calculate Effective Fetch for Wind Energy Transfer
#'
#' This function calculates the effective fetch (cosine-weighted fetch) from input fetch lines at a site,
#' by weighting each fetch distance by the cosine of its angular difference from a target direction.
#' Due to irregular shoreline in coastal waters the simple fetch length in a given
#' compass direction is not very effective since the width of fetch places a
#' substantial restriction on the fetch length (USCOE 1977). Therefore, the fetch
#' is modified by taking the cosine weighted average of all the rays within a
#' certain sector on either side of the fetch ray defined by angle θ. The cosine
#' weighted fetch is named as effective fetch.
#'
#' @param fetch An sf object containing fetch data with geometry `LINESTRING`,
#'   direction, fetch columns, and site identifiers. Must include
#'   columns for direction (in compass degrees), fetch distance, and site ID.
#' @param wind_energy_transfer_degrees Numeric. The half-angle (in degrees) of
#'   the wind energy transfer cone. Wind directions within ±this angle from the
#'   target direction contribute to the effective fetch calculation. Default is 45.
#'
#' @return An sf object (LINESTRING) containing effective fetch rays for each
#'   site and direction, with the following columns:
#'   \itemize{
#'     \item \code{direction} - Wind direction in compass degrees (0-360)
#'     \item \code{efetch} - Calculated effective fetch distance
#'     \item Additional columns from original fetch data
#'     \item \code{geometry} - LINESTRING geometry from center to effective fetch endpoint
#'   }
#'
#' @details
#' Effective fetch in `WEMo` is calculated by summing the product of fetch length
#' and cosine of the angle of departure from the ith heading over each of n number
#' of fetch rays and dividing by the sum of the cosine of all angles.
#' Effective fetch rays look trimmed down in the figure compared to fetch rays.
#' This method is based on assumption that wind moving over water surface transfers
#' energy to the water in the direction of the wind and in all directions within
#' 4/π radians (45°) on either side of the wind direction. (USCOE 1977)
#'
#' The function processes each site independently and calculates effective fetch
#' using a cosine-weighted average of all fetch rays within the specified angular
#' cone around each target direction. The calculation follows:
#'
#' \deqn{Eff F_{i} = \frac{\sum_{j=-n}^{n} F_j \cos(\theta_j)}{\sum_{j= -n}^{n} \cos(\theta_j)}}
#'
#' Where:
#' \itemize{
#'   \item \deqn{Eff F_i} Effective fetch for the ith direction fetch ray
#'   \item \deqn{F_j} length for j radiating fetch ray after clipping to shoreline
#'   \item \deqn{\theta_j} angle between the ith fetch ray and the jth ray
#'   \item \deqn{n} number of rays selected by user
#' }
#'
#' The function handles directional wrap-around by expanding the direction space
#' by ±360° to ensure proper calculation near 0°/360° boundary.
#'
#' @section Site Identification:
#' The function automatically detects site identifier columns by looking for:
#' \itemize{
#'   \item A column named exactly "site"
#'   \item The first column with a name starting with "site"
#' }
#'
#' @section Required Columns:
#' The input sf object must contain:
#' \itemize{
#'   \item \code{direction} - Wind direction in compass degrees
#'   \item \code{fetch} Fetch distance values
#'   \item Site identifier column (see Site Identification section)
#'   \item Valid `sf` geometry (typically `LINESTRING` or `POINT`)
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default 45-degree cone
#' effective_fetch <- effective_fetch_sf(fetch_data)
#'
#' # Use narrower energy transfer cone
#' effective_fetch <- effective_fetch_sf(fetch_data, wind_energy_transfer_degrees = 30)
#'
#' # Plot results
#' library(ggplot2)
#' ggplot(effective_fetch) +
#'   geom_sf(aes(color = efetch)) +
#'   facet_wrap(~site)
#' }
#'
effective_fetch <- function(fetch, wind_energy_transfer_degrees = 45) {
  # Check input type
  stopifnot("fetch must be an sf object" = inherits(fetch, "sf"))
  stopifnot("wind_energy_transfer_degrees must be positive" = wind_energy_transfer_degrees > 0)

  # Ensure fetch has site column for calculations
  if("site" %in% names(fetch)){
    site_var <- "site"
  } else {
    site_var <- names(fetch)[which(startsWith(names(fetch), "site"))]
    cat("using", site_var, "as site variable")
  }

  eff_fetch <- lapply(unique(fetch[[site_var]]), function(site){
    fetch_filtered <- fetch[fetch[[site_var]]==site, ]

    # Get center coordinates
    center_coords <- st_coordinates(fetch_filtered)[1,c("X", "Y")]

    # Create a dataframe without geometry for calculations
    fetch_df <- st_drop_geometry(fetch_filtered)

    # Expand directions by ±360° to handle edge cases
    fetch_df <- dplyr::bind_rows(fetch_df,
                                 fetch_df %>% dplyr::mutate(direction = .data$direction - 360),
                                 fetch_df %>% dplyr::mutate(direction = .data$direction + 360)) %>%
      dplyr::arrange(direction) %>%
      dplyr::distinct()

    # Calculate effective fetch for each direction
    eff_fetch <- lapply(unique(fetch_df$direction), function(theta) {
      # Filter for fetch rays within the defined wind energy transfer cone
      df <- fetch_df[which(fetch_df$direction >= theta - wind_energy_transfer_degrees &
                             fetch_df$direction <= theta + wind_energy_transfer_degrees), ]

      # Calculate weighted fetch
      df <- df %>%
        dplyr::mutate(a = fetch * cos((direction - theta) * pi / 180),
                      b = cos((direction - theta) * pi / 180)) %>%
        dplyr::summarize(suma = sum(a),
                         sumb = sum(b)) %>%
        dplyr::mutate(direction = theta,
                      efetch = suma/sumb,
                      suma = NULL, sumb = NULL)
      return(df)
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(direction)

    # Join with original data
    eff_fetch <- eff_fetch %>%
      dplyr::left_join(st_drop_geometry(fetch_filtered), ., by = "direction") %>%
      dplyr::mutate(X0 = center_coords[[1]],
                    Y0 = center_coords[[2]])

    # Convert compass degrees to radians and calculate endpoints
    eff_fetch <- eff_fetch %>%
      dplyr::mutate(radianDegrees = compassDegrees_to_radianDegrees(direction),
                    f = ifelse(fetch > efetch, efetch, fetch),
                    xdist = f * cos(radianDegrees * pi / 180),
                    ydist = f * sin(radianDegrees * pi / 180),
                    X1 = X0 + xdist,
                    Y1 = Y0 + ydist,
                    efetch = f,
                    f = NULL) %>%
      filter(direction>=0, direction<360)

    # Create list of LINESTRING geometries
    lines_list <- lapply(1:nrow(eff_fetch), function(i) {
      make_line(X0 = eff_fetch$X0[i], X1 = eff_fetch$X1[i], Y0 = eff_fetch$Y0[i], Y1 = eff_fetch$Y1[i])
    })

    eff_fetch_sf <- st_as_sf(
      eff_fetch,
      geometry = st_sfc(lines_list, crs = st_crs(fetch))
    )

    eff_fetch_sf %>%
      mutate(X0 = NULL,
             Y0 = NULL,
             X1 = NULL,
             Y1 = NULL,
             radianDegrees = NULL,
             xdist = NULL,
             ydist = NULL) %>%
      return()
  })

  # Combine all sites
  eff_fetch_sf <- do.call(rbind, eff_fetch)

  eff_fetch_sf[,which(names(eff_fetch_sf) %in% c(site_var, 'direction', 'fetch', 'efetch'))]
  return(eff_fetch_sf)
}
