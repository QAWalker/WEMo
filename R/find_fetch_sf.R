#' Calculate Wind Fetch Distances to Shoreline
#'
#' This function calculates wind fetch distances from water-based sites to shorelines
#' in multiple directions. Wind fetch represents the distance over which wind can act
#' on a water surface without obstruction by land, making it a critical parameter in
#' coastal and marine studies for wave modeling, habitat analysis, and environmental
#' assessments.
#'
#' @param site_layer An `sf` object containing point geometries representing the sites
#'   where fetch distances will be calculated. Sites located on land (intersecting
#'   with polygon_layer) are automatically filtered out.
#' @param polygon_layer An `sf` object containing polygon geometries representing
#'   land masses or shorelines. Fetch distances are calculated to the boundaries
#'   of these polygons.
#' @param directions A numeric vector of angular directions (in compass degrees)
#'   for which to calculate fetch distances. Default is 8 equally-spaced directions
#'   from 0° to 315° (every 45°). Directions follow compass convention where
#'   0° = North, 90° = East, 180° = South, 270° = West.
#' @param max_fetch A numeric value specifying the maximum fetch distance to consider
#'   (in map units). If no shoreline intersection is found within this distance,
#'   the fetch is recorded as the maximum value. Default is 10,000 units.
#'
#' @return An `sf` object containing linestring geometries representing fetch rays
#'   with the following columns:
#'   \describe{
#'     \item{geometry}{Linestring geometries showing fetch rays from each site}
#'     \item{direction}{Direction angle (in degrees) for each fetch ray}
#'     \item{fetch}{Calculated fetch distance (in map units)}
#'     \item{site_id}{Integer identifier linking each ray to its originating site}
#'   }
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Ensures coordinate reference systems (CRS) match between input layers
#'   \item Filters out sites located on land
#'   \item For each site and direction, creates a ray extending to max_fetch distance
#'   \item Identifies intersections between rays and shoreline polygons
#'   \item Calculates distances from sites to nearest intersection points
#'   \item Returns fetch rays as linestring geometries with associated measurements
#' }
#'
#' Direction angles are converted from compass degrees to mathematical radians
#' internally using a helper function [compassDegrees_to_radianDegrees].
#'
#' @examples
#' \dontrun{
#' # Basic usage with default 8 directions and 10,000 map units max fetch
#' fetch_rays <- find_fetch_sf(water_sites, land_polygons)
#'
#' # Custom directions - cardinal directions only
#' n_headings = 4
#' headings <- seq(0, 360 - 360 / n_headings, length.out = n_headings)
#'
#' fetch_rays <- find_fetch_sf(
#'   site_layer = sites,
#'   polygon_layer = coastline,
#'   directions = headings,
#'   max_fetch = 5000
#' )
#'
#' # High resolution analysis every 10 degrees
#' n_headings = 36
#' headings <- seq(0, 360 - 360 / n_headings, length.out = n_headings)
#' fetch_rays <- find_fetch_sf(
#'   site_layer = sites,
#'   polygon_layer = coastline,
#'   directions = headings,
#'   max_fetch = 1000
#' )
#'
#' }
#'
#' @export
find_fetch <- function(site_layer, polygon_layer,
                          directions = c(0, 45, 90, 135, 180, 225, 270, 315),
                          max_fetch = 10000) {
  # Ensure CRS match
  if (st_crs(site_layer) != st_crs(polygon_layer)) {
    warning("CRS mismatch: transforming site_layer to match polygon_layer")
    site_layer <- st_transform(site_layer, st_crs(polygon_layer))
  }

  # Ensure fetch has site column for calculations
  if("site" %in% names(site_layer)){
    site_var <- "site"
  } else {
    site_var <- names(site_layer)[which(startsWith(names(site_layer), "site"))]
    cat("using", site_var, "as site variable")
  }

  # Remove points on land
  sites_on_land <- st_intersects(site_layer, polygon_layer, sparse = FALSE)
  site_layer <- site_layer[!sites_on_land, ]

  # Store coordinates for speed
  coords <- st_coordinates(site_layer)

  # Setup progress bar
  # pb <- progress_bar$new(
  #   format = "  Processing [:bar] :percent eta: :eta",
  #   total = nrow(site_layer), clear = FALSE, width = 60
  # )

  # Generate rays from each site in all directions
  fetch_rays <- lapply(1:nrow(site_layer), function(i) {
    # pb$tick()
    site_geom <- st_geometry(site_layer[i,])
    site_coord <- coords[i, ]

    # heading = directions[1]
    site_fetch_rays <- lapply(directions, function(heading) {
      # for (heading in directions) {


      heading_rad <- compassDegrees_to_radianDegrees(heading) * (pi / 180)
      dest_x <- site_coord[1] + max_fetch * cos(heading_rad)
      dest_y <- site_coord[2] + max_fetch * sin(heading_rad)

      ray <- st_linestring(rbind(site_coord, c(dest_x, dest_y))) %>%
        st_sfc(crs = st_crs(site_layer)) %>%
        st_sf()

      # Find intersection with shoreline
      inter <- suppressWarnings(st_intersection(ray, polygon_layer))

      if (!inherits(inter, "sf") || nrow(inter) == 0 || all(st_is_empty(inter))) {
        final_ray <- st_linestring(rbind(site_coord, c(dest_x, dest_y)))
        dist_val <- max_fetch
      } else {
        # Cast to points
        inter_points <- suppressWarnings(tryCatch({
          inter_geom <- st_geometry(inter)
          if (!inherits(inter_geom, "sfc_POINT")) {
            st_cast(inter_geom, "POINT")
          } else {
            inter_geom
          }
        }, error = function(e) NULL))

        if (!is.null(inter_points) && length(inter_points) > 0 && !all(st_is_empty(inter_points))) {
          dists <- st_distance(site_geom, inter_points)
          min_idx <- which.min(dists[1,])
          closest_pt <- st_coordinates(inter_points[min_idx])

          final_ray <- st_linestring(rbind(site_coord, closest_pt))
          dist_val <- as.numeric(dists[min_idx])
        } else {
          final_ray <- st_linestring(rbind(site_coord, c(dest_x, dest_y)))
          dist_val <- max_fetch
        }
      }

      tibble(
        geometry = st_sfc(final_ray, crs = st_crs(site_layer)),
        direction = heading,
        fetch = dist_val,
        site = i
      )

    })
    site_fetch_rays <- do.call(rbind, site_fetch_rays)
    site_fetch_rays <- st_as_sf(site_fetch_rays)
    return(site_fetch_rays)
  })

  fetch_rays <- do.call(rbind, fetch_rays)
  fetch_rays <- st_as_sf(fetch_rays)
  return(fetch_rays)
}
