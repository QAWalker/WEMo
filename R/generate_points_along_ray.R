#' Generate Points Along Ray
#'
#' Generate equally spaced sampling points along a fetch ray geometry.
#'
#' @description
#' This function takes a linear geometry (fetch ray) and generates sampling points
#' at regular intervals along its length. The function handles the positioning of
#' points when the total length doesn't divide evenly by the sample distance.
#'
#' @param fetch_ray A spatial vector object (SpatVector) or sf object representing
#'   a linear geometry (fetch ray)
#' @param sample_dist Numeric. The distance between sampling points along the ray
#' @param extra_at_start Logical. If TRUE (default), when ray length doesn't
#'    divide evenly by sample_dist places extra point at the start. If FALSE,
#'    places extra point at the end
#'
#' @return A list containing:
#'   \item{points}{SpatVector of point geometries along the ray}
#'   \item{distances}{Numeric vector of distances from previous point}
#'
#' @details
#' When the total ray length doesn't divide evenly by sample_dist, there will be
#' a shorter distance between some points. The extra_at_start parameter controls
#' whether this shorter distance occurs near the site (TRUE) or away from it (FALSE).
#' For example, with a 90m ray and 25m sample distance, two points will be 15m apart.
#'
#' @export
generate_points_along_ray <- function(fetch_ray, sample_dist, extra_at_start = T) {
  # Check if input is SpatVector, convert if necessary
  if (inherits(fetch_ray, "SpatVector")) {
    ray_vect <- fetch_ray
  } else {
    ray_vect <- vect(fetch_ray)
  }

  # Get line coordinates
  coords <- geom(ray_vect)[, c("x", "y")]

  # Calculate total ray length
  total_length <- terra::perim(ray_vect)

  # Create target distances - equally spaced points and the final point
  # Allows the user to dictate if the short distance between two points is close
  # to the site point (extra_at_start = TRUE) or furthest away (extra_at_start = FALSE)
  # e.g. if the fetch ray is 90 meters long and sample_dist is 25 meters there will be
  # two points that are 15 meters away - default is to make them closest to the site
  if (extra_at_start) {
    target_distances <- unique(sort(c(seq(from = total_length, to = 0, by = -sample_dist), 0)))
  } else {
    target_distances <- unique(sort(c(seq(from = 0, to = total_length, by = sample_dist), total_length)))
  }

  # make the vector start with the far point
  distances <- rev(diff(target_distances))
  target_distances <- rev(target_distances)

  # Interpolate coordinates at target distances
  x_interp <- approx(c(0, total_length), coords[, 1], xout = target_distances)$y
  y_interp <- approx(c(0, total_length), coords[, 2], xout = target_distances)$y

  # Create SpatVector points
  sample_coords <- cbind(x_interp, y_interp)
  ray_points_vect <- vect(sample_coords, type = "points", crs = crs(ray_vect))

  return(list(points = ray_points_vect, distances = distances))
}
