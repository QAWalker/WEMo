#' Extract Bathymetry Along Fetch Ray
#'
#' Extract bathymetry values at regular intervals along a fetch ray from a raster.
#'
#' @description
#' This function generates sampling points along a fetch ray and extracts
#' corresponding bathymetry values from a raster dataset.
#'
#' @param bathy_raster A SpatRaster object containing bathymetry data
#' @param fetch_ray A spatial vector object (SpatVector) or sf object representing
#'   a linear geometry (fetch ray)
#' @param sample_dist Numeric. The distance between sampling points along the ray
#'
#' @return A list containing:
#'   \item{bathy}{Numeric vector of bathymetry values extracted at each sampling point}
#'   \item{distances}{Numeric vector of distances from ray start for each sampling point}
#'
#' @details
#' This function is a wrapper that combines point generation along the ray with
#' raster value extraction. It uses generate_points_along_ray() internally with
#' default settings (extra_at_start = TRUE).
#'
#' @export
extract_bathy_along_fetch <- function(bathy_raster, fetch_ray, sample_dist, extra_at_start = T){
  # Check if input is SpatVector, convert if necessary
  if (inherits(fetch_ray, "SpatVector")) {
    ray_vect <- fetch_ray
  } else{
    ray_vect <- vect(fetch_ray)
  }
  # Generate sampling points along the ray
  samp_points <- generate_points_along_ray(ray_vect, sample_dist, extra_at_start = extra_at_start)

  # Extract bathymetry values at sampling points
  bathy_values <- terra::extract(bathy_raster, samp_points$points, xy = F, ID = F, bind = T)

  # Return bathymetry values and distances
  return(
    list(
      bathy = terra::values(bathy_values)[[1]],
      distances = samp_points$distances
      )
  )
}
