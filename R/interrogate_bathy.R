#' Interrogate Bathymetry Along Fetch Rays
#'
#' Extract bathymetry data along multiple fetch rays and return updated input dataset.
#'
#' @description
#' This function processes multiple fetch rays, extracting bathymetry values along each
#' ray and storing the results as list columns in the input spatial dataframe.
#'
#' @param fetch A spatial dataframe (sf) containing fetch ray geometries
#' @param bathy_raster A SpatRaster object containing bathymetry data
#' @param sample_dist Numeric. The distance between sampling points along each ray (default: 10)
#' @param extra_at_start Logical. Controls point positioning when ray length doesn't
#'   divide evenly by sample_dist (default: TRUE). See generate_points_along_ray() for details
#'
#' @return A spatial dataframe with added list columns:
#'   \item{bathy}{List column containing bathymetry values for each fetch ray}
#'   \item{distances}{List column containing distance values for each fetch ray}
#'
#' @details
#' This function processes each row of the input fetch dataframe, extracting bathymetry
#' values along the corresponding geometry. Results are stored as list columns, allowing
#' each row to contain vectors of different lengths.
#'
#' @seealso
#' [extract_bathy_along_fetch] for single ray processing
#' [generate_points_along_ray] for point generation details
#'
#' @examples
#' # Process all fetch rays with 10m sampling
#' fetch_with_bathy <- interrogate_bathy(fetch_rays, bathy_raster, sample_dist = 10)
#'
#' @export
interrogate_bathy <- function(fetch, bathy_raster, sample_dist = 10, depths_or_elev = "elev", water_level = 0, extra_at_start = T){
  fetch_with_bathy <-
    lapply(seq_along(fetch$geometry), function(i){
      fetch_ray <- fetch[i, ]
      extracted_bathy <- extract_bathy_along_fetch(bathy, fetch_ray, sample_dist = sample_dist)

      if(depths_or_elev == 'depths') {
        depths <- extracted_bathy[["bathy"]] + water_level
      } else if(depths_or_elev == 'elev') {
        depths <- -1 * extracted_bathy[["bathy"]] + water_level
      }

      tibble::tibble(
        geometry = fetch_ray$geometry,
        bathy = list(extracted_bathy[["bathy"]]),
        distances = list(extracted_bathy[["distances"]]),
        depths = list(depths)
      ) %>%
        bind_cols(sf::st_drop_geometry(fetch_ray), .)
    }) %>%
    dplyr::bind_rows()
  return(fetch_with_bathy)
}
