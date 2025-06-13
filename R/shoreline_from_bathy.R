#' Extract shoreline from a bathymetric raster using a contour threshold
#'
#' This function identifies land areas from a bathymetric raster by applying a contour threshold,
#' converts the result into polygon geometry, and optionally saves it to a shapefile.
#'
#' @param bathy A `SpatRaster` object representing bathymetric or elevation data.
#' @param contour A numeric value representing the elevation threshold to define land.
#' @param save_output Logical. If `TRUE`, the output polygon will be saved to a shapefile. Default is `FALSE`.
#' @param filename A string giving the name of the shapefile to save if `save_output = TRUE`. Default is `"shoreline.shp"`.
#'
#' @return A `SpatVector` polygon object representing the extracted land areas (shoreline).
#'
#' @examples
#' \dontrun{
#' # Example with a bathymetric raster
#' shoreline <- shoreline_from_bathy(bathy = bathy_raster, contour = 0)
#'
#' # Save the result to a shapefile
#' shoreline <- shoreline_from_bathy(bathy_raster, contour = 0, save_output = TRUE, filename = "shoreline.shp")
#' }
#'
#' @export
shoreline_from_bathy <- function(bathy, contour, save_output = FALSE, filename = "shoreline.shp") {
  # Convert raster to binary: 1 for land (bathy >= contour), NA for water
  land_bathy <- terra::ifel(bathy >= contour, 1, NA)

  # Convert the binary raster to polygons
  land_poly <- terra::as.polygons(land_bathy)

  # Optionally write the polygon to file
  if (save_output) {
    terra::writeVector(land_poly, filename = filename, overwrite = TRUE)
  }

  return(land_poly)
}
