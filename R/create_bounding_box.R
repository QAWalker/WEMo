#' Create a bounding box around a point with an expansion distance
#'
#' @param site_point a SpatVector (or coercible to a SpatVector)
#' @param expansion_dist a vector of length 1 or 2 to add to the limits of the bounding box. If length 1, the value is used for both x and y expansion
#'
#' @return a SpatExtent object
#' @export
#'
#' @examples
#' site_point <- terra::vect(cbind(x = 406825, y = 4372214), crs = "EPSG:26918")
#' create_bounding_box(site_point, 1000)
create_bounding_box <- function(site_point, expansion_dist) {
  if(length(expansion_dist)==1){
    expansion_dist[2] <- expansion_dist[1]
  }
  bbox <- terra::ext(site_point)  # Get the bounding box of the site point
  xmin <- bbox$xmin - expansion_dist[1]  # Expand xmin by the expansion distance
  xmax <- bbox$xmax + expansion_dist[1]  # Expand xmax by the expansion distance
  ymin <- bbox$ymin - expansion_dist[2]  # Expand ymin by the expansion distance
  ymax <- bbox$ymax + expansion_dist[2]  # Expand ymax by the expansion distance
  return(terra::ext(xmin, xmax, ymin, ymax))  # Return the expanded bounding box
}
