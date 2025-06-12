#' Convert compass degrees (direction) to radian degrees
#'
#' @param compassDegrees compass direction or heading in degrees
#'
#' @return a numerical vector
#' @export
#'
#' @examples
#' # North on a compass is 0 deg
#' compassDegrees_to_radianDegrees(0)
#'
#' # East on a compass is 90 deg
#' compassDegrees_to_radianDegrees(90)
#'
#' # South on a compass is 180 deg
#' compassDegrees_to_radianDegrees(180)
#'
#' # West on a compass is 270 deg
#' compassDegrees_to_radianDegrees(270)
#'
#' # North on a compass is ALSO 360 deg
#' compassDegrees_to_radianDegrees(360)
compassDegrees_to_radianDegrees <- function(compassDegrees){
  (450 - compassDegrees) %% 360
}
