find_fetch_and_effective_fetch <- function(points, shoreline, directions = list(), max_fetch = 10000){
  fetch <- find_fetch(site_layer = points, polygon_layer = shoreline, directions = directions, max_fetch = max_fetch)

  eff_fetch <- effective_fetch(fetch = fetch)

  return(eff_fetch)
}
