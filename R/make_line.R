make_line <- function(X0, X1, Y0, Y1) {
  sf::st_linestring(matrix(c(X0, X1, Y0, Y1), 2, 2))
}

make_line_terra <- function(X0, X1, Y0, Y1) {
  terra::vect(matrix(c(X0, X1, Y0, Y1), 2, 2), type = "lines")
}
