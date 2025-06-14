% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_points_along_ray.R
\name{generate_points_along_ray}
\alias{generate_points_along_ray}
\title{Generate Points Along Ray}
\usage{
generate_points_along_ray(fetch_ray, sample_dist, extra_at_start = T)
}
\arguments{
\item{fetch_ray}{A spatial vector object (SpatVector) or sf object representing
a linear geometry (fetch ray)}

\item{sample_dist}{Numeric. The distance between sampling points along the ray}

\item{extra_at_start}{Logical. If TRUE (default), when ray length doesn't
divide evenly by sample_dist places extra point at the start. If FALSE,
places extra point at the end}
}
\value{
A list containing:
\item{points}{SpatVector of point geometries along the ray}
\item{distances}{Numeric vector of distances from previous point}
}
\description{
This function takes a linear geometry (fetch ray) and generates sampling points
at regular intervals along its length. The function handles the positioning of
points when the total length doesn't divide evenly by the sample distance.
}
\details{
Generate equally spaced sampling points along a fetch ray geometry.

When the total ray length doesn't divide evenly by sample_dist, there will be
a shorter distance between some points. The extra_at_start parameter controls
whether this shorter distance occurs near the site (TRUE) or away from it (FALSE).
For example, with a 90m ray and 25m sample distance, two points will be 15m apart.
}
