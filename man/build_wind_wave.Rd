% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_wind_wave.R
\name{build_wind_wave}
\alias{build_wind_wave}
\title{Generate and Propagate a Wind Wave Over Variable Bathymetry}
\usage{
build_wind_wave(fetch, depths = list(), distances = list(), wind_speed)
}
\arguments{
\item{fetch}{Numeric. Total wind fetch distance (in meters).}

\item{depths}{Numeric list. Water depths (in meters) at each segment; must be one element longer than \code{distances}.}

\item{distances}{Numeric list. Horizontal distances (in meters) between each depth point.}

\item{wind_speed}{Numeric. Wind speed (in m/s).}
}
\value{
A data frame with the following columns:
\describe{
\item{wave_height_final}{Final wave height at the end of the fetch (in meters).}
\item{WEI}{Wave Energy Index (proportional to energy flux).}
\item{wave_period}{Estimated wave period (in seconds).}
\item{wave_number}{Average wave number (radians per meter).}
\item{celerity_final}{Final wave celerity (phase speed, in m/s).}
\item{nnumber_final}{Final group velocity coefficient.}
}
}
\description{
This function builds a wind-wave over input bathymetry and distances returning
an estimate of final wave height and related wave characteristics
(e.g., wave energy index, wave period, wave number, celerity)
}
\examples{
fetch <- 90
depths <- list(1, 2, 3, 5, 5)
distances <- list(25, 25, 25, 15)
wind_speed <- 10
wind_proportion <- 0.25
result <- wind_wave_height_new(fetch, depths, distances, wind_speed)

}
