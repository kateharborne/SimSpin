% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cen_galaxy}
\alias{cen_galaxy}
\title{Centering the galaxy.}
\usage{
cen_galaxy(part_data)
}
\arguments{
\item{part_data}{A data.frame describing the particles ID, positions,
velocities and masses.}
}
\value{
A list containing the input data.frame with positions and velocities
 centered, as well as the median coordiantes and velocities of the original
 distribution
}
\description{
The purpose of this function is to centre the galaxy such that the origin of the
system lies at (0,0,0).
}
\examples{
  galaxy_data = data.frame("ID"=1:100, "x"=stats::rnorm(100),
                           "y"=stats::rnorm(100), "z"=stats::rnorm(100),
                           "vx"=stats::rnorm(100), "vy"=stats::rnorm(100),
                           "vz"=stats::rnorm(100), "Mass"=rep(1,100))
  centered_galaxy = cen_galaxy(galaxy_data)
}
