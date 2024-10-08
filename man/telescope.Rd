% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/telescope.R
\name{telescope}
\alias{telescope}
\title{A class to describe the basic properties of the observing telescope}
\usage{
telescope(
  type = "IFU",
  fov,
  aperture_shape = "circular",
  wave_range = c(3700, 5700),
  wave_centre,
  wave_res = 1.04,
  spatial_res,
  filter = "g",
  lsf_fwhm = 2.65,
  signal_to_noise = NA,
  method
)
}
\arguments{
\item{type}{String that describes the type of telescope you wish to create.
Current pre-loaded types include "SAMI", "MaNGA", "CALIFA", "MUSE" and
"Hector". Input is NOT case sensitive. If you wish to specify different
observing properties below, set \code{type = "IFU"} to define your own.}

\item{fov}{Numeric describing the field of view of the instrument in arcsec.}

\item{aperture_shape}{String to describe the shape of the IFU aperture.
Options include "circular", "hexagonal" or "square".}

\item{wave_range}{Numeric vector of length 2 describing the wave range of the
IFU (i.e. \code{c(wave_min, wave_max)}) in units of angstrom.}

\item{wave_centre}{Numeric describing the central wavelength of the
spectrograph used in the observation in units of angstrom. If unsupplied,
default is the exact centre of the provided `wave_range` parameter.}

\item{wave_res}{Numeric describing the wavelength resolution in angstrom.}

\item{spatial_res}{Numeric describing the size of spatial pixels in arcsec.}

\item{filter}{String describing the relevent filter through which luminosities
of individual particles are calculated.}

\item{lsf_fwhm}{Numeric describing the full-width half-maximum of the Gaussian
line spread function.}

\item{signal_to_noise}{Numeric describing the minimum signal-to-noise ratio per
angstrom. Default is NA, meaning that no noise will be added to the
observation.}

\item{method}{Providing backward compatibility for \code{method}
specification. String to describe whether cubes output are "spectral", "gas",
"sf gas" or "velocity" (as in SimSpin v1) along the z-axis. Default is
"spectral". For v2.1.6  onwards, \code{method} should be specified in
\link{build_datacube}. Support  for this input will remain in all
versions 2.X.X.}
}
\value{
Returns an object of class "telescope" that describes the properties
of the instrument doing the observation. Required to run
\code{build_datacube()}.
}
\description{
The purpose of this function is to generate the observing telescope whose
properties will be used to compute the properties of the observation.
Several telescope types are included and can be accessed using the
\code{type} parameter. Alternatively, you may describe the properties of
another telescope by specifying \code{type = "IFU"} and listing the other
required properties explicitly.
}
\examples{
telescope = telescope(type="SAMI")

}
