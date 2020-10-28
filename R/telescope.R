# Author: Kate Harborne
# Date: 28/10/2020
# Title: telescope - a class for describing telescope properties
#
#'A class to describe the basic properties of the observing telescope
#'
#'The purpose of this function is to generate the observing telescope whose
#' properties will be used to compute the properties of the observation.
#' Several telescope types are included and can be accessed using the
#' \code{type} parameter. Alternatively, you may describe the properties of
#' another telescope by specifying \code{type = "IFU"} and listing the other
#' required properties explicitly.
#'
#'@param type String that describes the type of telescope you wish to create.
#' Current pre-loaded types include "SAMI", "MaNGA", "CALIFA", and "Hector".
#' Input is NOT case sensitive. If you wish to specify different observing
#' properties below, set \code{type = "IFU"} to define your own telescope.
#'@param fov Numeric describing the field of view of the instrument in arcsec.
#'@param aperture_shape String to describe the shape of the IFU aperture.
#' Options include "circular", "hexagonal" or "square".
#'@param wave_range Numeric vector of length 2 describing the wave range of the
#' IFU (i.e. \code{c(wave_min, wave_max)}.
#'@param spatial_res Numeric describing the size of spatial pixels in arcsec.
#'@param wave_res Numeric describing the wavelength resolution in angstrom.
#'@param lsf_fwhm Numeric describing the full-width half-maximum of the Gaussian
#' line spread function.
#'@return Returns an object of class "telescope" that describes the properties
#' of the instrument doing the observation. Required to run
#' \code{build_datacube()}.
#'@examples
#'telescope = telescope(type="SAMI")
#'

telescope = function(type="IFU", fov=15, aperture_shape="circular", wave_range=c(3700,5700),
                     spatial_res=0.5, wave_res=1.04, lsf_fwhm=2.65){

  if (length(wave_range)!=2){
    stop("Error: length(wave_range) should be 2. \n Please specify wave_range=c(wave_min, wave_max) and try again.")
  }
  if (wave_range[1] > wave_range[2]){
    wave_range = c(wave_range[2], wave_range[1]) # ensure that the wave_range is c(wave_min, wave_max)
  }
  if (stringr::str_to_lower(aperture_shape) != "circular" &
      stringr::str_to_lower(aperture_shape) != "hexagonal" &
      stringr::str_to_lower(aperture_shape) != "square"){
    stop("Error: Invalid aperture_shape. \n Please specify aperture_shape='circular', 'hexagonal' or 'square' and try again.")
  }

  if(stringr::str_to_upper(type) == "SAMI"){
    output = list(type           = "SAMI",
                  fov            = 15,
                  aperture_shape = "circular",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.65,
                  sbin           = floor(15 / 0.5))
  }

  if(stringr::str_to_upper(type) == "MANGA"){
    output = list(type           = "MaNGA",
                  fov            = 22,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.8,
                  sbin           = floor(22 / 0.25))
  }

  if(stringr::str_to_upper(type) == "HECTOR"){
    output = list(type           = "Hector",
                  fov            = 30,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.05,
                  wave_res       = 1.6,
                  lsf_fwhm       = 1.3,
                  sbin           = floor(30 / 0.05))

  }

  if(stringr::str_to_upper(type) == "CALIFA"){
    output = list(type           = "CALIFA",
                  fov            = 30,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 1,
                  wave_res       = 2,
                  lsf_fwhm       = 5.65,
                  sbin           = floor(30 / 1))

  }

  if(stringr::str_to_upper(type) == "IFU"){
    output = list(type           = "IFU",
                  fov            = fov,
                  aperture_shape = stringr::str_to_lower(aperture_shape),
                  wave_range     = wave_range,
                  spatial_res    = spatial_res,
                  wave_res       = wave_res,
                  lsf_fwhm       = lsf_fwhm,
                  sbin           = floor(fov / spatial_res))
  }

  return(output)
}
