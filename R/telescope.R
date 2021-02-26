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
#' Current pre-loaded types include "SAMI", "MaNGA", "CALIFA", "MUSE" and
#' "Hector". Input is NOT case sensitive. If you wish to specify different
#' observing properties below, set \code{type = "IFU"} to define your own.
#'@param method String to describe whether cubes output are "spectral" or
#' "velocity" (as in SimSpin v1) along the z-axis. Default is "spectral".
#'@param fov Numeric describing the field of view of the instrument in arcsec.
#'@param aperture_shape String to describe the shape of the IFU aperture.
#' Options include "circular", "hexagonal" or "square".
#'@param wave_range Numeric vector of length 2 describing the wave range of the
#' IFU (i.e. \code{c(wave_min, wave_max)}.
#'@param wave_centre Numeric describing the central wavelength of the
#' spectrograph used in the observation. If unsupplied, default is the exact
#' centre of the provided `wave_range` parameter.
#'@param wave_res Numeric describing the wavelength resolution in angstrom.
#'@param spatial_res Numeric describing the size of spatial pixels in arcsec.
#'@param filter String describing the relevent filter through which luminosities
#' of individual particles are calculated.
#'@param lsf_fwhm Numeric describing the full-width half-maximum of the Gaussian
#' line spread function.
#'@param signal_to_noise Numeric describing the minimum signal-to-noise ratio per
#' angstrom.
#'@return Returns an object of class "telescope" that describes the properties
#' of the instrument doing the observation. Required to run
#' \code{build_datacube()}.
#'@examples
#'telescope = telescope(type="SAMI")
#'

telescope = function(type="IFU", method="spectral", fov=15, aperture_shape="circular", wave_range=c(3700,5700),
                     wave_centre, wave_res=1.04, spatial_res=0.5, filter="r", lsf_fwhm=2.65, signal_to_noise = 10){

  if (length(wave_range)!=2){
    stop("Error: length(wave_range) should be 2. \n Please specify wave_range = c(wave_min, wave_max) and try again.")
  }
  if (wave_range[1] > wave_range[2]){
    wave_range = c(wave_range[2], wave_range[1]) # ensure that the wave_range is c(wave_min, wave_max)
  }
  if (stringr::str_to_lower(aperture_shape) != "circular" &
      stringr::str_to_lower(aperture_shape) != "hexagonal" &
      stringr::str_to_lower(aperture_shape) != "square"){
    stop("Error: Invalid aperture_shape. \n Please specify aperture_shape = 'circular', 'hexagonal' or 'square' and try again.")
  }
  if (stringr::str_to_lower(method) != "spectral" &
      stringr::str_to_lower(method) != "velocity" &
      stringr::str_to_lower(method) != "gas"){
    stop("Error: Invalid method. \n Please specify method = 'spectral' or 'velocity' and try again.")
  }
  if (missing(wave_centre)){
    wave_centre = wave_range[1] + (diff(wave_range)/2)
  }

  if(stringr::str_to_upper(type)  == "SAMI"){
    fov = 15
    spatial_res = 0.5
    output = list(type            = "SAMI",
                  method          = method,
                  fov             = fov,
                  aperture_shape  = "circular",
                  wave_range      = c(3700,5700),
                  wave_centre     = 4800,
                  spatial_res     = spatial_res,
                  filter          = ProSpect::filt_g_SDSS,
                  wave_res        = 1.04,
                  lsf_fwhm        = 2.65,
                  signal_to_noise = 10,
                  sbin            =  floor(fov/spatial_res))
  }

  if(stringr::str_to_upper(type)  == "MANGA"){
    fov = fov
    spatial_res = 0.5
    output = list(type            = "MaNGA",
                  method          = method,
                  fov             = fov,
                  aperture_shape  = "hexagonal",
                  wave_range      = c(3700,5700),
                  wave_centre     = 4700,
                  spatial_res     = spatial_res,
                  filter          = ProSpect::filt_g_SDSS,
                  wave_res        = 1.04,
                  lsf_fwhm        = 2.8,
                  signal_to_noise = 10,
                  sbin            =  floor(fov/spatial_res))
  }

  if(stringr::str_to_upper(type)  == "MUSE"){
    fov = fov
    spatial_res = 0.2
    output = list(type            = "MUSE",
                  method          = method,
                  fov             = fov,
                  aperture_shape  = "square",
                  wave_range      = c(4650,9300),
                  wave_centre     = 6975,
                  spatial_res     = spatial_res,
                  filter          = ProSpect::filt_r_SDSS,
                  wave_res        = 1.25,
                  lsf_fwhm        = 2.63,
                  signal_to_noise = 10,
                  sbin            = floor(fov/spatial_res))

  }

  if(stringr::str_to_upper(type)  == "HECTOR"){
    fov = 30
    spatial_res = 0.1
    output = list(type            = "Hector",
                  method          = method,
                  fov             = fov,
                  aperture_shape  = "hexagonal",
                  wave_range      = c(3700,5700),
                  wave_centre     = 4700,
                  spatial_res     = spatial_res,
                  filter          = ProSpect::filt_g_SDSS,
                  wave_res        = 1.6,
                  lsf_fwhm        = 1.3,
                  signal_to_noise = 10,
                  sbin            =  floor(fov/spatial_res))

  }

  if(stringr::str_to_upper(type)  == "CALIFA"){
    fov = 30
    spatial_res = 0.8
    output = list(type            = "CALIFA",
                  method          = method,
                  fov             = fov,
                  aperture_shape  = "hexagonal",
                  wave_range      = c(3700,5700),
                  wave_centre     = 4700,
                  spatial_res     = spatial_res,
                  filter          = ProSpect::filt_g_SDSS,
                  wave_res        = 2,
                  lsf_fwhm        = 5.65,
                  signal_to_noise = 10,
                  sbin            = floor(fov/spatial_res))

  }

  if(stringr::str_to_upper(type)  == "IFU"){

    if (stringr::str_to_lower(filter) == "r"){
      filter = ProSpect::filt_r_SDSS
    } else if (stringr::str_to_lower(filter) == "u"){
      filter = ProSpect::filt_u_SDSS
    } else if (stringr::str_to_lower(filter) == "g"){
      filter = ProSpect::filt_g_SDSS
    } else if (stringr::str_to_lower(filter) == "i"){
      filter = ProSpect::filt_i_SDSS
    } else if (stringr::str_to_lower(filter) == "z"){
      filter = ProSpect::filt_z_SDSS
    } else {
      stop("Error: Invalid filter. \n Please specify filter = 'r', 'u' or 'g', 'i' or 'z' and try again.")
    }

    output = list(type            = "IFU",
                  fov             = fov,
                  method          = method,
                  aperture_shape  = stringr::str_to_lower(aperture_shape),
                  wave_range      = wave_range,
                  wave_centre     = wave_centre,
                  wave_res        = wave_res,
                  spatial_res     = spatial_res,
                  filter          = filter,
                  lsf_fwhm        = lsf_fwhm,
                  signal_to_noise = signal_to_noise,
                  sbin            = floor(fov / spatial_res))
  }

  return(output)
}
