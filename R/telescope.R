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
#'@param fov Numeric describing the field of view of the instrument in arcsec.
#'@param aperture_shape String to describe the shape of the IFU aperture.
#' Options include "circular", "hexagonal" or "square".
#'@param wave_range Numeric vector of length 2 describing the wave range of the
#' IFU (i.e. \code{c(wave_min, wave_max)}) in units of angstrom.
#'@param wave_centre Numeric describing the central wavelength of the
#' spectrograph used in the observation in units of angstrom. If unsupplied,
#' default is the exact centre of the provided `wave_range` parameter.
#'@param wave_res Numeric describing the wavelength resolution in angstrom.
#'@param spatial_res Numeric describing the size of spatial pixels in arcsec.
#'@param filter String describing the relevent filter through which luminosities
#' of individual particles are calculated.
#'@param lsf_fwhm Numeric describing the full-width half-maximum of the Gaussian
#' line spread function.
#'@param signal_to_noise Numeric describing the minimum signal-to-noise ratio per
#' angstrom.
#'@param method Providing backward compatibility for \code{method}
#' specification. String to describe whether cubes output are "spectral", "gas",
#' "sf gas" or "velocity" (as in SimSpin v1) along the z-axis. Default is
#' "spectral". For v2.1.6  onwards, \code{method} should be specified in
#' \link{build_datacube}. Support  for this input will remain in all
#' versions 2.X.X.
#'@return Returns an object of class "telescope" that describes the properties
#' of the instrument doing the observation. Required to run
#' \code{build_datacube()}.
#'@examples
#'telescope = telescope(type="SAMI")
#'

telescope = function(type="IFU", fov, aperture_shape="circular", wave_range=c(3700,5700),
                     wave_centre, wave_res=1.04, spatial_res, filter="r", lsf_fwhm=2.65,
                     signal_to_noise = 10, method){

  if (!missing(method)){

    if (method != "spectral" &
        method != "velocity" &
        method != "gas" &
        method != "sf gas" ){
      stop("Error: Invalid method. \n Please specify method = 'spectral', 'velocity', 'gas' or 'sf gas' and try again.")
    } else {
      warning(">>> WARNING! >>> \n
             `method` is now specified within the build_datacube function directly,
              rather than within the telescope() class. \n
              Support for this input will remain in versions 2.X.X, but please consider
              updating your code.")

      method = stringr::str_to_lower(method)
    }
  }

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
  if (missing(wave_centre)){
    wave_centre = wave_range[1] + (diff(wave_range)/2)
  }
  if (stringr::str_to_lower(filter) == "r"){
    filter = SimSpin::filt_r_SDSS
  } else if (stringr::str_to_lower(filter) == "u"){
    filter = SimSpin::filt_u_SDSS
  } else if (stringr::str_to_lower(filter) == "g"){
    filter = SimSpin::filt_g_SDSS
  } else if (stringr::str_to_lower(filter) == "i"){
    filter = SimSpin::filt_i_SDSS
  } else if (stringr::str_to_lower(filter) == "z"){
    filter = SimSpin::filt_z_SDSS
  } else {
    stop("Error: Invalid filter. \n Please specify filter = 'r', 'u' or 'g', 'i' or 'z' and try again.")
  }

  if (missing(method)){
    if(stringr::str_to_upper(type)  == "SAMI"){
      fov = 15
      spatial_res = 0.5
      output = list(type            = "SAMI",
                    fov             = fov,
                    aperture_shape  = "circular",
                    wave_range      = c(3750,5750),
                    wave_centre     = 4800,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.04,
                    lsf_fwhm        = 2.65,
                    signal_to_noise = signal_to_noise,
                    sbin            =  floor(fov/spatial_res))
    }

    if(stringr::str_to_upper(type)  == "MANGA"){

      spatial_res = 0.5

      if (missing(fov)){
        fov = 12
      } else {
        fov_pos = c(12, 17, 22, 27, 32)
        if (!fov %in% fov_pos){
          suggested_fov = fov_pos[which((fov_pos - fov) > 0)[1]]
          warning(paste0(">>> WARNING! >>> \n
                         `fov` specified is not an option for telescope type `MaNGA`. \n
                         `fov` options include 12, 17, 22, 27, or 32. \n
                         Setting 'fov' to the nearest available option. \n
                         fov = ", suggested_fov))
          fov = suggested_fov
        } # catching any inputs that are not suitable for MaNGA telescope
      }

      output = list(type            = "MaNGA",
                    fov             = fov,
                    aperture_shape  = "hexagonal",
                    wave_range      = c(3600,6350),
                    wave_centre     = 4700,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.04,
                    lsf_fwhm        = 2.85,
                    signal_to_noise = signal_to_noise,
                    sbin            = floor(fov/spatial_res))
    }

    if(stringr::str_to_upper(type)  == "MUSE"){

      if (missing(spatial_res)){
        spatial_res = 0.2
      } else {
        pos_spatial_res = c(0.025, 0.2)
        if (!spatial_res %in% pos_spatial_res){
          spatial_res = 0.2
          warning(">>> WARNING! >>> \n
                  `spatial_res` specified is not an option for telescope type `MUSE`. \n
                  `spatial_res` options include 0.025 (NFM) or 0.2 (WFM). \n
                  Setting `spatial_res` to the default. \n
                  spatial_res = 0.2 ")

        }
      }

      if (missing(fov)){
        fov = 5
      } else {
        if (fov > 60){
          warning(paste0(">>> WARNING! >>> \n
                         `fov` specified is not an option for telescope type 'MUSE'. \n
                         `fov` options must be less than 60. \n
                         Setting `fov` to the nearest available option. \n
                         fov = 60"))
          fov = 60
        } # catching any inputs that are not suitable for MaNGA telescope
      }

      output = list(type            = "MUSE",
                    fov             = fov,
                    aperture_shape  = "square",
                    wave_range      = c(4700.15,9351.4),
                    wave_centre     = 6975,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.25,
                    lsf_fwhm        = 2.51,
                    signal_to_noise = signal_to_noise,
                    sbin            = floor(fov/spatial_res))

    }

    if(stringr::str_to_upper(type)  == "HECTOR"){
      fov = 30
      spatial_res = 0.1
      output = list(type            = "Hector",
                    fov             = fov,
                    aperture_shape  = "hexagonal",
                    wave_range      = c(3720,5910),
                    wave_centre     = 4815,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.6,
                    lsf_fwhm        = 1.3,
                    signal_to_noise = signal_to_noise,
                    sbin            =  floor(fov/spatial_res))

    }

    if(stringr::str_to_upper(type)  == "CALIFA"){
      fov = 74
      spatial_res = 1
      output = list(type            = "CALIFA",
                    fov             = fov,
                    aperture_shape  = "hexagonal",
                    wave_range      = c(3700,4750),
                    wave_centre     = 4225,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 2.7,
                    lsf_fwhm        = 2.7,
                    signal_to_noise = signal_to_noise,
                    sbin            = floor(fov/spatial_res))

    }

    if(stringr::str_to_upper(type)  == "IFU"){

      if (missing(fov)){fov = 15} # setting to the default parameter expected
      if (missing(spatial_res)){spatial_res = 0.5}

      output = list(type            = "IFU",
                    fov             = fov,
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

  } else {

    if(stringr::str_to_upper(type)  == "SAMI"){
      fov = 15
      spatial_res = 0.5
      output = list(type            = "SAMI",
                    method          = method,
                    fov             = fov,
                    aperture_shape  = "circular",
                    wave_range      = c(3750,5750),
                    wave_centre     = 4800,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.04,
                    lsf_fwhm        = 2.65,
                    signal_to_noise = signal_to_noise,
                    sbin            =  floor(fov/spatial_res))
    }

    if(stringr::str_to_upper(type)  == "MANGA"){

      spatial_res = 0.5

      if (missing(fov)){
        fov = 12
      } else {
        fov_pos = c(12, 17, 22, 27, 32)
        if (!fov %in% fov_pos){
          suggested_fov = fov_pos[which((fov_pos - fov) > 0)[1]]
          warning(paste0(">>> WARNING! >>> \n
                         `fov` specified is not an option for telescope type `MaNGA`. \n
                         `fov` options include 12, 17, 22, 27, or 32. \n
                         Setting `fov` to the nearest available option. \n
                         fov = ", suggested_fov))
          fov = suggested_fov
        } # catching any inputs that are not suitable for MaNGA telescope
      }

      output = list(type            = "MaNGA",
                    method          = method,
                    fov             = fov,
                    aperture_shape  = "hexagonal",
                    wave_range      = c(3600,6350),
                    wave_centre     = 4700,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.04,
                    lsf_fwhm        = 2.85,
                    signal_to_noise = signal_to_noise,
                    sbin            =  floor(fov/spatial_res))
    }

    if(stringr::str_to_upper(type)  == "MUSE"){

      if (missing(spatial_res)){
        spatial_res = 0.2
      } else {
        pos_spatial_res = c(0.025, 0.2)
        if (!spatial_res %in% pos_spatial_res){
          spatial_res = 0.2
          warning(">>> WARNING! >>> \n
                  `spatial_res` specified is not an option for telescope type `MUSE`. \n
                  `spatial_res` options include 0.025 (NFM) or 0.2 (WFM). \n
                  Setting 'spatial_res' to the default. \n
                  spatial_res = 0.2 ")

        }
      }

      if (missing(fov)){
          fov = 5
        } else {
          if (fov > 60){
            warning(paste0(">>> WARNING! >>> \n
                         `fov` specified is not an option for telescope type `MUSE`. \n
                         `fov` options must be less than 60. \n
                         Setting `fov` to the nearest available option. \n
                         fov = 60"))
            fov = 60
          } # catching any inputs that are not suitable for MaNGA telescope
        }

      output = list(type            = "MUSE",
                    method          = method,
                    fov             = fov,
                    aperture_shape  = "square",
                    wave_range      = c(4700.15,9351.4),
                    wave_centre     = 6975,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.25,
                    lsf_fwhm        = 2.51,
                    signal_to_noise = signal_to_noise,
                    sbin            = floor(fov/spatial_res))

    }

    if(stringr::str_to_upper(type)  == "HECTOR"){
      fov = 30
      spatial_res = 0.1
      output = list(type            = "Hector",
                    method          = method,
                    fov             = fov,
                    aperture_shape  = "hexagonal",
                    wave_range      = c(3720,5910),
                    wave_centre     = 4815,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 1.6,
                    lsf_fwhm        = 1.3,
                    signal_to_noise = signal_to_noise,
                    sbin            =  floor(fov/spatial_res))

    }

    if(stringr::str_to_upper(type)  == "CALIFA"){
      fov = 74
      spatial_res = 1
      output = list(type            = "CALIFA",
                    method          = method,
                    fov             = fov,
                    aperture_shape  = "hexagonal",
                    wave_range      = c(3700,4750),
                    wave_centre     = 4225,
                    spatial_res     = spatial_res,
                    filter          = filter,
                    wave_res        = 2.7,
                    lsf_fwhm        = 2.7,
                    signal_to_noise = signal_to_noise,
                    sbin            = floor(fov/spatial_res))

    }

    if(stringr::str_to_upper(type)  == "IFU"){

      if (missing(fov)){fov = 15} # setting to the default parameter expected
      if (missing(spatial_res)){spatial_res = 0.5}

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
  }


  return(output)
}
