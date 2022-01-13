# Author: Kate Harborne
# Date: 28/10/2020
# Title: object_properties - a class for describing the properties of the object
# being observed
#
#'A class to describe the basic properties of the object being observed
#'
#'The purpose of this function is to generate a class that describes the
#'conditions of the observation for \code{build_datacube()}.
#'
#'@param distance Numeric describing the distance to which the observed galaxy
#' is projected. Default is 0.1 in units of redshift "z". Avaliable unit
#' options include `distance_unit = 'z'` (redshift distance),
#' `distance_unit = 'Mpc'`  (physical distance in units of mega-parsec) and
#' `distance_unit = 'kpc/pix'` (physical size of each pixel in kpc i.e.
#' kilo-parsec).
#'@param distance_units Character describing the type of units used for the
#' distance parameter. Default is `z` for distances given in units of redshift.
#' Avaliable options include `z` (redshift distance), `Mpc` (physical distance
#' in units of mega-parsec) and `kpc/pix` (physical size of each pixel in kpc
#' i.e. kilo-parsec). Input is NOT case sensitive.
#'@param inc_deg Numeric describing the projected inclination of the observed
#' galaxy relative to the z-axis - 0 deg places the galaxy face-on, 90 deg is
#' edge-on aligned with the horizontal axis. Default is 70.
#'@param twist_deg Numeric describing the viewer's orientation relative to the
#' x-axis - 0 deg places the galaxy face-on, 90 deg places the galaxy edge-on
#' aligned with the vertical axis. Default is 0.
#'@param pointing Numeric array `c(x,y)`. Two elements specifying the position
#' at which the observation is centred given as a shift relative to the centre.
#' Can be specified in either physical kilo-parsec offsets
#' (`pointing_unit = 'kpc'`) or in an angular reference frame using RA and DEC
#' (`pointing_unit = 'RA/DEC'`)
#'@param pointing_unit Character describing the units of the pointing. Options
#' include either physical kilo-parsec offsets `kpc` or an angular reference
#' frame using RA and DEC `RA/DEC`. All offsets are defined relative to the
#' centre of the galaxy. Input is NOT case sensitive.
#'@param blur Boolean describing whether seeing conditions should be applied.
#' Default is FALSE.
#'@param fwhm If \code{blur = TRUE}, a numeric that describes the full-width
#' half-maximum of the point spread function (PSF) used to blur the image,
#' given in arcsec.
#'@param psf A string to describe the shape of the PSF. Options include
#' "Gaussian" or "Moffat". Input is NOT case sensitive.
#'@return Returns an object of class "observing_strategy" that describes the
#' conditions in which the observation is made. Required to run
#' \code{build_datacube()}.
#'@examples
#'conditions = observing_strategy()
#'

objective = function(distance = 0.1, distance_unit = "z",
                     inc_deg = 70, twist_deg = 0,
                     pointing = c(0,0), pointing_unit ="kpc",
                     blur = F, fwhm=2, psf="Gaussian"){

  # checks for validity
  if (stringr::str_to_upper(distance_unit) != "Z" & stringr::str_to_upper(distance_unit) != "MPC" & stringr::str_to_upper(distance_unit) != "KPC/PIX"){
    stop("Error: Distance unit is not supported. \n Please specify distance_unit = 'z' OR 'Mpc' OR 'kpc/pix'.")
  }
  if (stringr::str_to_upper(distance_unit) == "Z" & distance <= 0){
    stop("Error: Unable to process supplied redshift, z <= 0. \n Please specify z > 0 and try again.")
  }
  if (stringr::str_to_upper(distance_unit) == "MPC" & distance < 0){
    stop("Error: Unable to process supplied distance to system < 0 Mpc. \n Please specify a distance >= 0 Mpc and try again.")
  }
  if (stringr::str_to_upper(distance_unit) == "KPC/PIX" & distance <= 0){
    stop("Error: Unable to process supplied pixel scale =< 0 kpc/pix. \n Please specify a pixel scale > 0 kpc/pix and try again.")
  }
  if (stringr::str_to_upper(pointing_unit) != "KPC" & stringr::str_to_upper(pointing_unit) != "RA/DEC"){
    stop("Error: Pointing unit is not supported. \n Please specify pointing_unit = 'kpc' OR 'RA/DEC'.")
  }


  if (blur){
    if(stringr::str_to_upper(psf) != "GAUSSIAN" & stringr::str_to_upper(psf) != "MOFFAT"){
      stop("Error: Unable to generate requested PSF. \n Please specify psf = 'Gaussian' or 'Moffat'.")
    }

    if (stringr::str_to_upper(psf) == "GAUSSIAN"){
      output = list(distance       = distance,
                    distance_unit  = distance_unit,
                    inc_deg        = inc_deg,
                    twist_deg      = twist_deg,
                    pointing       = pointing,
                    pointing_unit  = pointing_unit,
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Gaussian")
    }
    if (stringr::str_to_upper(psf) == "MOFFAT"){
      output = list(distance       = distance,
                    distance_unit  = distance_unit,
                    inc_deg        = inc_deg,
                    twist_deg      = twist_deg,
                    pointing       = pointing,
                    pointing_unit  = pointing_unit,
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Moffat")
    }
  } else {
    output = list(distance       = distance,
                  distance_unit  = distance_unit,
                  inc_deg        = inc_deg,
                  twist_deg      = twist_deg,
                  pointing       = pointing,
                  pointing_unit  = pointing_unit,
                  blur           = F)
  }

  return(output)


}
