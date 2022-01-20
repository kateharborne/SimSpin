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
#'@param distance An object of the class Distance. Includes the measured
#' distance to the object in a variety of possible units. See
#' \link[SimSpin::Distance]{Distance()} for further details.
#'@param inc_deg Numeric describing the projected inclination of the observed
#' galaxy relative to the z-axis - 0 deg places the galaxy face-on, 90 deg is
#' edge-on aligned with the horizontal axis. Default is 70.
#'@param twist_deg Numeric describing the viewer's orientation relative to the
#' x-axis - 0 deg places the galaxy face-on, 90 deg places the galaxy edge-on
#' aligned with the vertical axis. Default is 0.
#'@param pointing Numeric array `c(x,y)`. Two elements specifying the position
#' at which the observation is centred given as a shift relative to the centre.
#' Can be specified in either physical kilo-parsec offsets
#' (`pointing_unit = 'kpc'`) or in an angular reference frame using degrees
#' (`pointing_unit = 'deg'`)
#'@param pointing_unit Character describing the units of the pointing. Options
#' include either physical kilo-parsec offsets `kpc` or an angular reference
#' frame using angular offsets in degrees `deg`. All offsets are defined
#' relative to the centre of the galaxy. Input is NOT case sensitive.
#'@param blur Boolean describing whether seeing conditions should be applied.
#' Default is FALSE.
#'@param fwhm If \code{blur = TRUE}, a numeric that describes the full-width
#' half-maximum of the point spread function (PSF) used to blur the image,
#' given in arcsec.
#'@param psf A string to describe the shape of the PSF. Options include
#' "Gaussian" or "Moffat". Input is NOT case sensitive.
#'@return Returns an object of class "objective" that describes the
#' conditions in which the observation is made. Required to run
#' \code{build_datacube()}.
#'@examples
#'conditions = objective()
#'

objective = function(distance = Distance(z=0.3),
                     inc_deg = 70, twist_deg = 0,
                     pointing = c(0,0), pointing_unit = "kpc",
                     blur = F, fwhm=2, psf="Gaussian"){

  if (stringr::str_to_upper(pointing_unit) != "KPC" & stringr::str_to_upper(pointing_unit) != "DEG"){
    stop("Error: Unable to generate requested Pointing. \n
         Please specify ONE of the following: pointing_unit = 'kpc' OR 'deg'.")
  }

  if (blur){
    if(stringr::str_to_upper(psf) != "GAUSSIAN" & stringr::str_to_upper(psf) != "MOFFAT"){
      stop("Error: Unable to generate requested PSF. \n
           Please specify ONE of the following: psf = 'Gaussian' OR 'Moffat'.")
    }

    if (stringr::str_to_upper(psf) == "GAUSSIAN"){
      output = list(distance       = distance,
                    inc_deg        = inc_deg,
                    twist_deg      = twist_deg,
                    pointing       = if (stringr::str_to_upper(pointing_unit) == "KPC") {Pointing(xy_kpc = pointing, distance = distance)} else {Pointing(xy_deg = pointing, distance = distance)},
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Gaussian")
    }
    if (stringr::str_to_upper(psf) == "MOFFAT"){
      output = list(distance       = distance,
                    inc_deg        = inc_deg,
                    twist_deg      = twist_deg,
                    pointing       = if (stringr::str_to_upper(pointing_unit) == "KPC") {Pointing(xy_kpc = pointing, distance = distance)} else {Pointing(xy_deg = pointing, distance = distance)},
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Moffat")
    }
  } else {
    output = list(distance       = distance,
                  inc_deg        = inc_deg,
                  twist_deg      = twist_deg,
                  pointing       = if(stringr::str_to_upper(pointing_unit) == "KPC"){Pointing(xy_kpc = pointing, distance = distance)}else{Pointing(xy_deg = pointing, distance = distance)},
                  blur           = F)
  }

  return(output)


}
