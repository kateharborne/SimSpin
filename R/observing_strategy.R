# Author: Kate Harborne
# Date: 28/10/2020
# Title: observing_strategy - a class for describing the obsering conditions
#
#'A class to describe the basic properties of the observing conditions
#'
#'The purpose of this function is to generate a class that describes the
#'conditions of the observation for \code{build_datacube()}.
#'
#'@param z Numeric describing the redshift to which the observed galaxy is
#' projected. Default is 0.1.
#'@param inc_deg Numeric describing the projected inclination of the observed
#' galaxy - 0 deg places the galaxy face-on, 90 deg is edge-on. Default is 70.
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
observing_strategy = function(z = 0.1, inc_deg = 70, blur = F, fwhm=2, psf="Gaussian"){


  if(blur){
    if(stringr::str_to_upper(psf) != "GAUSSIAN" & stringr::str_to_upper(psf) != "MOFFAT"){
      stop("Error: Unable to generate requested PSF. \n Please specify psf = 'Gaussian' or 'Moffat'.")
    }

    if(stringr::str_to_upper(psf) == "GAUSSIAN"){
      output = list(z       = z,
                    inc_deg = inc_deg,
                    blur    = T,
                    fwhm    = fwhm,
                    psf     = "Gaussian")
    }
    if(stringr::str_to_upper(psf) == "MOFFAT"){
      output = list(z       = z,
                    inc_deg = inc_deg,
                    blur    = T,
                    fwhm    = fwhm,
                    psf     = "Moffat")
    }
  } else {
    output = list(z       = z,
                  inc_deg = inc_deg,
                  blur    = F)
  }

  return(output)

}