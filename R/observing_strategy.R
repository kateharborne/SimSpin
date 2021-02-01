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
#' galaxy relative to the z-axis - 0 deg places the galaxy face-on, 90 deg is
#' edge-on aligned with the horizontal axis. Default is 70.
#'@param twist_deg Numeric describing the viewer's orientation relative to the
#' x-axis - 0 deg places the galaxy face-on, 90 deg places the galaxy edge-on
#' aligned with the vertical axis. Default is 0.
#'@param particle_limit Numeric describing the number of particles that must be
#' present along the line of sight for the spaxel to be occupied.
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
observing_strategy = function(z = 0.1, inc_deg = 70, twist_deg = 0, particle_limit = 1, blur = F, fwhm=2, psf="Gaussian"){


  if(blur){
    if(stringr::str_to_upper(psf) != "GAUSSIAN" & stringr::str_to_upper(psf) != "MOFFAT"){
      stop("Error: Unable to generate requested PSF. \n Please specify psf = 'Gaussian' or 'Moffat'.")
    }

    if(stringr::str_to_upper(psf) == "GAUSSIAN"){
      output = list(z              = z,
                    inc_deg        = inc_deg,
                    twist_deg      = twist_deg,
                    particle_limit = particle_limit,
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Gaussian")
    }
    if(stringr::str_to_upper(psf) == "MOFFAT"){
      output = list(z              = z,
                    inc_deg        = inc_deg,
                    twist_deg      = twist_deg,
                    particle_limit = particle_limit,
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Moffat")
    }
  } else {
    output = list(z              = z,
                  inc_deg        = inc_deg,
                  twist_deg      = twist_deg,
                  particle_limit = particle_limit,
                  blur           = F)
  }

  return(output)

}
