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
#'@param dist_z The projected distance at which the galaxy is placed in terms of
#' redshift distance. This parameter is one of three optional ways of describing
#' the distance to the observed galaxy. Only ONE of these three must be
#' specified.
#'@param dist_Mpc The projected distance at which the galaxy is placed in terms
#' of physical distance in units of Mega-parsec. This parameter is the second of
#' three optional ways of describing the distance to the observed galaxy. Only
#' ONE of these three must be specified.
#'@param dist_kpc_per_arcsec The projected distance at which the galaxy is
#' placed in terms of angular distance in units of kilo-parsec per arcsecond.
#' This parameter is the third of three optional ways of describing the distance
#' to the observed galaxy. Only ONE of these three must be specified.
#'@param inc_deg Numeric describing the projected inclination of the observed
#' galaxy relative to the z-axis - 0 deg places the galaxy face-on, 90 deg is
#' edge-on aligned with the horizontal axis. Default is 70.
#'@param twist_deg Numeric describing the viewer's orientation relative to the
#' x-axis - 0 deg places the galaxy face-on, 90 deg places the galaxy edge-on
#' aligned with the vertical axis. Default is 0.
#'@param pointing_kpc Numeric array `c(x,y)`. Two elements specifying the
#' position at which the observation is centred given as a shift relative to the
#' centre at (0,0) in units of kilo-parsecs. This parameter is ONE of two
#' optional ways of describing this offset. Default is c(0,0) offset.
#'@param pointing_deg Numeric array `c(x,y)`. Two elements specifying the
#' position at which the observation is centred given as a shift relative to the
#' centre at (0,0) in units of degrees. This parameter is ONE of two
#' optional ways of describing this offset. Default is c(0,0) offset.
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

objective = function(dist_z = 0.05, dist_Mpc, dist_kpc_per_arcsec,
                     inc_deg = 70, twist_deg = 0,
                     pointing_kpc = c(0,0), pointing_deg,
                     blur = F, fwhm=2, psf="Gaussian"){

  if (missing(dist_Mpc) & missing(dist_kpc_per_arcsec)){
    distance = Distance(z = dist_z)
  } else if (missing(dist_z) & !missing(dist_Mpc) & missing(dist_kpc_per_arcsec)){
    distance = Distance(Mpc = dist_Mpc)
  } else if (missing(dist_z) & missing(dist_Mpc) & !missing(dist_kpc_per_arcsec)){
    distance = Distance(kpc_per_arcsec = dist_kpc_per_arcsec)
  } else {
    stop("Error: Please specify ONE method of projected distance to the observed galaxy. \n
         Please specify ONE of the following input parameters: 'dist_z' OR 'dist_Mpc' OR 'dist_kpc_per_arcsec'. ")
  }

  if (missing(pointing_deg)){
    pointing = Pointing(xy_kpc = pointing_kpc, distance = distance)
  } else if (missing(pointing_kpc) & !missing(pointing_deg)){
    pointing = Pointing(xy_deg = pointing_deg, distance = distance)
  } else {
    stop("Error: Unable to generate requested Pointing. \n
          Please specify ONE of the following: pointing_kpc OR pointing_deg.")
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
                    pointing       = pointing,
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Gaussian")
    }
    if (stringr::str_to_upper(psf) == "MOFFAT"){
      output = list(distance       = distance,
                    inc_deg        = inc_deg,
                    twist_deg      = twist_deg,
                    pointing       = pointing,
                    blur           = T,
                    fwhm           = fwhm,
                    psf            = "Moffat")
    }
  } else {
    output = list(distance       = distance,
                  inc_deg        = inc_deg,
                  twist_deg      = twist_deg,
                  pointing       = pointing,
                  blur           = F)
  }

  return(output)


}
