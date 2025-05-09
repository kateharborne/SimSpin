# Author: Kate Harborne
# Date: 13/01/2022
# Title: Pointing Class descriptors
#
#'A class to describe the center pointing of the telescope relative to the
#' centre of the observed object.
#
#'The purpose of this function is to generate an object that contains the
#' pointing of the centre of the observation relative to the centre of the
#' observed galaxy in a variety of units. Supply a single pair of inputs (i.e.
#' pointing in units of physical distance, "x_kpc" and "y_kpc" OR angular
#' distance, "x_deg" and "y_deg") and all others will be contained
#' by the object.
#'
#'@param xy_deg Offset from centre (0,0) in angular units of degrees along the
#' x-axis and y-axis respectively.
#'@param xy_kpc Offset from centre (0,0) in  physical units of kilo-parsec along
#' the x-axis and y-axis respectively.
#'@param distance The Distance object used to describe the projected distance
#' of the observed galaxy.
#'
#'@return Returns an object of class "Pointing" that summarises the
#' pointing position of the telescope relative to the centre of the object
#' being observed in all associated units.
#'
#'@examples
#'dist = Distance(z=0.3)
#'Pointing(xy_kpc = c(-1,0), distance = dist)
#'

Pointing <- function(xy_deg, xy_kpc, distance){

  if (missing(xy_deg) & missing(xy_kpc) & missing(distance)){
    stop("Error: Invalid pointing measure specified. \n",
         "Please specify inputs 'xy_deg' and 'distance' OR 'xy_kpc' & 'distance'. ")

  } else if (!missing(xy_deg) & missing(xy_kpc) & !missing(distance)){
    xy_deg = as.double(xy_deg)

    kpc_per_degree = kpc_per_arcsec(distance) * 3600

    xy_kpc = xy_deg * kpc_per_degree

    methods::new("Pointing", xy_deg = xy_deg, xy_kpc = xy_kpc)

  } else if (missing(xy_deg) & !missing(xy_kpc) & !missing(distance)) {

    xy_kpc = as.double(xy_kpc)

    kpc_per_degree = kpc_per_arcsec(distance) * 3600

    xy_deg = xy_kpc / kpc_per_degree

    methods::new("Pointing", xy_deg = xy_deg, xy_kpc = xy_kpc)

  } else if (missing(distance)) {
    stop("Error: Distance class object required to specify scale. \n",
         "Please specify 'distance' input. ")

  } else {
    stop("Error: Invalid pointing measure specified. \n",
         "Please specify just ONE of the following parameters: 'xy_deg' OR 'xy_kpc' ")
  }
}

# Description of the Pointing class --------------------------------------------

#' An S4 class to represent the telescope pointing relative to the galaxy centre
#'
#' @slot xy_deg Numeric of length 2. Describes the pointing relative to the
#'  galaxy centre in units of degrees.
#' @slot xy_kpc Numeric of length 2. Describes the pointing relative to the
#'  galaxy centre in units of kilo-parsecs.

setClass("Pointing",
         slots = c(
           xy_deg = "numeric",
           xy_kpc = "numeric"
         ),
         prototype = list(
           xy_deg = NA_real_,
           xy_kpc = NA_real_
         )
)

setValidity("Pointing", function(object){

  if (length(object@xy_deg) != 2) {
    "@xy_deg must be a numeric input of length 2"
  } else if (length(object@xy_kpc) != 2) {
    "@xy_kpc must be a numeric input of length 2"
  } else {
    TRUE
  }

} )

setMethod("show", "Pointing", function(object) {
  cat(methods::is(object)[[1]], " relative to observed galaxy centre at (0,0) is:", "\n",
      " (x,y):  (", object@xy_deg[1], ",", object@xy_deg[2], ") degrees -- ", "\n",
      " (x,y):  (", object@xy_kpc[1], ",", object@xy_kpc[2], ") kpc ------ ", "\n",
      sep = "")
})

# setting ability to query the slots of a pointing object
#' @describeIn Pointing Query the xy_deg of a Pointing object
#' @param x Pointing object
setGeneric("xy_deg", function(x) standardGeneric("xy_deg"))

#' Query the xy_deg of a Pointing object
#' @param x Pointing object
#' @return The (x,y) offset (in units of degrees) of the telescope pointing
#'  relative to the galaxy centre.
setMethod("xy_deg", "Pointing", function(x) x@xy_deg)

#' @describeIn Pointing Query the xy_kpc of a Pointing object
#' @param x Pointing object
setGeneric("xy_kpc", function(x) standardGeneric("xy_kpc"))

#' Query the xy_kpc of a Pointing object
#' @param x Pointing object
#' @return The (x,y) offset (in units of kilo-parsec) of the telescope pointing
#'  relative to the galaxy centre.
setMethod("xy_kpc", "Pointing", function(x) x@xy_kpc)
