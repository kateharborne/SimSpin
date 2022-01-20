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
#'@param x_deg Offset from centre (0,0) in angular units of degrees along the
#' x-axis.
#'@param y_deg Offset from centre (0,0) in angular units of degrees along the
#' y-axis.
#'@param x_kpc Offset from centre (0,0) in  physical units of kilo-parsec along
#' the x-axis.
#'@param y_kpc Offset from centre (0,0) in  physical units of kilo-parsec along
#' the y-axis.
#'@param distance The Distance object used to describe the projected distance
#' of the observed galaxy.
#'
#'@return Returns an object of class "Pointing" that summarises the
#' pointing position of the telescope relative to the centre of the object
#' being observed in all associated units.
#'
#'@examples
#'Pointing(x_kpc = -1.0, y_kpc = 0)
#'

Pointing <- function(x_deg, y_deg, x_kpc, y_kpc, distance){

  if (!missing(x_deg) & !missing(y_deg) & missing(x_kpc) & missing(y_kpc)){
    x_deg = as.double(x_deg)
    y_deg = as.double(y_deg)
    kpc_per_degree = kpc_per_arcsec(distance) * 3600

    x_kpc = x_deg * kpc_per_degree
    y_kpc = y_deg * kpc_per_degree

    new("Pointing", x_deg = x_deg, y_deg = y_deg, x_kpc = x_kpc, y_kpc = y_kpc)

  } else if (missing(x_deg) & missing(y_deg) & !missing(x_kpc) & !missing(y_kpc)) {

    x_kpc = as.double(x_kpc)
    y_kpc = as.double(y_kpc)
    kpc_per_degree = kpc_per_arcsec(distance) * 3600

    x_deg = x_kpc / kpc_per_degree
    y_deg = y_kpc / kpc_per_degree

    new("Pointing", x_deg = x_deg, y_deg = y_deg, x_kpc = x_kpc, y_kpc = y_kpc)

  } else {
    cat("Error: Invalid pointing measure specified. \n",
        "Please specify just ONE PAIR of the following parameters: 'x_deg', 'y_deg' OR 'x_kpc', 'y_kpc' ")
  }
}

# Description of the Pointing class --------------------------------------------
setClass("Pointing",
         slots = c(
           x_deg = "numeric",
           y_deg = "numeric",
           x_kpc = "numeric",
           y_kpc = "numeric"
         ),
         prototype = list(
           x_deg = NA_real_,
           y_deg = NA_real_,
           x_kpc = NA_real_,
           y_kpc = NA_real_
         )
)

setValidity("Pointing", function(object){

  if (length(object@x_deg) > 1) {
    "@x_deg must be a single numeric input"
  } else if (length(object@y_deg) > 1) {
    "@y_deg must be a single numeric input"
  } else if (length(object@x_kpc) > 1) {
    "@x_kpc must be a single numeric input"
  } else if (length(object@y_kpc) > 1) {
    "@y_kpc must be a single numeric input"
  } else {
    TRUE
  }

} )

setMethod("show", "Pointing", function(object) {
  cat(is(object)[[1]], " relative to observed galaxy centre at (0,0) is:", "\n",
      "    x_deg:", object@x_deg,  " degrees ---- ", "\n",
      "   y_deg:", object@y_deg, " degrees ---- ", "\n",
      " (x,y):  (", object@x_kpc, ",", object@y_kpc, ") kpc ---- ", "\n",
      sep = "")
})

# setting ability to query the slots of a pointing object
setGeneric("x_deg", function(x) standardGeneric("x_deg"))
setMethod("x_deg", "Pointing", function(x) x@x_deg)

setGeneric("y_deg", function(x) standardGeneric("y_deg"))
setMethod("y_deg", "Pointing", function(x) x@y_deg)

setGeneric("x_kpc", function(x) standardGeneric("x_kpc"))
setMethod("x_kpc", "Pointing", function(x) x@x_kpc)

setGeneric("y_kpc", function(x) standardGeneric("y_kpc"))
setMethod("y_kpc", "Pointing", function(x) x@y_kpc)
