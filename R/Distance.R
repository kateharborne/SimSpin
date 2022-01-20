# Author: Kate Harborne
# Date: 13/01/2022
# Title: Distance class descriptor
#
#'A class to describe the projected distance to the observed object
#'
#'The purpose of this function is to generate an object that contains the
#' projected distance to the observed galaxy in a variety of units. Supply a
#' single input (i.e. distance in units of redshift "z", physical distance
#' "Mpc", or angular distance "kpc_per_arcsec") and all others will be contained
#' by the object.
#'
#'@param z A numeric describing the distance measured to the observed
#' object in units of redshift, z.
#'@param Mpc A numeric describing the distance measured to the observed
#' object in units of physical distance in mega-parsec.
#'@param kpc_per_arcsec A numeric describing the distance measured to the
#' observed object in units of angular distance in kilo-parsec per arcsecond.
#'
#'@return Returns an object of class "Distance" that summarises the
#' projected distance to the object and the associated units.
#'
#'@examples
#'Distance(z = 0.3)
#'

Distance <- function(z, Mpc, kpc_per_arcsec){

  if (!missing(z) & missing(Mpc) & missing(kpc_per_arcsec)){

    z = as.double(z)  # redshift distance
    Mpc = celestial::cosdistLumDist(z, ref="Planck") # computing distance in units of Mpc
    kpc_per_arcsec = celestial::cosdistAngScale(z, ref="Planck") # angular size given z, kpc/"

    methods::new("Distance", z = z, Mpc = Mpc, kpc_per_arcsec = kpc_per_arcsec)

  } else if (missing(z) & missing(Mpc) & !missing(kpc_per_arcsec)){

    kpc_per_arcsec = as.double(kpc_per_arcsec)
    z = .get_z_from_ang_size(kpc_per_arcsec)
    Mpc = celestial::cosdistLumDist(z, ref="Planck") # computing distance in units of Mpc

    methods::new("Distance", z = z, Mpc = Mpc, kpc_per_arcsec = kpc_per_arcsec)

  } else if (missing(z) & !missing(Mpc) & missing(kpc_per_arcsec)){

    Mpc = as.double(Mpc)
    z = .get_z_from_phy_size(Mpc)
    kpc_per_arcsec = celestial::cosdistAngScale(z, ref="Planck") # angular size given z, kpc/"

    methods::new("Distance", z = z, Mpc = Mpc, kpc_per_arcsec = kpc_per_arcsec)

  } else {
    cat("Error: Invalid distance measure specified. \n",
        "Please specify just ONE of the following parameters: z, Mpc, kpc_per_arcsec")
  }
}

.possible_z        = seq(0, 1.2, length.out = 5000)
.possible_ang_size = celestial::cosdistAngScale(.possible_z, ref="Planck") # angular size given z, kpc/"
.possible_phy_size = celestial::cosdistLumDist(.possible_z, ref="Planck")

.get_z_from_ang_size = approxfun(x = .possible_ang_size, y = .possible_z, method = "linear", rule = 1, ties = mean)
.get_z_from_phy_size = approxfun(x = .possible_phy_size, y = .possible_z, method = "linear", rule = 1, ties = mean)

# Description of the Distance class --------------------------------------------
setClass("Distance",
         slots = c(
           z = "numeric",
           Mpc = "numeric",
           kpc_per_arcsec = "numeric"
         ),
         prototype = list(
           z = NA_real_,
           Mpc = NA_real_,
           kpc_per_arcsec = NA_real_
         )
)

setValidity("Distance", function(object){

  if (!(length(object@z) == length(object@Mpc) & length(object@Mpc) == length(object@kpc_per_arcsec))) {
    "@z, @Mpc and @kpc_per_arcsec must all be equal length"
  } else if (object@z <= 0) {
    "@z cannot be negative or zero"
  } else if (object@Mpc < 0) {
    "@Mpc cannot be negative"
  } else if (object@kpc_per_arcsec <= 0) {
    "@kpc_per_arcsec cannot be negative or zero"
  } else {
    TRUE
  }

} )

setMethod("show", "Distance", function(object) {
  cat(is(object)[[1]], " to observed galaxy is:", "\n",
      " ", object@Mpc,            " Mpc ---- physical distance,", "\n",
      " ", object@kpc_per_arcsec, " kpc/'' - angular scale on the sky,", "\n",
      " z = ", object@z,          " -------- projected redshift distance,", "\n",
      sep = "")
})

# setting ability to query the slots of a distance object
setGeneric("kpc_per_arcsec", function(x) standardGeneric("kpc_per_arcsec"))
setMethod("kpc_per_arcsec", "Distance", function(x) x@kpc_per_arcsec)

setGeneric("Mpc", function(x) standardGeneric("Mpc"))
setMethod("Mpc", "Distance", function(x) x@Mpc)

setGeneric("z", function(x) standardGeneric("z"))
setMethod("z", "Distance", function(x) x@z)
