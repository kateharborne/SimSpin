# Author: Kate Harborne
# Date: 13/01/2022
# Title: Distance class descriptor
#
#'A class helper function to describe the projected distance to the observed
#' object
#'
#'The purpose of this function is to generate an object that contains the
#' projected distance to the observed galaxy in a variety of units. Supply a
#' single input (i.e. distance in units of redshift "z", physical distance
#' "Mpc", or angular distance "kpc_per_arcsec") and all others will be contained
#' by the object.
#'
#'@param z A numeric describing the distance measured to the observed
#' object in units of redshift, z. Precision to +/- 1e-7.
#'@param Mpc A numeric describing the distance measured to the observed
#' object in units of physical distance in mega-parsec. Precision to +/- 1e-3.
#'@param kpc_per_arcsec A numeric describing the distance measured to the
#' observed object in units of angular distance in kilo-parsec per arcsecond.
#' Precision to +/- 1e-5.
#'@param H0 Hubble constant. Default as Planck 2018 results where H0 = 68.4.
#'@param OmegaM Omega matter, the relative component of energy in mass. Default
#' as Planck 2018 results where OmegaM = 0.301.
#'@param OmegaL Omega lambda, the relative component of energy in dark energy.
#' Default as Planck 2018 results where OmegaL = 0.699.
#'@param OmegaR Omega radiation, the relative component of the energy in
#' radiation (including neutrinos). Default as Planck 2018 results where OmegaR
#' = 8.985075e-5.
#'@return Returns an object of class "Distance" that summarises the
#' projected distance to the object and the associated units.
#'
#'@examples
#'Distance(z = 0.3)
#'

Distance <- function(z, Mpc, kpc_per_arcsec, H0 = 68.4, OmegaM = 0.301, OmegaL = 0.699, OmegaR = 8.985075e-05){

  .possible_z        = seq(0, 1.2, length.out = 5000)
  .possible_ang_size = celestial::cosdistAngScale(.possible_z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR) # angular size given z, kpc/"
  .possible_phy_size = celestial::cosdistLumDist(.possible_z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR)

  .get_z_from_ang_size = approxfun(x = .possible_ang_size, y = .possible_z, method = "linear", rule = 1, ties = mean)
  .get_z_from_phy_size = approxfun(x = .possible_phy_size, y = .possible_z, method = "linear", rule = 1, ties = mean)

  if (!missing(z) & missing(Mpc) & missing(kpc_per_arcsec)){

    z = as.double(z)  # redshift distance
    Mpc = celestial::cosdistLumDist(z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR) # computing distance in units of Mpc
    kpc_per_arcsec = celestial::cosdistAngScale(z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR) # angular size given z, kpc/"

    methods::new("Distance", z = z, Mpc = Mpc, kpc_per_arcsec = kpc_per_arcsec)

  } else if (missing(z) & missing(Mpc) & !missing(kpc_per_arcsec)){

    kpc_per_arcsec = as.double(kpc_per_arcsec)
    z = .get_z_from_ang_size(kpc_per_arcsec)
    Mpc = celestial::cosdistLumDist(z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR) # computing distance in units of Mpc

    methods::new("Distance", z = z, Mpc = Mpc, kpc_per_arcsec = kpc_per_arcsec)

  } else if (missing(z) & !missing(Mpc) & missing(kpc_per_arcsec)){

    Mpc = as.double(Mpc)
    z = .get_z_from_phy_size(Mpc)
    kpc_per_arcsec = celestial::cosdistAngScale(z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR) # angular size given z, kpc/"

    methods::new("Distance", z = z, Mpc = Mpc, kpc_per_arcsec = kpc_per_arcsec)

  } else {
    stop("Error: Invalid distance measure specified. \n",
         "Please specify just ONE of the following parameters: z, Mpc, kpc_per_arcsec")
  }
}


# Description of the Distance class --------------------------------------------

#' An S4 class to represent projected distance to a galaxy
#'
#' @slot z Describes the projected distance to a galaxy in units of redshift.
#' @slot Mpc Describes the projected distance to a galaxy in physical units of
#'  Mega-parsecs.
#' @slot kpc_per_arcsec Describes the projected distance to a galaxy in angular
#'  units of kilo-parsecs per arcsecond.
#' @slot ref Describes the cosmological parameters which which to compute
#'  distances.

setClass("Distance",
         slots = c(
           z = "numeric",
           Mpc = "numeric",
           kpc_per_arcsec = "numeric",
           H0 = "numeric",
           OmegaM = "numeric",
           OmegaL = "numeric",
           OmegaR = "numeric"
         ),
         prototype = list(
           z = NA_real_,
           Mpc = NA_real_,
           kpc_per_arcsec = NA_real_,
           H0 = NA_real_,
           OmegaM = NA_real_,
           OmegaL = NA_real_,
           OmegaR = NA_real_
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
  cat(methods::is(object)[[1]], " to observed galaxy is:", "\n",
      " ", object@Mpc,            " Mpc ---- physical distance,", "\n",
      " ", object@kpc_per_arcsec, " kpc/'' - angular scale on the sky,", "\n",
      " z = ", object@z,          " -------- projected redshift distance,", "\n",
      sep = "")
})

# setting ability to query the slots of a distance object
#' @describeIn Distance Query the kpc_per_arcsec of a Distance object
#' @param x Distance object
setGeneric("kpc_per_arcsec", function(x) standardGeneric("kpc_per_arcsec"))

#' Query the kpc_per_arcsec of a Distance object
#' @param x Distance object
#' @return The distance to the projected galaxy in units of kilo-parsec per
#'  arcsecond.
setMethod("kpc_per_arcsec", "Distance", function(x) x@kpc_per_arcsec)

#' @describeIn Distance Query the Mpc of a Distance object
#' @param x Distance object
setGeneric("Mpc", function(x) standardGeneric("Mpc"))

#' Query the Mpc of a Distance object
#' @inheritParams kpc_per_arcsec
#' @return The distance to the projected galaxy in units of Mpc.
setMethod("Mpc", "Distance", function(x) x@Mpc)

#' @describeIn Distance Query the z of a Distance object
#' @param x Distance object
setGeneric("z", function(x) standardGeneric("z"))

#' Distance Query the z of a Distance object
#' @inheritParams kpc_per_arcsec
#' @return The distance to the projected galaxy in units of redshift.
setMethod("z", "Distance", function(x) x@z)
