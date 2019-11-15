# Kate Harborne (last edit - 15/11/2019)
#'Constructing kinematic data cubes.
#'
#'The purpose of this basic function is to use the \code{SimSpin} package sub-functions to
#' construct an IFU kinematic data cube. This function will call each required sub-function
#' (\code{\link{obs_data_prep}}, \code{\link{flux_grid}}, \code{\link{ifu_cube}},
#' \code{\link{blur_cube}}) and produce a 3D array corresponding to the kinematic data cube
#' for that mock observation.
#'
#'@param simdata The simulation information data.frame output by \code{\link{sim_data}}.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param z The galaxy redshift.
#'@param fov The field of view of the IFU, diameter in arcseconds.
#'@param ap_shape The shape of the field of view, with options "circular", "square" or "hexagonal".
#'@param central_wvl The central filter wavelength used for the observation, given in angstroms.
#'@param lsf_fwhm The line spread function full-width half-max, given in angstroms.
#'@param pixel_sscale The corresponding spatial pixel scale associated with a given telescope
#' output in arcseconds.
#'@param pixel_vscale The corresponding velocity pixel scale associated with a given telescope
#' filter output in angstroms.
#'@param inc_deg The inclination at which to observe the galaxy in degrees.
#'@param filter If Age/Metallicity is supplied, the filter within which the SED is generated.
#'Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
#'@param blur \emph{Optional} Specified to apply observational seeing effects to the cube. A list
#' of the form \code{list("psf" = "Moffat", "fwhm" = 0.5)}. \code{"psf"} specifies the shape of the
#' PSF chosen and may be either \code{"Moffat"} or \code{"Gaussian"}. \code{"fwhm"} is a numeric
#' specifying the full-width half-maximum of the PSF given in units of arcseconds.
#'@return A list containing:
#' \item{\code{$datacube}}{A 3D array corresponding to the kinematic data cube.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' cube = build_datacube(simdata      = galaxy_data,
#'                       z            = 0.1,
#'                       fov          = 15,
#'                       ap_shape     = "circular",
#'                       central_wvl  = 4800,
#'                       lsf_fwhm     = 2.65,
#'                       pixel_sscale = 0.5,
#'                       pixel_vscale = 1.04,
#'                       inc_deg      = 0)
#'

build_datacube = function(simdata, r200 = 200, z, fov, ap_shape, central_wvl,
                          lsf_fwhm, pixel_sscale, pixel_vscale, inc_deg, filter="g",
                          blur){

  if (missing(blur)) {                                     # IF spatial blurring IS NOT requested

    observe_data = obs_data_prep(simdata, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg) # prep simulation data in observer units
    fluxes = flux_grid(obs_data = observe_data, filter = filter)
    ifu_imgs = ifu_cube(obs_data = observe_data, flux_data = fluxes) # construct IFU data cube

    output       = list("datacube"     = ifu_imgs$cube)

    return(output)

  } else {                                                 # IF spatial blurring IS requested

    observe_data = obs_data_prep(simdata, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg) # prep simulation data in observer units
    fluxes = flux_grid(obs_data = observe_data, filter = filter)
    ifu_imgs = ifu_cube(obs_data = observe_data, flux_data = fluxes) # construct IFU data cube
    blur_imgs = blur_cube(obs_data = observe_data, ifu_datacube = ifu_imgs, psf = blur$psf,
                             fwhm = blur$fwhm)
                                                           # blur IFU data cube

    output       = list("datacube"     = blur_imgs$cube)
    return(output)

  }

}
