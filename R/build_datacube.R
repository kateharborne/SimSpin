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
#'@param align Boolean indicating whether or not to align the semi-major axis with the x-axis.
#'@param blur \emph{Optional} Specified to apply observational seeing effects to the cube. A list
#' of the form \code{list("psf" = "Moffat", "fwhm" = 0.5)}. \code{"psf"} specifies the shape of the
#' PSF chosen and may be either \code{"Moffat"} or \code{"Gaussian"}. \code{"fwhm"} is a numeric
#' specifying the full-width half-maximum of the PSF given in units of arcseconds.
#'@param multi_thread A boolean specifying whether you would like to multi-thread the process.
#'@return A list containing:
#' \item{\code{$datacube}}{A 3D array corresponding to the kinematic data cube.}
#' \item{\code{$xbin_labels}}{Bin labels for the x-spatial dimension.}
#' \item{\code{$zbin_labels}}{Bin labels for the z-spatial dimension.}
#' \item{\code{$vbin_labels}}{Bin labels for the velocity dimension.}
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
#'                       inc_deg      = 0,
#'                       multi_thread = FALSE)
#'

build_datacube = function(simdata, r200 = 200, z, fov, ap_shape, central_wvl,
                          lsf_fwhm, pixel_sscale, pixel_vscale, inc_deg, filter="g",
                          blur, align=FALSE, multi_thread=TRUE){

  if (missing(blur)) {                                     # IF spatial blurring IS NOT requested

    observe_data = obs_data_prep(simdata = simdata, r200 = r200, z = z, fov = fov, ap_shape = ap_shape,
                                 central_wvl = central_wvl, lsf_fwhm = lsf_fwhm, pixel_sscale = pixel_sscale,
                                 pixel_vscale = pixel_vscale, inc_deg = inc_deg, align = align) # prep simulation data in observer units
    fluxes = flux_grid(obs_data = observe_data, filter = filter, multi_thread = multi_thread)
    ifu_imgs = ifu_cube(obs_data = observe_data, flux_data = fluxes) # construct IFU data cube

    output       = list("datacube"     = ifu_imgs$cube,
                        "xbin_labels"  = ifu_imgs$xbin_labels,
                        "zbin_labels"  = ifu_imgs$zbin_labels,
                        "vbin_labels"  = ifu_imgs$vbin_labels)

    return(output)

  } else {                                                 # IF spatial blurring IS requested

    observe_data = obs_data_prep(simdata = simdata, r200 = r200, z = z, fov = fov, ap_shape = ap_shape,
                                 central_wvl = central_wvl, lsf_fwhm = lsf_fwhm, pixel_sscale = pixel_sscale,
                                 pixel_vscale = pixel_vscale, inc_deg = inc_deg, align = align) # prep simulation data in observer units
    fluxes = flux_grid(obs_data = observe_data, filter = filter, multi_thread = multi_thread)
    ifu_imgs = ifu_cube(obs_data = observe_data, flux_data = fluxes) # construct IFU data cube
    blur_imgs = blur_cube(obs_data = observe_data, ifu_datacube = ifu_imgs, psf = blur$psf,
                             fwhm = blur$fwhm)
                                                           # blur IFU data cube

    output       = list("datacube"     = blur_imgs$cube,
                        "xbin_labels"  = ifu_imgs$xbin_labels,
                        "zbin_labels"  = ifu_imgs$zbin_labels,
                        "vbin_labels"  = ifu_imgs$vbin_labels)
    return(output)

  }

}
