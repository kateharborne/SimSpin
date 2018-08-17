# Kate Harborne (last edit - 23/04/2018)
#'Constructing kinematic data cubes.
#'
#'The purpose of this basic function is to use the \code{SimSpin} package sub-functions to
#' construct an IFU kinematic data cube. This function will call each required sub-function
#' (\code{\link{obs_data_prep}}, \code{\link{ifu_cube}}, \code{\link{blur_cube}}) and produce a 3D
#' array corresponding to the kinematic data cube for that mock observation, the corresponding axes
#' labels and the observed axis ratio of the observed galaxy.
#'
#'@param filename The Gadget file containing the particle information of the galaxy to be analysed.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the
#' simulation, 1 - gas, 2 - dark matter, 3 - disc, 4 - bulge, 5 - stars, 6 - boundary.
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
#'@param m2l_disc The mass-to-light ratio of the disc component in solar units.
#'@param m2l_bulge The mass-to-light ratio of the bulge component in solar units.
#'@param threshold The magnitude limit of the observation.
#'@param blur \emph{Optional} Specified to apply observational seeing effects to the cube. A list
#' of the form \code{list("psf" = "Moffat", "fwhm" = 0.5)}. \code{"psf"} specifies the shape of the
#' PSF chosen and may be either \code{"Moffat"} or \code{"Gaussian"}. \code{"fwhm"} is a numeric
#' specifying the full-width half-maximum of the PSF given in units of arcseconds.
#'@return A list containing:
#'\item{\code{$datacube}}{A 3D array corresponding to the kinematic data cube.}
#'\item{\code{$xbin_labels}}{Bin labels for the x-spatial dimension.}
#'\item{\code{$ybin_labels}}{Bin labels for the y-spatial dimension.}
#'\item{\code{$vbin_labels}}{Bin labels for the velocity dimension.}
#'\item{\code{$axis_ratio}}{The axis ratio of the observed galaxy in the form of a data frame where
#' \code{$a} is the semi-major axis and \code{$b} is the semi-minor axis given in kpc.}
#'@examples
#' cube = build_datacube(filename     = system.file("extdata", 'S0_vignette', package="SimSpin"),
#'                       r200         = 200,
#'                       z            = 0.1,
#'                       fov          = 15,
#'                       ap_shape     = "circular",
#'                       central_wvl  = 4800,
#'                       lsf_fwhm     = 2.65,
#'                       pixel_sscale = 0.5,
#'                       pixel_vscale = 1.04,
#'                       inc_deg      = 0,
#'                       m2l_disc     = 2,
#'                       m2l_bulge    = 1,
#'                       threshold    = 25)
#'

build_datacube = function(filename, ptype = NA, r200 = 200, z, fov, ap_shape, central_wvl,
                          lsf_fwhm, pixel_sscale, pixel_vscale, inc_deg, m2l_disc, m2l_bulge,
                          threshold, blur){

  if (missing(blur)) {                                     # IF spatial blurring IS NOT requested

    observe_data = obs_data_prep(filename, ptype, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg, m2l_disc, m2l_bulge)
                                                           # prep simulation data in observer units
    ifu_imgs     = ifu_cube(observe_data, threshold)       # construct IFU data cube
    output       = list("datacube"     = ifu_imgs$cube,
                        "xbin_labels"  = ifu_imgs$xbin_labels,
                        "ybin_labels"  = ifu_imgs$ybin_labels,
                        "vbin_labels"  = ifu_imgs$vbin_labels,
                        "axis_ratio"   = ifu_imgs$axis_ratio,
                        "angular_size" = observe_data$angular_size)
    return(output)

  } else {                                                 # IF spatial blurring IS requested

    observe_data = obs_data_prep(filename, ptype, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg, m2l_disc, m2l_bulge)
                                                           # prep simulation data in observer units
    ifu_imgs     = ifu_cube(observe_data, threshold)       # construct IFU data cube
    blur_imgs    = blur_cube(ifu_imgs, sbinsize = observe_data$sbinsize, psf = blur$psf,
                             fwhm = blur$fwhm, angular_size = observe_data$angular_size)
                                                           # blur IFU data cube
    output       = list("datacube"     = blur_imgs$cube,
                        "xbin_labels"  = blur_imgs$xbin_labels,
                        "ybin_labels"  = blur_imgs$ybin_labels,
                        "vbin_labels"  = blur_imgs$vbin_labels,
                        "axis_ratio"   = ifu_imgs$axis_ratio,
                        "angular_size" = observe_data$angular_size)
    return(output)

  }

}