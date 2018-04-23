# Kate Harborne (last edit - 23/04/2018)
#'Measuring observable galaxy kinematics.
#'
#'The purpose of this basic function is to use the \code{SimSpin} package sub-functions to measure
#' the observable ratio V/\eqn{\sigma} of a simulated galaxy model. This function will call
#' each required sub-function (\code{\link{obs_data_prep}}, \code{\link{ifu_cube}},
#' \code{\link{blur_cube}}, \code{\link{find_reff}}, \code{\link{obs_vsigma}} and
#' \code{\link{plot_ifu}}) and return the V/\eqn{\sigma} ratio within a user specified measurement
#' radius.
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
#'@param threshold The flux threshold of the observation.
#'@param measure_type A list specifying the radius within which \eqn{\lambda_R} is measured. If
#' \code{list(type = "fit", fac = 1)}, \eqn{\lambda_R} is measured in a specified multiple of
#' \eqn{R_{eff}}, where \eqn{R_{eff}} has been calculated from the unblurred galaxy  counts image.
#' If \code{list(type = "specified", fac = 1, axis_ratio = data.frame("a" = 2, "b" = 1))},
#' \eqn{\lambda_R} is calculated within some factor of the effective radius, \code{fac *}
#' \eqn{R_{eff}} with specified axis ratio externally as a data frame in units of kpc.
#'@param blur \emph{Optional} Specify if you wish to apply observational seeing effects to the
#' cube. A list of the form \code{list("psf" = "Moffat", "fwhm" = 0.5)}. \code{"psf"} specifies the
#' shape of the PSF chosen and may be either \code{"Moffat"} or \code{"Gaussian"}. \code{"fwhm"} is
#' a numeric specifying the full-width half-maximum of the PSF given in units of arcseconds.
#'@return A list containing:
#' \item{\code{$datacube}}{A 3D array corresponding to the kinematic data cube.}
#' \item{\code{$xbin_labels}}{Bin labels for the x-spatial dimension.}
#' \item{\code{$ybin_labels}}{Bin labels for the y-spatial dimension.}
#' \item{\code{$vbin_labels}}{Bin labels for the velocity dimension.}
#' \item{\code{$axis_ratio}}{The axis ratio of the observed galaxy in the form of a data frame
#'  where \code{$a} is the semi-major axis and \code{$b} is the semi-minor axis given in kpc.}
#' \item{\code{$vsigma}}{The observed V/\eqn{\sigma}.}
#' \item{\code{$counts_img}}{The observed flux image.}
#' \item{\code{$velocity_img}}{The observed line-of-sight velocity image.}
#' \item{\code{$dispersion_img}}{The observed line-of-sight velocity dispersion image.}
#' \item{\code{$reff_ellipse}}{The observed measurement radius.}
#' And optionally, specified by the \code{dispersion_analysis} parameter, the mean and median
#' values of the line-of-sight velocity dispersion (\code{$dispersion_analysis}). The observational
#' images will also be plotted.
#'@examples
#' vsigma = find_vsigma(filename     = system.file("extdata", 'S0_vignette', package="SimSpin"),
#'                      r200         = 200,
#'                      z            = 0.1,
#'                      fov          = 15,
#'                      ap_shape     = "circular",
#'                      central_wvl  = 4800,
#'                      lsf_fwhm     = 2.65,
#'                      pixel_sscale = 0.5,
#'                      pixel_vscale = 1.04,
#'                      inc_deg      = 0,
#'                      m2l_disc     = 2,
#'                      m2l_bulge    = 1,
#'                      threshold    = 25,
#'                      measure_type = list(type = "fit", fac = 1))
#'
#'
#' vsigma = find_vsigma(filename        = system.file("extdata", 'S0_vignette', package="SimSpin"),
#'                      r200            = 200,
#'                      z               = 0.1,
#'                      fov             = 15,
#'                      ap_shape        = "circular",
#'                      central_wvl     = 4800,
#'                      lsf_fwhm        = 2.65,
#'                      pixel_sscale    = 0.5,
#'                      pixel_vscale    = 1.04,
#'                      inc_deg         = 0,
#'                      m2l_disc        = 2,
#'                      m2l_bulge       = 1,
#'                      threshold       = 25,
#'                      measure_type    = list(type = "specified",
#'                                             axis_ratio = data.frame("a"=3.5, "b"=1.7),
#'                                             fac = 1),
#'                      blur            = list("psf" = "Moffat", "fwhm" = 2))
#'


find_vsigma = function(filename, ptype = NA, r200 = 200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                       pixel_sscale, pixel_vscale, inc_deg, m2l_disc, m2l_bulge, threshold,
                       measure_type = list(type="fit", fac=1), blur){

  if (missing(blur)) {                                     # IF spatial blurring IS NOT requested

    observe_data = obs_data_prep(filename, ptype, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg, m2l_disc, m2l_bulge)
                                                           # prep simulation data in observer units
    ifu_imgs     = ifu_cube(observe_data, threshold)       # construct IFU data cube

    if (measure_type$type == "fit"){                       # fit Reff from the unblurred counts_img
      reff_ar      = find_reff(filename, ptype = NA, r200, inc_deg,
                               axis_ratio = ifu_imgs$axis_ratio,
                               angular_size = observe_data$angular_size)
                                                           # Reff from data & measured axis ratio
      vsigma       = obs_vsigma(ifu_datacube = ifu_imgs,
                                reff_axisratio = measure_type$fac * reff_ar,
                                sbinsize = observe_data$sbinsize)
                                                           # measure lambdaR within specified Reff
    }

    if (measure_type$type == "specified"){                 # fitting Reff from specified axis_ratio
      reff_ar      = find_reff(filename, ptype = NA, r200, inc_deg,
                               axis_ratio = measure_type$axis_ratio,
                               angular_size = observe_data$angular_size)
                                                           # Reff from data & supplied axis ratio
      vsigma       = obs_vsigma(ifu_datacube = ifu_imgs,
                                reff_axisratio = measure_type$fac * reff_ar,
                                sbinsize = observe_data$sbinsize)
                                                           # measure lambdaR within specified Reff
    }

    plot_ifu(vsigma, appregion = observe_data$appregion)   # plot IFU images


    output       = list("datacube" = ifu_imgs$cube, "xbin_labels" = ifu_imgs$xbin_labels,
                        "ybin_labels" = ifu_imgs$ybin_labels, "vbin_labels" = ifu_imgs$vbin_labels,
                        "axis_ratio" = reff_ar, "vsigma" = vsigma$obs_vsigma,
                        "counts_img" = vsigma$counts_img, "velocity_img" = vsigma$velocity_img,
                        "dispersion_img" = vsigma$dispersion_img,
                        "reff_ellipse" = vsigma$reff_ellipse, "appregion" = observe_data$appregion,
                        "dispersion_analysis" = vsigma$dispersion_analysis)

    return(output)

  } else {                                                 # IF spatial blurring IS requested

    observe_data = obs_data_prep(filename, ptype, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg, m2l_disc, m2l_bulge)
                                                           # prep simulation data in observer units
    ifu_imgs     = ifu_cube(observe_data, threshold)       # construct IFU data cube
    blur_imgs    = blur_cube(ifu_imgs, sbinsize = observe_data$sbinsize, psf = blur$psf,
                             fwhm = blur$fwhm, angular_size = observe_data$angular_size)
                                                           # blur IFU cube

    if (measure_type$type == "fit"){                       # fit Reff from the unblurred counts_img
      reff_ar      = find_reff(filename, ptype = NA, r200, inc_deg,
                               axis_ratio = ifu_imgs$axis_ratio,
                               angular_size = observe_data$angular_size)
                                                           # Reff from data & measured axis ratio
      vsigma       = obs_vsigma(ifu_datacube = blur_imgs,
                                reff_axisratio = measure_type$fac * reff_ar,
                                sbinsize = observe_data$sbinsize)
                                                           # measure vsigma within number of Reff
    }

    if (measure_type$type == "specified"){                 # fitting Reff from specified axis_ratio
      reff_ar      = find_reff(filename, ptype = NA, r200, inc_deg,
                               axis_ratio = measure_type$axis_ratio,
                               angular_size = observe_data$angular_size)
                                                           # Reff from data and measured axis ratio
      vsigma       = obs_vsigma(ifu_datacube = blur_imgs,
                                reff_axisratio = measure_type$fac * reff_ar,
                                sbinsize = observe_data$sbinsize)
                                                           # measure vsigma within number of Reff
    }

    plot_ifu(vsigma, appregion = observe_data$appregion)   # plot IFU images

    output = list("datacube" = blur_imgs$cube, "xbin_labels" = blur_imgs$xbin_labels,
                  "ybin_labels" = blur_imgs$ybin_labels, "vbin_labels" = blur_imgs$vbin_labels,
                  "axis_ratio" = reff_ar, "vsigma" = vsigma$obs_vsigma,
                  "counts_img" = vsigma$counts_img, "velocity_img" = vsigma$velocity_img,
                  "dispersion_img" = vsigma$dispersion_img,
                  "reff_ellipse" = vsigma$reff_ellipse, "appregion" = observe_data$appregion)

    return(output)

  }
}
