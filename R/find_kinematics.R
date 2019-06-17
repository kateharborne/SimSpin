# Kate Harborne (last edit - 16/06/2019)
#'Measuring observable galaxy kinematics.
#'
#'The purpose of this basic function is to use the \code{SimSpin} package sub-functions to measure
#' both the observable spin parameter, \eqn{\lambda_R}, and ratio V/\eqn{\sigma} of a simulated
#' galaxy model. This function will call each required sub-function (\code{\link{obs_data_prep}},
#' \code{\link{ifu_cube}}, \code{\link{blur_cube}}, \code{\link{find_reff}},
#' \code{\link{obs_lambda}}, \code{\link{obs_vsigma}} and \code{\link{plot_ifu}}) and return the
#' \eqn{\lambda_R} and the V/\eqn{\sigma} ratio within a specified measurement radius.
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
#'@param threshold The flux threshold of the observation.
#'@param measure_type A list specifying the radius within which the kinematics are measured. There
#' are three options for this:
#' \enumerate{
#' \item If \code{list("type" = "fit", "fac" = 1)}, the kinematics are measured in a specified
#'  multiple of \code{$fac *} \eqn{R_{eff}}, where \eqn{R_{eff}} has been calculated from the
#'  unblurred galaxy counts image.
#' \item If \code{list("type" = "specified", "fract" = 0.5, "axis_ratio" = data.frame("a" = 2,
#'  "b" = 1))}, the effective radius of the galaxy is measured from the unblurred counts image as
#'  in "fit", but the axis ratio of the grown ellipse is kept at the supplied axis ratio and grown
#'  until it contains some fraction (\code{$fract}) of the total particles. The kinematics are
#'  calculated within this ellipse.
#' \item Finally, if \code{list("type" = "fixed", "fac" = 1, "axis_ratio" = data.frame("a" = 2,
#' "b" = 1, "ang" = 90))}, the kinematics are measured within an ellipse described by the supplied
#' "axis_ratio" in kpc (or a multiple of that ellipse size given by "fac") at the position angle
#' "ang" -  no measurement of the effective radius is made, assuming that the supplied values are
#' determined using another package.
#' }
#'@param blur \emph{Optional} Specify if you wish to apply observational seeing effects to the
#' cube. A list of the form \code{list("psf" = "Moffat", "fwhm" = 0.5)}. \code{"psf"} specifies
#' the shape of the PSF chosen and may be either \code{"Moffat"} or \code{"Gaussian"}.
#' \code{"fwhm"} is a numeric specifying the full-width half-maximum of the PSF given in units of
#' arcseconds.
#'@param dispersion_analysis \emph{Optional} If specified as \code{TRUE}, the code will output the
#' mean and median values of the LOS velocity dispersion. Default is \code{FALSE}.
#'@param IFU_plot \emph{Optional} If specified \code{FALSE}, the function will not output the IFU flux,
#'LOS velocity and LOS velocity dispersion images. Default is \code{TRUE}, where plots are output
#'automatically.
#'@return A list containing:
#' \item{\code{$datacube}}{The 3D array corresponding to the kinematic data cube.}
#' \item{\code{$xbin_labels}}{Bin labels for the x-spatial dimension.}
#' \item{\code{$ybin_labels}}{Bin labels for the y-spatial dimension.}
#' \item{\code{$vbin_labels}}{Bin labels for the velocity dimension.}
#' \item{\code{$axis_ratio}}{The axis ratio of the observed galaxy in the form of a data frame where
#'  \code{$a} is the semi-major axis and \code{$b} is the semi-minor axis given in kpc.}
#' \item{\code{$lambda_r}}{The observed spin parameter \eqn{\lambda_R}}
#' \item{\code{$vsigma}}{The observed V/\eqn{\sigma}.}
#' \item{\code{$counts_img}}{The observed flux image.}
#' \item{\code{$velocity_img}}{The observed line-of-sight velocity image.}
#' \item{\code{$dispersion_img}}{The observed line-of-sight velocity dispersion image.}
#' \item{\code{$angular_size}}{The angular size of the galaxy in kpc/arcecond at the provided
#'  redshift.}
#' \item{\code{$sbinsize}}{The size of the spatial bins in kpc.}
#' \item{\code{$vbinsize}}{The size of the velocity bins in km/s.}
#' \item{\code{$appregion}}{The aperture region mask used to remove flux outside of the specified
#'  aperture.}
#'
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' kinematics = find_kinematics(simdata      = galaxy_data,
#'                              r200         = 200,
#'                              z            = 0.1,
#'                              fov          = 15,
#'                              ap_shape     = "circular",
#'                              central_wvl  = 4800,
#'                              lsf_fwhm     = 2.65,
#'                              pixel_sscale = 0.5,
#'                              pixel_vscale = 1.04,
#'                              inc_deg      = 0,
#'                              threshold    = 25,
#'                              measure_type = list(type = "fixed",
#'                                                  axis_ratio = data.frame("a"=3.5, "b"=1.7, "angle"=90),
#'                                                  fac = 1),
#'                              IFU_plot     = FALSE)
#'

find_kinematics=function(simdata, r200 = 200, z=0.05, fov=15, ap_shape="circular", central_wvl=4800, lsf_fwhm=2.65,
                         pixel_sscale=0.5, pixel_vscale=1.04, inc_deg=70, threshold=25,
                         measure_type = list(type="fit", fac=1), blur,
                         dispersion_analysis = FALSE, IFU_plot = TRUE){

  if (missing(blur)) {                                     # IF spatial blurring IS NOT requested

    observe_data = obs_data_prep(simdata, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg)
    # prep simulation data in observer units
    ifu_imgs     = ifu_cube(observe_data, threshold)       # construct IFU data cube

    if (measure_type$type == "fit"){                       # fit Reff from the unblurred counts_img
      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = ifu_imgs$axis_ratio,
                               angular_size = observe_data$angular_size)
      # Reff from data & measured axis ratio
      reff_ar$a_kpc = reff_ar$a_kpc * measure_type$fac
      reff_ar$b_kpc = reff_ar$b_kpc * measure_type$fac
      reff_ar$a_arcsec = reff_ar$a_arcsec * measure_type$fac
      reff_ar$b_arcsec = reff_ar$b_arcsec * measure_type$fac

      kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                  reff_axisratio = reff_ar,
                                  sbinsize = observe_data$sbinsize, dispersion_analysis)
      # measure kinematics within specified Reff
    }

    if (measure_type$type == "specified"){                 # fitting Reff from specified axis_ratio
      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = measure_type$axis_ratio,
                               angular_size = observe_data$angular_size,
                               fract = measure_type$fract)
      # Reff from data & supplied axis ratio
      kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                  reff_axisratio = reff_ar,
                                  sbinsize = observe_data$sbinsize, dispersion_analysis)
      # measure kinematics within specified Reff
    }

    if (measure_type$type == "fixed"){                     # fitting Reff from specified axis_ratio
      reff_ar      = data.frame("a_kpc"    = measure_type$axis_ratio$a,
                                "b_kpc"    = measure_type$axis_ratio$b,
                                "a_arcsec" = measure_type$axis_ratio$a / observe_data$angular_size,
                                "b_arcsec" = measure_type$axis_ratio$b / observe_data$angular_size,
                                "angle"    = measure_type$axis_ratio$ang)
      # Reff at fixed specification
      reff_ar$a_kpc = reff_ar$a_kpc * measure_type$fac
      reff_ar$b_kpc = reff_ar$b_kpc * measure_type$fac
      reff_ar$a_arcsec = reff_ar$a_arcsec * measure_type$fac
      reff_ar$b_arcsec = reff_ar$b_arcsec * measure_type$fac

      kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                  reff_axisratio = reff_ar,
                                  sbinsize = observe_data$sbinsize, dispersion_analysis)
      # measure kinematics within specified Reff

    }

    if (dispersion_analysis == TRUE) {
      output       = list("datacube"=ifu_imgs$cube, "xbin_labels"=ifu_imgs$xbin_labels,
                          "ybin_labels"=ifu_imgs$ybin_labels, "vbin_labels"=ifu_imgs$vbin_labels,
                          "axis_ratio"=reff_ar, "lambda_r"=kinematics$obs_lambdar, "vsigma" = kinematics$obs_vsigma,
                          "counts_img"=kinematics$counts_img, "velocity_img"=kinematics$velocity_img,
                          "dispersion_img"=kinematics$dispersion_img,
                          "angular_size"=observe_data$angular_size,
                          "sbinsize"=observe_data$sbinsize,
                          "vbinsize"=observe_data$vbinsize,
                          "d_L"=observe_data$d_L,
                          "appregion"=observe_data$appregion,
                          "dispersion_analysis"=kinematics$dispersion_analysis)
    } else {
      output       = list("datacube"=ifu_imgs$cube, "xbin_labels"=ifu_imgs$xbin_labels,
                          "ybin_labels"=ifu_imgs$ybin_labels, "vbin_labels"=ifu_imgs$vbin_labels,
                          "axis_ratio"=reff_ar, "lambda_r"=kinematics$obs_lambdar, "vsigma" = kinematics$obs_vsigma,
                          "counts_img"=kinematics$counts_img, "velocity_img"=kinematics$velocity_img,
                          "dispersion_img"=kinematics$dispersion_img,
                          "angular_size"=observe_data$angular_size,
                          "sbinsize"=observe_data$sbinsize,
                          "vbinsize"=observe_data$vbinsize,
                          "d_L"=observe_data$d_L,
                          "appregion"=observe_data$appregion)
    }

    if (IFU_plot == TRUE){
      plot_ifu(output)   # plot IFU images
    }

    return(output)

  } else {                                                 # IF spatial blurring IS requested

    observe_data = obs_data_prep(simdata, r200, z, fov, ap_shape, central_wvl, lsf_fwhm,
                                 pixel_sscale, pixel_vscale, inc_deg)
    # prep simulation data in observer units
    ifu_imgs     = ifu_cube(observe_data, threshold)       # construct IFU data cube
    blur_imgs    = blur_cube(ifu_imgs, sbinsize = observe_data$sbinsize, psf = blur$psf,
                             fwhm = blur$fwhm, angular_size = observe_data$angular_size)
    # blur IFU cube

    if (measure_type$type == "fit"){                       # fit Reff from the unblurred counts_img
      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = ifu_imgs$axis_ratio,
                               angular_size = observe_data$angular_size)
      # Reff from data & measured axis ratio
      reff_ar$a_kpc = reff_ar$a_kpc * measure_type$fac
      reff_ar$b_kpc = reff_ar$b_kpc * measure_type$fac
      reff_ar$a_arcsec = reff_ar$a_arcsec * measure_type$fac
      reff_ar$b_arcsec = reff_ar$b_arcsec * measure_type$fac

      kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                  reff_axisratio = reff_ar,
                                  sbinsize = observe_data$sbinsize, dispersion_analysis)
      # measure kinematics within specified Reff

    }

    if (measure_type$type == "specified"){                 # fitting Reff from specified axis_ratio
      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = measure_type$axis_ratio,
                               angular_size = observe_data$angular_size,
                               fract = measure_type$fract)
      # Reff from data and measured axis ratio

      kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                  reff_axisratio = reff_ar,
                                  sbinsize = observe_data$sbinsize, dispersion_analysis)
      # measure kinematics within specified Reff
    }

    if (measure_type$type == "fixed"){                 # fitting Reff from specified axis_ratio
      reff_ar      = data.frame("a_kpc"    = measure_type$axis_ratio$a,
                                "b_kpc"    = measure_type$axis_ratio$b,
                                "a_arcsec" = measure_type$axis_ratio$a / observe_data$angular_size,
                                "b_arcsec" = measure_type$axis_ratio$b / observe_data$angular_size,
                                "angle"    = measure_type$axis_ratio$ang)
      # Reff at fixed specification
      reff_ar$a_kpc = reff_ar$a_kpc * measure_type$fac
      reff_ar$b_kpc = reff_ar$b_kpc * measure_type$fac
      reff_ar$a_arcsec = reff_ar$a_arcsec * measure_type$fac
      reff_ar$b_arcsec = reff_ar$b_arcsec * measure_type$fac

      kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                  reff_axisratio = reff_ar,
                                  sbinsize = observe_data$sbinsize, dispersion_analysis)
      # measure kinematics within specified Reff
    }

    if (dispersion_analysis == TRUE) {
      output = list("datacube" = blur_imgs$cube, "xbin_labels" = blur_imgs$xbin_labels,
                    "ybin_labels" = blur_imgs$ybin_labels, "vbin_labels" = blur_imgs$vbin_labels,
                    "axis_ratio" = reff_ar, "lambda_r" = kinematics$obs_lambdar, "vsigma" = kinematics$obs_vsigma,
                    "counts_img" = kinematics$counts_img, "velocity_img" = kinematics$velocity_img,
                    "dispersion_img" = kinematics$dispersion_img,
                    "angular_size" = observe_data$angular_size,
                    "sbinsize"= observe_data$sbinsize,
                    "vbinsize" = observe_data$vbinsize,
                    "d_L"=observe_data$d_L,
                    "appregion" = observe_data$appregion,
                    "dispersion_analysis" = kinematics$dispersion_analysis)
    } else {
      output = list("datacube" = blur_imgs$cube, "xbin_labels" = blur_imgs$xbin_labels,
                    "ybin_labels" = blur_imgs$ybin_labels, "vbin_labels" = blur_imgs$vbin_labels,
                    "axis_ratio" = reff_ar, "lambda_r" = kinematics$obs_lambdar, "vsigma" = kinematics$obs_vsigma,
                    "counts_img" = kinematics$counts_img, "velocity_img" = kinematics$velocity_img,
                    "dispersion_img" = kinematics$dispersion_img,
                    "angular_size" = observe_data$angular_size,
                    "sbinsize"= observe_data$sbinsize,
                    "vbinsize" = observe_data$vbinsize,
                    "d_L"=observe_data$d_L,
                    "appregion" = observe_data$appregion)
    }

    if (IFU_plot == TRUE){
      plot_ifu(output)   # plot IFU images
    }

    return(output)

  }

}
