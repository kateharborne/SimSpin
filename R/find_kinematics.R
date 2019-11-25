# Kate Harborne (last edit - 15/11/2019)
#'Measuring observable galaxy kinematics.
#'
#'The purpose of this basic function is to use the \code{SimSpin} package sub-functions to measure
#' both the observable spin parameter, \eqn{\lambda_R}, and ratio V/\eqn{\sigma} of a simulated
#' galaxy model. This function will call each required sub-function (\code{\link{obs_data_prep}},
#' \code{\link{flux_grid}}, \code{\link{ifu_cube}}, \code{\link{blur_cube}},
#' \code{\link{find_reff}}, \code{\link{kin_calc}}, and \code{\link{plot_ifu}}) and return the
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
#'@param filter If Age/Metallicity is supplied, the filter within which the SED is generated.
#'Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
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
#'@param radius_type The method of computing radii - "Circular" i.e. \eqn{r^{2} = x^{2} + y{2}} or
#'"Elliptical" where r is the semi-major axis of the ellipse having an axis ratio \eqn{b/a} on
#'which the pixel lies, i.e.
#'\eqn{r^{2} = \frac{x^{2} (1 - \epsilon)^{2} + y^{2}}{(1 - \epsilon)^2}}. Default is "Both" such
#'that both \eqn{\lambda_R} values are returned.
#'@param blur \emph{Optional} Specify if you wish to apply observational seeing effects to the
#' cube. A list of the form \code{list("psf" = "Moffat", "fwhm" = 0.5)}. \code{"psf"} specifies
#' the shape of the PSF chosen and may be either \code{"Moffat"} or \code{"Gaussian"}.
#' \code{"fwhm"} is a numeric specifying the full-width half-maximum of the PSF given in units of
#' arcseconds.
#'@param addSky A boolean to specify whether to add sky noise to the output images. Default is
#' FALSE. If TRUE, further parameters including \code{mag_threshold} and \code{mag_zero} described
#' below.
#'@param mag_zero The magnitude zero point with regards to the mangitude system being used (e.g.
#' AB or Vega).
#'@param IFU_plot \emph{Optional} If specified \code{FALSE}, the function will not output the IFU flux,
#'LOS velocity and LOS velocity dispersion images. Default is \code{TRUE}, where plots are output
#'automatically.
#'@return A list containing:
#' \item{\code{$datacube}}{The 3D array corresponding to the kinematic data cube.}
#' \item{\code{$axis_ratio}}{The axis ratio of the observed galaxy in the form of a data frame where
#'  \code{$a} is the semi-major axis and \code{$b} is the semi-minor axis given in kpc.}
#' \item{\code{$obs_lambdar}}{The observed spin parameter \eqn{\lambda_R} measured with circular
#' radii. \emph{(When \code{radius_type = "Both"} or \code{"Circular"}.)}}
#' \item{\code{$obs_elambdar}}{The observed spin parameter \eqn{\lambda_R} measured with elliptical
#' radii. \emph{(When \code{radius_type = "Both"} or \code{"Elliptical"}.)}}
#' \item{\code{$vsigma}}{The observed V/\eqn{\sigma}.}
#' \item{\code{$flux_img}}{The observed flux image.}
#' \item{\code{$velocity_img}}{The observed line-of-sight velocity image.}
#' \item{\code{$dispersion_img}}{The observed line-of-sight velocity dispersion image.}
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
#'                                                  axis_ratio = data.frame("a"=3.5,
#'                                                                          "b"=1.7,
#'                                                                          "angle"=90),
#'                                                  fac = 1),
#'                              IFU_plot     = FALSE)
#'

find_kinematics=function(simdata, r200 = 200, z=0.05, fov=15, ap_shape="circular", central_wvl=4800, lsf_fwhm=2.65,
                         pixel_sscale=0.5, pixel_vscale=1.04, inc_deg=70, threshold=25, filter="g",
                         measure_type = list(type="fit", fac=1), blur,
                         radius_type = "Both", addSky = FALSE, mag_zero = 8.9, IFU_plot = TRUE){
                         pixel_sscale=0.5, pixel_vscale=1.04, inc_deg=70, align=TRUE, filter="g", threshold=25,
                         measure_type = list(type="fit", fac=1), blur, radius_type = "Both", IFU_plot=FALSE){

  if (missing(blur)) {                                     # IF spatial blurring IS NOT requested

    observe_data = obs_data_prep(simdata = simdata, r200 = r200, z = z, fov = fov, ap_shape = ap_shape,
                                 central_wvl = central_wvl, lsf_fwhm = lsf_fwhm, pixel_sscale = pixel_sscale,
                                 pixel_vscale = pixel_vscale, inc_deg = inc_deg, align = TRUE) # prep simulation data in observer units
    fluxes = flux_grid(obs_data = observe_data, filter = filter)
    ifu_imgs = ifu_cube(obs_data = observe_data, flux_data = fluxes) # construct IFU data cube
    images = obs_imgs(obs_data = observe_data, ifu_datacube = ifu_imgs, threshold = threshold)

    if (measure_type$type == "fit"){                       # fit Reff from the unblurred counts_img
      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = images$axis_ratio)
      # Reff from data & measured axis ratio
      reff_ar$a = reff_ar$a * measure_type$fac
      reff_ar$b = reff_ar$b * measure_type$fac

      if (addSky){
        kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type,
                                    addSky = TRUE, mag_zero = mag_zero, threshold = threshold,
                                    pixel_sscale = pixel_sscale)
      } else {
        kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type)
        # measure kinematics within specified Reff
      }
      kinematics = kin_calc(obs_data = observe_data,
                            obs_images = images,axis_ratio = reff_ar,
                            radius_type = radius_type)
      # measure kinematics within specified Reff
    }

    if (measure_type$type == "specified"){                 # fitting Reff from specified axis_ratio
      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = measure_type$axis_ratio,
                               fract = measure_type$fract)
      # Reff from data & supplied axis ratio
      if (addSky){
        kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type,
                                    addSky = TRUE, mag_zero = mag_zero, threshold = threshold,
                                    pixel_sscale = pixel_sscale)

      } else {
        kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type)
        # measure kinematics within specified Reff
      }

      kinematics = kin_calc(obs_data = observe_data,
                            obs_images = images, axis_ratio = reff_ar,
                            radius_type = radius_type)
      # measure kinematics within specified Reff
    }

    if (measure_type$type == "fixed"){                     # fitting Reff from specified axis_ratio
      reff_ar      = data.frame("a"    = measure_type$axis_ratio$a,
                                "b"    = measure_type$axis_ratio$b,
                                "ang"  = measure_type$axis_ratio$ang)
      # Reff at fixed specification

      if (addSky){
        kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type,
                                    addSky = TRUE, mag_zero = mag_zero, threshold = threshold,
                                    pixel_sscale = pixel_sscale)

      } else {
        kinematics = obs_kinematics(ifu_datacube = ifu_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type)
        # measure kinematics within specified Reff
      }

    }

    if (radius_type == "Both" | radius_type == "both") {
      output       = list("datacube"=ifu_imgs$cube, "xbin_labels"=ifu_imgs$xbin_labels,
                          "ybin_labels"=ifu_imgs$ybin_labels, "vbin_labels"=ifu_imgs$vbin_labels,
                          "axis_ratio"=reff_ar, "lambda_r"=kinematics$obs_lambdar,
                          "elambda_r"=kinematics$obs_elambdar,  "vsigma" = kinematics$obs_vsigma,
                          "counts_img"=kinematics$counts_img, "velocity_img"=kinematics$velocity_img,
                          "dispersion_img"=kinematics$dispersion_img,
                          "angular_size"=observe_data$angular_size,
                          "sbinsize"=observe_data$sbinsize,
                          "vbinsize"=observe_data$vbinsize,
                          "d_L"=observe_data$d_L,
                          "appregion"=observe_data$appregion)
    } else if (radius_type == "Circular" | radius_type == "circular") {
      output       = list("datacube"=ifu_imgs$cube, "xbin_labels"=ifu_imgs$xbin_labels,
                          "ybin_labels"=ifu_imgs$ybin_labels, "vbin_labels"=ifu_imgs$vbin_labels,
                          "axis_ratio"=reff_ar, "lambda_r"=kinematics$obs_lambdar,
                          "vsigma" = kinematics$obs_vsigma,
                          "counts_img"=kinematics$counts_img, "velocity_img"=kinematics$velocity_img,
                          "dispersion_img"=kinematics$dispersion_img,
                          "angular_size"=observe_data$angular_size,
                          "sbinsize"=observe_data$sbinsize,
                          "vbinsize"=observe_data$vbinsize,
                          "d_L"=observe_data$d_L,
                          "appregion"=observe_data$appregion)
    } else if (radius_type == "Elliptical" | radius_type == "elliptical") {
      output       = list("datacube"=ifu_imgs$cube, "xbin_labels"=ifu_imgs$xbin_labels,
                          "ybin_labels"=ifu_imgs$ybin_labels, "vbin_labels"=ifu_imgs$vbin_labels,
                          "axis_ratio"=reff_ar,
                          "elambda_r"=kinematics$obs_elambdar,  "vsigma" = kinematics$obs_vsigma,
                          "counts_img"=kinematics$counts_img, "velocity_img"=kinematics$velocity_img,
                          "dispersion_img"=kinematics$dispersion_img,
                          "angular_size"=observe_data$angular_size,
                          "sbinsize"=observe_data$sbinsize,
                          "vbinsize"=observe_data$vbinsize,
                          "d_L"=observe_data$d_L,
                          "appregion"=observe_data$appregion)
    }
      reff_ar$a = reff_ar$a * measure_type$fac
      reff_ar$b = reff_ar$b * measure_type$fac

      kinematics = kin_calc(obs_data = observe_data,
                            obs_images = images, axis_ratio = reff_ar,
                            radius_type = radius_type)
      # measure kinematics within specified Reff

    }

      output       = list("datacube"=ifu_imgs$cube,
                          "axis_ratio"=reff_ar,
                          "obs_lambdar"=kinematics$obs_lambdar,
                          "obs_elambdar"=kinematics$obs_elambdar,
                          "obs_vsigma" = kinematics$obs_vsigma,
                          "flux_img"=images$flux_img,
                          "velocity_img"=images$velocity_img,
                          "dispersion_img"=images$dispersion_img,
                          "obs_data" = observe_data,
                          "obs_images" = images)


    if (IFU_plot != FALSE){
      plot_ifu(obs_data = observe_data, obs_images = images, reff=TRUE, axis_ratio=reff_ar, which_plots = IFU_plot)
      # plot IFU images
    }

    return(output)

  } else {                                                 # IF spatial blurring IS requested

    observe_data = obs_data_prep(simdata = simdata, r200 = r200, z = z, fov = fov, ap_shape = ap_shape,
                                 central_wvl = central_wvl, lsf_fwhm = lsf_fwhm, pixel_sscale = pixel_sscale,
                                 pixel_vscale = pixel_vscale, inc_deg = inc_deg, align = TRUE) # prep simulation data in observer units
    fluxes = flux_grid(obs_data = observe_data, filter = filter)
    ifu_imgs = ifu_cube(obs_data = observe_data, flux_data = fluxes) # construct IFU data cube
    blur_imgs = blur_cube(obs_data = observe_data, ifu_datacube = ifu_imgs, psf = blur$psf,
                             fwhm = blur$fwhm) # blur IFU cube
    blur_images = obs_imgs(obs_data = observe_data, ifu_datacube = blur_imgs, threshold = threshold)


    if (measure_type$type == "fit"){                       # fit Reff from the unblurred counts_img
      images = obs_imgs(obs_data = observe_data, ifu_datacube = ifu_imgs, threshold = threshold)
      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = images$axis_ratio)
      # Reff from data & measured axis ratio
      reff_ar$a = reff_ar$a * measure_type$fac
      reff_ar$b = reff_ar$b * measure_type$fac

        kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type,
                                    addSky = TRUE, mag_zero = mag_zero, threshold = threshold,
                                    pixel_sscale = pixel_sscale)
      } else {
        kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type)
        # measure kinematics within specified Reff
      }
      kinematics = kin_calc(obs_data = observe_data,
                            obs_images = blur_images, axis_ratio = reff_ar,
                            radius_type = radius_type)
      # measure kinematics within specified Reff

    }

      reff_ar      = find_reff(simdata, r200, inc_deg,
                               axis_ratio = measure_type$axis_ratio,
                               fract = measure_type$fract)
      # Reff from data and measured axis ratio

      if (addSky){
        kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type,
                                    addSky = TRUE, mag_zero = mag_zero, threshold = threshold,
                                    pixel_sscale = pixel_sscale)

      } else {
        kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type)
        # measure kinematics within specified Reff
      }
      # Reff from data & supplied axis ratio
      kinematics = kin_calc(obs_data = observe_data,
                            obs_images = blur_images, axis_ratio = reff_ar,
                            radius_type = radius_type)
      # measure kinematics within specified Reff
    }

    if (measure_type$type == "fixed"){                 # fitting Reff from specified axis_ratio
      reff_ar      = data.frame("a"    = measure_type$axis_ratio$a,
                                "b"    = measure_type$axis_ratio$b,
                                "ang"  = measure_type$axis_ratio$ang)
      # Reff at fixed specification

      if (addSky){
        kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type,
                                    addSky = TRUE, mag_zero = mag_zero, threshold = threshold,
                                    pixel_sscale = pixel_sscale)

      } else {
        kinematics = obs_kinematics(ifu_datacube = blur_imgs,
                                    reff_axisratio = reff_ar,
                                    sbinsize = observe_data$sbinsize, radius_type=radius_type)
        # measure kinematics within specified Reff
      }
    }

    if (radius_type == "Both" | radius_type == "both") {
      output = list("datacube" = blur_imgs$cube, "xbin_labels" = blur_imgs$xbin_labels,
                    "ybin_labels" = blur_imgs$ybin_labels, "vbin_labels" = blur_imgs$vbin_labels,
                    "axis_ratio" = reff_ar, "lambda_r" = kinematics$obs_lambdar,
                    "elambda_r" = kinematics$obs_elambdar, "vsigma" = kinematics$obs_vsigma,
                    "counts_img" = kinematics$counts_img, "velocity_img" = kinematics$velocity_img,
                    "dispersion_img" = kinematics$dispersion_img,
                    "angular_size" = observe_data$angular_size,
                    "sbinsize"= observe_data$sbinsize,
                    "vbinsize" = observe_data$vbinsize,
                    "d_L"=observe_data$d_L,
                    "appregion" = observe_data$appregion)
    } else if (radius_type == "Circular" | radius_type == "circular") {
      output = list("datacube" = blur_imgs$cube, "xbin_labels" = blur_imgs$xbin_labels,
                    "ybin_labels" = blur_imgs$ybin_labels, "vbin_labels" = blur_imgs$vbin_labels,
                    "axis_ratio" = reff_ar, "lambda_r" = kinematics$obs_lambdar,
                    "vsigma" = kinematics$obs_vsigma,
                    "counts_img" = kinematics$counts_img, "velocity_img" = kinematics$velocity_img,
                    "dispersion_img" = kinematics$dispersion_img,
                    "angular_size" = observe_data$angular_size,
                    "sbinsize"= observe_data$sbinsize,
                    "vbinsize" = observe_data$vbinsize,
                    "d_L"=observe_data$d_L,
                    "appregion" = observe_data$appregion)
     } else if (radius_type == "Elliptical" | radius_type == "elliptical") {
       output = list("datacube" = blur_imgs$cube, "xbin_labels" = blur_imgs$xbin_labels,
                     "ybin_labels" = blur_imgs$ybin_labels, "vbin_labels" = blur_imgs$vbin_labels,
                     "axis_ratio" = reff_ar,
                     "elambda_r" = kinematics$obs_elambdar, "vsigma" = kinematics$obs_vsigma,
                     "counts_img" = kinematics$counts_img, "velocity_img" = kinematics$velocity_img,
                     "dispersion_img" = kinematics$dispersion_img,
                     "angular_size" = observe_data$angular_size,
                     "sbinsize"= observe_data$sbinsize,
                     "vbinsize" = observe_data$vbinsize,
                     "d_L"=observe_data$d_L,
                     "appregion" = observe_data$appregion)
      }
      reff_ar$a = reff_ar$a * measure_type$fac
      reff_ar$b = reff_ar$b * measure_type$fac

      kinematics = kin_calc(obs_data = observe_data,
                            obs_images = blur_images, axis_ratio = reff_ar,
                            radius_type = radius_type)
    }

    output       = list("datacube"=ifu_imgs$cube,
                        "axis_ratio"=reff_ar,
                        "obs_lambdar"=kinematics$obs_lambdar,
                        "obs_elambdar"=kinematics$obs_elambdar,
                        "obs_vsigma" = kinematics$obs_vsigma,
                        "flux_img"=blur_images$flux_img,
                        "velocity_img"=blur_images$velocity_img,
                        "dispersion_img"=blur_images$dispersion_img,
                        "obs_data" = observe_data,
                        "obs_images" = blur_images)

    if (IFU_plot != FALSE){
      plot_ifu(obs_data = observe_data, obs_images = images, reff=TRUE, axis_ratio=reff_ar, which_plots = IFU_plot)
    }

    return(output)

  }

}
