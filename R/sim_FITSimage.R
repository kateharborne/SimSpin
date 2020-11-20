# Kate Harborne (last edit - 03/12/2018)
#'Generate a FITS file for the data cube produced in SimSpin
#'
#'The purpose of this function is to write out the data cube produced by SimSpin into the common
#' astronomy FITS data format.
#'
#'@param out_image The counts/velocity/dispersion image output from \code{\link{find_kinematics}}
#'  that you wish to save.
#'@param out_data The list output from \code{\link{build_datacube}}, or
#'\code{\link{find_kinematics}}.
#'@param pixel_sscale The corresponding spatial pixel scale associated with a given telescope
#' output in arcseconds.
#'@param z The redshift projected distance at which
#' \code{\link{find_kinematics}} was run.
#'@param psf_fwhm The FWHM of the PSF used for spatial blurring used when
#' \code{\link{find_kinematics}} was run.
#'@param out_file A string describing the path and file name of the FITS file to be written.
#'@param obs_name A string that describes the name of the observation. Default is "SimSpin
#' datacube".
#'@param addSky A boolean to specify whether to add sky noise to the output images. Default is
#' FALSE. If TRUE, further parameters including \code{mag_threshold} and \code{mag_zero} described
#' below.
#'@param threshold The magnitude limit of the observation.
#'@param mag_zero The magnitude zero point with regards to the mangitude system being used (e.g.
#' AB or Vega).
#'
#'@return Outputs a standard format FITS file.
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' kin = find_kinematics(simdata      = galaxy_data,
#'                       r200         = 200,
#'                       z            = 0.1,
#'                       fov          = 15,
#'                       ap_shape     = "circular",
#'                       central_wvl  = 4800,
#'                       lsf_fwhm     = 2.65,
#'                       pixel_sscale = 0.5,
#'                       pixel_vscale = 1.04,
#'                       inc_deg      = 0,
#'                       threshold    = 25,
#'                       measure_type = list(type = "fixed",
#'                                           axis_ratio = data.frame("a"=3.5,
#'                                                                   "b"=1.7,
#'                                                                   "angle"=90),
#'                                           fac = 1),
#'                       IFU_plot     = FALSE,
#'                       multi_thread = FALSE)
#' sim_FITSimage(out_image = kin$counts_img, out_data = kin, z = 0.1, pixel_sscale = 0.5,
#' out_file = "simdata_example.fits")
#' unlink("simdata_example.fits")


sim_FITSimage = function(out_image, out_data, z, pixel_sscale, psf_fwhm=0,
                         out_file, obs_name="SimSpin datacube",
                         addSky = FALSE, threshold=25, mag_zero=8.9){

  if (addSky){
    skyRMS = ProFound::profoundSB2Flux(threshold, mag_zero, pixel_sscale)
    out_image = out_image+rnorm(dim(out_image)[1]^2, sd=skyRMS)
  }

  crpix_sim = c(dim(out_data$cube)[1]/2, dim(out_data$cube)[2]/2)
  cdelt_sim = c(pixel_sscale, pixel_sscale)
  crval_sim = c(round(out_data$xbin_labels[length(out_data$xbin_labels)/2] + cdelt_sim[1]/2),
                round(out_data$xbin_labels[length(out_data$xbin_labels)/2] + cdelt_sim[2]/2))
  len_sim   = c(dim(out_data$cube)[1], dim(out_data$cube)[2])
  ctype_sim = c("X-Spaxel Size", "Y-Spaxel Size")
  cunit_sim = c("arcsec", "arcsec")
  head_sim = FITSio::newKwv("KEYWORD", "VALUE", "NOTE")
  head_sim = FITSio::addKwv("REDSHIFT", z, "redshift, z", header=head_sim)
  head_sim = FITSio::addKwv("SPIXSIZE", out_data$sbinsize, "spatial size, kpc/pixel", header=head_sim)
  head_sim = FITSio::addKwv("PSFFWHM", psf_fwhm, "FWHM, arcsec", header=head_sim)

  FITSio::writeFITSim(out_image, file = out_file, c1 = obs_name,
                      crpixn = crpix_sim,
                      crvaln = crval_sim,
                      cdeltn = cdelt_sim,
                      ctypen = ctype_sim,
                      cunitn = cunit_sim,
                      header = head_sim)


  return(sprintf("New FITS file written."))

}
