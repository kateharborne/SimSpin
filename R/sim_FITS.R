# Kate Harborne (last edit - 03/12/2018)
#'Generate a FITS file for the data cube produced in SimSpin
#'
#'The purpose of this function is to write out the data cube produced by SimSpin into the common
#' astronomy FITS data format.
#'
#'@param out_data The list output from \code{\link{build_datacube}}, \code{\link{find_lambda}} or
#' \code{\link{find_vsigma}}.
#'@param out_file A string describing the path and file name of the FITS file to be written.
#'@param obs_name A string that describes the name of the observation. Default is "SimSpin datacube".
#'
#'@return Outputs a standard format FITS file.
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' lambdar = find_lambda(simdata      = galaxy_data,
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
#'                                           axis_ratio = data.frame("a"=3.5, "b"=1.7, "angle"=90),
#'                                           fac = 1),
#'                       IFU_plot     = FALSE)
#' sim_FITS(out_data = lambdar, out_file = "simdata_example.fits")
#' unlink("simdata_example.fits")


sim_FITS = function(out_data, out_file, obs_name="SimSpin datacube"){

  crpix_sim = c(length(out_data$xbin_labels)/2,length(out_data$ybin_labels)/2,length(out_data$vbin_labels)/2)
  cdelt_sim = c(out_data$sbinsize, out_data$sbinsize, out_data$vbinsize)
  crval_sim = c(round(out_data$xbin_labels[crpix_sim[1]] + cdelt_sim[1]/2),
                round(out_data$xbin_labels[crpix_sim[2]] + cdelt_sim[2]/2),
                out_data$vbin_labels[crpix_sim[3]] + cdelt_sim[3]/2)
  len_sim   = c(length(out_data$xbin_labels), length(out_data$ybin_labels), length(out_data$vbin_labels))
  ctype_sim = c("X-Spaxel Size", "Y-Spaxel Size", "Voxel Size")
  cunit_sim = c("kpc", "kpc", "km s-1")

  FITSio::writeFITSim(out_data$datacube, file = out_file, c1 = obs_name,
                      crpixn = crpix_sim,
                      crvaln = crval_sim,
                      cdeltn = cdelt_sim,
                      ctypen = ctype_sim,
                      cunitn = cunit_sim)

  return(sprintf("New FITS file written."))

}
