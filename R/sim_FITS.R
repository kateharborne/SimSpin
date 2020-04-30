# Kate Harborne (last edit - 03/12/2018)
#'Generate a FITS file for the data cube produced in SimSpin
#'
#'The purpose of this function is to write out the data cube produced by SimSpin into the common
#' astronomy FITS data format.
#'
#'@param out_data The list output from \code{\link{build_datacube}} or
#' \code{\link{find_kinematics}}.
#'@param pixel_sscale The corresponding spatial pixel scale associated with a given telescope
#' output in arcseconds.
#'@param pixel_vscale The corresponding velocity pixel scale associated with a given telescope
#' filter output in angstroms.
#'@param z The redshift projected distance at which
#'\code{\link{find_kinematics}} was run.
#'@param psf_fwhm The FWHM of the PSF used for spatial blurring used when
#' \code{\link{find_kinematics}} was run.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param r50 The half mass radius specified by the simulation, kpc.
#'@param Hdisk The scale length of the disk component in the simulation, kpc.
#'@param Ahalo The scale height of the dark matter halo component in the simulation, kpc.
#'@param Abulge The scale height of the bulge component in the simulation, kpc.
#'@param out_file A string describing the path and file name of the FITS file to be written.
#'@param obs_name A string that describes the name of the observation. Default is "SimSpin
#'datacube".
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
#'                                           axis_ratio = data.frame("a"=3.5, "b"=1.7, "angle"=90),
#'                                           fac = 1),
#'                       IFU_plot     = FALSE,
#'                       multi_thread = FALSE)
#' sim_FITS(out_data = kin, pixel_sscale = 0.5, pixel_vscale = 1.04, z=0.1, r200=200, r50=10,
#' Hdisk=5, Ahalo=34.5, Abulge=3.45, out_file = "simdata_example.fits")
#' unlink("simdata_example.fits")


sim_FITS = function(out_data, z, pixel_sscale, pixel_vscale, psf_fwhm=0, r200=200, r50=10, Hdisk=5.64,
                    Ahalo=34.5, Abulge=3.45, out_file, obs_name="SimSpin datacube"){

  crpix_sim = c(dim(out_data$datacube)[1]/2, dim(out_data$datacube)[2]/2, dim(out_data$datacube)[3]/2)
  cdelt_sim = c(pixel_sscale, pixel_sscale, pixel_vscale)
  crval_sim = c(round(out_data$xbin_labels[length(out_data$xbin_labels)/2] + cdelt_sim[1]/2),
                round(out_data$xbin_labels[length(out_data$xbin_labels)/2] + cdelt_sim[2]/2),
                out_data$vbin_labels[crpix_sim[3]] + cdelt_sim[3]/2)
  len_sim   = c(dim(out_data$datacube)[1], dim(out_data$datacube)[2], dim(out_data$datacube)[3])
  ctype_sim = c("X-Spaxel Size", "Y-Spaxel Size", "Voxel Size")
  cunit_sim = c("arcsec", "arcsec", "angstrom")
  head_sim = FITSio::newKwv("KEYWORD", "VALUE", "NOTE")
  head_sim = FITSio::addKwv("REDSHIFT", z, "redshift, z", header=head_sim)
  head_sim = FITSio::addKwv("SPIXSIZE", out_data$sbinsize, "spatial size, kpc/pixel", header=head_sim)
  head_sim = FITSio::addKwv("VPIXSIZE", out_data$vbinsize, "velcoity size, kms-1/pixel", header=head_sim)
  head_sim = FITSio::addKwv("PSFFWHM", psf_fwhm, "FWHM, arcsec", header=head_sim)
  head_sim = FITSio::addKwv("SIMR200", r200, "virial radius from sim, kpc", header=head_sim)
  head_sim = FITSio::addKwv("SIMHDISK", Hdisk, "disk scale length from sim, kpc", header=head_sim)
  head_sim = FITSio::addKwv("SIMAHALO", Ahalo, "halo scale height from sim, kpc", header=head_sim)
  head_sim = FITSio::addKwv("SIMABULG", Abulge, "bulge scale height from sim, kpc", header=head_sim)
  head_sim = FITSio::addKwv("HALFMASS", r50, "half mass radius of sim, kpc", header=head_sim)

  FITSio::writeFITSim(out_data$datacube, file = out_file, c1 = obs_name,
                      crpixn = crpix_sim,
                      crvaln = crval_sim,
                      cdeltn = cdelt_sim,
                      ctypen = ctype_sim,
                      cunitn = cunit_sim,
                      header = head_sim)


  return(sprintf("New FITS file written."))

}
