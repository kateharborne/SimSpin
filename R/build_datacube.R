# Author: Kate Harborne
# Date: 26/10/2020
# Title: build_datacube - a function for generating a data cube from the observation
#
#'A function for making a mock spectral cube
#'
#'The purpose of this function is to generate a mock spectral IFU data cube from
#' a SimSpin file, as generated using the \code{make_simspin_file()} function.
#'
#'@param simspin_file The path to the location of the SimSpin .Rdata file OR
#' output list from \code{make_simspin_file()}.
#'@param telescope An object of the telescope class describing the
#' specifications of the observing telescope (i.e. field of view, spatial
#' resolution, wavelength resolution, etc.).
#'@param observing_strategy An object of the observing_strategy class that
#' describes the properties of the observed simulation (i.e. redshift,
#' inclination, seeing conditions).
#'@param verbose Default is \code{FALSE}. If you would like the code to give
#' updates about its progress, change this parameter to \code{TRUE}.
#'@param write_fits Default is \code{FALSE}. If you would like the code to
#' output a FITS file, change this parameter to \code{TRUE}.
#'@param output_location Optional parameter that describes the path to the FITS
#' file output if \code{write_fits = TRUE}. If \code{write_fits = TRUE} and no
#' \code{output_location} is specified, the FITS file will be written to the
#' same directory as the input \code{simspin_file}.
#'@return Returns an .fits file that contains a the generated spectral cube and
#' relevant header describing the mock observation.
#'@examples
#'ss_gadget = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata",
#'                         package = "SimSpin")
#'cube = build_datacube(simspin_file = ss_gadget,
#'                      telescope = telescope(type="SAMI"),
#'                      observing_strategy = observing_strategy())
#'

build_datacube = function(simspin_file, telescope, observing_strategy,
                          verbose = F, write_fits = F, output_location){

  if (verbose){cat("Computing observation parameters... \n")}
  observation = observation(telescope = telescope, observing_strategy = observing_strategy)

  if (typeof(simspin_file) == "character"){ # if provided with path to file
    simspin_data = readRDS(simspin_file)
  }
  if (typeof(simspin_file) == "list"){ # if provided with output list
    simspin_data = simspin_file
  }

  twisted_data = twist_galaxy(simspin_data$star_part, twist_rad = observation$twist_rad)

  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad) # projecting the galaxy to given inclination

  if (verbose){cat("Assigning particles to spaxels... \n")}
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z_obs, breaks=observation$sbin_seq, labels=F)) - (observation$sbin) # assigning particles to positions in cube

  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$pixel_region[!is.na(observation$pixel_region)],] # trimming particles that lie outside the aperture of the telescope

  original_wave  = simspin_data$wave # read original wavelengths
  wavelength = (observation$z * original_wave) + original_wave # and then applying a shift due to redshift, z
  lsf_fwhm   = (observation$z * observation$lsf_fwhm) + observation$lsf_fwhm # adjusting the LSF for the resolution at z

  spec_res_sigma_sq = lsf_fwhm^2 - (min(diff(wavelength)))^2
  if (spec_res_sigma_sq < 0){
    warning(cat("WARNING! - Wavelength resolution of provided spectra is lower than the requested telescope resolution.\n"))
    cat("LSF = ", observation$lsf_fwhm,  " A < wavelength resolution ", min(diff(wavelength)), " A. \n")
    cat("No LSF will be applied in this case.\n")
    LSF_conv = FALSE
  } else {
    LSF_conv = TRUE
    lsf_sigma = (sqrt(spec_res_sigma_sq) / (2 * sqrt(2*log(2))))
  }

  spectra = matrix(data = NA, ncol = observation$wave_bin, nrow = observation$sbin^2)
  vel_los = array(data = NA, dim = observation$sbin^2)
  dis_los = array(data = NA, dim = observation$sbin^2)

  if (verbose){cat("Generating spectra per spaxel... \n")}
  for (i in sort(unique(galaxy_data$pixel_pos))){ # computing the spectra at each occupied spatial pixel position
     particle_IDs = which(galaxy_data$pixel_pos == i)
     galaxy_sample = galaxy_data[particle_IDs,]
     intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]), nrow = length(particle_IDs), byrow = T) * (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra
     velocity_los = galaxy_sample$vy_obs # the LOS velocities of each particle
     vel_los[i] = mean(velocity_los)
     dis_los[i] = sd(velocity_los)
     wave = matrix(data = rep(wavelength, length(particle_IDs)), nrow = length(particle_IDs), byrow=T)
     wave_shift = ((velocity_los / .speed_of_light) * wave) + wave # using doppler formula to compute the shift in wavelengths cause by LOS velocity
     luminosity = .interpolate_spectra(shifted_wave = wave_shift, spectra = intrinsic_spectra, wave_seq = observation$wave_seq)
     if (LSF_conv){
       luminosity = .lsf_convolution(observation=observation, luminosity=luminosity, lsf_sigma=lsf_sigma)
     }
     if (!is.na(observation$signal_to_noise) | observation$signal_to_noise == 0){
       luminosity = .add_noise(luminosity, observation$signal_to_noise)
     }
     spectra[i,] = (luminosity*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) # flux in units erg/s/cm^2/Ang
     if (verbose){cat(i, "... ", sep = "")}
  }

  cube = array(data = spectra, dim = c(observation$sbin, observation$sbin, observation$wave_bin))
  vel_image = array(data = vel_los, dim = c(observation$sbin, observation$sbin))
  dis_image = array(data = dis_los, dim = c(observation$sbin, observation$sbin))

  if (observation$psf_fwhm > 0){
    if (verbose){cat("Convolving cube with PSF... \n")    }
    cube = blur_datacube(cube = list("spectral_cube" = cube,
                                     "observation" = observation)) # apply psf convolution to each cube plane
  }

  if (verbose){cat("Done! \n")}

  if (write_fits){
    if (verbose){cat("Writing FITS... \n")}
    if (missing(output_location)){
      out_file_name = stringr::str_remove(simspin_file, ".Rdata")
      output_location = paste(out_file_name, "_inc", observation$inc_deg, "deg_seeing",
                              observation$psf_fwhm,"fwhm.FITS", sep="")
    }

    write_simspin_FITS(output_file = output_location, simspin_cube = cube, observation = observation)
  }

  output = list("spectral_cube"    = cube,
                "observation"      = observation,
                "velocity_image"   = vel_image,
                "dispersion_image" = dis_image)

  return(output)
}

