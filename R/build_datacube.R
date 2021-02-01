# Author: Kate Harborne
# Date: 26/10/2020
# Title: build_datacube - a function for generating a data cube from the observation
#
#'A function for making a mock spectral/velocity data cube
#'
#'The purpose of this function is to generate a mock IFU data cube from
#' a SimSpin file, as generated using the \code{make_simspin_file()} function.
#' The data cube produced can be either a spectral cube or a velocity cube.
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
#'@param cores Float describing the number of cores to run the interpolation
#' and velocity gridding on. Default is 1.
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
                          verbose = F, write_fits = F, output_location,
                          cores=1){

  if (verbose){cat("Computing observation parameters... \n")}
  observation = observation(telescope = telescope, observing_strategy = observing_strategy)

  # Reading in SimSpin file data
  if (typeof(simspin_file) == "character"){ # if provided with path to file
    simspin_data = readRDS(simspin_file)
  }
  if (typeof(simspin_file) == "list"){ # if provided with output list
    simspin_data = simspin_file
  }

  # Twisting galaxy about the z-axis to look from an angle
  twisted_data = twist_galaxy(simspin_data$star_part, twist_rad = observation$twist_rad)

  # Projecting the galaxy to given inclination
  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)

  if (verbose){cat("Assigning particles to spaxels... \n")}
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z_obs, breaks=observation$sbin_seq, labels=F)) - (observation$sbin)

  # Trimming particles that lie outside the aperture of the telescope
  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$pixel_region[!is.na(observation$pixel_region)],]



  if (observation$method == "spectral"){

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
    lum_map = array(data = NA, dim = observation$sbin^2)

    if (verbose){cat("Sorting spaxels... \n")}
    occupied = sort(unique(galaxy_data$pixel_pos))
    particle_IDs = .particles_to_pixels(galaxy_data, occupied, cores = cores)
    particles_in_spaxel = lapply(particle_IDs, length)

    if (verbose){cat("Generating spectra per spaxel... \n")}
    for (i in 1:length(particle_IDs)){ # computing the spectra at each occupied spatial pixel position
      if (particles_in_spaxel[[i]] > observation$particle_limit){ # if the number of particles in the spaxel is greater than the particle limit
        galaxy_sample = galaxy_data[particle_IDs[[i]],]
        intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]), nrow = particles_in_spaxel[[i]], byrow = T) * (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra
        vel_los[occupied[i]] = mean(galaxy_sample$vy_obs)
        dis_los[occupied[i]] = sd(galaxy_sample$vy_obs)
        wave = matrix(data = rep(wavelength, particles_in_spaxel[[i]]), nrow = particles_in_spaxel[[i]], byrow=T)
        wave_shift = ((galaxy_sample$vy_obs / .speed_of_light) * wave) + wave # using doppler formula to compute the shift in wavelengths cause by LOS velocity
        luminosity = .interpolate_spectra(shifted_wave = wave_shift, spectra = intrinsic_spectra, wave_seq = observation$wave_seq)
        if (LSF_conv){
          luminosity = .lsf_convolution(observation=observation, luminosity=luminosity, lsf_sigma=lsf_sigma)
        }
        if (!is.na(observation$signal_to_noise) | observation$signal_to_noise == 0){
          luminosity = .add_noise(luminosity, observation$signal_to_noise)
        }
        spectra[occupied[i],] = (luminosity*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) # flux in units erg/s/cm^2/Ang
        lum_map[occupied[i]] = ProSpect::bandpass(wave = observation$wave_seq,
                                                  flux = spectra[occupied[i],],
                                                  filter = observation$filter, flux_in = "wave", flux_out = "wave")
      }
      if (verbose){cat(i, "... ", sep = "")}
    }

    cube = array(data = spectra, dim = c(observation$sbin, observation$sbin, observation$wave_bin))
    vel_image = array(data = vel_los, dim = c(observation$sbin, observation$sbin))
    dis_image = array(data = dis_los, dim = c(observation$sbin, observation$sbin))
    lum_image = array(data = lum_map, dim = c(observation$sbin, observation$sbin))

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
                  "dispersion_image" = dis_image,
                  "flux_image"       = lum_image)

  }

  if (observation$method == "velocity"){

    observation$vbin = ceiling((max(abs(galaxy_data$vy_obs))*2) / observation$vbin_size) # the number of velocity bins in the cube
    if (observation$vbin <= 2){observation$vbin = 3}

    observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
    observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2

    vel_spec = matrix(data = 0, ncol = observation$vbin, nrow = observation$sbin^2)
    vel_los  = array(data = NA, dim = observation$sbin^2)
    dis_los  = array(data = NA, dim = observation$sbin^2)
    lum_map  = array(data = NA, dim = observation$sbin^2)

    if (verbose){cat("Sorting spaxels... \n")}
    occupied = sort(unique(galaxy_data$pixel_pos))
    particle_IDs = .particles_to_pixels(galaxy_data, occupied, cores = cores)
    particles_in_spaxel = lapply(particle_IDs, length)

    if (verbose){cat("Generating velocity distributions per spaxel... \n")}
    for (i in 1:length(particle_IDs)){
      if (particles_in_spaxel[[i]] > observation$particle_limit){ # if the number of particles in the spaxel is greater than the particle limit
        galaxy_sample = galaxy_data[particle_IDs[[i]],]

        intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]), nrow = particles_in_spaxel[[i]], byrow = T) * (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra
        spectral_flux = (intrinsic_spectra*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) # flux in units erg/s/cm^2/Ang

        galaxy_sample$luminosity = apply(spectral_flux, 1, ProSpect::bandpass,
                                         wave = simspin_data$wave,
                                         filter = observation$filter, flux_in = "wave", flux_out = "wave") #computing the r-band luminosity per particle from spectra

        vel_los[occupied[i]] = mean(galaxy_sample$vy_obs)
        dis_los[occupied[i]] = sd(galaxy_sample$vy_obs)
        lum_map[occupied[i]] = sum(galaxy_sample$luminosity)
        vel_spec[occupied[i],] = .sum_velocities(galaxy_sample = galaxy_sample, observation = observation, cores = 1)
      } # adding the "gaussians" of each particle to the velocity bins

      if (verbose){cat(occupied[i], "... ", sep = "")}

    }

    cube = array(data = vel_spec, dim = c(observation$sbin, observation$sbin, observation$vbin))
    vel_image = array(data = vel_los, dim = c(observation$sbin, observation$sbin))
    dis_image = array(data = dis_los, dim = c(observation$sbin, observation$sbin))
    lum_image = array(data = lum_map, dim = c(observation$sbin, observation$sbin))

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

    output = list("velocity_cube"    = cube,
                  "observation"      = observation,
                  "velocity_image"   = vel_image,
                  "dispersion_image" = dis_image,
                  "flux_image"       = lum_image)

  }

  return(output)
}

