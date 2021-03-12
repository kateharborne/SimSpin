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
#'@param object_name Optional string used in \code{write_simspin_FITS} to
#' describe the name of the object observed in FITS header.
#'@param telescope_name Optional string used in \code{write_simspin_FITS} to
#' describe the name of the telescope used for observation in FITS header.
#'@param observer_name Optional string used in \code{write_simspin_FITS} to
#' describe the name of the person who ran the observation in FITS header.
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
                          object_name="GalaxyID_unknown",
                          telescope_name="SimSpin",
                          observer_name="Anonymous",
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

  if (observation$method == "spectral" | observation$method == "velocity"){
    galaxy_data = simspin_data$star_part
  } else {
    galaxy_data = simspin_data$gas_part[simspin_data$gas_part$SFR > 0,]
  }

  # Twisting galaxy about the z-axis to look from an angle
  twisted_data = twist_galaxy(galaxy_data, twist_rad = observation$twist_rad)

  # Projecting the galaxy to given inclination
  obs_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)

  galaxy_data$x = obs_data$x;   galaxy_data$y = obs_data$y;   galaxy_data$z = obs_data$z
  galaxy_data$vx = obs_data$vx; galaxy_data$vy = obs_data$vy; galaxy_data$vz = obs_data$vz

  if (verbose){cat("Assigning particles to spaxels... \n")}
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z, breaks=observation$sbin_seq, labels=F)) - (observation$sbin)

  # Trimming particles that lie outside the aperture of the telescope
  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$pixel_region[!is.na(observation$pixel_region)],]

  if (verbose){cat("Sorting spaxels... \n")}

  particle_ID = NULL # initiallising varible to avoid CRAN error in checks

  galaxy_data_table = data.table::data.table("particle_ID" = seq(1, length(galaxy_data$x)), "spaxel_ID"=galaxy_data$pixel_pos)

  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data_table[, list(val=list(particle_ID)), by = "spaxel_ID"]

  if (observation$method == "spectral"){

    original_wave  = simspin_data$wave # read original wavelengths
    wavelength = (observation$z * original_wave) + original_wave # and then applying a shift due to redshift, z
    lsf_fwhm   = (observation$z * observation$lsf_fwhm) + observation$lsf_fwhm # adjusting the LSF for the resolution at z

    spec_res_sigma_sq = lsf_fwhm^2 - (min(diff(wavelength)))^2
    if (spec_res_sigma_sq < 0){ # if the lsf is smaller than the wavelength resolution of the spectra
      warning(cat("WARNING! - Wavelength resolution of provided spectra is lower than the requested telescope resolution.\n"))
      cat("LSF = ", observation$lsf_fwhm,  " A < wavelength resolution ", min(diff(wavelength)), " A. \n")
      cat("No LSF will be applied in this case.\n")
      observation$LSF_conv = FALSE
    } else {
      observation$LSF_conv = TRUE
      observation$lsf_sigma = (sqrt(spec_res_sigma_sq) / (2 * sqrt(2*log(2))))
    }

    if (verbose){cat("Generating spectra per spaxel... \n")}

    if (cores == 1){
      output = .spectral_spaxels(part_in_spaxel, wavelength, observation, galaxy_data, simspin_data, verbose)
    }
    if (cores > 1){
      output = .spectral_spaxels_mc(part_in_spaxel, wavelength, observation, galaxy_data, simspin_data, verbose, cores)
    }

    cube = array(data = output[[1]], dim = c(observation$sbin, observation$sbin, observation$wave_bin))
    lum_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin))
    vel_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin))
    dis_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin))

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      cube = blur_datacube(cube = list("spectral_cube" = cube,
                                       "observation" = observation)) # apply psf convolution to each cube plane
    }

    if (verbose){cat("Done! \n")}

    output = list("spectral_cube"    = cube,
                  "observation"      = observation,
                  "velocity_image"   = vel_image,
                  "dispersion_image" = dis_image,
                  "flux_image"       = lum_image)

    if (write_fits){
      if (verbose){cat("Writing FITS... \n")}
      if (missing(output_location)){
        tryCatch(out_file_name = {stringr::str_remove(simspin_file, ".Rdata")},
                 error = function(e){"./"})
        output_location = paste(out_file_name, "_inc", observation$inc_deg, "deg_seeing",
                                observation$psf_fwhm,"fwhm_spectral.FITS", sep="")
      }

      write_simspin_FITS(output_file = output_location,
                         simspin_data = output, object_name = object_name,
                         telescope_name = telescope_name, instrument_name = telescope$type,
                         observer_name = observer_name,
                         input_simspin_file = rev(stringr::str_split(simspin_file, "/")[[1]])[1])
    }

  }

  if (observation$method == "velocity"){

    observation$vbin = ceiling((max(abs(galaxy_data$vy))*2) / observation$vbin_size) # the number of velocity bins in the cube
    if (observation$vbin <= 2){observation$vbin = 3}

    observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
    observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2

    if (verbose){cat("Generating stellar velocity distributions per spaxel... \n")}
    if (cores == 1){
      output = .velocity_spaxels(part_in_spaxel, observation, galaxy_data, simspin_data, verbose)
    }
    if (cores > 1){
      output = .velocity_spaxels_mc(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, cores)
    }

    cube = array(data = output[[1]], dim = c(observation$sbin, observation$sbin, observation$vbin))
    lum_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin))
    vel_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin))
    dis_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin))
    age_image = array(data = output[[5]], dim = c(observation$sbin, observation$sbin))
    met_image = array(data = output[[6]], dim = c(observation$sbin, observation$sbin))

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      cube = blur_datacube(cube = list("spectral_cube" = cube,
                                       "observation" = observation)) # apply psf convolution to each cube plane
    }

    if (verbose){cat("Done! \n")}

    output = list("velocity_cube"    = cube,
                  "observation"      = observation,
                  "velocity_image"   = vel_image,
                  "dispersion_image" = dis_image,
                  "flux_image"       = lum_image,
                  "age_image"        = age_image,
                  "metallicity_image"= met_image)

    if (write_fits){
      if (verbose){cat("Writing FITS... \n")}
      if (missing(output_location)){
        tryCatch(out_file_name = {stringr::str_remove(simspin_file, ".Rdata")},
                 error = function(e){"./"})
        output_location = paste(out_file_name, "_inc", observation$inc_deg, "deg_seeing",
                                observation$psf_fwhm,"fwhm_images.FITS", sep="")
      }

      write_simspin_FITS(output_file = output_location,
                         simspin_data = output, object_name = object_name,
                         telescope_name = telescope_name , instrument_name = telescope$type,
                         observer_name = observer_name,
                         input_simspin_file = rev(stringr::str_split(simspin_file, "/")[[1]])[1])
    }

  }

  if (observation$method == "gas"){

    observation$vbin = ceiling((max(abs(galaxy_data$vy))*2) / observation$vbin_size) # the number of velocity bins in the cube
    if (observation$vbin <= 2){observation$vbin = 3}

    observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
    observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2

    if (verbose){cat("Generating gas velocity distributions per spaxel... \n")}
    if (cores == 1){
      output = .gas_velocity_spaxels(part_in_spaxel, observation, galaxy_data, simspin_data, verbose)
    }
    if (cores > 1){
      output = .gas_velocity_spaxels_mc(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, cores)
    }

    cube = array(data = output[[1]], dim = c(observation$sbin, observation$sbin, observation$vbin))
    mass_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin))
    vel_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin))
    dis_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin))
    met_image = array(data = output[[5]], dim = c(observation$sbin, observation$sbin))
    OH_image  = array(data = output[[6]], dim = c(observation$sbin, observation$sbin))

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      cube = blur_datacube(cube = list("spectral_cube" = cube,
                                       "observation" = observation)) # apply psf convolution to each cube plane
    }

    if (verbose){cat("Done! \n")}

    output = list("velocity_cube"    = cube,
                  "observation"      = observation,
                  "velocity_image"   = vel_image,
                  "dispersion_image" = dis_image,
                  "mass_image"       = mass_image,
                  "metallicity_image"= met_image,
                  "OH_image"         = OH_image)

    if (write_fits){
      if (verbose){cat("Writing FITS... \n")}
      if (missing(output_location)){
        tryCatch(out_file_name = {stringr::str_remove(simspin_file, ".Rdata")},
                 error = function(e){"./"})
        output_location = paste(out_file_name, "_inc", observation$inc_deg, "deg_seeing",
                                observation$psf_fwhm,"fwhm_images.FITS", sep="")
      }

      write_simspin_FITS(output_file = output_location,
                         simspin_data = output, object_name = object_name,
                         telescope_name = telescope_name , instrument_name = telescope$type,
                         observer_name = observer_name,
                         input_simspin_file = rev(stringr::str_split(simspin_file, "/")[[1]])[1])
    }

  }

  return(output)
}

