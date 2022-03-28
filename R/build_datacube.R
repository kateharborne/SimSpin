# Author: Kate Harborne
# Date: 27/08/21
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
#' resolution, wavelength resolution, etc.). See
#' \code{\link{telescope}} help for more details.
#'@param observing_strategy An object of the observing_strategy class that
#' describes the properties of the observed simulation (i.e. redshift,
#' inclination, seeing conditions). See \code{\link{observing_strategy}}
#' help for more details.
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
#'@param split_save Boolean describing whether to split the output from
#' \code{build_datacube()} into multiple files while saving to FITS. If TRUE,
#' several FITS files will be saved with file names that reflect their content
#' (i.e."_spectral_cube.FITS", "_velocity_image.FITS", "_flux_images.FITS",
#' etc.). Default option is FALSE.
#'@param cores Float describing the number of cores to run the interpolation
#' and velocity gridding on. Default is 1.
#'@param mass_flag Boolean flag that, when set to TRUE, will compute properties
#' using a mass weighting rather than a luminosity weighting. Default is FALSE.
#'@return Returns a list containing four elements:
#'\enumerate{
#' \item \code{spectral_cube} or \code{velocity_cube} - a 3D array containing
#' either a spectral cube or a velocity cube, where the output type is
#' determined by the \code{method} selected in the \code{telescope} function.
#' \item \code{observation} - a list containing a summary of the details of the
#' observation (i.e. the output from the function \code{observation()}).
#' \item \code{raw_images} - a list of 2D arrays, where each 2D array represents
#' a raw particle image gridded as for the observation details.
#' \item \code{observed_images} - NULL or a list of 2D arrays (again, where the
#' output type is determined by the \code{method} selected in the
#' \code{telescope} function.) containing kinematic images of the collapsed
#' cube. If \code{blur=T}, these images will be blurred to the specified amount.
#'}
#' If \code{write_fits = T}, a .fits file that contains a the generated cube and
#' relevant header describing the mock observation will also be produced at the
#' specified \code{output_location}.
#'
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
                          split_save=F,
                          cores=1, mass_flag = F){

  if (verbose){cat("Computing observation parameters... \n")}
  observation = observation(telescope = telescope, observing_strategy = observing_strategy)

  # Reading in SimSpin file data
  if (typeof(simspin_file) == "character"){ # if provided with path to file
    simspin_data = readRDS(simspin_file)
  }
  if (typeof(simspin_file) == "list"){ # if provided with output list
    simspin_data = simspin_file
    simspin_file = paste0("./", object_name)
  }

  if (observation$method == "spectral" | observation$method == "velocity"){
    galaxy_data = simspin_data$star_part
  } else if (observation$method == "sf gas"){
    galaxy_data = simspin_data$gas_part[simspin_data$gas_part$SFR > 0,]
  } else if (observation$method == "gas"){
    galaxy_data = simspin_data$gas_part
  } else {
    stop("Error: Invalid method. \n Please specify observation$method = 'spectral', 'velocity', 'sf gas', or 'gas' and try again.")
  }

  if (!data.table::is.data.table(galaxy_data)){
    # if we are working with a simspin file from before v2.1.5
    warning(cat("WARNING! - You are using an old SimSpin file (< v2.1.5). \n"))
    cat("For quicker processing, consider re-making your SimSpin files using the make_simspin_file() function. \n")
    # convert from a data.frame() to a data.table.
    galaxy_data = data.table::as.data.table(galaxy_data)
    simspin_data$spectra = data.table::as.data.table(simspin_data$spectra)
  }

  # Twisting galaxy about the z-axis to look from an angle
  twisted_data = twist_galaxy(galaxy_data, twist_rad = observation$twist_rad)

  # Projecting the galaxy to given inclination
  obs_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)

  galaxy_data$x = (obs_data$x + observation$pointing_kpc[1])   # adjusting pointing of the aperture by x_kpc
  galaxy_data$y = obs_data$y                                   #   and y_kpc
  galaxy_data$z = (obs_data$z + observation$pointing_kpc[2])
  galaxy_data$vx = obs_data$vx; galaxy_data$vy = obs_data$vy; galaxy_data$vz = obs_data$vz

  remove(obs_data, twisted_data)

  if (verbose){cat("Assigning particles to spaxels... \n")}
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z, breaks=observation$sbin_seq, labels=F)) - (observation$sbin)

  # Trimming particles that lie outside the aperture of the telescope
  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$pixel_region[!is.na(observation$pixel_region)],]

  if (verbose){cat("Sorting spaxels... \n")}

  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data[, list(val=list(ID), .N), by = "pixel_pos"]

  if (observation$method == "spectral"){

    wavelength = (observation$z * simspin_data$wave) + simspin_data$wave # applying a shift due to redshift, z, to original wavelength
    lsf_fwhm   = observation$lsf_fwhm
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
    raw_images = list(
      flux_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin)),
      particle_image = array(data = output[[5]], dim = c(observation$sbin, observation$sbin))
      )

    output = list("spectral_cube"    = cube,
                  "observation"      = observation,
                  "raw_images"       = raw_images,
                  "observed_images"  = NULL)

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      output = blur_datacube(output) # apply psf convolution to each cube plane
    }

    if (verbose){cat("Done! \n")}

  }

  if (observation$method == "velocity"){

    observation$vbin = ceiling((max(abs(galaxy_data$vy))*2) / observation$vbin_size) # the number of velocity bins in the cube
    if (observation$vbin <= 2){observation$vbin = 3}

    observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
    observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2

    if (verbose){cat("Generating stellar velocity distributions per spaxel... \n")}
    if (cores == 1){
      output = .velocity_spaxels(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, mass_flag)
    }
    if (cores > 1){
      output = .velocity_spaxels_mc(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, cores, mass_flag)
    }

    cube = array(data = output[[1]], dim = c(observation$sbin, observation$sbin, observation$vbin))
    raw_images = list(
      flux_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin)),
      age_image = array(data = output[[5]], dim = c(observation$sbin, observation$sbin)),
      metallicity_image = array(data = output[[6]], dim = c(observation$sbin, observation$sbin)),
      particle_image = array(data = output[[7]], dim = c(observation$sbin, observation$sbin))
      )

    output = list("velocity_cube"   = cube,
                  "observation"     = observation,
                  "raw_images"      = raw_images,
                  "observed_images"  = vector(mode = "list", length=3))

    if (mass_flag){ # if mass flag is T, the flux image is really just a mass image
      names(output$raw_images)[which(names(output$raw_images) == "flux_image")] = "mass_image"
    }

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      output = blur_datacube(output) # apply psf convolution to each cube plane
    } else {
      dims = dim(output$raw_images$velocity_image)

      names(output$observed_images) = c("flux_image", "velocity_image", "dispersion_image") # default calling flux/mass as flux_image
      output$observed_images$flux_image       = array(0.0, dim = dims[c(1,2)])
      output$observed_images$velocity_image   = array(0.0, dim = dims[c(1,2)])
      output$observed_images$dispersion_image = array(0.0, dim = dims[c(1,2)])

      for (c in 1:dims[1]){
        for (d in 1:dims[2]){
          output$observed_images$flux_image[c,d]       = sum(output$velocity_cube[c,d,])
          output$observed_images$velocity_image[c,d]   = .meanwt(observation$vbin_seq, output$velocity_cube[c,d,])
          output$observed_images$dispersion_image[c,d] = sqrt(.varwt(observation$vbin_seq, output$velocity_cube[c,d,], output$observed_images$velocity_image[c,d]))
        }
      }

      if (mass_flag){ # if mass flag is T, renaming the flux image as mass image
        names(output$observed_images)[which(names(output$observed_images) == "flux_image")] = "mass_image"
      }
    }

    if (verbose){cat("Done! \n")}

  }

  if (observation$method == "gas" | observation$method == "sf gas"){

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
    raw_images = list(
      mass_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin)),
      SFR_image = array(data = output[[5]], dim = c(observation$sbin, observation$sbin)),
      metallicity_image = array(data = output[[6]], dim = c(observation$sbin, observation$sbin)),
      OH_image  = array(data = output[[7]], dim = c(observation$sbin, observation$sbin))
      )

    output = list("velocity_cube"   = cube,
                  "observation"     = observation,
                  "raw_images"      = raw_images,
                  "observed_images" = vector(mode = "list", length=3))

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      output = blur_datacube(output) # apply psf convolution to each cube plane
    } else {
      dims = dim(output$raw_images$velocity_image)

      names(output$observed_images) = c("mass_image", "velocity_image", "dispersion_image")
      output$observed_images$mass_image       = array(0.0, dim = dims[c(1,2)])
      output$observed_images$velocity_image   = array(0.0, dim = dims[c(1,2)])
      output$observed_images$dispersion_image = array(0.0, dim = dims[c(1,2)])

      for (c in 1:dims[1]){
        for (d in 1:dims[2]){
          output$observed_images$mass_image[c,d]       = sum(output$velocity_cube[c,d,])
          output$observed_images$velocity_image[c,d]   = .meanwt(observation$vbin_seq, output$velocity_cube[c,d,])
          output$observed_images$dispersion_image[c,d] = sqrt(.varwt(observation$vbin_seq, output$velocity_cube[c,d,], output$observed_images$velocity_image[c,d]))
        }
      }
    }

    if (verbose){cat("Done! \n")}

  }

  # Trimming off extra zeros from images outside the aperture of the telescope
  aperture_region = matrix(data = observation$aperture_region, nrow = observation$sbin, ncol = observation$sbin)

  for (raw_image in names(output$raw_images)){
    output$raw_images[[raw_image]] = output$raw_images[[raw_image]] * aperture_region
  }
  for (obs_image in names(output$observed_images)){
    output$observed_images[[obs_image]] = output$observed_images[[obs_image]] * aperture_region
  }

  if (write_fits){
    if (verbose){cat("Writing FITS... \n")}
    if (missing(output_location)){
      out_file_name = character(1)
      out_file_name = tryCatch({stringr::str_remove(simspin_file, ".Rdata")},
                               error = function(e){"./"})
      output_location = paste(out_file_name, "_inc", observation$inc_deg, "deg_seeing",
                              observation$psf_fwhm,"fwhm.FITS", sep="")
    }

    write_simspin_FITS(output_file = output_location,
                       simspin_datacube = output, object_name = object_name,
                       telescope_name = telescope_name, instrument_name = telescope$type,
                       observer_name = observer_name, split_save=split_save,
                       input_simspin_file = rev(stringr::str_split(simspin_file, "/")[[1]])[1])
  }

  return(output)
}

