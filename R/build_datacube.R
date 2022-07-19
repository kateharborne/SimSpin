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
#'@param method String to describe whether cubes output are "spectral", "gas",
#' "sf gas" or "velocity" (as in SimSpin v1) along the z-axis. Default is
#' "spectral".
#'@param verbose Default is \code{FALSE}. If you would like the code to give
#' updates about its progress, change this parameter to \code{TRUE}.
#'@param write_fits Default is \code{FALSE}. If you would like the code to
#' output a FITS file, change this parameter to \code{TRUE}.
#'@param output_location Optional parameter that describes the path and file
#' name of the FITS file output if \code{write_fits = TRUE}.
#' If \code{output_location} is specified as just a path, the file name will be
#' auto-generated based on the name of the input \code{simspin_file} and the
#' observing conditions and written to the specified directory.
#' If \code{write_fits = TRUE} and no \code{output_location} is specified, the
#' FITS file name will be auto-generated and written to the same directory as
#' the input \code{simspin_file}.
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
                          method, verbose = F, write_fits = F,
                          output_location, object_name="GalaxyID_unknown",
                          telescope_name="SimSpin",
                          observer_name="Anonymous",
                          split_save=F,
                          cores=1, mass_flag = F){

  if (missing(method)){
    if ("method" %in% names(telescope)){
      warning(">>> WARNING! >>> \n
              `method` is now specified within the build_datacube function directly,
              rather than within the telescope() class. \n
              Support for this input will remain in versions 2.X.X, but please consider
              updating your code.")
      method = telescope$method
    } else {
      method = "spectral"
    }
  }

  if (!missing(method) & ("method" %in% names(telescope))){
    warning(">>> WARNING >>> \n
            `method` has been specified in BOTH build_datacube() and telescope(). \n
            Using the `method` specified in build_datacube() and ignoring telescope(method). \n
            Please remove the `method` specified in telescope() to suppress this warning.")
  }

  method = stringr::str_to_lower(method)

  if (method != "spectral" &
      method != "velocity" &
      method != "gas" &
      method != "sf gas" ){
    stop("Error: Invalid method. \n Please specify method = 'spectral', 'velocity', 'gas' or 'sf gas' and try again.")
  }

  if (verbose){cat("Computing observation parameters... \n")}
  observation = observation(telescope = telescope, observing_strategy = observing_strategy, method = method)

  # Reading in SimSpin file data
  if (typeof(simspin_file) == "character"){ # if provided with path to file
    simspin_data = readRDS(simspin_file)
  }
  if (typeof(simspin_file) == "list"){ # if provided with output list
    simspin_data = simspin_file
    simspin_file = paste0("./", object_name)
  }

  if (length(simspin_data) == 4){
    # if we are working with a simspin file from before v2.3.0
    simspin_data$header = list("InputFile" = "Unknown",
                               "OutputFile" = simspin_file,
                               "Type" = "Unknown",
                               "Template" = "",
                               "Template_LSF" = "",
                               "Template_waveres" = "",
                               "Origin" = "< SimSpin v2.3.0",
                               "Date"   = "Unknown")

    if (length(simspin_data$wave) == 1221 | length(simspin_data$wave) == 842 ){
      simspin_data$header$Template = "BC03lr"
      simspin_data$header$Template_LSF = 3 # as according to Bruzual & Charlot (2003) MNRAS 344, pg 1000-1028
      simspin_data$header$Template_waveres = min(diff(simspin_data$wave))
    } else if (length(simspin_data$wave) == 6900 | length(simspin_data$wave) == 6521 ){
      simspin_data$header$Template = "BC03hr"
      simspin_data$header$Template_LSF = 3 # as according to Bruzual & Charlot (2003) MNRAS 344, pg 1000-1028
      simspin_data$header$Template_waveres = min(diff(simspin_data$wave))
    } else if (length(simspin_data$wave) == 53689 | length(simspin_data$wave) == 20356 ) {
      simspin_data$header$Template = "EMILES"
      simspin_data$header$Template_LSF = 2.51 # as according to Vazdekis et al (2016) MNRAS 463, pg 3409-3436
      simspin_data$header$Template_waveres = min(diff(simspin_data$wave))
    } else {
      stop("Error: Unknown spectral templates with no header information available. \n
           Please remake your SimSpin file with make_simspin_file > v2.3.0 to use custom templates.")
    }
    warning(cat("WARNING! - You are using an old SimSpin file (< v2.3.0). \n"))
    cat(paste0("Assuming that the ", simspin_data$header$Template, " has been used to build this SimSpin file. \n",
               "Consider re-making your SimSpin files using the make_simspin_file() function. \n"))
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

  if (length(galaxy_data$ID) == 0){
    stop(paste0("Error: There are no simulation particles within the aperture of the telescope. \n
         Please check that the method, `", method, "` is suitable for your input simulation file. \n
         Else, consider increasing your aperture size or adjusting the pointing of the telescope."))
  }

  if (verbose){cat("Sorting spaxels... \n")}

  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data[, list(val=list(ID), .N), by = "pixel_pos"]

  # SPECTRAL mode method =======================================================
  if (observation$method == "spectral"){

    # read original wavelengths of the template spectra
    wavelength = simspin_data$wave * (observation$z + 1) # and then applying a shift to those spectra due to redshift, z

    # If the requested wavelength resolution of the telescope is a smaller number than the intrinsic wavelength
    # resolution of the template spectra, the interpolation onto a finer grid can cause errors that pPXF cannot
    # account for in extreme cases. Issue a warning to the user in this case.

    if (observation$wave_res < min(diff(wavelength))){
      warning(cat("WARNING! - Wavelength resolution of provided template spectra at this redshift is too coarse for the requested telescope wavelength resolution.\n"))
      cat("Dlambda_telescope = ", observation$wave_res,  " A < Dlambda_templates ", min(diff(wavelength)), " A. \n")
      cat("This will cause some interpolation that may make spectral fitting techniques fail. \n")
      }

    # Similarly, the template spectra for each particle have some intrinsic spectral resolution
    # These vary dependent on the template. If the requested spectral resolution of the telescope
    # is lower than the spectral resolution of the templates, we can't convolve them.

    lsf_fwhm      = observation$lsf_fwhm
    lsf_fwhm_temp = simspin_data$header$Template_LSF * (observation$z + 1)
    # applying a shift to that intrinsic template LSF due to redshift, z

    spec_res_sigma_sq = ((lsf_fwhm^2) - (lsf_fwhm_temp^2))

    if (spec_res_sigma_sq < 0){ # if the lsf is smaller than the wavelength resolution of the spectra
      warning(cat("WARNING! - Spectral resolution of provided template spectra is greater than the requested telescope spectral resolution.\n"))
      cat("LSF_telescope = ", lsf_fwhm,  " A < LSF_templates (at redshift z) ", lsf_fwhm_temp, " A. \n")
      cat("No LSF convolution will be applied in this case. \n")
      cat("Intrinsic LSF of observation = ", lsf_fwhm_temp, " A for comparison with kinematic cubes. \n")
      observation$LSF_conv = FALSE
    } else {
      observation$LSF_conv = TRUE
      observation$lsf_sigma = (sqrt(spec_res_sigma_sq) / (2 * sqrt(2*log(2)))) / (simspin_data$header$Template_waveres * (1 + observation$z))
      # To get to the telescope's LSF, we only need to convolve with a Gaussian the width of the additional
      # difference between the redshifted template and the intrinsic telescope LSF.
      # This is the scaled for the wavelength pixel size at redshift "z".

      for (spectrum in 1:length(simspin_data$spectra)){
         convolved_spectrum = .lsf_convolution(observation, simspin_data$spectra[[spectrum]], observation$lsf_sigma)
         simspin_data$spectra[[spectrum]] = convolved_spectrum
      } # convolving the intrinsic spectra with the convolution kernel sized for the LSF

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

  # VELOCITY mode method =======================================================
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

      names(output$observed_images) = c("flux_image", "velocity_image", "dispersion_image", "h3_image", "h4_image") # default calling flux/mass as flux_image
      output$observed_images$flux_image       = array(0.0, dim = dims[c(1,2)])
      output$observed_images$velocity_image   = array(0.0, dim = dims[c(1,2)])
      output$observed_images$dispersion_image = array(0.0, dim = dims[c(1,2)])
      output$observed_images$h3_image         = array(0.0, dim = dims[c(1,2)])
      output$observed_images$h4_image         = array(0.0, dim = dims[c(1,2)])

      for (c in 1:dims[1]){
        for (d in 1:dims[2]){
          output$observed_images$flux_image[c,d]       = sum(output$velocity_cube[c,d,])
          output$observed_images$velocity_image[c,d]   = .meanwt(observation$vbin_seq, output$velocity_cube[c,d,])
          output$observed_images$dispersion_image[c,d] = sqrt(.varwt(observation$vbin_seq, output$velocity_cube[c,d,], output$observed_images$velocity_image[c,d]))
          h3h4 = optim(par = c(0,0),
                       fn = .losvd_fit,
                       x = observation$vbin_seq,
                       losvd = (output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T))),
                       vel = output$observed_images$velocity_image[c,d],
                       sig = output$observed_images$dispersion_image[c,d],
                       method="BFGS", control=list(reltol=1e-9))$par
          output$observed_images$h3_image[c,d]       = h3h4[1]
          output$observed_images$h4_image[c,d]       = h3h4[2]
        }
      }

      if (mass_flag){ # if mass flag is T, renaming the flux image as mass image
        names(output$observed_images)[which(names(output$observed_images) == "flux_image")] = "mass_image"
      }
    }

    if (verbose){cat("Done! \n")}

  }

  # GAS mode method =======================================================
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
    if (length(grep(".fits", output_location)) == 0 & length(grep(".FITS", output_location)) == 0 ){
      # if no filename has been specified, assume that the output location is just a path
      out_file_name = character(1)
      output_name = rev(stringr::str_split(simspin_file, "/")[[1]])[1]
      out_file_name = tryCatch({stringr::str_remove(output_name, ".Rdata")},
                               error = function(e){"./"})
      output_location = paste(output_location, "/", out_file_name, "_inc", observation$inc_deg, "deg_seeing",
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

