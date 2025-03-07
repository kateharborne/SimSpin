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
#'@param moments Integer "2" or "4". Default is "4". When method="velocity",
#'  "gas" or "sf_gas", this parameter dictates if the kinematics are fit with
#'  a Gaussian (moments="2") or Gauss-Hermite polynomial (moments="4").
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
#'@param voronoi_bin Boolean flag that, when set to TRUE, will bin pixels into
#' voronoi tessellated cells that contain a minimum number of particles per
#' pixel, specified by \code{vorbin_limit}. Default is FALSE.
#'@param vorbin_limit Integer float that describes the minimum number of
#' particles per pixel within a given bin, only used if \code{voronoi_bin = T}.
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
                          method, moments = 4, verbose = F, write_fits = F,
                          output_location, object_name="GalaxyID_unknown",
                          telescope_name="SimSpin",
                          observer_name="Anonymous",
                          split_save=F,
                          cores=1, mass_flag = F,
                          voronoi_bin=F, vorbin_limit=10){

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

  if (moments != 2 & moments != 4){
    stop("ERROR: Invalid `moments` requested. Please specify `2` or `4` for the Gaussian or Gauss-Hermite polynomial fit.")
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

  if (!"header" %in% names(simspin_data)){
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
    if (length(galaxy_data)==0){
      stop(c("Error: No stellar particles exist in this SimSpin file. \n",
             "Please specify a different method ('gas' or 'sf gas') and try again. \n"))
    }
  } else if (observation$method == "sf gas"){
    galaxy_data = simspin_data$gas_part[simspin_data$gas_part$SFR > 0,]
    if (length(galaxy_data)==0){
      stop(c("Error: No star forming gas particles exist in this SimSpin file. \n",
             "Please specify a different method ('gas', 'velocity' or 'spectral') and try again. \n"))
    }
  } else if (observation$method == "gas"){
    galaxy_data = simspin_data$gas_part
    if (length(galaxy_data)==0){
      stop(c("Error: No gas particles exist in this SimSpin file. \n",
             "Please specify a different method ('velocity' or 'spectral') and try again. \n"))
    }
  }

  if (!data.table::is.data.table(galaxy_data)){
    # if we are working with a simspin file from before v2.1.5
    warning(cat("WARNING! - You are using an old SimSpin file (< v2.1.5). \n"))
    cat("For quicker processing, consider re-making your SimSpin files using the make_simspin_file() function. \n")
    # convert from a data.frame() to a data.table.
    galaxy_data = data.table::as.data.table(galaxy_data)
    simspin_data$spectra = data.table::as.data.table(simspin_data$spectra)
  }

  if (!"spectral_weights" %in% names(simspin_data)){
    # if we are working with a simspin file from before v2.6.0
    warning(cat("WARNING! - You are using an old SimSpin file (< v2.6.0). \n"))
    cat("For smaller SimSpin file sizes, consider re-making your SimSpin files using the make_simspin_file() function. \n")
    spectra_flag = 1
  } else {spectra_flag = 0}

  if (simspin_data$header$Template == "EMILES"){
      temp = SimSpin::EMILES
    } else if (simspin_data$header$Template == "BC03hr"){
      temp = SimSpin::BC03hr
    } else {
      temp = SimSpin::BC03lr
    }

  observation$moments = moments

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

  # function for adding flux information to stellar methods
  if (observation$method == "spectral" | observation$method == "velocity"){
    if (cores > 1){
      galaxy_data = .compute_flux_mc(observation, galaxy_data, simspin_data,
                                  template=temp, verbose, spectra_flag, cores)
    } else {
      galaxy_data = .compute_flux(observation, galaxy_data, simspin_data,
                                  template=temp, verbose, spectra_flag)
    }
  }

  if (length(galaxy_data$ID) == 0){
    stop(paste0("Error: There are no simulation particles within the aperture of the telescope. \n
         Please check that the method, `", method, "` is suitable for your input simulation file. \n
         Else, consider increasing your aperture size or adjusting the pointing of the telescope."))
  }

  if (verbose){cat("Sorting spaxels... \n")}

  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data[, list(val=list(ID), .N), by = "pixel_pos"]


  if (method == "spectral" | method == "velocity"){
    summed_images = galaxy_data[, list(.N,
                                       luminosity = sum(luminosity),
                                       filter_luminosity = sum(filter_luminosity),
                                       mass = sum(Mass)),
                                by = "pixel_pos"]


    empty_pixels = data.table::data.table("pixel_pos" = which(!(seq(1:(observation$sbin*observation$sbin))
                                                                %in% summed_images$pixel_pos)),
                                          "N" = 0,
                                          "luminosity" = 0.,
                                          "filter_luminosity" = 0.,
                                          "mass" = 0.)

    summed_images = data.table::rbindlist(list(summed_images, empty_pixels))
    summed_images = summed_images[order(pixel_pos)]
    remove(empty_pixels)

  } else if (method == "gas" | method == "sf gas"){
    summed_images = galaxy_data[, list(.N,
                                       sfr = sum(SFR),
                                       mass = sum(Mass)),
                                by = "pixel_pos"]

    empty_pixels = data.table::data.table("pixel_pos" = which(!(seq(1:(observation$sbin*observation$sbin))
                                                                %in% summed_images$pixel_pos)),
                                          "N" = 0,
                                          "sfr" = 0.,
                                          "mass" = 0.)

    summed_images = data.table::rbindlist(list(summed_images, empty_pixels))
    summed_images = summed_images[order(pixel_pos)]
    remove(empty_pixels)
  }

  if (voronoi_bin){ # Returning the binned pixels based on some voronoi limit
    if (verbose){cat("Binning spaxels into voronoi bins... \n")}

    part_in_spaxel = voronoi(part_in_spaxel=part_in_spaxel, obs=observation,
                             particle_limit = as.numeric(vorbin_limit),
                             roundness_limit = 0.3, uniform_limit = 0.8,
                             verbose = verbose)

    observation$particle_limit = as.numeric(vorbin_limit)

  }


  # SPECTRAL mode method =======================================================
  if (observation$method == "spectral"){

    # read original wavelengths of the template spectra
    wavelength = temp$Wave * (observation$z + 1) # and then applying a shift to those spectra due to redshift, z

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

    spec_res_fwhm_sq = ((lsf_fwhm^2) - (lsf_fwhm_temp^2))

    if (spec_res_fwhm_sq < 0){ # if the lsf is smaller than the wavelength resolution of the spectra
      warning(cat("WARNING! - Spectral resolution of provided template spectra is greater than the requested telescope spectral resolution.\n"))
      cat("LSF_telescope = ", lsf_fwhm,  " A < LSF_templates (at redshift z) ", lsf_fwhm_temp, " A. \n")
      cat("No LSF convolution will be applied in this case. \n")
      observation$LSF_conv = FALSE
    } else {
      observation$LSF_conv = TRUE
      observation$lsf_sigma = (sqrt(spec_res_fwhm_sq) / (2 * sqrt(2*log(2))) / observation$wave_res)
      # To get to the telescope's LSF, we only need to convolve with a Gaussian the width of the additional
      # difference between the redshifted template and the intrinsic telescope LSF.
      # This is the scaled for the wavelength pixel size of the observation.
    }

    if (verbose){cat("Generating spectra per spaxel... \n")}

    if (cores == 1){
      output = .spectral_spaxels(part_in_spaxel, wavelength, observation, galaxy_data, simspin_data, temp, verbose, spectra_flag)
    }
    if (cores > 1){
      output = .spectral_spaxels_mc(part_in_spaxel, wavelength, observation, galaxy_data, simspin_data, temp, verbose, cores, spectra_flag)
    }

    cube = array(data = output[[1]], dim = c(observation$sbin, observation$sbin, observation$wave_bin))
    raw_images = list(
      flux_image = array(data = summed_images$luminosity, dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin)),
      age_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin)),
      metallicity_image = array(data = output[[5]], dim = c(observation$sbin, observation$sbin)),
      particle_image = array(data = summed_images$N, dim = c(observation$sbin, observation$sbin)),
      voronoi_bins = array(data = output[[6]], dim = c(observation$sbin, observation$sbin)),
      mass_image = array(data = summed_images$mass, dim = c(observation$sbin, observation$sbin))
      )

    output = list("spectral_cube"    = cube,
                  "observation"      = observation,
                  "raw_images"       = raw_images,
                  "observed_images"  = NULL,
                  "variance_cube"    = NULL)

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      output = blur_datacube(output) # apply psf convolution to each cube plane
    }

    if (observation$LSF_conv){
      if (verbose){cat("Convolving spectra with LSF... \n")    }
      for (a in 1:dim(output$spectral_cube)[1]){
        for (b in 1:dim(output$spectral_cube)[2]){
          output$spectral_cube[a,b,] = .lsf_convolution(observation,
                                                        output$spectral_cube[a,b,],
                                                        observation$lsf_sigma)
        }
      }
    }

    if (!is.na(observation$signal_to_noise)){ # should we add noise?
      if (verbose){cat("Adding noise... \n")    }
      noise_cube = array(data = 0, dim = dim(output$spectral_cube))
      output$variance_cube = noise_cube # initialising empty arrays

      noise_cube = .add_noise(output$spectral_cube,
                              sqrt(median(raw_images$flux_image[raw_images$particle_image > 0], na.rm=T))/
                                (observation$signal_to_noise*sqrt(raw_images$flux_image)))

      output$spectral_cube = output$spectral_cube + noise_cube
      output$variance_cube = 1/(noise_cube)^2
      output$variance_cube[is.infinite(output$variance_cube)] = 0
      remove(noise_cube)
    }

    if (verbose){cat("Done! \n")}

  }

  # VELOCITY mode method =======================================================
  if (observation$method == "velocity"){

    observation$vbin = 5*ceiling((max(abs(galaxy_data$vy))*2) / observation$vbin_size) # the number of velocity bins in the cube
    if (observation$vbin <= 2){observation$vbin = 3}

    observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
    observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2
    if (mass_flag){observation$mass_flag = TRUE}

    if (verbose){cat("Generating stellar velocity distributions per spaxel... \n")}
    if (cores == 1){
      output = .velocity_spaxels(part_in_spaxel, observation, galaxy_data, simspin_data, temp, verbose, mass_flag, spectra_flag)
    }
    if (cores > 1){
      output = .velocity_spaxels_mc(part_in_spaxel, observation, galaxy_data, simspin_data, temp, verbose, cores, mass_flag, spectra_flag)
    }

    cube = array(data = output[[1]], dim = c(observation$sbin, observation$sbin, observation$vbin))
    raw_images = list(
      flux_image = array(data = summed_images$luminosity, dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin)),
      age_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin)),
      metallicity_image = array(data = output[[5]], dim = c(observation$sbin, observation$sbin)),
      mass_image = array(data = summed_images$mass, dim = c(observation$sbin, observation$sbin)),
      particle_image = array(data = summed_images$N, dim = c(observation$sbin, observation$sbin)),
      voronoi_bins = array(data = output[[6]], dim = c(observation$sbin, observation$sbin))
      )
    observed_images = list(
      flux_image = array(data = summed_images$filter_luminosity, dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      h3_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      h4_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      residuals = array(0.0, dim = c(observation$sbin, observation$sbin)),
      mass_image = array(data = summed_images$mass, dim = c(observation$sbin, observation$sbin))
      )

    output = list("velocity_cube"   = cube,
                  "observation"     = observation,
                  "raw_images"      = raw_images,
                  "observed_images" = observed_images,
                  "variance_cube"   = NULL)

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      output = blur_datacube(output) # apply psf convolution to each cube plane
    }

    if (!is.na(observation$signal_to_noise)){ # should we add noise?
      if (verbose){cat("Adding noise... \n")    }
      noise_cube = array(data = 0, dim = dim(output$velocity_cube))
      output$variance_cube = noise_cube # initialising empty arrays

      noise_cube = .add_noise(output$velocity_cube,
                              sqrt(median(raw_images$flux_image[raw_images$particle_image > 0], na.rm=T))/
                                (observation$signal_to_noise*sqrt(raw_images$flux_image)))
      noise_image = output$observed_images$flux_image*(rowSums(noise_cube, dims=2)/rowSums(output$velocity_cube, dims=2))

      output$velocity_cube = output$velocity_cube + noise_cube
      output$observed_images$flux_image = output$observed_images$flux_image + noise_image

      if (any(output$velocity_cube<0)){
        output$velocity_cube[which(output$velocity_cube<0)] = 0
      } # removing negative values introduced at the edges of the LOSVD
      output$variance_cube = 1/(noise_cube)^2
      output$variance_cube[is.infinite(output$variance_cube)] = 0
      remove(noise_cube)
    }

    dims = dim(output$raw_images$velocity_image)

    if (moments == 4){
      for (c in 1:dims[1]){
        for (d in 1:dims[2]){
          vel_ini = .meanwt(observation$vbin_seq, output$velocity_cube[c,d,])
          sd_ini  = sqrt(.varwt(observation$vbin_seq, output$velocity_cube[c,d,], vel_ini))

          kin   = tryCatch({stats::optim(par   = c(vel_ini,sd_ini,0,0),
                                         fn    = .losvd_fit_h3h4,
                                         x     = observation$vbin_seq,
                                         losvd = (output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T))),
                                         method="L-BFGS-B", lower = c(NA, 1e-10, NA, NA),
                                         control=list(pgtol=1e-9))$par},
                           error = function(e){c(0,0,0,0)})

          output$observed_images$velocity_image[c,d]   = kin[1]
          output$observed_images$dispersion_image[c,d] = kin[2]
          output$observed_images$h3_image[c,d]         = kin[3]
          output$observed_images$h4_image[c,d]         = kin[4]
          output$observed_images$residuals[c,d]        = mean(abs(.losvd_out_h3h4(x=observation$vbin_seq, vel=kin[1], sig=kin[2], h3=kin[3], h4=kin[4]) -
                                                                    output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T)) ), na.rm=T)
        }
      }
    }

    if (moments == 2){
      for (c in 1:dims[1]){
        for (d in 1:dims[2]){
          vel_ini = .meanwt(observation$vbin_seq, output$velocity_cube[c,d,])
          sd_ini  = sqrt(.varwt(observation$vbin_seq, output$velocity_cube[c,d,], vel_ini))

          kin   = tryCatch({stats::optim(par   = c(vel_ini,sd_ini),
                                         fn    = .losvd_fit_vsig,
                                         x     = observation$vbin_seq,
                                         losvd = (output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T))),
                                         method="L-BFGS-B", lower = c(NA, 1e-10),
                                         control=list(pgtol=1e-9))$par},
                           error = function(e){c(0,0)})

          output$observed_images$velocity_image[c,d]   = kin[1]
          output$observed_images$dispersion_image[c,d] = kin[2]
          output$observed_images$h3_image[c,d]         = 0
          output$observed_images$h4_image[c,d]         = 0
          output$observed_images$residuals[c,d]        = mean(abs(.losvd_out_vsig(x=observation$vbin_seq, vel=kin[1], sig=kin[2]) -
                                                                    output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T)) ), na.rm=T)
        }
      }
    }

    if (verbose){cat("Done! \n")}

  }

  # GAS mode method =======================================================
  if (observation$method == "gas" | observation$method == "sf gas"){

    # Check if SimSpin file is prior to version 2.3.15, re-format SFR units
    if (as.numeric(stringr::str_split(simspin_data$header$Origin, pattern = "2.")[[1]][2]) < 3.16){
      warning(c("In SimSpin files built with < v2.3.16, gas star formation rates are stored in g/s. \n",
                "Re-formatting to display SFR in units of Msol/yr."))
      galaxy_data$SFR = galaxy_data$SFR*(.g_to_msol/.s_to_yr)
    }

    observation$vbin = 5*ceiling((max(abs(galaxy_data$vy))*2) / observation$vbin_size) # the number of velocity bins in the cube
    if (observation$vbin <= 2){observation$vbin = 3}

    observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
    observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2

    observation$mass_flag = TRUE

    if (verbose){cat("Generating gas velocity distributions per spaxel... \n")}
    if (cores == 1){
      output = .gas_velocity_spaxels(part_in_spaxel, observation, galaxy_data, simspin_data, verbose)
    }
    if (cores > 1){
      output = .gas_velocity_spaxels_mc(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, cores)
    }

    cube = array(data = output[[1]], dim = c(observation$sbin, observation$sbin, observation$vbin))
    raw_images = list(
      mass_image = array(data = summed_images$mass, dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(data = output[[2]], dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(data = output[[3]], dim = c(observation$sbin, observation$sbin)),
      SFR_image = array(data = summed_images$sfr, dim = c(observation$sbin, observation$sbin)),
      metallicity_image = array(data = output[[4]], dim = c(observation$sbin, observation$sbin)),
      OH_image  = array(data = output[[5]], dim = c(observation$sbin, observation$sbin)),
      particle_image = array(data = summed_images$N, dim = c(observation$sbin, observation$sbin)),
      voronoi_bins = array(data = output[[6]], dim = c(observation$sbin, observation$sbin))
      )
    observed_images = list(
      mass_image = array(data = summed_images$mass, dim = c(observation$sbin, observation$sbin)),
      velocity_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      dispersion_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      h3_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      h4_image = array(0.0, dim = c(observation$sbin, observation$sbin)),
      residuals = array(0.0, dim = c(observation$sbin, observation$sbin)),
      SFR_image = array(data = summed_images$sfr, dim = c(observation$sbin, observation$sbin))
    )

    output = list("velocity_cube"   = cube,
                  "observation"     = observation,
                  "raw_images"      = raw_images,
                  "observed_images" = observed_images,
                  "variance_cube"   = NULL)

    if (observation$psf_fwhm > 0){
      if (verbose){cat("Convolving cube with PSF... \n")    }
      output = blur_datacube(output) # apply psf convolution to each cube plane
    }

    if (!is.na(observation$signal_to_noise)){ # should we add noise?
      if (verbose){cat("Adding noise... \n")    }
      noise_cube = array(data = 0, dim = dim(output$velocity_cube))
      output$variance_cube = noise_cube # initialising empty arrays

      noise_cube = .add_noise(output$velocity_cube,
                              sqrt(median(raw_images$mass_image[raw_images$particle_image > 0], na.rm=T))/
                                (observation$signal_to_noise*sqrt(raw_images$mass_image)))
      noise_image = output$observed_images$flux_image*(rowSums(noise_cube, dims=2)/rowSums(output$velocity_cube, dims=2))

      output$velocity_cube = output$velocity_cube + noise_cube
      output$observed_images$flux_image = output$observed_images$flux_image + noise_image

      if (any(output$velocity_cube<0)){
        output$velocity_cube[which(output$velocity_cube<0)] = 0
      } # removing negative values introduced at the edges of the LOSVD
      output$variance_cube = 1/(noise_cube)^2
      output$variance_cube[is.infinite(output$variance_cube)] = 0
      remove(noise_cube)
    }

    dims = dim(output$raw_images$velocity_image)

    if (moments == 4){
      for (c in 1:dims[1]){
        for (d in 1:dims[2]){
          vel_ini = .meanwt(observation$vbin_seq, output$velocity_cube[c,d,])
          sd_ini  = sqrt(.varwt(observation$vbin_seq, output$velocity_cube[c,d,], vel_ini))

          kin  = tryCatch({stats::optim(par   = c(vel_ini,sd_ini,0,0),
                                        fn    = .losvd_fit_h3h4,
                                        x     = observation$vbin_seq,
                                        losvd = (output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T))),
                                        method="L-BFGS-B", lower = c(NA, 0, NA, NA),
                                        control=list(pgtol=1e-9))$par},
                          error = function(e){c(0,0,0,0)})

          output$observed_images$velocity_image[c,d]   = kin[1]
          output$observed_images$dispersion_image[c,d] = kin[2]
          output$observed_images$h3_image[c,d]       = kin[3]
          output$observed_images$h4_image[c,d]       = kin[4]
          output$observed_images$residuals[c,d]      = mean(abs(.losvd_out_h3h4(x=observation$vbin_seq, vel=kin[1], sig=kin[2], h3=kin[3], h4=kin[4]) -
                                                                  output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T)) ), na.rm=T)

        }
      }
    }

    if (moments == 2){
      for (c in 1:dims[1]){
        for (d in 1:dims[2]){
          vel_ini = .meanwt(observation$vbin_seq, output$velocity_cube[c,d,])
          sd_ini  = sqrt(.varwt(observation$vbin_seq, output$velocity_cube[c,d,], vel_ini))

          kin  = tryCatch({stats::optim(par   = c(vel_ini,sd_ini),
                                        fn    = .losvd_fit_vsig,
                                        x     = observation$vbin_seq,
                                        losvd = (output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T))),
                                        method="L-BFGS-B", lower = c(NA, 0),
                                        control=list(pgtol=1e-9))$par},
                          error = function(e){c(0,0)})

          output$observed_images$velocity_image[c,d]   = kin[1]
          output$observed_images$dispersion_image[c,d] = kin[2]
          output$observed_images$h3_image[c,d]       = 0
          output$observed_images$h4_image[c,d]       = 0
          output$observed_images$residuals[c,d]      = mean(abs(.losvd_out_vsig(x=observation$vbin_seq, vel=kin[1], sig=kin[2]) -
                                                                  output$velocity_cube[c,d,]/(max(output$velocity_cube[c,d,], na.rm=T)) ), na.rm=T)

        }
      }
    }


    if (verbose){cat("Done! \n")}

  }

  if (!voronoi_bin){ # if vorbin has not been requested, don't supply it in the output
    output$raw_images = output$raw_images[-which(names(output$raw_images) == "voronoi_bins")]
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
                       input_simspin_file_path = rev(stringr::str_split(simspin_file, "/")[[1]])[1])
  }

  return(output)
}

