# Author: Kate Harborne
# Date: 26/10/2020
# Title: write_simspin_FITS - a function for making a spectral FITS file output
#
#'A function for making a mock spectral cube FITS file output
#'
#'The purpose of this function is to write the spectral cube generated using
#' the \code{build_datacube()} function to a FITS file. This code will also
#' write the images produced by \code{build_datacube()} to a further series of
#' FITS files automatically. Images are either added to subsequent HDU within
#' the same FITS file, or saved seperately to a series of files. The number of
#' FITS HDU or files written is dependent on the input simspin_datacube. Files
#' and HDU will be named suitably to reflect their content automatically.
#'
#'@param output_file The path and file name to the location where the FITS file
#' should be written.
#'@param simspin_datacube The list output from \code{build_datacube()}.
#'@param object_name A string that described the name of the observed object.
#'@param telescope_name A string that describes the name of the telescope used.
#'@param instrument_name A string that describes the used instrument on that
#' telescope.
#'@param observer_name A string that describes the name of the observer.
#'@param input_simspin_file_path A string describing the path to the original
#' SimSpin file (i.e. the location of the file output from make_simspin_file).
#'@param split_save Boolean describing whether to split the output from
#' \code{build_datacube()} into multiple files. If TRUE, several FITS files
#' will be saved with file names that reflect their content (i.e.
#' "_spectral_cube.FITS", "_velocity_image.FITS", "_dispersion_images.FITS",
#' etc.). Default option is FALSE.
#'@param mask (Optional) A binary array describing the masked regions of the
#' cube/images.
#'@param galaxy_centre (Optional) A numeric array (x,y,z) describing the centre
#' of potential of the observed galaxy within it's simulation.
#'@return Returns an .FITS file that contains a the generated data and
#' relevant header describing the mock observation.
#'@examples
#'\dontrun{
#'ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
#'temp_loc = tempdir()
#'make_simspin_file(ss_eagle, output = paste(temp_loc, "spectra.Rdata", sep=""))
#'cube = build_datacube(simspin_file = paste(temp_loc, "spectra.Rdata", sep=""),
#'                      telescope = telescope(type="SAMI", lsf_fwhm = 4, wave_res = 1.6),
#'                      observing_strategy = observing_strategy())
#'write_simspin_FITS(output_file = paste(temp_loc, "cube.FITS", sep=""),
#'                   simspin_datacube = cube, object_name = "SimSpin EAGLE example",
#'                   telescope_name = "AAO", instrument_name = "SAMI",
#'                   observer_name = "K.E.Harborne",
#'                   input_simspin_file_path = paste(temp_loc, "spectra.Rdata", sep=""))
#'unlink(paste(temp_loc, "spectra.fst", sep=""))
#'unlink(paste(temp_loc, "cube.FITS", sep=""))
#'}
#'

write_simspin_FITS = function(output_file, simspin_datacube, object_name,
                              telescope_name, instrument_name, observer_name,
                              input_simspin_file_path, split_save=F, mask=NA,
                              galaxy_centre = c(0,0,0), ...){

  # Running checks on inputs ---------------------------------------------------
  if (missing(input_simspin_file_path)){
    if (exists(input_simspin_file)){
      warning("Warning: No provided `input_simspin_file_path`. \n  You have provided `input_simspin_file` instead, which will unavailable in future SimSpin versions (>v2). \n Please adjust and try again.")
      input_simspin_file_path = input_simspin_file
      remove(input_simspin_file)
    }
  }

  if (!is.character(input_simspin_file_path)){
    stop("Error: `input_simspin_file_path` is expected as a string, describing the path to the input SimSpin file. \n Please adjust and try again.")
  }

  if (!is.list(simspin_datacube)){
    stop("Error: `simspin_datacube` is expected as an R object list, as output by `build_datacube`. \n Please adjust and try again.")
  }

  if (!all(is.character(object_name), is.character(telescope_name),
           is.character(instrument_name), is.character(observer_name),
           is.character(output_file))){
    stop("Error: One of `output_file`, `object_name`, `telescope_name`, `instrument_name`, `observer_name` is expected as a string, but was not provided in the correct format. \n Please adjust and try again.")
  }

  if (!is.logical(split_save)){
    stop("Error: `split_save` should be provided as a logical TRUE or FALSE.\n Please adjust and try again.")
  }

  if (!all(is.na(mask)) & !is.double(mask)){
    stop("Error: `mask` should be provided either as a logical NA, or as a 2-dimensional numeric array.\n Please adjust and try again.")
  }


  # Sorting names for each file ------------------------------------------------
  output_name = rev(stringr::str_split(output_file, "/")[[1]])[1]
  output_dir  = stringr::str_remove(output_file, output_name)
  output_file_root = stringr::str_remove(stringr::str_remove(output_name, ".fits"), ".FITS")

  observation = simspin_datacube$observation

  voronoi = FALSE # flag to determine behaviour of writing FITS files of voronoi binned data
  if ("voronoi_bins" %in% names(simspin_datacube$raw_images)){voronoi=TRUE}

  mass_flag = observation$mass_flag # flag to determine whether units of output cubes

  simspin_cube = simspin_datacube[[1]] # getting either the "velocity cube" or the "spectral cube"

  galaxy_centre_norm = galaxy_centre / sqrt((galaxy_centre[1]^2) + (galaxy_centre[2]^2) + (galaxy_centre[3]^2))

  if (any(is.na(galaxy_centre_norm))){
    galaxy_centre_norm = c(0,0,0)
  }

  gal_DEC = asin(galaxy_centre_norm[3])
  gal_RA  = asin(galaxy_centre_norm[2]/cos(gal_DEC))

  # FITS files always have the first HDU containing the header information.

  # Basic header information is initiallised in this function to begin.
  # This is then modified in the following "if"s depending on the telescope
  # mode employed.

  # Making header for all files ----

  header_keyvalues = list("SIMPLE"=TRUE, "BITPIX"=8, "NAXIS"=0, "EXTEND"=TRUE,
                          "DATE"=Sys.time(), "ORIGIN"=observation$origin,
                          "TELESCOP"=telescope_name,
                          "INSTRUME"=instrument_name, "RA"=(observation$pointing_deg[1]+(gal_RA*(180/pi))),
                          "DEC"=(observation$pointing_deg[2]+(gal_DEC*(180/pi))), "EQINOX"=2000,
                          "RADECSYS"="FK5", "EXPTIME"=1320, "MJD-OBS"=58906.11,
                          "DATE-OBS"=observation$date, "UTC"=9654, "LST"=30295.18,
                          "PI-COI"="UNKNOWN", "OBSERVER"=observer_name, "REDSHIFT"=observation$z,
                          "PIPEFILE"=output_name,
                          "BUNIT"=if(mass_flag){"Msol"}else{"erg/s/cm**2"},
                          "ARCFILE"=input_simspin_file_path,
                          "DATAMD5"="4aece79473a5c88f6533382655e948bf",
                          "OBJECT"=object_name)

  header_keycomments = list("SIMPLE"="file does conform to FITS standard",
                            "BITPIX"="number of bits per data pixel",
                            "NAXIS"="number of data axes",
                            "EXTEND"="FITS dataset may contain extensions",
                            "DATE"="file creation date (YYYY-MM-DD hh:mm:ss timezone)",
                            "ORIGIN"="Mock observation",
                            "TELESCOP"="Telescope name",
                            "INSTRUME"="Instrument used",
                            "RA"="[deg] 11:41:21.1 RA (J2000) pointing",
                            "DEC"="[deg] -01:34:59.0 DEC (J2000) pointing",
                            "EQINOX"="Standard FK5",
                            "RADECSYS"="Coordinate system",
                            "EXPTIME"="Integration time",
                            "MJD-OBS"="Obs start",
                            "DATE-OBS"="Observing date",
                            "UTC"="[s] 02:40:54.000 UTC",
                            "LST"="[s] 08:24:55.178 LST",
                            "PI-COI"="PI-COI name.", "OBSERVER"="Name of observer.",
                            "REDSHIFT"="Observed redshift.",
                            "PIPEFILE"="Filename of data product",
                            "BUNIT"="units of values in image or cube",
                            "ARCFILE"="Archive File Name",
                            "DATAMD5"="MD5 checksum",
                            "OBJECT"="Original target.")

  # Data header to be modified and added for each HDU ----

  data_keyvalues = list("XTENSION" = "IMAGE", "BITPIX"=-32, "NAXIS"=3,
                        "NAXIS1"=dim(simspin_cube)[1], "NAXIS2"=dim(simspin_cube)[2],
                        "NAXIS3"=dim(simspin_cube)[3], "PCOUNT"=0, "GCOUNT"=1,
                        "EXTNAME"="DATA", "HDUCLASS"="ESO",
                        "HDUDOC"="DICD", "HDUVERS"="DCID version 6",
                        "HDUCLAS1"="IMAGE", "HDUCLAS2"="DATA", "ERRDATA"="STAT",
                        "OBJECT"=object_name, "BUNIT"=if(mass_flag){"Msol"}else{"erg/s/cm**2"},
                        "CRPIX1"=1, "CRPIX2"=1,
                        "CDELT1"=observation$sbin_size/observation$ang_size/3600,
                        "CDELT2"=observation$sbin_size/observation$ang_size/3600,
                        "CUNIT1"="deg", "CUNIT2"="deg",
                        "CTYPE1"="RA---TAN", "CTYPE2"="DEC--TAN",
                        "CSYER1"=character(1), "CSYER2"=character(1),
                        "CRVAL1"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                        "CRVAL2"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                        "CTYPE3"=character(1), "CUNIT3"=character(1),
                        "CDELT3"=numeric(1), "CRPIX3"=1,
                        "CRVAL3"=numeric(1))

  data_keycomments = list("XTENSION"="IMAGE extension",
                          "BITPIX"="number of bits per data pixel",
                          "NAXIS"="number of data axes",
                          "NAXIS1"="axis length", "NAXIS2"="axis length",
                          "NAXIS3"="axis length",
                          "PCOUNT"="required keyword; must = 0",
                          "GCOUNT"="required keyword; must = 1",
                          "EXTNAME"="This extension contains data values",
                          "HDUCLASS"="class name (ESO format)",
                          "HDUDOC"="document with class description",
                          "HDUVERS"="version number (according to spec v2.5.1)",
                          "HDUCLAS1"="Image data format",
                          "HDUCLAS2"="this extension contains the data itself",
                          "ERRDATA"="pointer to the variance extension",
                          "OBJECT"="simulation and galaxy ID number",
                          "BUNIT"="units of values in image or cube",
                          "CRPIX1"="Pixel coordinate of reference point",
                          "CRPIX2"="Pixel coordinate of reference point",
                          "CDELT1"="Coordinate transformation matrix element",
                          "CDELT2"="Coordinate transformation matrix element",
                          "CUNIT1"="Units of coordinate increment and value",
                          "CUNIT2"="Units of coordinate increment and value",
                          "CTYPE1"="Right ascension, gnomonic projection",
                          "CTYPE2"="Declination, gnomonic projection",
                          "CSYER1"="[deg] Systematic error in coordinate",
                          "CSYER2"="[deg] Systematic error in coordinate",
                          "CRVAL1"="Pixel value at reference point",
                          "CRVAL2"="Pixel value at reference point",
                          "CTYPE3"=character(1),
                          "CUNIT3"="Units of coordinate increment and value",
                          "CDELT3"="Coordinate transformation matrix element",
                          "CRPIX3"="Pixel coordinate of reference point",
                          "CRVAL3"="Pixel value at reference point")

  # Building a table for the observation summary for all files -------------------
  # Need to begin by coercing the values that are arrays into comma seperated
  # character vectors of length 1.

  obs_summary = stats::setNames(data.table::data.table(matrix("",
                                                              ncol = 3,
                                                              nrow = (length(observation)-4))),
                                c("Name", "Value", "Units"))

  obs_names = names(observation)
  obs_names = obs_names[-c(which(obs_names %in% c("aperture_region","psf_kernel", "pixel_region")))]

  obs_units = data.table::data.table("ang_size" = "num: scale at given distance in kpc/arcsec",
                                     "aperture_shape" = "str: shape of aperture",
                                     "aperture_size" = "num: field of view diameter width in kpc",
                                     "date" = "str: date and time of mock observation",
                                     "fov" = "num: field of view diameter in arcsec",
                                     "filter" = "str: filter name",
                                     "inc_deg" = "num: projected inclination of object in degrees about the horizontal axis",
                                     "inc_rad" = "num: projected inclination of object in radians about the horizontal axis",
                                     "twist_deg" = "num: projected inclination of object in degrees about the vertical axis",
                                     "twist_rad" = "num: projected inclination of object in radians about the vertical axis",
                                     "lsf_fwhm" = "num: line-spread function of telescope given as full-width half-maximum in Angstrom",
                                     "LSF_conv" = "bool: has line spread function convolution been applied?",
                                     "lsf_sigma" = "num: line-spread function of telescope given as a sigma width in Angstrom",
                                     "lum_dist" = "num: distance to object in Mpc",
                                     "mass_flag" = "bool: kinematics are mass-weighted if true",
                                     "method" = "str: name of observing method employed",
                                     "moments" = "num: Whether a 2-moment (Gaussian) or 4-moment (Gauss-Hermite) function has been fitted to the LOSVD",
                                     "origin" = "str: version of SimSpin used for observing",
                                     "particle_limit" = "int: minimum number of particles per pixel. If 0, model has not been voronoi binned",
                                     "pointing_kpc" = "num: x-y position of field of view centre relative to object centre in units of kpc",
                                     "pointing_deg" = "num: x-y position of field of view centre relative to object centre in units of degrees",
                                     "psf_fwhm" = "num: the full-width half-maximum of the point spread function kernel in arcsec",
                                     "sbin" =  "num: the number of spatial pixels across the diameter of the field of view",
                                     "sbin_seq" = "num: the min and max spatial bin centres in kpc",
                                     "sbin_size" = "num: the size of each pixel in kpc",
                                     "spatial_res" = "num: the size of each pixel in arcsec",
                                     "signal_to_noise" = "num: the signal-to-noise ratio for observed spectrum",
                                     "wave_bin" = "num: the number of wavelength bins for a given telescope",
                                     "wave_centre" = "num: the central wavelength for a given telescope in Angstrom",
                                     "wave_res" = "num: the width of each wavelength bin in Angstrom",
                                     "wave_seq" = "num: the min and max wavelength bin centres in Angstrom",
                                     "wave_edges" = "num: the wavelength bin edges in Angstrom",
                                     "vbin" = "num: the number of velocity bins for a given telescope resolution",
                                     "vbin_seq" = "num: the min and max velocity bin centres in km/s",
                                     "vbin_edges" = "num: the min and max velocity bin edges in km/s",
                                     "vbin_size" = "num: the size of each velocity bin in km/s",
                                     "vbin_error" = "num: the velocity uncertainty given the telescope LSF in km/s",
                                     "z" = "num: the redshift distance of the object observed")

  if ("vbin_seq" %in% obs_names){
    observation[["vbin_seq"]] = range(observation[["vbin_seq"]])
    observation[["vbin_edges"]] = range(observation[["vbin_edges"]])
  }
  observation[["wave_seq"]] = range(observation[["wave_seq"]])
  observation[["wave_edges"]] = range(observation[["wave_edges"]])
  observation[["sbin_seq"]] = range(observation[["sbin_seq"]])

  if ("filter_name" %in% names(observation)){
    observation[["filter"]] = observation[["filter_name"]]
  } else {
    if (min(range(observation$filter$wave)) == 2980){
      observation[["filter_name"]] = "u_SDSS"
    }
    if (min(range(observation$filter$wave)) == 3630){
      observation[["filter_name"]] = "g_SDSS"
    }
    if (min(range(observation$filter$wave)) == 5380){
      observation[["filter_name"]] = "r_SDSS"
    }
    if (min(range(observation$filter$wave)) == 6430){
      observation[["filter_name"]] = "i_SDSS"
    }
    if (min(range(observation$filter$wave)) == 7730){
      observation[["filter_name"]] = "z_SDSS"
    }
  }
  observation[["filter"]] = observation[["filter_name"]]
  obs_names = obs_names[-c(which(obs_names %in% c("filter_name")))]

  obs_summary[, "Name":= obs_names]

  for (val in 1:(length(observation)-4)){

    if ( length(observation[[obs_names[val]]]) > 1){
      rounded_val = signif(observation[[obs_names[val]]], digits = 6)
      add_value = paste(as.character(rounded_val), collapse = ",")
    } else {
      if (is.na(observation[[obs_names[val]]])){
        add_value = "None"
      } else {
        rounded_val = if(is.numeric(observation[[obs_names[val]]])){
          signif(observation[[obs_names[val]]], digits = 6)
        } else {
            observation[[obs_names[val]]]
          }
        add_value = as.character(rounded_val)
      }
    }

    data.table::set(obs_summary, i =val, j="Value", value=add_value)
    data.table::set(obs_summary, i =val, j="Units", value=obs_units[[obs_names[val]]])
  }


  # Writing data to HDUs based on the observation method employed ----

  # SPECTRAL mode method =======================================================
  if (observation$method == "spectral"){

    if (split_save){
      cube_file_name = paste0(output_dir, "/", output_file_root, "_spectral_cube.FITS") # adding "cube" to end of file name
    } else {
      cube_file_name = output_file
    }

    Rfits::Rfits_write_header(filename = cube_file_name, keyvalues = header_keyvalues,
                              keycomments = header_keycomments, ext=1, create_file = T,
                              overwrite_file = TRUE)

    data_keyvalues$CTYPE3 = "WAVE"
    data_keyvalues$CUNIT3 = "Angstrom"
    data_keyvalues$CDELT3 = observation$wave_res
    data_keyvalues$CRVAL3 = observation$wave_seq[1]

    Rfits::Rfits_write_cube(data = simspin_cube, filename = cube_file_name, ext=2,
                            keyvalues = data_keyvalues, keycomments = data_keycomments,
                            create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

    # Adding observation summary to FITS file ----

    if (split_save){
      obs_summary_file_name = paste0(output_dir, "/", output_file_root, "_observation_summary.FITS")
      Rfits::Rfits_write_header(filename = obs_summary_file_name,
                                keyvalues = header_keyvalues,
                                keycomments = header_keycomments, ext=1, create_file = T,
                                overwrite_file = TRUE)

      Rfits::Rfits_write_table(obs_summary, filename = obs_summary_file_name,
                               ext = 2, extname = "OB_TABLE",
                               create_ext = TRUE, create_file = FALSE,
                               overwrite_file = FALSE)
    } else {
      Rfits::Rfits_write_table(obs_summary, filename = cube_file_name,
                               ext = 3, extname = "OB_TABLE",
                               create_ext = TRUE, create_file = FALSE,
                               overwrite_file = FALSE)

    }

    # Adding images to FITS file ----

    image_keyvalues = list("XTENSION"="IMAGE", "BITPIX"=-64, "NAXIS"=2,
                           "NAXIS1"=dim(simspin_cube)[1], "NAXIS2"=dim(simspin_cube)[2],
                           "PCOUNT"=0, "GCOUNT"=1, "BUNIT"=character(1),
                           "CRPIX1"=1,
                           "CRVAL1"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                           "CDELT1"=observation$sbin_size/observation$ang_size/3600,
                           "CTYPE1"="RA---TAN", "CUNIT1"="deg",
                           "CRPIX2"=1,
                           "CRVAL2"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                           "CDELT2"=observation$sbin_size/observation$ang_size/3600,
                           "CTYPE2"="DEC--TAN", "CUNIT2"="deg", "EXTNAME"=character(1))

    image_keycomments = list("XTENSION"="IMAGE extension",
                             "BITPIX"="number of bits per data pixel",
                             "NAXIS"="number of data axes",
                             "NAXIS1"="axis length", "NAXIS2"="axis length",
                             "PCOUNT"="required keyword; must = 0",
                             "GCOUNT"="required keyword; must = 1",
                             "BUNIT"="units of image values",
                             "CRPIX1"="Pixel coordinate of reference point",
                             "CRVAL1"="Pixel value at reference point",
                             "CDELT1"="Coordinate transformation matrix element",
                             "CTYPE1"="Right ascension, gnomonic projection",
                             "CUNIT1"="Units of coordinate increment and value",
                             "CRPIX2"="Pixel coordinate of reference point",
                             "CRVAL2"="Pixel value at reference point",
                             "CDELT2"="Coordinate transformation matrix element",
                             "CTYPE2"="Declination, gnomonic projection",
                             "CUNIT2"="Units of coordinate increment and value",
                             "EXTNAME"="Image extension name")

    extnames =
      if (voronoi){
        c("RAW_FLUX", "RAW_VEL", "RAW_DISP", "RAW_AGEM", "RAW_AGEL", "RAW_Z", "NPART", "RAW_MASS", "VORONOI")
      } else {
        c("RAW_FLUX", "RAW_VEL", "RAW_DISP", "RAW_AGEM", "RAW_AGEL", "RAW_Z", "NPART", "RAW_MASS")
      }

    bunits =
      if (voronoi){
        c("erg/s/cm**2", "km/s", "km/s", "mass-weighted Gyr", "light-weighted Gyr", "Z_solar", "Particle number", "Msol", "Bin ID")
        } else {
        c("erg/s/cm**2", "km/s", "km/s", "mass-weighted Gyr", "light-weighted Gyr", "Z_solar", "Particle number", "Msol")
        }

    image_names =
      if (voronoi){
        c("flux_image", "velocity_image", "dispersion_image", "ageM_image", "ageL_image",
          "metallicity_image", "particle_image", "mass_image", "voronoi_bins")
      } else {
        c("flux_image", "velocity_image", "dispersion_image", "ageM_image", "ageL_image",
          "metallicity_image", "particle_image", "mass_image")
      }

    extnum =
      if (voronoi){
        c(4,5,6,7,8,9,10,11,12)
      } else {
        c(4,5,6,7,8,9,10,11)
      }

    output_image_file_names = paste0(output_dir, "/", output_file_root, "_raw_", image_names, ".FITS")

    if (split_save){ # if writing each image to a seperate file
      for (i in 1:length(extnum)){ # for each image in the build_datacube output,

        # 1. Write the header again to the new file to HDU 1
        Rfits::Rfits_write_header(filename = output_image_file_names[i], keyvalues = header_keyvalues,
                                  keycomments = header_keycomments, ext=1, create_file = T,
                                  overwrite_file = TRUE)

        image_keyvalues$BUNIT = bunits[i]
        image_keyvalues$EXTNAME = extnames[i]

        # 2. Write the image to this new file HDU 2
        Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                 filename = output_image_file_names[i], ext=2,
                                 keyvalues = image_keyvalues, keycomments = image_keycomments,
                                 create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
      }
    } else { # if writing all images to the same file
      for (i in 1:length(extnum)){ # for each image in the build_datacube output,

        image_keyvalues$BUNIT = bunits[i]
        image_keyvalues$EXTNAME = extnames[i]

        # Write each subsequent image to the next HDU
        Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                 filename = cube_file_name, ext=extnum[i],
                                 keyvalues = image_keyvalues, keycomments = image_keycomments,
                                 create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
      }

    }

    if (!is.na(observation$signal_to_noise)){
      # Adding variance cube to FITS file ----
      data_keyvalues$EXTNAME = "STAT"
      data_keyvalues$BUNIT = "1**40 (erg/s/cm**2)**-2"

      if (split_save){
        stat_summary_file_name = paste0(output_dir, "/", output_file_root, "_inv_variance_cube.FITS")

        Rfits::Rfits_write_header(filename = stat_summary_file_name,
                                  keyvalues = header_keyvalues,
                                  keycomments = header_keycomments, ext=1, create_file = T,
                                  overwrite_file = TRUE)
        Rfits::Rfits_write_cube(data = simspin_datacube$variance_cube/1e40,
                                filename = stat_summary_file_name, ext=2,
                                keyvalues = data_keyvalues, keycomments = data_keycomments,
                                create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

      } else {
        Rfits::Rfits_write_cube(data = simspin_datacube$variance_cube/1e40,
                                filename = cube_file_name, ext=(max(extnum)+1),
                                keyvalues = data_keyvalues,
                                keycomments = data_keycomments,
                                create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
      }
    }

  }

  # VELOCITY mode method =======================================================
  if (observation$method == "velocity"){

    if (split_save){
      cube_file_name = paste0(output_dir, "/", output_file_root, "_velocity_cube.FITS")
    } else {
      cube_file_name = output_file
    }

    Rfits::Rfits_write_header(filename = cube_file_name, keyvalues = header_keyvalues,
                              keycomments = header_keycomments, ext=1, create_file = T,
                              overwrite_file = TRUE)

    data_keyvalues$CTYPE3 = "VELOCITY"
    data_keyvalues$CUNIT3 = "km/s"
    data_keyvalues$CDELT3  = observation$vbin_size
    data_keyvalues$CRVAL3 = observation$vbin_seq[1]

    Rfits::Rfits_write_cube(data = simspin_cube, filename = cube_file_name, ext=2,
                            keyvalues = data_keyvalues, keycomments = data_keycomments,
                            create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

    # Adding observation summary to FITS file ----
    if (split_save){
      obs_summary_file_name = paste0(output_dir, "/", output_file_root, "_observation_summary.FITS")
      Rfits::Rfits_write_header(filename = obs_summary_file_name,
                                keyvalues = header_keyvalues,
                                keycomments = header_keycomments, ext=1, create_file = T,
                                overwrite_file = TRUE)

      Rfits::Rfits_write_table(obs_summary, filename = obs_summary_file_name,
                               ext = 2, extname = "OB_TABLE",
                               create_ext = TRUE, create_file = FALSE,
                               overwrite_file = FALSE)
    } else {
      Rfits::Rfits_write_table(obs_summary, filename = cube_file_name,
                               ext = 3, extname = "OB_TABLE",
                               create_ext = TRUE, create_file = FALSE,
                               overwrite_file = FALSE)

    }

    # Adding images to FITS file ----

    image_keyvalues = list("XTENSION"="IMAGE", "BITPIX"=-64, "NAXIS"=2,
                           "NAXIS1"=dim(simspin_cube)[1], "NAXIS2"=dim(simspin_cube)[2],
                           "PCOUNT"=0, "GCOUNT"=1, "BUNIT"=character(1),
                           "CRPIX1"=1,
                           "CRVAL1"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                           "CDELT1"=observation$sbin_size/observation$ang_size/3600,
                           "CTYPE1"="RA---TAN", "CUNIT1"="deg",
                           "CRPIX2"=1,
                           "CRVAL2"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                           "CDELT2"=observation$sbin_size/observation$ang_size/3600,
                           "CTYPE2"="DEC--TAN", "CUNIT2"="deg", "EXTNAME"=character(1))

    image_keycomments = list("XTENSION"="IMAGE extension",
                             "BITPIX"="number of bits per data pixel",
                             "NAXIS"="number of data axes",
                             "NAXIS1"="axis length", "NAXIS2"="axis length",
                             "PCOUNT"="required keyword; must = 0",
                             "GCOUNT"="required keyword; must = 1",
                             "BUNIT"="units of image values",
                             "CRPIX1"="Pixel coordinate of reference point",
                             "CRVAL1"="Pixel value at reference point",
                             "CDELT1"="Coordinate transformation matrix element",
                             "CTYPE1"="Right ascension, gnomonic projection",
                             "CUNIT1"="Units of coordinate increment and value",
                             "CRPIX2"="Pixel coordinate of reference point",
                             "CRVAL2"="Pixel value at reference point",
                             "CDELT2"="Coordinate transformation matrix element",
                             "CTYPE2"="Declination, gnomonic projection",
                             "CUNIT2"="Units of coordinate increment and value",
                             "EXTNAME"="Image extension name")

    extnames    = if (voronoi){
                      c("OBS_FLUX", "OBS_VEL", "OBS_DISP", "OBS_H3", "OBS_H4", "RESIDUAL", "OBS_MASS", "RAW_FLUX", "RAW_MASS", "RAW_VEL", "RAW_DISP", "RAW_AGEM", "RAW_AGEL", "RAW_Z", "NPART", "VORONOI")
                    } else {
                      c("OBS_FLUX", "OBS_VEL", "OBS_DISP", "OBS_H3", "OBS_H4", "RESIDUAL", "OBS_MASS", "RAW_FLUX", "RAW_MASS", "RAW_VEL", "RAW_DISP", "RAW_AGEM", "RAW_AGEL", "RAW_Z", "NPART")
                    }

    bunits      = if (voronoi){
                      c("erg/s/cm**2", "km/s", "km/s", "unitless", "unitless", "percentage", "Msol", "erg/s/cm**2", "Msol", "km/s", "km/s", "mass-weighted Gyr", "light-weighted Gyr", "Z_solar", "Particle number", "Bin ID")
                    } else {
                      c("erg/s/cm**2", "km/s", "km/s", "unitless", "unitless", "percentage", "Msol", "erg/s/cm**2", "Msol", "km/s", "km/s", "mass-weighted Gyr", "light-weighted Gyr", "Particle number")
                    }

    image_names = if (voronoi){
                      c("flux_image", "velocity_image", "dispersion_image", "h3_image", "h4_image",
                        "residuals", "mass_image",
                        "flux_image", "mass_image", "velocity_image", "dispersion_image", "ageM_image", "ageL_image",
                        "metallicity_image", "particle_image", "voronoi_bins")
                    } else {
                      c("flux_image", "velocity_image", "dispersion_image", "h3_image", "h4_image",
                        "residuals", "mass_image",
                        "flux_image", "mass_image", "velocity_image", "dispersion_image", "ageM_image", "ageL_image",
                        "metallicity_image", "particle_image")
                    }

    rawobs = if (voronoi){
                 c("obs", "obs", "obs", "obs", "obs", "obs",  "obs",
                   "raw", "raw", "raw", "raw", "raw", "raw", "raw", "raw", "raw")
               } else {
                 c("obs", "obs", "obs", "obs", "obs", "obs",  "obs",
                   "raw", "raw", "raw", "raw", "raw", "raw", "raw", "raw")
               }

    output_image_file_names = paste0(output_dir, "/", output_file_root, "_", rawobs, "_", image_names, ".FITS")

    extnum = if (voronoi){
                  c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
                } else {
                  c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
                }

    if (split_save){ # if writing each image to a seperate file
      for (i in 1:length(extnum)){ # for each image in the build_datacube output,

        # 1. Write the header again to the new file to HDU 1
        Rfits::Rfits_write_header(filename = output_image_file_names[i], keyvalues = header_keyvalues,
                                  keycomments = header_keycomments, ext=1, create_file = T,
                                  overwrite_file = TRUE)

        image_keyvalues$BUNIT = bunits[i]
        image_keyvalues$EXTNAME = extnames[i]

        if (i < 8){ # 2. Write the image to this new file HDU 2
          Rfits::Rfits_write_image(data = simspin_datacube$observed_images[[which(names(simspin_datacube$observed_images) == image_names[i])]],
                                   filename = output_image_file_names[i], ext=2,
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
        } else {
          Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                   filename = output_image_file_names[i], ext=2,
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

        }

      }
    } else { # if writing all images to the same file
      for (i in 1:length(extnum)){ # for each image in the build_datacube output,

        image_keyvalues$BUNIT = bunits[i]
        image_keyvalues$EXTNAME = extnames[i]

        if (i < 8){ # Write each subsequent image to the next HDU
          Rfits::Rfits_write_image(data = simspin_datacube$observed_images[[which(names(simspin_datacube$observed_images) == image_names[i])]],
                                   filename = cube_file_name, ext=extnum[i],
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
        } else {
          Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                   filename = cube_file_name, ext=extnum[i],
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

        }

      }

    }

    if (!is.na(observation$signal_to_noise)){
      # Adding variance cube to FITS file ----
      data_keyvalues$EXTNAME = "STAT"
      data_keyvalues$BUNIT   = "1**40 (erg/s/cm**2)**-2"
      if (split_save){
        stat_summary_file_name = paste0(output_dir, "/", output_file_root, "_inv_variance_cube.FITS")

        Rfits::Rfits_write_header(filename = stat_summary_file_name,
                                  keyvalues = header_keyvalues,
                                  keycomments = header_keycomments, ext=1, create_file = T,
                                  overwrite_file = TRUE)
        Rfits::Rfits_write_cube(data = simspin_datacube$variance_cube/1e40,
                                filename = stat_summary_file_name, ext=2,
                                keyvalues = data_keyvalues, keycomments = data_keycomments,
                                create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

      } else {
        Rfits::Rfits_write_cube(data = simspin_datacube$variance_cube/1e40,
                                filename = cube_file_name, ext=(max(extnum)+1),
                                keyvalues = data_keyvalues,
                                keycomments = data_keycomments,
                                create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
      }
    }
  }

  # GAS mode method ============================================================
  if (observation$method == "gas" | observation$method == "sf gas"){

    if (split_save){
      cube_file_name = paste0(output_dir, "/", output_file_root, "_gas_velocity_cube.FITS")
    } else {
      cube_file_name = output_file
    }

    Rfits::Rfits_write_header(filename = cube_file_name, keyvalues = header_keyvalues,
                              keycomments = header_keycomments, ext=1, create_file = T,
                              overwrite_file = TRUE)

    data_keyvalues$CTYPE3 = "GAS_VELO"
    data_keyvalues$CUNIT3 = "km/s"
    data_keyvalues$CDELT3 = observation$vbin_size
    data_keyvalues$CRVAL3 = observation$vbin_seq[1]

    Rfits::Rfits_write_cube(data = simspin_cube, filename = cube_file_name, ext=2,
                            keyvalues = data_keyvalues, keycomments = data_keycomments,
                            create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

    # Adding observation summary to FITS file ----
    if (split_save){
      obs_summary_file_name = paste0(output_dir, "/", output_file_root, "_observation_summary.FITS")

      Rfits::Rfits_write_header(filename = obs_summary_file_name,
                                keyvalues = header_keyvalues,
                                keycomments = header_keycomments, ext=1, create_file = T,
                                overwrite_file = TRUE)
      Rfits::Rfits_write_table(obs_summary, filename = obs_summary_file_name,
                               ext = 2, extname = "OB_TABLE",
                               create_ext = TRUE, create_file = FALSE,
                               overwrite_file = FALSE)
    } else {
      Rfits::Rfits_write_table(obs_summary, filename = cube_file_name,
                               ext = 3, extname = "OB_TABLE",
                               create_ext = TRUE, create_file = FALSE,
                               overwrite_file = FALSE)

    }

    # Adding images to FITS file ----

    image_keyvalues = list("XTENSION"="IMAGE", "BITPIX"=-64, "NAXIS"=2,
                           "NAXIS1"=dim(simspin_cube)[1], "NAXIS2"=dim(simspin_cube)[2],
                           "PCOUNT"=0, "GCOUNT"=1, "BUNIT"=character(1),
                           "CRPIX1"=1,
                           "CRVAL1"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                           "CDELT1"=-observation$sbin_size/observation$ang_size/3600,
                           "CTYPE1"="RA---TAN", "CUNIT1"="deg",
                           "CRPIX2"=1,
                           "CRVAL2"=(diff(observation$sbin_seq[1:2])/2 + observation$sbin_seq[1])/observation$ang_size/3600,
                           "CDELT2"=observation$sbin_size/observation$ang_size/3600,
                           "CTYPE2"="DEC--TAN", "CUNIT2"="deg", "EXTNAME"=character(1))

    image_keycomments = list("XTENSION"="IMAGE extension",
                             "BITPIX"="number of bits per data pixel",
                             "NAXIS"="number of data axes",
                             "NAXIS1"="axis length", "NAXIS2"="axis length",
                             "PCOUNT"="required keyword; must = 0",
                             "GCOUNT"="required keyword; must = 1",
                             "BUNIT"="units of image values",
                             "CRPIX1"="Pixel coordinate of reference point",
                             "CRVAL1"="Pixel value at reference point",
                             "CDELT1"="Coordinate transformation matrix element",
                             "CTYPE1"="Right ascension, gnomonic projection",
                             "CUNIT1"="Units of coordinate increment and value",
                             "CRPIX2"="Pixel coordinate of reference point",
                             "CRVAL2"="Pixel value at reference point",
                             "CDELT2"="Coordinate transformation matrix element",
                             "CTYPE2"="Declination, gnomonic projection",
                             "CUNIT2"="Units of coordinate increment and value",
                             "EXTNAME"="Image extension name")

    extnames =
      if (voronoi){
        c("OBS_MASS", "OBS_VEL", "OBS_DISP", "OBS_H3", "OBS_H4", "RESIDUAL", "OBS_SFR",
          "RAW_MASS", "RAW_VEL", "RAW_DISP", "RAW_Z", "RAW_OH", "RAW_SFR", "NPART", "VORONOI")
        } else {
          c("OBS_MASS", "OBS_VEL", "OBS_DISP", "OBS_H3", "OBS_H4", "RESIDUAL", "OBS_SFR",
            "RAW_MASS", "RAW_VEL", "RAW_DISP", "RAW_Z", "RAW_OH", "RAW_SFR", "NPART")
        }

    bunits =
      if (voronoi){
        c("Msol", "km/s", "km/s", "unitless", "unitless", "percentage", "Msol/year",
          "Msol", "km/s", "km/s", "log10(Z/Z_solar)", "log10(O/H)+12", "Msol/year",
          "Particle number", "Bin ID")
      } else {
        c("Msol", "km/s", "km/s", "unitless", "unitless", "percentage", "Msol/year",
          "Msol", "km/s", "km/s", "log10(Z/Z_solar)", "log10(O/H)+12", "Msol/year",
          "Particle number")
      }

    image_names =
      if (voronoi){
        c("mass_image", "velocity_image", "dispersion_image", "h3_image", "h4_image",
          "residuals", "SFR_image",
          "mass_image", "velocity_image", "dispersion_image", "metallicity_image",
          "OH_image", "SFR_image", "particle_image", "voronoi_bins")
        } else {
        c("mass_image", "velocity_image", "dispersion_image", "h3_image", "h4_image",
          "residuals", "SFR_image",
          "mass_image", "velocity_image", "dispersion_image", "metallicity_image",
          "OH_image", "SFR_image", "particle_image")
        }

    rawobs =
      if (voronoi){
        c("obs", "obs", "obs", "obs", "obs", "obs", "obs",
          "raw", "raw", "raw", "raw", "raw", "raw", "raw", "raw")
      } else {
        c("obs", "obs", "obs", "obs", "obs", "obs", "obs",
          "raw", "raw", "raw", "raw", "raw", "raw", "raw")
      }

    output_image_file_names = paste0(output_dir, "/", output_file_root, "_", rawobs, "_", image_names, ".FITS")

    extnum =
      if (voronoi){
        c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
      } else {
        c(4,5,6,7,8,9,10,11,12,13,14,15,16,17)
      }

    if (split_save){  # if writing each image to a seperate file
      for (i in 1:length(extnum)){ # for each image in the build_datacube output,

        Rfits::Rfits_write_header(filename = output_image_file_names[i], keyvalues = header_keyvalues,
                                  keycomments = header_keycomments, ext=1, create_file = T,
                                  overwrite_file = TRUE)

        image_keyvalues$BUNIT = bunits[i]
        image_keyvalues$EXTNAME = extnames[i]

        if (i < 8){ # write observed mass, velocity and dispersion images to the file
          Rfits::Rfits_write_image(data = simspin_datacube$observed_images[[which(names(simspin_datacube$observed_images) == image_names[i])]],
                                   filename = output_image_file_names[i], ext=2,
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
        } else { # and the raw particle data for metallicity and OH
          Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                   filename = output_image_file_names[i], ext=2,
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

        }

      }

    } else {  # if writing all images to the same file
      for (i in 1:length(extnum)){ # for each image in the build_datacube output,

        image_keyvalues$BUNIT = bunits[i]
        image_keyvalues$EXTNAME = extnames[i]

        # Write each subsequent image to the next HDU
        if (i < 8){ # write observed mass, velocity and dispersion images to the file
          Rfits::Rfits_write_image(data = simspin_datacube$observed_images[[which(names(simspin_datacube$observed_images) == image_names[i])]],
                                   filename = cube_file_name, ext=extnum[i],
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
        } else { # and the raw particle data for metallicity and OH
          Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                   filename = cube_file_name, ext=extnum[i],
                                   keyvalues = image_keyvalues, keycomments = image_keycomments,
                                   create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

        }

      }

    }

    if (!is.na(observation$signal_to_noise)){
      # Adding variance cube to FITS file ----
      data_keyvalues$EXTNAME = "STAT"
      data_keyvalues$BUNIT   = "1**40 (erg/s/cm**2)**-2"

      if (split_save){
        stat_summary_file_name = paste0(output_dir, "/", output_file_root, "_inv_variance_cube.FITS")

        Rfits::Rfits_write_header(filename = stat_summary_file_name,
                                  keyvalues = header_keyvalues,
                                  keycomments = header_keycomments, ext=1, create_file = T,
                                  overwrite_file = TRUE)
        Rfits::Rfits_write_cube(data = simspin_datacube$variance_cube/1e40,
                                filename = stat_summary_file_name, ext=2,
                                keyvalues = data_keyvalues, keycomments = data_keycomments,
                                create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

      } else {
        Rfits::Rfits_write_cube(data = simspin_datacube$variance_cube/1e40,
                                filename = cube_file_name, ext=(max(extnum)+1),
                                keyvalues = data_keyvalues,
                                keycomments = data_keycomments,
                                create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
      }
    }
  }


  if (!all(is.na(mask))){

    mask_keyvalues = list("XTENSION"="IMAGE", "BITPIX"=-64, "NAXIS"=2,
                          "NAXIS1"=dim(mask)[1], "NAXIS2"=dim(mask)[2],
                          "PCOUNT"=0, "GCOUNT"=1, "CRPIX1"=1, "CRVAL1"=1,
                          "CDELT1"=1, "CTYPE1"="", "CUNIT1"="", "CRPIX2"=1,
                          "CRVAL2"=1, "CDELT2"=1, "CTYPE2"="", "CUNIT2"="",
                          "EXTNAME"="MASK")

    mask_keycomments = list("XTENSION"="Image extension",
                            "BITPIX"="array data type",
                            "NAXIS"="number of array dimensions",
                            "NAXIS1"=NA, "NAXIS2"=NA,
                            "PCOUNT"="number of parameters",
                            "GCOUNT"="number of groups", "CRPIX1"=NA, "CRVAL1"=NA,
                            "CDELT1"=NA, "CTYPE1"=NA, "CUNIT1"=NA, "CRPIX2"=NA,
                            "CRVAL2"=NA, "CDELT2"=NA, "CTYPE2"=NA, "CUNIT2"=NA,
                            "EXTNAME"="Extension name")

    if (split_save){
      mask_file_name = paste0(output_dir, "/", output_file_root, "_mask.FITS")

      Rfits::Rfits_write_header(filename = mask_file_name, keyvalues = header_keyvalues,
                                keycomments = header_keycomments, ext=1, create_file = T,
                                overwrite_file = TRUE)

      Rfits::Rfits_write_image(data = mask,
                               filename = mask_file_name,
                               keyvalues = mask_keyvalues, keycomments = mask_keycomments,
                               create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
    } else {
      Rfits::Rfits_write_image(data = mask,
                               filename = cube_file_name,
                               keyvalues = mask_keyvalues, keycomments = mask_keycomments,
                               create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

    }


  }

  message("FITS files written to directory: ", output_dir)

}
