# Author: Kate Harborne
# Date: 26/10/2020
# Title: write_simspin_FITS - a function for making a spectral FITS file output
#
#'A function for making a mock spectral cube FITS file output
#'
#'The purpose of this function is to write the spectral cube generated using
#' the \code{build_datacube()} function to a FITS file.
#'
#'@param output_file The path to the location where the FITS file should be
#' written.
#'@param simspin_datacube The list output from \code{build_datacube()}.
#'@param object_name A string that described the name of the observed object.
#'@param telescope_name A string that describes the name of the telescope used.
#'@param instrument_name A string that describes the used instrument on that
#' telescope.
#'@param observer_name A string that describes the name of the observer.
#'@param input_simspin_file A string describing the original SimSpin file
#' (i.e. the file output from make_simspin_file).
#'@param mask (Optional) A binary array describing the masked regions of the
#' cube/images.
#'@return Returns an .fits file that contains a the generated data and
#' relevant header describing the mock observation.
#'@examples
#'\dontrun{
#'ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
#'temp_loc = tempdir()
#'make_simspin_file(ss_eagle, output = paste(temp_loc, "spectra.Rdata", sep=""))
#'cube = build_datacube(simspin_file = paste(temp_loc, "spectra.Rdata", sep=""),
#'                      telescope = telescope(type="SAMI"),
#'                      observing_strategy = observing_strategy())
#'write_simspin_FITS(output_file = paste(temp_loc, "cube.fits", sep=""),
#'                   simspin_datacube = cube, object_name = "SimSpin EAGLE example",
#'                   telescope_name = "AAO", instrument_name = "SAMI",
#'                   observer_name = "K.E.Harborne",
#'                   input_simspin_file = paste(temp_loc, "spectra.Rdata", sep=""))
#'unlink(paste(temp_loc, "spectra.fst", sep=""))
#'unlink(paste(temp_loc, "cube.fits", sep=""))
#'}
#'

write_simspin_FITS = function(output_file, simspin_datacube, object_name,
                              telescope_name, instrument_name, observer_name,
                              input_simspin_file, mask=NA){

  observation = simspin_datacube$observation
  simspin_cube = simspin_datacube[[1]] # getting either the "velocity cube" or the "spectral cube"

  # FITS files always have the first HDU containing the header information.
  #
  # Then, depending on the type of data cube that has been constructed, images
  # and their respective header information is stored in subsequent HDUs.
  #
  # Basic header information is initiallised in this function to begin.
  # This is then modified in the following "if"s depending on the telescope
  # mode employed.

  # Making header for file and saving to HDU 1 ----
  output_name = rev(stringr::str_split(output_file, "/")[[1]])[1]

  header_keyvalues = list("SIMPLE"=TRUE, "BITPIX"=8, "NAXIS"=0, "EXTEND"=TRUE,
                          "DATE"=Sys.time(), "ORIGIN"="SimSpin", "TELESCOP"=telescope_name,
                          "INSTRUME"=instrument_name, "RA"=observation$pointing_deg[1],
                          "DEC"=observation$pointing_deg[2], "EQINOX"=2000,
                          "RADECSYS"="FK5", "EXPTIME"=1320, "MJD-OBS"=58906.11,
                          "DATE-OBS"=observation$date, "UTC"=9654, "LST"=30295.18,
                          "PI-COI"="UNKNOWN", "OBSERVER"=observer_name, "REDSHIFT"=observation$z,
                          "PIPEFILE"=output_name,
                          "BUNIT"="erg/s/cm**2",
                          "ARCFILE"=input_simspin_file,
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
                            "BUNIT"="Angstrom",
                            "ARCFILE"="Archive File Name",
                            "DATAMD5"="MD5 checksum",
                            "OBJECT"="Original target.")

  Rfits::Rfits_write_header(filename = output_file, keyvalues = header_keyvalues,
                            keycomments = header_keycomments, ext=1, create_file = T,
                            overwrite_file = TRUE)

  # Data header to be modified and added for each HDU ----

  data_keyvalues = list("XTENSION" = "IMAGE", "BITPIX"=-32, "NAXIS"=3,
                        "NAXIS1"=dim(simspin_cube)[1], "NAXIS2"=dim(simspin_cube)[2],
                        "NAXIS3"=dim(simspin_cube)[3], "PCOUNT"=0, "GCOUNT"=1,
                        "EXTNAME"="DATA", "HDUCLASS"="ESO",
                        "HDUDOC"="DICD", "HDUVERS"="DCID version 6",
                        "HDUCLAS1"="IMAGE", "HDUCLAS2"="DATA", "ERRDATA"="STAT",
                        "OBJECT"=object_name, "BUNIT"="erg/s/cm**2",
                        "CRPIX1"=1, "CRPIX2"=1,
                        "CDELT1"=-observation$sbin_size/observation$ang_size/3600,
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
                          "BUNIT"="Angstrom",
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

  # Writing data to HDUs based on the observation method employed ----

  if (observation$method == "spectral"){

    data_keyvalues$CTYPE3 = "WAVE"
    data_keyvalues$CUNIT3 = "Angstrom"
    data_keyvalues$CDELT3 = observation$wave_res
    data_keyvalues$CRVAL3 = observation$wave_seq[1]

    Rfits::Rfits_write_cube(data = simspin_cube, filename = output_file, ext=2,
                            keyvalues = data_keyvalues, keycomments = data_keycomments,
                            create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

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

    extnames = c("FLUX", "LOS_VEL", "LOS_DISP", "NPART")
    bunits = c("erg/s/cm**2", "km/s", "km/s", "Particle number")
    extnum = c(3,4,5,6)
    image_names = c("flux_image", "velocity_image", "dispersion_image", "particle_image")

    for (i in 1:4){
      image_keyvalues$BUNIT = bunits[i]
      image_keyvalues$EXTNAME = extnames[i]

      Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                               filename = output_file, ext=extnum[i],
                               keyvalues = image_keyvalues, keycomments = image_keycomments,
                               create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
    }

  }

  if (observation$method == "velocity"){

    data_keyvalues$CTYPE3 = "VELOCITY"
    data_keyvalues$CUNIT3 = "km/s"
    data_keyvalues$CDELT3  = observation$vbin_size
    data_keyvalues$CRVAL3 = observation$vbin_seq[1]

    Rfits::Rfits_write_cube(data = simspin_cube, filename = output_file, ext=2,
                            keyvalues = data_keyvalues, keycomments = data_keycomments,
                            create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

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

    extnames = if("flux_image" %in% names(simspin_datacube$raw_images)){c("FLUX", "LOS_VEL", "LOS_DISP", "AGE", "METALS", "NPART")}else{c("MASS", "LOS_VEL", "LOS_DISP", "AGE", "METALS", "NPART")}
    bunits = if("flux_image" %in% names(simspin_datacube$raw_images)){c("erg/s/cm**2", "km/s", "km/s", "Gyr", "Z_solar", "Particle number")}else{c("Msol", "km/s", "km/s", "Gyr", "Z_solar", "Particle number")}
    extnum = c(3,4,5,6,7,8)
    image_names = if("flux_image" %in% names(simspin_datacube$raw_images)){c("flux_image", "velocity_image", "dispersion_image", "age_image", "metallicity_image", "particle_image")}else{c("mass_image", "velocity_image", "dispersion_image", "age_image", "metallicity_image", "particle_image")}

    for (i in 1:6){
      image_keyvalues$BUNIT = bunits[i]
      image_keyvalues$EXTNAME = extnames[i]

      if (i < 4){
        Rfits::Rfits_write_image(data = simspin_datacube$observed_images[[which(names(simspin_datacube$observed_images) == image_names[i])]],
                                 filename = output_file, ext=extnum[i],
                                 keyvalues = image_keyvalues, keycomments = image_keycomments,
                                 create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
      } else {
        Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                 filename = output_file, ext=extnum[i],
                                 keyvalues = image_keyvalues, keycomments = image_keycomments,
                                 create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

      }

    }

  }

  if (observation$method == "gas" | observation$method == "sf gas"){

    data_keyvalues$CTYPE3 = "GAS_VELO"
    data_keyvalues$CUNIT3 = "km/s"
    data_keyvalues$CDELT3  = observation$vbin_size
    data_keyvalues$CRVAL3 = observation$vbin_seq[1]

    Rfits::Rfits_write_cube(data = simspin_cube, filename = output_file, ext=2,
                            keyvalues = data_keyvalues, keycomments = data_keycomments,
                            create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

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

    extnames = c("MASS", "LOS_VEL", "LOS_DISP", "METALS", "OH_ABUND", "SFR")
    bunits = c("Msol", "km/s", "km/s", "log10(Z/Z_solar)", "log10(O/H)+12", "Msol/year")
    extnum = c(3,4,5,6,7,8)
    image_names = c("mass_image", "velocity_image", "dispersion_image", "metallicity_image", "OH_image", "SFR_image")

    for (i in 1:6){
      image_keyvalues$BUNIT = bunits[i]
      image_keyvalues$EXTNAME = extnames[i]

      if (i < 4){ # write observed mass, velocity and dispersion images to the file
        Rfits::Rfits_write_image(data = simspin_datacube$observed_images[[which(names(simspin_datacube$observed_images) == image_names[i])]],
                                 filename = output_file, ext=extnum[i],
                                 keyvalues = image_keyvalues, keycomments = image_keycomments,
                                 create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)
      } else { # and the raw particle data for metallicity and OH
        Rfits::Rfits_write_image(data = simspin_datacube$raw_images[[which(names(simspin_datacube$raw_images) == image_names[i])]],
                                 filename = output_file, ext=extnum[i],
                                 keyvalues = image_keyvalues, keycomments = image_keycomments,
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

    Rfits::Rfits_write_image(data = mask,
                             filename = output_file,
                             keyvalues = mask_keyvalues, keycomments = mask_keycomments,
                             create_ext = TRUE, create_file = FALSE, overwrite_file = FALSE)

  }

  message("FITS file written to: ", output_file)

}
