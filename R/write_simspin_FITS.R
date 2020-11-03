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
#'@param simspin_cube The spectral data cube output when \code{build_datacube()}
#' is run.
#'@param observation An object of the observation class, as generated when
#' \code{build_datacube()} is run.
#'@return Returns an .fits file that contains a the generated spectral cube and
#' relevant header describing the mock observation.
#'@examples
#'\dontrun{
#'ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
#'temp_loc = tempdir()
#'make_simspin_file(ss_eagle, output = paste(temp_loc, "spectra.fst", sep=""))
#'cube = build_datacube(simspin_file = paste(temp_loc, "spectra.fst", sep=""),
#'                      telescope = telescope(type="SAMI"),
#'                      observing_strategy = observing_strategy())
#'write_simspin_FITS(output_file = paste(temp_loc, "cube.fits", sep=""),
#'                   simspin_cube = cube$spectral_cube,
#'                   observation = cube$observation)
#'unlink(paste(temp_loc, "spectra.fst", sep=""))
#'unlink(paste(temp_loc, "cube.fits", sep=""))
#'}
#'



write_simspin_FITS = function(output_file, simspin_cube, observation){

  crpixn = dim(simspin_cube)%/%2
  crvaln = c(observation$sbin_seq[crpixn[1]]/observation$ang_size/3600,
             observation$sbin_seq[crpixn[2]]/observation$ang_size/3600,
             observation$wave_seq[crpixn[3]])
  cdeltn = c(observation$sbin_size/observation$ang_size/3600,
             observation$sbin_size/observation$ang_size/3600,
             observation$wave_res)
  ctypen = c("RA---TAN", "DEC--TAN", "AWAV")
  cunitn = c("deg", "deg", "Angstrom")
  axDat  = data.frame("crpix" = crpixn, "crval" = crvaln, "cdelt" = cdeltn,
                      "len" = dim(simspin_cube), "ctype" = ctypen,
                      "cunit" = cunitn)

  FITSio::writeFITSim(X = simspin_cube, file = output_file, crpixn = crpixn,
                      crvaln = crvaln, cdeltn = cdeltn, ctypen = ctypen,
                      cunitn = cunitn, axDat = axDat)

  message("FITS file written to: ", output_file)

}
