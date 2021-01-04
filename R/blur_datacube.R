# Author: Kate Harborne
# Date: 29/10/2020
# Title: Blurring the cube
#
#'A function for applying seeing conditions to the mock cube
#'
#'The purpose of this function is to apply a convolution kernel to each spatial
#' plane of the cube produced in \code{\link{build_datacube}}.
#'
#'@param cube The list output from \code{\link{build_datacube}}.
#'@return Returns a new list that contains the spectral cube produced by
#' \code{\link{build_datacube}}, but with the seeing conditions applied.
#'@examples
#'\dontrun{
#'ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
#'temp_loc = tempdir()
#'make_simspin_file(ss_eagle, output = paste(temp_loc, "spectra.fst", sep=""))
#'cube = build_datacube(simspin_file = paste(temp_loc, "spectra.fst", sep=""),
#'                      telescope = telescope(type="SAMI"),
#'                      observing_strategy = observing_strategy())
#'unlink(paste(temp_loc, "spectra.fst", sep=""))
#'}
#'

blur_datacube = function(cube){

  spectral_cube = cube$spectral_cube
  spectral_cube[is.na(spectral_cube)] = 0
  cube_dims = dim(spectral_cube)
  observation = cube$observation

  calc_region = observation$aperture_region
  aperture_region = matrix(data = calc_region, nrow = cube_dims[1], ncol=cube_dims[2])

  blur_cube = array(data = NA, dim = cube_dims)

  for (spatial_plane in seq(1, dim(spectral_cube)[3])){
    blur_cube[,,spatial_plane] = ProFit::profitBruteConv(spectral_cube[,,spatial_plane], observation$psf_kernel) * aperture_region
  }

  return(blur_cube)

}
