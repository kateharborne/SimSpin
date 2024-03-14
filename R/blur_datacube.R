# Author: Kate Harborne
# Date: 29/10/2020
# Title: Blurring the cube
#
#'A function for applying seeing conditions to the mock cube
#'
#'The purpose of this function is to apply a convolution kernel to each spatial
#' plane of the cube produced in \code{\link{build_datacube}}.
#'
#'@param datacube_output The list output from \code{\link{build_datacube}}.
#'@return Returns a new list that contains the spectral cube produced by
#' \code{\link{build_datacube}}, but with the seeing conditions applied.
#'@examples
#'ss_gadget = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata",
#'                         package = "SimSpin")
#'cube = build_datacube(simspin_file = ss_gadget,
#'                      telescope = telescope(type="SAMI"),
#'                      observing_strategy = observing_strategy())
#'blurred_cube = blur_datacube(datacube_output = cube)
#'

blur_datacube = function(datacube_output){

  # Behavior of function changes with input cube -
  # If a `spectral` cube is input:
  #   1. blur each layer of the full cube
  # If a `velocity` cube is input:
  #   1. blur each layer of the full cube
  #   2. recompute the velocity and dispersion images from the blurred cube.
  #   3. convolve other images (age/metallicity/particles) with kernel.
  # If a `gas/sf gas` cube is input:
  #   1. blur each layer of the full cube
  #   2. recompute the velocity and dispersion images from the blurred cube.
  #   3. convolve other images (SFR/metallicity/particles) with kernel.


  if ("spectral_cube" %in% names(datacube_output)){

    # Pulling out the 3D spectral datacube
    cube = datacube_output$spectral_cube
    cube[is.na(cube)] = 0
    cube_dims = dim(cube)

    # Pulling out the observation details and describing aperture edge
    observation = datacube_output$observation
    calc_region = observation$aperture_region
    aperture_region = matrix(data = calc_region, nrow = cube_dims[1], ncol=cube_dims[2])

    # 1. Blurring each plane of the spectral datacube
    blur_cube = array(data = 0.0, dim = cube_dims)
    for (spatial_plane in seq(1, dim(cube)[3])){
      blur_cube[,,spatial_plane] = ProFit::profitBruteConv(cube[,,spatial_plane], observation$psf_kernel) * aperture_region
    }

    # Returning output in same format as input
    blur_output = list("spectral_cube"    = blur_cube,
                       "observation"      = observation,
                       "raw_images"       = datacube_output$raw_images,
                       "observed_images"  = NULL,
                       "variance_cube"    = NULL)
  }

  if ("velocity_cube" %in% names(datacube_output)){

    # Pulling out the 3D velocity datacube
    cube = datacube_output$velocity_cube
    cube[is.na(cube)] = 0
    cube_dims = dim(cube)

    # Pulling out the observation details and describing aperture edge
    observation = datacube_output$observation
    calc_region = observation$aperture_region
    aperture_region = matrix(data = calc_region, nrow = cube_dims[1], ncol=cube_dims[2])

    # 1. Blurring each plane of the velocity datacube
    blur_cube = array(data = 0.0, dim = cube_dims)
    for (spatial_plane in seq(1, cube_dims[3])){
      blur_cube[,,spatial_plane] = ProFit::profitBruteConv(cube[,,spatial_plane], observation$psf_kernel) * aperture_region
    }

    # 2. Blurring the observed flux map
    if (observation$method == "velocity"){  # if your observation is of stars, blur the flux map
      blur_image = array(data = 0.0, dim = cube_dims[c(1,2)])
      blur_image = ProFit::profitBruteConv(datacube_output$observed_images$flux_image, observation$psf_kernel) * aperture_region
      datacube_output$observed_images$flux_image = blur_image
    }

    blur_mass_image = array(data = 0.0, dim = cube_dims[c(1,2)])
    blur_mass_image = ProFit::profitBruteConv(datacube_output$observed_images$mass_image, observation$psf_kernel) * aperture_region
    datacube_output$observed_images$mass_image = blur_mass_image

    if (observation$method == "gas" | observation$method == "sf gas"){  # if your observation is of gas, blur the sfr map
      blur_sfr_image = array(data = 0.0, dim = cube_dims[c(1,2)])
      blur_sfr_image = ProFit::profitBruteConv(datacube_output$observed_images$SFR_image, observation$psf_kernel) * aperture_region
      datacube_output$observed_images$SFR_image = blur_sfr_image
    }


    # Returning output in same format as input
    blur_output = list("velocity_cube"    = blur_cube,
                       "observation"      = observation,
                       "raw_images"       = datacube_output$raw_images,
                       "observed_images"  = datacube_output$observed_images,
                       "variance_cube"    = NULL)

  }

  return(blur_output)

}
