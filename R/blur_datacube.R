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
                       "observed_images"  = NULL)
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

    # Returning output in same format as input
    blur_output = list("velocity_cube"    = blur_cube,
                       "observation"      = observation,
                       "raw_images"       = datacube_output$raw_images,
                       "observed_images"  = vector(mode = "list", length=3))

    names(blur_output$observed_images) = c("flux_image", "velocity_image", "dispersion_image")

    # 2. Recompute the velocity and dispersion images from the blurred cube.
    # Initializing empty arrays
    blur_flux = array(0.0, dim = c(cube_dims[c(1,2)]))
    blur_velocity = array(0.0, dim = c(cube_dims[c(1,2)]))
    blur_dispersion = array(0.0, dim = c(cube_dims[c(1,2)]))

    # Filling array based on blurred cubes
    for (c in 1:cube_dims[1]){
      for (d in 1:cube_dims[2]){
        blur_flux[c,d]       = sum(blur_cube[c,d,])
        blur_velocity[c,d]   = .meanwt(observation$vbin_seq, blur_cube[c,d,])
        blur_dispersion[c,d] = sqrt(.varwt(observation$vbin_seq, blur_cube[c,d,], blur_velocity[c,d]))
      }
    }

    # Trimming any data blurred out of the aperture region
    blur_flux       = blur_flux * aperture_region
    blur_velocity   = blur_velocity * aperture_region
    blur_dispersion = blur_dispersion * aperture_region

    blur_output$observed_images$flux_image = blur_flux
    blur_output$observed_images$velocity_image = blur_velocity
    blur_output$observed_images$dispersion_image = blur_dispersion

    if ("mass_image" %in% names(blur_output$raw_images)){
      names(blur_output$observed_images)[which(names(blur_output$observed_images) == "flux_image")] = "mass_image"
    }
  }

  return(blur_output)

}
