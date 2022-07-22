# Date: 22/07/2022
# Title: Testing the plot_images.R code

library(testthat)
context("Testing plotting functions.\n")

ss_pd_gadget   = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")
ss_pd_hdf5  = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_pd_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
ss_pd_magneticum = system.file("extdata", "SimSpin_example_Magneticum.hdf5", package = "SimSpin")
ss_pd_horizon = system.file("extdata", "SimSpin_example_HorizonAGN.hdf5", package = "SimSpin")

ss_gadget   = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
ss_hdf5     = make_simspin_file(ss_pd_hdf5, write_to_file = FALSE)
ss_eagle    = make_simspin_file(ss_pd_eagle, write_to_file = FALSE, template = "EMILES")
ss_magneticum = make_simspin_file(ss_pd_magneticum, write_to_file = FALSE, template = "BC03hr", sph_spawn_n=10)
ss_horizon    = make_simspin_file(ss_pd_horizon, write_to_file = FALSE, template = "BC03lr")

gadget_cube = build_datacube(simspin_file = ss_gadget,
                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                             observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                             method="velocity",
                             verbose = F)

test_that("Each image type can be plotted for build_datacube images - Gadget", {

  expect_invisible(plot_flux(gadget_cube$raw_images$flux_image))
  expect_invisible(plot_velocity(gadget_cube$raw_images$velocity_image))
  expect_invisible(plot_dispersion(gadget_cube$raw_images$dispersion_image))
  expect_invisible(plot_age(gadget_cube$raw_images$age_image))
  expect_invisible(plot_metallicity(gadget_cube$raw_images$metallicity_image))
  expect_invisible(plot_particles(gadget_cube$raw_images$particle_image))

  expect_invisible(plot_flux(gadget_cube$observed_images$flux_image))
  expect_invisible(plot_velocity(gadget_cube$observed_images$velocity_image))
  expect_invisible(plot_dispersion(gadget_cube$observed_images$dispersion_image))
  expect_invisible(plot_h3(gadget_cube$observed_images$h3_image))
  expect_invisible(plot_h4(gadget_cube$observed_images$h4_image))

})

test_that("Each image type can be plotted for build_datacube images - HDF5", {
  hdf5_cube = build_datacube(simspin_file = ss_hdf5,
                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                             observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 70, blur = T),
                             method="velocity", mass_flag=T,
                             verbose = F)

  expect_invisible(plot_mass(hdf5_cube$raw_images$mass_image))
  expect_invisible(plot_velocity(hdf5_cube$raw_images$velocity_image))
  expect_invisible(plot_dispersion(hdf5_cube$raw_images$dispersion_image))
  expect_invisible(plot_age(hdf5_cube$raw_images$age_image))
  expect_invisible(plot_metallicity(hdf5_cube$raw_images$metallicity_image))
  expect_invisible(plot_particles(hdf5_cube$raw_images$particle_image))

  expect_invisible(plot_mass(hdf5_cube$observed_images$mass_image))
  expect_invisible(plot_velocity(hdf5_cube$observed_images$velocity_image))
  expect_invisible(plot_dispersion(hdf5_cube$observed_images$dispersion_image))
  expect_invisible(plot_h3(hdf5_cube$observed_images$h3_image))
  expect_invisible(plot_h4(hdf5_cube$observed_images$h4_image))

})

test_that("Each image type can be plotted for build_datacube images - EAGLE", {
  eagle_cube = build_datacube(simspin_file = ss_eagle,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 30, blur = T),
                              method="velocity",
                              verbose = F)

  expect_invisible(plot_flux(eagle_cube$raw_images$flux_image))
  expect_invisible(plot_velocity(eagle_cube$raw_images$velocity_image))
  expect_invisible(plot_dispersion(eagle_cube$raw_images$dispersion_image))
  expect_invisible(plot_age(eagle_cube$raw_images$age_image))
  expect_invisible(plot_metallicity(eagle_cube$raw_images$metallicity_image))
  expect_invisible(plot_particles(eagle_cube$raw_images$particle_image))

  expect_invisible(plot_flux(eagle_cube$observed_images$flux_image))
  expect_invisible(plot_velocity(eagle_cube$observed_images$velocity_image))
  expect_invisible(plot_dispersion(eagle_cube$observed_images$dispersion_image))
  expect_invisible(plot_h3(eagle_cube$observed_images$h3_image))
  expect_invisible(plot_h4(eagle_cube$observed_images$h4_image))

})

test_that("Each image type can be plotted for build_datacube images - Magneticum", {
  magneticum_cube = build_datacube(simspin_file = ss_magneticum,
                                   telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                   observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 30, blur = T),
                                   method="velocity", mass_flag=T,
                                   verbose = F)

  expect_invisible(plot_mass(magneticum_cube$raw_images$mass_image))
  expect_invisible(plot_velocity(magneticum_cube$raw_images$velocity_image))
  expect_invisible(plot_dispersion(magneticum_cube$raw_images$dispersion_image))
  expect_invisible(plot_age(magneticum_cube$raw_images$age_image))
  expect_invisible(plot_metallicity(magneticum_cube$raw_images$metallicity_image))
  expect_invisible(plot_particles(magneticum_cube$raw_images$particle_image))

  expect_invisible(plot_mass(magneticum_cube$observed_images$mass_image))
  expect_invisible(plot_velocity(magneticum_cube$observed_images$velocity_image))
  expect_invisible(plot_dispersion(magneticum_cube$observed_images$dispersion_image))
  expect_invisible(plot_h3(magneticum_cube$observed_images$h3_image))
  expect_invisible(plot_h4(magneticum_cube$observed_images$h4_image))

})

test_that("Each image type can be plotted for build_datacube images - HorizonAGN", {
  horizon_cube = build_datacube(simspin_file = ss_horizon,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 80, blur = T),
                                method="spectral",
                                verbose = F)

  expect_invisible(plot_flux(horizon_cube$raw_images$flux_image))
  expect_invisible(plot_velocity(horizon_cube$raw_images$velocity_image))
  expect_invisible(plot_dispersion(horizon_cube$raw_images$dispersion_image))
  expect_invisible(plot_particles(horizon_cube$raw_images$particle_image))
})

test_that("Can we add a radii symbol to the plot without error?", {
  expect_invisible(plot_flux(gadget_cube$raw_images$flux_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_mass(gadget_cube$raw_images$flux_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_velocity(gadget_cube$raw_images$velocity_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_dispersion(gadget_cube$raw_images$dispersion_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_age(gadget_cube$raw_images$age_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_metallicity(gadget_cube$raw_images$metallicity_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_particles(gadget_cube$raw_images$particle_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_h3(gadget_cube$observed_images$h3_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
  expect_invisible(plot_h4(gadget_cube$observed_images$h4_image, radii = list("a" = 2, "b" = 4, "ang" = 90)))
})

test_that("Can we move the title on the plot without error?", {
  expect_invisible(plot_flux(gadget_cube$raw_images$flux_image, titleshift = -2))
  expect_invisible(plot_mass(gadget_cube$raw_images$flux_image, titleshift = -2))
  expect_invisible(plot_velocity(gadget_cube$raw_images$velocity_image, titleshift = -2))
  expect_invisible(plot_dispersion(gadget_cube$raw_images$dispersion_image, titleshift = -2))
  expect_invisible(plot_age(gadget_cube$raw_images$age_image, titleshift = -2))
  expect_invisible(plot_metallicity(gadget_cube$raw_images$metallicity_image, titleshift = -2))
  expect_invisible(plot_particles(gadget_cube$raw_images$particle_image, titleshift = -2))
  expect_invisible(plot_h3(gadget_cube$observed_images$h3_image, titleshift = -2))
  expect_invisible(plot_h4(gadget_cube$observed_images$h4_image, titleshift = -2))
})

test_that("Can we change the number of values on the plot without error?", {
  expect_invisible(plot_flux(gadget_cube$raw_images$flux_image, labN = 3))
  expect_invisible(plot_mass(gadget_cube$raw_images$flux_image, labN = 3))
  expect_invisible(plot_velocity(gadget_cube$raw_images$velocity_image, labN = 3))
  expect_invisible(plot_dispersion(gadget_cube$raw_images$dispersion_image, labN = 3))
  expect_invisible(plot_age(gadget_cube$raw_images$age_image, labN = 3))
  expect_invisible(plot_metallicity(gadget_cube$raw_images$metallicity_image, labN = 3))
  expect_invisible(plot_particles(gadget_cube$raw_images$particle_image, labN = 3))
  expect_invisible(plot_h3(gadget_cube$observed_images$h3_image, labN = 3))
  expect_invisible(plot_h4(gadget_cube$observed_images$h4_image, labN = 3))
})

test_that("Can we change the number of values on the plot without error?", {
  expect_invisible(plot_flux(gadget_cube$raw_images$flux_image, zlim = c(5e-10, 2e-9)))
  expect_invisible(plot_mass(gadget_cube$raw_images$flux_image, zlim = c(5e-10, 2e-9)))
  expect_invisible(plot_velocity(gadget_cube$raw_images$velocity_image, zlim = c(-100,100)))
  expect_invisible(plot_dispersion(gadget_cube$raw_images$dispersion_image, zlim = c(0,100)))
  expect_invisible(plot_age(gadget_cube$raw_images$age_image, zlim = c(0,6)))
  expect_invisible(plot_metallicity(gadget_cube$raw_images$metallicity_image, zlim = c(0.01,0.02)))
  expect_invisible(plot_particles(gadget_cube$raw_images$particle_image, zlim = c(5,10)))
  expect_invisible(plot_h3(gadget_cube$observed_images$h3_image, zlim = c(-0.2,0.1)))
  expect_invisible(plot_h4(gadget_cube$observed_images$h4_image, zlim = c(-0.1,0.15)))
})

test_that("Can we change the units on the plot without error?", {
  expect_invisible(plot_flux(gadget_cube$raw_images$flux_image, units = "test"))
  expect_invisible(plot_mass(gadget_cube$raw_images$flux_image, units = "test"))
  expect_invisible(plot_velocity(gadget_cube$raw_images$velocity_image, units = "test"))
  expect_invisible(plot_dispersion(gadget_cube$raw_images$dispersion_image, units = "test"))
  expect_invisible(plot_age(gadget_cube$raw_images$age_image, units = "test"))
  expect_invisible(plot_metallicity(gadget_cube$raw_images$metallicity_image, units = "test"))
  expect_invisible(plot_particles(gadget_cube$raw_images$particle_image, units = "test"))
  expect_invisible(plot_h3(gadget_cube$observed_images$h3_image, units = "test"))
  expect_invisible(plot_h4(gadget_cube$observed_images$h4_image, units = "test"))
})

test_that("Does a meaningful error occur if the input image is all 0's?", {
  zero_map = gadget_cube$raw_images$flux_image
  zero_map[,] = 0

  expect_error(plot_flux(zero_map))
  expect_error(plot_mass(zero_map))
  expect_error(plot_velocity(zero_map))
  expect_error(plot_dispersion(zero_map))
  expect_error(plot_age(zero_map))
  expect_error(plot_metallicity(zero_map))
  expect_error(plot_particles(zero_map))
  expect_error(plot_h3(zero_map))
  expect_error(plot_h4(zero_map))
})

unlink("Rplots.pdf")
