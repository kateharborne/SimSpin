# Author: Kate Harborne
# Co-author: Alice Serene
# Date: 13/01/2023
# Title: Testing the build_datacube.R code

library(testthat)
context("Testing build_datacube function.\n")

ss_pd_gadget     = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")
ss_pd_hdf5       = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_pd_eagle      = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
ss_pd_magneticum = system.file("extdata", "SimSpin_example_Magneticum.hdf5", package = "SimSpin")
ss_pd_horizon    = system.file("extdata", "SimSpin_example_HorizonAGN.hdf5", package = "SimSpin")
ss_pd_illustris  = system.file("extdata", "SimSpin_example_IllustrisTNG.hdf5", package = "SimSpin")

ss_gadget_old = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata", package = "SimSpin")
ss_gadget     = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
ss_hdf5       = make_simspin_file(ss_pd_hdf5, template = "BC03hr", write_to_file = FALSE)
ss_eagle      = make_simspin_file(ss_pd_eagle, write_to_file = FALSE, template = "EMILES")
ss_magneticum = make_simspin_file(ss_pd_magneticum, write_to_file = FALSE, template = "BC03hr", sph_spawn_n=10)
ss_horizon    = make_simspin_file(ss_pd_horizon, write_to_file = FALSE, template = "BC03lr")
ss_illustris  = make_simspin_file(ss_pd_illustris, write_to_file = FALSE, template = "BC03lr")

temp_loc = tempdir()

built_cube_size = 5
ob_table_loc = 3
variance_loc_spectral = 11

spectra_raw_images_size = 7
spectra_observed_images_size = NULL
spectra_number_of_hdu_sntrue = 11
spectra_number_of_hdu_snfalse = 10
spectral_raw_vel_loc = 5

velocity_raw_images_size = 7
velocity_observed_images_size = 7

velocity_number_of_hdu_sntrue = 18
velocity_number_of_hdu_snfalse = 17
velocity_obs_vel_loc = 5
velocity_obs_h3_loc = 7
velocity_obs_h4_loc = 8
velocity_obs_res_loc = 9
velocity_obs_mass_loc = 10
velocity_raw_mass_loc = 12
velocity_variance_loc = 18

gas_raw_images_size = 7
gas_observed_images_size = 7
gas_number_of_hdu_snfalse = 17
gas_number_of_hdu_sntrue = 18
gas_obs_res_loc = 9
gas_obs_sfr_loc = 10

# Testing that build_datacube works in spectral mode ----

test_that("Gadget files can be built - spectral mode", {
  expect_warning(build_datacube(simspin_file = ss_gadget_old,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                verbose = F))
  gadget_spectra = build_datacube(simspin_file = ss_gadget,
                                  telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 30),
                                  observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                  verbose = T)
  expect_length(gadget_spectra, built_cube_size)
  expect_length(gadget_spectra$raw_images, spectra_raw_images_size)
  expect_null(gadget_spectra$observed_images)
  expect_false(is.null(gadget_spectra$variance_cube))
})

test_that("HDF5 files can be built - spectral mode", {
  hdf5_spectra = build_datacube(simspin_file = ss_hdf5,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T))
  expect_length(hdf5_spectra, built_cube_size)
  expect_length(hdf5_spectra$raw_images, spectra_raw_images_size)
  expect_null(hdf5_spectra$observed_images)
})

test_that("EAGLE files can be built - spectral mode and be identical in series and parallel", {
  eagle_spectra = build_datacube(simspin_file = ss_eagle,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                 observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T))
  expect_length(eagle_spectra, built_cube_size)
  expect_length(eagle_spectra$raw_images, spectra_raw_images_size)
  expect_null(eagle_spectra$observed_images)
  expect_null(eagle_spectra$variance_cube)

  eagle_parallel_spectra = build_datacube(simspin_file = ss_eagle,
                                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                          observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                          cores = 2)
  expect_length(eagle_parallel_spectra, built_cube_size)
  expect_length(eagle_parallel_spectra$raw_images, spectra_raw_images_size)
  expect_null(eagle_parallel_spectra$observed_images)
  expect_null(eagle_parallel_spectra$variance_cube)

  expect_true(all.equal(eagle_spectra$spectral_cube, eagle_parallel_spectra$spectral_cube))
  expect_true(all.equal(eagle_spectra$raw_images$flux_image, eagle_parallel_spectra$raw_images$flux_image))
  expect_true(all.equal(eagle_spectra$raw_images$velocity_image, eagle_parallel_spectra$raw_images$velocity_image))
  expect_true(all.equal(eagle_spectra$raw_images$dispersion_image, eagle_parallel_spectra$raw_images$dispersion_image))

})

test_that("Magneticum files can be built - spectral mode and be identical in series and parallel", {
  magneticum_spectra = build_datacube(simspin_file = ss_magneticum,
                                      telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                      observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T))
  expect_length(magneticum_spectra, built_cube_size)
  expect_length(magneticum_spectra$raw_images, spectra_raw_images_size)
  expect_null(magneticum_spectra$observed_images)
  expect_null(magneticum_spectra$variance_cube)

  magneticum_parallel_spectra = build_datacube(simspin_file = ss_magneticum,
                                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                               cores = 2)
  expect_length(magneticum_parallel_spectra, built_cube_size)
  expect_length(magneticum_parallel_spectra$raw_images, spectra_raw_images_size)
  expect_null(magneticum_parallel_spectra$observed_images)

  expect_true(all.equal(magneticum_spectra$spectral_cube, magneticum_parallel_spectra$spectral_cube))
  expect_true(all.equal(magneticum_spectra$raw_images$flux_image, magneticum_parallel_spectra$raw_images$flux_image))
  expect_true(all.equal(magneticum_spectra$raw_images$velocity_image, magneticum_parallel_spectra$raw_images$velocity_image))
  expect_true(all.equal(magneticum_spectra$raw_images$dispersion_image, magneticum_parallel_spectra$raw_images$dispersion_image))

})

test_that("HorizonAGN files can be built - spectral mode and be identical in series and parallel", {
  horizon_spectra = build_datacube(simspin_file = ss_horizon,
                                   telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                   observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T))
  expect_length(horizon_spectra, built_cube_size)
  expect_length(horizon_spectra$raw_images, spectra_raw_images_size)
  expect_null(horizon_spectra$observed_images)

  horizon_parallel_spectra = build_datacube(simspin_file = ss_horizon,
                                            telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                            observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                            cores = 2)
  expect_length(horizon_parallel_spectra, built_cube_size)
  expect_length(horizon_parallel_spectra$raw_images, spectra_raw_images_size)
  expect_null(horizon_parallel_spectra$observed_images)

  expect_true(all.equal(horizon_spectra$spectral_cube, horizon_parallel_spectra$spectral_cube))
  expect_true(all.equal(horizon_spectra$raw_images$flux_image, horizon_parallel_spectra$raw_images$flux_image))
  expect_true(all.equal(horizon_spectra$raw_images$velocity_image, horizon_parallel_spectra$raw_images$velocity_image))
  expect_true(all.equal(horizon_spectra$raw_images$dispersion_image, horizon_parallel_spectra$raw_images$dispersion_image))

})

test_that("IllustrisTNG files can be built - spectral mode and be identical in series and parallel", {
  illustris_spectra = build_datacube(simspin_file = ss_illustris,
                                     telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                     observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T))
  expect_length(illustris_spectra, built_cube_size)
  expect_length(illustris_spectra$raw_images, spectra_raw_images_size)
  expect_null(illustris_spectra$observed_images)

  illustris_parallel_spectra = build_datacube(simspin_file = ss_illustris,
                                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                              cores = 2)
  expect_length(illustris_parallel_spectra, built_cube_size)
  expect_length(illustris_parallel_spectra$raw_images, spectra_raw_images_size)
  expect_null(illustris_parallel_spectra$observed_images)

  expect_true(all.equal(illustris_spectra$spectral_cube, illustris_parallel_spectra$spectral_cube))
  expect_true(all.equal(illustris_spectra$raw_images$flux_image, illustris_parallel_spectra$raw_images$flux_image))
  expect_true(all.equal(illustris_spectra$raw_images$velocity_image, illustris_parallel_spectra$raw_images$velocity_image))
  expect_true(all.equal(illustris_spectra$raw_images$dispersion_image, illustris_parallel_spectra$raw_images$dispersion_image))

})

# Testing that build_datacube works in velocity mode ----
test_that("Gadget files can be built - velocity mode.", {
  gadget_velocity = build_datacube(simspin_file = ss_gadget,
                                   telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                   observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                   method = "velocity",
                                   verbose = T)
  expect_length(gadget_velocity, built_cube_size)
  expect_length(gadget_velocity$raw_images, velocity_raw_images_size)
  expect_length(gadget_velocity$observed_images, velocity_observed_images_size)

})

test_that("HDF5 files can be built - velocity mode.", {
  hdf5_velocity = build_datacube(simspin_file = ss_hdf5,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                 method = "velocity")
  expect_length(hdf5_velocity, built_cube_size)
  expect_length(hdf5_velocity$raw_images, velocity_raw_images_size)
  expect_length(hdf5_velocity$observed_images, velocity_observed_images_size)

})

test_that("EAGLE files can be built - velocity mode and be identical in series and parallel.", {
  eagle_velocity = build_datacube(simspin_file = ss_eagle,
                                  telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                  observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                  method = "velocity")
  expect_length(eagle_velocity, built_cube_size)
  expect_length(eagle_velocity$raw_images, velocity_raw_images_size)
  expect_length(eagle_velocity$observed_images, velocity_observed_images_size)

  eagle_parallel_velocity = build_datacube(simspin_file = ss_eagle,
                                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                           observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                           method = "velocity",
                                           cores = 2)
  expect_length(eagle_parallel_velocity, built_cube_size)
  expect_length(eagle_parallel_velocity$raw_images, velocity_raw_images_size)
  expect_length(eagle_parallel_velocity$observed_images, velocity_observed_images_size)

  expect_true(all.equal(eagle_velocity$velocity_cube, eagle_parallel_velocity$velocity_cube))
  expect_true(all.equal(eagle_velocity$raw_images$flux_image, eagle_parallel_velocity$raw_images$flux_image))
  expect_true(all.equal(eagle_velocity$raw_images$velocity_image, eagle_parallel_velocity$raw_images$velocity_image))
  expect_true(all.equal(eagle_velocity$raw_images$dispersion_image, eagle_parallel_velocity$raw_images$dispersion_image))
  expect_true(all.equal(eagle_velocity$observed_images$flux_image, eagle_parallel_velocity$observed_images$flux_image))
  expect_true(all.equal(eagle_velocity$observed_images$velocity_image, eagle_parallel_velocity$observed_images$velocity_image))
  expect_true(all.equal(eagle_velocity$observed_images$dispersion_image, eagle_parallel_velocity$observed_images$dispersion_image))

})

test_that("Magneticum files can be built - velocity mode and be identical in series and parallel.", {
  magneticum_velocity = build_datacube(simspin_file = ss_magneticum,
                                       telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                       observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                       method = "velocity")
  expect_length(magneticum_velocity, built_cube_size)
  expect_length(magneticum_velocity$raw_images, velocity_raw_images_size)
  expect_length(magneticum_velocity$observed_images, velocity_observed_images_size)

  magneticum_parallel_velocity = build_datacube(simspin_file = ss_magneticum,
                                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                                observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                                method = "velocity",
                                                cores = 2)
  expect_length(magneticum_parallel_velocity, built_cube_size)
  expect_length(magneticum_parallel_velocity$raw_images, velocity_raw_images_size)
  expect_length(magneticum_parallel_velocity$observed_images, velocity_observed_images_size)

  expect_true(all.equal(magneticum_velocity$velocity_cube, magneticum_parallel_velocity$velocity_cube))
  expect_true(all.equal(magneticum_velocity$raw_images$flux_image, magneticum_parallel_velocity$raw_images$flux_image))
  expect_true(all.equal(magneticum_velocity$raw_images$velocity_image, magneticum_parallel_velocity$raw_images$velocity_image))
  expect_true(all.equal(magneticum_velocity$raw_images$dispersion_image, magneticum_parallel_velocity$raw_images$dispersion_image))
  expect_true(all.equal(magneticum_velocity$observed_images$flux_image, magneticum_parallel_velocity$observed_images$flux_image))
  expect_true(all.equal(magneticum_velocity$observed_images$velocity_image, magneticum_parallel_velocity$observed_images$velocity_image))
  expect_true(all.equal(magneticum_velocity$observed_images$dispersion_image, magneticum_parallel_velocity$observed_images$dispersion_image))

})

test_that("HorizonAGN files can be built - velocity mode and be identical in series and parallel.", {
  horizon_velocity = build_datacube(simspin_file = ss_horizon,
                                    telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                    observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                    method = "velocity")
  expect_length(horizon_velocity, built_cube_size)
  expect_length(horizon_velocity$raw_images, velocity_raw_images_size)
  expect_length(horizon_velocity$observed_images, velocity_observed_images_size)
  expect_false(horizon_velocity$observation$mass_flag)

  horizon_parallel_velocity = build_datacube(simspin_file = ss_horizon,
                                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                             observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                             method = "velocity",
                                             cores = 2)
  expect_length(horizon_parallel_velocity, built_cube_size)
  expect_length(horizon_parallel_velocity$raw_images, velocity_raw_images_size)
  expect_length(horizon_parallel_velocity$observed_images, velocity_observed_images_size)
  expect_false(horizon_parallel_velocity$observation$mass_flag)

  expect_true(all.equal(horizon_velocity$velocity_cube, horizon_parallel_velocity$velocity_cube))
  expect_true(all.equal(horizon_velocity$raw_images$flux_image, horizon_parallel_velocity$raw_images$flux_image))
  expect_true(all.equal(horizon_velocity$raw_images$velocity_image, horizon_parallel_velocity$raw_images$velocity_image))
  expect_true(all.equal(horizon_velocity$raw_images$dispersion_image, horizon_parallel_velocity$raw_images$dispersion_image))
  expect_true(all.equal(horizon_velocity$observed_images$flux_image, horizon_parallel_velocity$observed_images$flux_image))
  expect_true(all.equal(horizon_velocity$observed_images$velocity_image, horizon_parallel_velocity$observed_images$velocity_image))
  expect_true(all.equal(horizon_velocity$observed_images$dispersion_image, horizon_parallel_velocity$observed_images$dispersion_image))

})

test_that("IllustrisTNG files can be built - velocity mode and be identical in series and parallel.", {
  illustris_velocity = build_datacube(simspin_file = ss_illustris,
                                      telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                      observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                      method = "velocity", mass_flag = T)
  expect_length(illustris_velocity, built_cube_size)
  expect_length(illustris_velocity$raw_images, velocity_raw_images_size)
  expect_length(illustris_velocity$observed_images, velocity_observed_images_size)
  expect_true(illustris_velocity$observation$mass_flag)

  illustris_parallel_velocity = build_datacube(simspin_file = ss_illustris,
                                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                               method = "velocity", mass_flag = T,
                                               cores = 2)
  expect_length(illustris_parallel_velocity, built_cube_size)
  expect_length(illustris_parallel_velocity$raw_images, velocity_raw_images_size)
  expect_length(illustris_parallel_velocity$observed_images, velocity_observed_images_size)
  expect_true(illustris_parallel_velocity$observation$mass_flag)

  expect_true(all.equal(illustris_velocity$velocity_cube, illustris_parallel_velocity$velocity_cube))
  expect_true(all.equal(illustris_velocity$raw_images$flux_image, illustris_parallel_velocity$raw_images$flux_image))
  expect_true(all.equal(illustris_velocity$raw_images$velocity_image, illustris_parallel_velocity$raw_images$velocity_image))
  expect_true(all.equal(illustris_velocity$raw_images$dispersion_image, illustris_parallel_velocity$raw_images$dispersion_image))
  expect_true(all.equal(illustris_velocity$observed_images$flux_image, illustris_parallel_velocity$observed_images$flux_image))
  expect_true(all.equal(illustris_velocity$observed_images$velocity_image, illustris_parallel_velocity$observed_images$velocity_image))
  expect_true(all.equal(illustris_velocity$observed_images$dispersion_image, illustris_parallel_velocity$observed_images$dispersion_image))

})

# Testing that build_datacube works in gas mode ----

test_that("EAGLE files can be built - gas mode and be identical in series and parallel.", {
  eagle_gas = build_datacube(simspin_file = ss_eagle,
                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                             observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                             method = "gas")
  expect_length(eagle_gas, built_cube_size)
  expect_length(eagle_gas$raw_images, gas_raw_images_size)
  expect_length(eagle_gas$observed_images, gas_observed_images_size)
  expect_true(eagle_gas$observation$mass_flag)

  eagle_parallel_gas = build_datacube(simspin_file = ss_eagle,
                                      telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                      observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                      method = "gas",
                                      cores = 2)
  expect_length(eagle_parallel_gas, built_cube_size)
  expect_length(eagle_parallel_gas$raw_images, gas_raw_images_size)
  expect_length(eagle_parallel_gas$observed_images, gas_observed_images_size)

  expect_true(all.equal(eagle_gas$velocity_cube, eagle_parallel_gas$velocity_cube))
  expect_true(all.equal(eagle_gas$raw_images$mass_image, eagle_parallel_gas$raw_images$mass_image))
  expect_true(all.equal(eagle_gas$raw_images$velocity_image, eagle_parallel_gas$raw_images$velocity_image))
  expect_true(all.equal(eagle_gas$raw_images$dispersion_image, eagle_parallel_gas$raw_images$dispersion_image))
  expect_true(all.equal(eagle_gas$raw_images$SFR_image, eagle_parallel_gas$raw_images$SFR_image))
  expect_true(all.equal(eagle_gas$raw_images$metallicity_image, eagle_parallel_gas$raw_images$metallicity_image))
  expect_true(all.equal(eagle_gas$raw_images$OH_image, eagle_parallel_gas$raw_images$OH_image))
  expect_true(all.equal(eagle_gas$raw_images$particle_image, eagle_parallel_gas$raw_images$particle_image))
  expect_true(all.equal(eagle_gas$observed_images$mass_image, eagle_parallel_gas$observed_images$mass_image))
  expect_true(all.equal(eagle_gas$observed_images$velocity_image, eagle_parallel_gas$observed_images$velocity_image))
  expect_true(all.equal(eagle_gas$observed_images$dispersion_image, eagle_parallel_gas$observed_images$dispersion_image))
  expect_true(all.equal(eagle_gas$observed_images$h3_image, eagle_parallel_gas$observed_images$h3_image))
  expect_true(all.equal(eagle_gas$observed_images$h4_image, eagle_parallel_gas$observed_images$h4_image))

})

test_that("Magneticum files can be built - gas mode and be identical in series and parallel.", {
  magneticum_gas = build_datacube(simspin_file = ss_magneticum,
                                  telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                  observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 0, twist_deg = 90, blur = T),
                                  method = "gas")
  expect_length(magneticum_gas, built_cube_size)
  expect_length(magneticum_gas$raw_images, gas_raw_images_size)
  expect_length(magneticum_gas$observed_images, gas_observed_images_size)

  magneticum_parallel_gas = build_datacube(simspin_file = ss_magneticum,
                                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                           observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 0, twist_deg = 90, blur = T),
                                           method = "gas",
                                           cores = 2)
  expect_length(magneticum_parallel_gas, built_cube_size)
  expect_length(magneticum_parallel_gas$raw_images, gas_raw_images_size)
  expect_length(magneticum_parallel_gas$observed_images, gas_observed_images_size)

  expect_true(all.equal(magneticum_gas$velocity_cube, magneticum_parallel_gas$velocity_cube))
  expect_true(all.equal(magneticum_gas$raw_images$mass_image, magneticum_parallel_gas$raw_images$mass_image))
  expect_true(all.equal(magneticum_gas$raw_images$velocity_image, magneticum_parallel_gas$raw_images$velocity_image))
  expect_true(all.equal(magneticum_gas$raw_images$dispersion_image, magneticum_parallel_gas$raw_images$dispersion_image))
  expect_true(all.equal(magneticum_gas$raw_images$SFR_image, magneticum_parallel_gas$raw_images$SFR_image))
  expect_true(all.equal(magneticum_gas$raw_images$metallicity_image, magneticum_parallel_gas$raw_images$metallicity_image))
  expect_true(all.equal(magneticum_gas$raw_images$OH_image, magneticum_parallel_gas$raw_images$OH_image))
  expect_true(all.equal(magneticum_gas$raw_images$particle_image, magneticum_parallel_gas$raw_images$particle_image))
  expect_true(all.equal(magneticum_gas$observed_images$mass_image, magneticum_parallel_gas$observed_images$mass_image))
  expect_true(all.equal(magneticum_gas$observed_images$velocity_image, magneticum_parallel_gas$observed_images$velocity_image))
  expect_true(all.equal(magneticum_gas$observed_images$dispersion_image, magneticum_parallel_gas$observed_images$dispersion_image))
  expect_true(all.equal(magneticum_gas$observed_images$h3_image, magneticum_parallel_gas$observed_images$h3_image))
  expect_true(all.equal(magneticum_gas$observed_images$h4_image, magneticum_parallel_gas$observed_images$h4_image))

})

test_that("HorizonAGN files can be built - gas mode and be identical in series and parallel.", {
  horizon_gas = build_datacube(simspin_file = ss_horizon,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 0, twist_deg = 90, blur = F),
                               method = "gas")
  expect_length(horizon_gas, built_cube_size)
  expect_length(horizon_gas$raw_images, gas_raw_images_size)
  expect_length(horizon_gas$observed_images, gas_observed_images_size)

  horizon_parallel_gas = build_datacube(simspin_file = ss_horizon,
                                        telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                        observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 0, twist_deg = 90, blur = F),
                                        method = "gas",
                                        cores = 2)
  expect_length(horizon_parallel_gas, built_cube_size)
  expect_length(horizon_parallel_gas$raw_images, gas_raw_images_size)
  expect_length(horizon_parallel_gas$observed_images, gas_observed_images_size)

  expect_true(all.equal(horizon_gas$velocity_cube, horizon_parallel_gas$velocity_cube))
  expect_true(all.equal(horizon_gas$raw_images$mass_image, horizon_parallel_gas$raw_images$mass_image))
  expect_true(all.equal(horizon_gas$raw_images$velocity_image, horizon_parallel_gas$raw_images$velocity_image))
  expect_true(all.equal(horizon_gas$raw_images$dispersion_image, horizon_parallel_gas$raw_images$dispersion_image))
  expect_true(all.equal(horizon_gas$raw_images$SFR_image, horizon_parallel_gas$raw_images$SFR_image))
  expect_true(all.equal(horizon_gas$raw_images$metallicity_image, horizon_parallel_gas$raw_images$metallicity_image))
  expect_true(all.equal(horizon_gas$raw_images$OH_image, horizon_parallel_gas$raw_images$OH_image))
  expect_true(all.equal(horizon_gas$raw_images$particle_image, horizon_parallel_gas$raw_images$particle_image))
  expect_true(all.equal(horizon_gas$observed_images$mass_image, horizon_parallel_gas$observed_images$mass_image))
  expect_true(all.equal(horizon_gas$observed_images$velocity_image, horizon_parallel_gas$observed_images$velocity_image))
  expect_true(all.equal(horizon_gas$observed_images$dispersion_image, horizon_parallel_gas$observed_images$dispersion_image))
  expect_true(all.equal(horizon_gas$observed_images$h3_image, horizon_parallel_gas$observed_images$h3_image))
  expect_true(all.equal(horizon_gas$observed_images$h4_image, horizon_parallel_gas$observed_images$h4_image))

})

test_that("IllustrisTNG files can be built - gas mode and be identical in series and parallel.", {
  illustris_gas = build_datacube(simspin_file = ss_illustris,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                 observing_strategy = observing_strategy(dist_z = 0.1, inc_deg = 0, twist_deg = 90, blur = F),
                                 method = "gas")
  expect_length(illustris_gas, built_cube_size)
  expect_length(illustris_gas$raw_images, gas_raw_images_size)
  expect_length(illustris_gas$observed_images, gas_observed_images_size)

  illustris_parallel_gas = build_datacube(simspin_file = ss_illustris,
                                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                          observing_strategy = observing_strategy(dist_z = 0.1, inc_deg = 0, twist_deg = 90, blur = F),
                                          method = "gas",
                                          cores = 2)
  expect_length(illustris_parallel_gas, built_cube_size)
  expect_length(illustris_parallel_gas$raw_images, gas_raw_images_size)
  expect_length(illustris_parallel_gas$observed_images, gas_observed_images_size)

  expect_true(all.equal(illustris_gas$velocity_cube, illustris_parallel_gas$velocity_cube))
  expect_true(all.equal(illustris_gas$raw_images$mass_image, illustris_parallel_gas$raw_images$mass_image))
  expect_true(all.equal(illustris_gas$raw_images$velocity_image, illustris_parallel_gas$raw_images$velocity_image))
  expect_true(all.equal(illustris_gas$raw_images$dispersion_image, illustris_parallel_gas$raw_images$dispersion_image))
  expect_true(all.equal(illustris_gas$raw_images$SFR_image, illustris_parallel_gas$raw_images$SFR_image))
  expect_true(all.equal(illustris_gas$raw_images$metallicity_image, illustris_parallel_gas$raw_images$metallicity_image))
  expect_true(all.equal(illustris_gas$raw_images$OH_image, illustris_parallel_gas$raw_images$OH_image))
  expect_true(all.equal(illustris_gas$raw_images$particle_image, illustris_parallel_gas$raw_images$particle_image))
  expect_true(all.equal(illustris_gas$observed_images$mass_image, illustris_parallel_gas$observed_images$mass_image))
  expect_true(all.equal(illustris_gas$observed_images$velocity_image, illustris_parallel_gas$observed_images$velocity_image))
  expect_true(all.equal(illustris_gas$observed_images$dispersion_image, illustris_parallel_gas$observed_images$dispersion_image))
  expect_true(all.equal(illustris_gas$observed_images$h3_image, illustris_parallel_gas$observed_images$h3_image))
  expect_true(all.equal(illustris_gas$observed_images$h4_image, illustris_parallel_gas$observed_images$h4_image))

})

# Testing that build_datacube works in sf gas mode ----

test_that("EAGLE files can be built - sf gas mode and be identical in series and parallel.", {
  eagle_sf_gas = build_datacube(simspin_file = ss_eagle,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                method = "sf gas")
  expect_length(eagle_sf_gas, built_cube_size)
  expect_length(eagle_sf_gas$raw_images, gas_raw_images_size)
  expect_length(eagle_sf_gas$observed_images, gas_observed_images_size)

  eagle_parallel_sf_gas = build_datacube(simspin_file = ss_eagle,
                                         telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                         observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                         method = "sf gas",
                                         cores = 2)
  expect_length(eagle_parallel_sf_gas, built_cube_size)
  expect_length(eagle_parallel_sf_gas$raw_images, gas_raw_images_size)
  expect_length(eagle_parallel_sf_gas$observed_images, gas_observed_images_size)

  expect_true(all.equal(eagle_sf_gas$velocity_cube, eagle_parallel_sf_gas$velocity_cube))
  expect_true(all.equal(eagle_sf_gas$raw_images$mass_image, eagle_parallel_sf_gas$raw_images$mass_image))
  expect_true(all.equal(eagle_sf_gas$raw_images$velocity_image, eagle_parallel_sf_gas$raw_images$velocity_image))
  expect_true(all.equal(eagle_sf_gas$raw_images$dispersion_image, eagle_parallel_sf_gas$raw_images$dispersion_image))
  expect_true(all.equal(eagle_sf_gas$raw_images$SFR_image, eagle_parallel_sf_gas$raw_images$SFR_image))
  expect_true(all.equal(eagle_sf_gas$raw_images$metallicity_image, eagle_parallel_sf_gas$raw_images$metallicity_image))
  expect_true(all.equal(eagle_sf_gas$raw_images$OH_image, eagle_parallel_sf_gas$raw_images$OH_image))
  expect_true(all.equal(eagle_sf_gas$raw_images$particle_image, eagle_parallel_sf_gas$raw_images$particle_image))
  expect_true(all.equal(eagle_sf_gas$observed_images$mass_image, eagle_parallel_sf_gas$observed_images$mass_image))
  expect_true(all.equal(eagle_sf_gas$observed_images$velocity_image, eagle_parallel_sf_gas$observed_images$velocity_image))
  expect_true(all.equal(eagle_sf_gas$observed_images$dispersion_image, eagle_parallel_sf_gas$observed_images$dispersion_image))
  expect_true(all.equal(eagle_sf_gas$observed_images$h3_image, eagle_parallel_sf_gas$observed_images$h3_image))
  expect_true(all.equal(eagle_sf_gas$observed_images$h4_image, eagle_parallel_sf_gas$observed_images$h4_image))

})

test_that("Magneticum files can be built - sf gas mode and be identical in series and parallel.", {
  magneticum_sf_gas = build_datacube(simspin_file = ss_magneticum,
                                     telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                     observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                     method = "sf gas")

  expect_length(magneticum_sf_gas, built_cube_size)
  expect_length(magneticum_sf_gas$raw_images, gas_raw_images_size)
  expect_length(magneticum_sf_gas$observed_images, gas_observed_images_size)

  magneticum_parallel_sf_gas = build_datacube(simspin_file = ss_magneticum,
                                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                              method = "sf gas",
                                              cores = 2)
  expect_length(magneticum_parallel_sf_gas, built_cube_size)
  expect_length(magneticum_parallel_sf_gas$raw_images, gas_raw_images_size)
  expect_length(magneticum_parallel_sf_gas$observed_images, gas_observed_images_size)

  expect_true(all.equal(magneticum_sf_gas$velocity_cube, magneticum_parallel_sf_gas$velocity_cube))
  expect_true(all.equal(magneticum_sf_gas$raw_images$mass_image, magneticum_parallel_sf_gas$raw_images$mass_image))
  expect_true(all.equal(magneticum_sf_gas$raw_images$velocity_image, magneticum_parallel_sf_gas$raw_images$velocity_image))
  expect_true(all.equal(magneticum_sf_gas$raw_images$dispersion_image, magneticum_parallel_sf_gas$raw_images$dispersion_image))
  expect_true(all.equal(magneticum_sf_gas$raw_images$SFR_image, magneticum_parallel_sf_gas$raw_images$SFR_image))
  expect_true(all.equal(magneticum_sf_gas$raw_images$metallicity_image, magneticum_parallel_sf_gas$raw_images$metallicity_image))
  expect_true(all.equal(magneticum_sf_gas$raw_images$OH_image, magneticum_parallel_sf_gas$raw_images$OH_image))
  expect_true(all.equal(magneticum_sf_gas$raw_images$particle_image, magneticum_parallel_sf_gas$raw_images$particle_image))
  expect_true(all.equal(magneticum_sf_gas$observed_images$mass_image, magneticum_parallel_sf_gas$observed_images$mass_image))
  expect_true(all.equal(magneticum_sf_gas$observed_images$velocity_image, magneticum_parallel_sf_gas$observed_images$velocity_image))
  expect_true(all.equal(magneticum_sf_gas$observed_images$dispersion_image, magneticum_parallel_sf_gas$observed_images$dispersion_image))
  expect_true(all.equal(magneticum_sf_gas$observed_images$h3_image, magneticum_parallel_sf_gas$observed_images$h3_image))
  expect_true(all.equal(magneticum_sf_gas$observed_images$h4_image, magneticum_parallel_sf_gas$observed_images$h4_image))

})

test_that("HorizonAGN files error to be built due to insufficient particle number - sf gas mode.", {
  expect_error(build_datacube(simspin_file = ss_horizon,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                              method = "sf gas",
                              cores = 1))
  expect_error(build_datacube(simspin_file = ss_horizon,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                              method = "sf gas",
                              cores = 2))
})

test_that("IllustrisTNG files error to be built due to insufficient particle number - sf gas mode.", {
  expect_error(build_datacube(simspin_file = ss_illustris,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                              method = "sf gas",
                              cores = 1))
  expect_error(build_datacube(simspin_file = ss_illustris,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                              method = "sf gas",
                              cores = 2))
})

# Testing that build_datacube errors when invalid method given ----
test_that("Error occurs when invalid method given.", {
  expect_error(build_datacube(simspin_file = ss_eagle,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                              method = "magpie"))
})

test_that("Warnings occurs when method is specified in telescope rather than build_datacube, but still makes a cube.", {
  expect_warning(build_datacube(simspin_file = ss_eagle,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, method = "velocity"),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T)))

  warn_cube = suppressWarnings(build_datacube(simspin_file = ss_eagle,
                                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, method = "velocity"),
                                              observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T)))

  expect_length(warn_cube, built_cube_size)
  remove(warn_cube)
})

test_that("Warning occurs if method is specified in BOTH telescope and build_datacube inputs.", {
  expect_warning(build_datacube(simspin_file = ss_eagle,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, method = "spectral"),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                method="velocity"))

  warn_cube = suppressWarnings(build_datacube(simspin_file = ss_eagle,
                                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, method = "spectral"),
                                              observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                              method="velocity"))

  expect_true(warn_cube$observation$method == "velocity")

})

# Testing the mass flag functionalilty ----
test_that("Data cubes can be generated using mass rather than luminosity weighting", {
  eagle_mass = build_datacube(simspin_file = ss_eagle,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                              method = "velocity",
                              write_fits = F, mass_flag = T)

  eagle_flux = build_datacube(simspin_file = ss_eagle,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                              method = "velocity",
                              write_fits = F, mass_flag = F)

  expect_length(eagle_mass, built_cube_size)
  expect_true("velocity_cube" %in% names(eagle_mass))
  expect_false("spectral_cube" %in% names(eagle_mass))
  expect_true(eagle_mass$observation$mass_flag)
  expect_false(eagle_flux$observation$mass_flag)
  expect_false(all(eagle_flux$velocity_cube == eagle_mass$velocity_cube))

})

# Testing that build_datacube works to write to FITS file ----
test_that("Data cubes can be written to a single files", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T), built_cube_size)

  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS"))

  spectral_fits = Rfits::Rfits_read_all("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS", header = T)
  expect_true(length(spectral_fits) == spectra_number_of_hdu_snfalse)
  expect_true(spectral_fits[[2]]$keyvalues$CTYPE3 == "WAVE")
  expect_true(all(dim(spectral_fits[[2]]$imDat) == c(spectral_fits[[2]]$keyvalues$NAXIS1, spectral_fits[[2]]$keyvalues$NAXIS2, spectral_fits[[2]]$keyvalues$NAXIS3)))
  expect_true(spectral_fits[[spectral_raw_vel_loc]]$keyvalues$EXTNAME == "RAW_VEL")
  expect_true(names(spectral_fits)[ob_table_loc] == "OB_TABLE")

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T), built_cube_size)

  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS"))

  spectral_fits = Rfits::Rfits_read_all("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS", header = T)
  expect_true(length(spectral_fits) == spectra_number_of_hdu_sntrue)
  expect_true(spectral_fits[[2]]$keyvalues$CTYPE3 == "WAVE")
  expect_true(all(dim(spectral_fits[[2]]$imDat) == c(spectral_fits[[2]]$keyvalues$NAXIS1, spectral_fits[[2]]$keyvalues$NAXIS2, spectral_fits[[2]]$keyvalues$NAXIS3)))
  expect_true(spectral_fits[[spectral_raw_vel_loc]]$keyvalues$EXTNAME == "RAW_VEL")
  expect_true(names(spectral_fits)[ob_table_loc] == "OB_TABLE")
  expect_true(names(spectral_fits)[variance_loc_spectral] == "STAT")

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="velocity",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))) == velocity_number_of_hdu_snfalse)
  expect_true(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))[[velocity_obs_vel_loc]]$keyvalues$EXTNAME == "OBS_VEL")
  expect_true(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))[[velocity_obs_h3_loc]]$keyvalues$EXTNAME == "OBS_H3")
  expect_true(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))[[velocity_obs_h4_loc]]$keyvalues$EXTNAME == "OBS_H4")
  expect_true(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))[[velocity_obs_res_loc]]$keyvalues$EXTNAME == "RESIDUAL")
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS")))[ob_table_loc] == "OB_TABLE")
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS")))[velocity_raw_mass_loc] == "RAW_MASS")

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 30),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="velocity", mass_flag = T,
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget_mft.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_mft.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft.FITS"))) == velocity_number_of_hdu_sntrue)
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft.FITS")))[velocity_raw_mass_loc] == "RAW_MASS")
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft.FITS")))[velocity_obs_mass_loc] == "OBS_MASS")
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft.FITS")))[velocity_variance_loc] == "STAT")
  expect_true(stringr::str_detect(Rfits::Rfits_read_header_raw(paste0(temp_loc, "/ss_gadget_mft.FITS")), "Msol"))

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 30),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="velocity", mass_flag = T,
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget_mft_voronoi.FITS"),
                               split_save=F, voronoi_bin = T, vorbin_limit = 10), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_mft_voronoi.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_voronoi.FITS"))) == velocity_number_of_hdu_sntrue+1)
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_voronoi.FITS")))[velocity_raw_mass_loc] == "RAW_MASS")
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_voronoi.FITS")))[(velocity_variance_loc+1)] == "STAT")
  expect_true(stringr::str_detect(Rfits::Rfits_read_header_raw(paste0(temp_loc, "/ss_gadget_mft_voronoi.FITS")), "Msol"))


  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 30),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="velocity", mass_flag = F,
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget_mff_snt.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_mff_snt.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mff_snt.FITS"))) == velocity_number_of_hdu_sntrue)
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mff_snt.FITS")))[velocity_raw_mass_loc] == "RAW_MASS")
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mff_snt.FITS")))[velocity_variance_loc] == "STAT")
  expect_false(stringr::str_detect(Rfits::Rfits_read_header_raw(paste0(temp_loc, "/ss_gadget_mff_snt.FITS")), "Msol"))


  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="velocity", mass_flag = T,
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget_mft_snf.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_mft_snf.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf.FITS"))) == velocity_number_of_hdu_snfalse)
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf.FITS")))[velocity_raw_mass_loc] == "RAW_MASS")
  expect_true(is.na(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf.FITS")))[velocity_variance_loc] == "STAT"))

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="velocity", mass_flag = T,
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS"),
                               split_save=F, voronoi_bin = T, vorbin_limit = 10), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS"))) == velocity_number_of_hdu_snfalse+1)
  expect_true(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS")))[velocity_raw_mass_loc] == "RAW_MASS")
  expect_true(is.na(names(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS")))[(velocity_variance_loc+1)] == "STAT"))
  expect_false(all(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS"))[["OBS_MASS"]]$imDat ==
                     Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS"))[["RAW_MASS"]]$imDat))

  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="MUSE", fov = 10, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.27, inc_deg = 60, blur = T, fwhm = 0.6),
                               method="gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_eagle.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_eagle.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_eagle.FITS"))) == gas_number_of_hdu_snfalse)
  expect_true(stringr::str_detect(Rfits::Rfits_read_header_raw(paste0(temp_loc, "/ss_eagle.FITS")), "Msol"))

  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="MUSE", fov = 10, signal_to_noise = 30),
                               observing_strategy = observing_strategy(dist_z = 0.27, inc_deg = 60, blur = T, fwhm = 0.6),
                               method="gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_eagle_snt.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_snt.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_eagle_snt.FITS"))) == gas_number_of_hdu_sntrue)

  expect_length(build_datacube(simspin_file = ss_magneticum,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="sf gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_magneticum.FITS"),
                               split_save=F, voronoi_bin = T, vorbin_limit = 10), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_magneticum.FITS"))) == (gas_number_of_hdu_snfalse+1))
  expect_true(Rfits::Rfits_read(paste0(temp_loc, "/ss_magneticum.FITS"))[[gas_obs_res_loc]]$keyvalues$EXTNAME == "RESIDUAL")
  expect_true(Rfits::Rfits_read(paste0(temp_loc, "/ss_magneticum.FITS"))[[gas_obs_sfr_loc]]$keyvalues$EXTNAME == "OBS_SFR")
  expect_false(all(Rfits::Rfits_read(paste0(temp_loc, "/ss_magneticum.FITS"))[["OBS_SFR"]]$imDat ==
                     Rfits::Rfits_read(paste0(temp_loc, "/ss_magneticum.FITS"))[["RAW_SFR"]]$imDat))

  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 10, spatial_res = 0.25, fov=50),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T, output_location = paste0(temp_loc, "/ss_hdf5.FITS")), built_cube_size)

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T, voronoi_bin = T, vorbin_limit = 10,
                               output_location = paste0(temp_loc, "/ss_gadget_voronoi.FITS")), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_voronoi.FITS")))

  vorbin_fits = Rfits::Rfits_read_all(paste0(temp_loc, "/ss_gadget_voronoi.FITS"), header = T)
  expect_true(length(vorbin_fits) == (spectra_number_of_hdu_snfalse+1))
  expect_true(vorbin_fits[[2]]$keyvalues$CTYPE3 == "WAVE")
  expect_true(all(dim(vorbin_fits[[2]]$imDat) == c(vorbin_fits[[2]]$keyvalues$NAXIS1, vorbin_fits[[2]]$keyvalues$NAXIS2, vorbin_fits[[2]]$keyvalues$NAXIS3)))
  expect_true(vorbin_fits[[spectral_raw_vel_loc]]$keyvalues$EXTNAME == "RAW_VEL")
  expect_true(names(vorbin_fits)[ob_table_loc] == "OB_TABLE")
  expect_true(vorbin_fits[[(spectra_number_of_hdu_snfalse+1)]]$keyvalues$EXTNAME == "VORONOI")

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 30),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T, voronoi_bin = T, vorbin_limit = 10,
                               output_location = paste0(temp_loc, "/ss_gadget_voronoi_sntrue.FITS")), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_voronoi_sntrue.FITS")))

  vorbin_fits_sntrue = Rfits::Rfits_read_all(paste0(temp_loc, "/ss_gadget_voronoi_sntrue.FITS"), header = T)
  expect_true(length(vorbin_fits_sntrue) == (spectra_number_of_hdu_sntrue+1))
  expect_true(vorbin_fits_sntrue[[2]]$keyvalues$CTYPE3 == "WAVE")
  expect_true(all(dim(vorbin_fits_sntrue[[2]]$imDat) == c(vorbin_fits_sntrue[[2]]$keyvalues$NAXIS1,
                                                          vorbin_fits_sntrue[[2]]$keyvalues$NAXIS2,
                                                          vorbin_fits_sntrue[[2]]$keyvalues$NAXIS3)))
  expect_true(vorbin_fits_sntrue[[spectral_raw_vel_loc]]$keyvalues$EXTNAME == "RAW_VEL")
  expect_true(names(vorbin_fits_sntrue)[ob_table_loc] == "OB_TABLE")
  expect_true(vorbin_fits_sntrue[[spectra_number_of_hdu_sntrue]]$keyvalues$EXTNAME == "VORONOI")

})

unlink(c("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS",
         paste0(temp_loc, "/ss_gadget.FITS"),
         paste0(temp_loc, "/ss_gadget_mft.FITS"),
         paste0(temp_loc, "/ss_gadget_mft_voronoi.FITS"),
         paste0(temp_loc, "/ss_gadget_mff_snt.FITS"),
         paste0(temp_loc, "/ss_gadget_mft_snf.FITS"),
         paste0(temp_loc, "/ss_gadget_mft_snf_voronoi.FITS"),
         paste0(temp_loc, "/ss_eagle.FITS"),
         paste0(temp_loc, "/ss_eagle_snt.FITS"),
         paste0(temp_loc, "/ss_magenticum.FITS"),
         paste0(temp_loc, "/ss_hdf5.FITS"),
         paste0(temp_loc, "/ss_gadget_voronoi.FITS"),
         paste0(temp_loc, "/ss_gadget_voronoi_sntrue.FITS")))

test_that("Data cubes can be written to multiple files", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T, split_save=T, voronoi = T, vorbin_limit = 10), built_cube_size)

  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_observation_summary.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_flux_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_mass_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_velocity_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_dispersion_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_age_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_metallicity_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_particle_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_raw_voronoi_bins.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_inv_variance_cube.FITS"))

  spectral_fits = Rfits::Rfits_read("GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS")
  expect_true(length(spectral_fits) == 2)
  expect_true(spectral_fits[[2]]$keyvalues$CTYPE3 == "WAVE")
  expect_true(all(dim(spectral_fits[[2]]$imDat) == c(spectral_fits[[2]]$keyvalues$NAXIS1, spectral_fits[[2]]$keyvalues$NAXIS2, spectral_fits[[2]]$keyvalues$NAXIS3)))

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="velocity",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_velocity_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_observation_summary.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_obs_flux_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_obs_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_obs_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_obs_h3_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_obs_h4_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_raw_flux_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_raw_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_raw_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_raw_age_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_raw_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_raw_particle_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_inv_variance_cube.FITS")))

  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_eagle.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_gas_velocity_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_observation_summary.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_obs_mass_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_obs_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_obs_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_obs_h3_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_obs_h4_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_raw_mass_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_raw_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_raw_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_raw_OH_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_raw_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_raw_SFR_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_raw_particle_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_inv_variance_cube.FITS")))

  expect_length(build_datacube(simspin_file = ss_magneticum,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method="sf gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_magneticum.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_gas_velocity_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_observation_summary.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_obs_mass_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_obs_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_obs_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_obs_h3_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_obs_h4_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_raw_mass_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_raw_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_raw_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_raw_OH_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_raw_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_raw_SFR_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_raw_particle_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_inv_variance_cube.FITS")))

  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T, output_location = paste0(temp_loc, "/ss_hdf5.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_spectral_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_observation_summary.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_raw_flux_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_raw_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_raw_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_raw_age_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_raw_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_raw_particle_image.FITS")))
  expect_false(file.exists(paste0(temp_loc, "/ss_hdf5_inv_variance_cube.FITS")))

})

unlink(c("GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_observation_summary.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_flux_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_mass_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_velocity_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_dispersion_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_age_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_metallicity_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_particle_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_raw_voronoi_bins.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_inv_variance_cube.FITS",

         paste0(temp_loc, "/ss_gadget_velocity_cube.FITS"),
         paste0(temp_loc, "/ss_gadget_observation_summary.FITS"),
         paste0(temp_loc, "/ss_gadget_obs_flux_image.FITS"),
         paste0(temp_loc, "/ss_gadget_obs_velocity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_obs_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_gadget_obs_h3_image.FITS"),
         paste0(temp_loc, "/ss_gadget_obs_h4_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_flux_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_mass_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_velocity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_age_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_particle_image.FITS"),
         paste0(temp_loc, "/ss_gadget_inv_variance_cube.FITS"),

         paste0(temp_loc, "/ss_eagle_gas_velocity_cube.FITS"),
         paste0(temp_loc, "/ss_eagle_observation_summary.FITS"),
         paste0(temp_loc, "/ss_eagle_obs_mass_image.FITS"),
         paste0(temp_loc, "/ss_eagle_obs_velocity_image.FITS"),
         paste0(temp_loc, "/ss_eagle_obs_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_eagle_obs_h3_image.FITS"),
         paste0(temp_loc, "/ss_eagle_obs_h4_image.FITS"),
         paste0(temp_loc, "/ss_eagle_raw_mass_image.FITS"),
         paste0(temp_loc, "/ss_eagle_raw_velocity_image.FITS"),
         paste0(temp_loc, "/ss_eagle_raw_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_eagle_raw_OH_image.FITS"),
         paste0(temp_loc, "/ss_eagle_raw_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_eagle_raw_SFR_image.FITS"),
         paste0(temp_loc, "/ss_eagle_raw_particle_image.FITS"),
         paste0(temp_loc, "/ss_eagle_inv_variance_cube.FITS"),

         paste0(temp_loc, "/ss_magenticum_gas_velocity_cube.FITS"),
         paste0(temp_loc, "/ss_magenticum_observation_summary.FITS"),
         paste0(temp_loc, "/ss_magneticum_obs_mass_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_obs_velocity_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_obs_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_obs_h3_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_obs_h4_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_raw_mass_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_raw_velocity_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_raw_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_raw_OH_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_raw_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_raw_SFR_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_raw_particle_image.FITS"),
         paste0(temp_loc, "/ss_magneticum_inv_variance_cube.FITS"),

         paste0(temp_loc, "/ss_hdf5_spectral_cube.FITS"),
         paste0(temp_loc, "/ss_hdf5_observation_summary.FITS"),
         paste0(temp_loc, "/ss_hdf5_raw_flux_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_raw_mass_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_raw_velocity_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_raw_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_raw_age_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_raw_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_raw_particle_image.FITS")
))

test_that("Mask can be included in FITS files correctly", {
  cube = build_datacube(simspin_file = ss_gadget,
                        telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                        observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                        write_fits = F)

  mask = cube$raw_images$particle_image
  mask[which(cube$raw_images$particle_image <= 1)] = NA

  write_simspin_FITS(output_file = paste0(temp_loc, "/ss_gadget.FITS"),
                     simspin_datacube = cube, split_save = F,
                     mask = mask, object_name = "ss_gadget",
                     telescope_name = "SimSpin",
                     instrument_name = "SAMI", observer_name = "K Harborne",
                     input_simspin_file = "ss_gadget.Rdata")

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget.FITS")))
  fits_w_mask = Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))
  expect_true(length(fits_w_mask) == (spectra_number_of_hdu_sntrue+1))
  expect_true(fits_w_mask[[(spectra_number_of_hdu_sntrue+1)]]$keyvalues$EXTNAME == "MASK")

  write_simspin_FITS(output_file = paste0(temp_loc, "/ss_gadget.FITS"),
                     simspin_datacube = cube, split_save = T,
                     mask = mask, object_name = "ss_gadget",
                     telescope_name = "SimSpin",
                     instrument_name = "SAMI", observer_name = "K Harborne",
                     input_simspin_file = "ss_gadget.Rdata")

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_mask.FITS")))
})

unlink(c(paste0(temp_loc, "/ss_gadget.FITS"),
         paste0(temp_loc, "/ss_gadget_spectral_cube.FITS"),
         paste0(temp_loc, "/ss_gadget_observation_summary.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_flux_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_mass_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_velocity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_age_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_raw_particle_image.FITS"),
         paste0(temp_loc, "/ss_gadget_mask.FITS")))

test_that("FITS files will be written with automatic names at directory given by `output_location` if only PATH is specified", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               write_fits = T, split_save=T,
                               output_location = temp_loc), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_observation_summary.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_flux_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_mass_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_age_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_particle_image.FITS")))
})

unlink(c(paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_observation_summary.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_flux_image.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_mass_image.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_velocity_image.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_dispersion_image.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_age_image.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_metallicity_image.FITS"),
         paste0(temp_loc,"/GalaxyID_unknown_inc45deg_seeing2fwhm_raw_particle_image.FITS")))


# Testing that build_datacube will give warning if the spectra given is low res ----
test_that("build_datacube issues warning when spectral resolution < LSF fwhm.", {
  expect_warning(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="IFU", lsf_fwhm = 0.9),
                                observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45)))
  expect_warning(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="IFU", lsf_fwhm = 3.08),
                                observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45)))
})

test_that("build_datacube issues warning when wavelength resolution of telescope is < templates.", {
  expect_warning(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="IFU", lsf_fwhm = 3.6, wave_res = 0.5),
                                observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45)))
})

# Testing that the velocity shift functions work as expected -------------------
test_that("velocity shift for wavelengths work correctly", {

  wavelength = SimSpin::EMILES$Wave
  velocity_los = c(27.04932, 40.94573)
  wave = matrix(data = rep(wavelength, length(velocity_los)), nrow = length(velocity_los), byrow=T)
  wave_shift = ((velocity_los / .speed_of_light) * wave) + wave # using doppler formula to compute the shift in wavelengths cause by LOS velocity

  wave_shift_comp = matrix(data=NA, nrow=2, ncol=length(wavelength))

  for (i in 1:2){
    wave_shift_comp[i,] = ((velocity_los[i] / .speed_of_light) * wave[i,]) + wave[i,]
  }

  expect_equal(wave_shift, wave_shift_comp)

})

# Test that the spectra pulled for each particle are correct -------------------
test_that("Repeated spectra are included in intrinsic spectra", {
  simspin_data = ss_gadget
  galaxy_data = simspin_data$star_part

  obs_data = obs_galaxy(part_data = galaxy_data, inc_rad = 0.7853982) # projecting the galaxy to given inclination
  galaxy_data$x = obs_data$x;   galaxy_data$y = obs_data$y;   galaxy_data$z = obs_data$z
  galaxy_data$vx = obs_data$vx; galaxy_data$vy = obs_data$vy; galaxy_data$vz = obs_data$vz

  sbin_seq = c(-7.502401, -7.002241, -6.502080, -6.001920, -5.501760, -5.001600,
               -4.501440, -4.001280, -3.501120, -3.000960, -2.500800, -2.000640,
               -1.500480, -1.000320, -0.500160,  0.000000,  0.500160,  1.000320,
               1.500480, 2.000640, 2.500800, 3.000960, 3.501120, 4.001280, 4.501440,
               5.001600, 5.501760, 6.001920, 6.502080, 7.002241, 7.502401)
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=sbin_seq, labels=F) +
    (30 * cut(galaxy_data$z, breaks=sbin_seq, labels=F)) - (30) # assigning particles to positions in cube

  i = 411

  particle_IDs = which(galaxy_data$pixel_pos == i)
  galaxy_sample = galaxy_data[particle_IDs,]

  intrinsic_spectra = array(data = 0.0, dim = c(5, 842))
  for (n in 1:length(particle_IDs)){
    intrinsic_spectra[n,] = .spectra(SW=simspin_data$spectral_weights[[galaxy_sample$sed_id[n]]], Template = SimSpin::BC03lr)
  }

  spectra = intrinsic_spectra * (galaxy_sample$Initial_Mass * 1e10) # reading relevant spectra

  expect_equal((intrinsic_spectra[1,] * galaxy_sample$Initial_Mass[1] * 1e10), spectra[1,], tolerance = 0.001)
  expect_equal((intrinsic_spectra[2,] * galaxy_sample$Initial_Mass[2] * 1e10), spectra[2,], tolerance = 0.001)
})

# Test that the output of blurred and un-blurred format is the same ------------
test_that("Format of blurring output is the same as unblurred output", {

  unblurred = build_datacube(simspin_file = ss_hdf5,
                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                             observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = F))

  blurred = build_datacube(simspin_file = ss_hdf5,
                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                           observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T))

  expect_equal(typeof(blurred$spectral_cube), "double")
  expect_equal(typeof(blurred$observation), "list")
  expect_equal(typeof(blurred$raw_images$velocity_image), "double")
  expect_equal(typeof(blurred$raw_images$dispersion_image), "double")
  expect_equal(typeof(unblurred$spectral_cube), "double")
  expect_equal(typeof(unblurred$observation), "list")
  expect_equal(typeof(unblurred$raw_images$velocity_image), "double")
  expect_equal(typeof(unblurred$raw_images$dispersion_image), "double")

})

# Testing the twisting and inclination changes make sense ----------------------
test_that("Twisting and inclination work as expected", {
  SAMI = telescope(type="SAMI")
  SAMI$lsf_fwhm = 3.6
  strategy = SimSpin::observing_strategy(dist_z = 0.03, inc_deg = 90, twist_deg = 0) # viewing from the front
  observation = SimSpin::observation(SAMI, strategy, method = "spectral")
  twisted_data = twist_galaxy(ss_eagle$star_part, twist_rad = observation$twist_rad)
  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  front = galaxy_data$vy

  strategy = SimSpin::observing_strategy(dist_z = 0.03, inc_deg = 90, twist_deg = 180) # viewing from the back
  observation = SimSpin::observation(SAMI, strategy, method = "spectral")
  twisted_data = twist_galaxy(ss_eagle$star_part, twist_rad = observation$twist_rad)
  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  back = galaxy_data$vy

  expect_equal(front, -1*(back)) # velocities should be equal but opposite signs
})

# Test the observations get dimmer with distance added -------------------------
test_that("Observations get dimmer with increasing redshift", {
  cube_near = build_datacube(simspin_file = ss_gadget,
                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                             observing_strategy = observing_strategy(dist_z = 0.01, inc_deg = 45, blur = F),
                             verbose = F)

  cube_far  = build_datacube(simspin_file = ss_gadget,
                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                             observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = F),
                             verbose = F)

  expect_true(sum(cube_near$raw_images$flux_image, na.rm=T) > sum(cube_far$raw_images$flux_image, na.rm=T))

})

# Test that velocity cubes can be built in mass mode ---------------------------
test_that("EAGLE cubes can be built with mass weighting rather than luminosity", {
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method = "velocity",
                               mass_flag = T), built_cube_size)
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method = "velocity",
                               mass_flag = T, cores=2), built_cube_size)

  expect_true(all(build_datacube(simspin_file = ss_eagle,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                 method = "velocity",
                                 mass_flag = T)$flux_image ==
                    build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                   observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                   method = "velocity",
                                   mass_flag = T, cores=2)$flux_image, na.rm=T))

})

test_that("Mass/flux images are different for the same observing conditions.", {
  expect_true(all((build_datacube(simspin_file = ss_eagle,
                                  telescope = telescope(type="IFU",
                                                        lsf_fwhm = 3.6, signal_to_noise = NA),
                                  observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = F),
                                  method = "velocity",
                                  mass_flag = T)$flux_image) !=
                    (build_datacube(simspin_file = ss_eagle,
                                    telescope = telescope(type="IFU",
                                                          lsf_fwhm = 3.6, signal_to_noise = NA),
                                    observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = F),
                                    method = "velocity",
                                    mass_flag = F)$flux_image), na.rm =T))

})

# Testing that the pointing element in build_datacube works effectively --------
test_that("Pointing description works effectively", {

  centred_gal = build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type = "SAMI"),
                               observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(0,0)),
                               method="velocity",
                               mass_flag = T)

  centre_index = which(centred_gal$raw_images$mass_image == max(centred_gal$raw_images$mass_image, na.rm=T))

  shifted_up_gal = build_datacube(simspin_file = ss_gadget,
                                  telescope = telescope(type = "SAMI"),
                                  observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(0,5)),
                                  method="velocity",
                                  mass_flag = T)

  up_index = which(shifted_up_gal$raw_images$mass_image == max(shifted_up_gal$raw_images$mass_image, na.rm=T))

  shifted_up_gal_deg = build_datacube(simspin_file = ss_gadget,
                                      telescope = telescope(type = "SAMI"),
                                      observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_deg = c(0,0.00139)),
                                      method="velocity",
                                      mass_flag = T)

  up_index_deg = which(shifted_up_gal_deg$raw_images$mass_image == max(shifted_up_gal_deg$raw_images$mass_image, na.rm=T))

  shifted_down_gal = build_datacube(simspin_file = ss_gadget,
                                    telescope = telescope(type = "SAMI"),
                                    observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(0,-5)),
                                    method="velocity",
                                    mass_flag = T)

  down_index = which(shifted_down_gal$raw_images$mass_image == max(shifted_down_gal$raw_images$mass_image, na.rm=T))

  shifted_left_gal = build_datacube(simspin_file = ss_gadget,
                                    telescope = telescope(type = "SAMI"),
                                    observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(-5,0)),
                                    method="velocity",
                                    mass_flag = T)

  left_index =which(shifted_left_gal$raw_images$mass_image == max(shifted_left_gal$raw_images$mass_image, na.rm=T))

  shifted_right_gal = build_datacube(simspin_file = ss_gadget,
                                     telescope = telescope(type = "SAMI"),
                                     observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(5,0)),
                                     method="velocity",
                                     mass_flag = T)

  right_index = which(shifted_right_gal$raw_images$mass_image == max(shifted_right_gal$raw_images$mass_image, na.rm=T))

  expect_true((centre_index - up_index) == -300) # check all centres have shifted by the expected number of pixels
  expect_true((centre_index - down_index) == 300)
  expect_true((centre_index - right_index) == -10)
  expect_true((centre_index - left_index) == 10)
  expect_true(up_index == up_index_deg) # check that pointing is the same if spefified in degrees instead
})

# Testing that all images are produced with 0's where no particles are present (rather than NAs!)
test_that("No NAs or NaNs are included in mock observed images", {

  horizon_vel = build_datacube(simspin_file = ss_horizon,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 0, twist_deg = 90, blur = T),
                               method = "velocity")

  expect_false(any(is.na(horizon_vel$observed_images$velocity_image)))
  expect_false(any(is.na(horizon_vel$observed_images$dispersion_image)))
})

# Testing that the build_datacube function will work with old SimSpin files -----
test_that("SimSpin files < v2.3.0 do not stop the code from working!!!", {

  ss_gadget_v230 = ss_gadget[2:4]
  ss_hdf5_v230 = ss_hdf5[2:4]
  ss_eagle_v230 = ss_eagle[2:4]

  ss_gadget_v230$wave = SimSpin::BC03lr$Wave[1:842]
  ss_hdf5_v230$wave = SimSpin::BC03lr$Wave[1:6521]
  ss_eagle_v230$wave = SimSpin::EMILES$Wave[1:20356]


  expect_warning(build_datacube(simspin_file = ss_gadget_v230,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                verbose = F))

  BC03lr_cube = suppressWarnings(build_datacube(simspin_file = ss_gadget_v230,
                                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                                verbose = F))
  expect_length(BC03lr_cube, built_cube_size)


  expect_warning(build_datacube(simspin_file = ss_hdf5_v230,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                verbose = F))

  BC03hr_cube = suppressWarnings(build_datacube(simspin_file = ss_hdf5_v230,
                                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                                verbose = F))
  expect_length(BC03hr_cube, built_cube_size)

  expect_warning(build_datacube(simspin_file = ss_eagle_v230,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                verbose = F))

  EMILES_cube = suppressWarnings(build_datacube(simspin_file = ss_eagle_v230,
                                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                                verbose = F))
  expect_length(EMILES_cube, built_cube_size)


  ss_gadget_v230$wave = ss_gadget_v230$wave[1:840]

  expect_error(build_datacube(simspin_file = ss_gadget_v230,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                              observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                              verbose = F))

})

test_that("SimSpin files < 2.6.0 do not stop the code from working!!!", {
  # SimSpin files before v2.6.0 include full spectra in the file - using
  # .spectra() function to recreate this and test for a hydro and N-body model

  ss_gadget_v260 = ss_gadget
  ss_gadget_v260$spectra = vector(mode="list", length = length(ss_gadget_v260$spectral_weights))
  for (part in 1:length(ss_gadget_v260$spectral_weights)){
    ss_gadget_v260$spectra[[part]] = .spectra(ss_gadget_v260$spectral_weights[[part]], SimSpin::BC03lr)
  }
  ss_gadget_v260 = ss_gadget_v260[names(ss_gadget_v260) != "spectral_weights"]
  ss_gadget_v260$wave = SimSpin::BC03lr$Wave

  ss_eagle_v260 = ss_eagle
  ss_eagle_v260$spectra = vector(mode="list", length = length(ss_eagle_v260$spectral_weights))
  for (part in 1:length(ss_eagle_v260$spectral_weights)){
    ss_eagle_v260$spectra[[part]] = .spectra(ss_eagle_v260$spectral_weights[[part]], SimSpin::EMILES)
  }
  ss_eagle_v260 = ss_eagle_v260[names(ss_eagle_v260) != "spectral_weights"]
  ss_eagle_v260$wave = SimSpin::BC03lr$EMILES

  expect_warning(build_datacube(simspin_file = ss_gadget_v260,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                verbose = F))

  nbody_cube  = suppressWarnings(build_datacube(simspin_file = ss_gadget_v260,
                                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                                verbose = F))
  expect_length(nbody_cube, built_cube_size)


  expect_warning(build_datacube(simspin_file = ss_eagle_v260,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                verbose = F))

  hydro_cube  = suppressWarnings(build_datacube(simspin_file = ss_eagle_v260,
                                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3, wave_res = 1.06),
                                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                                verbose = F))
  expect_length(hydro_cube, built_cube_size)

})

# Testing that flux conservation is effective and works as expected ------------
test_that("Flux conservation works as expected", {

  BC03_test = make_simspin_file(filename = ss_pd_gadget, disk_age = 5, disk_Z = 0.004, template = "BC03", write_to_file = F)
  EMILES_test = make_simspin_file(filename = ss_pd_gadget, disk_age = 5, disk_Z = 0.004, template = "EMILES", write_to_file = F)

  BC03_temp = SimSpin::BC03lr
  BC03_test$star_part$sed_id = rep(1, length(BC03_test$star_part$sed_id))
  BC03_test$star_part$Initial_Mass = rep(200006, length(BC03_test$star_part$sed_id))
  BC03_test$spectra = .spectra(BC03_test$spectral_weights[[1]], BC03_temp)

  EMILES_temp = SimSpin::EMILES
  EMILES_test$star_part$sed_id = rep(1, length(EMILES_test$star_part$sed_id))
  EMILES_test$star_part$Initial_Mass = rep(200006, length(EMILES_test$star_part$sed_id))
  EMILES_test$spectra = .spectra(EMILES_test$spectral_weights[[1]], EMILES_temp)

  wave_range = c(3700,5700)
  BC03_wave_int = which(BC03_temp$Wave > wave_range[1] & BC03_temp$Wave < wave_range[2])
  EMILES_wave_int = which(EMILES_temp$Wave > wave_range[1] & EMILES_temp$Wave < wave_range[2])
  raw_diff_perc = (sum(BC03_test$spectra[BC03_wave_int] * .qdiff(BC03_temp$Wave)[BC03_wave_int])/
                     sum(EMILES_test$spectra[EMILES_wave_int] * .qdiff(EMILES_temp$Wave)[EMILES_wave_int]))*100

  #magplot(EMILES_temp$Wave[EMILES_wave_int], EMILES_test$spectra[EMILES_wave_int], type="l", col = "blue", lwd=2)
  #lines(BC03_temp$Wave[BC03_wave_int], BC03_test$spectra[BC03_wave_int], col = "red", lwd=2)

  BC03_obs   = suppressWarnings(build_datacube(simspin_file = BC03_test,
                                               telescope = telescope(type="IFU", lsf_fwhm = 0, signal_to_noise = NA, wave_res = 1.06),
                                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = F),
                                               method = "spectral", write_fits = F))

  EMILES_obs = suppressWarnings(build_datacube(simspin_file = EMILES_test,
                                               telescope = telescope(type="IFU", lsf_fwhm = 0, signal_to_noise = NA, wave_res = 1.06),
                                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = F),
                                               method = "spectral", write_fits = F))

  BC03_obs_range   = which(BC03_obs$observation$wave_seq > wave_range[1] & BC03_obs$observation$wave_seq < wave_range[2])
  EMILES_obs_range = which(EMILES_obs$observation$wave_seq > wave_range[1] & EMILES_obs$observation$wave_seq < wave_range[2])

  #magplot(EMILES_obs$observation$wave_seq[EMILES_obs_range], EMILES_obs$spectral_cube[15,15,][EMILES_obs_range], type="l", col = "blue", lwd=2)
  #lines(BC03_obs$observation$wave_seq[BC03_obs_range], BC03_obs$spectral_cube[15,15,][BC03_obs_range], col = "red", lwd=2)

  obs_diff_perc = matrix(NA, nrow=30, ncol=30)
  for (a in 1:30){
    for (b in 1:30){
      obs_diff_perc[a,b] = (sum(BC03_obs$spectral_cube[a,b,][BC03_obs_range] * .qdiff(BC03_obs$observation$wave_seq[BC03_obs_range]))/
                              sum(EMILES_obs$spectral_cube[a,b,][EMILES_obs_range] * .qdiff(EMILES_obs$observation$wave_seq[EMILES_obs_range])))*100

    }
  }

  frac = raw_diff_perc - mean(obs_diff_perc, na.rm=T)

  expect_true(frac < 1) # fluxes consistent within 1%

})

# Warning is issued when using SimSpin files older than 2.3.16 in method = "gas" mode
test_that("Warning is issued when using SimSpin files older than 2.3.16 in method = 'gas' or 'sf gas' mode", {

  ss_eagle$header$Origin = "SimSpin_v2.3.14"

  expect_warning(build_datacube(simspin_file = ss_eagle,
                                telescope = telescope(type = "IFU", fov = 3, lsf_fwhm = 3),
                                observing_strategy = observing_strategy(),
                                method = "gas"))

  expect_warning(build_datacube(simspin_file = ss_eagle,
                                telescope = telescope(type = "IFU", fov = 3, lsf_fwhm = 3),
                                observing_strategy = observing_strategy(),
                                method = "sf gas"))

})

# Error is given if you try to build an observation with a simspin file that doesn't contain those particles
test_that("Error is given if you try to build an observation with a simspin file that doesn't contain those particles", {

  ss_eagle_gasonly = list("header"    = ss_eagle$header,
                          "star_part" = NULL,
                          "gas_part"  = ss_eagle$gas_part,
                          "spectral_weights" = ss_eagle$spectral_weights)

  ss_eagle_nsfgasonly = list("header"    = ss_eagle$header,
                             "star_part" = NULL,
                             "gas_part"  = ss_eagle$gas_part[ss_eagle$gas_part$SFR == 0,],
                             "spectral_weights" = ss_eagle$spectral_weights)

  expect_error(build_datacube(simspin_file = ss_eagle_gasonly,
                              telescope = telescope(type = "IFU", fov = 3, lsf_fwhm = 3),
                              observing_strategy = observing_strategy(),
                              method = "velocity"))

  expect_error(build_datacube(simspin_file = ss_eagle_gasonly,
                              telescope = telescope(type = "IFU", fov = 3, lsf_fwhm = 3),
                              observing_strategy = observing_strategy(),
                              method = "spectral"))

  expect_error(build_datacube(simspin_file = ss_gadget,
                              telescope = telescope(type = "IFU", fov = 3, lsf_fwhm = 3),
                              observing_strategy = observing_strategy(),
                              method = "gas"))

  expect_error(build_datacube(simspin_file = ss_gadget,
                              telescope = telescope(type = "IFU", fov = 3, lsf_fwhm = 3),
                              observing_strategy = observing_strategy(),
                              method = "sf gas"))

  expect_error(build_datacube(simspin_file = ss_eagle_nsfgasonly,
                              telescope = telescope(type = "IFU", fov = 3, lsf_fwhm = 3),
                              observing_strategy = observing_strategy(),
                              method = "sf gas"))
})

# Test that the LOSVD reaches zero at all spatial bins either end of the velocity scale
test_that("The LOSVD is fully sampled by the velocity bins at all spaxels", {

  eagle_velocity = build_datacube(simspin_file = ss_eagle,
                                  telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                  observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                  method = "velocity")

  expect_equal(sum(eagle_velocity$velocity_cube[,,1]) + sum(eagle_velocity$velocity_cube[,,eagle_velocity$observation$vbin]), 0)

})

test_that("Noise increases as expected in velocity cubes", {

  gadget_velocity_nonoise = build_datacube(simspin_file = ss_gadget,
                                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                           observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                           method = "velocity")

  gadget_velocity_sn30    = build_datacube(simspin_file = ss_gadget,
                                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 30),
                                           observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                           method = "velocity")

  gadget_velocity_sn5     = build_datacube(simspin_file = ss_gadget,
                                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 5),
                                           observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                           method = "velocity")

  expect_true(sd(gadget_velocity_sn30$observed_images$flux_image - gadget_velocity_nonoise$observed_images$flux_image, na.rm = T) <
                sd(gadget_velocity_sn5$observed_images$flux_image - gadget_velocity_nonoise$observed_images$flux_image, na.rm = T))

})

test_that("Noise increases as expected in spectral cubes", {

  gadget_spectra_nonoise = build_datacube(simspin_file = ss_gadget,
                                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA),
                                          observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                          method = "spectral")

  gadget_spectra_sn30    = build_datacube(simspin_file = ss_gadget,
                                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 30),
                                          observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                          method = "spectral")

  gadget_spectra_sn5     = build_datacube(simspin_file = ss_gadget,
                                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 5),
                                          observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                                          method = "spectral")

  expect_true(sd(gadget_spectra_sn30$spectral_cube[15,15,] - gadget_spectra_nonoise$spectral_cube[15,15,]) <
                sd(gadget_spectra_sn5$spectral_cube[15,15,] - gadget_spectra_nonoise$spectral_cube[15,15,]))

})

# Testing that spawning conserves different properties -----------------------
test_that("Gas spawning creates even images despite different n-values",{

  ss_eagle_500  = make_simspin_file(ss_pd_eagle, write_to_file = FALSE,
                                    template = "BC03lr", sph_spawn_n = 500)

  ss_eagle_100 = ss_illustris  = make_simspin_file(ss_pd_eagle,
                                                   write_to_file = FALSE,
                                                   template = "BC03lr",
                                                   sph_spawn_n = 100)

  n500 =  build_datacube(simspin_file = ss_eagle_500,
                         telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA, fov = 25),
                         observing_strategy = observing_strategy(dist_z = 0.1, inc_deg = 45, blur = F),
                         method = "gas")

  n100 = build_datacube(simspin_file = ss_eagle_100,
                        telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = NA, fov = 25),
                        observing_strategy = observing_strategy(dist_z = 0.1, inc_deg = 45, blur = F),
                        method = "gas")


  mask100 = n100$raw_images$particle_image; mask100[mask100 >= 1] = 1; mask100[mask100!=1]=NA
  mask500 = n500$raw_images$particle_image; mask500[mask500 >= 1] = 1; mask500[mask500!=1]=NA

  expect_equal(sum(n500$raw_images$mass_image*mask500, na.rm=T), sum(n100$raw_images$mass_image*mask100, na.rm=T))
  expect_equal(sum(n500$raw_images$SFR_image), sum(n100$raw_images$SFR_image))

})

# Checking that the voronoi binning works as we expect -------------------------
test_that("Voronoi bin maps are produced and saved to the raw_images output", {

  vorbin_spectral = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "spectral", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_spectral$raw_images))
  expect_length(vorbin_spectral$raw_images, (spectra_raw_images_size + 1))

  vorbin_velocity = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "velocity", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_velocity$raw_images))
  expect_length(vorbin_velocity$raw_images, (velocity_raw_images_size + 1))

  vorbin_velocity_mf = build_datacube(simspin_file = ss_eagle,
                                      telescope = telescope(type = "SAMI"),
                                      observing_strategy = observing_strategy(dist_z = 0.03),
                                      method = "velocity", verbose = F, mass_flag = T,
                                      voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_velocity_mf$raw_images))
  expect_length(vorbin_velocity_mf$raw_images, (velocity_raw_images_size + 1))

  vorbin_gas      = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "gas", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_gas$raw_images))

  vorbin_sfgas    = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "sf gas", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_sfgas$raw_images))

})

test_that("Summed images are consistent in the binned and unbinned images - stellar", {

  vorbin_velocity = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "velocity", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  notbin_velocity = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "velocity", verbose = F,
                                   voronoi_bin = F)

  expect_equal(vorbin_velocity$raw_images$particle_image, notbin_velocity$raw_images$particle_image)
  expect_equal(vorbin_velocity$raw_images$flux_image, notbin_velocity$raw_images$flux_image)
  expect_equal(vorbin_velocity$raw_images$mass_image, notbin_velocity$raw_images$mass_image)
  expect_equal(vorbin_velocity$observed_images$flux_image, notbin_velocity$observed_images$flux_image)
  expect_false(all(vorbin_velocity$observed_images$velocity_image == notbin_velocity$observed_images$velocity_image))
  expect_false(all(vorbin_velocity$raw_images$velocity_image == notbin_velocity$raw_images$velocity_image))

})

test_that("Summed images are consistent in the binned and unbinned images - gas", {

  vorbin_gas      = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "gas", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  notbin_gas      = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "gas", verbose = F,
                                   voronoi_bin = F)

  expect_equal(vorbin_gas$raw_images$particle_image, notbin_gas$raw_images$particle_image)
  expect_equal(vorbin_gas$raw_images$mass_image, notbin_gas$raw_images$mass_image)
  expect_equal(vorbin_gas$raw_images$SFR_image, notbin_gas$raw_images$SFR_image)
  expect_false(all(vorbin_gas$observed_images$velocity_image == notbin_gas$observed_images$velocity_image))
  expect_false(all(vorbin_gas$raw_images$velocity_image == notbin_gas$raw_images$velocity_image))

})

test_that("Vorbin images made with single and multi-core methods are the same", {

  vorbin_velocity = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03),
                                   method = "velocity", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  vorbin_velocity_mc = build_datacube(simspin_file = ss_eagle,
                                      telescope = telescope(type = "SAMI"),
                                      observing_strategy = observing_strategy(dist_z = 0.03),
                                      method = "velocity", verbose = F,
                                      voronoi_bin = T, vorbin_limit = 10, cores = 2)

  expect_equal(vorbin_velocity$raw_images$voronoi_bins, vorbin_velocity_mc$raw_images$voronoi_bins)
})

test_that("Voronoi binned images can be blurred with the PSF", {

  vorbin_spectral = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03, blur = T),
                                   method = "spectral", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_spectral$raw_images))
  expect_length(vorbin_spectral$raw_images, (spectra_raw_images_size + 1))

  vorbin_velocity = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03, blur = T),
                                   method = "velocity", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_velocity$raw_images))
  expect_length(vorbin_velocity$raw_images, (velocity_raw_images_size + 1))

  vorbin_velocity_mf = build_datacube(simspin_file = ss_eagle,
                                      telescope = telescope(type = "SAMI"),
                                      observing_strategy = observing_strategy(dist_z = 0.03, blur = T),
                                      method = "velocity", verbose = F, mass_flag = T,
                                      voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_velocity_mf$raw_images))
  expect_length(vorbin_velocity_mf$raw_images, (velocity_raw_images_size + 1))

  vorbin_gas      = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03, blur = T),
                                   method = "gas", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_gas$raw_images))

  vorbin_sfgas    = build_datacube(simspin_file = ss_eagle,
                                   telescope = telescope(type = "SAMI"),
                                   observing_strategy = observing_strategy(dist_z = 0.03, blur = T),
                                   method = "sf gas", verbose = F,
                                   voronoi_bin = T, vorbin_limit = 10)

  expect_true("voronoi_bins" %in% names(vorbin_sfgas$raw_images))

})
