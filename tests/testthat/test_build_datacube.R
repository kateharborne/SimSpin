# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the build_datacube.R code

library(testthat)
context("Testing build_datacube function.\n")

ss_pd_gadget   = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")
ss_pd_hdf5  = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_pd_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
ss_pd_magneticum = system.file("extdata", "SimSpin_example_Magneticum.hdf5", package = "SimSpin")

ss_gadget_old = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata", package = "SimSpin")
ss_gadget   = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
ss_hdf5     = make_simspin_file(ss_pd_hdf5, write_to_file = FALSE)
ss_eagle    = make_simspin_file(ss_pd_eagle, write_to_file = FALSE)
ss_magneticum = make_simspin_file(ss_pd_magneticum, write_to_file = FALSE)

temp_loc = tempdir()

# Testing that build_datacube works in spectral mode ----
built_cube_size = 4
spectra_raw_images_size = 4
spectra_observed_images_size = NULL
velocity_raw_images_size = 6
velocity_observed_images_size = 3

test_that("Gadget files can be built - spectral mode", {
  expect_warning(build_datacube(simspin_file = ss_gadget_old,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                verbose = F))
  gadget_spectra = build_datacube(simspin_file = ss_gadget,
                                  telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                  observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                  verbose = T)
  expect_length(gadget_spectra, built_cube_size)
  expect_length(gadget_spectra$raw_images, spectra_raw_images_size)
  expect_null(gadget_spectra$observed_images)
})

test_that("HDF5 files can be built - spectral mode", {
  hdf5_spectra = build_datacube(simspin_file = ss_hdf5,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T))
  expect_length(hdf5_spectra, built_cube_size)
  expect_length(hdf5_spectra$raw_images, spectra_raw_images_size)
  expect_null(hdf5_spectra$observed_images)
})

test_that("EAGLE files can be built - spectral mode", {
  eagle_spectra = build_datacube(simspin_file = ss_eagle,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T))
  expect_length(eagle_spectra, built_cube_size)
  expect_length(eagle_spectra$raw_images, spectra_raw_images_size)
  expect_null(eagle_spectra$observed_images)
})

test_that("EAGLE files can be built in parallel - spectral mode", {
  eagle_parallel_spectra = build_datacube(simspin_file = ss_eagle,
                                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                          observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                          cores = 2)
  expect_length(eagle_parallel_spectra, built_cube_size)
  expect_length(eagle_parallel_spectra$raw_images, spectra_raw_images_size)
  expect_null(eagle_parallel_spectra$observed_images)

})

test_that("Magneticum files can be built - spectral mode", {
  magneticum_spectra = build_datacube(simspin_file = ss_magneticum,
                                      telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                      observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T))
  expect_length(magneticum_spectra, built_cube_size)
  expect_length(magneticum_spectra$raw_images, spectra_raw_images_size)
  expect_null(magneticum_spectra$observed_images)
})

test_that("Magneticum files can be built in parallel - spectral mode", {
  magneticum_parallel_spectra = build_datacube(simspin_file = ss_magneticum,
                                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                               cores = 2)
  expect_length(magneticum_parallel_spectra, built_cube_size)
  expect_length(magneticum_parallel_spectra$raw_images, spectra_raw_images_size)
  expect_null(magneticum_parallel_spectra$observed_images)
})

# Testing that build_datacube works in velocity mode ----
test_that("Gadget files can be built - velocity mode.", {
  gadget_velocity = build_datacube(simspin_file = ss_gadget,
                                   telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                   observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                   method = "velocity",
                                   verbose = T)

  expect_length(gadget_velocity, built_cube_size)
  expect_length(gadget_velocity$raw_images, velocity_raw_images_size)
  expect_length(gadget_velocity$observed_images, velocity_observed_images_size)

})

test_that("HDF5 files can be built - velocity mode.", {
  hdf5_velocity = build_datacube(simspin_file = ss_hdf5,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                 method = "velocity")
  expect_length(hdf5_velocity, built_cube_size)
  expect_length(hdf5_velocity$raw_images, velocity_raw_images_size)
  expect_length(hdf5_velocity$observed_images, velocity_observed_images_size)

})

test_that("EAGLE files can be built - velocity mode.", {
  eagle_velocity = build_datacube(simspin_file = ss_eagle,
                                  telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                  observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                  method = "velocity")
  expect_length(eagle_velocity, built_cube_size)
  expect_length(eagle_velocity$raw_images, velocity_raw_images_size)
  expect_length(eagle_velocity$observed_images, velocity_observed_images_size)

})

test_that("EAGLE files can be built in parallel - velocity mode.", {
  eagle_parallel_velocity = build_datacube(simspin_file = ss_eagle,
                                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                           observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                           method = "velocity",
                                           cores = 2)
  expect_length(eagle_parallel_velocity, built_cube_size)
  expect_length(eagle_parallel_velocity$raw_images, velocity_raw_images_size)
  expect_length(eagle_parallel_velocity$observed_images, velocity_observed_images_size)
})

test_that("Magneticum files can be built - velocity mode.", {
  magneticum_velocity = build_datacube(simspin_file = ss_magneticum,
                                       telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                       observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                       method = "velocity")
  expect_length(magneticum_velocity, built_cube_size)
  expect_length(magneticum_velocity$raw_images, velocity_raw_images_size)
  expect_length(magneticum_velocity$observed_images, velocity_observed_images_size)

})

test_that("Magneticum files can be built in parallel - velocity mode.", {
  magneticum_parallel_velocity = build_datacube(simspin_file = ss_magneticum,
                                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                                method = "velocity",
                                                cores = 2)
  expect_length(magneticum_parallel_velocity, built_cube_size)
  expect_length(magneticum_parallel_velocity$raw_images, velocity_raw_images_size)
  expect_length(magneticum_parallel_velocity$observed_images, velocity_observed_images_size)

})

# Testing that build_datacube works in gas mode ----

test_that("EAGLE files can be built - gas mode.", {
  eagle_gas = build_datacube(simspin_file = ss_eagle,
                             telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                             observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                             method = "gas")
  expect_length(eagle_gas, built_cube_size)
  expect_length(eagle_gas$raw_images, velocity_raw_images_size)
  expect_length(eagle_gas$observed_images, velocity_observed_images_size)
})

test_that("EAGLE files can be built in parallel - gas mode.", {
  eagle_parallel_gas = build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method = "gas",
                               cores = 2)
  expect_length(eagle_parallel_gas, built_cube_size)
  expect_length(eagle_parallel_gas$raw_images, velocity_raw_images_size)
  expect_length(eagle_parallel_gas$observed_images, velocity_observed_images_size)

})

test_that("Magneticum files can be built - gas mode.", {
  magneticum_gas = build_datacube(simspin_file = ss_magneticum,
                                  telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                  observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 0, twist_deg = 90, blur = T),
                                  method = "gas")
  expect_length(magneticum_gas, built_cube_size)
  expect_length(magneticum_gas$raw_images, velocity_raw_images_size)
  expect_length(magneticum_gas$observed_images, velocity_observed_images_size)
})

test_that("Magneticum files can be built in parallel - gas mode.", {
  magneticum_parallel_gas = build_datacube(simspin_file = ss_magneticum,
                                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                           observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                           method = "gas",
                                           cores = 2)
  expect_length(magneticum_parallel_gas, built_cube_size)
  expect_length(magneticum_parallel_gas$raw_images, velocity_raw_images_size)
  expect_length(magneticum_parallel_gas$observed_images, velocity_observed_images_size)
})

# Testing that build_datacube works in sf gas mode ----

test_that("EAGLE files can be built - sf gas mode.", {
  eagle_sf_gas = build_datacube(simspin_file = ss_eagle,
                                telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                method = "sf gas")
  expect_length(eagle_sf_gas, built_cube_size)
  expect_length(eagle_sf_gas$raw_images, velocity_raw_images_size)
  expect_length(eagle_sf_gas$observed_images, velocity_observed_images_size)
})

test_that("EAGLE files can be built in parallel - sf gas mode.", {
  eagle_parallel_sf_gas = build_datacube(simspin_file = ss_eagle,
                                         telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                         observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                         method = "sf gas",
                                         cores = 2)
  expect_length(eagle_parallel_sf_gas, built_cube_size)
  expect_length(eagle_parallel_sf_gas$raw_images, velocity_raw_images_size)
  expect_length(eagle_parallel_sf_gas$observed_images, velocity_observed_images_size)
})

test_that("Magneticum files can be built - sf gas mode.", {
  magneticum_sf_gas = build_datacube(simspin_file = ss_magneticum,
                                     telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                     observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                     method = "sf gas")

  expect_length(magneticum_sf_gas, built_cube_size)
  expect_length(magneticum_sf_gas$raw_images, velocity_raw_images_size)
  expect_length(magneticum_sf_gas$observed_images, velocity_observed_images_size)
})


test_that("Magneticum files can be built in parallel - sf gas mode.", {
  magneticum_parallel_sf_gas = build_datacube(simspin_file = ss_magneticum,
                                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                              observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                              method = "sf gas",
                                              cores = 2)
  expect_length(magneticum_parallel_sf_gas, built_cube_size)
  expect_length(magneticum_parallel_sf_gas$raw_images, velocity_raw_images_size)
  expect_length(magneticum_parallel_sf_gas$observed_images, velocity_observed_images_size)
})

# Testing that build_datacube errors when invalid method given ----
test_that("Error occurs when invalid method given.", {
  expect_error(build_datacube(simspin_file = ss_eagle,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                              method = "magpie"))
})

# Testing the mass flag functionalilty
test_that("Data cubes can be generated using mass rather than luminosity weighting", {
  eagle_mass = build_datacube(simspin_file = ss_eagle,
                              telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                              observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                              method = "velocity",
                              write_fits = F, mass_flag = T)
  expect_length(eagle_mass, built_cube_size)
  expect_true("velocity_cube" %in% names(eagle_mass))
  expect_false("spectral_cube" %in% names(eagle_mass))
  expect_true("mass_image" %in% names(eagle_mass$raw_images))
  expect_true("mass_image" %in% names(eagle_mass$observed_images))
  expect_false("flux_image" %in% names(eagle_mass$raw_images))
  expect_false("flux_image" %in% names(eagle_mass$observed_images))
  })


# Testing that build_datacube works to write to FITS file
test_that("Data cubes can be written to a single files", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               write_fits = T), built_cube_size)

  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS"))

  spectral_fits = Rfits::Rfits_read_all("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS", header = T)
  expect_true(length(spectral_fits) == 6)
  expect_true(spectral_fits[[2]]$keyvalues$CTYPE3 == "WAVE")
  expect_true(all(dim(spectral_fits[[2]]$imDat) == c(spectral_fits[[2]]$keyvalues$NAXIS1, spectral_fits[[2]]$keyvalues$NAXIS2, spectral_fits[[2]]$keyvalues$NAXIS3)))
  expect_true(spectral_fits[[4]]$keyvalues$EXTNAME == "RAW_VEL")

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method="velocity",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))) == 8)
  expect_true(Rfits::Rfits_read(paste0(temp_loc, "/ss_gadget.FITS"))[[4]]$keyvalues$EXTNAME == "OBS_VEL")

  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method="gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_eagle.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_eagle.FITS")))

  expect_length(build_datacube(simspin_file = ss_magneticum,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method="sf gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_magneticum.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum.FITS")))

  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               write_fits = T, output_location = paste0(temp_loc, "/ss_hdf5.FITS"),
                               split_save=F), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5.FITS")))
  expect_true(length(Rfits::Rfits_read(paste0(temp_loc, "/ss_hdf5.FITS"))) == 6)

})

unlink(c("GalaxyID_unknown_inc45deg_seeing2fwhm.FITS",
         paste0(temp_loc, "/ss_gadget.FITS"),
         paste0(temp_loc, "/ss_eagle.FITS"),
         paste0(temp_loc, "/ss_magenticum.FITS"),
         paste0(temp_loc, "/ss_hdf5.FITS")))

test_that("Data cubes can be written to multiple files", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               write_fits = T, split_save=T), built_cube_size)

  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_flux_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_velocity_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_dispersion_image.FITS"))
  expect_true(file.exists("GalaxyID_unknown_inc45deg_seeing2fwhm_particle_image.FITS"))

  spectral_fits = Rfits::Rfits_read("GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS")
  expect_true(length(spectral_fits) == 2)
  expect_true(spectral_fits[[2]]$keyvalues$CTYPE3 == "WAVE")
  expect_true(all(dim(spectral_fits[[2]]$imDat) == c(spectral_fits[[2]]$keyvalues$NAXIS1, spectral_fits[[2]]$keyvalues$NAXIS2, spectral_fits[[2]]$keyvalues$NAXIS3)))

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method="velocity",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_gadget.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_velocity_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_flux_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_age_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_gadget_particle_image.FITS")))

  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method="gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_eagle.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_gas_velocity_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_mass_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_OH_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_eagle_SFR_image.FITS")))

  expect_length(build_datacube(simspin_file = ss_magneticum,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method="sf gas",
                               write_fits = T, output_location = paste0(temp_loc, "/ss_magneticum.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_gas_velocity_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_mass_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_OH_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_metallicity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_magneticum_SFR_image.FITS")))

  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               write_fits = T, output_location = paste0(temp_loc, "/ss_hdf5.FITS"),
                               split_save=T), built_cube_size)

  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_spectral_cube.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_flux_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_velocity_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_dispersion_image.FITS")))
  expect_true(file.exists(paste0(temp_loc, "/ss_hdf5_particle_image.FITS")))

})

unlink(c("GalaxyID_unknown_inc45deg_seeing2fwhm_spectral_cube.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_flux_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_velocity_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_dispersion_image.FITS",
         "GalaxyID_unknown_inc45deg_seeing2fwhm_particle_image.FITS",

         paste0(temp_loc, "/ss_gadget_velocity_cube.FITS"),
         paste0(temp_loc, "/ss_gadget_flux_image.FITS"),
         paste0(temp_loc, "/ss_gadget_velocity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_gadget_age_image.FITS"),
         paste0(temp_loc, "/ss_gadget_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_particle_image.FITS"),

         paste0(temp_loc, "/ss_eagle_gas_velocity_cube.FITS"),
         paste0(temp_loc, "/ss_eagle_mass_image.FITS"),
         paste0(temp_loc, "/ss_eagle_velocity_image.FITS"),
         paste0(temp_loc, "/ss_eagle_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_eagle_OH_image.FITS"),
         paste0(temp_loc, "/ss_eagle_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_eagle_SFR_image.FITS"),

         paste0(temp_loc, "/ss_magenticum_gas_velocity_cube.FITS"),
         paste0(temp_loc, "/ss_magenticum_mass_image.FITS"),
         paste0(temp_loc, "/ss_magenticum_velocity_image.FITS"),
         paste0(temp_loc, "/ss_magenticum_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_magenticum_OH_image.FITS"),
         paste0(temp_loc, "/ss_magenticum_metallicity_image.FITS"),
         paste0(temp_loc, "/ss_magenticum_SFR_image.FITS"),

         paste0(temp_loc, "/ss_hdf5_spectral_cube.FITS"),
         paste0(temp_loc, "/ss_hdf5_flux_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_velocity_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_hdf5_particle_image.FITS")
         ))

test_that("Mask can be included in FITS files correctly", {
  cube = build_datacube(simspin_file = ss_gadget,
                        telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                        observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
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
  expect_true(length(fits_w_mask) == 7)
  expect_true(fits_w_mask[[7]]$keyvalues$EXTNAME == "MASK")

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
         paste0(temp_loc, "/ss_gadget_flux_image.FITS"),
         paste0(temp_loc, "/ss_gadget_velocity_image.FITS"),
         paste0(temp_loc, "/ss_gadget_dispersion_image.FITS"),
         paste0(temp_loc, "/ss_gadget_particle_image.FITS"),
         paste0(temp_loc, "/ss_gadget_mask.FITS")))

# Testing that build_datacube will give warning if the spectra given is low res
test_that("build_datacube issues warning when spectral resolution < LSF fwhm.", {
  expect_warning(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="IFU", lsf_fwhm = 0.9),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45)))
})

# Testing that the velocity shift functions work as expected -------------------
test_that("velocity shift for wavelengths work correctly", {

  wavelength = SimSpin::EMILES$Wave
  velocity_los = c(27.04932, 40.94573)
  wave = matrix(data = rep(wavelength, length(velocity_los)), nrow = length(velocity_los), byrow=T)
  wave_shift = ((velocity_los / .speed_of_light) * wave) + wave # using doppler formula to compute the shift in wavelengths cause by LOS velocity

  wave_shift_comp = matrix(data=NA, nrow=2, ncol=53689)

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

  intrinsic_spectra = simspin_data$spectra[ , galaxy_sample$sed_id, with=FALSE]
  spectra = intrinsic_spectra * (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra

  expect_true(all(intrinsic_spectra[,c(3),] == simspin_data$spectra[["V2"]]))
  expect_equal((intrinsic_spectra[,c(1),] * galaxy_sample$Initial_Mass[1] * 1e10), spectra[,c(1),], tolerance = 0.001)
  expect_equal((intrinsic_spectra[,c(2),] * galaxy_sample$Initial_Mass[2] * 1e10), spectra[,c(2),], tolerance = 0.001)
})

# Test that the output of blurred and un-blurred format is the same ------------
test_that("Format of blurring output is the same as unblurred output", {

  unblurred = build_datacube(simspin_file = ss_hdf5,
                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                          observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = F))

  blurred = build_datacube(simspin_file = ss_hdf5,
                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                           observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T))

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
  strategy = SimSpin::observing_strategy(dist_z = 0.05, inc_deg = 90, twist_deg = 0) # viewing from the front
  observation = SimSpin::observation(SAMI, strategy, method = "spectral")
  twisted_data = twist_galaxy(ss_eagle$star_part, twist_rad = observation$twist_rad)
  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  front = galaxy_data$vy

  strategy = SimSpin::observing_strategy(dist_z = 0.05, inc_deg = 90, twist_deg = 180) # viewing from the back
  observation = SimSpin::observation(SAMI, strategy, method = "spectral")
  twisted_data = twist_galaxy(ss_eagle$star_part, twist_rad = observation$twist_rad)
  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  back = galaxy_data$vy

  expect_equal(front, -1*(back)) # velocities should be equal but opposite signs
})

# Test the observations get dimmer with distance added -------------------------
test_that("Observations get dimmer with increasing redshift", {
  cube_near = build_datacube(simspin_file = ss_gadget,
                             telescope = telescope(type="IFU", signal_to_noise = NA),
                             observing_strategy = observing_strategy(dist_z = 0.01, inc_deg = 45, blur = F),
                             verbose = F)

  cube_far  = build_datacube(simspin_file = ss_gadget,
                             telescope = telescope(type="IFU", signal_to_noise = NA),
                             observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = F),
                             verbose = F)

  expect_true(sum(cube_near$raw_images$flux_image, na.rm=T) > sum(cube_far$raw_images$flux_image, na.rm=T))

})

# Test that velocity cubes can be built in mass mode ---------------------------
test_that("EAGLE cubes can be built with mass weighting rather than luminosity", {
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method = "velocity",
                               mass_flag = T), built_cube_size)
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                               method = "velocity",
                               mass_flag = T, cores=2), built_cube_size)

  expect_true(all(build_datacube(simspin_file = ss_eagle,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                 method = "velocity",
                                 mass_flag = T)$flux_image ==
                  build_datacube(simspin_file = ss_eagle,
                                 telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = T),
                                 method = "velocity",
                                 mass_flag = T, cores=2)$flux_image, na.rm=T))

})

test_that("Mass/flux images are different for the same observing conditions.", {
  expect_true(all((build_datacube(simspin_file = ss_eagle,
                                  telescope = telescope(type="IFU",
                                                        lsf_fwhm = 3.6, signal_to_noise = NA),
                                  observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = F),
                                  method = "velocity",
                                  mass_flag = T)$flux_image) !=
                    (build_datacube(simspin_file = ss_eagle,
                                    telescope = telescope(type="IFU",
                                                          lsf_fwhm = 3.6, signal_to_noise = NA),
                                    observing_strategy = observing_strategy(dist_z = 0.05, inc_deg = 45, blur = F),
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

  centre_index = which(centred_gal$observed_images$mass_image == max(centred_gal$observed_images$mass_image, na.rm=T))

  shifted_up_gal = build_datacube(simspin_file = ss_gadget,
                                  telescope = telescope(type = "SAMI"),
                                  observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(0,5)),
                                  method="velocity",
                                  mass_flag = T)

  up_index = which(shifted_up_gal$observed_images$mass_image == max(shifted_up_gal$observed_images$mass_image, na.rm=T))

  shifted_up_gal_deg = build_datacube(simspin_file = ss_gadget,
                                  telescope = telescope(type = "SAMI"),
                                  observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_deg = c(0,0.00139)),
                                  method="velocity",
                                  mass_flag = T)

  up_index_deg = which(shifted_up_gal_deg$observed_images$mass_image == max(shifted_up_gal_deg$observed_images$mass_image, na.rm=T))

  shifted_down_gal = build_datacube(simspin_file = ss_gadget,
                                    telescope = telescope(type = "SAMI"),
                                    observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(0,-5)),
                                    method="velocity",
                                    mass_flag = T)

  down_index = which(shifted_down_gal$observed_images$mass_image == max(shifted_down_gal$observed_images$mass_image, na.rm=T))

  shifted_left_gal = build_datacube(simspin_file = ss_gadget,
                                    telescope = telescope(type = "SAMI"),
                                    observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(-5,0)),
                                    method="velocity",
                                    mass_flag = T)

  left_index =which(shifted_left_gal$observed_images$mass_image == max(shifted_left_gal$observed_images$mass_image, na.rm=T))

  shifted_right_gal = build_datacube(simspin_file = ss_gadget,
                                     telescope = telescope(type = "SAMI"),
                                     observing_strategy = observing_strategy(dist_kpc_per_arcsec = 1, inc_deg = 90, pointing_kpc = c(5,0)),
                                     method="velocity",
                                     mass_flag = T)

  right_index = which(shifted_right_gal$observed_images$mass_image == max(shifted_right_gal$observed_images$mass_image, na.rm=T))

  expect_true((centre_index - up_index) == -300) # check all centres have shifted by the expected number of pixels
  expect_true((centre_index - down_index) == 300)
  expect_true((centre_index - right_index) == -10)
  expect_true((centre_index - left_index) == 10)
  expect_true(up_index == up_index_deg) # check that pointing is the same if spefified in degrees instead
})
