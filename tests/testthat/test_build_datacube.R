# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the build_datacube.R code

library(testthat)
context("Testing build_datacube function.\n")

ss_pd_hdf5  = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_pd_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")

ss_gadget   = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata", package = "SimSpin")
ss_hdf5     = make_simspin_file(ss_pd_hdf5, write_to_file = FALSE)
ss_eagle    = make_simspin_file(ss_pd_eagle, write_to_file = FALSE)

temp_loc = tempdir()

# Testing that build_datacube works in spectral mode
test_that("Gadget files can be built.", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", method = "spectral", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               verbose = T), 5)
})

test_that("HDF5 files can be built.", {
  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", method = "spectral", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 5)
})

test_that("EAGLE files can be built.", {
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", method = "spectral", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 5)
})

test_that("EAGLE files can be built in parallel.", {
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", method = "spectral", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               cores = 2), 5)
})

# Testing that build_datacube works in velocity mode
test_that("Gadget files can be built.", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               verbose = T), 7)
})

test_that("HDF5 files can be built.", {
  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 7)
})

test_that("EAGLE files can be built.", {
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 7)
})

test_that("EAGLE files can be built in parallel.", {
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               cores = 2), 7)
})

# Testing that build_datacube works to write to FITS file
test_that("Data cubes can be written to file", {
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               write_fits = T), 5)

  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", method="velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               write_fits = T, output_location = paste0(temp_loc, "velocity_cube.FITS")), 7)

  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               write_fits = T, output_location = paste0(temp_loc, "cube.FITS")), 5)
})

unlink(c(paste(temp_loc, "cube.FITS", sep=""), paste(temp_loc, "velocity_cube.FITS", sep=""), paste(stringr::str_remove(ss_gadget, ".Rdata"), "_inc45deg_seeing2fwhm.FITS", sep="")))

# Testing that build_datacube will give warning if the spectra given is low res
test_that("build_datacube issues warning when spectral resolution < LSF fwhm.", {
  expect_warning(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="IFU", lsf_fwhm = 0.9),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45)))
})

# Testing that the velocity shift functions work as expected
test_that("velocity shift for wavelengths work correctly", {

  wavelength = ProSpect::EMILES$Wave
  velocity_los = c(27.04932, 40.94573)
  wave = matrix(data = rep(wavelength, length(velocity_los)), nrow = length(velocity_los), byrow=T)
  wave_shift = ((velocity_los / .speed_of_light) * wave) + wave # using doppler formula to compute the shift in wavelengths cause by LOS velocity

  wave_shift_comp = matrix(data=NA, nrow=2, ncol=53689)

  for (i in 1:2){
    wave_shift_comp[i,] = ((velocity_los[i] / .speed_of_light) * wave[i,]) + wave[i,]
  }

  expect_equal(wave_shift, wave_shift_comp)

})

# Test that the spectra pulled for each particle are correct
test_that("Repeated spectra are included in intrinsic spectra", {
  simspin_data = readRDS(ss_gadget)
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

  intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]), nrow = length(particle_IDs), byrow = T)
  spectra = intrinsic_spectra * (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra

  expect_true(all(intrinsic_spectra[3,] == simspin_data$spectra[[2]]))
  expect_equal((intrinsic_spectra[1,] * galaxy_sample$Initial_Mass[1] * 1e10), spectra[1,])
  expect_equal((intrinsic_spectra[2,] * galaxy_sample$Initial_Mass[2] * 1e10), spectra[2,])
})

# Test that the output of blurred and un-blurred format is the same
test_that("Format of blurring output is the same as unblurred output", {

  unblurred = build_datacube(simspin_file = ss_hdf5,
                          telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                          observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = F))

  blurred = build_datacube(simspin_file = ss_hdf5,
                           telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                           observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T))

  expect_equal(typeof(blurred$spectral_cube), "double")
  expect_equal(typeof(blurred$observation), "list")
  expect_equal(typeof(blurred$velocity_image), "double")
  expect_equal(typeof(blurred$dispersion_image), "double")
  expect_equal(typeof(unblurred$spectral_cube), "double")
  expect_equal(typeof(unblurred$observation), "list")
  expect_equal(typeof(unblurred$velocity_image), "double")
  expect_equal(typeof(unblurred$dispersion_image), "double")

})

# Testing the twisting and inclination changes make sense
test_that("Twisting and inclination work as expected", {
  SAMI = telescope(type="SAMI")
  strategy = SimSpin::observing_strategy(z = 0.05, inc_deg = 90, twist_deg = 0) # viewing from the front
  observation = observation(SAMI, strategy)
  twisted_data = twist_galaxy(ss_eagle$star_part, twist_rad = observation$twist_rad)
  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  front = galaxy_data$vy

  strategy = SimSpin::observing_strategy(z = 0.05, inc_deg = 90, twist_deg = 180) # viewing from the back
  observation = SimSpin::observation(SAMI, strategy)
  twisted_data = twist_galaxy(ss_eagle$star_part, twist_rad = observation$twist_rad)
  galaxy_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  back = galaxy_data$vy

  expect_equal(front, -1*(back)) # velocities should be equal but opposite signs
})

# Test the observations get dimmer with distance added
test_that("Observations get dimmer with increasing redshift", {
  cube_near = build_datacube(simspin_file = ss_gadget,
                             telescope = telescope(type="IFU", method = "spectral", signal_to_noise = NA),
                             observing_strategy = observing_strategy(z = 0.01, inc_deg = 45, blur = F),
                             verbose = F)

  cube_far  = build_datacube(simspin_file = ss_gadget,
                             telescope = telescope(type="IFU", method = "spectral", signal_to_noise = NA),
                             observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = F),
                             verbose = F)

  expect_true(sum(cube_near$flux_image, na.rm=T) > sum(cube_far$flux_image, na.rm=T))

})

# Test that velocity cubes can be built in mass mode
test_that("EAGLE cubes can be built with mass weighting rather than luminosity", {
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               mass_flag = T), 7)
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                               mass_flag = T, cores=2), 7)

  expect_true(all(build_datacube(simspin_file = ss_eagle,
                                 telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                                 mass_flag = T)$flux_image ==
                  build_datacube(simspin_file = ss_eagle,
                                 telescope = telescope(type="IFU", method = "velocity", lsf_fwhm = 3.6, signal_to_noise = 3),
                                 observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T),
                                 mass_flag = T, cores=2)$flux_image, na.rm=T))

})

test_that("Mass/flux images are different for the same observing conditions.", {
  expect_true(all((build_datacube(simspin_file = ss_eagle,
                                  telescope = telescope(type="IFU", method = "velocity",
                                                        lsf_fwhm = 3.6, signal_to_noise = NA),
                                  observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = F),
                                  mass_flag = T)$flux_image) !=
                    (build_datacube(simspin_file = ss_eagle,
                                    telescope = telescope(type="IFU", method = "velocity",
                                                          lsf_fwhm = 3.6, signal_to_noise = NA),
                                    observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = F),
                                    mass_flag = F)$flux_image), na.rm =T))

})
