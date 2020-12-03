# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the build_datacube.R code

library(testthat)
context("Testing build_datacube function.\n")

ss_gadget = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata", package = "SimSpin")
ss_hdf5   = system.file("extdata", "SimSpin_example_HDF5_spectra.Rdata", package = "SimSpin")
ss_eagle  = system.file("extdata", "SimSpin_example_EAGLE_spectra.Rdata", package = "SimSpin")

# Testing that build_datacube works
test_that("Gadget files can be built.", {
  skip_on_travis()
  expect_length(build_datacube(simspin_file = ss_gadget,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 4)
})

test_that("HDF5 files can be built.", {
  skip_on_travis()
  expect_length(build_datacube(simspin_file = ss_hdf5,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 4)
})

test_that("EAGLE files can be built.", {
  skip_on_travis()
  expect_length(build_datacube(simspin_file = ss_eagle,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 4)
})

# Testing that build_datacube will give warning if the spectra given is low res
test_that("build_datacube issues warning when spectral resolution < LSF fwhm.", {
  expect_warning(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="SAMI"),
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
  simspin_data = readRDS(ss_eagle)

  galaxy_data = obs_galaxy(part_data = simspin_data$star_part, inc_rad = 0.7853982) # projecting the galaxy to given inclination
  sbin_seq = c(-7.502401, -7.002241, -6.502080, -6.001920, -5.501760, -5.001600,
               -4.501440, -4.001280, -3.501120, -3.000960, -2.500800, -2.000640,
               -1.500480, -1.000320, -0.500160,  0.000000,  0.500160,  1.000320,
               1.500480, 2.000640, 2.500800, 3.000960, 3.501120, 4.001280, 4.501440,
               5.001600, 5.501760, 6.001920, 6.502080, 7.002241, 7.502401)
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=sbin_seq, labels=F) +
    (30 * cut(galaxy_data$z_obs, breaks=sbin_seq, labels=F)) - (30) # assigning particles to positions in cube

  i = 411

  particle_IDs = which(galaxy_data$pixel_pos == i)
  galaxy_sample = galaxy_data[particle_IDs,]

  intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]), nrow = length(particle_IDs), byrow = T)
  spectra = intrinsic_spectra * (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra

  expect_true(all(intrinsic_spectra[1,] == simspin_data$spectra[[5]]))
  expect_equal((intrinsic_spectra[1,] * galaxy_sample$Initial_Mass[1] * 1e10), spectra[1,])
  expect_equal((intrinsic_spectra[2,] * galaxy_sample$Initial_Mass[2] * 1e10), spectra[2,])
})
