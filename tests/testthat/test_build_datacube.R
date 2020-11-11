# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the build_datacube.R code

library(testthat)
context("Testing build_datacube function.\n")

ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
temp_loc = tempdir()

test_that("simspin file can be made", {
  skip_on_travis()
  expect_null(make_simspin_file(ss_eagle, output = paste(temp_loc, "spectra.fst", sep="")))
})

# Testing that build_datacube works
test_that("Initial run of build_datacube function #1.", {
  skip_on_travis()
  expect_length(build_datacube(simspin_file = paste(temp_loc, "spectra.fst", sep=""),
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), 2)
})

# Testing that build_datacube will give warning if the spectra given is low res
test_that("build_datacube issues warning when spectral resolution < LSF fwhm.", {
  skip_on_travis()
  expect_warning(build_datacube(simspin_file = paste(temp_loc, "spectra.fst", sep=""), telescope = telescope(type="SAMI"),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45)))
})

unlink(paste(temp_loc, "spectra.fst", sep=""))

# Testing that the spectral_weights functions work as expected
test_that("spectra are weighted by mass correctly", {
  weights = list(V5389 = c(7, 8, 0.5117805, 0.4882195, 39, 40, 0.6736590, 0.3263410),
                 V5464 = c(6, 7, 0.09855354,  0.90144646, 40, 41, 0.79954516, 0.20045484))
  mass = c(0.0001274539, 0.0001326889)
  spectra = matrix(data=NA, nrow=2, ncol=53689)
  Template = ProSpect::EMILES
  for (i in 1:2){
    w = weights[[i]]; m = mass[i]
    spectra[i,] = ((Template$Zspec[[w[2]]][w[6],] * w[4]*w[8]) +
                   (Template$Zspec[[w[2]]][w[5],] * w[4]*w[7]) +
                   (Template$Zspec[[w[1]]][w[6],] * w[3]*w[8]) +
                   (Template$Zspec[[w[1]]][w[5],] * w[3]*w[7])) * 1e10 * m

  }
  expect_equal(.compute_spectra(weights, mass, Template), spectra, tolerance = 1e-6)
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
