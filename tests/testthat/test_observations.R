# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the telescope, observing strategy and observation functions

library(testthat)
context("Testing observation description functions.\n")

telescope_length = 11
obs_strategy_wo_blur_length = 5
obs_strategy_w_blur_length = 7
obs_length = 33

# Testing that the telescope() function works with each of the possible default
#  "types".
test_that("Initial run of telescope() function with default types - SAMI.", {
  expect_length(telescope(type="SAMI"), telescope_length)
  })

test_that("Initial run of telescope() function with default types - MaNGA", {
  expect_length(telescope(type="MaNGA"), telescope_length)
  expect_warning(telescope(type="MaNGA", fov = 15))

  manga15 = suppressWarnings(telescope(type="MaNGA", fov = 15))
  expect_true(manga15$fov == 17)
  expect_length(manga15, telescope_length)
  remove(manga15)

  manga33 = suppressWarnings(telescope(type="MaNGA", fov = 33))
  expect_true(manga33$fov == 32)
  expect_length(manga33, telescope_length)
  remove(manga33)
})

test_that("Initial run of telescope() function with default types - MUSE", {
  expect_length(telescope(type="MUSE"), telescope_length)
  expect_length(telescope(type="MUSE", fov = 16), telescope_length)
  expect_length(telescope(type="MUSE", spatial_res = 0.2), telescope_length)
  expect_length(telescope(type="MUSE", spatial_res = 0.025), telescope_length)
  expect_warning(telescope(type="MUSE", spatial_res = 0.1))
  expect_warning(telescope(type="MUSE", fov = 61))

  # checking that fov is automatically rounded down to 60
  muse60 = suppressWarnings(telescope(type="MUSE", fov = 61))
  expect_true(muse60$fov == 60)
  expect_length(muse60, telescope_length)
  remove(muse60)

  # checking that spatial res is automatically set to 0.2
  muse02 = suppressWarnings(telescope(type="MUSE", spatial_res = 0.1))
  expect_true(muse02$spatial_res == 0.2)
  expect_length(muse02, telescope_length)
  remove(muse02)

  # checking that spatial res can be set to NFM
  museNFM = telescope(type="MUSE", spatial_res = 0.025)
  expect_true(museNFM$spatial_res == 0.025)
  expect_length(museNFM, telescope_length)
  remove(museNFM)


})

test_that("Initial run of telescope() function with default types - Hector", {
  expect_length(telescope(type="Hector"), telescope_length)
})

test_that("Initial run of telescope() function with default types - CALIFA", {
  expect_length(telescope(type="CALIFA"), telescope_length)
})

test_that("Initial run of telescope() function with default types - IFU", {
  expect_length(telescope(type="IFU"), telescope_length)
})

# Testing case sensitivity of "type" parameter.
test_that("Checking case sensitivity #1", {
  expect_length(telescope(type="SAMI"), telescope_length)
})

test_that("Checking case sensitivity #2", {
  expect_length(telescope(type="sami"), telescope_length)
})

test_that("Checking case sensitivity #3", {
  expect_length(telescope(type="Sami"), telescope_length)
})

test_that("Checking case sensitivity #4", {
  expect_length(telescope(type="SaMi"), telescope_length)
  expect_length(telescope(type="sAmI"), telescope_length)
})

# Testing that the telescope() function works with each of the possible default
#  "filters".
test_that("Initial run of telescope() function with default filters - r", {
  expect_length(telescope(type="IFU", filter = "r", wave_range = c(3000,9000)), telescope_length)
})

test_that("Initial run of telescope() function with default filters - u", {
  expect_length(telescope(type="IFU", filter = "u", wave_range = c(3000,9000)), telescope_length)
})

test_that("Initial run of telescope() function with default filters - g", {
  expect_length(telescope(type="IFU", filter = "g", wave_range = c(3000,9000)), telescope_length)
})

test_that("Initial run of telescope() function with default filters - i", {
  expect_length(telescope(type="IFU", filter = "i", wave_range = c(3000,9000)), telescope_length)
})

test_that("Initial run of telescope() function with default filters - z", {
  expect_length(telescope(type="IFU", filter = "z", wave_range = c(3000,9000)), telescope_length)
})

test_that("Error when requestedfilter is outside range of wavlengths", {
  expect_error(telescope(type="IFU", filter = "z", wave_range = c(3000,6000)))
  expect_error(telescope(type="IFU", filter = "u", wave_range = c(5000,6000)))
  expect_error(telescope(type="SAMI", filter = "z", wave_range = c(3000,6000)))
})

# Testing that an error triggers when you give it a range in the wrong format,
#  or an unsupported aperture shape.
test_that("telescope() issues error when incompatible parameters are given #1.", {
  expect_error(telescope(wave_range = seq(3700,5700)))
})

test_that("telescope() issues error when incompatible parameters are given #2.", {
  expect_error(telescope(aperture_shape = "octogon"))
})

test_that("telescope() issues error when incompatible parameters are given #3.", {
  expect_error(telescope(filter = "p"))
})

test_that("telescope() fixes the inordered wave-range.", {
  expect_equal(telescope(wave_range = c(5700,3700))$wave_range, c(3700,5700))
})

# Testing that warnings are issued if the method parameter is specified.
test_that("telescope() issues warning when method parameter is given.", {
  expect_warning(telescope(method = "spectral"))
  scope = suppressWarnings(telescope(method = "spectral"))
  expect_length(scope, telescope_length+1)
  remove(scope)

  expect_warning(telescope(type="SAMI", method = "spectral"))
  expect_warning(telescope(type="MaNGA", method = "spectral"))
  expect_warning(telescope(type="MaNGA", method = "spectral", fov = 15))
  expect_warning(telescope(type="MUSE", method = "spectral"))
  expect_warning(telescope(type="MUSE", method = "spectral", fov = 61))
  expect_warning(telescope(type="MUSE", method = "spectral", spatial_res = 0.1))
  expect_warning(telescope(type="Hector", method = "spectral"))
  expect_warning(telescope(type="CALIFA", method = "spectral"))
})

test_that("telescope() issues error when incompatible parameters are given.", {
  expect_error(telescope(method = "flat"))
})

# Testing that the observing_strategy() function works with each of the possible
#  default "types".
test_that("Initial run of observing_strategy() function - w/o blur", {
  expect_length(observing_strategy(), obs_strategy_wo_blur_length)
})

test_that("Initial run of observing_strategy() function - w/ blur", {
  expect_length(observing_strategy(blur=T), obs_strategy_w_blur_length)
})

# Testing that an error triggers if incompatible parameters are given
test_that("observing_strategy() issues error when incompatible parameters are given.", {
  expect_error(observing_strategy(blur = T, psf = "round"))
  expect_error(observing_strategy(dist_z=0))
  expect_error(observing_strategy(dist_z=-0.2))
})

# Testing that you can generate an observation with each of the telescope types
test_that("Initial run of observation() function with default types #1.", {
  expect_length(observation(telescope(type="SAMI"), observing_strategy(), method="spectral"), obs_length)
})

test_that("Initial run of observation() function with default types #2.", {
  expect_length(observation(telescope(type="MaNGA"), observing_strategy(), method="spectral"), obs_length)
})

test_that("Initial run of observation() function with default types #3.", {
  expect_length(observation(telescope(type="Hector"), observing_strategy(), method="spectral"), obs_length)
})

test_that("Initial run of observation() function with default types #4.", {
  expect_length(observation(telescope(type="CALIFA"), observing_strategy(), method="spectral"), obs_length)
})

test_that("Initial run of observation() function with default types #5 - w/o blur.", {
  expect_length(observation(telescope(type="IFU"), observing_strategy(), method="spectral"), obs_length)
})

test_that("Initial run of observation() function with default types #6 w/ blur.", {
  expect_length(observation(telescope(type="IFU"), observing_strategy(blur=T), method="spectral"), obs_length)
})

# Testing that the psf produced is symmetrical
test_that("The PSF shape produced is symmetrical - Gaussian", {
  expect_true(isSymmetric(observation(telescope(type="SAMI"), observing_strategy(blur=T), method="spectral")$psf_kernel))
  expect_true(isSymmetric(observation(telescope(type="SAMI"), observing_strategy(blur=T, psf="Moffat"), method="spectral")$psf_kernel))
})

# Testing uncovered features
test_that("Aperture shape = 'square' works", {
  expect_length(telescope(type="IFU", aperture_shape = "square"), telescope_length)
  expect_length(observation(telescope(type="IFU", aperture_shape = "square"), observing_strategy(), method="spectral"), obs_length)
})

test_that("Check even field-of-view leads to psf scaled", {
  expect_length(telescope(type="IFU", fov = 10), telescope_length)
  expect_length(observation(telescope(type="IFU", fov = 10), observing_strategy(blur=T), method="spectral")$psf_kernel, 19*19)
})

test_that("Initial run of telescope() function with non-default types - check fov leads to psf scaled", {
  expect_length(telescope(type="IFU", fov = 11.5), telescope_length)
  expect_length(observation(telescope(type="IFU", fov = 11.5), observing_strategy(blur=T), method="spectral")$psf_kernel, 23*23)

})

# Testing Distance class
test_that("Check that when input z, or Mpc, or kpc_per_arcsec to Distance, other answers are consistent", {
  set_z = Distance(z=0.3)
  set_Mpc = Distance(Mpc=1588.662)
  set_kpc = Distance(kpc_per_arcsec = 4.557426)

  z_tol = 1e-7
  Mpc_tol = 1e-3
  kpc_tol = 1e-5

  expect_equal(z(set_z), z(set_Mpc), tolerance = z_tol)
  expect_equal(z(set_z), z(set_kpc), tolerance = z_tol)
  expect_equal(Mpc(set_z), Mpc(set_Mpc), tolerance = Mpc_tol)
  expect_equal(Mpc(set_z), Mpc(set_kpc), tolerance = Mpc_tol)
  expect_equal(kpc_per_arcsec(set_z), kpc_per_arcsec(set_Mpc), tolerance = kpc_tol)
  expect_equal(kpc_per_arcsec(set_z), kpc_per_arcsec(set_kpc), tolerance = kpc_tol)
})

test_that("Distance class errors if wrong parameter specified", {
  expect_error(Distance(k = 22))
})

test_that("Distance class errors if multiple parameters specified", {
  expect_error(Distance(z = 0.2, kpc_per_arcsec = 1))
})

# Testing Distance class still performs as expected within observing_strategy() function
test_that("Check that when input z, or Mpc, or kpc_per_arcsec to observing_strategy, other answers are consistent", {
  set_z = observing_strategy(dist_z=0.3)
  set_Mpc = observing_strategy(dist_Mpc=1588.662)
  set_kpc = observing_strategy(dist_kpc_per_arcsec = 4.557426)

  z_tol = 1e-7
  Mpc_tol = 1e-3
  kpc_tol = 1e-5

  expect_equal(z(set_z$distance), z(set_Mpc$distance), tolerance = z_tol)
  expect_equal(z(set_z$distance), z(set_kpc$distance), tolerance = z_tol)
  expect_equal(Mpc(set_z$distance), Mpc(set_Mpc$distance), tolerance = Mpc_tol)
  expect_equal(Mpc(set_z$distance), Mpc(set_kpc$distance), tolerance = Mpc_tol)
  expect_equal(kpc_per_arcsec(set_z$distance), kpc_per_arcsec(set_Mpc$distance), tolerance = kpc_tol)
  expect_equal(kpc_per_arcsec(set_z$distance), kpc_per_arcsec(set_kpc$distance), tolerance = kpc_tol)
})

test_that("observing_strategy() errors if wrong parameter specified for Distance class", {
  expect_error(observing_strategy(k = 22))
})

test_that("observing_strategy errors if multiple parameters specified for Distance class", {
  expect_error(observing_strategy(dist_z = 0.2, dist_kpc_per_arcsec = 1))
})

# Testing Pointing class
test_that("Check that when input xy_deg or xy_kpc to Pointing, other answers are consistent", {
  distance = Distance(z=0.3)
  set_deg = Pointing(xy_deg = c(0.1,0.05), distance = distance)
  set_kpc = Pointing(xy_kpc = c(1640.673,820.3367), distance = distance)

  deg_tol = 1e-7
  kpc_tol = 1e-3

  expect_equal(xy_deg(set_kpc), xy_deg(set_deg), tolerance = deg_tol)
  expect_equal(xy_kpc(set_kpc), xy_kpc(set_deg), tolerance = kpc_tol)
  expect_output(show(distance))
  expect_output(show(set_deg))

})

test_that("Pointing class errors if wrong parameter specified", {
  expect_error(Pointing(k = 22))
})

test_that("Pointing class errors if multiple parameters specified", {
  expect_error(Pointing(xy_deg = c(0.1,0.05), xy_kpc = c(1640.673,820.3367)))
  expect_error(Pointing(xy_deg = c(0.1,0.05), xy_kpc = c(1640.673,820.3367), distance = Distance(z=0.3)))
})

# Testing Pointing class still performs as expected within observing_strategy() function
test_that("Check that when input xy_deg or xy_kpc to Pointing, other answers are consistent", {
  set_deg = observing_strategy(dist_z=0.3, pointing_deg = c(0.1,0.05))
  set_kpc = observing_strategy(dist_z=0.3, pointing_kpc = c(1640.673,820.3367))

  deg_tol = 1e-7
  kpc_tol = 1e-3

  expect_equal(xy_deg(set_kpc$pointing), xy_deg(set_deg$pointing), tolerance = deg_tol)
  expect_equal(xy_kpc(set_kpc$pointing), xy_kpc(set_deg$pointing), tolerance = kpc_tol)

})

test_that("Pointing class errors if wrong parameter specified", {
  expect_error(observing_strategy(dist_z=0.3, pointing_k = c(0.1)))
})

test_that("Pointing class errors if multiple parameters specified", {
  expect_error(observing_strategy(dist_z=0.3, pointing_kpc = c(0.1), pointing_deg = c(123,125)))
})

# Testing that a warning is issued if "z" input is specified -------------------
test_that("Warning is issued if old input parameter is specified", {
  expect_warning(observing_strategy(z = 0.3))
})
