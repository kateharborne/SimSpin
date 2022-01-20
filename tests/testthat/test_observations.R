# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the telescope, observing strategy and observation functions

library(testthat)
context("Testing observation description functions.\n")

# Testing that the telescope() function works with each of the possible default
#  "types".
test_that("Initial run of telescope() function with default types - SAMI.", {
  expect_length(telescope(type="SAMI"), 12)
  })

test_that("Initial run of telescope() function with default types - MaNGA", {
  expect_length(telescope(type="MaNGA"), 12)
})

test_that("Initial run of telescope() function with default types - MaNGA", {
  expect_length(telescope(type="MUSE"), 12)
})

test_that("Initial run of telescope() function with default types - Hector", {
  expect_length(telescope(type="Hector"), 12)
})

test_that("Initial run of telescope() function with default types - CALIFA", {
  expect_length(telescope(type="CALIFA"), 12)
})

test_that("Initial run of telescope() function with default types - IFU", {
  expect_length(telescope(type="IFU"), 12)
})

# Testing case sensitivity of "type" parameter.
test_that("Checking case sensitivity #1", {
  expect_length(telescope(type="SAMI"), 12)
})

test_that("Checking case sensitivity #2", {
  expect_length(telescope(type="sami"), 12)
})

test_that("Checking case sensitivity #3", {
  expect_length(telescope(type="Sami"), 12)
})

test_that("Checking case sensitivity #4", {
  expect_length(telescope(type="SaMi"), 12)
  expect_length(telescope(type="sAmI"), 12)
})

# Testing that the telescope() function works with each of the possible default
#  "filters".
test_that("Initial run of telescope() function with default filters - r", {
  expect_length(telescope(type="IFU", filter = "r"), 12)
})

test_that("Initial run of telescope() function with default filters - u", {
  expect_length(telescope(type="IFU", filter = "u"), 12)
})

test_that("Initial run of telescope() function with default filters - g", {
  expect_length(telescope(type="IFU", filter = "g"), 12)
})

test_that("Initial run of telescope() function with default filters - i", {
  expect_length(telescope(type="IFU", filter = "i"), 12)
})

test_that("Initial run of telescope() function with default filters - z", {
  expect_length(telescope(type="IFU", filter = "z"), 12)
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
  expect_error(telescope(method = "flat"))
})

test_that("telescope() issues error when incompatible parameters are given #4.", {
  expect_error(telescope(filter = "p"))
})

test_that("telescope() fixes the inordered wave-range.", {
  expect_equal(telescope(wave_range = c(5700,3700))$wave_range, c(3700,5700))
})

# Testing that the objective() function works with each of the possible
#  default "types".
test_that("Initial run of objective() function - w/o blur", {
  expect_length(objective(), 5)
})

test_that("Initial run of objective() function - w/ blur", {
  expect_length(objective(blur=T), 7)
})

# Testing that an error triggers if incompatible parameters are given
test_that("objective() issues error when incompatible parameters are given.", {
  expect_error(objective(blur = T, psf = "round"))
  expect_error(objective(distance = Distance(z=0)))
  expect_error(objective(distance = Distance(z=-0.2)))
})

# Testing that you can generate an observation with each of the telescope types
test_that("Initial run of observation() function with default types #1.", {
  expect_length(observation(telescope(type="SAMI"), objective()), 32)
})

test_that("Initial run of observation() function with default types #2.", {
  expect_length(observation(telescope(type="MaNGA"), objective()), 32)
})

test_that("Initial run of observation() function with default types #3.", {
  expect_length(observation(telescope(type="Hector"), objective()), 32)
})

test_that("Initial run of observation() function with default types #4.", {
  expect_length(observation(telescope(type="CALIFA"), objective()), 32)
})

test_that("Initial run of observation() function with default types #5 - w/o blur.", {
  expect_length(observation(telescope(type="IFU"), objective()), 32)
})

test_that("Initial run of observation() function with default types #6 w/ blur.", {
  expect_length(observation(telescope(type="IFU"), objective(blur=T)), 32)
})

# Testing that the psf produced is symmetrical
test_that("The PSF shape produced is symmetrical - Gaussian", {
  expect_true(isSymmetric(observation(telescope(type="SAMI"), objective(blur=T))$psf_kernel))
  expect_true(isSymmetric(observation(telescope(type="SAMI"), objective(blur=T, psf="Moffat"))$psf_kernel))
})

# Testing uncovered features
test_that("Aperture shape = 'square' works", {
  expect_length(telescope(type="IFU", aperture_shape = "square"), 12)
  expect_length(observation(telescope(type="IFU", aperture_shape = "square"), objective()), 32)
})

test_that("Check even field-of-view leads to psf scaled", {
  expect_length(telescope(type="IFU", fov = 10), 12)
  expect_length(observation(telescope(type="IFU", fov = 10), objective(blur=T))$psf_kernel, 19*19)
})

test_that("Initial run of telescope() function with non-default types - check fov leads to psf scaled", {
  expect_length(telescope(type="IFU", fov = 11.5), 12)
  expect_length(observation(telescope(type="IFU", fov = 11.5), objective(blur=T))$psf_kernel, 23*23)

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

# Testing Distance class still performs as expected within objective() function
test_that("Check that when input z, or Mpc, or kpc_per_arcsec to objective, other answers are consistent", {
  set_z = objective(distance = Distance(z=0.3))
  set_Mpc = objective(distance = Distance(Mpc=1588.662))
  set_kpc = objective(distance = Distance(kpc_per_arcsec = 4.557426))

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

test_that("objective() errors if wrong parameter specified for Distance class", {
  expect_error(objective(distance = Distance(k = 22)))
})

test_that("objective errors if multiple parameters specified for Distance class", {
  expect_error(objective(distance = Distance(z = 0.2, kpc_per_arcsec = 1)))
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

})

test_that("Pointing class errors if wrong parameter specified", {
  expect_error(Pointing(k = 22))
})

test_that("Pointing class errors if multiple parameters specified", {
  expect_error(Pointing(xy_deg = c(0.1,0.05), xy_kpc = c(1640.673,820.3367)))
  expect_error(Pointing(xy_deg = c(0.1,0.05), xy_kpc = c(1640.673,820.3367), distance = Distance(z=0.3)))
})

# Testing Pointing class still performs as expected within objective() function
test_that("Check that when input xy_deg or xy_kpc to Pointing, other answers are consistent", {
  distance = Distance(z=0.3)
  set_deg = objective(distance = distance, pointing = c(0.1,0.05), pointing_unit = "deg")
  set_kpc = objective(distance = distance, pointing = c(1640.673,820.3367), pointing_unit = "kpc")

  deg_tol = 1e-7
  kpc_tol = 1e-3

  expect_equal(xy_deg(set_kpc$pointing), xy_deg(set_deg$pointing), tolerance = deg_tol)
  expect_equal(xy_kpc(set_kpc$pointing), xy_kpc(set_deg$pointing), tolerance = kpc_tol)

})

test_that("Pointing class errors if wrong parameter specified", {
  expect_error(objective(distance = distance, pointing = c(0.1), pointing_unit = "k"))
})

test_that("Pointing class errors if multiple parameters specified", {
  expect_error(objective(distance = distance, pointing = c(0.1), pointing_unit = "kpc"))
})
