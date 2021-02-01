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

# Testing that the observing_strategy() function works with each of the possible
#  default "types".
test_that("Initial run of observing_strategy() function - w/o blur", {
  expect_length(observing_strategy(), 5)
})

test_that("Initial run of observing_strategy() function - w/ blur", {
  expect_length(observing_strategy(blur=T), 7)
})

# Testing that an error triggers if an incompatible psf shape is given.
test_that("observing_strategy() issues error when incompatible parameters are given.", {
  expect_error(observing_strategy(blur = T, psf = "round"))
})

# Testing that you can generate an observation with each of the telescope types
test_that("Initial run of observation() function with default types #1.", {
  expect_length(observation(telescope(type="SAMI"), observing_strategy = observing_strategy()), 30)
})

test_that("Initial run of observation() function with default types #2.", {
  expect_length(observation(telescope(type="MaNGA"), observing_strategy = observing_strategy()), 30)
})

test_that("Initial run of observation() function with default types #3.", {
  expect_length(observation(telescope(type="Hector"), observing_strategy = observing_strategy()), 30)
})

test_that("Initial run of observation() function with default types #4.", {
  expect_length(observation(telescope(type="CALIFA"), observing_strategy = observing_strategy()), 30)
})

test_that("Initial run of observation() function with default types #5 - w/o blur.", {
  expect_length(observation(telescope(type="IFU"), observing_strategy = observing_strategy()), 30)
})

test_that("Initial run of observation() function with default types #6 w/ blur.", {
  expect_length(observation(telescope(type="IFU"), observing_strategy = observing_strategy(blur=T)), 30)
})

# Testing that the psf produced is symmetrical
test_that("The PSF shape produced is symmetrical - Gaussian", {
  expect_true(isSymmetric(observation(telescope(type="SAMI"), observing_strategy = observing_strategy(blur=T))$psf_kernel))
  expect_true(isSymmetric(observation(telescope(type="SAMI"), observing_strategy = observing_strategy(blur=T, psf="Moffat"))$psf_kernel))
})
