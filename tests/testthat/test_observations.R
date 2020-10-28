# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the telescope, observing strategy and observation functions

library(testthat)
context("Testing observation description functions.\n")

# Testing that the telescope() function works with each of the possible default
#  "types".
test_that("Initial run of telescope() function with default types.", {
  expect_vector(telescope(type="SAMI"), ptype=list(), size = 8)
  expect_vector(telescope(type="MaNGA"), ptype=list(), size = 8)
  expect_vector(telescope(type="Hector"), ptype=list(), size = 8)
  expect_vector(telescope(type="CALIFA"), ptype=list(), size = 8)
  expect_vector(telescope(type="IFU"), ptype=list(), size = 8)
  })

# Testing case sensitivity of "type" parameter.
test_that("Initial run of telescope() function with default types.", {
  expect_vector(telescope(type="SAMI"), ptype=list(), size = 8)
  expect_vector(telescope(type="sami"), ptype=list(), size = 8)
  expect_vector(telescope(type="Sami"), ptype=list(), size = 8)
  expect_vector(telescope(type="SaMi"), ptype=list(), size = 8)
  expect_vector(telescope(type="sAmI"), ptype=list(), size = 8)
})

# Testing that an error triggers when you give it a range in the wrong format,
#  or an unsupported aperture shape.
test_that("telescope() issues error when incompatible parameters are given.", {
  expect_error(telescope(wave_range = seq(3700,5700)))
  expect_error(telescope(aperture_shape = "octogon"))
})

# Testing that the observing_strategy() function works with each of the possible
#  default "types".
test_that("Initial run of observing_strategy() function", {
  expect_vector(observing_strategy(), ptype=list(), size = 3)
  expect_vector(observing_strategy(blur=T), ptype = list(), size=5)
})

# Testing that an error triggers if an incompatible psf shape is given.
test_that("observing_strategy() issues error when incompatible parameters are given.", {
  expect_error(observing_strategy(blur = T, psf = "round"))
})

# Testing that you can generate an observation with each of the telescope types
test_that("Initial run of telescope() function with default types.", {
  expect_vector(observation(telescope(type="SAMI"), observing_strategy = observing_strategy()), ptype=list(), size = 18)
  expect_vector(observation(telescope(type="MaNGA"), observing_strategy = observing_strategy()), ptype=list(), size = 18)
  expect_vector(observation(telescope(type="Hector"), observing_strategy = observing_strategy()), ptype=list(), size = 18)
  expect_vector(observation(telescope(type="CALIFA"), observing_strategy = observing_strategy()), ptype=list(), size = 18)
  expect_vector(observation(telescope(type="IFU"), observing_strategy = observing_strategy()), ptype=list(), size = 18)
  expect_vector(observation(telescope(type="IFU"), observing_strategy = observing_strategy(blur=T)), ptype=list(), size = 18)
})
