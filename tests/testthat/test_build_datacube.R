# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the build_datacube.R code

library(testthat)
context("Testing build_datacube function.\n")

ss_gadget = system.file("extdata", "SimSpin_example_Gadget_spectra.fst", package = "SimSpin")
temp_loc = tempdir()


# Testing that build_datacube will run
test_that("Initial run of build_datacube function with defaults.", {
  expect_vector(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="SAMI"),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45)), ptype=list(), size = 2)
  expect_vector(build_datacube(simspin_file = ss_gadget, telescope = telescope(type="MaNGA"),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)), ptype=list(), size = 2)
})
