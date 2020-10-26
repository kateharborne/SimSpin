# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the make_simspin_file.R code

library(testthat)
context("Testing make_simspin_file function.\n")

ss_gadget = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")
ss_hdf5 = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")

temp_loc = tempdir()

# Test that they the function runs successfully without error
test_that("Initial run of each simulation type.", {
  expect_null(make_simspin_file(ss_gadget, output = paste(temp_loc, "gadget_test_spectra.fst", sep="")))
  expect_null(make_simspin_file(ss_hdf5, output = paste(temp_loc, "hdf5_test_spectra.fst", sep="")))
  expect_null(make_simspin_file(ss_eagle, output = paste(temp_loc, "eagle_test_spectra.fst", sep="")))
})

# Test that the function fails when the file already exists
test_that("Error when output file already exists and overwrite = F",{
  expect_error(make_simspin_file(ss_gadget, output = paste(temp_loc, "gadget_test_spectra.fst", sep="")))
  expect_error(make_simspin_file(ss_hdf5, output = paste(temp_loc, "hdf5_test_spectra.fst", sep="")))
  expect_error(make_simspin_file(ss_eagle, output = paste(temp_loc, "eagle_test_spectra.fst", sep="")))
})

# Test that the function will overwrite if overwrite = T
test_that("Test that old files can be overwritten", {
  expect_null(make_simspin_file(ss_eagle, output = paste(temp_loc, "eagle_test_spectra.fst", sep=""),
                                 overwrite = T))
})

unlink(c(paste(temp_loc, "gadget_test_spectra.fst", sep=""), paste(temp_loc, "hdf5_test_spectra.fst", sep=""),
         paste(temp_loc, "eagle_test_spectra.fst", sep="")))
