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
test_that("Initial run of each simulation type - Gadget.", {
  expect_null(make_simspin_file(ss_gadget, template = "BC03hr", output = paste(temp_loc, "gadget_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "gadget_test", sep="")), 4)
  expect_true(length(readRDS(paste(temp_loc, "gadget_test", sep=""))$gas_part) == 0)
})

test_that("Initial run of each simulation type - HDF5", {
  expect_null(make_simspin_file(ss_hdf5, output = paste(temp_loc, "hdf5_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "hdf5_test", sep="")), 4)
  expect_true(length(readRDS(paste(temp_loc, "hdf5_test", sep=""))$gas_part) == 0)
})

test_that("Initial run of each simulation type - EAGLE", {
  expect_null(make_simspin_file(ss_eagle, output = paste(temp_loc, "eagle_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "eagle_test", sep="")), 4)
  expect_true(length(readRDS(paste(temp_loc, "eagle_test", sep=""))$gas_part) == 14)
})

# Test that the function fails when the file already exists
test_that("Error when output file already exists and overwrite = F - Gadget",{
  expect_error(make_simspin_file(ss_gadget, output = paste(temp_loc, "gadget_test", sep="")))
  })

test_that("Error when output file already exists and overwrite = F - HDF5",{
  expect_error(make_simspin_file(ss_hdf5, output = paste(temp_loc, "hdf5_test", sep="")))
})

test_that("Error when output file already exists and overwrite = F - EAGLE",{
  expect_error(make_simspin_file(ss_eagle, output = paste(temp_loc, "eagle_test", sep="")))
})

# Test that function can output to environment
test_that("Function works to output List to environment", {
  ss_file = make_simspin_file(ss_gadget, write_to_file = FALSE)
  expect_true("ss_file" %in% ls())
})

# Test that the function fails when the template provided is unsupported
test_that("Error when template is unsupported", {
  expect_error(make_simspin_file(ss_eagle, template = "abcd"))
})

# Test that the function will overwrite if overwrite = T
test_that("Test that old files can be overwritten", {
  expect_null(make_simspin_file(ss_gadget, template = "EMILES", output = paste(temp_loc, "gadget_test", sep=""),
                                 overwrite = T))
})

test_that("Test that function works on multiple cores", {
  expect_null(make_simspin_file(ss_eagle, template = "EMILES", output = paste(temp_loc, "gadget_test", sep=""),
                                overwrite = T, cores = 2))
})

# Test that none of the output arrays have NAs in
test_that("Values are successfully associated with variables", {
  gadget = readRDS(paste(temp_loc, "gadget_test", sep=""))
  hdf5   = readRDS(paste(temp_loc, "hdf5_test", sep=""))
  eagle  = readRDS(paste(temp_loc, "eagle_test", sep=""))

  expect_true(all(!is.na(gadget$star_part$x)))
  expect_true(all(!is.na(gadget$star_part$y)))
  expect_true(all(!is.na(gadget$star_part$z)))
  expect_true(all(!is.na(gadget$star_part$vx)))
  expect_true(all(!is.na(gadget$star_part$vy)))
  expect_true(all(!is.na(gadget$star_part$vz)))
  expect_true(all(!is.na(gadget$star_part$Mass)))
  expect_true(all(!is.na(gadget$star_part$sed_id)))
  expect_true(all(!is.na(gadget$star_part$Metallicity)))
  expect_true(all(!is.na(gadget$star_part$Age)))
  expect_true(all(!is.na(gadget$star_part$Initial_Mass)))

  expect_true(all(!is.na(hdf5$star_part$x)))
  expect_true(all(!is.na(hdf5$star_part$y)))
  expect_true(all(!is.na(hdf5$star_part$z)))
  expect_true(all(!is.na(hdf5$star_part$vx)))
  expect_true(all(!is.na(hdf5$star_part$vy)))
  expect_true(all(!is.na(hdf5$star_part$vz)))
  expect_true(all(!is.na(hdf5$star_part$Mass)))
  expect_true(all(!is.na(hdf5$star_part$sed_id)))
  expect_true(all(!is.na(hdf5$star_part$Metallicity)))
  expect_true(all(!is.na(hdf5$star_part$Age)))
  expect_true(all(!is.na(hdf5$star_part$Initial_Mass)))

  expect_true(all(!is.na(eagle$star_part$x)))
  expect_true(all(!is.na(eagle$star_part$y)))
  expect_true(all(!is.na(eagle$star_part$z)))
  expect_true(all(!is.na(eagle$star_part$vx)))
  expect_true(all(!is.na(eagle$star_part$vy)))
  expect_true(all(!is.na(eagle$star_part$vz)))
  expect_true(all(!is.na(eagle$star_part$Mass)))
  expect_true(all(!is.na(eagle$star_part$sed_id)))
  expect_true(all(!is.na(eagle$star_part$Metallicity)))
  expect_true(all(!is.na(eagle$star_part$Age)))
  expect_true(all(!is.na(eagle$star_part$Initial_Mass)))

  expect_true(all(!is.na(eagle$gas_part$x)))
  expect_true(all(!is.na(eagle$gas_part$y)))
  expect_true(all(!is.na(eagle$gas_part$z)))
  expect_true(all(!is.na(eagle$gas_part$vx)))
  expect_true(all(!is.na(eagle$gas_part$vy)))
  expect_true(all(!is.na(eagle$gas_part$vz)))
  expect_true(all(!is.na(eagle$gas_part$Mass)))
  expect_true(all(!is.na(eagle$gas_part$Z)))
  expect_true(all(!is.na(eagle$gas_part$Density)))
  expect_true(all(!is.na(eagle$gas_part$Temp)))
  expect_true(all(!is.na(eagle$gas_part$SFR)))
  expect_true(all(!is.na(eagle$gas_part$Metallicity)))
  expect_true(all(!is.na(eagle$gas_part$OEOS)))

})

unlink(c(paste(temp_loc, "gadget_test", sep=""), paste(temp_loc, "hdf5_test", sep=""),
         paste(temp_loc, "eagle_test", sep="")))
