# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the build_datacube.R code

library(testthat)
context("Testing build_datacube function.\n")

ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
temp_loc = tempdir()
make_simspin_file(ss_eagle, output = paste(temp_loc, "spectra.fst", sep=""))

# Testing that build_datacube works
test_that("Initial run of build_datacube function #1.", {
  expect_vector(build_datacube(simspin_file = paste(temp_loc, "spectra.fst", sep=""),
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45, blur = T)),
                ptype=list(), size = 2)
})

# Testing that build_datacube will give warning if the spectra given is low res
test_that("build_datacube issues warning when spectral resolution < LSF fwhm.", {
  expect_warning(build_datacube(simspin_file = paste(temp_loc, "spectra.fst", sep=""), telescope = telescope(type="SAMI"),
                               observing_strategy = observing_strategy(z = 0.05, inc_deg = 45)))
})

unlink(paste(temp_loc, "spectra.fst", sep=""))
