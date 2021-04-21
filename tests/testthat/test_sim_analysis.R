# Date: 21/04/21
# Title: Testing the sim_analysis.R code

library(testthat)
context("Testing sim_analysis functions.\n")

ss_pd_hdf5  = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_pd_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
ss_pd_magneticum = system.file("extdata", "SimSpin_example_Magneticum.hdf5", package = "SimSpin")

ss_gadget   = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata", package = "SimSpin")
ss_hdf5     = make_simspin_file(ss_pd_hdf5, write_to_file = FALSE)
ss_eagle    = make_simspin_file(ss_pd_eagle, write_to_file = FALSE)
ss_magneticum = make_simspin_file(ss_pd_magneticum, write_to_file = FALSE)

# Testing that sim_analysis works with each galaxy type - type="stars"----
test_that("Gadget files can be analysed", {
  expect_length(sim_analysis(ss_gadget), 3)
})

test_that("HDF5 files can be analysed", {
  expect_length(sim_analysis(ss_hdf5), 3)
})

test_that("EAGLE files can be analysed", {
  expect_length(sim_analysis(ss_eagle), 3)
})

test_that("Magenticum files can be analysed", {
  expect_length(sim_analysis(ss_magneticum), 3)
})

# Testing that sim_analysis works with each galaxy type - type="gas"----

test_that("Gadget files in gas mode fail", {
  expect_error(sim_analysis(ss_gadget, type="gas"))
})

test_that("HDF5 files in gas mode fail", {
  expect_error(sim_analysis(ss_hdf5, type = "gas"))
})

test_that("EAGLE files can be analysed in gas mode", {
  expect_length(sim_analysis(ss_eagle, type = "gas"), 3)
})

test_that("Magenticum files can be analysed in gas mode", {
  expect_length(sim_analysis(ss_magneticum, type = "gas"), 3)
})
