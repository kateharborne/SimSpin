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

sa_out = 4

# Testing that sim_analysis works with each galaxy type - type="stars"----
test_that("Gadget files can be analysed", {
  expect_length(sim_analysis(ss_gadget), sa_out)
})

test_that("HDF5 files can be analysed", {
  expect_length(sim_analysis(ss_hdf5), sa_out)
})

test_that("EAGLE files can be analysed", {
  eagle_stars = sim_analysis(ss_eagle)
  expect_length(eagle_stars, sa_out)
  expect_true(!"MeanSFR" %in% names(eagle_stars$Properties))
  expect_true("MeanAge" %in% names(eagle_stars$Properties))
})

test_that("Magenticum files can be analysed", {
  expect_length(sim_analysis(ss_magneticum), sa_out)
})

# Testing that sim_analysis works with each galaxy type - type="gas"----

test_that("Gadget files in gas mode fail", {
  expect_error(sim_analysis(ss_gadget, type="gas"))
})

test_that("HDF5 files in gas mode fail", {
  expect_error(sim_analysis(ss_hdf5, type = "gas"))
})

test_that("EAGLE files can be analysed in gas mode", {
  eagle_gas = sim_analysis(ss_eagle, type = "gas")
  expect_length(eagle_gas, sa_out)
  expect_true("MeanSFR" %in% names(eagle_gas$Properties))
  expect_true("Temperature" %in% names(eagle_gas$RadialTrends_Spherical))
  expect_true("Temperature" %in% names(eagle_gas$RadialTrends_Cylindrical))
  expect_true(!"MeanAge" %in% names(eagle_gas$Properties))
})

test_that("Magenticum files can be analysed in gas mode", {
  expect_length(sim_analysis(ss_magneticum, type = "gas"), sa_out)
})
