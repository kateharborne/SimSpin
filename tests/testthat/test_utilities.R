# Author: Kate Harborne
# Co-author: Alice Serene
# Date: 13/01/2023
# Title: Testing the utilities.R code

library(testthat)
context("Testing utilities functions.\n")

gadget_data = .read_gadget(system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin"))
hdf5_data   = .read_hdf5(system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin"))
eagle_data  = .read_hdf5(system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin"))
#illustris_data = .read_hdf5(system.file("extdata", "SimSpin_example_IllustrisTNG.hdf5", package = "SimSpin"))

test_that("Gadget file can be read", {
  expect_length(gadget_data, 3)
})

test_that("EAGLE HDF5 file can be read", {
  expect_length(eagle_data, 4)
})

test_that("The number of stellar particles is equal to the number of SSP particles", {
  expect_equal(dim(eagle_data$star_part)[1], length(eagle_data$ssp$Initial_Mass))
})

test_that("Standard HDF5 file can be read", {
  expect_length(hdf5_data, 3)
})

test_that("Length of particle data and info in header are the same", {
  expect_equal(dim(hdf5_data$star_part)[1], sum(hdf5_data$head$Npart))
})

test_that("Gas information returns NULL for galaxy without gas", {
  expect_null(gadget_data$gas_part)
})

test_that("Gas information returns df for galaxy with gas", {
  expect_true(is.data.frame(eagle_data$gas_part))
})

test_that("Galaxy centering works as expected for Gadget w/o gas", {
  cen_data = .centre_galaxy(gadget_data)
  expect_length(cen_data, length(gadget_data)) # produces a list of the same size as input
  # centers the data such that the new medians are zero.
  expect_equal(c(median(cen_data$star_part$x),median(cen_data$star_part$y), median(cen_data$star_part$z)), c(0,0,0))
  expect_equal(c(median(cen_data$star_part$vx),median(cen_data$star_part$vy), median(cen_data$star_part$vz)), c(0,0,0))
})

test_that("Galaxy centering works as expected for EAGLE w/ gas", {
  new_data = data.table::copy(eagle_data)
  cen_data = .centre_galaxy(new_data)
  expect_length(cen_data, length(eagle_data)) # produces a list of the same size as input
  # centers the data such that the new medians are zero.
  expect_equal(c(median(cen_data$star_part$x),median(cen_data$star_part$y), median(cen_data$star_part$z)), c(0,0,0))
  expect_equal(c(median(cen_data$star_part$vx),median(cen_data$star_part$vy), median(cen_data$star_part$vz)), c(0,0,0))

  expect_equal((median(eagle_data$gas_part$x) - median(eagle_data$star_part$x)), median(cen_data$gas_part$x))
  expect_equal((median(eagle_data$gas_part$y) - median(eagle_data$star_part$y)), median(cen_data$gas_part$y))
  expect_equal((median(eagle_data$gas_part$z) - median(eagle_data$star_part$z)), median(cen_data$gas_part$z))
  expect_equal((median(eagle_data$gas_part$vx) - median(eagle_data$star_part$vx)), median(cen_data$gas_part$vx))
  expect_equal((median(eagle_data$gas_part$vy) - median(eagle_data$star_part$vy)), median(cen_data$gas_part$vy))
  expect_equal((median(eagle_data$gas_part$vz) - median(eagle_data$star_part$vz)), median(cen_data$gas_part$vz))

})

test_that("Galaxy alignment is aligning the gas correctly", {
  new_data = data.table::copy(eagle_data)
  new_data = .centre_galaxy(new_data) # centre galaxy
  # measure angle between J vectors before align
  J_ang_init = as.numeric(.vector_angle(.vector_unit(angmom_galaxy(new_data$star_part)),
                                        .vector_unit(angmom_galaxy(new_data$gas_part))))

  aligned_data = .align_galaxy(new_data)
  # measure angle between J vectors after align
  J_ang_end  = as.numeric(.vector_angle(.vector_unit(angmom_galaxy(aligned_data$star_part)),
                                        .vector_unit(angmom_galaxy(aligned_data$gas_part))))

  expect_equal(J_ang_init, J_ang_end)

})

test_that("SED is returned as a list and data.table when 1 or multiple cores are used", {
  ages = sample(seq(3e7, 1.4e10, length.out = 100), 100)/1e9
  metallicity = sample(seq(1e-4, 4e-2, length.out = 100), 100)
  temp = SimSpin::BC03lr
  expect_type(.spectral_weights(Metallicity = metallicity, Age = ages, Template = temp, cores = 1), "list")
  expect_type(.spectral_weights(Metallicity = metallicity, Age = ages, Template = temp, cores = 2), "list")
  expect_true(data.table::is.data.table(.spectral_weights(Metallicity = metallicity, Age = ages, Template = temp, cores = 1)))
  expect_true(data.table::is.data.table(.spectral_weights(Metallicity = metallicity, Age = ages, Template = temp, cores = 2)))
})

test_that("SED is returned as a list when only one age and metallicity are given", {
  temp = SimSpin::BC03lr
  expect_type(.spectral_weights(Metallicity = 1e-4, Age = 5, Template = temp, cores = 1), "list")
  expect_true(data.table::is.data.table(.spectral_weights(Metallicity = 1e-4, Age = 5, Template = temp, cores = 1)))
})

test_that("Filters and templates can be loaded successfully", {
  BC03lr = SimSpin::BC03lr
  BC03hr = SimSpin::BC03hr
  EMILES = SimSpin::EMILES

  expect_type(BC03lr, "list")
  expect_type(BC03hr, "list")
  expect_type(EMILES, "list")

  g_filter = SimSpin::filt_g_SDSS
  r_filter = SimSpin::filt_r_SDSS
  u_filter = SimSpin::filt_u_SDSS
  i_filter = SimSpin::filt_i_SDSS
  z_filter = SimSpin::filt_z_SDSS

  expect_type(g_filter, "list")
  expect_type(r_filter, "list")
  expect_type(u_filter, "list")
  expect_type(i_filter, "list")
  expect_type(z_filter, "list")
})

test_that("interp_quick fails if given more than a scalar", {
 expect_error(.interp_quick(x = c(1,2,3), params = c(seq(1,10), log = T)))
})

test_that("kernel weights sum to 1", {
  expect_equal(sum(.generate_uniform_sphere(100, kernel = "WC2")$weight), 1)
  expect_equal(sum(.generate_uniform_sphere(100, kernel = "WC6")$weight), 1)
  expect_equal(sum(.generate_uniform_sphere(100, kernel = "M4")$weight), 1)
})
