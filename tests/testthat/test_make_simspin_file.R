# Author: Kate Harborne
# Co-author: Alice Serene
# Date: 13/01/2023
# Title: Testing the make_simspin_file.R code

library(testthat)
context("Testing make_simspin_file function.\n")

ss_gadget     = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")
ss_hdf5       = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_eagle      = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
ss_magneticum = system.file("extdata", "SimSpin_example_Magneticum.hdf5", package = "SimSpin")
ss_horizon    = system.file("extdata", "SimSpin_example_HorizonAGN.hdf5", package = "SimSpin")
ss_illustris  = system.file("extdata", "SimSpin_example_IllustrisTNG.hdf5", package = "SimSpin")

ss_file_length = 5
ss_file_header_length = 15

temp_loc = tempdir()

# Test that the function runs successfully without error
test_that("Initial run of each simulation type - Gadget.", {
  expect_null(make_simspin_file(ss_gadget, template = "BC03hr", output = paste(temp_loc, "/gadget_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "/gadget_test", sep="")), ss_file_length)
  expect_true(length(readRDS(paste(temp_loc, "/gadget_test", sep=""))$gas_part) == 0)
})

test_that("Initial run of each simulation type - HDF5", {
  expect_null(make_simspin_file(ss_hdf5, output = paste(temp_loc, "/hdf5_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "/hdf5_test", sep="")), ss_file_length)
  expect_true(length(readRDS(paste(temp_loc, "/hdf5_test", sep=""))$gas_part) == 0)
})

test_that("Initial run of each simulation type - EAGLE", {
  expect_null(make_simspin_file(ss_eagle, output = paste(temp_loc, "/eagle_test", sep=""), centre = c(18318,61583,38667), half_mass = 1483809589))
  expect_length(readRDS(paste(temp_loc, "/eagle_test", sep="")), ss_file_length)
  expect_true(length(readRDS(paste(temp_loc, "/eagle_test", sep=""))$gas_part) == 16)
})

test_that("Initial run of each simulation type - Magneticum", {
  expect_null(make_simspin_file(ss_magneticum, output = paste(temp_loc, "/magneticum_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "/magneticum_test", sep="")), ss_file_length)
  expect_true(length(readRDS(paste(temp_loc, "/magneticum_test", sep=""))$gas_part) == 16)
})

test_that("Initial run of each simulation type - HorizonAGN", {
  expect_null(make_simspin_file(ss_horizon, output = paste(temp_loc, "/horizon_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "/horizon_test", sep="")), ss_file_length)
  expect_true(length(readRDS(paste(temp_loc, "/horizon_test", sep=""))$gas_part) == 16)
})

test_that("Initial run of each simulation type - IllustrisTNG", {
  expect_null(make_simspin_file(ss_illustris, output = paste(temp_loc, "/illustris_test", sep="")))
  expect_length(readRDS(paste(temp_loc, "/illustris_test", sep="")), ss_file_length)
  expect_true(length(readRDS(paste(temp_loc, "/illustris_test", sep=""))$gas_part) == 16)
})

# Test that the function fails when the file already exists
test_that("Error when output file already exists and overwrite = F - Gadget",{
  expect_error(make_simspin_file(ss_gadget, output = paste(temp_loc, "/gadget_test", sep="")))
  })

test_that("Error when output file already exists and overwrite = F - HDF5",{
  expect_error(make_simspin_file(ss_hdf5, output = paste(temp_loc, "/hdf5_test", sep="")))
})

test_that("Error when output file already exists and overwrite = F - EAGLE",{
  expect_error(make_simspin_file(ss_eagle, output = paste(temp_loc, "/eagle_test", sep="")))
})

test_that("Error when output file already exists and overwrite = F - Magneticum",{
  expect_error(make_simspin_file(ss_magneticum, output = paste(temp_loc, "/magneticum_test", sep="")))
})

test_that("Error when output file already exists and overwrite = F - HorizonAGN",{
  expect_error(make_simspin_file(ss_horizon, output = paste(temp_loc, "/horizon_test", sep="")))
})

test_that("Error when output file already exists and overwrite = F - IllustrisTNG",{
  expect_error(make_simspin_file(ss_illustris, output = paste(temp_loc, "/illustris_test", sep="")))
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
  expect_null(make_simspin_file(ss_gadget, template = "EMILES", output = paste(temp_loc, "/gadget_test", sep=""),
                                 overwrite = T))
})

# Test that the function works on multiple cores
test_that("Test that function works on multiple cores", {
  expect_null(make_simspin_file(ss_eagle, template = "EMILES", output = paste(temp_loc, "/eagle_test", sep=""),
                                overwrite = T, cores = 2))
})

# Test that none of the output arrays have NAs in
test_that("Values are successfully associated with variables", {
  gadget    = readRDS(paste(temp_loc, "/gadget_test", sep=""))
  hdf5      = readRDS(paste(temp_loc, "/hdf5_test", sep=""))
  eagle     = readRDS(paste(temp_loc, "/eagle_test", sep=""))
  illustris = readRDS(paste(temp_loc, "/illustris_test", sep=""))

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

  expect_true(all(!is.na(illustris$star_part$x)))
  expect_true(all(!is.na(illustris$star_part$y)))
  expect_true(all(!is.na(illustris$star_part$z)))
  expect_true(all(!is.na(illustris$star_part$vx)))
  expect_true(all(!is.na(illustris$star_part$vy)))
  expect_true(all(!is.na(illustris$star_part$vz)))
  expect_true(all(!is.na(illustris$star_part$Mass)))
  expect_true(all(!is.na(illustris$star_part$sed_id)))
  expect_true(all(!is.na(illustris$star_part$Metallicity)))
  expect_true(all(!is.na(illustris$star_part$Age)))
  expect_true(all(!is.na(illustris$star_part$Initial_Mass)))

  expect_true(all(!is.na(illustris$gas_part$x)))
  expect_true(all(!is.na(illustris$gas_part$y)))
  expect_true(all(!is.na(illustris$gas_part$z)))
  expect_true(all(!is.na(illustris$gas_part$vx)))
  expect_true(all(!is.na(illustris$gas_part$vy)))
  expect_true(all(!is.na(illustris$gas_part$vz)))
  expect_true(all(!is.na(illustris$gas_part$Mass)))
  expect_true(all(!is.na(illustris$gas_part$SFR)))
  expect_true(all(!is.na(illustris$gas_part$Density)))
  expect_true(all(!is.na(illustris$gas_part$Temperature)))
  expect_true(all(!is.na(illustris$gas_part$Metallicity)))
  # SmoothingLength, Carbon/Hydrogen/Oxygen?
})

# Testing the sph_spawn functionality ----
test_that("Test that sph_spawn functionality works", {
  expect_null(make_simspin_file(ss_eagle, template = "EMILES", output = paste(temp_loc, "/eagle_test", sep=""),
                                overwrite = T, cores = 1, sph_spawn_n = 10))
  expect_null(make_simspin_file(ss_eagle, template = "EMILES", output = paste(temp_loc, "/eagle_test", sep=""),
                                overwrite = T, cores = 2, sph_spawn_n = 10))
  expect_null(make_simspin_file(ss_magneticum, template = "BC03lr", output = paste(temp_loc, "/magneticum_test", sep=""),
                                overwrite = T, cores = 1, sph_spawn_n = 10))
  expect_null(make_simspin_file(ss_magneticum, template = "BC03lr", output = paste(temp_loc, "/magneticum_test", sep=""),
                                overwrite = T, cores = 2, sph_spawn_n = 10))
  expect_null(make_simspin_file(ss_horizon, template = "BC03hr", output = paste(temp_loc, "/horizon_test", sep=""),
                                overwrite = T, cores = 1, sph_spawn_n = 10))
  expect_null(make_simspin_file(ss_horizon, template = "BC03hr", output = paste(temp_loc, "/horizon_test", sep=""),
                                overwrite = T, cores = 2, sph_spawn_n = 10))
  expect_null(make_simspin_file(ss_illustris, template = "BC03lr", output = paste(temp_loc, "/illustris_test", sep=""),
                               overwrite = T, cores = 1, sph_spawn_n = 10))
  expect_null(make_simspin_file(ss_illustris, template = "BC03lr", output = paste(temp_loc, "/illustris_test", sep=""),
                               overwrite = T, cores = 2, sph_spawn_n = 10))

  test_horizon = readRDS(paste(temp_loc, "/horizon_test", sep=""))
  expect_true(data.table::is.data.table(test_horizon$gas_part))
  remove(test_horizon)
})

test_that("Test that sph_spawn errors when N is not an integer", {
  expect_error(make_simspin_file(ss_eagle, template = "EMILES", output = paste(temp_loc, "/eagle_test", sep=""),
                                overwrite = T, cores = 1, sph_spawn_n = 10.2))
})

test_that("Test that sph_spawn functionality works on multiple cores - EAGLE", {
  gas_data_c1 = make_simspin_file(ss_eagle, template = "EMILES", write_to_file = FALSE, cores = 1, sph_spawn_n = 10)
  gas_data_c2 = make_simspin_file(ss_eagle, template = "EMILES", write_to_file = FALSE, cores = 2, sph_spawn_n = 10)
  expect_equal(gas_data_c1$gas_part$ID, gas_data_c2$gas_part$ID)
  expect_length(gas_data_c1$gas_part$ID, 1000) # sph_spawn_n = 10, original file contains 100 gas particles
  expect_length(gas_data_c2$gas_part$ID, 1000)
})

test_that("Test that sph_spawn functionality works on multiple cores - Magneticum", {
  gas_data_c1 = make_simspin_file(ss_magneticum, template = "EMILES", write_to_file = FALSE, cores = 1, sph_spawn_n = 10)
  gas_data_c2 = make_simspin_file(ss_magneticum, template = "EMILES", write_to_file = FALSE, cores = 2, sph_spawn_n = 10)
  expect_equal(gas_data_c1$gas_part$ID, gas_data_c2$gas_part$ID)
  expect_length(gas_data_c1$gas_part$ID, 1000) # sph_spawn_n = 10, original file contains 100 gas particles
  expect_length(gas_data_c2$gas_part$ID, 1000)
})

test_that("Test that sph_spawn functionality works on multiple cores - HorizonAGN", {
  gas_data_c1 = make_simspin_file(ss_horizon, template = "BC03", write_to_file = FALSE, cores = 1, sph_spawn_n = 10)
  gas_data_c2 = make_simspin_file(ss_horizon, template = "BC03", write_to_file = FALSE, cores = 2, sph_spawn_n = 10)
  expect_equal(gas_data_c1$gas_part$ID, gas_data_c2$gas_part$ID)
  expect_length(gas_data_c1$gas_part$ID, 1000) # sph_spawn_n = 10, original file contains 100 gas particles
  expect_length(gas_data_c2$gas_part$ID, 1000)
})

test_that("Test that sph_spawn functionality works on multiple cores - IllustrisTNG", {
  gas_data_c1 = make_simspin_file(ss_illustris, template = "BC03", write_to_file = FALSE, cores = 1, sph_spawn_n = 10)
  gas_data_c2 = make_simspin_file(ss_illustris, template = "BC03", write_to_file = FALSE, cores = 2, sph_spawn_n = 10)
  expect_equal(gas_data_c1$gas_part$ID, gas_data_c2$gas_part$ID)
  expect_length(gas_data_c1$gas_part$ID, 1000) # sph_spawn_n = 10, original file contains 100 gas particles
  expect_length(gas_data_c2$gas_part$ID, 1000)
})

# Test that the added header information works as expected ---------------------
test_that("Testing that the header data works as expected", {
  gadget     = readRDS(paste(temp_loc, "/gadget_test", sep=""))
  hdf5       = readRDS(paste(temp_loc, "/hdf5_test", sep=""))
  eagle      = readRDS(paste(temp_loc, "/eagle_test", sep=""))
  magneticum = readRDS(paste(temp_loc, "/magneticum_test", sep=""))
  horizon    = readRDS(paste(temp_loc, "/horizon_test", sep=""))
  illustris  = readRDS(paste(temp_loc, "/illustris_test", sep=""))

  expect_equal(gadget$header$Type, "nbody")
  expect_equal(hdf5$header$Type, "nbody")
  expect_equal(eagle$header$Type, "EAGLE")
  expect_equal(magneticum$header$Type, "Magneticum")
  expect_equal(horizon$header$Type, "Horizon-AGN")
  expect_equal(illustris$header$Type, "Illustris-TNG")

  expect_equal(gadget$header$Template, "EMILES")
  expect_equal(hdf5$header$Template, "BC03lr")
  expect_equal(eagle$header$Template, "EMILES")
  expect_equal(magneticum$header$Template, "BC03lr")
  expect_equal(horizon$header$Template, "BC03hr")
  expect_equal(illustris$header$Template, "BC03lr")                             # correct?

  expect_equal(gadget$header$Template_LSF, 2.51)
  expect_equal(hdf5$header$Template_LSF, 3)
  expect_equal(eagle$header$Template_LSF, 2.51)
  expect_equal(magneticum$header$Template_LSF, 3)
  expect_equal(horizon$header$Template_LSF, 3)
  expect_equal(illustris$header$Template_LSF, 3)                                # 3?

  expect_equal(gadget$header$Template_waveres, eagle$header$Template_waveres)
  expect_equal(hdf5$header$Template_waveres, 1)
  expect_equal(eagle$header$Template_waveres, 0.9)
  expect_equal(magneticum$header$Template_waveres, 1)
  expect_equal(horizon$header$Template_waveres, 1)
  expect_equal(illustris$header$Template_waveres, 1)                            # 1?

  remove(gadget, hdf5, eagle, magneticum, horizon, illustris)
 })

# Illustris-specific tests -----------------------------------------------------
# CGS conversion factors that are inherently 0 should be converted to 1
test_that("CGS conversion values of 0 get converted to 1",{
  illustris = readRDS(paste(temp_loc, "/illustris_test", sep=""))
  expect_true(all(illustris$star_part$Metallicity!=0))
})

# Stellar wind particles should be removed from the dataset
# (There are 4 wind particles in the testing set, these should be removed.)
test_that("Negative stellar wind particles are removed",{
  illustris = readRDS(paste(temp_loc, "/illustris_test", sep=""))
  expect_length(illustris$star_part$Metallicity, 1996)
})

# Gas temperature should fall within a reasonable range
test_that("Temperature does not go outside a reasonable range",{
  illustris = readRDS(paste(temp_loc, "/illustris_test", sep=""))
  expect_true(min(illustris$gas_part$Temperature)>8000)
  expect_true(max(illustris$gas_part$Temperature)<50000000)
})

unlink(c(paste(temp_loc, "/gadget_test", sep=""), paste(temp_loc, "/hdf5_test", sep=""),
         paste(temp_loc, "/eagle_test", sep=""), paste(temp_loc, "/magneticum_test", sep=""),
         paste(temp_loc, "/horizon_test", sep=""), paste(temp_loc, "/illustris_test", sep="")))

# Testing that the centre parameter works as expected ------------
test_that("Objects are centered correctly based on the specified central coordinates", {
  expect_error(make_simspin_file(filename = ss_horizon, template = "EMILES", write_to_file = F, centre = c(0,0,0)))

  centre_right = make_simspin_file(filename = ss_horizon, template = "BC03lr", write_to_file = F, centre = c(12328,72177.97,32388.04))
  centre_left  = make_simspin_file(filename = ss_horizon, template = "BC03lr", write_to_file = F, centre = c(12332,72177.97,32388.04))

  # Checking values exist in fields for LEFT
  expect_true(all(!is.na(centre_left$star_part$x)))
  expect_true(all(!is.na(centre_left$star_part$y)))
  expect_true(all(!is.na(centre_left$star_part$z)))
  expect_true(all(!is.na(centre_left$star_part$vx)))
  expect_true(all(!is.na(centre_left$star_part$vy)))
  expect_true(all(!is.na(centre_left$star_part$vz)))
  expect_true(all(!is.na(centre_left$star_part$Mass)))
  expect_true(all(!is.na(centre_left$star_part$sed_id)))
  expect_true(all(!is.na(centre_left$star_part$Metallicity)))
  expect_true(all(!is.na(centre_left$star_part$Age)))
  expect_true(all(!is.na(centre_left$star_part$Initial_Mass)))

  expect_true(all(!is.na(centre_left$gas_part$x)))
  expect_true(all(!is.na(centre_left$gas_part$y)))
  expect_true(all(!is.na(centre_left$gas_part$z)))
  expect_true(all(!is.na(centre_left$gas_part$vx)))
  expect_true(all(!is.na(centre_left$gas_part$vy)))
  expect_true(all(!is.na(centre_left$gas_part$vz)))
  expect_true(all(!is.na(centre_left$gas_part$Mass)))
  expect_true(all(!is.na(centre_left$gas_part$SFR)))
  expect_true(all(!is.na(centre_left$gas_part$Density)))
  expect_true(all(!is.na(centre_left$gas_part$Temperature)))
  expect_true(all(!is.na(centre_left$gas_part$CellSize)))
  expect_true(all(!is.na(centre_left$gas_part$Metallicity)))
  expect_true(all(!is.na(centre_left$gas_part$Carbon)))
  expect_true(all(!is.na(centre_left$gas_part$Oxygen)))
  expect_true(all(!is.na(centre_left$gas_part$Hydrogen)))

  # Checking values exist in fields for RIGHT
  expect_true(all(!is.na(centre_right$star_part$x)))
  expect_true(all(!is.na(centre_right$star_part$y)))
  expect_true(all(!is.na(centre_right$star_part$z)))
  expect_true(all(!is.na(centre_right$star_part$vx)))
  expect_true(all(!is.na(centre_right$star_part$vy)))
  expect_true(all(!is.na(centre_right$star_part$vz)))
  expect_true(all(!is.na(centre_right$star_part$Mass)))
  expect_true(all(!is.na(centre_right$star_part$sed_id)))
  expect_true(all(!is.na(centre_right$star_part$Metallicity)))
  expect_true(all(!is.na(centre_right$star_part$Age)))
  expect_true(all(!is.na(centre_right$star_part$Initial_Mass)))

  expect_true(all(!is.na(centre_right$gas_part$x)))
  expect_true(all(!is.na(centre_right$gas_part$y)))
  expect_true(all(!is.na(centre_right$gas_part$z)))
  expect_true(all(!is.na(centre_right$gas_part$vx)))
  expect_true(all(!is.na(centre_right$gas_part$vy)))
  expect_true(all(!is.na(centre_right$gas_part$vz)))
  expect_true(all(!is.na(centre_right$gas_part$Mass)))
  expect_true(all(!is.na(centre_right$gas_part$SFR)))
  expect_true(all(!is.na(centre_right$gas_part$Density)))
  expect_true(all(!is.na(centre_right$gas_part$Temperature)))
  expect_true(all(!is.na(centre_right$gas_part$CellSize)))
  expect_true(all(!is.na(centre_right$gas_part$Metallicity)))
  expect_true(all(!is.na(centre_right$gas_part$Carbon)))
  expect_true(all(!is.na(centre_right$gas_part$Oxygen)))
  expect_true(all(!is.na(centre_right$gas_part$Hydrogen)))

  galaxy_data = tryCatch(expr = {.read_gadget(ss_horizon)},
                         error = function(e){.read_hdf5(ss_horizon, cores=1)})
  galaxy_data_right = .centre_galaxy(galaxy_data, centre=c(12328,72177.97,32388.04)) # centering the galaxy based on stellar particles
  galaxy_data_left = .centre_galaxy(galaxy_data, centre=c(12332,72177.97,32388.04))

  # Need to run without align! This causes variations in the end output.
  expect_true(median(galaxy_data_left$star_part$x) < median(galaxy_data_right$star_part$x))

})

# Testing that the half_mass parameter works as expected ------------
test_that("Half-mass parameter can be specified", {

  eagle_default = make_simspin_file(filename = ss_eagle, half_mass = NA, write_to_file = F)

  eagle_hdf5 = .read_hdf5(ss_eagle, cores =1)
  total_mass = sum(eagle_hdf5$star_part$Mass)
  half_mass = (total_mass/2)

  eagle_mass_specified = make_simspin_file(filename = ss_eagle, half_mass = half_mass, write_to_file = F)
  expect_true(all(eagle_default$star_part$x == eagle_mass_specified$star_part$x))
  expect_true(all(eagle_default$star_part$y == eagle_mass_specified$star_part$y))
  expect_true(all(eagle_default$star_part$z == eagle_mass_specified$star_part$z))

})

test_that("Half-mass parameter errors when specified outside of reasonable ranges", {

  eagle_default = make_simspin_file(filename = ss_eagle, half_mass = NA, write_to_file = F)
  eagle_default = .read_hdf5(ss_eagle, cores =1)
  total_mass = sum(eagle_default$star_part$Mass)

  expect_error(make_simspin_file(filename = ss_eagle, half_mass = (total_mass+1), write_to_file = F))
  expect_error(make_simspin_file(filename = ss_eagle, half_mass = min(eagle_default$star_part$Mass), write_to_file = F))
})

# Warning if asking for half-mass with nbody model
test_that("Warning if asking for half-mass with nbody model", {

  expect_warning(make_simspin_file(ss_gadget, half_mass = 350, write_to_file = F))

})

# Test that we can still mock observe an IllustrisTNG file that has attributes with the original output names
test_that("An IllustrisTNG file with traditionally named attributes (i.e. aexp-scaling-exponent vs. a_scaling) will process successfully.",{
  file.copy(from = ss_illustris, to = paste0(temp_loc, "/SimSpin_example_IllustrisTNG_copy.hdf5"), overwrite = T)
  modified_illustris = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_IllustrisTNG_copy.hdf5"), mode = "r+") # read in the horizonagn file and rename the RubLabel
  modified_illustris[["PartType0/TestDataset"]] <- numeric(length=100) # rename an the attribute name
  hdf5r::h5attr(modified_illustris[["PartType0/TestDataset"]], "aexp-scale-exponent") <- 1
  hdf5r::h5attr(modified_illustris[["PartType0/TestDataset"]], "h-scale-exponent") <- 1
  hdf5r::h5attr(modified_illustris[["PartType0/TestDataset"]], "CGSConversionFactor") <- 20
  hdf5r::h5close(modified_illustris) #close

  expect_length(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_IllustrisTNG_copy.hdf5"), write_to_file = F), ss_file_length)

  unlink(paste0(temp_loc, "/SimSpin_example_IllustrisTNG_copy.hdf5"))
})

test_that("A file with wrongly named attributes will error.",{
  # Testing also that the code errors if the file does not contain the correct attributes.
  file.copy(from = ss_illustris, to = paste0(temp_loc, "/SimSpin_example_IllustrisTNG_copy.hdf5"), overwrite = T)
  modified_illustris = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_IllustrisTNG_copy.hdf5"), mode = "r+") # read in the horizonagn file and rename the RubLabel
  modified_illustris[["PartType0/MisnamedDataset"]] <- numeric(length=100) # rename an the attribute name
  hdf5r::h5attr(modified_illustris[["PartType0/MisnamedDataset"]], "a_misnamed_attribute") <- 7620
  hdf5r::h5attr(modified_illustris[["PartType0/MisnamedDataset"]], "h_scaling") <- 0
  hdf5r::h5attr(modified_illustris[["PartType0/MisnamedDataset"]], "to_cgs") <- 0
  hdf5r::h5close(modified_illustris) #close

  expect_error(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_IllustrisTNG_copy.hdf5"), write_to_file = F))

  unlink(paste0(temp_loc, "/SimSpin_example_IllustrisTNG.hdf5"))

})

test_that("An EAGLE file without a SmoothingLength field will error.",{
  # Testing also that the code errors if the file does not contain the correct attributes.
  file.copy(from = ss_eagle, to = paste0(temp_loc, "/SimSpin_example_EAGLE_copy.hdf5"), overwrite = T)
  modified_eagle = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_EAGLE_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel
  modified_eagle.part0 <- modified_eagle[["PartType0"]]
  modified_eagle.part0$link_delete("SmoothingLength")
  hdf5r::h5close(modified_eagle) #close

  expect_error(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_EAGLE_copy.hdf5"), write_to_file = F))

  unlink(paste0(temp_loc, "/SimSpin_example_EAGLE_copy.hdf5"))
})

test_that("A Magneticum file without a Temperature field will error.",{
  # Testing also that the code errors if the file does not contain the correct attributes.
  file.copy(from = ss_magneticum, to = paste0(temp_loc, "/SimSpin_example_Magneticum_copy.hdf5"), overwrite = T)
  modified_mag = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_Magneticum_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel
  modified_mag.part0 <- modified_mag[["PartType0"]]
  modified_mag.part0$link_delete("Temperature")
  hdf5r::h5close(modified_mag) #close

  expect_error(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_Magneticum_copy.hdf5"), write_to_file = F))

  unlink(paste0(temp_loc, "/SimSpin_example_Magneticum_copy.hdf5"))
})

test_that("A HorizonAGN file without a Temperature field will error.",{
  # Testing also that the code errors if the file does not contain the correct attributes.
  file.copy(from = ss_horizon, to = paste0(temp_loc, "/SimSpin_example_HorizonAGN_copy.hdf5"), overwrite = T)
  modified_hoz = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_HorizonAGN_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel
  modified_hoz.part0 <- modified_hoz[["PartType0"]]
  modified_hoz.part0$link_delete("Temperature")
  hdf5r::h5close(modified_hoz) #close

  expect_error(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_HorizonAGN_copy.hdf5"), write_to_file = F))

  unlink(paste0(temp_loc, "/SimSpin_example_HorizonAGN_copy.hdf5"))
})

test_that("A IllustrisTNG file without an ElectronAbundance field will error.",{
  # Testing also that the code errors if the file does not contain the correct attributes.
  file.copy(from = ss_illustris, to = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), overwrite = T)
  modified_ill = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel
  modified_ill.part0 <- modified_ill[["PartType0"]]
  modified_ill.part0$link_delete("ElectronAbundance")
  hdf5r::h5close(modified_ill) #close

  expect_error(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), write_to_file = F))

  unlink(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"))
})

test_that("A hydro sim will error without neccesary header field.",{
  # Testing also that the code errors if the file does not contain the correct attributes.
  file.copy(from = ss_illustris, to = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), overwrite = T)
  modified_tng = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel
  modified_tng.head <- modified_tng[["Header"]]
  modified_tng.head$attr_delete("Redshift")
  hdf5r::h5close(modified_tng) #close

  expect_error(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), write_to_file = F))

  unlink(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"))
})

test_that("An HDF5 file for a hydro sim will give warning without neccesary RunLabel header field.",{
  # Testing also that the code errors if the file does not contain the correct attributes.
  file.copy(from = ss_illustris, to = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), overwrite = T)
  modified_tng = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel
  modified_tng.head <- modified_tng[["Header"]]
  modified_tng.head$attr_delete("RunLabel")
  hdf5r::h5close(modified_tng) #close

  expect_error(make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), write_to_file = F))

  unlink(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"))
})

test_that("A galaxy with a star with age of 0 will NOT raise an error.", {
  # Testing that we can build a SimSpin file if we have a particle with Age = 0
  built_cube_size = 4

  file.copy(from = ss_illustris, to = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), overwrite = T)
  modified_age = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel

  sft = hdf5r::readDataSet(modified_age[["PartType4/GFM_StellarFormationTime"]])
  aexp = hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "a_scaling")
  h = hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "h_scaling")
  cgs  = hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "to_cgs")
  sft[1] = 1

  modified_age[["PartType4/"]]$link_delete("GFM_StellarFormationTime")

  modified_age[["PartType4/GFM_StellarFormationTime"]] = sft
  hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "a_scaling") = aexp
  hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "h_scaling") = h
  hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "to_cgs") = cgs

  hdf5r::h5close(modified_age)

  ss_file = make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), write_to_file = F)

  expect_length(ss_file, ss_file_length)
  expect_length(build_datacube(simspin_file = ss_file,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method = "velocity"), built_cube_size)

  unlink(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"))
})


test_that("A galaxy with a star with a very small age will NOT raise an error.", {
  # Testing that we can build a SimSpin file if we have a particle with Age = 0
  built_cube_size = 4

  file.copy(from = ss_illustris, to = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), overwrite = T)
  modified_age = hdf5r::h5file(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), mode = "r+") # read in the eagle file and rename the RubLabel

  sft = hdf5r::readDataSet(modified_age[["PartType4/GFM_StellarFormationTime"]])
  aexp = hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "a_scaling")
  h = hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "h_scaling")
  cgs  = hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "to_cgs")
  sft[1] = 0.99999999999

  modified_age[["PartType4/"]]$link_delete("GFM_StellarFormationTime")

  modified_age[["PartType4/GFM_StellarFormationTime"]] = sft
  hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "a_scaling") = aexp
  hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "h_scaling") = h
  hdf5r::h5attr(modified_age[["PartType4/GFM_StellarFormationTime"]], "to_cgs") = cgs

  hdf5r::h5close(modified_age)

  ss_file = make_simspin_file(filename = paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"), write_to_file = F)

  expect_length(ss_file, ss_file_length)
  expect_length(build_datacube(simspin_file = ss_file,
                               telescope = telescope(type="IFU", lsf_fwhm = 3.6, signal_to_noise = 3),
                               observing_strategy = observing_strategy(dist_z = 0.03, inc_deg = 45, blur = T),
                               method = "velocity"), built_cube_size)

  unlink(paste0(temp_loc, "/SimSpin_example_illustris_copy.hdf5"))
})

# Testing that the added header parameters are included in the new headers
test_that("Expect header information in each of the built simspin files.", {
  gadget_ss_file     = make_simspin_file(ss_gadget, template = "BC03hr", write_to_file = F)
  hdf5_ss_file       = make_simspin_file(ss_hdf5, template = "BC03lr", write_to_file = F)
  eagle_ss_file      = make_simspin_file(ss_eagle, template = "BC03hr", write_to_file = F, centre = c(18316,61583,38667))
  magneticum_ss_file = make_simspin_file(ss_magneticum, template = "EMILES", write_to_file = F, half_mass = 1007671144)
  horizon_ss_file    = make_simspin_file(ss_horizon, template = "BC03hr", write_to_file = F, sph_spawn_n = 3)
  illustris_ss_file  = make_simspin_file(ss_illustris, template = "BC03lr", write_to_file = F)

  expect_equal(gadget_ss_file$header$Centre, NA)
  expect_equal(gadget_ss_file$header$HalfMass, NA)
  expect_equal(gadget_ss_file$header$SmoothingN, 1)
  expect_equal(length(names(gadget_ss_file$header)), ss_file_header_length)

  expect_equal(hdf5_ss_file$header$Centre, NA)
  expect_equal(hdf5_ss_file$header$HalfMass, NA)
  expect_equal(hdf5_ss_file$header$SmoothingN, 1)
  expect_equal(length(names(hdf5_ss_file$header)), ss_file_header_length)

  expect_equal(eagle_ss_file$header$Centre, c(18316,61583,38667))
  expect_equal(eagle_ss_file$header$HalfMass, NA)
  expect_equal(eagle_ss_file$header$SmoothingN, 1)
  expect_equal(length(names(eagle_ss_file$header)), ss_file_header_length)

  expect_equal(magneticum_ss_file$header$Centre, NA)
  expect_equal(magneticum_ss_file$header$HalfMass, 1007671144)
  expect_equal(magneticum_ss_file$header$SmoothingN, 1)
  expect_equal(magneticum_ss_file$header$Alignment, "Specified")
  expect_equal(length(names(magneticum_ss_file$header)), ss_file_header_length)

  expect_equal(horizon_ss_file$header$Centre, NA)
  expect_equal(horizon_ss_file$header$HalfMass, NA)
  expect_equal(horizon_ss_file$header$SmoothingN, 3)
  expect_equal(length(names(horizon_ss_file$header)), ss_file_header_length)

  expect_equal(illustris_ss_file$header$Centre, NA)
  expect_equal(illustris_ss_file$header$HalfMass, NA)
  expect_equal(illustris_ss_file$header$SmoothingN, 1)
  expect_equal(length(names(illustris_ss_file$header)), ss_file_header_length)

})

