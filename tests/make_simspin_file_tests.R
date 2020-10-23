# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the make_simspin_file.R code

ss_gadget = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")
ss_hdf5 = system.file("extdata", "SimSpin_example_HDF5.hdf5", package = "SimSpin")
ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")

# Test that they run successfully without error
gadget_test = make_simspin_file(ss_gadget)
hdf5_test   = make_simspin_file(ss_hdf5)
eagle_test  = make_simspin_file(ss_eagle)
