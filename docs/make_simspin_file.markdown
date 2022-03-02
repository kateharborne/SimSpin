---
layout: default
title: make_simspin_file
parent: Documentation
nav_order: 1
last_modified_date: "Wed, 16 February 2022 17:10:00 AWST"
---

# Making the input file for SimSpin

In order to build a data cube with `SimSpin`, we need to construct an input file.
{: .fs-5 .fw-300 .pb-2 }

The purpose of this input file is to prepare all the information of the supplied simulation in a consistent format.
It also aims to prepare the computationally expensive steps only once (i.e. the generation of spectra for each star).
Once we have generated a `SimSpin` file, that file can be used as input to `build_datacube` many times. 
{: .fw-300 }


[See an example](#example){: .btn .btn-purple }

---

The following code shows the default parameters used in the `make_simspin_file` function. Calling the function without specifying anything other than the required input `filename` will produce a SimSpin file saved at the same directory location as the input simulation file with the following defaults. 

```R
make_simspin_file(filename,                         # REQUIRED input file 
                  disk_age = 5, disk_Z = 0.024,     # for N-body disk particles
                  bulge_age = 10, bulge_Z = 0.001,  # for N-body bulge particles
                  cores = 1, write_to_file = TRUE, 
                  output, overwrite = F,
                  template = "BC03lr",              # template choice for spectra
                  centre = NA, half_mass = NA,      # alignment choice
                  sph_spawn_n = 1)                  # gas smoothing choice

```

---

## Parameters

| **filename**      	| The path to the snapshot file containing the galaxy of interest. Note that this file can be a Gadget binary file or an HDF5 file directly from a simulation. In the case that your simulation outputs data in an alternative format, either see [this example](examples/generating_hdf5) for more information about how to set up an HDF5 file for `SimSpin`, or [file an issue on GitHub](https://github.com/kateharborne/SimSpin/issues/new/choose) to see the inclusion of your format in the `SimSpin` code directly.                                                                                                                                                                                                                                                                                                                                               	|
| **cores**         	| The number of cores across which to multi-thread the problem.                                                                                                                                                                                                                                                                                                                                                     	|
| **disk_age**      	| The age of the disk particles in Gyr, used if the input file is an N-body model.                                                                                                                                                                                                                                                                                                                                  	|
| **bulge_age**     	| The age of the bulge particles in Gyr, used if the input file is an N-body model.                                                                                                                                                                                                                                                                                                                                 	|
| **disk_Z**        	| The metallicity of the disk particles in units of solar, used if the input file is an N-body model.                                                                                                                                                                                                                                                                                                               	|
| **bulge_Z**       	| The metallicity of the bulge particles in units of solar, used if the input file is an N-body model.                                                                                                                                                                                                                                                                                                              	|
| **template**      	| The stellar templates from which to derive the SEDs for each stellar particle. Options include "BC03lr" (GALEXEV low resolution, Bruzual & Charlot 2003), "BC03hr" (GALEXEV high resolution, Bruzual & Charlot 2003) and "EMILES" (Vazdekis et al 2016).                                                                                                                                                          	|
| **write_to_file** 	| Boolean to specify whether the list produced should be written to a ".Rdata" file or output to the environment. Default is TRUE, so that files can be re-observed without having the generate spectra each time.                                                                                                                                                                                                  	|
| **output**        	| The path at which the output file is written. If not provided, file will be written at the location of the input filename with the addition of "_spectra.Rdata".                                                                                                                                                                                                                                                  	|
| **overwrite**     	| If true, and the file already exists at the output location, a new file will be written over the old one. Default is FALSE.                                                                                                                                                                                                                                                                                       	|
| **centre**        	| If simulation file contains all particles cutout from a box (rather than just particles from a single galaxy), you can specify the point around which the view should be centred. Numeric length = 3. Default is NA, in which case the system is centred around the median position of stellar particles.                                                                                                         	|
| **half_mass**     	| If simulation file contains all particles cutout from a box (rather than just particles from a single galaxy), you can the half-mass value at which the alignment function is run. Numeric length = 1. Default is NA, in which case half the total mass of the supplied simulation data is used.                                                                                                                  	|
| **sph_spawn_n**   	| Numeric describing the number of gas particles with which to sample the SPH smoothing length. Default is 1, which will not spawn additional gas particles. Increasing this value increases the number of particles used to model the gas distribution. This value may need to be tested for convergence depending on the resolution of the grid used to image the gas properties at the `build_datacube()` stage. 	|

## Example

```R
# Using a Gadget model stored within the SimSpin software:
> simulation_file = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")

# Writing a SimSpin file directly to the R environment with all other defaults:
> simspin_file = make_simspin_file(filename = simulation_file,    
                                   write_to_file = FALSE)

# Examining the SimSpin file structure, we can see that it contains:
> summary(simspin_file)
           Length Class      Mode
 header       8   -none-     list    
 star_part   12   data.frame list   
 gas_part     0   -none-     NULL   
 spectra      2   -none-     list   
 wave      1221   -none-     numeric

# `header` is a list that contains details of the file you've just created. This allows you to 
# recreate this file in the future using the same information. 
> names(simspin_file$header)
[1] "InputFile"        "OutputFile"       "Type"             "Template"         "Template_LSF"    
[6] "Template_waveres" "Origin"           "Date"        

# `star_part` is a data.frame that contains the stellar particle information for that simulation.
# As the input model was an N-body model, this is a summary of both the disk (`PartType2`) and 
# bulge (`PartType3`) particle properties. If the input model was a hydrodynamic simulation, this
# would just contain the stellar particle information (`PartType4`).
> names(simspin_file$star_part)
 [1] "ID"           "x"            "y"            "z"            "vx"           "vy"          
 [7] "vz"           "Mass"         "sed_id"       "Metallicity"  "Age"          "Initial_Mass"

# `gas_part` in this file is NULL as the input simulation was an N-body model. However, in the case
# that the input was a hydrodynamical simulation with gas particles, this would be a data.frame 
# similar to `star_part` but summarising the gas particle properties instead. 
> simspin_file[["gas_part"]]
NULL

# `spectra` is a list where each element describes the flux at a given wavelength for one group of 
# stellar particles. These have been generated within the function for a distinct group of stellar
# particle ages and metallicities using the specified `template` spectra. Mapping between stellar 
# particles and their respective spectra is given by the `sed_id` within the `star_part` 
# data.frame. Each number `sed_id` corresponds to the element within `spectra` list. 
> simspin_file[["spectra"]]
[[1]]
   [1] 7.659288e-10 1.176103e-09 1.558386e-09 2.104659e-09 2.769032e-09 3.559435e-09 4.178723e-09
   [8] 5.006696e-09 6.134861e-09 8.460731e-09 1.130473e-08 1.429757e-08 5.147714e-08 6.069533e-08
  [15] ... [1221]
[[2]]
   [1] 1.089325e-08 1.397532e-08 1.675843e-08 2.010100e-08 2.321242e-08 2.587942e-08 2.695232e-08
   [8] 2.921529e-08 3.182460e-08 4.229105e-08 4.810211e-08 5.303458e-08 1.067975e-07 1.186373e-07
  [15] ... [1221]

# Because the input simulation was an N-body model, we only have two associated spectra - one for
# the disk particles (with `disk_age = 5` and `disk_Z = 0.024`) and one for the bulge particles 
# (with `bulge_age = 10` and `bulge_Z = 0.001`). 

# Finally `wave` contains the wavelength values in units of angstrom for all of the spectra 
# contained in the `spectra` list.
> simspin_file[["wave"]]
   [1]    91    94    96    98   100   102   104   106   108   110   114   118   121   125   127   128
  [17]   131   132   134   137   140   143   147   151   155   159   162   166   170   173   177   180
  [33] ... [1221]
```

