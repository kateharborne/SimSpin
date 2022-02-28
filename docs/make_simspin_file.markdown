---
layout: default
title: make_simspin_file
parent: Documentation
nav_order: 1
last_modified_date: "Wed, 16 February 2022 15:57:00 AWST"
---

# Making the input file for SimSpin

In order to build a data cube with `SimSpin`, we need to construct an input file. 
The purpose of this file is to prepare all the information of the supplied simulation into a consistent format and pree
It also aims to run the computationally expensive steps only once (i.e. the generation of spectra for each star).
Once we have generated a `SimSpin` file, that file can be used as input to `build_datacube` many times. 

```R
make_simspin_file(filename,
                  cores = 1,
                  disk_age = 5,
                  bulge_age = 10,
                  disk_Z = 0.024,
                  bulge_Z = 0.001,
                  template = "BC03lr",
                  write_to_file = TRUE,
                  output,
                  overwrite = F,
                  centre = NA,
                  half_mass = NA,
                  sph_spawn_n = 1)

```

## Parameters

| **filename**      	| The path to the snapshot file containing the galaxy of interest. Note                                                                                                                                                                                                                                                                                                                                             	|
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