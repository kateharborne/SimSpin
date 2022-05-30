---
layout: default
title: make_simspin_file
parent: Documentation
nav_order: 1
last_modified_date: "Fri, 27 May 2022 15:57:00 AWST"
---

# Making the input file for SimSpin

In order to build a data cube with `SimSpin`, we need to construct an input file.
{: .fs-5 .fw-300 .pb-2 }

The purpose of this input file is to prepare all the information of the supplied simulation in a consistent format.
It also aims to prepare the computationally expensive steps only once (i.e. the generation of spectra for each star).
Once we have generated a `SimSpin` file, that file can be used as input to `build_datacube` many times. 
{: .fw-300 }


[See an example](#example){: .btn .btn-purple }
[See the source code](https://github.com/kateharborne/SimSpin/blob/d020398fb66274443bb2f70ea1fdd8346c4476ae/R/make_simspin_file.R#L57){: .btn .btn-purple }

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

## Input Parameters

| `filename`     	| The path to the snapshot file containing the galaxy of interest. Note that this file can be a Gadget binary file or an HDF5 file directly from a simulation. In the case that your simulation outputs data in an alternative format, either see [this example](examples/generating_hdf5) for more information about how to set up an HDF5 file for `SimSpin`, or [file an issue on GitHub](https://github.com/kateharborne/SimSpin/issues/new/choose) to see the inclusion of your format in the `SimSpin` code directly.                                                                                                                                                                                                                                                                                                                                               	|
| `cores`         	| The number of cores across which to multi-thread the problem.                                                                                                                                                                                                                                                                                                                                                     	|
| `disk_age`      	| The age of the disk particles in Gyr, used only if the input file is an N-body model.                                                                                                                                                                                                                                                                                                                                  	|
| `bulge_age`     	| The age of the bulge particles in Gyr, used only if the input file is an N-body model.                                                                                                                                                                                                                                                                                                                                 	|
| `disk_Z`        	| The metallicity of the disk particles in units of solar, used only if the input file is an N-body model.                                                                                                                                                                                                                                                                                                               	|
| `bulge_Z`       	| The metallicity of the bulge particles in units of solar, used only if the input file is an N-body model.                                                                                                                                                                                                                                                                                                              	|
| `template`      	| The stellar templates from which to derive the SEDs for each stellar particle. Options include `"BC03lr"` (GALEXEV low resolution, Bruzual & Charlot 2003), `"BC03hr"` (GALEXEV high resolution, Bruzual & Charlot 2003) and `"EMILES"` (Vazdekis et al 2016). The details of these templates are shown in the [additional information](#notes) below.                                                                                                                                                         	|
| `write_to_file` 	| Boolean to specify whether the list produced should be written to a ".Rdata" file or output to the environment. Default is TRUE, so that files can be re-observed without having the generate spectra each time.                                                                                                                                                                                                  	|
| `output`        	| The path at which the output file is written. If not provided, file will be written at the location of the input filename with the addition of "_spectra.Rdata".                                                                                                                                                                                                                                                  	|
| `overwrite`     	| If true, and the file already exists at the output location, a new file will be written over the old one. Default is FALSE.                                                                                                                                                                                                                                                                                       	|
| `centre`        	| If simulation file contains all particles cutout from a box (rather than just particles from a single galaxy), you can specify the point around which the view should be centred. Numeric length = 3. Default is NA, in which case the system is centred around the median position of stellar particles.                                                                                                         	|
| `half_mass`     	| If simulation file contains all particles cutout from a box (rather than just particles from a single galaxy), you can the half-mass value at which the alignment function is run. Numeric length = 1. Default is NA, in which case half the total mass of the supplied simulation data is used.                                                                                                                  	|
| `sph_spawn_n`   	| Numeric describing the number of gas particles with which to sample the SPH smoothing length. Default is 1, which will not spawn additional gas particles. Increasing this value increases the number of particles used to model the gas distribution. This value may need to be tested for convergence depending on the resolution of the grid used to image the gas properties at the `build_datacube()` stage. 	|

---

## Output Value
The output of `make_simspin_file` is a *List* element.
If `write_to_file = T`, this list will be written to an .Rdata file at the location specified by `output`. Otherwise, it will be written as a variable to the environment.

The list will contain the following 5 elements:

1. `header` -  *List* element containing eight labelled properties for data transparency and reproducability.

    | `InputFile`        | the file path to the original simulation file. |
    | `OutputFile`       | the original location that the SimSpin file was written to. |
    | `Type`             | the classification of the input simulation file, one of `nbody`, `EAGLE`, or `Magneticum` (though this list is ever growing). |
    | `Template`         | the spectral template used to construct spectra for each stellar particle. |
    | `Template_LSF`     | the line-spread-function assoicated with the template spectra in angstrom. |
    | `Template_waveres` | the wavelength resolution given by the binning size in angstrom.  |
    | `Origin`           | the version of SimSpin used to build the file. |
    | `Date`             | the date and time at which the SimSpin file was built. |

1. `star_part` - A *data.table* element containing properties of each stellar particle contained in the simulation. (`NULL` if no star particles are present in the simulation.)

    | `ID` | the unique stellar particle ID number. |
    | `x`, `y` and `z` | the stellar particle positions in each of the 3 dimensions in kpc. |
    | `vx`, `vy`, `vz` | the stellar particle velocities in each of the 3 dimensions in km/s. |
    | `Mass` | the stellar particle mass in solar masses. |
    | `sed_id` | the index of the spectrum associated with each stellar particle (i.e. the respective row in the `spectra` element below.) |
    | `Metallicity` | the stellar metallicity of a given particle in units of solar metallicity. For `nbody` models, these values are assigned by the function. For any hydrodynamical simulations, these values are pulled from the simulation file. |
    | `Age` | the stellar age of each particle given in Gyr. For `nbody` models, these values are assigned by the function. For any hydrodynamical simulations, these values are pulled from the simulation file. |
    | `Initial_Mass` | the stellar mass of each star particle at the birth of the star (used for scaling spectra). For `nbody` models, these values are assigned by the function. We assume that the star has doubled in mass since its birth by default. For any hydrodynamical simulations, these values are pulled from the simulation file. |

1. `gas_part` - A *data.table* element containing properties of each gas particle contained within the simulation. (`NULL` if no gas particles are present in the simulation.)

    | `ID` | the unique gas particle ID number. |
    | `x`, `y` and `z` | the gas particle positions in each of the 3 dimensions in kpc. |
    | `vx`, `vy`, `vz` | the gas particle velocities in each of the 3 dimensions in km/s. |
    | `Mass` | the gas particle mass in solar masses. |
    | `SFR` | the instananeous star formation rate of the gas particle. |
    | `Density` | the density of the gas within the simulation given in solar masses per kpc^3. |
    | `Temperature` | the temperature of that gas particle in K. |
    | `SmoothingLength` | the size of the the SPH smoothing kernel in kpc. |
    | `Metallicity` | the smoothed gas mass of elements heavier than Helium divided by the mass of the particle. |
    | `Carbon` | the fraction of carbon within the particle, given as a fraction of its total mass. |
    | `Hydrogen` | the fraction of hydrogen within the particle, given as a fraction of its total mass. |
    | `Oxygen` | the fraction of oxygen within the particle, given as a fraction of its total mass. |

1. `spectra` - A *data.table* element containing the spectral flux associated with each particle, linked to each particle by the `sed_id` within the `star_part` element. Each row corresponds to the spectrum associated with an individual spectral ID. 

1. `wave` - *Numeric* element containing the wavelengths in angstrom for the spectra in the element above. 

---

## Example

Using a Gadget model stored within the SimSpin software, we can write a SimSpin file directly to the R environment with all other defaults:
```R
# Load a Gadget model...
simulation_file = system.file("extdata", "SimSpin_example_Gadget", package = "SimSpin")
# ... use to build a SimSpin file. 
simspin_file = make_simspin_file(filename = simulation_file,
                                 write_to_file = FALSE)
```
Examining the SimSpin file produced, we can see that it contains:
```R
summary(simspin_file)
#           Length Class      Mode
# header       8   -none-     list    
# star_part   12   data.frame list   
# gas_part     0   -none-     NULL   
# spectra      2   -none-     list   
# wave       842   -none-     numeric
```
- `header` is a list that contains details of the file you've just created. This allows you to recreate this file in the future using the same information. 
```R
names(simspin_file$header)
# [1] "InputFile"        "OutputFile"       "Type"             "Template"            
# [5] "Template_LSF"     "Template_waveres" "Origin"           "Date"        
```
- `star_part` is a data.frame that contains the stellar particle information for that simulation. As the input model was an N-body model, this is a summary of both the disk (`PartType2`) and bulge (`PartType3`) particle properties. If the input model was a hydrodynamic simulation, this would just contain the stellar particle information (`PartType4`).
```R
names(simspin_file$star_part)
#  [1] "ID"           "x"            "y"            "z"            "vx"           
#  [6] "vy"           "vz"           "Mass"         "sed_id"       "Metallicity"
# [11] "Age"          "Initial_Mass"
```
- `gas_part` in this file is NULL as the input simulation was an N-body model. However, in the case that the input was a hydrodynamical simulation with gas particles, this would be a data.frame similar to `star_part` but summarising the gas particle properties instead. 
```R
simspin_file[["gas_part"]]
# NULL
```
- `spectra` is a list where each element describes the flux at a given wavelength for one group of stellar particles. These have been generated within the function for a distinct group of stellar particle ages and metallicities using the specified `template` spectra. Mapping between stellar particles and their respective spectra is given by the `sed_id` within the `star_part` data.frame. Each number `sed_id` corresponds to the element within `spectra` list. 
Because the input simulation was an N-body model, we only have two associated spectra - one for the disk particles (with `disk_age = 5` and `disk_Z = 0.024`) and one for the bulge particles (with `bulge_age = 10` and `bulge_Z = 0.001`). 
```R
simspin_file[["spectra"]]
#[[1]]
#   [1] 7.659288e-10 1.176103e-09 1.558386e-09 2.104659e-09 2.769032e-09 3.559435e-09 
#   [7] 4.178723e-09 5.006696e-09 6.134861e-09 8.460731e-09 1.130473e-08 1.429757e-08 
#  [13] ... [1221]
#[[2]]
#   [1] 1.089325e-08 1.397532e-08 1.675843e-08 2.010100e-08 2.321242e-08 2.587942e-08 
#   [7] 2.695232e-08 2.921529e-08 3.182460e-08 4.229105e-08 4.810211e-08 5.303458e-08 
#  [13] ... [1221]
```
- Finally `wave` contains the wavelength values in units of angstrom for all of the spectra contained in the `spectra` list.
```R
simspin_file[["wave"]]
#   [1]    91    94    96    98   100   102   104   106   108   110   114   118   121   125   
#  [15]   127   128   131   132   134   137   140   143   147   151   155   159   162   166   
#  [29] ... [1221]
```

---
## Notes

### Spectral Template Choice

Spectral templates provide a grid of modelled spectra for a known stellar population. The grid associates a stellar age and metallicity to each spectrum in the template. These are used throughout integral field spectroscopy to determine the stellar population properties of the observed galaxy.

In SimSpin, we follow this methodology in reverse. As opposed to fitting the spectrum to find the age and metallicity of the underlying stellar population, we use the stellar age and metallicity to assign a spectrum using the chosen template library. While the output is consistent between hydrodynamical and *N*-body models, we discuss the differences in methodology between these inputs below:

- **Hydrodynamical simulations** track the age and metallicity of each stellar particle. We use these parameters to select the spectrum that should be associated with each star. Individual particles are initially binned into groups of age and metallicity with width 20 Myr and 0.1 Z<sub>&#9737;</sub> respectively. Stellar particles that fall within each age-metallicity bin will be associated with a single spectrum. This simplification is made to reduce the number of spectra that  need to be stored and hence optimise memory usage within SimSpin.
- **N-body models** consider stellar particles as collisionless bulge and disk components, but do not model their the stellar properties of age and metallicity. These properties must be assigned before SimSpin can generate a spectrum.  Of course, the values chosen will change the shape of the associated spectrum and the relative brightness of the individual constituents. Commonly, we assign different age and luminosity values to each population of particles to represent the younger, metal-rich disk and older, metal-poor bulge. This choice is arbitrary and depends on the science in question. To avoid systematic differences, the age and metallicity of all components can be set equal. 

Available spectral templates include: 

| Name     | Age Steps | Age Range (Gyr) | Z Steps | Z Range (Z<sub>&#9737;</sub>)| &lambda; Steps (&#8491;) | &lambda; Range (&#8491;) | 
|----------|-----------|-----------------|---------|------------------------------|--------------------------|--------------------------|
| `BC03lr` | 221       | 0 - 20          | 6       | 0.0001 - 0.05                | 842                      | 91 - 20000               |
| `BC03hr` | 221       | 0 - 20          | 6       | 0.0001 - 0.05                | 6521                     | 91 - 20000               |
| `EMILES` | 53        | 0.03 - 14       | 12      | 0.0001 - 0.04                | 20356                    | 1680 - 20000             |

---

Now that you've built a SimSpin file, you can go one to set up the observations with a mock telescope...

[Next Step: telescope](telescope.markdown){: .btn .btn-purple }