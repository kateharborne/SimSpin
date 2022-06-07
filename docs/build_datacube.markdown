---
layout: default
title: build_datacube
parent: Documentation
nav_order: 4
last_modified_date: "Tue, 7 June 2022 15:57:00 AWST"
---

# Constructing your data cube

Using the `build_datacube` function, we can now make a mock observation of a simulated galaxy. 
{: .fs-5 .fw-300 .pb-2 }

The output of this function will be a mock IFS data cube. The details of the observation will vary depending on the requested `method`, `telescope` and `observing_strategy`, but the format will be consistent across the board. In all cases, the output can be written to a FITS file, with headers that allow the user to reproduce these observations in the future. 
{: .fw-300 }


[See an example](#example){: .btn .btn-purple }
[See the source code](https://github.com/kateharborne/SimSpin/blob/bb371d9e4d981d1ebaba3aa07978bb61a2d434f0/R/build_datacube.R#L77){: .btn .btn-purple }

---

The following code shows the default parameters used in the `make_simspin_file` function. Calling the function without specifying anything other than the required input `filename` will produce a SimSpin file saved at the same directory location as the input simulation file with the following defaults. 

```R
build_datacube(simspin_file,                     # REQUIRED input SimSpin file      
               telescope,                        # REQUIRED telescope class description
               observing_strategy,               # REQUIRED observing_strategy class description
               method = "spectral",              # "spectral", "velocity", "gas" or "sf_gas"
               verbose = F, 
               write_fits = F, output_location,  # write to a FITS file? If so, where?
               object_name = "GalaxyID_unknown", # information for FITS file headers
               telescope_name = "SimSpin",
               observer_name = "Anonymous",
               split_save = F, cores = 1,        
               mass_flag = F)                    # mass- or luminosity-weighted?

```

---

## Input Parameters

| **simspin_file**       | The path and filename of the SimSpin file OR the output environment list element constructed using [`make_simspin_file`](make_simspin_file.markdown). |
| **telescope**          | An object of the [`telescope`](telescope.markdown) class that describes the properties of the observing telescope (e.g. the field of view, the spatial resolution, the wavelength resolution, etc.).   |
| **observing_strategy** | An object of the [`observing_strategy`](observing_strategy.markdown) class that defines the properties of the observed simulation (e.g. the projected distance, the relative pointing, the seeing conditions, etc.). |
| **method**             | A character string that describes whether to generate a spectral data cube or a kinematic data cube. Options include `"spectral"` and `"velocity"` to produce these for the stellar component of the simulation. If wishing to analyse the kinematics of the gas, you can also specify `"gas"` or `"sf_gas"` (i.e. star forming gas) to produce a kinematic data cube of this component alone.   |
| **verbose**            | Should the code output text to describe its progress through the build? By default, this will be `FALSE`. |
| **write_fits**         | Should the code write the output data cube to a FITS file? By default, this is `FALSE` and output will be stored in an environment variable. |
| **output_location**    | This is an optional parameter that will only be used if `write_fits = TRUE`. Code behavior will depend on the text input. If location is given as just a path to a directory location, the file name will be automatically generated based on the input `simspin_file` name and the observing conditions. If the full path and filename are given (ending in ".fits" or ".FITS") the generated file will be written to that location with that name. Else, in the case that `write_fits = TRUE` but no `output_location` is specified, the output FITS file will be written to the same directory as the input SimSpin file. |
| **object_name**        | String used to describe the name of the object observed in the FITS file header. |
| **telescope_name**     | String used to describe the name of the telescope used for the observation in the FITS file header. |
| **observer_name**      | String used to describe the name of the person who generated the observation in the FITS file header. |
| **split_save**         | Only used when `write_fits = TRUE`. Should the output FITS be saved as one file with multiple HDUs to describe the output cube and observed images (`split_save = FALSE`)?  Or should each cube/image be saved to a seperate file (`split_save = TRUE`)? In this case, the file name root will be taken from the state of `output_location` and descriptive names will automatically be appended to individual files (i.e. "_spectral_cube.FITS", "_velocity_image.FITS", etc.). |
| **cores**              | Specifiying the number of cores to run the code with parallellism turned on. |
| **mass_flag**          | When `mass_flag = TRUE` and `method = "velocity"`, the output kinematic properties will be mass weighted rather than luminosity weighted. By default, this parameter is set to `FALSE`. |

---

## Example 