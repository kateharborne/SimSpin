---
title: SimSpin
layout: default
filename: basic_usage.md
---

*Following this [installation](https://kateharborne.github.io/SimSpin/installation), you will have access to the range of functions within the `SimSpin` package. In this section, we take you through a simple mock observation routine using the default parameters for each function.*

Once you have installed `SimSpin` successfully, you can build your first mock data cube.
There are four steps to this process: 
    1. [making the input file](#step-1),
    1. [defining the conditions of the observing telescope and target galaxy](#step-2),
    1. generating the observation, 
    1. writing the output to FITS. 

We take you through each of these steps in the simplest form below.
More extensive examples can be found in ... 

## Step 1
### Making the input file

The main input for `SimSpin` is a galaxy simulation. 
This may be either an isolated galaxy model (*N*-body or hydro-dynamical) or a region cut out from a cosmological simulation. 
We define groups of particles or cells in each simulation according to the convention defined by Gadget [(Springel, 2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract):

* \code{PartType0} - gas,
* \code{PartType1} - dark matter, 
* \code{PartType2} - $N$-body disk, 
* \code{PartType3} - $N$-body bulge, 
* \code{PartType4} - stars, 
* \code{PartType5} - boundary/black holes. 

We begin by running the [`make_simspin_file`](https://kateharborne.github.io/SimSpin/docs/make_simspin_file) function. 
This will organise the provided data into a consistent format for `SimSpin` processing and assign SEDs to stellar particles given their individual ages and metallicities. 
We will discuss the variety of options available for each of the supported simulation input types within ... 

For this basic usage example, we will use one of the simulations included in the package:

```R
simulation_data = system.file("extdata", "SimSpin_example_Gadget", 
                              package = "SimSpin")
                              
simspin_data    = make_simspin_file(filename = simulation_data, 
                                    disk_age = 5, # ages are assigned in Gyr
                                    bulge_age = 10, 
                                    disk_Z = 0.024, # metallicities are wrt solar 
                                    bulge_Z = 0.001,
                                    template = "BC03lr", # SSP template used
                                    write_to_file = FALSE) 
```

The output of this function will be a list that can either be written to an \code{.Rdata} file or written to an environment variable within your R session. 
In the case above, we assign the output to the environment variable \code{`simspin\_data'}. 

Despite accepting a variety of different simulation inputs, the output of the \makesimspinfile{} function will always be the same. 
The format will be a list containing 4 elements: 

```R
typeof(simspin_file) 
# [1] "list"

summary(simspin_file) 
#           Length Class      Mode   
# star_part   12   data.frame list   
# gas_part     0   -none-     NULL   
# spectra      2   -none-     list   
# wave      1221   -none-     numeric
```

While the length of each element may change depending on the input simulation type or observation properties requested, the element structure of the output will remain consistent. 
In this case, because the input file is an *N*-body model containing both bulge and disk particles, the resulting element `spectra` is a `list` with *2* entries (with one spectra assigned for the bulge particles and the second to the disk particles).
For a hydrodynamical model, the  `spectra` element may have many more entries due to the broader variety in stellar ages and metallicities, but will remain an element of mode `list`. 
Once a file has been produced for a single simulated galaxy, the same file can be used to produce any number of different observations.
*Bear in mind that the resolution and wavelength range of the spectral templates used to generate the file need to overlap with the resolution and wavelength range of the observing telescope.*

## Step 2
### Defining the observation properties

`SimSpin` acts as a virtual telescope wrapper. 
You can choose to observe your galaxy model in a variety of different ways with any integral field unit (IFU) instrument. 
This requires you to set two distinct groups of properties - the properties of the instrument used to take the observation i.e. the `telescope()`, and the properties of the object under scrutiny i.e. the `observing_strategy()`. 

The properties are split in this way to enable a suite of observations to be generated in a straightforward manner. 
It is common that an observer will wish to observe a suite of galaxies using the same telescope, but may want to iterate over a number of projected inclinations, distances or seeing conditions. 
Hence, we have split the description classes for the observing telescope and observed object into two to enable this method more smoothly. 

#### telescope()
This class describes the properties of the instrument used for the observation.
Several current instruments are included within the class (named by either the telescope or the survey that used them i.e. SAMI, MaNGA, CALIFA and MUSE) but you can also define your own instrument by specifying the `type = 'IFU'`. 
An example of a user defined instrument is shown below, with the default parameters specified:

```R
ifu = telescope(type = "IFU", # other options include pre-defined instruments
                method = "spectral", # to generate spectral data cube
                fov = 15, # field-of-view in arcsec
                aperture_shape = "circular", # shape of fov
                wave_range = c(3700,5700), # given in angstroms
                wave_centre = 4700, # central wavelength given in angstroms
                wave_res = 1.04, # wavelength resolution in angstroms
                spatial_res = 0.5, # size of spatial pixels in arcsec
                filter = "r", # filter for calc of particle luminosity
                lsf_fwhm = 2.65, # full-width half-max of line-spread-function
                signal_to_noise = 10) # s/n ratio per angstrom
```

This will return a `list` that contains 12 elements. Here we have precomputed properties that will be necessary for further computation steps. The units of each returned element have been added in the final column.
```R
typeof(ifu) 
# [1] "list"

summary(ifu) 
#                 Length Class      Mode        
# type            1      -none-     character   # describes the telescope
# fov             1      -none-     numeric     # arcsec
# method          1      -none-     character   # describes the type of output cube
# aperture_shape  1      -none-     character   # shape of fov
# wave_range      2      -none-     numeric     # angstrom
# wave_centre     1      -none-     numeric     # angstrom
# wave_res        1      -none-     numeric     # angstrom
# spatial_res     1      -none-     numeric     # arcsec per pixel
# filter          2      data.frame list        # filter wavelengths and response
# lsf_fwhm        1      -none-     numeric     # angstrom 
# signal_to_noise 1      -none-     numeric     # dimensionless fraction
# sbin            1      -none-     numeric     # number of spatial bins across fov
```

#### observing_strategy()

We describe the conditions of the object begin observed using this function. 
For example, the level of atmospheric blurring can be modified, the angle from which the object is observed can be changed and the distance at which the galaxy is projected can be adjusted. An example with the default parameters specified is shown here:

```R
strategy = observing_strategy(z = 0.1, # redshift distance the galaxy is placed
                              inc_deg = 70, # projection (0 - face-on, 90 - edge-on)
                              twist_deg = 0, # angle of observation (longitude)
                              blur = F, # whether seeing conditions are included
                              fwhm = 2, # the fwhm of the blurring kernel in arcsec
                              psf = "Gaussian") # the shape of the blurring kernel
``` 

This will return a `list` that contains 4 elements. If `blur = T`, an additional two elements will be included to describe the point-spread-function shape and size. Similar to the previous output, this class just stores precomputed properties relevant to the observed object that will be necessary in further processing steps. The units of each returned element are listed in the final column.

```R
typeof(strategy) 
# [1] "list"

summary(strategy) 
#                Length Class  Mode   
# z              1      -none- numeric # redshift distance 
# inc_deg        1      -none- numeric # degrees
# twist_deg      1      -none- numeric # degrees
# blur           1      -none- logical # describes whether seeing is included
```

With these properties specified, we have fully described the parameters necessary to construct the observation. 
