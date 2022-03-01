---
layout: default
title: Basic Usage
permalink: /basic_usage/
description: "A walkthrough demonstrating how to generate your first mock observation."
nav_order: 3
last_modified_date: "Wed, 16 February 2022 13:57:00 AWST"
---

# A walkthrough SimSpin's functions
{: .no_toc }

Here, we demonstrate how to create your first mock observation. 
At each step, we take you through the simplest routine using the default parameters for each function.
{: .fs-5 .fw-300 .pb-2 }

---

There are four steps to build your first `SimSpin` mock data cube: 
{: .fw-300 }

1. TOC
{:toc .pb-1}

We take you through each of these steps in the simplest form below.
More extensive examples can be found in [Examples](https://kateharborne.github.io/SimSpin/examples).
{: .fw-300 }

---

## Making the input file

The main input for `SimSpin` is a galaxy simulation. 
This may be either an isolated galaxy model (*N*-body or hydro-dynamical) or a region cut out from a cosmological simulation. 
We define groups of particles or cells in each simulation according to the convention defined by Gadget [(Springel, 2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract):

* `PartType0` - gas,
* `PartType1` - dark matter, 
* `PartType2` - *N*-body disk, 
* `PartType3` - *N*-body bulge, 
* `PartType4` - stars, 
* `PartType5` - boundary/black holes. 
{: .lh-tight }

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

The output of this function will be a list that can either be written to an {.Rdata} file or written to an environment variable within your R session. 
In the case above, we assign the output to the environment variable {`simspin_data'}. 

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

---

## Defining the observation properties

`SimSpin` acts as a virtual telescope wrapper. 
You can choose to observe your galaxy model in a variety of different ways with any integral field unit (IFU) instrument. 
This requires you to set two distinct groups of properties - the properties of the instrument used to take the observation i.e. the `telescope()`, and the properties of the object under scrutiny i.e. the `observing_strategy()`. 

The properties are split in this way to enable a suite of observations to be generated in a straightforward manner. 
It is common that an observer will wish to observe a suite of galaxies using the same telescope, but may want to iterate over a number of projected inclinations, distances or seeing conditions. 
Hence, we have split the description classes for the observing telescope and observed object into two to enable this method more smoothly. 

### telescope()
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

### observing_strategy()
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

---

## Generating the observation

Once the particulars of the observation have been defined, we can go ahead and make the observation. 
This means feeding our prepared parameters into the `build_datacube()` function. 
Here, we choose not to write the output directly to file (`write_fits = F`) and output to a variable called `gadget_cube` instead. 

```r
gadget_cube = build_datacube(simspin_file       = simspin_data, 
                             telescope          = ifu,
                             observing_strategy = strategy,
                             write_fits         = F)
```

The returned variable will always be a list containing 4 elements:

1. the first always being the data cube,
2. followed by the observation summary,
3. raw simulation images and
4. mock observation images generated from collapsing the data cube. 

```r
> typeof(gadget_cube)
# [1] "list"
> summary(gadget_cube)
#                 Length  Class  Mode   
# spectral_cube   1731600 -none- numeric
# observation          33 -none- list   
# raw_images            3 -none- list   
# observed_images       0 -none- NULL   
```

**add text about the structure of each element in the list**

The overall format of this output will be consistent.
However, the names of individual elements change to reflect the variety in the requested properties specified by [`telescope`]() and [`observing_strategy`]().
* For example,  `method = 'spectral'` will return a variable `spectral_cube` as its first element; specifying instead `method = 'velocity'` will return a `velocity_cube`.
* Similarly, if we are working in velocity mode with `mass_flag = T`, the images within the `raw_images` and `observed_images` elements will include an array called `mass_image`, rather than `flux_image`.

---

## Writing the output to FITS

Using the `write_simspin_FITS()` function, the new observation can be written to a file that can be shared or re-analysed using observational pipelines.
`SimSpin` can write such a FITS file using the variable produced by `build_datacube()`.
A few further descriptive details are suggested as input for users of the file later on, but the function defaults to sensible suggestions, as shown in the code below.

```r
write_simspin_FITS(simspin_datacube   = gadget_cube,        # build_datacube() output
                   output_file        = "gadget_cube.FITS", # filename of output file
                   input_simspin_file = simulation_data,    # filename of input sim     
                   object_name        = "GalaxyID_unknown", # name of observed object
                   telescope_name     = "Telescope",        # name of telescope
                   instrument_name    = "IFU",              # name of instrument 
                   observer_name      = "Anonymous")        # name of the observer
```

The output of running this code will be a FITS file written at the location specified by `output_file`. 
To explore this file within your R session, we suggest using the package `Rfits`, which can be downloaded using the instructions at the link provided [here](https://github.com/asgr/Rfits).
We use the `Rfits_info()` function to examine the structure of the file that has been written:

```r
library("Rfits")
fits_summary = Rfits_info(filename = "gadget_cube.FITS")
fits_summary$summary

# [1] "SIMPLE  =                    T / file does conform to FITS standard"
# [2] "XTENSION= 'IMAGE   '           / IMAGE extension"                   
```

Here, we can clearly see that the file contains six HDU: a header table and a number of images. The exact number of HDU can vary depending on the type of observation made by `build_datacube()`. In all cases, HDU `ext = 1` will always be the header information and HDU `ext = 2` will always be the data cube produced. In this case, for a spectral data cube there are a series of `raw_images` output into subsequent HDUs. In the case of a velocity data cube, we would also see several `observed_images` included in the output.

The first HDU (`ext = 1`) is a header element describes the properties of the observation. This acts as a record for reproducing that identical observation again in the future. Information that is recorded in the header is listed below. These are consistent across all FITS files written using this function, regardless of variety in `build_datacube` options.

```r
head = Rfits_read_header("~/Desktop/gadget_cube.FITS")
head$header
# [1] "SIMPLE  =                    T / file does conform to FITS standard"             
# [2] "BITPIX  =                    8 / number of bits per data pixel"                  
# [3] "NAXIS   =                    0 / number of data axes"                            
# [4] "EXTEND  =                    T / FITS dataset may contain extensions"            
# [5] "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy"
# [6] "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H" 
# [7] "ORIGIN  = 'SimSpin '           / Mock observation"                               
# [8] "TELESCOP= 'Telescope '         / Telescope name"                                 
# [9] "INSTRUME= 'IFU     '           / Instrument used"                                
#[10] "RA      =                    0 / [deg] 11:41:21.1 RA (J2000) pointing"           
#[11] "DEC     =                    0 / [deg] -01:34:59.0 DEC (J2000) pointing"         
#[12] "EQINOX  =                 2000 / Standard FK5"                                   
#[13] "RADECSYS= 'FK5     '           / Coordinate system"                              
#[14] "EXPTIME =                 1320 / Integration time"                               
#[15] "MJD-OBS =             58906.11 / Obs start"                                      
#[16] "DATE-OBS= '2021-11-16 14:47:20' / Observing date"                                
#[17] "UTC     =                 9654 / [s] 02:40:54.000 UTC"                           
#[18] "LST     =             30295.18 / [s] 08:24:55.178 LST"                           
#[19] "PI-COI  = 'UNKNOWN '           / PI-COI name."                                   
#[20] "OBSERVER= 'Anonymous'          / Name of observer."                              
#[21] "REDSHIFT=                  0.1 / Observed redshift."                             
#[22] "PIPEFILE= 'gadget_cube.FITS'   / Filename of data product"                       
#[23] "BUNIT   = 'erg/s/cm**2'        / Angstrom"                                       
#[24] "ARCFILE = '/Users/path/Library/R/4.0/library/SimSpin/extdata/SimSpin_exampl'"
#[25] "DATAMD5 = '4aece79473a5c88f6533382655e948bf' / MD5 checksum"                     
#[26] "OBJECT  = 'GalaxyID_unknown'   / Original target."            

```

Observing details such as the name of the person who ran the file, the name of the object and the redshift at which is was observed are all included in this header. 
Importantly, we also record the name and location of the `SimSpin` file from which this observation was made (denoted as `ARCFILE`), enabling files to be reproduced in the future. 
Admittedly, several aspects to this header have to be invented as they do not have any physical meaning for mock observations (`RA` and `DEC`, for example), but are required for consistency with their observational counterparts. 
**We have made these *invented* properties clear to the user by stating `mocked` in the comment of the header parameters.**

Finally, we can go about reading in and working with the data cube itself. 
Each HDU comes along with its own header to describe the properties of the image contained i.e. the units and dimensions of the pixels. 
These details can be used to examine and process the data.

**code demonstration of reading in the spectral data cube**

Similar methodology can be used to examine both the 3-dimensional data cubes and the 2-dimensional images contained within sequential HDUs. 
The only difference is the number of CVAL units (three for cubes versus two for images).

*These steps can be used to generate a mock observation of a simulated galaxy. In this Basic Usage section, we have run the simplest recipe using the inbuilt defaults of each function. In the next section, we go into details of additional functionality that may be relevant for your work.*

