---
layout: default
title: telescope
parent: Documentation
nav_order: 3
last_modified_date: "Wed 24 Aug 2022 15:57:00 AWST"
---

# Defining the properties of the telescope

Recent Updates
{: .label .label-purple } 

In order to build a data cube, we need to describe the properties of the telescope used to observe the model. 
{: .fs-5 .fw-300 .pb-2 }

Using the `telecope` function, a single telescope object is generated with a set number of expected properties. These properties are used with the `observing_strategy` to describe a specific observation. You can fully specify the particulars of your chosen telescope or you can select from a number of inbuilt telescopes designed to mimic current IFS surveys. 
{: .fw-300 }

---

The following code shows the default parameters used in the `telescope` function. Calling the function without specifying any input will produce a telescope object with the following properties:

```R
telescope(type="IFU",                # specify a type or define your own using "IFU"
          fov=15,                    # diameter of the field-of-view in arcsec          
          aperture_shape="circular", 
          wave_range=c(3700,5700),   # in angstrom
          wave_centre,               # wave_centre is auto-generated if not supplied
          wave_res=1.5,              # wavelength resolution in angstrom
          spatial_res=0.5,           # spatial resolution in arcsec
          filter="g",                # luminosities output in this band  
          lsf_fwhm=3.2,              # spectral uncertainty due to instrument in angstrom
          signal_to_noise = NA)      # target signal-to-noise ratio, default is no noise
```

[Input parameters](#input-parameters){: .btn .btn-purple }
[Output parameters](#output-parameters){: .btn .btn-purple }
[See an example](#example){: .btn .btn-purple }
[See the source code](https://github.com/kateharborne/SimSpin/blob/d020398fb66274443bb2f70ea1fdd8346c4476ae/R/telescope.R#L44){: .btn .btn-purple }

*As of version 2.2.0, the `method` input parameter has been moved directly to the [`build_datacube`](build_datacube.markdown) function. For backwards compatibility, this parameter can still be specified here, but a warning will be issued.*
{: .fw-150 }

---

## Input Parameters

| `type`              | A character string describing the type of telescope you wish to create. Default is `"IFU"`, in which case the telescope properties are described by the input parameters to this function. Other input types include: `"SAMI"`, `"MaNGA"`, `"CALIFA"` and `"MUSE"`. With these, the telescope properties are set to the default parameters specified by the respective instrument.    |
| `fov`               | A numeric describing the diameter of the field-of-view of the telescope in arcsec.    |
| `aperture_shape`    | A character string describing the shape of the field-of-view. Options include: `"circular"`, `"hexagonal"`, and `"square"`.   |
| `wave_range`        | Two numeric parameters describing the beginning and end of the wavelength range visible to the telescope, given in angstrom. Format is expected `c(wave_start, wave_end)`, and errors will be issued if not in this format.   |
| `wave_centre`       | A numeric parameter describing the centre of the wavelength range in angstrom. If no `wave_centre` is provided, SimSpin assumes that the central wavelength is at the centre of the wavelength range provided above.   |
| `wave_res`         | A numeric that gives the resolution of the wavelength axis of the telescope in angstrom. For a given `wave_range`, this resolution defines the width of the bins along that range.    |
| `spatial_res`       |  A numeric that gives the resolution of the spatial axis of the telescope in arcsec. For a given `fov`, this resolution defines the width of the pixels within the observation produced.  |
| `filter`            | A character string that describes which band to use to calculate luminoisity at the telescope. Options include: `"r"`, `"g"`, `"u"`, `"i"`, and `"z"` and use the SDSS filter bandpass functions to compute the luminosity that would be detected in that band.   |
| `lsf_fwhm`          | A numeric describing the full-width half-maximum of the Gaussian line-spread-function of the telescope in units of angstrom, i.e. the spectral uncertainty of the underlying spectrograph used in the observation.    |
| `signal_to_noise`   | A numeric describing the target signal-to-noise ratio at any pixel.   |

---

## Output Parameters

The output of `telescope` is a *List* element that will be stored as a variable to the environment. 

The list will contain the following 12 elements:
{: .fw-300 }

1. `type` - *Character* element recording the requested telescope type. Must be one of `"IFU"`, `"SAMI"`, `"MaNGA"`, `"MUSE"` or `"Hector"`.

1. `fov` - *Numeric* element describing the diameter of the field-of-view in arcseconds. 

1. `aperture_shape` - *Character* element recording the shape of the field-of-view. Must be one of `"circular"`, `"hexagonal"` or `"square"`.

1. `wave_range` - *Numeric* element of length 2 describing the beginning and end of the wavelength axis in angstrom. 

1. `wave_centre` - *Numeric* element describing the centre of the wavelength axis in angstrom. 

1. `wave_res` - *Numeric* element describing the width of each bin on the wavelength axis in angstrom. 

1. `spatial_res` - *Numeric* element describing the width of each spatial pixel in arcseconds.

1. `filter_name` - A *character* string element describing the name of the filter used. Could be "r_SDSS", "g_SDSS", "u_SDSS", "i_SDSS" or "z_SDSS". 

1. `filter` - A *data.table* element with two columns:

    | `wave` | *Numeric* element describing the wavelengths at which the spectral reponse has been computed for a given spectral filter. |
    | `response` | *Numeric* element describing the spectral response of that requested filter, given as a fraction. |

1. `lsf_fwhm` - *Numeric* element describing the line-spread-function full-width-half-maximum of the spectrograph in the specified telescope. Given in units of angstrom. 

1. `signal_to_noise` - *Numeric* element describing the maximum signal-to-noise ratio per pixel. 

1. `sbin` - *Numeric* element describing the number of spatial pixels across the diameter of the aperture. 
{: .bg-grey-lt-000 }
---

## Example

Start by building a telescope using the default parameters:
```R
scope = telescope()
```
We can inspect the object produced and see that the default parameters have been stored as a list.
```R
summary(scope)
#                 Length Class      Mode     
# type            1      -none-     character
# fov             1      -none-     numeric  
# aperture_shape  1      -none-     character
# wave_range      2      -none-     numeric  
# wave_centre     1      -none-     numeric  
# wave_res        1      -none-     numeric  
# spatial_res     1      -none-     numeric  
# filter          2      data.table list     
# lsf_fwhm        1      -none-     numeric  
# signal_to_noise 1      -none-     numeric  
# sbin            1      -none-     numeric  
``` 
To inspect individual properties of our telescope, use the named elements in one of two ways to achieve the same result:

```R
scope$type
# [1] "IFU"

scope[["type"]]
# [1] "IFU"
```
For pre-defined `type` classes, certain parameters cannot be modified. Values highlighted in purple below (for example, the field-of-view for a MaNGA telescope) can be selected by adding additional parameters to the inputs. 

| *telescope parameter* | *units* | **SAMI** | **MaNGA** | **CALIFA** | **MUSE** | **Hector** |
| --------------------- | ------- | -------- | --------- | ---------- | -------- | ---------- |
| `fov` | arcsec | 15 | <mark style="background-color: #e3d0fe"><i>n</i> = 12, 17, 22, 27 or 32</mark> | 74 | <mark style="background-color: #e3d0fe"><i>n</i> < 60</mark> | 30 |
| `aperture_shape` |   | "circular" | "hexagonal" | "hexagonal" | "square" | "hexagonal" |
| `wave_range` | &#8491;<sup>1</sup> | 3750 - 5750 | 3600 - 6350 | 3700 – 4750 | 4700.15 - 9351.4 | 3720 - 5910 |
| `wave_centre` | &#8491;<sup>1</sup> | 4800 | 4700 | 4225 | 6975 | 4815 |
| `wave_res` | &#8491; | 1.04 | 1.04 | 2.7 | 1.25 | 1.60 |
| `spatial_res` | arcsec/pixel | 0.5 | 0.5 | 1 | <mark style="background-color: #e3d0fe">0.2 (WFM) or 0.025 (NFM)</mark> | 0.2 |
| `lsf_fwhm` | &#8491; | 2.65 | 2.85 | 2.7 | 2.51 | 1.3 |

<sup>1</sup> *Wave ranges and central wave lengths are quoted for the blue arm of each spectrograph as this is the wavelength range across which the kinematics are commonly measured.* 

For example, we can specify a variety of different field-of-views for the MaNGA telescope:

```r
# MaNGA telescope with 22" field-of-view
manga_22 = telescope(type="MaNGA", fov = 22)

# MaNGA telescope with 12" field-of-view
manga_12 = telescope(type="MaNGA", fov = 12)

```
The field-of-view (and hence the number of spatial bins, `sbin`) in each output list will change, but all other telescope properties remain consistent. We can check this with a verbose little loop:

```r
for (each in names(manga_12)){ 
    # for each element in the telescope list
  if (!all(manga_12[[each]] == manga_22[[each]])){  
    # if the parameters do not match, let the user know.
    cat("NOT EQUAL!", each, "\n")
  } else {
    # otherwise, all good!
    cat("EQUAL:", each, "\n") 
  }
} 

# EQUAL: type 
# NOT EQUAL! fov  
# EQUAL: aperture_shape 
# EQUAL: wave_range 
# EQUAL: wave_centre 
# EQUAL: spatial_res 
# EQUAL: filter 
# EQUAL: wave_res 
# EQUAL: lsf_fwhm 
# EQUAL: signal_to_noise 
# NOT EQUAL! sbin 

```

If a parameter is requested that is not compatible with the requested telescope, the code will issue a warning. A telescope object will be produced, but the code will assume a compatible answer by either selecting the closest available option, or replacing with the most common request.

```r
# Trying to make a MaNGA instrument with an incompatible field-of-view
manga_15 = telescope(type = "MaNGA", fov = 15)

# Warning message:
# In telescope(type = "MaNGA", fov = 15) :
# >>> WARNING! >>> 
# `fov` specified is not an option for telescope type `MaNGA`. 
# `fov` options include 12, 17, 22, 27, or 32. 
# Setting 'fov' to the nearest available option. 
# fov = 17

# Trying to make a MUSE instrument with an incompatible spatial resolution
muse_nfm = telescope(type = "MUSE", spatial_res=0.001)

# Warning message:
# In telescope(type = "MUSE", spatial_res = 0.001) :
# >>> WARNING! >>> 
# `spatial_res` specified is not an option for telescope type `MUSE`. 
# `spatial_res` options include 0.025 (NFM) or 0.2 (WFM). 
# Setting `spatial_res` to the default. 
# spatial_res = 0.2 

```

Or, if you request a `filter` whose wavelength range does not overlap with the wavelength range of the telescope, you will also receive an error:

```r
# Trying to build an instrument with incompatible wave-length coverage and filters
telescope(type="IFU", wave_range = c(1000,3000), filter = "g")

#  Error in telescope(type = "IFU", wave_range = c(1000, 3000), filter = "g") : 
#  Error: Requested filter will not overlap with the telescope wavelength range. 
#         Please select a different filter or extend your telescope wavelength range.
```

Further notes on filter choice and wavelength range compatibility can be found below. 

---

## Notes

### Filter Choice

Recent Updates
{: .label .label-purple } 

With the changes implemented in the upgrade to v2.3.8, cubes built with `method = "velocity"` and `mass_flag = FALSE` will have velocity maps weighted by the observed flux in a given filter band. Similarly, the `flux_image` produced will be the flux in a given filter. 

Below, we show the wavelength coverage for each of the available filters within SimSpin. The filter selected in this function will need to overlap with the wavelength range of the requested telescope. If they do not overlap, SimSpin will issue an error. 

| **Filter Name** | **Wavelength Range,** &#8491; | **Data Shape** | **Description**                                                                                    |
|-----------------|--------------------------------|----------------|----------------------------------------------------------------------------------------------------|
| `filt_u_SDSS`   | 2980 - 4130                    | (47, 2)        | The relative response of the **SDSS u-band** filter across the relevant wavelength range in Angstroms. |
| `filt_g_SDSS`   | 3630 - 5830                    | (89, 2)        | The relative response of the **SDSS g-band** filter across the relevant wavelength range in Angstroms. |
| `filt_r_SDSS`   | 5380 - 7230                    | (75, 2)        | The relative response of the **SDSS r-band** filter across the relevant wavelength range in Angstroms. |
| `filt_i_SDSS`   | 6430 - 8630                    | (89, 2)        | The relative response of the **SDSS i-band** filter across the relevant wavelength range in Angstroms. |
| `filt_z_SDSS`   | 7730 - 11230                   | (141, 2)       | The relative response of the **SDSS z-band** filter across the relevant wavelength range in Angstroms. |

---

Now that you've initiallised the telescope, you can move on to describing the strategy of your observing run...

[Next Step: observing_strategy](observing_strategy.markdown){: .btn .btn-purple }