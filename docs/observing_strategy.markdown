---
layout: default
title: observing_strategy
parent: Documentation
nav_order: 4
last_modified_date: "Fri, 27 May 2022 15:57:00 AWST"
---

# Defining the properties of the object being observed

In order to build a data cube, we need to describe the conditions in which the observation has been made. 
{: .fs-5 .fw-300 .pb-2 }

Using the ``observing_strategy`` function, a single object is generated to describe the object specific properties such as the projected distance, the projected inclination, the pointing location and the level of seeing. This is used, along with the `telescope` object, to describe the specific observing conditions input to `build_datacube`. 
{: .fw-300 }

[See an example](#example){: .btn .btn-purple }
[See the source code](https://github.com/kateharborne/SimSpin/blob/d020398fb66274443bb2f70ea1fdd8346c4476ae/R/observing_strategy.R#L56){: .btn .btn-purple }

---

The following code shows the default parameters used in the `observing_strategy` function. Calling the function without specifying any input will produce an observing strategy object with the following properties:

```R
observing_strategy(dist_z       = 0.05,   # projected distance to galaxy model
                   inc_deg      = 70,     # projected angle, 0 deg (o) face on, 90 deg (-) edge on
                   twist_deg    = 0,      # projected angle, 0 deg (o) face on, 90 deg (|) edge on 
                   pointing_kpc = c(0,0), # pointing of telescope centre relative to galaxy centre
                   blur         = F)      # seeing conditions
```

---

## Input Parameters

| `dist_Mpc`, `dist_kpc_per_arcsec`, `dist_z` | The projected distance at which the galaxy is placed. Only one of these three input parameters should be defined at once. Each input describes the projected distance in different units (1) distance in mega-parsec (Mpc) (2) distance as a spatial scale in kilo-parsec per arcsecond (kpc/'') or (3) redshift distance.  |
| `inc_deg`  | The projected inclination of the observed galaxy relative to the z-axis - 0 deg places the galaxy face-on (o), 90 deg is edge-on aligned with the horizontal axis (-). Default is 70. |
| `twist_deg` | The projected inclination of the observed galaxy relative to the x-axis - 0 deg places the galaxy face-on (o), 90 deg is edge-on aligned with the vertical axis (\|). Default is 0.   |
| `pointing_deg`, `pointing_kpc`|  A numeric array `c(x,y)` specifying the position at which the field-of-view is centred given as a shift relative to the centre at (0,0) in units of (1) degrees or (2) kilo-parsecs. Only one of these two parameters should be supplied. Default is c(0,0) offset. |
| `blur` | Should the observation be convolved with a PSF to mimic the effects of seeing conditions? If `TRUE`, the parameters below will describe the size and shape of this convolution kernel. |
| `psf` | The shape of the point-spread function with which the cube is convolved to generate a seeing-convolved observation. Options include "Gaussian" and "Moffat". Only used if `blur = TRUE `. |
| `fwhm` | The size of the point-spread function kernel, as measured by the full-width half-maximum in arcseconds. Only used if `blur = TRUE `. |

---

## Output Value

The output of `observing_strategy` is a *List* element that will be stored as a variable to the environment. 

The list will contain either 5 or 7 elements, dependent on whether the observation will be blurred for seeing conditions or not:

1. `distance` - A "distance" object that summarises the different distance measures to the projected galaxy. This will include three elements, determined from the one input above:
    - `z` which describes the projected redshift distance to the simulation.
    - `Mpc` which describes the physcial luminosity distance in Mpc to the simulation.
    - `kpc_per_arcsec` which describes the angular distance on the sky given in kpc/arcseconds to the simulation.
1. `inc_deg` - The projected inclination of the galaxy about the horizontal axis in degrees.
1. `twist_deg` - The projected inclination (or "twist") of the galaxy about the vertical axis in degrees.
1. `pointing` - A "pointing" object that summarises the distance between the centre of the telescope field of view and the centre of the galaxy model. This includes two elements determined from the one input above:
    - `xy_deg` which describes the offset in x and y in units of degrees (for the given projected distance).
    - `xy_kpc` which describes the offset in x and y in units of kpc. 
1. `blur` - A boolean that describes whether or not seeing conditions will be applied during the mock cube building process. 
1. `fwhm` - (*Optional*) If `blur=T`, this value will describe the full-width half-maximum of the point-spread function (PSF) used to blur the resulting data cube. 
1. `psf` - (*Optional*) If `blur=T`, this value will name the sape of the point-spread function (PSF). Currently, this can be either "Gaussian" or "Moffat". 

---

## Example

Start by initiallising an observing strategy using the default parameters:

```R
strategy = observing_strategy()
```

We can inspect the object produced and see that the default parameters have been stored as a list.

```R
summary(strategy)
#           Length Class    Mode   
# distance  1      Distance S4     
# inc_deg   1      -none-   numeric
# twist_deg 1      -none-   numeric
# pointing  1      Pointing S4     
# blur      1      -none-   logical
``` 

To inspect individual properties of our observing strategy, use the named elements in one of two ways to achieve the same result:

```R
strategy$inc_deg
# [1] 70

strategy[["inc_deg"]]
# [1] 70
```

The `distance` and `pointing` parameters are listed as *S4 objects*. When we examine each of these parameters, they will be formatted to include all definitions of distance/pointing, for example:

```R
strategy$distance
# Distance to observed galaxy is:
#  227.4797 Mpc ---- physical distance,
#  1.00032 kpc/'' - angular scale on the sky,
#  z = 0.05 -------- projected redshift distance.

strategy$pointing
# Pointing relative to observed galaxy centre at (0,0) is:
#  (x,y):  (0,0) degrees -- 
#  (x,y):  (0,0) kpc ------ 

```

In order to access individual components of these parameters, such as the projected distance to the galaxy in units of Mpc, we can use the S4 methods shown below:

```R
z(strategy$distance)
# [1] 0.05

Mpc(strategy$distance)
# [1] 227.4797

kpc_per_arcsec(strategy$distance)
# [1] 1.00032

xy_deg(strategy$pointing)
# [1] 0 0

xy_kpc(strategy$pointing)
# [1] 0 0

```

---

Now that you've described the observing strategy, you can move on to building a mock observation...

[Next Step: build_datacube](build_datacube.markdown){: .btn .btn-purple }