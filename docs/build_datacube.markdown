---
layout: default
title: build_datacube
parent: Documentation
nav_order: 5
last_modified_date: "Tue, 7 June 2022 15:57:00 AWST"
---

# Constructing your data cube

Recent Updates
{: .label .label-purple } 

Using the `build_datacube` function, we can now make a mock observation of a simulated galaxy. 
{: .fs-5 .fw-300 .pb-2 }

The output of this function will be a mock IFS data cube. The details of the observation will vary depending on the requested `method`, `telescope` and `observing_strategy`, but the format will be consistent. In all cases, the output can be written to a FITS file.
{: .fw-300 }

*As of version 2.2.0, the `method` input parameter has been moved directly to the [`build_datacube`](build_datacube.markdown) function. For backwards compatibility, this parameter can still be specified within [`telescope`](telescope.markdown), but a warning will be issued.*
{: .fw-150 }

---


The following code shows the default parameters used in the `build_datacube` function. Calling the function without specifying any input will produce an observation with the following properties:

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
               mass_flag = F)                    # mass or luminosity-weighted? (if method = "velocity").

```

{: pb-4 }

[Input parameters](#input-parameters){: .btn .btn-purple }
[Output parameters](#output-parameters){: .btn .btn-purple }
[See an example](#example){: .btn .btn-purple }
[See the source code](https://github.com/kateharborne/SimSpin/blob/bb371d9e4d981d1ebaba3aa07978bb61a2d434f0/R/build_datacube.R#L77){: .btn .btn-purple }

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
## Output parameters

The output of `build_datacube` is either:
 -  a *List* element that will be stored as a variable to the environment, or
 -  a *FITS file* that is written to the specified `output_location`. 

Here, we give details about the *List* written to the environment, but the same properties are stored within the output FITS file. An example of how to extract the necessary information from the saved FITS can be found [here](/SimSpin/examples/working_with_FITS).


The list will always contain 4 elements, though the contents of this list will change dependent on the mode in which the observation has been constructed. These elements are numbered below:
{: .fw-300 }

1.  `spectral_cube` or `velocity_cube` - A three-dimensional numeric array with spatial coordinates in x-y and either wavelength (in the case of a `spectral_cube`) or velocity (in the case of a `velocity_cube`) along the z-axis. 

1.  `observation` - A *List* element that summarises the methods used to construct the observation. 

    <details>
        <summary> <i>Expand to see list details</i> </summary>

        In the elements below, <sup>1</sup> denotes output included only in `method = "spectral"`, and <sup>2</sup> denotes output included only in `method = "velocity"`.
        
        <table class="tg">
        <thead>
        <tr>
            <th class="tg-0pky"><span style="font-weight:bold">List element</span></th>
            <th class="tg-0pky"><span style="font-weight:bold">Type</span></th>
            <th class="tg-0pky"><span style="font-weight:bold">Description</span></th>
        </tr>
        </thead>
        <tbody>
        <tr>
            <td class="tg-0pky">ang_size</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">scale at given distance in kpc/arcsec</td>
        </tr>
        <tr>
            <td class="tg-0pky">aperture_region</td>
            <td class="tg-0pky">bool</td>
            <td class="tg-0pky">pixels within the aperture</td>
        </tr>
        <tr>
            <td class="tg-0pky">aperture_shape</td>
            <td class="tg-0pky">str</td>
            <td class="tg-0pky">shape of aperture</td>
        </tr>
        <tr>
            <td class="tg-0pky">aperture_size</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">field of view diameter width in kpc</td>
        </tr>
        <tr>
            <td class="tg-0pky">date</td>
            <td class="tg-0pky">str</td>
            <td class="tg-0pky">date and time of mock observation</td>
        </tr>
        <tr>
            <td class="tg-0pky">fov</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">field of view diameter in arcsec</td>
        </tr>
        <tr>
            <td class="tg-0pky">filter</td>
            <td class="tg-0pky">str</td>
            <td class="tg-0pky">filter name</td>
        </tr>
        <tr>
            <td class="tg-0pky">inc_deg</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">inclination of object in degrees about the horizontal </td>
        </tr>
        <tr>
            <td class="tg-0pky">inc_rad</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">inclination of object in radians about the horizontal</td>
        </tr>
        <tr>
            <td class="tg-0pky">LSF_conv<sup>1</sup></td>
            <td class="tg-0pky">bool</td>
            <td class="tg-0pky">has line spread function convolution been applied?</td>
        </tr>
        <tr>
            <td class="tg-0pky">lsf_fwhm</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">full-width half-maximum of line-spread function in Angstrom</td>
        </tr>
        <tr>
            <td class="tg-0pky">lsf_sigma<sup>1</sup></td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">sigma width of line-spread function in Angstrom</td>
        </tr>
        <tr>
            <td class="tg-0pky">lum_dist</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">distance to object in Mpc</td>
        </tr>
        <tr>
            <td class="tg-0pky">method</td>
            <td class="tg-0pky">str</td>
            <td class="tg-0pky">name of observing method employed</td>
        </tr>
        <tr>
            <td class="tg-0pky">origin</td>
            <td class="tg-0pky">str</td>
            <td class="tg-0pky">version of SimSpin used for observing</td>
        </tr>
        <tr>
            <td class="tg-0pky">pointing_kpc</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">x-y position of field of view centre relative to object centre in units of kpc</td>
        </tr>
        <tr>
            <td class="tg-0pky">pointing_deg</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">x-y position of field of view centre relative to object centre in units of degrees</td>
        </tr>
        <tr>
            <td class="tg-0pky">psf_fwhm</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the full-width half-maximum of the point spread function kernel in arcsec</td>
        </tr>
        <tr>
            <td class="tg-0pky">sbin</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the number of spatial pixels across the diameter of the field of view</td>
        </tr>
        <tr>
            <td class="tg-0pky">sbin_seq</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the spatial bin centres in kpc</td>
        </tr>
        <tr>
            <td class="tg-0pky">sbin_size</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the size of each pixel in kpc</td>
        </tr>
        <tr>
            <td class="tg-0pky">spatial_res</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the size of each pixel in arcsec</td>
        </tr>
        <tr>
            <td class="tg-0pky">signal_to_noise</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the signal-to-noise ratio for observed spectrum</td>
        </tr>
        <tr>
            <td class="tg-0pky">twist_deg</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">inclination of object in degrees about the vertical</td>
        </tr>
        <tr>
            <td class="tg-0pky">twist_rad</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">inclination of object in radians about the vertical</td>
        </tr>
        <tr>
            <td class="tg-0pky">wave_bin</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the number of wavelength bins for a given telescope</td>
        </tr>
        <tr>
            <td class="tg-0pky">wave_centre</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the central wavelength for a given telescope in Angstrom</td>
        </tr>
        <tr>
            <td class="tg-0pky">wave_res</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the width of each wavelength bin in Angstrom</td>
        </tr>
        <tr>
            <td class="tg-0pky">wave_seq</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the wavelength bin centres in Angstrom</td>
        </tr>
        <tr>
            <td class="tg-0pky">wave_edges </td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the wavelength bin edges in Angstrom</td>
        </tr>
        <tr>
            <td class="tg-0pky">vbin<sup>2</sup></td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the number of velocity bins for a given telescope resolution</td>
        </tr>
        <tr>
            <td class="tg-0pky">vbin_seq<sup>2</sup></td>
            <td class="tg-0pky">num </td>
            <td class="tg-0pky">the velocity bin centres in km/s</td>
        </tr>
        <tr>
            <td class="tg-0pky">vbin_edges<sup>2</sup></td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the velocity bin edges in km/s</td>
        </tr>
        <tr>
            <td class="tg-0pky">vbin_size</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the size of each velocity bin in km/s</td>
        </tr>
        <tr>
            <td class="tg-0pky">vbin_error</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the velocity uncertainty given the telescope LSF in km/s</td>
        </tr>
        <tr>
            <td class="tg-0pky">z</td>
            <td class="tg-0pky">num</td>
            <td class="tg-0pky">the redshift distance of the object observed</td>
        </tr>
        </tbody>
        </table>
        
    </details>

    {: pb-3 }


1.  `raw_images` - This will be a list of a variable number of images.
    -   In `method = "spectral"` or `"velocity"`, there will be 6 raw images listed:

        | `flux_image` or `mass_image`  | a 2D numeric array giving the total flux in the requested band described by `filter` or the total mass within a given pixel. |
        | `velocity_image`   | a 2D numeric array giving the particle mass-weighted mean LOS velocity in a given pixel. |
        | `dispersion_image` | a 2D numeric array giving the particle mass-weighted standard deviation of the LOS velocity in a given pixel. |
        | `age_image` | a 2D numeric array giving the mass-weighted particle mean age in a given pixel. |
        | `metallicity_image` | a 2D numeric array giving the mass-weighted particle mean metallicity in a given pixel. |
        | `particle_image` | a 2D numeric array giving the number of particles per pixel. |

    -   In `method = "gas"` or `"sf_gas"`, there will be 7 raw images listed:

        | `mass_image`  | a 2D numeric array giving the total mass within a given pixel. |
        | `velocity_image`   | a 2D numeric array giving the particle mass-weighted mean LOS velocity in a given pixel. |
        | `dispersion_image` | a 2D numeric array giving the particle mass-weighted standard deviation of the LOS velocity in a given pixel. |
        | `SFR_image` | a 2D numeric array giving the mass-weighted particle mean star formation rate (SFR) in a given pixel. |
        | `metallicity_image` | a 2D numeric array giving the log10 mean gas metallicity in a given pixel, i.e. [log10(Z / 0.0127)]. |
        | `OH_image` | a 2D numeric array giving the log10 ratio of oxygen to hydrogen abundance in a given pixel, i.e. [log10(O/H) + 12]. |
        | `particle_image` | a 2D numeric array giving the number of particles per pixel. |

1.  `observed_images` - This will be a list of a variable number of images.
    -   In `method = "spectral"`, the returned element will be `NULL`. This is because observed images will need to be fitted using external software such as pPXF. 

    -   In all other methods, (`method = "velocity"`, `"gas"` or `"sf_gas"`), the returned element will be a list of length 5 containing:

        | `flux_image` or `mass_image`  | a 2D numeric array giving the total flux in the requested band described by `filter` or the total mass within a given pixel (always returned for gas observations). |
        | `velocity_image`   | a 2D numeric array giving the observed mean LOS velocity in a given pixel. |
        | `dispersion_image` | a 2D numeric array giving the observed standard deviation of the LOS velocity in a given pixel. |
        | `h3_image` | a 2D numeric array giving the observed h3 higher-order moment in a given pixel. |
        | `h4_image` | a 2D numeric array giving the observed h4 higher-order moment in a given pixel. |
{: .bg-grey-lt-000 }
---

## Example 