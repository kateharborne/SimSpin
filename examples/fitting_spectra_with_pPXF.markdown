---
layout: default
title: Fitting spectra with pPXF
parent: Examples
nav_order: 3
last_modified_date: "Thu, 4th Aug 2022 14:08:00 AWST"
---

# Fitting stellar kinematics using pPXF

In order to produce observables from spectral data cubes, we must first process the outputs using external software. In this example, we demonstrate how an IFS cube produced using SimSpin can be fit using the industry standard, pPXF. 
{: .fs-5 .fw-300 .pb-2 }
---

When making an observation using an IFU, the observed spectra at each pixel must be fit to determine the observed kinematics. Most commonly, the penalised pixel fitting software, [pPXF](https://pypi.org/project/ppxf/) (Cappellari & Emsellem 2004) is used by IFU surveys to produce observable kinematics.

The spectral cubes built by SimSpin are directly compatible with this software. In the example below, we demonstrate how to run a spectral SimSpin cube through pPXF in order to measure the observed line-of-sight velocity distribution. As pPXF is a Python code, we will present some python functions for reading and fitting the SimSpin FITS files produced in the previous example [Working with FITS](/SimSpin/examples/working_with_FITS).

We break this down into a series of steps:
1. [Installing pPXF](#installing-ppxf)
2. [Reading in the observation](#reading-in-the-observation)
3. [Reading in the pPXF templates](#reading-in-the-ppxf-templates)
4. [Running the fit](#running-the-fit)

---

## Installing pPXF

First, let's ensure we have pPXF installed. This can be done simply using the command:

```
pip install ppxf
```

This should also download the necessary dependencies. To check that this has installed correctly, run:

```
pip show ppxf

# Name: ppxf
# Version: 8.1.0
# Summary: pPXF: Full Spectrum and SED Fitting of Galactic and Stellar Spectra
# Home-page: https://purl.org/cappellari/software
# Author: Michele Cappellari
# Author-email: michele.cappellari@physics.ox.ac.uk
# License: Other/Proprietary License
# Location: /path/to/lib/python3.9/site-packages
# Requires: astropy, matplotlib, numpy, scipy
# Required-by: 
```

This tells us that we have the software installed and the location at which the package has been installed. 
Once this is successful, we can begin the process of fitting a SimSpin observation. We will be using the spectral FITS created for the last example. That was done using the R code below:

```R
# Load a model...
simulation_file = system.file("extdata", "SimSpin_example_EAGLE", package = "SimSpin")
# ... use to build a SimSpin file. 
simspin_eagle  = make_simspin_file(filename = simulation_file,
                                   template = "EMILES",
                                   write_to_file = FALSE)

# Building a datacube with default spectral parameters 
cube   = build_datacube(simspin_file = simspin_eagle,
                        telescope = telescope(),
                        observing_strategy = observing_strategy(),
                        method = "spectral",
                        write_fits = T, 
                        output_location = "SimSpin_spectral_example_EAGLE.FITS",  
                        object_name = "SimSpin_example_EAGLE",
                        telescope_name = "SimSpin",
                        observer_name = "Anonymous",
                        split_save = F)
```

The `SimSpin_spectral_example_EAGLE.FITS` file produced by this function will be used in the fitting procedure. For ease running as a function with parameters, we have organised the pPXF fitting code with ArgParse, such that the arguments for each pPXF fit can be passed in from the command line. The full code for this automation can be found [here]() on GitHub, but we take you through each of the steps below.

Let's begin by importing the necessary packages and functions:

```python 
# Importing the maths and astronomy packages required for performing pPXF fits.

import numpy as np
from astropy.io import fits
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
from scipy import ndimage,signal
import glob
import os
import sys

from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

from ppxf.ppxf import ppxf, rebin
import ppxf.ppxf_util as util
from ppxf import miles_util as lib
```

Once these steps are done, we are ready to start fitting with pPXF!

---

## Reading in the observation

The methodology behind pPXF is that a series of template spectra are broadened by a LOSVD and iteratively fit to the observed spectrum to find the best fitting combination of parameters. There are several preparatory steps before we begin this process, however. Let's begin by declaring the paths and variables required for running the fit. We would like to fit for *V*, $$\sigma$$, $$h_3$$ and $$h_4$$ (i.e. 4 moments of the LOSVD).

```python
#### INPUT VARIABLES ####
c = 299792.458 # speed of light in km/s
template_loc = "path/to/lib/python3.9/site-packages/ppxf/miles_models/Eun1.30*.fits" # path to the template spectra used for pPXF fitting
FWHM_temp = 2.51 # Angstrom (resolution of the template EMILES spectra used for fitting)
file_input = "SimSpin_example_EAGLE.FITS" # spectral cube generated using SimSpin
n_mom = 4 # number of LOSVD moments to compute
############################
```

First, we will bring all of the observed spectra back to the rest-frame. This is done because our pPXF fitting template spectra are at rest-frame. This requires us to know the observational redshift of the spectra (as given in the `OB_TABLE` under `z`). In this example, the simulated galaxy has been put at z = 0.05, so we use this to de-shift the observed spectrum. 

It is also necessary to convolve the templates with a Gaussian such that the line-spread-functions (LSF) of the mock observation and the pPXF fitting templates match. As the cubes have been built from SimSpin template spectra themselves which have an underlying LSF intrinsically (before the addition of the observing telescope LSF), we need to take this into account. In this case, we used `EMILES` as our templates from which the SimSpin file was built. These templates have an intrinsic LSF = 2.51 A. 

```python
#---------------------------------#
# Step 1: Read in the observation #
#---------------------------------#

z = 0.05  # observed redshift of the model galaxy
telescope_FWHM = 3.2  # as specified in the telescope function
SimSpin_template_FWHM  = 2.51 # as dictated by the spectra from which the SimSpin file was built

# If the observing telescope had a requested LSF *less* than the templates themselves, we just use the reshifted SimSpin intrinsic template LSF as the value associated with the observed spectra. 
if telescope_FWHM <= SimSpin_template_FWHM:
    FWHM_gal = (SimSpin_template_FWHM * (1 + z)) # Angstroms (resolution for this "observed" template)
else:
    FWHM_gal = telescope_FWHM

hdu = fits.open(f"{file_input}")

cube     = hdu["DATA"].data    # accessing the spectral datacube within the file
head_obs = hdu["DATA"].header  #    and the associated axes labels
npart    = hdu["NPART"].data   # number of particles per pixel (used to describe the level of noise per pixel). 

print(f"Reading spectral cube of size {cube.shape}.")

# Reshape the cube into an array
npix = cube.shape[0]
nspec = cube.shape[1] * cube.shape[2]
spectrum = cube.reshape(npix, -1) # create an array of spectra [npix, (sbin_x * sbin_y)]
npart_flat = npart.reshape(-1)
noise_flat = 1/np.sqrt(npart_flat)

print(f"Re-shaping into array of shape {spectrum.shape}, with {spectrum.shape[0]} wavelengths and {spectrum.shape[1]} pixels.")

# Begin by pulling out the wavelength axis of the cube -------------------------------------------------------------
# These properties are necessary for logarithmically re-binning the spectra before
# fitting with pPXF.
wave_axis = head_obs['CRVAL3'] + (np.arange(0., (head_obs['NAXIS3']))*head_obs['CDELT3'])
wave_range = head_obs["CRVAL3"] + np.array([0., head_obs["CDELT3"]*(head_obs["NAXIS3"]-1)])
velscale = c * np.diff(np.log(wave_axis[-2:])) # velocity step
velscale = velscale[0]

print(f"The velocity scale of the `observed` spectrum is {velscale} km/s per pixel.")

# We also need to bring the observed spectrum back to the rest-frame ------------------------------------
if z >= 0.0:
    wave_range = wave_range / (1 + z)
    FWHM_gal = FWHM_gal / (1 + z)
    z = 0
    print(f"The wavelength range of the 'observed' galaxy (returned to rest-frame) is {wave_range} Angstrom.")
    print(f"The FWHM at this restframe is {FWHM_gal}.")

else:
    print(f"The wavelength range of the 'observed' galaxy is {wave_range} Angstrom.")
    print(f"The FWHM of the observation is {FWHM_gal}.")

```
With these parameters ready, we next need to read in the templates we'll use to fit the observed spectra. 


## Reading in the pPXF templates

