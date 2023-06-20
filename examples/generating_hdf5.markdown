---
layout: default
title: HDF5 simulation format
parent: Examples
nav_order: 1
last_modified_date: "Tue, 20 June 2023 11:11:00 AWST"
---

# Make your own HDF5 formatted simulation data 

SimSpin attempts to be agnostic to the type of simulation file input (as the aim is to produce a comparable data product for all simulations). 
{: .fs-5 .fw-300 .pb-2 }
---

Currently, SimSpin has ability to support:
* Gadget2 simulation outputs written in binary class 1.
* Gadget2 simulation outputs written in HDF5.
* Galaxies cut from the EAGLE simulation and written to HDF5 using the [`read_eagle`](https://github.com/jchelly/read_eagle) routines.
* Galaxies cut from the Magneticum simulation and written to HDF5 using the [`g3read`](https://github.com/aragagnin/g3read) routines.
{: .lh-tight }

The are two possible options available if you would like to use **a different simulation output**. Either raise an issue on GitHub using the button below, or read on to learn how to write your own HDF5 file. 

[Raise an issue](https://github.com/kateharborne/SimSpin/issues/new/choose){: .btn .btn-purple }
[Write your own HDF5](#write-your-own){: .btn .btn-purple }

---

## Write your own

Hierarchical Data Format version 5 files (or HDF5 files) are designed to store large amounts of data in a structured way. 
For astronomers, the format is very similar to FITS files.

A single file contains `Groups` and `Datasets`, as well as `Metadata` information that records details about each group or dataset. 
A `Group` is like a folder or directory. It can contain `Datasets` of arbitrary size, or other `Groups`. 

For simulation outputs in HDF5 format, we have individual `Groups` for each particle type, as well as a `Header` that contains the metadata summary for the full file. For each particle type, we use the `PartType` definitions from Gadget [(Springel, 2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract), with each contained group listed here:

* `Header`
* `PartType0` - gas,
* `PartType1` - dark matter, 
* `PartType2` - *N*-body disk, 
* `PartType3` - *N*-body bulge, 
* `PartType4` - stars, 
* `PartType5` - boundary/black holes. 
{: .lh-tight }

Within each of these groups, properties such as positions and velocities are stored for particles of that type. The specific properties contained depend on the type of simulation input. We provide examples here for a galaxy cut out from a hydrodynamical simulation and from an N-body simulation.

[Hydrodynamic Example](#hydrodynamical-simulations){: .btn  }
[N-body Example](#n-body-simulations){: .btn  }

----

### Hydrodynamical Simulations

Here we discuss the tables included for each particle type within an input hydrodynamical simulation HDF5 file.
{: .fs-5 .fw-300 .pb-2 }

The following simulations are already supported for automatic read in, with example files available within the R-package.
<img align="centre" src="assets/images/supported_hydro_sims.png" width="750" height="237" />

These files must contain the following groups:
* [`Header`](#header)
* [`PartType0`](#parttype0) (i.e. gas particles/cells)
* [`PartType4`](#parttype4) (i.e. stellar particles)
{: .lh-tight }

Any additional groups (such as those describing the number of dark matter and black hole properties) will not be used within SimSpin at this time and as a result will be ignored if provided within the input simulation file. 

If `PartType0` or `PartType4` are mising entirely from the file, it is assumed that there are no particles of this type within the simulation cut-out. If a `Header` is not provided, the code will return an error. 

See example files for specific EAGLE, Magneticum and HorizonAGN formatting below.

[EAGLE file](https://github.com/kateharborne/SimSpin/blob/master/inst/extdata/SimSpin_example_EAGLE.hdf5){: .btn  }
[Magneticum file](https://github.com/kateharborne/SimSpin/blob/master/inst/extdata/SimSpin_example_Magneticum.hdf5){: .btn  }
[HorizonAGN file](https://github.com/kateharborne/SimSpin/blob/main/inst/extdata/SimSpin_example_HorizonAGN.hdf5){: .btn }
{: .pb-6 }

#### Header 
{: .text-purple-200 .fs-4 }

The header for each hydrodynamic simulation contains atrributes that record the information necessary for SimSpin to distinguish the input simulation type. It also contains information about the properties of the simulation from which the galaxy has been extracted. In particular, there are a number of attribute elements that are important to include within the `Header` group:

| **Attribute Name** | **Shape** | **Description**                                                                                                   |
|--------------------|----------------|-------------------------------------------------------------------------------------------------------------------|
| `BoxSize`          | 1              | The linear extent of the simulation cube in Mpc/h.                                                                |
| `Redshift`         | 1              | The redshift of the simulation at this snapshot.                                                                  |
| `HubbleParam`      | 1              | The Hubble Parameter (H0 / 100 km/s/Mpc).                                                                         |
| `MassTable`        | 6              | The mass of particles for a given particle type. If 0, see the `Mass` table within the respective PartType group. |
| `NumPart_Total`    | 6              | The number of particles of each particle type within the simulation file.                                              |
| `RunLabel`<sup>1</sup> | 1              | Name of the simulation.                                                                                           |

Other properties are often included as default. These will be ignored by `make_simspin_file`, so can be left as is, or trimmed to reduce the HDF5 file size. Below we show examples of the header information for each of our example files within SimSpin.

<sup>1</sup> If `RunLabel` is missing from the header of the file, SimSpin will assume that the input file is an N-body simulation and look for stellar components in the bulge and disk Groups (`PartType2` and `PartType3`). Ensure that a value is listed under this attribute if the underlying simulation is a hydrodynamical simulation with stellar particles listed in `PartType4`. 

<figure>
<img align="centre" src="assets/images/example_HDF5_header.png" width="750" height="403" />
</figure>
{: .pb-6 }

#### Attributes
{: .text-purple-200 .fs-4 }

Each property in the following PartType groups have associated attributes that define the units of that property. These allow you convert from co-moving properties output by the simulation into the physical cordinates in which we would like to work for individual galaxies. The following attributes are expected and used to convert the properties into CGS (centimeter-grams-seconds) units for SimSpin processing:

| **Attribute Name**    | **Shape** | **Description**                                                                                     |
|-----------------------|----------------|-----------------------------------------------------------------------------------------------------|
| `CGSConversionFactor` | 1              | The constant parameter required to convert the value from internal simulation units to CGS          |
| `aexp-scale-exponent` | 1              | The expansion factor (*a*) exponent necessary to convert between co-moving and physical coordinates |
| `h-scale-exponent`    | 1              | The Hubble parameter (*h*) exponent necessary to convert between co-moving and physical coordinates |

Using the following formula, we can then convert the values within the input file into CGS:

<img align="centre" src="assets/images/attributes_equation.png" width="750" height="36" />

where *a* is the scale factor at the snapshot (i.e. *a = 1/(1+z)), *h* is the Hubble parameter, the *simulation parameter* is the value returned from that table for a given particle, and the *physical parameter* is the output value in CGS units. This is done automatically in the `make_simspin_file` function and so it is important that these attributes are listed against each table correctly.

If your parameters are already in **physical units**, these values can simply be set to:
* `CGSConversionFactor` = Constant required to convert into CGS units
* `aexp-scale-exponent` = 0
* `h-scale-exponent`    = 0 
{: .lh-tight }

*It has been noticed that some IllustrisTNG has different names for these conversion factors. SimSpin will automatically accept either the default outputs from IllustrisTNG (`a_scaling`, `h_scaling` and `to_cgs`) or the named parameters above. The `to_cgs` factor output by IllustrisTNG is sometimes 0 within their HDF5 outputs. This will result in all values read in by SimSpin to become 0. To avoid this, SimSpin will take `to_cgs` = 0 and assume this indicates that the values are already in CGS units (i.e. CGSConversionFactor = 1).*

This will then read simply within the current architecture of the code. 
{: .pb-6 }

#### PartType0
{: .text-purple-200 .fs-4 }

This group contains a series of Datasets describing the gas particle/cell properties. 
Each of these Datasets should have the attributes described [here](#attributes).

For a mock observation to be constructed with SimSpin, the following Datasets will be used for processing.
These Datasets will be in table format with shape described relative to the number of gas particles in that simulation file, *N<sub>gas</sub>*. 

Necessary Datasets for functionality are indicated by the &#9733; symbol. 
For Datasets without a &#9733;, some may be required for *cell*-based codes and others for *particle*-based codes. 
These are further elaborated in the connected footnotes. 

| **Table Name**                         | **Shape** | **Description**                                           |
|----------------------------------------|-----------|-----------------------------------------------------------|
| `Coordinates`  &#9733;                 | N<sub>gas</sub> x 3 | The (x,y,z) coordinates of each particle/cell.            |
| `Density`      &#9733;                 | N<sub>gas</sub>     | The gas density.                                          |
| `Mass`         &#9733;                 | N<sub>gas</sub>     | Particle/cell mass.                                       |
| `ParticleIDs`  &#9733;                 | N<sub>gas</sub>     | Unique particle identification number.                    |
| `ElementAbundance/Carbon`<sup>1,2</sup>| N<sub>gas</sub>     | Mass fraction of Carbon for a given particle/cell mass.   |
| `ElementAbundance/Oxygen`<sup>1,2</sup>| N<sub>gas</sub>     | Mass fraction of Oxygen for a given particle/cell mass.   |
| `ElementAbundance/Hydrogen`<sup>1,2</sup>| N<sub>gas</sub>     | Mass fraction of Hydrogen for a given particle/cell mass.          |
| `Metallicity`<sup>2</sup> &#9733;      | N<sub>gas</sub>     | Mass fraction of all elements heavier than Helium.        |
| `SmoothingLength` <sup>3</sup>         | N<sub>gas</sub>     | SPH smoothing kernel radius (SPH simulations only)        |
| `StarFormationRate` &#9733;            | N<sub>gas</sub>     | Instantaneous star formation rate                         |
| `Temperature`<sup>4</sup>              | N<sub>gas</sub>     | Temperature.                                              |
| `Pressure`<sup>4, 5</sup>                 | N<sub>gas</sub>     | Gas pressure.                                              |
| `InternalEnergy`<sup>4, 5</sup>           | N<sub>gas</sub>     | Internal thermal energy per unit mass for each gas *cell*.|
| `ElectronAbundance`<sup>4</sup>        | N<sub>gas</sub>     | Fractional electron number density with respect to the total hydrogen number density. |
| `Velocity`  &#9733;                    | N<sub>gas</sub> x 3 | The (vx, vy, vz) coordinates of each particle/cell.       |

Further properties can be passed into SimSpin without error, but in the case that one of these necessary groups is missing the code will return an error.
{: .pb-6 }

<sup>1</sup> This element of the HDF5 file is a group `ElementAbundance` containing a further three groups (`Carbon`, `Oxygen`, and `Hydrogen`) that contain data tables. 

<sup>2</sup> Some simulations use different default names for individual datasets. For EAGLE outputs, these Group names will be prefixed with the name "Smoothed" to indicate that the parameters have been smoothed across the SPH smoothing kernel volume, as this is suggested to be the more appropriate output value to use when working with EAGLE outputs as in [Wiersma et al. 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.399..574W/abstract). In IllustrisTNG, these names will be prefixed with "GFM_". SimSpin is equipped to interpret these fields as intended such that you can choose to use the default field names rather than the ones listed here. 

<sup>3</sup> For *cell*-based codes, the `SmoothingLength` value will not be present in output files. To approximate the variable smoothing of cells across the pixelated image, SimSpin will compute an equivalent `SmoothingLength` for cell-based codes using the `Mass` and `Density` fields. 

<sup>4</sup> For *cell*-based gas simulations, the temperature for each cell can be approximated using the `InternalEnergy` and `ElectronAbundance` properties. In this case, the `Temperature` Group is not required as long as the `InternalEnergy` and `ElectronAbundance` Groups are present. Equally, the  `ElectronAbundance` Group is not required if `Temperature` is already present. We note here that **caution should be taken when using output gas temperatures for the simulations** as there are cooling floors employed in most hydrodynamical simulations that cause an artificial lower temperture limit. This is particularly important for low density and low temperature gas, as well as star-forming gas.  

<sup>5</sup> Gas within a simulation has a thermal contribution to the observable velocity dispersion that must be included to give a realistic representation of the true gas dispersion. To quantify the thermal contribution to this, we use $\sigma_{\text{thermal}}^2 =  P / \rho =  u (1 - \gamma)$, where the thermal dispersion contribution can be computed using either the pressure and denisty, or the internal energy. For particles with temperature below the cooling floor of the simulation, this vlaue is fixed to an isotropic thermal value of 11 km/s. 
 

#### PartType4
{: .text-purple-200 .fs-4 }

This group contains Datasets describing the stellar particle properties. 
Each of these Datasets should have the attributes described [here](#attributes).

These Datasets will be in table format with shape described relative to the number of gas particles in that simulation file, *N<sub>star</sub>*. 
All listed Datasets in the table below are required for functionality. 

| **Table Name**         | **Shape** | **Description**                                           |
|------------------------|----------------|-----------------------------------------------------------|
| `Coordinates`          | N<sub>star</sub> x 3  | The (x,y,z) coordinates of each particle/cell.            |
| `InitialMass`          | N<sub>star</sub>      | The particle mass at formation time.                      |
| `Mass`                 | N<sub>star</sub>      | Particle/cell mass.                                       |
| `ParticleIDs`          | N<sub>star</sub>      | Unique particle identification number.                    |
| `Metallicity`<sup>1</sup>          | N<sub>star</sub>      | Mass fraction of all elements heavier than Helium.        |
| `StellarFormationTime`<sup>2</sup> | N<sub>star</sub>      | Expansion factor, a, when the star was born.              |
| `Velocity`             | N<sub>star</sub> x 3  | The (vx, vy, vz) coordinates of each particle/cell.       |

Further properties can be passed into SimSpin without error, but in the case that one of these groups is missing the code will return an error.
{: .pb-6 }

<sup>1</sup> The Magneticum simulations give `Metallicity` as a table of raw metal masses rather than a metallicity in the form *Z = M<sub>Z</sub> / M<sub>total</sub>*. This is accounted for when read into SimSpin. 

<sup>2</sup> In the case that `StellarFormationTime` is negative (as is used to denote stellar wind particles in IllustrisTNG), SimSpin will remove these rows in all Dataset tables.

----

### N-body Simulations

SimSpin was originally designed to work with N-body models produced using Gadget2. 
{: .fs-5 .fw-300 .pb-2 }

When using binary or HDF5 N-body outputs from Gadget2, SimSpin will work with the raw output simulation files without modification. 
We list the parameters recovered for each particle flavour below in case a user would like to apply SimSpin to non-Gadget N-body models.
This is only described for HDF5 files, but for interested parties, the Gadget2 binary format is quite clearly deconstructed in the lines of code [here](https://github.com/kateharborne/SimSpin/blob/ef171e40626edacc350c1e119cc6f807c6335529/R/utilities.R#L67).

#### Header 
{: .text-purple-200 .fs-4 }

The Header block for output Gadget2 N-body models contain the following attributes, as listed in Table 4 of the [User's guide for GADGET-2](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf):

| **Attribute Name** | **Shape** | **Description**                                                                                                   |
|--------------------|-----------|-------------------------------------------------------------------------------------------------------------------|
| `BoxSize`          | 1         | Gives the box size if periodic boundary conditions are used.                                                      |
| `Flag_Cooling`     | 1         | Flag for cooling (unused in public version of Gadget-2). |
| `Flag_Entropy_ICs` | 1         | Flags that initial conditions contain entropy instead of internal energy in the u block. |
| `Flag_Feedback`    | 1         | Flag for feedback (unused in public version of Gadget-2). |
| `Flag_Metals`      | 1         | Flag for metallicity values (unused in public version of Gadget-2). |
| `Flag_Sfr`         | 1         | Flag for star formation (unused in public version of Gadget-2). |
| `Flag_StellarAge`  | 1         | Flag for creation times of stars (unused in public version of Gadget-2). |
| `HubbleParam`      | 1         | Gives "h", the Hubble constant in units of 100 km/s/Mpc (for cosmological integrations).                          |
| `MassTable`        | 6         | The mass of particles for a given particle type. If 0, see the `Mass` table within the respective PartType group. |
| `NumFilesPerSnapshot`| 1       | Number of files in each snapshot. |
| `NumPart_ThisFile` | 6         | The number of particles of each type in the present file. |
| `NumPart_Total`    | 6         | The number of particles of each particle type within the simulation (which may be split across multiple files).   |
| `NumPart_Total_HighWord` | 6   | For simulations that use more than 2<sup>32</sup> particles, these fields hold the most significant word of 64-bit total particle numbers. Otherwise, zero. |
| `Omega0`           | 1         | Matter density at z = 0 in units of the critical density (only relevant in cosmological integrations). |
| `OmegaLambda`      | 1         | Vacuum energy denisty at z = 0 in units of the criical density. |
| `Redshift`         | 1         | The redshift of the simulation at this snapshot.                |
| `Time`             | 1         | Time of output in Gyr, or expansion factor for cosmological simulations. |

Note that the `RunLabel` attribute is missing in this Header. 
This is used as an indication to SimSpin that the input HDF5 file is from an N-body code, rather than a hydrodynamical code. 
SimSpin then looks for Groups under the heading `PartType2` and `PartType3` for the stellar component of the galaxy. 

<figure>
<img align="centre" src="assets/images/nbody_example_HDF5_header.png" width="299" height="321" />
</figure>

#### PartType2
{: .text-purple-200 .fs-4 }

When run in N-body mode, the `PartType2` Group is reserved for **disk** particles, i.e. particles initialised with disk-like distribution functions. 

Within SimSpin, we assume that these particles trace a species of stars and will combine any particles in this group with any in `PartType3` bulge stars to form the stellar component table for analysis. 

Metallicities and ages are assigned to N-body disk particles in the [make_simspin_file](/SimSpin/docs/make_simspin_file.markdown) function seperately to the bulge particles.  Hence, the necessary Datasets within the Group only include:

| **Table Name**         | **Shape** | **Description**                                           |
|------------------------|----------------|-----------------------------------------------------------|
| `Coordinates`          | N<sub>disk</sub> x 3  | The (x,y,z) coordinates of each particle/cell.            |
| `Mass`                 | N<sub>disk</sub>      | Particle/cell mass.                                       |
| `ParticleIDs`          | N<sub>disk</sub>      | Unique particle identification number.                    |
| `Velocity`             | N<sub>disk</sub> x 3  | The (vx, vy, vz) coordinates of each particle/cell.       |

Further properties can be passed into SimSpin without error, but in the case that one of these groups is missing the code will return an error.
{: .pb-6 }

#### PartType3
{: .text-purple-200 .fs-4 }

When run in N-body mode, the `PartType3` Group is reserved for **bulge** particles, i.e. particles initialised with bulge-like distribution functions. 

Within SimSpin, we assume that these particles trace a species of stars and will combine any particles in this group with any in `PartType2` disk stars to form the stellar component table for analysis. 

Metallicities and ages are assigned to N-body bulge particles in the [make_simspin_file](/SimSpin/docs/make_simspin_file.markdown) function seperately to the disk particles.  Hence, the necessary Datasets within the Group include only: 

| **Table Name**         | **Shape** | **Description**                                           |
|------------------------|----------------|-----------------------------------------------------------|
| `Coordinates`          | N<sub>bulge</sub> x 3  | The (x,y,z) coordinates of each particle/cell.            |
| `Mass`                 | N<sub>bulge</sub>      | Particle/cell mass.                                       |
| `ParticleIDs`          | N<sub>bulge</sub>      | Unique particle identification number.                    |
| `Velocity`             | N<sub>bulge</sub> x 3  | The (vx, vy, vz) coordinates of each particle/cell.       |

Further properties can be passed into SimSpin without error, but in the case that one of these groups is missing the code will return an error.
{: .pb-6 }

---

This concludes the description of input file types supported in SimSpin at the current time. 
Check back here in the future for further simulation support developments.

If you would like us to formally integrate your simulation format into SimSpin, please get in contact by raising an issue:

[Raise an issue!](https://github.com/kateharborne/SimSpin/issues/new/choose){: .btn .btn-purple }