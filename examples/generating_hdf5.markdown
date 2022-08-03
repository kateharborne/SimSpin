---
layout: default
title: HDF5 simulation format
parent: Examples
nav_order: 1
last_modified_date: "Tues, 14 June 2022 14:08:00 AWST"
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

For simulation outputs in HDF5 format, we have individual `Groups` for the metadata of the file, as well as each particle type. We use the `PartType` definitions from Gadget [(Springel, 2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract), with each contained group listed here:

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

If `PartType0` or `PartType4` are mising entirely from the file, it is assumed that there are no particles of this type within the simulation cut-out. If `header` is not provided, the code will return an error. 

See example files for specific EAGLE and Magneticum formatting below.

[EAGLE file](https://github.com/kateharborne/SimSpin/blob/master/inst/extdata/SimSpin_example_EAGLE.hdf5){: .btn  }
[Magneticum file](https://github.com/kateharborne/SimSpin/blob/master/inst/extdata/SimSpin_example_Magneticum.hdf5){: .btn  }
{: .pb-6 }

#### Header 
{: .text-purple-200 .fs-4 }

The header for each hydrodynamic simulation contains atrributes that record the information necessary for SimSpin to distinguish the input simulation type. It also contains information about the properties of the simulation from which the galaxy has been extracted. In particular, there are a number of attribute elements that are important to include within the `Header` group:

| **Attribute Name** | **Data Shape** | **Description**                                                                                                   |
|--------------------|----------------|-------------------------------------------------------------------------------------------------------------------|
| `BoxSize`          | 1              | The linear extent of the simulation cube in Mpc/h.                                                                |
| `Redshift`         | 1              | The redshift of the simulation at this snapshot.                                                                  |
| `HubbleParam`      | 1              | The Hubble Parameter (H0 / 100 km/s/Mpc).                                                                         |
| `MassTable`        | 6              | The mass of particles for a given particle type. If 0, see the `Mass` table within the respective PartType group. |
| `NumPart_Total`    | 6              | The number of particles of each particle type within the simulation.                                              |
| `RunLabel`         | 1              | Name of the simulation.                                                                                           |


Other properties are often included as default. These will be ignored by `make_simspin_file`, so can be left as is, or trimmed to reduce the HDF5 file size. Below we show examples of the header information for each of our example files within SimSpin.

<figure>
<img align="centre" src="assets/images/example_HDF5_header.png" width="750" height="403" />
</figure>
{: .pb-6 }

#### Attributes
{: .text-purple-200 .fs-4 }

Each property in the following PartType groups have associated attributes that define the units of that property. These allow you know convert between physical and co-moving coordinates for a given property. These expected attributes are used to convert the properties within each table into CGS units for SimSpin processing:

| **Attribute Name**    | **Data Shape** | **Description**                                                                                     |
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

This will then read simply within the current architecture of the code. 
{: .pb-6 }

#### PartType0
{: .text-purple-200 .fs-4 }

This group contains tables describing the gas particle/cell properties. Each of these tables should have the attributes described [here](#attributes).
For a mock observation to be constructed with SimSpin, the following tables are necessary for processing:

| **Table Name**                         | **Data Shape** | **Description**                                           |
|----------------------------------------|----------------|-----------------------------------------------------------|
| `Coordinates`                          | N x 3          | The (x,y,z) coordinates of each particle/cell.            |
| `Density`                              | N              | The gas density.                                          |
| `Mass`                                 | N              | Particle/cell mass.                                       |
| `ParticleIDs`                          | N              | Unique particle identification number.                    |
| `SmoothedElementAbundance/Carbon`[^1]  | N              | Mass fraction of Carbon for a given particle/cell mass.   |
| `SmoothedElementAbundance/Oxygen`[^1]  | N              | Mass fraction of Oxygen for a given particle/cell mass.   |
| `SmoothedElementAbundance/Hydrogen`[^1]| N              | Mass fraction of Hydrogen for a given particle/cell mass. |
| `SmoothedMetallicity`                  | N              | Mass fraction of all elements heavier than Helium.        |
| `SmoothingLength`                      | N              | SPH smoothing kernel radius (SPH simulations only)        |
| `StarFormationRate`                    | N              | Instantaneous star formation rate                         |
| `Temperature`                          | N              | Temperature                                               |
| `Velocity`                             | N x 3          | The (vx, vy, vz) coordinates of each particle/cell.       |

[^1]: this element of the HDF5 file is a group `SmoothedElementAbundance` containing a further three groups (`Carbon`, `Oxygen`, and `Hydrogen`) that contain data tables. 

where *N* denotes the number of particles within the simulation cut-out. 

Further properties can be passed into SimSpin without error, but in the case that one of these groups is missing the code will return an error.
{: .pb-6 }

#### PartType4
{: .text-purple-200 .fs-4 }

This group contains tables describing the stellar particle properties. Each of these tables should have the attributes described [here](#attributes).
For a mock observation to be constructed with SimSpin, the following tables are necessary for processing:

| **Table Name**                         | **Data Shape** | **Description**                                           |
|----------------------------------------|----------------|-----------------------------------------------------------|
| `Coordinates`                          | N x 3          | The (x,y,z) coordinates of each particle/cell.            |
| `InitialMass`                          | N              | The particle mass at formation time.                      |
| `Mass`                                 | N              | Particle/cell mass.                                       |
| `ParticleIDs`                          | N              | Unique particle identification number.                    |
| `SmoothedMetallicity`                  | N              | Mass fraction of all elements heavier than Helium.        |
| `StellarFormationTime`                 | N              | Expansion factor, a, when the star was born.              |
| `Velocity`                             | N x 3          | The (vx, vy, vz) coordinates of each particle/cell.       |

where *N* denotes the number of particles within the simulation cut-out. 

Further properties can be passed into SimSpin without error, but in the case that one of these groups is missing the code will return an error.
{: .pb-6 }

----

### N-body Simulations

---