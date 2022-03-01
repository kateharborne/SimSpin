---
layout: default
title: HDF5 simulation format
parent: Examples
nav_order: 1
last_modified_date: "Tues, 1 March 2022 14:08:00 AWST"
---

# Make your own HDF5 formatted simulation data 

`SimSpin` attempts to be agnostic to the type of simulation file input (as the aim is to produce a comparable data product for all simulations). 
{: .fs-5 .fw-300 .pb-2 }
---

Currently, `SimSpin` has ability to support:
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

For simulation outputs in HDF5 format, we have individual `Groups` for each particle type. We use the `PartType` definitions from Gadget [(Springel, 2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract):

* `PartType0` - gas,
* `PartType1` - dark matter, 
* `PartType2` - *N*-body disk, 
* `PartType3` - *N*-body bulge, 
* `PartType4` - stars, 
* `PartType5` - boundary/black holes. 
{: .lh-tight }

Within each of these groups, properties such as positions and velocities are stored for particles of that type. The specific properties contained depend on the type of simulation input. For example, a hydrodynamical simulation 