---
layout: default
title: build_datacube
parent: Documentation
nav_order: 4
last_modified_date: "Wed, 16 February 2022 15:57:00 AWST"
---

# Constructing your data cube

`build_datacube()`


## Input Parameters

| **type**              | A character string describing the type of telescope you wish to create. Default is `"IFU"`, in which case the telescope properties are described by the input parameters to this function. Other input types include: `"SAMI"`, `"MaNGA"`, `"CALIFA"` and `"MUSE"`. With these, the telescope properties are set to the default parameters specified by the respective instrument.    |
| **method**            |  A character string that describes whether to generate a spectral data cube or a kinematic data cube. Options include `"spectral"` and `"velocity"` to produce these for the stellar component of the simulation. If wishing to analyse the kinematics of the gas, you can also specify `"gas"` or `"sf_gas"` (i.e. star forming gas) to produce a kinematic data cube of this component alone.   |
