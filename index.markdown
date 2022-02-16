---
layout: default
title: Home
nav_order: 1
description: "SimSpin is a package for producing mock IFS observations of galaxy simulations."
permalink: /SimSpin
---

# SimSpin v2.1.1 
{: .fs-9 }

A package for producing mock observations of particle simulations
{: .fs-6 .fw-300 }

The purpose of the [SimSpin](https://github.com/kateharborne/SimSpin) package is to take a particle simulation of a galaxy and produce an integral field spectroscopy (IFS) observation.

This produces a **data cube** - i.e. spatial information in projection (*xy*) with spectral or kinematic information along the line-of-sight (*z*). 

A mock data cube can be produced using this package. 
This is a simple process comprised of four steps:

<img class="center" src="assets/images/SimSpin_methodology.png" /> 

We incorporate some of the limitations encountered by observers so that more consistent comparisons can be made between observations and theory.

From this data cube, "observables" can be measured using observational pipelines. 
This package, once installed, is fully documented and tested.

Another implementation of this code (SimSpin v1.1.3) written in Julia is also available at [SimSpin.jl](https://github.com/kateharborne/SimSpin.jl) developed by [Gerry Gralton](https://github.com/gerrygralton). 

The purpose of [this guide](https://kateharborne.github.io/SimSpin/) is to act as a reference for users. 
It outlines the intended functionality and processes available in the `SimSpin` code. 
This is both a form of documentation and a collection of examples. 
For a basic walk-through of installing `SimSpin` and running a simple observation, please see the links on the left. 
