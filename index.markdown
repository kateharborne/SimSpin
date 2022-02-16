---
layout: default
title: Home
nav_order: 1
description: "SimSpin is a package for producing mock IFS observations of galaxy simulations."
permalink: /
last_modified_date: "Wed, 16 February 2022 13:57:00 AWST"
---

# SimSpin v2.1.1 
{: .fs-9 }

A package for producing mock observations of particle simulations
{: .fs-6 .fw-300 .mb-5 }

[View on GitHub](https://github.com/kateharborne/SimSpin){: .btn .btn-purple }
[Download Citation](https://github.com/kateharborne/SimSpin/blob/master/CITATION.cff){: .btn .btn-purple }
[View on ADS](https://ui.adsabs.harvard.edu/abs/2019ascl.soft03006H/abstract){: .btn .btn-purple }
{: .mb-1 }

<img align="right" src="assets/images/logo.png" width="175" height="175" />
{: .pl-3 .pb-1 } 

[`SimSpin`](https://github.com/kateharborne/SimSpin) is an R-package designed to take a simulation of a galaxy and produce a "mock" integral field spectroscopy (IFS) observation.

This produces a **data cube** - i.e. spatial information in projection (*xy*) with spectral or kinematic information along the line-of-sight (*z*). 

A mock data cube can be produced using this package. 
This is a simple process comprised of three steps:

  1. Read in your particle data and produce the relevant spectra using the `make_simspin_file` function.
  1. Setup the observation by defining your `telescope` and `observing_strategy`.
  1. Build your data cube using the `build_datacube`.

We incorporate some of the limitations encountered by observers so that more consistent comparisons can be made between observations and theory.

From this data cube, "observables" can be measured using observational pipelines. 
This package, once installed, is fully documented and tested.

<!-- badges: start -->
<a href="https://github.com/kateharborne/SimSpin/actions"><img src="https://github.com/kateharborne/SimSpin/actions/workflows/r.yml/badge.svg" alt="R-CMD-check"/></a>
<a href="https://app.codecov.io/gh/kateharborne/SimSpin"><img src="https://codecov.io/gh/kateharborne/SimSpin/branch/master/graph/badge.svg?token=2T1BDWZYSV" alt="codecov"/></a>
<a href="https://ascl.net/1903.006"><img src="https://img.shields.io/badge/ascl-1903.006-blue.svg?colorB=262255" alt="ascl:1903.006" /></a>
<!-- badges: end -->

Another implementation of this code (SimSpin v1.1.3) written in Julia is also available at [SimSpin.jl](https://github.com/kateharborne/SimSpin.jl) developed by [Gerry Gralton](https://github.com/gerrygralton). 

The purpose of [this guide](https://kateharborne.github.io/SimSpin/) is to act as a reference for users. 
It outlines the intended functionality and processes available in the `SimSpin` code. 
This is both a form of documentation and a collection of examples. 
For a basic walk-through of installing `SimSpin` and running a simple observation, please see the links on the left. 
