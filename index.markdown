---
layout: default
title: Home
nav_order: 1
description: "SimSpin is a package for producing mock IFS observations of galaxy simulations."
permalink: /
last_modified_date: "Friday 28 October 2022 13:57:00 AWST"
---

# SimSpin v2.4.1
{: .fs-9 }

A package for producing mock observations of simulations
{: .fs-6 .fw-300 .mb-3 .lh-tight }

[View on GitHub](https://github.com/kateharborne/SimSpin){: .btn .btn-purple }
[Download Citation](https://github.com/kateharborne/SimSpin/blob/master/CITATION.cff){: .btn .btn-purple }
[View on ADS](https://ui.adsabs.harvard.edu/abs/2019ascl.soft03006H/abstract){: .btn .btn-purple }
{: .lh-tight }

<img align="right" src="assets/images/logo.png" width="175" height="175" />
{: .pl-4 .pb-1 } 

[SimSpin](https://github.com/kateharborne/SimSpin) is an R-package designed to take a simulation of a galaxy and produce a "mock" integral field spectroscopy (IFS) observation.
{: .fs-5 .fw-300 }

This software can be used to produce a synthetic **data cube** - i.e. spatial information in projection (*xy*) with spectral or kinematic information along the line-of-sight (*z*). 

A mock data cube can be produced in three simple steps:

  1. Read in your particle data and produce the relevant spectra using the [`make_simspin_file`](docs/make_simspin_file) function.
  1. Setup the observation by defining your [`telescope`](docs/telescope.markdown) and [`observing_strategy`](docs/observing_strategy.markdown).
  1. Build your data cube using the [`build_datacube`](docs/build_datacube.markdown).
{: .lh-tight }

<img align="centre" src="assets/images/simspin_v2_wo_logo.png" width="600" height="438" />
{: .pt-4 .pb-1 } 

We incorporate limitations encountered by observers so that more consistent comparisons can be made between observations and theory.

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
