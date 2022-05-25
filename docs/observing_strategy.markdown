---
layout: default
title: observing_strategy
parent: Documentation
nav_order: 3
last_modified_date: "Wed, 16 February 2022 15:57:00 AWST"
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
observing_strategy(dist_z       = 0.05,  #
                   inc_deg      = 70,
                   twist_deg    = 0,
                   pointing_kpc = c(0,0),
                   blur         = F)
```