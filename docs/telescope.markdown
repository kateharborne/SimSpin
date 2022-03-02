---
layout: default
title: telescope
parent: Documentation
nav_order: 2
last_modified_date: "Wed, 16 February 2022 15:57:00 AWST"
---

# Defining the properties of the telescope

In order to build a data cube, we need to describe the properties of the telescope used to observe the model. 
{: .fs-5 .fw-300 .pb-2 }

Using the `telecope` function, a single telescope object is generated with a set number of expected properties. These 
{: .fw-300 }


[See an example](#example){: .btn .btn-purple }

---

The following code shows the default parameters used in the `make_simspin_file` function. Calling the function without specifying anything other than the required input `filename` will produce a SimSpin file saved at the same directory location as the input simulation file with the following defaults. 

```R
make_simspin_file(filename,                         # REQUIRED input file 
                  disk_age = 5, disk_Z = 0.024,     # for N-body disk particles
                  bulge_age = 10, bulge_Z = 0.001,  # for N-body bulge particles
                  cores = 1, write_to_file = TRUE, 
                  output, overwrite = F,
                  template = "BC03lr",              # template choice for spectra
                  centre = NA, half_mass = NA,      # alignment choice
                  sph_spawn_n = 1)                  # gas smoothing choice

```

---

## Parameters
