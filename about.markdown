---
layout: default
title: About
nav_order: 2
description: "Some background about SimSpin."
permalink: /about/
last_modified_date: "Tue, 20 June 2023 14:50:00 AWST"
---

# SimSpin v2.5.x
{: .fs-9 }

A package for producing mock observations of simulations
{: .fs-6 .fw-300 .mb-3 .lh-tight }

<img align="right" src="/SimSpin/assets/images/logo.png" width="175" height="175" />
{: .pl-4 .pb-1 } 

SimSpin v1.1.0 was written as part of the PhD Thesis "From Particles to Pixels: Using numerical simulations to investigate observable galaxy kinematics" by K.E. Harborne (2020). 
{: .fs-5 .fw-300 }

This code has since evolved to both construct *kinematic* and *spectral* data cubes, as denoted by the upgrade to version 2.X.X. Further details about these version numbers is described below. 

SimSpin has been used for a number of research publications including:

-  ["The diverse nature and formation paths of slow rotator galaxies in the EAGLE simulations"](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.4372L/abstract)(2022) *C.D.P. Lagos, E. Emsellem, J. van de Sande, **K.E. Harborne**, L. Cortese, T. Davison, C. Foster and R. Wright*, 
-  ["The MAGPI survey: Science goals, design, observing strategy, early results and theoretical framework"](https://ui.adsabs.harvard.edu/abs/2021PASA...38...31F/abstract)(2021) *C. Foster, J.T. Mendel, C.D.P. Lagos, E. Wisnioski, T. Yuan, F. D'Eugenio, T. M. Barone, **K.E. Harborne** and others* 
-  ["Recovering λR and V/σ from seeing-dominated IFS data"](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.2018H/abstract)(2020) **K.E. Harborne**, *J. van de Sande, L. Cortese, C. Power, A.S.G. Robotham, C.D.P. Lagos and S. Croom*  
-  ["A numerical twist on the spin parameter, λR"](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483..249H/abstract)(2019) **K.E. Harborne**, *C. Power, A.S.G. Robotham, L. Cortese and D.S. Taranu* 

[Click here to see on NASA ADS](https://ui.adsabs.harvard.edu/abs/2020PASA...37...16H/citations){: .btn }

---

#### Referencing the code
{: .fs-4 .pb-4 } 

If you use this code for your research, please make sure to include a citation to the source code paper. This reference will soon be updated for the newer v2.X.X. You can do so using the CITATION.cff file within the GitHub repo, or by including the following within your bibliography:
{: .fs-5 .fw-300 }

**NEW CITATION** With the release of v2.5.0, we have a new publication to record the significant upgrades since v1.1.1. This is under review, but available on the arXiv [here](https://arxiv.org/abs/2307.02618).
{: .note }

[Download Citation](https://github.com/kateharborne/SimSpin/blob/main/CITATION.cff){: .btn .btn-purple }
[See the paper](https://arxiv.org/abs/2307.02618){: .btn .btn-purple }

- K.E. Harborne,  A. Serene, E.J.A. Davies,  C. Derkenne, S. Vaughan, A. Burdon, C. del P. Lagos,  R. McDermid,  S. O'Toole, C. Power, A.S.G. Robotham, G. Santucci, and R. Tobar (2023), [SimSpin v2.5.0 -- Constructing synthetic spectral IFU cubes for comparison with observational surveys](https://arxiv.org/abs/2307.02618)

-   K.E. Harborne, C. Power and A.S.G. Robotham, (2020), ["SIMSPIN - Constructing mock IFS kinematic data cubes"](https://ui.adsabs.harvard.edu/abs/2020PASA...37...16H/abstract), Publications of the Astronomical Society of Australia, Volume 37, article id. e016

-   K.E. Harborne, (2019), ["SimSpin: Kinematic analysis of galaxy simulations"](https://ui.adsabs.harvard.edu/abs/2019ascl.soft03006H/abstract), Astrophysics Source Code Library, record ascl:1903.006

---
#### SimSpin design and versioning
{: .fs-4 .pb-4 } 

SimSpin is a well-documented and tested code that follows Semantic versioning standards. Further details of the design specifications are described for interested readers below. 

Developments in SimSpin always occur on Git branches with the prefix "dev-". 
{: .note }

In these branches, the code will not be stable as they keep a record of ongoing development changes to functionality in the code. Do not pull SimSpin from a "dev-" branch without being aware that there may be ongoing adjustments and breaking changes within these iterations. Once the development of a new feature has been completed, a series of tests and checks will be performed to ensure that this development branch will not introduce any breaking changes to the stable Main branch of the code. If these tests all pass without fault, this development branch will be merged with the Main branch of the code. 

Every time a development branch is merged with the Main branch, SimSpin will receive an updated version number at its smallest digit (i.e. v2.5.X). At this point, a summary of changes will be added to the NEWS table within the top directory of the code. This table also notes the new tagged version and the commit reference at which the branch was merged.

In the scenario that a large number of small changes have been made, we may choose to release a new complete version of the code. Code releases occur when the second digit of the semantic versioning has been updated (i.e. v2.X.0). In these scenarios, the code will always be backwards compatible i.e. existing analysis written will not error, though you may receive some warnings regarding the use of older syntax and the withdrawal of support for a given feature in the future. 

When a large or breaking change is incorporated into the code, the largest digit in the series (i.e. vX.0.0) will be updated. An update of this sort may introduce changes that cause current analysis codes to error. These changes will be made clear within the documentation and NEWS, as well as being communicated via useful error messages to the user. 

---

#### SimSpin on the web
{: .fs-4 }

<img align="right" src="/SimSpin/assets/images/ADACSlogo_LR.png" width="240" height="100" />
{: .pl-6 .pb-3 } 

In 2022, SimSpin was awarded an Astronomy Data and Computing Services (ADACS) allocation for the development of a web application with the full capability of the code accessible through a graphical user interface. 
{: .fs-5 .fw-300 }

[Launch the app!](https://simspin.datacentral.org.au/app/){: .btn .btn-purple }
{: .lh-tight }

<img align="centre" src="/SimSpin/assets/images/simspin_webapp.png" width="600" height="438" />
{: .pt-4 .pb-1 } 

This interactive tool can be used for visualising galaxy simulations in a browser, enabling the dynamic exploration and visualisation of simulated galaxies and downloading of generated files for offline use.

The application is a React Single Page App communicating asynchronously with a RESTful API to access the SimSpin package. The app allows for instant data exploration via a dedicated viewer, where authenticated users can re-visit previous queries and share results with others. All services (web, db, redis-cache, celery-workers, celery-node, celery-beat) are containerised and managed by docker compose, such that the project is easily re-deployable. The API is fully documented, and comes with an API Schema (adhering to the OpenAPI Specification) to aid users in calling the API from other services.

The SimSpin app removes the barrier of entry for novice astronomers (no R installation required, minimal tool understanding, instant data visualisation), providing an accessible and time-saving tool for simulated galaxy visualisations. 


*KH acknowledges the funding support of the Australian Research Council Centre of Excellence for All-Sky Astrophysics in 3 Dimensions (ASTRO3D), through project number CE170100013. This work has been made possible through the Astronomy Data and Computing Services (ADACS) with direct support from Elizabeth Mannering, Felipe Jimenez-Ibarra and Simon O’Toole. This work has also been supported by resources provided by the Pawsey Supercomputing Centre with funding from the Australian Government and the Government of Western Australia*
{: .fs-3 .fw-300 }