# SimSpin

<!-- badges: start -->

<a href="https://github.com/kateharborne/SimSpin/actions"><img src="https://github.com/kateharborne/SimSpin/actions/workflows/r.yml/badge.svg" alt="R-CMD-check"/></a> <a href="https://app.codecov.io/gh/kateharborne/SimSpin"><img src="https://codecov.io/gh/kateharborne/SimSpin/branch/master/graph/badge.svg?token=2T1BDWZYSV" alt="codecov"/></a> <a href="https://ascl.net/1903.006"><img src="https://img.shields.io/badge/ascl-1903.006-blue.svg?colorB=262255" alt="ascl:1903.006"/></a>

<!-- badges: end -->

<img src="https://raw.githubusercontent.com/kateharborne/SimSpin.jl/master/docs/src/assets/logo.png" align="right" width="200" height="200" style="padding-left:10px"/>

<p>&nbsp;</p>

v2.3.17 - A package for producing mock observations:

SimSpin allows you to take a simulation of a galaxy and produce a data cube in the style of an Integral Field Spectroscopy (IFS) instrument. A mock observation is produced using three simple steps:

1.  Read in your particle data and produce the relevant spectra using the `make_simspin_file` function.
2.  Setup the observation by defining your `telescope` and `observing_strategy`.
3.  Build your data cube using the `build_datacube`.

The data cube produced can be saved as a FITS file for ease of comparison with real observations. This package, once installed, is fully documented and tested.

Another implementation of this code (SimSpin v1.1.3) written in Julia is also available at [SimSpin.jl](https://github.com/kateharborne/SimSpin.jl) developed by [Gerry Gralton](https://github.com/gerrygralton).

## Package installation

To install directly into R:

    install.packages("devtools")
    library(devtools)
    install_github("kateharborne/SimSpin")

Else, please `fork` a copy of this repository for your own development using the `code` button to the upper right.

## Package examples

The following lines will take you from a particle simulation of a galaxy, an example of which is included in this package, through to the production of a mock-IFS data cube.

    ss_file = make_simspin_file(filename = system.file("extdata","SimSpin_example_Gadget",
                                                        package = "SimSpin"),
                                write_to_file = FALSE) # generate a spectra file
                                
    SAMI = telescope(type="SAMI")                      # initialise a telescope

    strategy = observing_strategy(dist_z = 0.05)       # initialise obsering conditions

    cube = build_datacube(simspin_file = ss_file,      # build data cube
                          telescope = SAMI,
                          observing_strategy = strategy,
                          method = "spectral")
                                

For a longer example of how each function can be used, please take a look at the documentation for the package. Short examples for each function are provided in each function's documentation, as well as an explanation of each of the possible input variables.

From within R, you can display the package documentation by typing `?SimSpin` and select "Index" at the bottom of the page to view all available functions. Alternatively, type `?` followed by the function name to see function specific documentation.

Longer examples are published [here](https://rpubs.com/kateharborne) and demonstrate a walk-through the basic code operation for SimSpin v2.1.1 (and v1.1.3).

If you have any further questions or requests for features in the code, report an issue or drop an email to [katherine.harborne\@uwa.edu.au](mailto:katherine.harborne@uwa.edu.au).

### Citation

If you use this code in any of your own published research, please make sure to include the following citation in your bibliography:

K.E. Harborne, C.Power and A.S.G. Robotham, (2020), ["SIMSPIN - Constructing mock IFS kinematic data cubes"](https://ui.adsabs.harvard.edu/abs/2020PASA...37...16H/abstract), Publications of the Astronomical Society of Australia, Volume 37, article id. e016

K.E. Harborne, (2019), ["SimSpin: Kinematic analysis of galaxy simulations"](https://ui.adsabs.harvard.edu/abs/2019ascl.soft03006H/abstract), Astrophysics Source Code Library, record ascl:1903.006
