# SimSpin
SimSpin - A package for the kinematic analysis of galaxy simulations

The purpose of the Simspin R-package is to take a galaxy simulation and measure an the kinematics of that model as if it had been observed using an IFU. A kinematic data cube can be produced using the functions in this package; from this cube, "observables" can be measured. Specifically, the observable spin parameter, &#955;_r. This package, once installed, is fully documented and tested.

While it is possible to use each function in the R-package in order and examine the output at each stage, there are three basic analysis functions designed to give the data in a user friendly format. We suggest first using these functions:

1. sim_analysis() - This function is designed to output the kinematic properties of the galaxy model to be observed. This provides the comparison to the kinematic observables produced in the following functions. 

2. build_datacube() - This function produces the kinematic data cube prior to kinematic analysis. This allows the user to take the cubes to use in some other form of analysis without having to calculate &#955;_r.

3. find_lambda() - This function produces a kinematic data cube and calculates the observed spin parameter, ellipticity, inclination and the corresponding flux, line-of-sight velocity and line-of-sight velocity dispersion images. 

By varying the effects of observational seeing, the measurement radius, projected inclination and distance, and the telescope parameters within the find_lambda() function, we can begin to understand how inherent limitations of observing galaxies can effect the measurement of &#955;_r by comparing to the true spin parameter than is measured in the sim_analysis() function.
