# SimSpin
SimSpin - A package for the kinematic analysis of galaxy simulations

The purpose of the Simspin R-package is to take a galaxy simulation and measure an the kinematics of that model as if it had been observed using an IFU. A kinematic data cube can be produced using the functions in this package; from this cube, "observables" can be measured. Specifically, the observable spin parameter, &#955;<sub>r</sub>. This package, once installed, is fully documented and tested.

To install directly into R:
```
> install.packages("devtools")
> library(devtools)
> install_github("kateharborne/SimSpin")
```
Simulation data needs to be laid out in HDF5 format for processing with SimSpin. If using a Gadget binary or HDF5 simulation output, see [create_SimSpinFile](https://github.com/kateharborne/create_SimSpinFile) for automatic generation of SimSpin compatible files.  If you would rather generate the SimSpin file independently, the expected format is outlined below. This SimSpin file can then be read into SimSpin using the function:

```
galaxy_data = SimSpin::sim_data(filename = SimSpin_example.hdf5)
```

This function produces a list that can be accessed by each of the basic SimSpin analysis functions listed below. While it is possible to use each function in the R-package in order and examine the output at each stage, there are three basic analysis functions designed to give the data in a user friendly format. We suggest first using these functions:

1. sim_analysis() - This function is designed to output the kinematic properties of the galaxy model to be observed. This provides the comparison to the kinematic observables produced in the following functions. 

2. build_datacube() - This function produces the kinematic data cube prior to kinematic analysis. This allows the user to take the cubes to use in some other form of analysis without having to calculate &#955;<sub>r</sub>.

3. find_kinematics() - This function produces a kinematic data cube and calculates the observed spin parameter, ellipticity, inclination and the corresponding flux, line-of-sight velocity and line-of-sight velocity dispersion images. 

By varying the effects of observational seeing, the measurement radius, projected inclination and distance, and the telescope parameters within the find_lambda() function, we can begin to understand how inherent limitations of observing galaxies can effect the measurement of &#955;<sub>r</sub> by comparing to the true spin parameter than is measured in the sim_analysis() function.

For more detailed examples of using each of the analysis functions above, please see the SimSpin vignettes published on [Rpubs](http://rpubs.com/kateharborne). This [web app](http://simspin.icrar.org) can also be used to explore the N-body analysis capabilities of SimSpin. 

## SimSpin format

Here we outline the expected file format accepted by SimSpin.  If you would like to generate this file automatically, a short Python function has been written that uses the [pynbody](https://github.com/pynbody/pynbody) package to read in various simulation data types and generate a SimSpin compatible HDF5 file. See [create_SimSpinFile](https://github.com/kateharborne/create_SimSpinFile).

If you would rather generate the SimSpin file independently, the expected file format is outlined below.

```
> SimSpin_example.hdf5

>> /PartType0           # Each particle type included in the simulation has its own group.
>>> /PartType0/Mass     # Each group then has a series of data sets assocaited,
>>> /PartType0/vx       #   including the position, velocity and Mass of each particle. 
>>> /PartType0/vy
>>> /PartType0/vz
>>> /PartType0/x
>>> /PartType0/y
>>> /PartType0/z

>> /PartType1 
>>> ...
```
We use the same PartType definition as Gadget: PartTypeX where 0 - gas, 1 - dark matter, 2 - disc, 3 - bulge, 4 - stars. For PartType0-3, each PartType group contains the same data sets as above. If the simulation contains stars, the Age and Metallicity information for each particle is also included:

```
> SimSpin_example.hdf5
>> /PartType4
>>> /PartType4/Age
>>> /PartType4/Mass
>>> /PartType4/Metallicity
>>> /PartType4/vx        
>>> /PartType4/vy
>>> /PartType4/vz
>>> /PartType4/x
>>> /PartType4/y
>>> /PartType4/z
```
If the file is set up in this way, the simulation data can easily be read into the SimSpin package. 

### References
A. Pontzen, R Roskar, G. Stinson and R. Woods, (2013), "pynbody: N-Body/SPH analysis for python",  Astrophysics Source Code Library, record ascl:1305.002
