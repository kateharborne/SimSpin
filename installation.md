---
title: SimSpin
layout: default
filename: installation.md
navigation:
  - title: Home
    url: /
  - title: Installation
    url: installation
  - title: Basic Usage
    url: basic_usage
  - title: Documentation
    url: documentation
  - title: Examples
    url: examples
---

Here, we explain how to install `SimSpin` on your machine, listing specific instructions for both MacOSX and Linux. 
C library dependencies, including FFTW and HDF5, need to be installed prior to `SimSpin` as outlined for individual operating systems below.  
If you already have the C library dependencies installed, skip ahead to [Installing SimSpin](#installing-simspin) to get started.

`SimSpin` is an open source R-package, registered with the Astrophysics Source Code Library, DOI: [1903.006](https://ascl.net/1903.006).
This package build has been tested on OSX and Linux operating systems using [Travis CI](https://travis-ci.org/github/kateharborne/SimSpin) and should install simply from within R using the instructions in [Installing SimSpin](#installing-simspin).

You will also (obviously) need an installation of [R](https://www.r-project.org/), which can be downloaded from your local mirror for your operating system. 
We suggest downloading the [RStudio IDE](https://rstudio.com/) for a friendly user environment. 
You can download these programs using the instructions at the links embedded above. 

Beyond base R itself, there are some dependencies that will be required that may not exist on your machine.
`SimSpin` requires FFTW and HDF5 libraries. 
Instructions for downloading these libraries are provided below for Mac OSX and Linux machines. 
If you encounter issues installing directly from within R, please check the guide specific for your system below. 

### Installing dependencies - Mac OSX
In order to install `SimSpin` you will need a copy of  [XCode](https://apps.apple.com/us/app/xcode/id497799835?mt=12) (11 or greater), which can be downloaded and installed from the Apple App Store for free. 
Once complete, you will also need to ensure that the Command Line Tools have also been installed by running the following line from within a Terminal window. 
*This is only necessary for Xcode versions prior to v13.1.*

```
xcode-select --install 
```

The FFTW and HDF5 packages can then be installed via Homebrew using the commands:

```
brew install fftw 
brew install hdf5
```

### Installing dependencies - Linux
To install the FFTW and HDF5 packages on a Linux machine:

```
sudo apt-get update
sudo apt-get libhdf5-dev
sudo apt-get libfftw3-dev
```

### Installing SimSpin
Assuming that you have successfully installed a copy of R on your machine, the most recent release of `SimSpin` can be installed from GitHub from within your R session using the following commands:

```R
install.packages("devtools")
library("devtools")
install_github("kateharborne/SimSpin")
```

Required R-package dependencies will be installed at this time. 
If you encounter any errors, please check that all dependency libraries are successfully installed according to your operating system using the instructions above. 
If all else fails, [report it as an issue on GitHub](https://github.com/kateharborne/SimSpin/issues). 

Once this installation is complete, load the package into your R session by typing:

```R
library("SimSpin")
```

Following this, the functionality and documentation of the package should all be available within your R session. 
To check the installation has been successful, try typing: 

```R
?SimSpin
```

to see the cover page of the `SimSpin` package documentation. 