---
layout: default
title: Installation
permalink: /installation/
description: "How to install SimSpin on your machine."
nav_order: 2
last_modified_date: "Wed, 16 February 2022 13:57:00 AWST"
---

# Installing SimSpin
{: .no_toc }

Here, we explain how to install SimSpin on your machine, listing specific instructions for both MacOSX and Linux. 
{: .fs-5 .fw-300 }
---
## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc .pb-3}

<!-- badges: start -->
<a href="https://github.com/kateharborne/SimSpin/actions"><img src="https://github.com/kateharborne/SimSpin/actions/workflows/r.yml/badge.svg" alt="R-CMD-check"/></a>
<a href="https://app.codecov.io/gh/kateharborne/SimSpin"><img src="https://codecov.io/gh/kateharborne/SimSpin/branch/master/graph/badge.svg?token=2T1BDWZYSV" alt="codecov"/></a>
<a href="https://ascl.net/1903.006"><img src="https://img.shields.io/badge/ascl-1903.006-blue.svg?colorB=262255" alt="ascl:1903.006" /></a>
<!-- badges: end -->

SimSpin is an open source R-package, registered with the Astrophysics Source Code Library, DOI: [1903.006](https://ascl.net/1903.006).
This package has been built and tested on OSX and Linux operating systems using [R-CMD-check](https://github.com/kateharborne/SimSpin/actions/workflows/r.yml) and coverage of tests has been measured using [CodeCov](https://codecov.io/gh/kateharborne/SimSpin/branch/master/graph/badge.svg?token=2T1BDWZYSV).
If you are familiar with R and astronomy software, you should be able to install simply from within R using the instructions in [Installing SimSpin](#installing-simspin-1) by clicking the button below. 
If you encounter problems, check for missing dependencies on your operating system using the instructions below. 
{: .fw-300 }

[Get started now](#installing-simspin-1){: .btn .btn-purple }

---

## Installing dependencies 

You will need an installation of **R** to run this software. 
This can be downloaded from your local mirror for your operating system.
We also suggest downloading the **RStudio IDE** for a friendly user environment. 
You can download these programs using the instructions at the buttons below. 

[Install R](https://www.r-project.org/){: .btn  }
[Install RStudio](https://rstudio.com/){: .btn  }

Beyond base R itself, there are some dependencies that will be required that may not exist on your machine.
`SimSpin` requires the C libraries **FFTW** and **HDF5**. 
Instructions for downloading these libraries are provided below for Mac OSX and Linux machines. 

---

### Mac OSX
In order to install `SimSpin` you will need a copy of  **XCode** (11 or greater), which can be downloaded and installed from the Apple App Store for free. 

[Install XCode](https://apps.apple.com/us/app/xcode/id497799835?mt=12){: .btn }

Once complete, you will also need to ensure that the Command Line Tools have also been installed by running the following line from within a Terminal window. 
*This is only necessary for Xcode versions prior to v13.1.*

```
xcode-select --install 
```

The FFTW and HDF5 packages can then be installed via the command line using Homebrew. Open a Terminal window and use the following commands:

```
brew install fftw 
brew install hdf5
```

---

### Linux
To install the FFTW and HDF5 packages on a Linux machine, open a Terminal window and use the commands:

```
sudo apt-get update
sudo apt-get libhdf5-dev
sudo apt-get libfftw3-dev
```

---

## Installing SimSpin

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