\name{SimSpin-package}
\alias{SimSpin-package}
\alias{SimSpin}
\docType{package}

\title{
\packageTitle{SimSpin}
}
\description{
\packageDescription{SimSpin}
}
\details{
\tabular{ll}{
Package: \tab SimSpin\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2020-10-28\cr
License: \tab LGPL-3\cr
Depends: \tab R (>= 3.0), Rcpp\cr
Suggests: \tab testthat\cr
}
}
\author{
\packageAuthor{SimSpin}
Maintainer: \packageMaintainer{SimSpin}
}

\examples{
\dontrun{
ss_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5", package = "SimSpin")
temp_loc = tempdir()
make_simspin_file(ss_eagle, output = paste(temp_loc, "spectra.fst", sep=""))
cube = build_datacube(simspin_file = paste(temp_loc, "spectra.fst", sep=""),
                     telescope = telescope(type="SAMI"),
                     observing_strategy = observing_strategy())
unlink(paste(temp_loc, "spectra.fst", sep=""))
}
}
\keyword{data}
