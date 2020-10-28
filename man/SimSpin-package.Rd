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
telescope = telescope(type="SAMI")
observing_strategy = observing_strategy(z = 0.05, inc_deg = 90)
cube = build_datacube(simspin_file = system.file("extdata",
                                                 "SimSpin_example_Gadget_spectra.fst",
                                                  package = "SimSpin"),
                      telescope = telescope,
                      observing_strategy = observing_strategy)

}
\keyword{data}
