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
Version: \tab 0.1.0\cr
Date: \tab 2018-04-18\cr
License: \tab LGPL-3\cr
Depends: \tab R (>= 3.0), Rcpp, snapshot, magicaxis, celestial, ProFit, fields, RColorBrewer\cr
Suggests: \tab knitr, rmarkdown \cr
}
}
\author{
\packageAuthor{SimSpin}
Maintainer: \packageMaintainer{SimSpin}
}

\examples{
S0_lambda = find_lambda(filename     = system.file("extdata", 'S0_vignette', package="SimSpin"),
                        ptype        = NA,
                        z            = 0.05,
                        fov          = 15,
                        ap_shape     = "circular",
                        central_wvl  = 4800,
                        lsf_fwhm     = 2.65,
                        pixel_sscale = 0.5,
                        pixel_vscale = 1.04,
                        inc_deg      = 70,
                        m2l_disc     = 1,
                        m2l_bulge    = 1,
                        threshold    = 25,
                        measure_type = list(type="fit", fac=1))
}
\keyword{data}
