#' u-band SDSS filter responses at a range of wavelengths
#'
#' A dataset containing the relative response of the SDSS u-band filter across
#' the relevent wavelength range in Angstroms.
#'
#' @format A data table with 2 variables and 47 observations:
#' \describe{
#' \item{wave}{A numeric representing the wavelength at which the relative response is recorded, given in Angstroms}
#' \item{response}{A numeric representing the relative response of the filter at that wavelength}
#' }
#' @source \url{https://classic.sdss.org/dr7/instruments/imager/filters/u.dat}
#'
"filt_u_SDSS"

#' g-band SDSS filter responses at a range of wavelengths
#'
#' A dataset containing the relative response of the SDSS g-band filter across
#' the relevent wavelength range in Angstroms.
#'
#' @format A data table with 2 variables and 89 observations:
#' \describe{
#' \item{wave}{A numeric representing the wavelength at which the relative response is recorded, given in Angstroms}
#' \item{response}{A numeric representing the relative response of the filter at that wavelength}
#' }
#' @source \url{https://classic.sdss.org/dr7/instruments/imager/filters/g.dat}
#'
"filt_g_SDSS"

#' r-band SDSS filter responses at a range of wavelengths
#'
#' A dataset containing the relative response of the SDSS r-band filter across
#' the relevent wavelength range in Angstroms.
#'
#' @format A data table with 2 variables and 75 observations:
#' \describe{
#' \item{wave}{A numeric representing the wavelength at which the relative response is recorded, given in Angstroms}
#' \item{response}{A numeric representing the relative response of the filter at that wavelength}
#' }
#' @source \url{https://classic.sdss.org/dr7/instruments/imager/filters/r.dat}
#'
"filt_r_SDSS"

#' i-band SDSS filter responses at a range of wavelengths
#'
#' A dataset containing the relative response of the SDSS i-band filter across
#' the relevent wavelength range in Angstroms.
#'
#' @format A data table with 2 variables and 89 observations:
#' \describe{
#' \item{wave}{A numeric representing the wavelength at which the relative response is recorded, given in Angstroms}
#' \item{response}{A numeric representing the relative response of the filter at that wavelength}
#' }
#' @source \url{https://classic.sdss.org/dr7/instruments/imager/filters/i.dat}
#'
"filt_i_SDSS"

#' z-band SDSS filter responses at a range of wavelengths
#'
#' A dataset containing the relative response of the SDSS z-band filter across
#' the relevent wavelength range in Angstroms.
#'
#' @format A data table with 2 variables and 141 observations:
#' \describe{
#' \item{wave}{A numeric representing the wavelength at which the relative response is recorded, given in Angstroms}
#' \item{response}{A numeric representing the relative response of the filter at that wavelength}
#' }
#' @source \url{https://classic.sdss.org/dr7/instruments/imager/filters/z.dat}
#'
"filt_z_SDSS"

#' BC03lr - Bruzual & Charlot 2003 Low-Resolution Spectral Templates for Stellar Populations
#'
#' A dataset containing the grid of spectral templates for a stellar population
#' with a given stellar age and metallicity. At each age/metallicity bin, a
#' template spectrum is given in units of solar luminoisities per angstrom for
#' 1 solar mass of star formation. Values are given at a series of wavelengths,
#' contained within the dataset. Note that while the first age bin is a
#' template for Age = 0 stars, we have re-labelled this bin with a non-zero small
#' number to avoid log errors.
#'
#' @format A list with 4 elements:
#' \describe{
#' \item{Z}{A numeric array that describes the 6 possible metallicity bins provided by the template spectra.}
#' \item{Age}{A numeric array that describes the 221 possible age bins provided by the template spectra.}
#' \item{Wave}{A numeric array that describes the 1221 wavelengths (in Angstrom) at which the spectrum is provided.}
#' \item{Zspec}{A list containing the template spectrum for a galaxy. The first level of the list corresponds to
#' the metallicity of the population. Each list element then contains a numeric matrix, with each row describing the
#' spectrum that would be associated with a population of a given age.}
#' }
#' @source \url{http://www.bruzual.org/bc03}
#'
"BC03lr"

#' BC03hr - Bruzual & Charlot 2003 High-Resolution Spectral Templates for Stellar Populations
#'
#' A dataset containing the grid of spectral templates for a stellar population
#' with a given stellar age and metallicity. At each age/metallicity bin, a
#' template spectrum is given in units of solar luminoisities per angstrom for
#' 1 solar mass of star formation. Values are given at a series of wavelengths,
#' contained within the dataset. Note that while the first age bin is a
#' template for Age = 0 stars, we have re-labelled this bin with a non-zero small
#' number to avoid log errors.
#'
#' @format A list with 4 elements:
#' \describe{
#' \item{Z}{A numeric array that describes the 6 possible metallicity bins provided by the template spectra.}
#' \item{Age}{A numeric array that describes the 221 possible age bins provided by the template spectra.}
#' \item{Wave}{A numeric array that describes the 6900 wavelengths (in Angstrom) at which the spectrum is provided.}
#' \item{Zspec}{A list containing the template spectrum for a galaxy. The first level of the list corresponds to
#' the metallicity of the population. Each list element then contains a numeric matrix, with each row describing the
#' spectrum that would be associated with a population of a given age.}
#' }
#' @source \url{http://www.bruzual.org/bc03}
#'
"BC03hr"

#' EMILES - UV-extended E-MILES stellar population models
#'
#' A dataset containing the grid of spectral templates for a stellar population
#' with a given stellar age and metallicity. At each age/metallicity bin, a
#' template spectrum is given in units of solar luminoisities per angstrom for
#' 1 solar mass of star formation. Values are given at a series of wavelengths,
#' contained within the dataset.
#'
#' @format A list with 4 elements:
#' \describe{
#' \item{Z}{A numeric array that describes the 12 possible metallicity bins provided by the template spectra.}
#' \item{Age}{A numeric array that describes the 53 possible age bins provided by the template spectra.}
#' \item{Wave}{A numeric array that describes the 53689 wavelengths (in Angstrom) at which the spectrum is provided.}
#' \item{Zspec}{A list containing the template spectrum for a galaxy. The first level of the list corresponds to
#' the metallicity of the population. Each list element then contains a numeric matrix, with each row describing the
#' spectrum that would be associated with a population of a given age.}
#' }
#' @source \url{http://miles.iac.es/}
#'
"EMILES"
