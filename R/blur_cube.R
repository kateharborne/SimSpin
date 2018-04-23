# Kate Harborne (last edit - 23/04/2018)
#'Mimicking seeing conditions for the IFU data cube.
#'
#'The purpose of this function is to mimic the spatial blurring caused by seeing and beam smearing.
#'
#'@param ifu_datacube The list output from the function \code{\link{ifu_cube}} containing the mock
#' IFU cube and the aperture region image (\code{$appregion}).
#'@param sbinsize The size of each spatial bin in kpc, output from the function
#' \code{\link{obs_data_prep}}.
#'@param psf The type of PSF with choices of "Gaussian" or "Moffat".
#'@param fwhm The full-width half-maximum size of the chosen PSF specified here in arcseconds.
#'@param angular_size The conversion factor for a given redshift between kpc and arcseconds as
#' produced by the \code{\link{obs_data_prep}} function.
#'@return Returns a list that contains:
#'\item{\code{$blurcube}}{The original IFU cube produced by the function \code{\link{ifu_cube}},
#' but blurred with the specified PSF.}
#'\item{\code{$xbin_labels}}{Bin labels for the x-spatial dimension.}
#'\item{\code{$ybin_labels}}{Bin labels for the y-spatial dimension.}
#'\item{\code{$vbin_labels}}{Bin labels for the velocity dimension.}
#'@examples
#'  data      = obs_data_prep(filename = system.file("extdata", 'S0_vignette', package="SimSpin"))
#'  ifucube   = ifu_cube(obs_data = data)
#'
#'  blurcube  = blur_cube(ifu_datacube   = ifucube,
#'                        sbinsize       = data$sbinsize,
#'                        psf            = "Moffat",
#'                        fwhm           = 0.5,
#'                        angular_size   = data$angular_size)


blur_cube = function(ifu_datacube, sbinsize, psf, fwhm, angular_size){

  fwhm_scaled     = (fwhm * angular_size)/ sbinsize  # the fwhm scaled to image pixel dimensions
  sbin            = length(ifu_datacube$xbin_labels) # number of spatial bins in data cube
  vbin            = length(ifu_datacube$vbin_labels) # number of velocity bins in data cube

  if (psf == "Gaussian"){
    psf_k = ProFit::profitMakeGaussianPSF(fwhm = fwhm_scaled, dim = c(25,25))
  }
  if (psf == "Moffat"){
    psf_k = ProFit::profitCubaMoffat(fwhm = fwhm_scaled, mag = 1, con = 5, dim = c(25,25))
  }
                                                     # creating each PSF kernel
  blurcube        = array(data = 0, dim = c(sbin, sbin, vbin))
  for (c in 1:vbin){
    blurcube[,,c] = ProFit::profitBruteConv(ifu_datacube$cube[,,c], psf_k)
  }
                                                     # apply PSF to spatial plane via convolution
  blurcube    = blurcube * ifu_datacube$appregion    # trimming cube to data within the aperture

  output = list("cube"        = blurcube,
                "xbin_labels" = ifu_datacube$xbin_labels,
                "ybin_labels" = ifu_datacube$ybin_labels,
                "vbin_labels" = ifu_datacube$vbin_labels)

  return(output)
}




