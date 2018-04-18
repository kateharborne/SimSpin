# Kate Harborne (last edit - 16/04/2018)
#'Mimicking seeing conditions for the IFU data cube.
#'
#'The purpose of this function is to mimic the spatial blurring caused by seeing and beam smearing.
#'
#'@param ifu_datacube The list output from the function \code{ifu_cube()} containing the mock IFU cube and the apperture region
#' image (\code{$appregion}).
#'@param sbinsize The size of each spatial bin in kpc, output from the function \code{obs_data_prep()}.
#'@param psf This parameter gives the user choice between a "Gaussian" or "Moffat" PSF.
#'@param fwhm This parameter will blur the observation using the specified PSF with full width half maximum as
#' specified here in arcseconds, in order to mimic the effects caused by beam smearing and seeing.
#'@param angular_size The conversion factor for a given redshift between kpc and arcseconds as produced by the \code{obs_data_prep()} function.
#'@return Returns a list that contains the original IFU cube produced by the function \code{ifu_cube()}, but blurred with the specified PSF,
#' \code{$blurcube}, and the bin labels for each dimension \code{$xbin_labels}, \code{$ybin_labels}, \code{$vbin_labels}.
#'@examples
#' \dontrun{
#'  data      = obs_data_prep()
#'  ifucube   = ifu_cube()
#'
#'  blur_cube(ifu_datacube   = ifucube,
#'            sbinsize       = data$sbinsize,
#'            psf            = "Moffat",
#'            fwhm           = 0.5,
#'            angular_size   = data$angular_size)
#' }

blur_cube = function(ifu_datacube, sbinsize, psf, fwhm, angular_size){

  fwhm_scaled     = (fwhm * angular_size)/ sbinsize                          # the fwhm of the beam smearing scaled to the image pixel dimensions
  sbin            = length(ifu_datacube$xbin_labels)                         # number of spatial bins in data cube
  vbin            = length(ifu_datacube$vbin_labels)                         # number of velocity bins in data cube

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
                                                                             # applying blurring PSF to each spatial plane via convolution
  blurcube    = blurcube * ifu_datacube$appregion                            # trimming blurred cube to data within the aperture

  output = list("cube"        = blurcube,
                "xbin_labels" = ifu_datacube$xbin_labels,
                "ybin_labels" = ifu_datacube$ybin_labels,
                "vbin_labels" = ifu_datacube$vbin_labels)

  return(output)
}




