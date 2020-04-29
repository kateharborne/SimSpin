# Kate Harborne (last edit - 15/11/2019)
#'Mimicking seeing conditions for the IFU data cube.
#'
#'The purpose of this function is to mimic the spatial blurring caused by seeing and beam smearing.
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param ifu_datacube The list output from the function \code{\link{ifu_cube}} containing the mock
#' IFU cube and the aperture region image (\code{$ap_region}).
#'@param psf The type of PSF with choices of "Gaussian" or "Moffat".
#'@param fwhm The full-width half-maximum size of the chosen PSF specified here in arcseconds.
#'@param threshold The magnitude limit of the observation in AB mag.
#'@return Returns a list that contains:
#' \item{\code{$blurcube}}{The original IFU cube produced by the function \code{\link{ifu_cube}},
#'  but blurred with the specified PSF.}
#' \item{\code{$ap_region}}{An image that describes the shape of the aperture such that any further
#'  convolutions that are applied to mimic beam smearing or seeing can be trimmed to the
#'  appropriate aperture shape.}
#' \item{\code{$xbin_labels}}{Bin labels for the x-spatial dimension.}
#' \item{\code{$ybin_labels}}{Bin labels for the y-spatial dimension.}
#' \item{\code{$vbin_labels}}{Bin labels for the velocity dimension.}
#'@examples
#'  galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#'  data        = obs_data_prep(galaxy_data)
#'  fluxes      = flux_grid(obs_data = data)
#'  cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#'
#'  blurcube  = blur_cube(ifu_datacube   = ifucube,
#'                        sbinsize       = data$sbinsize,
#'                        psf            = "Moffat",
#'                        fwhm           = 0.5,
#'                        angular_size   = data$angular_size)
#'


blur_cube = function(obs_data, ifu_datacube, psf, fwhm, threshold){

  fwhm_scaled     = (fwhm * obs_data$angular_size)/ obs_data$sbinsize  # the fwhm scaled to image pixel dimensions
  sbin            = obs_data$sbin # number of spatial bins in data cube
  vbin            = obs_data$vbin # number of velocity bins in data cube

  if (sbin < 25 && sbin %% 2 != 0){psf_dim = sbin} else if (sbin < 25 && sbin %% 2 == 0){psf_dim = sbin-1} else {psf_dim = 25}

  if (psf == "Gaussian"){
    psf_k = ProFit::profitMakeGaussianPSF(fwhm = fwhm_scaled, dim = c(psf_dim,psf_dim))
  }
  if (psf == "Moffat"){
    psf_k = ProFit::profitCubaMoffat(fwhm = fwhm_scaled, mag = 1, con = 5, dim = c(psf_dim,psf_dim))
  }
                                                     # creating each PSF kernel
  blurcube        = array(data = 0, dim = c(sbin, sbin, vbin))
  for (c in 1:vbin){
    blurcube[,,c] = ProFit::profitBruteConv(ifu_datacube$cube[,,c], psf_k)
  }

  threshold_flux = ProSpect::magAB2Jansky(threshold)

  for (i in 1:vbin){
    below_threshold = which(blurcube[,,i]<threshold_flux)
    blurcube[,,i][below_threshold] = 0
  }
                                                     # apply PSF to spatial plane via convolution
  blurcube    = blurcube * ifu_datacube$ap_region    # trimming cube to data within the aperture

  output = list("cube"        = blurcube,
                "ap_region"   = ifu_datacube$ap_region,
                "xbin_labels" = ifu_datacube$xbin_labels,
                "zbin_labels" = ifu_datacube$zbin_labels,
                "vbin_labels" = ifu_datacube$vbin_labels)

  return(output)
}




