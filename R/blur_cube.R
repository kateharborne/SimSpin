# Kate Harborne (last edit - 13/09/2017)
#'Mimicking seeing conditions for the IFU data cube.
#'
#'The purpose of this function is to mimic the spatial blurring caused by seeing and beam smearing.
#'
#'@param ifu_datacube The list output from the function \code{ifu_cube()} containing the mock IFU cube and the apperture region image (\code{$appregion}).
#'@param sbinsize The size of each spatial bin in kpc, output from the function \code{obs_data_prep()}.
#'@param psf This parameter gives the user choice between a "Gaussian" or "Moffat" PSF. The default is Gaussian.
#'@param fwhm This parameter will blur the observation using the specified PSF with full width half maximum as
#' specified here in arcseconds, in order to mimic the effects caused by beam smearing and seeing.
#'@param angular_size This parameter only needs to be specified if a fwhm is given. The conversion factor for a given redshift between
#' kpc and arcseconds as produced by the \code{obs_data_prep()} function.
#'@return Returns a list that contains the observed \eqn{\lambda}_R value (\code{$obs_lambdar}), and three matricies reflecting the images produced
#' in IFU surveys - a luminosity counts image (\code{$counts_img}), a velocity image (\code{$velocity_img}) and a velocity dispersion image
#' (\code{$dispersion_img}) - and the coordinates of the effective radius ellipse within which \eqn{\lambda}_R is measured.
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

  fwhm_scaled     = (fwhm * angular_size)/ sbinsize  # the fwhm of the beam smearing scaled to the image pixel dimensions
  sbin            = length(ifu_datacube$xbin_labels) # number of spatial bins in data cube
  vbin            = length(ifu_datacube$vbin_labels) # number of velocity bins in data cube
  xbin_labels     = ifu_datacube$blur_info[1]
  zbin_labels     = ifu_datacube$blur_info[2]

  if (psf == "Gaussian"){
    psf_k         = ProFit::profitMakeGaussianPSF(fwhm = fwhm_scaled, dim = c(25,25))
  }
  if (psf == "Moffat"){
    psf_k         = ProFit::profitCubaMoffat(fwhm = fwhm_scaled, mag = 1, con = 5, dim = c(25,25))
  }

  blurcube        = array(data = 0, dim = c(sbin, sbin, vbin))
  for (c in 1:vbin){
    blurcube[,,c] = ProFit::profitBruteConv(ifu_datacube$cube[,,c], psf_k)
  }

  blurcube    = blurcube * ifu_datacube$appregion
  counts_img  = apply(blurcube, c(1,2), sum)
  counts_flat = c(counts_img)

  xcen       = .meanwt(xbin_labels, counts_flat) # calculating the centre of the galaxy
  zcen       = .meanwt(zbin_labels, counts_flat)
  sx         = sqrt(.varwt(xbin_labels, counts_flat, xcen)) # calculating the standard deviation along each direction
  sz         = sqrt(.varwt(zbin_labels, counts_flat, zcen))
  sxz        = .covarwt(xbin_labels, zbin_labels, counts_flat, xcen, zcen) # calculating the covariance
  temprad    = .cov2eigval(sx, sz, sxz) # solving for the eigenvalues to give the major and minor axes lengths
  major      = sqrt(abs(temprad$hi))
  minor      = sqrt(abs(temprad$lo))
  axis_ratio = data.frame("a" = major, "b" = minor)

  output = list("cube" = blurcube, "xbin_labels" = ifu_datacube$xbin_labels, "ybin_labels" = ifu_datacube$ybin_labels, "vbin_labels" = ifu_datacube$vbin_labels, "axis_ratio" = axis_ratio)

  return(output)
}




