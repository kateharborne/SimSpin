# Kate Harborne (last edit - 13/09/2017)
#'Calculating the observable spin parameter, \eqn{\lambda}_R.
#'
#'The purpose of this function is to calculate the spin parameter that would be observed given an IFU data cube.
#'This function will also mimic the spatial blurring caused by seeing and beam smearing if the \code{fwhm} parameter is supplied.
#'
#'@param ifu_images The list output from the function \code{ifu_img()} containing the mock IFU images (\code{$counts_img},
#'\code{$velocity_img}, \code{$dispersion_img}) and the apperture region image (\code{$appregion}).
#'@param reff_axisratio The semi-major and semi-minor axes output from the \code{find_reff()} function.
#'@param sbinsize The size of each spatial bin in kpc, output from the function \code{obs_data_prep()}.
#'@param psf \emph{Optional} This parameter gives the user choice between a "Gaussian" or "Moffat" PSF. The default is Gaussian.
#'@param fwhm \emph{Optional} This parameter will blur the observation using the specified PSF with full width half maximum as
#' specified here in arcseconds, in order to mimic the effects caused by beam smearing and seeing.
#'@param angular_size \emph{Optional} This parameter only needs to be specified if a fwhm is given. The conversion factor for a given redshift between
#' kpc and arcseconds as produced by the \code{obs_data_prep()} function.
#'@return Returns a list that contains the observed \eqn{\lambda}_R value (\code{$obs_lambdar}), and three matricies reflecting the images produced
#' in IFU surveys - a luminosity counts image (\code{$counts_img}), a velocity image (\code{$velocity_img}) and a velocity dispersion image
#' (\code{$dispersion_img}) and the coordinates of the effective radius ellipse within which \eqn{\lambda}_R is measured.
#'@examples
#' \dontrun{
#'  data      = obs_data_prep()
#'  ifu_imgs  = ifu_img()
#'  reff_data = find_reff()
#'
#'  obs_lambda(ifu_images     = ifu_imgs,
#'             reff_axisratio = reff_data,
#'             sbinsize       = data$sbinsize)
#'
#'  obs_lambda(ifu_images     = ifu_imgs,
#'             reff_axisratio = reff_data,
#'             sbinsize       = data$sbinsize,
#'             psf            = "Moffat",
#'             fwhm           = 0.5,
#'             angular_size   = data$angular_size)
#' }
#'

obs_lambda = function(ifu_images, reff_axisratio, sbinsize, psf = "Gaussian", fwhm, angular_size = NULL){
  if (missing(fwhm)) {
    sbin            = nrow(ifu_images$appregion) # number of spatial bins in data cube
    calcregion_reff = array(data = rep(0,(sbin*sbin)), dim=c(sbin,sbin)) # for calculating lambdaR within the specified reff
    radius          = matrix(data = 0, nrow=sbin, ncol=sbin) # creating the radial positions within the calcregion

    xcentre = sbin/2 + 0.5
    ycentre = sbin/2 + 0.5 # finding the centre pixel of the image
    a = reff_axisratio$a / sbinsize
    b = reff_axisratio$b / sbinsize

    for (x in 1:sbin){
      for (y in 1:sbin){
        xx = x - xcentre
        yy = y - ycentre
        rr = (xx^2 / a^2) + (yy^2 / b^2)
        if (rr <= 1){
          calcregion_reff[x,y] = 1
          radius[x,y]          = sqrt(xx^2 + yy^2) * sbinsize
        }
      }
    } # creating two arrays - one of a multiplier to consider lambdaR within Reff,
    # and one of the radial values within Reff

    counts     = ifu_images$counts_img * calcregion_reff
    velocity   = ifu_images$velocity_img * calcregion_reff
    dispersion = ifu_images$dispersion_img * calcregion_reff
    lambda     = sum(counts * radius * abs(velocity))/sum(counts * radius * (sqrt((velocity * velocity) + (dispersion * dispersion))))

    elli_x = seq(-a, a, length.out = 500)
    elli_y = (b / a) * sqrt(a^2 - elli_x^2)
    reff_elli = matrix(data=NA, nrow = 1000, ncol = 2)
    reff_elli[1:1000,1] = c(elli_x, rev(elli_x)) + rep(sbin/2, 1000)
    reff_elli[1:1000,2] = c(elli_y, rev(elli_y)*-1) + rep(sbin/2, 1000)

    if (max(reff_elli)>sbin){cat("WARNING: reff > aperture, the value of $obs_lambdar produced will not be the true value evaluated at reff.")}

    output = list("obs_lambdar"    = lambda,
                  "counts_img"     = ifu_images$counts_img,
                  "velocity_img"   = ifu_images$velocity_img,
                  "dispersion_img" = ifu_images$dispersion_img,
                  "reff_ellipse"   = reff_elli)

  } else {
    if (is.null(angular_size)){stop("Please enter an angular_size and re-run obs_lambda()")}
    sbin            = nrow(ifu_images$appregion) # number of spatial bins in data cube
    fwhm_scaled     = (fwhm * angular_size)/ sbinsize # the fwhm of the beam smearing scaled to the image pixel dimensions
    calcregion_reff = array(data = rep(0,(sbin*sbin)), dim=c(sbin,sbin)) # for calculating lambdaR within the specified reff
    radius          = matrix(data = 0, nrow=sbin, ncol=sbin) # creating the radial positions within the calcregion

    xcentre = sbin/2 + 0.5
    ycentre = sbin/2 + 0.5 # finding the centre pixel of the image
    a = reff_axisratio$a / sbinsize
    b = reff_axisratio$b / sbinsize

    for (x in 1:sbin){
      for (y in 1:sbin){
        xx = x - xcentre
        yy = y - ycentre
        rr = (xx^2 / a^2) + (yy^2 / b^2)
        if (rr <= 1){
          calcregion_reff[x,y] = 1
          radius[x,y]          = sqrt(xx^2 + yy^2) * sbinsize
        }
      }
    } # creating two arrays - one of a multiplier to consider lambdaR within Reff,
    # and one of the radial values within Reff

    if (psf == "Gaussian"){
      psf_k         = ProFit::profitMakeGaussianPSF(fwhm = fwhm_scaled, dim = c(25,25))
    }
    if (psf == "Moffat"){
      psf_k         = ProFit::profitCubaMoffat(fwhm = fwhm_scaled, mag = 1, con = 5, dim = c(25,25))
    }

    counts_img     = ProFit::profitBruteConv(ifu_images$counts_img, psf_k) * ifu_images$appregion
    velocity_img   = ProFit::profitBruteConv(ifu_images$velocity_img, psf_k) * ifu_images$appregion
    dispersion_img = ProFit::profitBruteConv(ifu_images$dispersion_img, psf_k) * ifu_images$appregion
    counts         = counts_img * calcregion_reff
    velocity       = velocity_img * calcregion_reff
    dispersion     = dispersion_img * calcregion_reff
    lambda         = sum(counts * radius * abs(velocity))/sum(counts * radius * (sqrt((velocity * velocity) + (dispersion * dispersion))))

    elli_x = seq(-a, a, length.out = 500)
    elli_y = (b / a) * sqrt(a^2 - elli_x^2)
    reff_elli = matrix(data=NA, nrow = 1000, ncol = 2)
    reff_elli[1:1000,1] = c(elli_x, rev(elli_x)) + rep(sbin/2, 1000)
    reff_elli[1:1000,2] = c(elli_y, rev(elli_y)*-1) + rep(sbin/2, 1000)

    if (max(reff_elli)>sbin){cat("WARNING: reff > aperture, the value of $obs_lambdar produced will not be the true value evaluated at reff.")}

    output = list("obs_lambdar"    = lambda,
                  "counts_img"     = counts_img,
                  "velocity_img"   = velocity_img,
                  "dispersion_img" = dispersion_img,
                  "reff_ellipse"      = reff_elli)
  }

  return(output)

}
