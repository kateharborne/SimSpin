# Kate Harborne (last edit - 23/04/2018)
#'Creating images from a mock IFU kinematic data cube.
#'
#'The purpose of this function is to construct the observational images from a mock IFU data cube.
#' It accepts output parameters from \code{obs_data_prep()} and \code{ifu_cube()} and returns three
#' images (flux, line-of-sight velocity and line-of-sight velocity dispersion).
#'
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param ifu_datacube The list output from the function \code{\link{ifu_cube}} containing the mock
#' IFU cube and the apperture region image (\code{$appregion}).
#'@param addSky A boolean to specify whether to add sky noise to the output images. Default is
#' FALSE. If TRUE, further parameters including \code{mag_threshold} and \code{mag_zero} described
#' below.
#'@param mag_zero The magnitude zero point with regards to the mangitude system being used (e.g.
#' AB or Vega).
#'@param threshold The magnitude limit of the observation in AB mag.
#'@param pixel_sscale The corresponding spatial pixel scale associated with a given telescope
#' output in arcseconds.
#'@return Returns a list containing:
#'\item{\code{$flux_img}}{The flux image produced from flattening the cube along the velocity
#'domain.}
#'\item{\code{$velocity_img}}{The line-of-sight velocity image produced from taking the flux-
#'weighted velocities along the velocity domain.}
#'\item{\code{$dispersion_img}}{The line-of-sight velocity dispersion image produced from taking
#'the standard deviation of the flux-weighted velocities along the velocity domain.}
#'\item{\code{$axis_ratio}}{A list containing the semi-major (\code{$axis_ratio$a}) and semi-minor
#' (\code{$axis_ratio$b}) axis lengths of the galaxy image in kpc.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' fluxes      = flux_grid(obs_data = data, multi_thread=FALSE)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#' images      = obs_imgs(obs_data = data, ifu_datacube = cube)
#'

obs_imgs = function(obs_data, ifu_datacube, threshold=25, addSky=FALSE,
                    mag_zero=8.9, pixel_sscale=0.5){

  sbin = obs_data$sbin # dimensions of the final image = sbin*sbin
  vbin = obs_data$vbin # depth of cube
  xbin_ls      = ifu_datacube$xbin_labels
  zbin_ls      = ifu_datacube$zbin_labels
  vbin_ls      = ifu_datacube$vbin_labels

  counts_img = apply(ifu_datacube$cube, c(1,2), sum)
  velocity_img   = matrix(data=0, nrow=sbin, ncol=sbin)
  dispersion_img = matrix(data=0, nrow=sbin, ncol=sbin)
  for (c in 1:sbin){
    for (d in 1:sbin){
      velocity_img[c,d]   = .meanwt(vbin_ls, ifu_datacube$cube[c,d,])
      dispersion_img[c,d] = sqrt(.varwt(vbin_ls, ifu_datacube$cube[c,d,], velocity_img[c,d]))
    }
  }

  if (addSky){
    skyRMS = ProFound::profoundSB2Flux(threshold, mag_zero, pixel_sscale)
    noise = rnorm(dim(counts_img)[1]^2, sd=skyRMS)
    counts_img = counts_img+noise
    velocity_img = velocity_img+noise
    dispersion_img = dispersion_img+noise
  }

  threshold_flux = ProSpect::magAB2Jansky(threshold)

  below_threshold = which(counts_img<threshold_flux)
  counts_img[below_threshold] = 0;
  velocity_img[(is.na(velocity_img))] = 0; velocity_img[below_threshold] = 0
  dispersion_img[(is.na(dispersion_img))] = 0; dispersion_img[below_threshold] = 0

  counts_img = counts_img*obs_data$ap_region
  velocity_img = velocity_img*obs_data$ap_region
  dispersion_img = dispersion_img*obs_data$ap_region

  xbin_labels     = expand.grid(matrix(data = xbin_ls,
                                       nrow = sbin, ncol = sbin))
  zbin_labels     = expand.grid(t(matrix(data = zbin_ls,
                                         nrow = sbin, ncol = sbin)))
  #spatial coordinates in each spaxel in a column for ellipticity calculation

  xcen       = .meanwt(xbin_labels, as.vector(counts_img)) # calculating the centre of the galaxy
  zcen       = .meanwt(zbin_labels, as.vector(counts_img))
  sx         = sqrt(.varwt(xbin_labels, as.vector(counts_img), xcen)) # standard deviation along
  sz         = sqrt(.varwt(zbin_labels, as.vector(counts_img), zcen)) #   each direction
  sxz        = .covarwt(xbin_labels, zbin_labels, as.vector(counts_img), xcen, zcen) # covariance
  temprad    = .cov2eigval(sx, sz, sxz) # solving for the eigenvalues to give the major and
  major      = sqrt(abs(temprad$hi))    #    minor axes lengths
  minor      = sqrt(abs(temprad$lo))
  tempgrad   = .cov2eigvec(sx, sz, sxz)
  ang        = (180-atan(tempgrad)*180/pi)%%180
  axis_ratio = data.frame("a" = major, "b" = minor, "ang" = ang) # axis ratio output, kpc

  output = list("flux_img"       = counts_img,
                "velocity_img"   = velocity_img,
                "dispersion_img" = dispersion_img,
                "axis_ratio"     = axis_ratio)

  return(output)

}
