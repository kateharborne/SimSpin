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
#'@param threshold The magnitude limit of the observation in AB mag.
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
#' fluxes      = flux_grid(obs_data = data)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#' images      = obs_images(obs_data = data, ifu_datacube = cube)
#'

obs_imgs = function(obs_data, ifu_datacube, threshold=25){

  sbin = obs_data$sbin # dimensions of the final image = sbin*sbin
  vbin = obs_data$vbin # depth of cube
  vbinsize = obs_data$vbinsize
  sbinsize = obs_data$sbinsize
  vseq         = seq(-(vbin * vbinsize) / 2, (vbin * vbinsize) / 2, by=vbinsize)
  # velocity bin break positions
  sseq         = seq(-(sbin * sbinsize) / 2, (sbin * sbinsize) / 2, by=sbinsize)
  # spatial bin break positions

  xbins        = levels(cut(obs_data$galaxy_obs$x, breaks=sseq))
  zbins        = levels(cut(obs_data$galaxy_obs$z_obs, breaks=sseq))
  # spatial/velocity bins boundaries
  vbins        = levels(cut(obs_data$galaxy_obs$vy_obs, breaks=vseq))

  xbin_ls      = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", xbins)),
                                as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", xbins))))
  zbin_ls      = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", zbins)),
                                as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", zbins))))
  vbin_ls      = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", vbins)),
                                as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", vbins))))

  threshold_flux = ProSpect::magAB2Jansky(threshold)

  counts_img = apply(ifu_datacube$cube, c(1,2), sum)

  below_threshold = which(counts_img<threshold_flux)
  counts_img[below_threshold] = 0

  velocity_img   = matrix(data=0, nrow=sbin, ncol=sbin)
  dispersion_img = matrix(data=0, nrow=sbin, ncol=sbin)
  for (c in 1:sbin){
    for (d in 1:sbin){
      velocity_img[c,d]   = .meanwt(vbin_ls, ifu_datacube$cube[c,d,])
      dispersion_img[c,d] = sqrt(.varwt(vbin_ls, ifu_datacube$cube[c,d,], velocity_img[c,d]))
    }
  }

  velocity_img[(is.na(velocity_img))] = 0; velocity_img[below_threshold] = 0
  dispersion_img[(is.na(dispersion_img))] = 0; dispersion_img[below_threshold] = 0

  xbin_labels     = expand.grid(matrix(data =
                                         rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", xbins)),
                                                        as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", xbins)))),
                                       nrow = sbin, ncol = sbin))
  zbin_labels     = expand.grid(t(matrix(data =
                                           rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", zbins)),
                                                          as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", zbins)))),
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
  axis_ratio = data.frame("a" = major, "b" = minor) # axis ratio output, kpc

  output = list("flux_img"       = counts_img,
                "velocity_img"   = velocity_img,
                "dispersion_img" = dispersion_img,
                "axis_ratio"     = axis_ratio,
                "xbin_labels"    = xbin_ls)

  return(output)

}
