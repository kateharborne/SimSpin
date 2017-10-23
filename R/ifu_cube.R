# Kate Harborne (last edit - 13/09/2017)
#'Creating a set of mock IFU observational images.
#'
#'The purpose of this function is to construct an IFU data cube. It accepts output parameters from \code{obs_data_prep()} and returns
#' a 3D array containing the spatial and velocity corrdinates of each particle in the galaxy simulation in order to mimic the arrangement
#' of an IFU data cube.
#'
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param threshold The minimum number of counts in the image.
#'@return Returns a list containing a mock IFU data cube including a counts image (\code{$counts_img}), a velocity image (\code{$velocity_img})
#' and a dispersion image (\code{$dispersion_img}) as required for calculating the observed spin parameter. Also contained is an image that
#' describes the shape of the apperture (\code{$appregion}) such that any further convolutions that are applied to mimic beam smearing or seeing
#' can be trimmed to the appropriate apperture shape.
#'@examples
#' \dontrun{
#' data = obs_data_prep()
#'
#' ifu_cube(obs_data  = data,
#'         threshold = 0)
#'
#' ifu_cube(obs_data  = data,
#'         threshold = 20)
#' }
#'

ifu_cube = function(obs_data, threshold) {

  galaxy_obs      = obs_data$galaxy_obs # data to populate the aperture
  sbin            = obs_data$sbin # number of spatial bins along the x and z axes
  sbin_breaks     = seq(-(sbin * obs_data$sbinsize) / 2, (sbin * obs_data$sbinsize) / 2, by=obs_data$sbinsize)
  vbin            = obs_data$vbin
  vbinsize        = obs_data$vbinsize
  vseq            = seq(-(vbin * vbinsize) / 2, (vbin * vbinsize) / 2, by=vbinsize)
  lsf_size        = obs_data$lsf_size
  galaxy_obs$binx = cut(galaxy_obs$x, breaks=sbin_breaks, labels=F)
  galaxy_obs$binz = cut(galaxy_obs$z_obs, breaks=sbin_breaks, labels=F)
  galaxy_obs$binn = galaxy_obs$binx + (sbin * galaxy_obs$binz) - sbin # labelling data based on it's position in the aperture

  xbins           = levels(cut(galaxy_obs$x, breaks=sbin_breaks))
  zbins           = levels(cut(galaxy_obs$z_obs, breaks=sbin_breaks))
  vbins           = levels(cut(galaxy_obs$vy_obs, breaks=vseq))
  xbin_labels     = expand.grid(matrix(data = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", xbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", xbins)))), nrow = sbin, ncol = sbin))
  zbin_labels     = expand.grid(t(matrix(data = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", zbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", zbins)))), nrow = sbin, ncol = sbin)))
  xbin_ls         = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", xbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", xbins))))
  zbin_ls         = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", zbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", zbins))))
  vbin_ls         = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", vbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", vbins))))

  counts      = with(galaxy_obs, as.integer(binn))
  counts_df   = as.data.frame(table(counts))
  counts_flat = data.frame("counts"=rep(0,sbin*sbin)) # a 1D array containing the number of particles in each spatial bin
  counts_flat[as.integer(as.vector(counts_df$counts)),] = as.integer(as.vector(counts_df$Freq))
  nr          = nrow(counts_flat)

  xcen       = .meanwt(xbin_labels, counts_flat) # calculating the centre of the galaxy
  zcen       = .meanwt(zbin_labels, counts_flat)
  sx         = sqrt(.varwt(xbin_labels, counts_flat, xcen)) # calculating the standard deviation along each direction
  sz         = sqrt(.varwt(zbin_labels, counts_flat, zcen))
  sxz        = .covarwt(xbin_labels, zbin_labels, counts_flat, xcen, zcen) # calculating the covariance
  temprad    = .cov2eigval(sx, sz, sxz) # solving for the eigenvalues to give the major and minor axes lengths
  major      = sqrt(abs(temprad$hi))
  minor      = sqrt(abs(temprad$lo))
  axis_ratio = data.frame("a" = major, "b" = minor)

  vel_obs = data.frame("binn" = galaxy_obs$binn, "vz" = galaxy_obs$vy_obs)
  vel_obs = vel_obs[order(vel_obs$binn),]
  cum_counts = cumsum(counts_flat)
  cell_counts = array(data = 0, dim = 1)
  cube = array(data = 0, dim = c(sbin, sbin, vbin))

  for (i in 1:nr){
    if (counts_flat$counts[i] == 0){
      i = i+1 # if there are no counts in the cell, add one to the sequence and start repeat of loop again
    } else {
      coord = c((i%/%sbin + 1), i%%sbin) # determining coordinates to fill
      if (i%%sbin == 0){coord = c((i%/%sbin), sbin)} # if at the end of a row, setting the coordinate correctly
      cell_counts = vel_obs$vz[(cum_counts$counts[i-1] + 1): cum_counts$counts[i]] # all the velocities in that spatial bin
      if (length(cell_counts) < threshold){
        cube[coord[2], coord[1],] = 0 # if there are too few counts in that cell, set the cube cell to zero and skipping to the next cell
      } else {
        for (j in 1:length(cell_counts)){
          cube[coord[2], coord[1],] = cube[coord[2], coord[1],] + diff(pnorm(vseq, mean=cell_counts[j], sd=lsf_size)) # adding the "gaussians" of each particle to the velocity bins
        }
      }
      cell_counts = array(data = 0, dim = 1) # reset for the next iteration
    }
  }

  appcube = array(data = obs_data$appregion, dim = c(sbin,sbin,vbin))
  cube = cube * appcube

  output = list("cube" = cube, "vbin_labels" = vbin_ls, "xbin_labels" = xbin_ls, "zbin_labels" = zbin_ls, "appregion" = appcube, "axis_ratio" = axis_ratio)

  return(output)

}
