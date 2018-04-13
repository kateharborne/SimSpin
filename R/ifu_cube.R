# Kate Harborne (last edit - 13/09/2017)
#'Creating a set of mock IFU observational images.
#'
#'The purpose of this function is to construct an IFU data cube. It accepts output parameters from \code{obs_data_prep()} and returns
#' a 3D array containing the spatial and velocity corrdinates of each particle in the galaxy simulation in order to mimic the arrangement
#' of an IFU data cube.
#'
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param threshold The flux threshold of the observation.
#'@return Returns a list containing a mock IFU data cube as required for calculating the observed spin parameter. Also contained is an image that
#' describes the shape of the apperture (\code{$appregion}) such that any further convolutions that are applied to mimic beam smearing or seeing
#' can be trimmed to the appropriate apperture shape. The bin labels for each dimension are supplied as \code{$xbin_labels}, \code{$ybin_labels},
#' \code{$vbin_labels}. The axis ratio of the galaxy is calculated in this fuction and output in the list as \code{$axis_ratio$a} and
#' \code{$axis_ratio$b}.
#'@examples
#' \dontrun{
#' data = obs_data_prep()
#'
#' ifu_cube(obs_data = data,
#'         threshold = 0)
#'
#' ifu_cube(obs_data = data,
#'         threshold = 20)
#' }
#'

ifu_cube = function(obs_data, threshold) {

  galaxy_obs      = obs_data$galaxy_obs                                                         # data to populate the aperture
  sbin            = obs_data$sbin                                                               # number of spatial bins along the x and z axes
  pixel_sscale    = obs_data$pixel_sscale                                                       # size of spatial bins
  sbin_breaks     = seq(-(sbin * obs_data$sbinsize) / 2,
                        (sbin * obs_data$sbinsize) / 2, by=obs_data$sbinsize)                   # spatial bin break positions
  vbin            = obs_data$vbin                                                               # number of velocity bins
  vbinsize        = obs_data$vbinsize                                                           # size of velocity bins
  vseq            = seq(-(vbin * vbinsize) / 2, (vbin * vbinsize) / 2, by=vbinsize)             # velocity bin break positions
  lsf_size        = obs_data$lsf_size                                                           # size of line spread function fwhm
  threshold_flux  = 1.3608e22 * 10^((threshold + 26.832) / -2.5) * (pixel_sscale^2)             # calculating threshold flux in units for ifu cube
  galaxy_obs$binx = cut(galaxy_obs$x, breaks=sbin_breaks, labels=F)
  galaxy_obs$binz = cut(galaxy_obs$z_obs, breaks=sbin_breaks, labels=F)
  galaxy_obs$binn = galaxy_obs$binx + (sbin * galaxy_obs$binz) - sbin                           # assigning particles to positions in apperture
  xbins           = levels(cut(galaxy_obs$x, breaks=sbin_breaks))
  zbins           = levels(cut(galaxy_obs$z_obs, breaks=sbin_breaks))                           # spatial/velocity bins boundaries required for the
  vbins           = levels(cut(galaxy_obs$vy_obs, breaks=vseq))                                 #    galaxy observation
  xbin_ls         = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", xbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", xbins))))
  zbin_ls         = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", zbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", zbins))))
  vbin_ls         = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", vbins)), as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", vbins))))
                                                                                                # calculating the centres of spatial/velocity bins
  counts      = with(galaxy_obs, as.integer(binn))                                              # preliminary spatial bin positions of each particle
  counts_df   = as.data.frame(table(counts))                                                    # numbers of particles in each occupied spatial bin
  counts_flat = data.frame("counts"=rep(0,sbin*sbin))                                           # a 1D array containing the number of particles in
  counts_flat[as.integer(as.vector(counts_df$counts)),] = as.integer(as.vector(counts_df$Freq)) #    each spatial bin
  nr          = nrow(counts_flat)                                                               # number of bins in the ifu cube
  vel_obs     = data.frame("binn" = galaxy_obs$binn,
                           "vz"   = galaxy_obs$vy_obs,
                           "flux" = galaxy_obs$flux)                                            # particle data needed for cube
  vel_obs     = vel_obs[order(vel_obs$binn),]                                                   # order by bin number
  cum_counts  = cumsum(counts_flat)                                                             # culumative number of particles up to bin nr
  cell_counts = array(data = 0, dim = 1)
  cell_flux   = array(data = 0, dim = 1)
  cube        = array(data = 0, dim = c(sbin, sbin, vbin))                                      # empty arrays to fill for ifu cube

  for (i in 1:nr){
    if (counts_flat$counts[i] == 0){                                                            # if there are no counts in the cell, add one to
      i = i+1                                                                                   #    the sequence and start loop again
    } else {
      coord = c((i%/%sbin + 1), i%%sbin)                                                        # determining coordinates to fill
      if (i%%sbin == 0){coord = c((i%/%sbin), sbin)}                                            # if at the end of a row, setting the coordinate correctly
      cell_counts = vel_obs$vz[(cum_counts$counts[i-1] + 1): cum_counts$counts[i]]              # all the velocities in that spatial bin
      cell_flux = vel_obs$flux[(cum_counts$counts[i-1] + 1): cum_counts$counts[i]]              # all the flux in that spatial bin
      if (sum(cell_flux) < threshold_flux){                                                     # if there are too few counts in that cell, set
        cube[coord[2], coord[1],] = 0                                                           #    the cell to zero and skip to the next cell
      } else {
        for (j in 1:length(cell_counts)){
          cube[coord[2], coord[1],] = cube[coord[2], coord[1],] +
            diff(cell_flux[j] * pnorm(vseq, mean=cell_counts[j], sd=lsf_size))                  # adding the "gaussians" of each particle to the velocity bins
        }
      }
      cell_counts = array(data = 0, dim = 1)                                                    # reset for the next iteration
      flux_counts = array(data = 0, dim = 1)                                                    # reset for the next iteration
    }
  }

  appcube = array(data = obs_data$appregion, dim = c(sbin,sbin,vbin))                           # reformatting the 1D apperture array into 3D cube
  cube = cube * appcube                                                                         # filtering out counts that lie outside the apperture

  ## ellipticity calculation ## -------------------------------------------------------------------------------------------------------------------------------
  xbin_labels     = expand.grid(matrix(data = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", xbins)),
                                                             as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", xbins)))), nrow = sbin, ncol = sbin))
  zbin_labels     = expand.grid(t(matrix(data = rowMeans(cbind(as.numeric(sub("\\((.+),.*", "\\1", zbins)),
                                                               as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", zbins)))), nrow = sbin, ncol = sbin)))
                                                                                                # spatial coordinates in each spaxel in a column for
                                                                                                #    elliptcity calculation
  xcen       = .meanwt(xbin_labels, counts_flat)                                                # calculating the centre of the galaxy
  zcen       = .meanwt(zbin_labels, counts_flat)
  sx         = sqrt(.varwt(xbin_labels, counts_flat, xcen))                                     # calculating standard deviation along each direction
  sz         = sqrt(.varwt(zbin_labels, counts_flat, zcen))
  sxz        = .covarwt(xbin_labels, zbin_labels, counts_flat, xcen, zcen)                      # calculating the covariance
  temprad    = .cov2eigval(sx, sz, sxz)                                                         # solving for the eigenvalues to give the major and
  major      = sqrt(abs(temprad$hi))                                                            #    minor axes lengths
  minor      = sqrt(abs(temprad$lo))
  axis_ratio = data.frame("a" = major, "b" = minor)                                             # axis ratio output, kpc
  blur_info = data.frame("xbins" = xbin_labels, "ybins" = zbin_labels)                          # dimension info for blur_cube() function

  output = list("cube"        = cube,
                "xbin_labels" = xbin_ls,
                "ybin_labels" = zbin_ls,
                "vbin_labels" = vbin_ls,
                "appregion"   = appcube,
                "axis_ratio"  = axis_ratio,
                "blur_info"   = blur_info)

  return(output)

}
