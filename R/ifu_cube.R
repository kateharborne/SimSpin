# Kate Harborne (last edit - 15/11/2019)
#'Creating a mock IFU kinematic data cube.
#'
#'The purpose of this function is to construct an IFU data cube. It accepts output parameters from
#' \code{obs_data_prep()} and \code{flux_grid()} and returns a 3D array containing the spatial and
#'  velocity corrdinates of each particle in the galaxy simulation in order to mimic the
#'  arrangement of an IFU data cube.
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param flux_data The list output from the \code{flux_grid()} function.
#'@param threshold The magnitude limit of the observation in AB mag.
#'@return Returns a list that contains:
#' \item{\code{$cube}}{A mock IFU data cube as required for calculating the observed spin parameter}
#' \item{\code{$ap_region}}{An image that describes the shape of the aperture such that any further
#'  convolutions that are applied to mimic beam smearing or seeing can be trimmed to the
#'  appropriate aperture shape.}
#' \item{\code{$xbin_labels}}{Bin labels for the x-spatial dimension.}
#' \item{\code{$ybin_labels}}{Bin labels for the y-spatial dimension.}
#' \item{\code{$vbin_labels}}{Bin labels for the velocity dimension.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' fluxes      = flux_grid(obs_data = data)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#'

ifu_cube = function(obs_data, flux_data, threshold=25) {

  numCores     = parallel::detectCores()
  image_grid   = flux_data$image_grid
  lengths_grid = lapply(image_grid, length) # how many particles in each cell?

  galaxy_obs   = obs_data$galaxy_obs

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
  lsf_size     = obs_data$lsf_size
  cube         = array(data = 0, dim = c(sbin, sbin, vbin))
  ap_cube      = array(data = obs_data$ap_region, dim = c(sbin,sbin,vbin))

  for (cell in seq(1, sbin*sbin*vbin)){ # for each cell in cube
    if (lengths_grid[[cell]] > 0){
      coord = c(cell%%sbin, (cell%%(sbin*sbin)%/%sbin + 1), cell%/%(sbin*sbin)+1)
      if (cell%%sbin == 0 & cell%%(sbin*sbin) == 0){coord = c(sbin, sbin, cell%/%(sbin*sbin))} # case at the end of row & column
      if (cell%%sbin == 0 & cell%%(sbin*sbin) != 0){coord = c(sbin, (cell%%(sbin*sbin)%/%sbin), cell%/%(sbin*sbin)+1)} # case at the end of row
      for (j in 1:lengths_grid[[cell]]){
        cell_mass = sum(obs_data$galaxy_obs$Mass[image_grid[[cell]]])
        cell_flux = flux_data$flux_grid[coord[1], coord[2], coord[3]]*(obs_data$galaxy_obs$Mass[image_grid[[cell]]][j]/cell_mass)
        cube[coord[1], coord[2],] = cube[coord[1], coord[2],] +
          diff((cell_flux*
               pnorm(vseq, mean=obs_data$galaxy_obs$vy_obs[image_grid[[cell]][j]], sd=lsf_size)))
          # adding the "gaussians" of each particle to the velocity bins
        }
      }
  }
  threshold_flux = ProSpect::magAB2Jansky(threshold)

  for (i in 1:vbin){
    below_threshold = which(cube[,,i]<threshold_flux)
    cube[,,i][below_threshold] = 0
  }

  output = list("cube"        = cube,
                "ap_region"   = ap_cube,
                "xbin_labels" = xbin_ls,
                "zbin_labels" = zbin_ls,
                "vbin_labels" = vbin_ls)

  return(output)

}
