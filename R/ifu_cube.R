# Kate Harborne (last edit - 15/11/2019)
#'Creating a mock IFU kinematic data cube.
#'
#'The purpose of this function is to construct an IFU data cube. It accepts output parameters from
#' \code{obs_data_prep()} and \code{flux_grid()} and returns a 3D array containing the spatial and
#'  velocity corrdinates of each particle in the galaxy simulation in order to mimic the
#'  arrangement of an IFU data cube.
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param flux_data The list output from the \code{flux_grid()} function.
#'@return Returns a list that contains:
#' \item{\code{$cube}}{A mock IFU data cube as required for calculating the observed spin parameter}
#' \item{\code{$ap_region}}{An image that describes the shape of the aperture such that any further
#'  convolutions that are applied to mimic beam smearing or seeing can be trimmed to the
#'  appropriate aperture shape.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' fluxes      = flux_grid(obs_data = data)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#'

ifu_cube = function(obs_data, flux_data) {

  numCores     = parallel::detectCores()
  image_grid   = flux_data$image_grid
  lengths_grid = lapply(image_grid, length) # how many particles in each cell?

  galaxy_obs   = obs_data$galaxy_obs
  sbin         = obs_data$sbin
  sbinsize     = obs_data$sbinsize
  vbin         = obs_data$vbin
  vbinsize     = obs_data$vbinsize
  vseq         = seq(-(vbin * vbinsize) / 2, (vbin * vbinsize) / 2, by=vbinsize) # velocity bin break positions
  lsf_size     = obs_data$lsf_size
  cube         = array(data = 0, dim = c(sbin, sbin, vbin))
  ap_cube      = array(data = obs_data$ap_region, dim = c(sbin,sbin,vbin))

  for (cell in seq(1, sbin*sbin*vbin)){ # for each cell in cube
    if (lengths_grid[[cell]] > 1){
      coord = c(cell%%sbin, (cell%%(sbin*sbin)%/%sbin + 1), cell%/%(sbin*sbin)+1)
      if (cell%%sbin == 0 & cell%%(sbin*sbin) == 0){coord = c(sbin, sbin, cell%/%(sbin*sbin))} # case at the end of row & column
      if (cell%%sbin == 0 & cell%%(sbin*sbin) != 0){coord = c(sbin, (cell%%(sbin*sbin)%/%sbin), cell%/%(sbin*sbin)+1)} # case at the end of row
      for (j in 1:lengths_grid[[cell]]){
          cube[coord[1], coord[2],] = cube[coord[1], coord[2],] +
            diff(flux_data$flux_grid[coord[1], coord[2], coord[3]] *
                   pnorm(vseq, mean=obs_data$galaxy_obs$vy_obs[image_grid[[cell]][j]], sd=lsf_size))
          # adding the "gaussians" of each particle to the velocity bins
        }
      }
  }

  output = list("cube"        = cube,
                "ap_region"   = ap_cube)

  return(output)

}
