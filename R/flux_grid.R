# Kate Harborne (last edit - 15/11/2019)
#'Generating the fluxes for each element of the IFU data-cube.
#'
#'The purpose of this function is to construct the mock flux values within each cell of the IFU
#' cube. It accepts output parameters from \code{obs_data_prep()} and returns a 3D array containing
#' the flux at each cell position due to contributions from the particles. If particle ages and
#' metallicities are supplied, an SED is generated in each cell using ProSpect. Else, the
#' luminosity in each cell is converted to flux.
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param filter If Age/Metallicity is supplied, the filter within which the SED is generated.
#'Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
#'@return Returns a list containing:
#'\item{\code{$flux_grid}}{A 3D array of flux in Jansky.}
#'\item{\code{$image_grid}}{A list containing the indices of particles in each spaxel position.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' flux_img    = flux_grid(obs_data  = data)
#'


flux_grid = function(obs_data, filter="g"){

  numCores = parallel::detectCores()
  z = obs_data$z
  galIDs = as.data.frame(obs_data$galaxy_obs$binn) # cell that each particle sits in
  sbin = obs_data$sbin # dimensions of the final image = sbin*sbin
  vbin = obs_data$vbin # depth of cube
  cube_size = sbin*sbin*vbin
  image_grid = vector("list", cube_size)
  image_grid = parallel::mclapply(seq(1,cube_size), function(x) which(galIDs == x),
                     mc.cores=numCores) # list of particles that exist in each spaxel
  lengths_grid = lapply(image_grid, length) # number of particles in each spaxel

  if ("Metallicity" %in% names(obs_data$galaxy_obs) & "Age" %in% names(obs_data$galaxy_obs)){ # if SSP is required
    if (filter=="g"){tempfilt=list(ProSpect::filt_g_SDSS)} # loading the filter
    if (filter=="r"){tempfilt=list(ProSpect::filt_r_SDSS)}
        flux = parallel::mclapply(seq(1, cube_size), .assign_flux, obs__data=obs_data, image__grid=image_grid,
                              lengths__grid=lengths_grid, temp_filt=tempfilt, redshift=z, mc.cores=numCores)

  } else if ("Lum" %in% names(obs_data$galaxy_obs)){ # if just M2L conversion
    flux = parallel::mclapply(seq(1, cube_size), .masslum2flux, obs__data=obs_data, image__grid=image_grid,
                              lengths__grid=lengths_grid, redshift=z, mc.cores=numCores)
  }
  flux = as.numeric(flux)
  ap_region = obs_data$ap_region; ap_region[ap_region == 0] = NA
  outside_region = which(is.na(flux*as.vector(ap_region)))

  if (length(outside_region) != 0){
    flux[outside_region] = 0 # and any cells that are outside the aperture
    for (mm in 1:length(outside_region)){
      image_grid[[outside_region[mm]]] = integer(0)
    } # remove any entries in image grid that are outside the aperture
  }

  flux = array(data = flux, dim=c(sbin,sbin,vbin))

  output = list("flux_grid"  = flux,
                "image_grid" = image_grid)
  return(output)

}
