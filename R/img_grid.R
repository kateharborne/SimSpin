# Kate Harborne (last edit - 10/09/2019)
#'Creating a mock flux image.
#'
#'The purpose of this function is to construct a mock flux image. It accepts output parameters from
#' \code{obs_data_prep()} and returns a 2D array containing the flux at each spaxel position due to
#' contributions from the particles. If particle ages and metallicities are supplied, an SED is
#' generated in each spaxel using ProSpect. Else, the luminosity in each cell is converted to flux.
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param z The projected redshift of the observation, as used in the \code{obs_data_prep()}
#' function.
#'@param filter If Age/Metallicity is supplied, the filter within which the SED is generated.
#'Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
#'@param threshold The magnitude limit of the observation in AB mag.
#'@return Returns a list containing:
#'\item{\code{$flux}}{A 2D array of flux in Jansky.}
#'\item{\code{$image_grid}}{A list containing the indices of particles in each spaxel position.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' flux_img    = img_grid(obs_data  = data)
#'


img_grid = function(obs_data, z=0.05, filter="g", threshold=25){

  galIDs = as.data.frame(obs_data$galaxy_obs$binn) # spaxel that each particle sits in
  sbin = obs_data$sbin # dimensions of the final image = sbin*sbin
  image_grid = lapply(seq(1,sbin*sbin), function(x) which(galIDs == x)) # list of particles that exist in each spaxel
  lengths_grid = lapply(image_grid, length) # number of particles in each spaxel
  numCores = parallel::detectCores()

  if ("Metallicity" %in% names(obs_data$galaxy_obs) & "Age" %in% names(obs_data$galaxy_obs)){ # if SSP is required
    if (filter=="g"){tempfilt=list(ProSpect::filt_g_SDSS)} # loading the filter
    if (filter=="r"){tempfilt=list(ProSpect::filt_r_SDSS)}
    flux = parallel::mclapply(seq(1, sbin*sbin), .assign_flux, obs__data=obs_data, image__grid=image_grid,
                              lengths__grid=lengths_grid, temp_filt=tempfilt, redshift=z, mc.cores=numCores)

  } else if ("Lum" %in% names(obs_data$galaxy_obs)){ # if just M2L conversion
    flux = parallel::mclapply(seq(1, sbin*sbin), .masslum2flux, obs__data=obs_data, image__grid=image_grid,
                              lengths__grid=lengths_grid, redshift=z, mc.cores=numCores)
  }
  flux = as.numeric(flux)
  threshold_flux  = ProSpect::magAB2Jansky(threshold)
  below_threshold = which(flux < threshold_flux)
  ap_region = obs_data$ap_region; ap_region[ap_region == 0] = NA
  outside_region = which(is.na(flux*as.vector(ap_region)))

  flux[below_threshold] = 0 # remove any cells that aren't above the threshold
  flux[outside_region] = 0 # and any cells that are outside the aperture
  for (nn in 1:length(below_threshold)){
    image_grid[[below_threshold[nn]]] = integer(0)
  } # remove any entries in image grid that are below the flux threshold
  for (mm in 1:length(outside_region)){
    image_grid[[outside_region[mm]]] = integer(0)
  } # remove any entries in image grid that are outside the aperture

  flux = matrix(data = flux, nrow=sbin, ncol=sbin, byrow=FALSE)

  output = list("flux_img"   = flux,
                "image_grid" = image_grid)
  return(output)

}
