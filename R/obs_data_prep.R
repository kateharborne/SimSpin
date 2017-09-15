# Kate Harborne (last edit - 13/09/2017)
#'Prepare simulation data for observational kinematic analysis.
#'
#'The purpose of this function is to calculate the factors necessary for constructing an IFU observation data cube of a simulated galaxy at a user
#'specified inclination.
#'
#'@param filename The Gadget file containing the particle information of the galaxy to be analysed.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the simulation, 1 - gas, 2 - dark matter, 3 - disc, 4 - bulge, 5 - stars, 6 - boundary.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param z The galaxy redshift.
#'@param fov The field of view of the telescope in arcseconds.
#'@param ap_shape The shape of the field of view, with options "circular", "square" or "hexagonal".
#'@param fibres The number of observing fibres active in the telescope.
#'@param pivot_wvl The central filter wavelength used for the observation, given in angstroms.
#'@param pixel_res The corresponding velocity pixel resolution associated with a given telescope filter in angstroms.
#'@param inc_deg The inclination at which to observe the galaxy in degrees.
#'@return Returned is a list that contains a data frame of the observed particle information (\code{$galaxy_obs} containing the galaxy particle info
#' \code{$ID}, \code{$x}, \code{$z_obs}, \code{$r_obs}, \code{$vy_obs} and observed luminosity \code{$L}), and then numerical factors including the number of spatial bins (\code{$sbin}),
#' the size of those spatial bins in kpc (\code{$sbinsize}) and arcseconds (\code{pixsize}), the number of velocity bins (\code{$vbin}), the size of those velocity bins in km/s
#' (\code{$vbinsize}), the angular size in kpc/arcecond at the provided redshift (\code{$angular_size}) and the luminosity distance of the galaxy given in Mpc (\code{lum_dist}).
#'@examples
#' \dontrun{
#' obs_data_prep(filename  = "path/to/some/snapshot_XXX",
#'               r200      = 200,
#'               z         = 0.1,
#'               fov       = 15,
#'               ap_shape  = "circular",
#'               fibres    = 793,
#'               pivot_wvl = 4500,
#'               pixel_res = 1.04,
#'               inc_deg   = 0)
#'
#' obs_data_prep(filename  = "path/to/some/snapshot_XXX",
#'               ptype     = c(3,4),
#'               r200      = 200,
#'               z         = 0.1,
#'               fov       = 15,
#'               ap_shape  = "hexagonal",
#'               fibres    = 793,
#'               pivot_wvl = 4500,
#'               pixel_res = 1.04,
#'               inc_deg   = 0)
#' }
#'
obs_data_prep = function(filename, ptype=NA, r200, z, fov, ap_shape, fibres, pivot_wvl, pixel_res, inc_deg){

  set.seed(42);

  galaxy_data  = snapshot::snapread(filename) # reading in the snapshot data into large list
  galaxy_data$part$part_type = rep(0, nrow(galaxy_data$part))
  p = seq(1,6) # all possible particle values
  ppart = galaxy_data$head$Npart[which(galaxy_data$head$Nall[p] != 0)] # present particles

  for (i in 1:length(ppart)){
    if (i == 1){
      galaxy_data$part$part_type[1:ppart[i]] =  which(galaxy_data$head$Nall[p] != 0)[i]
    } else {
      galaxy_data$part$part_type[ppart[i-1]+1:ppart[i]] = which(galaxy_data$head$Nall[p] != 0)[i]
    }
  } # labelling the data frame with particle types

  if (is.na(ptype[1])){ptype = which(galaxy_data$head$Nall[p] != 0)} # for all particles, leave ptype = NA
  if (0 %in% galaxy_data$head$Nall[ptype]){
    cat("Particles of ptype =", which(galaxy_data$head$Nall[ptype] %in% 0), " are missing in this model. \n")
    stop("Npart Error")
  } # error returned if a requested ptype is not present in the simulation

  galaxy_data$part = galaxy_data$part[galaxy_data$part$part_type %in% ptype,] # leaving only particles of requested ptype
  galaxy_data$head$Npart[p[!p %in% ptype]] = as.integer(0)
  galaxy_data$head$Nall[p[!p %in% ptype]] = as.integer(0) # removing record of the removed particles in the header

  inc_rad      = inc_deg * (pi / 180) # the galaxy inclination in radians
  ang_size     = celestial::cosdistAngScale(z, ref="Planck15") # the angular size of the galaxy given the redshift, kpc/arcsecond
  ap_size      = ang_size * fov # aperture size of the telescope in kpc

  if (ap_shape == "circular"){
    sbin       = floor(sqrt(fibres * (4 / pi))) # number of spatial bins in the x- and y/z_obs- directions
  }
  if (ap_shape == "square"){
    sbin       = floor(sqrt(fibres))
  }
  if (ap_shape == "hexagonal"){
    sbin       = floor(2 * sqrt((fibres * 2) / (3 * sqrt(3))))
  }
  sbinsize     = (ap_size * 2) / sbin # kpc per bin
  vbinsize     = (pixel_res / pivot_wvl) * (3e8 / 1e3)
  pixsize      = sbinsize / ang_size # arcscond per bin

  appregion    = matrix(data = 0, ncol = sbin, nrow = sbin)
  xcentre = sbin/2 + 0.5
  ycentre = sbin/2 + 0.5
  if (ap_shape == "circular"){
    for (x in 1:sbin){
      for (y in 1:sbin){
        xx = x - xcentre
        yy = y - ycentre
        rr = sqrt(xx^2 + yy^2)
        if (rr <= fov+0.5){
          appregion[x,y] = 1
        }
      }
    }
  }
  if (ap_shape == "square"){
    appregion = matrix(data = 1, ncol = sbin, nrow = sbin)
  }
  if (ap_shape == "hexagonal"){
    for (x in 1:sbin){
      for (y in 1:sbin){
        xx = x - xcentre
        yy = y - ycentre
        rr = (2 * (sbin / 4) * (sbin * sqrt(3) / 4)) - ((sbin / 4) ) * abs(yy) - ((sbin * sqrt(3) / 4)) * abs(xx)
        if ((rr >= 0) && (abs(xx) < sbin/2) && (abs(yy) < (sbin  * sqrt(3) / 4))){
          appregion[x,y] = 1
        }
      }
    }
  }

  galaxy_df   = obs_galaxy(galaxy_data$part, centre=TRUE, inc_rad) # extracting the position and line of sight of galaxy particles for a given inclination
  galaxy_df   = galaxy_df[(galaxy_df$r < r200),]

  if (ap_shape == "circular"){
    galaxy_cdf  = galaxy_df[galaxy_df$r_obs < ap_size,]
  }
  if (ap_shape == "square"){
    galaxy_cdf  = galaxy_df[abs(galaxy_df$x) < (sbin / 2) * sbinsize,]
    galaxy_cdf  = galaxy_cdf[abs(galaxy_cdf$z_obs) < (sbin / 2) * sbinsize,]
  }
  if (ap_shape == "hexagonal"){
    galaxy_cdf  = galaxy_df[abs(galaxy_df$x) < (sbin / 2) * sbinsize,]
    galaxy_cdf  = galaxy_cdf[abs(galaxy_cdf$z_obs) < (sbin  * sqrt(3) / 4) * sbinsize,]
    dotprod     = (2 * (sbin / 4) * sbinsize * (sbin * sqrt(3) / 4) * sbinsize) - ((sbin / 4) * sbinsize) * abs(galaxy_cdf$z_obs) - ((sbin * sqrt(3) / 4) * sbinsize) * abs(galaxy_cdf$x)
    galaxy_cdf  = galaxy_cdf[dotprod >= 0,]
  } # cutting the number of particles to only contain those within the telescope aperture

  galaxy_cdf$vy_obs = galaxy_cdf$vy_obs + rnorm(galaxy_cdf$vy_obs, mean = 0, sd = vbinsize) # adding noise to the velocities to constraints due to velocity resolution

  output = list("galaxy_obs"  = galaxy_cdf,
                "sbin"        = sbin,
                "sbinsize"    = sbinsize,
                "pixsize"     = pixsize,
                "angular_size"= ang_size,
                "vbinsize"    = vbinsize,
                "appregion"   = appregion)

  return(output)
}
