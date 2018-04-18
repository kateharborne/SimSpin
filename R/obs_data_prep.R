# Kate Harborne (last edit - 16/04/2018)
#'Prepare simulation data for observational kinematic analysis.
#'
#'The purpose of this function is to calculate the properties necessary for constructing an IFU observation data cube of a simulated galaxy at a user
#' specified inclination and redshift.
#'
#'@param filename The Gadget output file containing the particle information of the galaxy to be analysed.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the simulation, 1 - gas, 2 - dark matter, 3 - disc,
#'4 - bulge, 5 - stars, 6 - boundary.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param z The galaxy redshift.
#'@param fov The field of view of the IFU, diameter in arcseconds.
#'@param ap_shape The shape of the field of view, with options "circular", "square" or "hexagonal".
#'@param central_wvl The central filter wavelength used for the observation, given in angstroms.
#'@param lsf_fwhm The line spread function full-width half-max, given in angstroms.
#'@param pixel_sscale The corresponding spatial pixel scale associated with a given telescope output in arcseconds.
#'@param pixel_vscale The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.
#'@param inc_deg The inclination at which to observe the galaxy in degrees.
#'@param m2l_disc The mass-to-light ratio of the disc component in solar units.
#'@param m2l_bulge The mass-to-light ratio of the bulge component in solar units.
#'@return Returned is a list that contains a data frame of the observed particle information (\code{$galaxy_obs} (which contains the galaxy particle
#' info from the GADGET file, with added coordinates \code{$z_obs}, \code{$r_obs}, and \code{$vy_obs} and particle fluxes, \code{$flux}, in units of
#' 1e-16 erg s-1 cm-2) and numerical factors including the number of spatial bins (\code{$sbin}), the size of those spatial bins in kpc and arcseconds
#' (\code{$sbinsize} and \code{pixsize}), the number of velocity bins (\code{$vbin}), the size of those velocity bins in km/s (\code{$vbinsize}), the
#' gaussian standard deviation of the line spread function in km/s (\code{lsf_size}) and the angular size of the galaxy in kpc/arcecond at the provided
#' redshift (\code{$angular_size}).
#'@examples
#' \dontrun{
#' obs_data_prep(filename     = "path/to/some/snapshot_XXX",
#'               r200         = 200,
#'               z            = 0.1,
#'               fov          = 15,
#'               ap_shape     = "circular",
#'               central_wvl  = 4800,
#'               lsf_fwhm     = 2.65,
#'               pixel_sscale = 0.5,
#'               pixel_vscale = 1.04,
#'               inc_deg      = 0,
#'               m2l_disc     = 2,
#'               m2l_bulge    = 1)
#'
#' obs_data_prep(filename     = "path/to/some/snapshot_XXX",
#'               ptype        = c(3,4),
#'               r200         = 200,
#'               z            = 0.1,
#'               fov          = 15,
#'               ap_shape     = "hexagonal",
#'               central_wvl  = 4800,
#'               lsf_fwhm     = 2.65,
#'               pixel_sscale = 0.5,
#'               pixel_vscale = 1.04,
#'               inc_deg      = 0,
#'               m2l_disc     = 1,
#'               m2l_bulge    = 1)
#' }
#'
obs_data_prep = function(filename, ptype=NA, r200=200, z, fov, ap_shape, central_wvl, lsf_fwhm, pixel_sscale, pixel_vscale, inc_deg, m2l_disc, m2l_bulge){

  set.seed(42);
  sol_lum = 3.827e33;                                                            # solar luminosity in erg s-1
  galaxy_data  = snapshot::snapread(filename)                                    # reading in the snapshot data
  galaxy_data$part$part_type = rep(0, nrow(galaxy_data$part))                    # add a "particle type" column
  p = seq(1,6)                                                                   # all possible particle values
  ppart = cumsum(galaxy_data$head$Npart[which(galaxy_data$head$Nall[p] != 0)])   # present particles
  for (i in 1:length(ppart)){
    if (i == 1){
      galaxy_data$part[1:as.integer(ppart[i]),]$part_type =  which(galaxy_data$head$Nall[p] != 0)[i]
    } else {
      galaxy_data$part[as.integer(ppart[i-1]+1):as.integer(ppart[i]),]$part_type = which(galaxy_data$head$Nall[p] != 0)[i]
    }
  }
                                                                                 # labelling the data frame with particle types
  if (is.na(ptype[1])){ptype = which(galaxy_data$head$Nall[p] != 0)}             # for all particles, leave ptype = NA
  if (0 %in% galaxy_data$head$Nall[ptype]){
    cat("Particles of ptype =", which(galaxy_data$head$Nall[ptype] %in% 0), " are missing in this model. \n")
    stop("Npart Error")
  }
                                                                                 # error returned if a requested ptype is not present
  galaxy_data$part = galaxy_data$part[galaxy_data$part$part_type %in% ptype,]    # leaving only particles of requested ptype
  galaxy_data$head$Npart[p[!p %in% ptype]] = as.integer(0)
  galaxy_data$head$Nall[p[!p %in% ptype]] = as.integer(0)                        # removing record of the removed particles in the header

  inc_rad      = inc_deg * (pi / 180)                                            # galaxy inclination in radians
  ang_size     = celestial::cosdistAngScale(z, ref="Planck15")                   # angular size of the galaxy given the redshift, kpc/"
  lum_dist     = celestial::cosdistLumDist(z, ref="Planck15")                    # the luminosity distance of the galaxy, Mpc
  ap_size      = ang_size * fov                                                  # aperture diameter size of the telescope, kpc
  sbin         = floor(fov / pixel_sscale)                                       # number of spatial bins in the x- and y/z_obs- directions
  sbinsize     = ap_size / sbin                                                  # kpc per bin
  vbinsize     = (pixel_vscale / central_wvl) * (3e8 / 1e3)                      # km/s per velocity bin
  lsf_size     = ((lsf_fwhm / central_wvl) * (3e8 / 1e3)) / (2 * sqrt(2*log(2))) # velocity uncertainty (standard deviation, not FWHM) due to LSF

  appregion    = matrix(data = 0, ncol = sbin, nrow = sbin)                      # empty matrix for aperture mask
  xcentre = sbin/2 + 0.5
  ycentre = sbin/2 + 0.5
  if (ap_shape == "circular"){
    x = matrix(data = rep(seq(1,sbin), each=sbin), nrow = sbin, ncol = sbin)
    y = matrix(data = rep(seq(sbin,1), sbin), nrow = sbin, ncol = sbin)
    xx = x - xcentre
    yy = y - ycentre
    rr = sqrt(xx^2 + yy^2)
    appregion[rr<= sbin/2] = 1
  }
                                                                                 # circular apperture mask
  if (ap_shape == "square"){
    appregion = matrix(data = 1, ncol = sbin, nrow = sbin)
  }
                                                                                 # square apperture mask
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
                                                                                 # hexagonal apperture mask

  galaxy_df   = obs_galaxy(galaxy_data$part, centre=TRUE, inc_rad)               # extracting the position and line of sight data for inc
  galaxy_df   = galaxy_df[(galaxy_df$r < r200),]                                 # removing particles beyond r200

  if (ap_shape == "circular"){
    galaxy_cdf  = galaxy_df[galaxy_df$r_obs < (ap_size / 2),]
  }
                                                                                 # removing particles along LOS outside circular apperture
  if (ap_shape == "square"){
    galaxy_cdf  = galaxy_df[abs(galaxy_df$x) < (sbin / 2) * sbinsize,]
    galaxy_cdf  = galaxy_cdf[abs(galaxy_cdf$z_obs) < (sbin / 2) * sbinsize,]
  }
                                                                                 # removing particles along LOS outside square apperture
  if (ap_shape == "hexagonal"){
    galaxy_cdf  = galaxy_df[abs(galaxy_df$x) < (sbin / 2) * sbinsize,]
    galaxy_cdf  = galaxy_cdf[abs(galaxy_cdf$z_obs) < (sbin  * sqrt(3) / 4) * sbinsize,]
    dotprod     = (2 * (sbin / 4) * sbinsize * (sbin * sqrt(3) / 4) * sbinsize) - ((sbin / 4) * sbinsize) * abs(galaxy_cdf$z_obs) - ((sbin * sqrt(3) / 4) * sbinsize) * abs(galaxy_cdf$x)
    galaxy_cdf  = galaxy_cdf[dotprod >= 0,]
  }
                                                                                 # removing particles along LOS outside hexagonal apperture

  vbin = ceiling((max(galaxy_cdf$vy_obs) - min(galaxy_cdf$vy_obs)) / vbinsize)   # the number of velocity bins
  galaxy_cdf$flux = rep(0, nrow(galaxy_cdf))
  galaxy_cdf[(galaxy_cdf$part_type== 3),]$flux =
    ((galaxy_cdf[(galaxy_cdf$part_type== 3),]$Mass * 1e10 * sol_lum * (pixel_sscale^2)) / (m2l_disc)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2)
                                                                                 # calculating the flux from each disc particle
  galaxy_cdf[(galaxy_cdf$part_type== 4),]$flux =
    ((galaxy_cdf[(galaxy_cdf$part_type== 4),]$Mass * 1e10 * sol_lum * (pixel_sscale^2)) / (m2l_bulge)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2)
                                                                                 # calculating the flux from each bulge particle
                                                                                 # fluxes in units of (1e-16 erg s-1 cm-2 arcsecond-2)

  output = list("galaxy_obs"  = galaxy_cdf,
                "sbin"        = sbin,
                "sbinsize"    = sbinsize,
                "angular_size"= ang_size,
                "vbin"        = vbin,
                "vbinsize"    = vbinsize,
                "lsf_size"    = lsf_size,
                "pixel_sscale"= pixel_sscale,
                "appregion"   = appregion)

  return(output)
}
