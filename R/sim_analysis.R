# Kate Harborne (last edit - 13/09/2017)
#'Kinematic analysis of simulation data.
#'
#'The purpose of this function is to calculate the kinematic properties of a simulated galaxy, specifically the velocity and dispersion of particles within
#'certain user defined bins. The user must specify which direction they wish to study the kinematics using \code{bin_dir} (where \code{= "r"} specifies out
#'radially in 3D spherical bins, \code{= "cr"} specifies radially in 2D circular bins, and \code{= "z"} specifies directly in 1D out of the plane of the galaxy).
#'The default, \code{bin_dir = "r"} will return a comprehensive list of the simulation's kinematic properties (i.e. contained radial mass
#'and densities, velocities, velocity dispersions, anisotropy, rotational and circular velocities and spin parameter).
#'Other options for \code{bin_dir} will not contain the anisotropy, rotational and circular velocities or spin parameter as these are only physical when defined
#'using 3D spherical shells. If a simulation contains a vary large number of particles, it is possible to extract and analyse a sample of these to reduce
#'computational expense. This can be dome by specifying \code{samplerate = #a reduced number of particles}.
#'
#'@param filename The Gadget file containing the particle information of the galaxy to be analysed.
#'@param bin_type The direction in which to bin the simulation model - "r" (default) bins radially in 3D spherical shells, "cr" bins radially in 2D circular rings,
#' "z" bins in 1D off the plane of the galaxy.
#'@param rmax The maximum radial coordinate considered within the simulated galaxy in kpc.
#'@param rbin The number of radial bins considered.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the simulation, 1 - gas, 2 - dark matter, 3 - disc, 4 - bulge, 5 - stars,
#'6 - boundary.
#'@param samplerate \emph{Optional} If specified, the particle data will be randomly sampled with the number specified here.
#'@return A data frame containing kinematic features of radial bins including (at least) the outer radius of each bin (\code{$r}/\code{$cr}/\code{$z}), contained mass
#'(\code{$Mass}), logarithmic density (\code{$logp}), velocities and velocity dispersions (\code{$vr}/\code{$vcr}/\code{$vz} and \code{$sigma_vr}/\code{$sigma_vcr}/
#'\code{$sigma_vz}) and angular momentum magnitude (\code{$J}). In the case of \code{bin_dir = "r"} additionally the circular velocity (\code{$vc}), the velocity
#'anisotropy (\code{$B}), rotational velocity ($vrot) and the spin parameter (\code{$lambda}) are included in the output data frame.
#'@examples
#' \dontrun{
#'  sim_analysis(filename = "path/to/some/snapshot_XXX")
#'
#'  sim_analysis(filename    = "path/to/some/snapshot_XXX",
#'               pype        = c(3,4),
#'               rmax        = 300,
#'               rbin        = 100,
#'               samplerate  = 1000000)
#' }
#'

sim_analysis = function(filename, bin_type="r", rmax=200, rbin=200, ptype=NA, samplerate=NA){

  galaxy_data = snapshot::snapread(filename) # reading in the snapshot data into large list
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
    cat("Particles of ptype = c(", ptype[which(galaxy_data$head$Nall[ptype] == 0)], ") are missing in this model. \n")
    stop("Npart Error")
  } # error returned if a requested ptype is not present in the simulation

  galaxy_data$part = galaxy_data$part[galaxy_data$part$part_type %in% ptype,] # leaving only particles of requested ptype
  galaxy_data$head$Npart[p[!p %in% ptype]] = as.integer(0)
  galaxy_data$head$Nall[p[!p %in% ptype]] = as.integer(0) # removing record of the removed particles in the header

  if (is.na(samplerate) == FALSE){
    galaxy_data$part = galaxy_data$part[sample(nrow(galaxy_data$part), size = samplerate),]
  } # taking a sample of particles according to samplerate

  galaxy_df  = sim_galaxy(galaxy_data$part, centre=TRUE)
  grp        = matrix()
  grp_mass   = 0
  grp_J      = 0
  grp_vtheta = 0
  grp_vphi   = 0 # empty placeholders for loops
  G = 4.516e-29 # gravitational constant in units of kpc^3/([1e10 Msolar]s^2)

  if (bin_type == "r"){
    galaxy_cdf       = galaxy_df[galaxy_df$r < rmax,]
    galaxy_cdf$group = as.integer(cut(galaxy_cdf$r, breaks=seq(0,rmax,by=rmax/rbin), labels=seq(1,rbin)))
    grp_num          = data.frame("rbin" = seq(1, rmax, length.out=rbin), "Freq" = integer(rbin))
    grp_obins        = as.data.frame(table(with(galaxy_cdf, group)))
    grp_num[as.integer(levels(grp_obins$Var1)),2] = grp_obins$Freq
    grp_num$cumsum   = cumsum(grp_num$Freq)
    galaxy_odf       = galaxy_cdf[order(galaxy_cdf$group),]
    rbin_labels      = seq(0,rmax,rmax/rbin)
    profile = data.frame("r"        = numeric(rbin),
                         "Mass"     = numeric(rbin),
                         "logp"     = numeric(rbin),
                         "vc"       = numeric(rbin),
                         "J"        = numeric(rbin),
                         "vr"       = numeric(rbin),
                         "sigma_vr" = numeric(rbin),
                         "sigma_vt" = numeric(rbin),
                         "B"        = numeric(rbin),
                         "vrot"     = numeric(rbin),
                         "lambda"   = numeric(rbin))
    for (j in 1:rbin){
      if (j == 1){
        grp = galaxy_odf[1:grp_num$cumsum[j],]
      } else {
        grp = galaxy_odf[grp_num$cumsum[j-1]:grp_num$cumsum[j],]
      } # all data in an individual radius bin
      profile$r[j] = rbin_labels[j+1] # the outer edge of each radial bin
      grp_mass = sum(grp$Mass)
      if (j == 1){
        profile$Mass[j] = grp_mass
      } else {
        profile$Mass[j] = grp_mass + profile$Mass[j-1]
      } # finding the enclosed mass
      profile$logp[j] = log10(grp_mass / ((4 / 3) * pi * ((profile$r[j] * profile$r[j] * profile$r[j]) - (rbin_labels[j] * rbin_labels[j] * rbin_labels[j])))) # finding the log10 of shell density
      profile$vc[j] = (sqrt((G  * profile$Mass[j]) / profile$r[j])) * 3.086e16 # the circular velocity of particles at this radius, km/s
      grp_J = sum(grp$J)
      if (j == 1){
        profile$J[j] = grp_J
      } else {
        profile$J[j] = grp_J + profile$J[j-1]
      } # finding the enclosed angular momentum
      profile$vr[j] = mean(grp$vr) # finding the mean radial velocity in the shell, km/s
      grp_vtheta = mean(grp$vtheta)
      grp_vphi   = mean(grp$vphi)
      profile$sigma_vr[j] = sqrt(mean(grp$vr * grp$vr) - (profile$vr[j] * profile$vr[j])) # radial velocity dispersion, (km/s)
      profile$sigma_vt[j] = sqrt((mean(grp$vtheta * grp$vtheta) - (grp_vtheta * grp_vtheta)) +
                                   (mean(grp$vphi * grp$vphi) - (grp_vphi * grp_vphi))) # tangential velocity dispersion, (km/s)
      profile$B[j] = 1 - ((profile$sigma_vt[j] * profile$sigma_vt[j]) / (2 * profile$sigma_vr[j] * profile$sigma_vr[j])) # velocity anisotropy, unitless
      profile$vrot[j] = grp_J / (grp_mass * profile$r[j] * 3.086e16) # rotational velocity, km/s
      profile$lambda[j] = profile$J[j] / (1.414214 * profile$Mass[j] * profile$vc[j] * profile$r[j] * 3.086e16) # spin parameter, unitless
    }
  }
  if (bin_type == "cr"){
    galaxy_cdf       = galaxy_df[(galaxy_df$cr < rmax) & (abs(galaxy_df$z) < rmax),]
    galaxy_cdf$group = as.integer(cut(galaxy_cdf$cr, breaks=seq(0,rmax,rmax/rbin), labels=seq(1,rbin)))
    grp_num          = data.frame("rbin" = seq(1, rmax, length.out=rbin), "Freq" = integer(rbin))
    grp_obins        = as.data.frame(table(with(galaxy_cdf, group)))
    grp_num[as.integer(levels(grp_obins$Var1)),2] = grp_obins$Freq
    grp_num$cumsum   = cumsum(grp_num$Freq)
    galaxy_odf       = galaxy_cdf[order(galaxy_cdf$group),]
    rbin_labels      = seq(0,rmax,rmax/rbin)
    profile = data.frame("cr"        = numeric(rbin),
                         "Mass"     = numeric(rbin),
                         "logp"     = numeric(rbin),
                         "J"        = numeric(rbin),
                         "vcr"      = numeric(rbin),
                         "sigma_vcr"= numeric(rbin))
    for (j in 1:rbin){
      if (j == 1){
        grp = galaxy_odf[1:grp_num$cumsum[j],]
      } else {
        grp = galaxy_odf[grp_num$cumsum[j-1]:grp_num$cumsum[j],]
      } # all data in an individual radius bin
      profile$cr[j] = rbin_labels[j+1] # the outer edge of each radial bin
      grp_mass = sum(grp$Mass)
      if (j == 1){
        profile$Mass[j] = grp_mass
      } else {
        profile$Mass[j] = grp_mass + profile$Mass[j-1]
      } # finding the enclosed mass
      profile$logp[j] = log10(grp_mass / ((4 / 3) * pi * ((profile$cr[j] * profile$cr[j] * profile$cr[j]) - (rbin_labels[j] * rbin_labels[j] * rbin_labels[j])))) # finding the log10 of shell density
      if (j == 1){
        profile$J[j] = sum(grp$J)
      } else {
        profile$J[j] = sum(grp$J) + profile$J[j-1]
      } # finding the enclosed angular momentum
      profile$vcr[j] = mean(grp$vcr) # finding the mean radial velocity in the shell, km/s
      profile$sigma_vcr[j] = sqrt(mean(grp$vcr * grp$vcr) - (profile$vcr[j] * profile$vcr[j])) # radial velocity dispersion, (km/s)
    }
  }
  if (bin_type == "z"){
    galaxy_cdf       = galaxy_df[(galaxy_df$cr < rmax) & ((galaxy_df$z) < rmax) & ((galaxy_df$z) > 0),]
    galaxy_cdf$group = as.integer(cut(galaxy_cdf$z, breaks=seq(0,rmax,rmax/rbin), labels=seq(1,rbin)))
    grp_num          = data.frame("rbin" = seq(1, rmax, length.out=rbin), "Freq" = integer(rbin))
    grp_obins        = as.data.frame(table(with(galaxy_cdf, group)))
    grp_num[as.integer(levels(grp_obins$Var1)),2] = grp_obins$Freq
    grp_num$cumsum   = cumsum(grp_num$Freq)
    galaxy_odf       = galaxy_cdf[order(galaxy_cdf$group),]
    rbin_labels      = seq(0,rmax,rmax/rbin)
    profile = data.frame("z"        = numeric(rbin),
                         "Mass"     = numeric(rbin),
                         "logp"     = numeric(rbin),
                         "J"        = numeric(rbin),
                         "vz"      = numeric(rbin),
                         "sigma_vz"= numeric(rbin))
    for (j in 1:rbin){
      if (j == 1){
        grp = galaxy_odf[1:grp_num$cumsum[j],]
      } else {
        grp = galaxy_odf[grp_num$cumsum[j-1]:grp_num$cumsum[j],]
      } # all data in an individual radius bin
      profile$z[j] = rbin_labels[j+1] # the outer edge of each radial bin
      grp_mass = sum(grp$Mass)
      if (j == 1){
        profile$Mass[j] = grp_mass
      } else {
        profile$Mass[j] = grp_mass + profile$Mass[j-1]
      } # finding the enclosed mass
      profile$logp[j] = log10(grp_mass / ((4 / 3) * pi * ((profile$z[j] * profile$z[j] * profile$z[j]) - (rbin_labels[j] * rbin_labels[j] * rbin_labels[j])))) # finding the log10 of shell density
      if (j == 1){
        profile$J[j] = sum(grp$J)
      } else {
        profile$J[j] = sum(grp$J) + profile$J[j-1]
      } # finding the enclosed angular momentum
      profile$vz[j] = mean(grp$vz) # finding the mean radial velocity in the shell, km/s
      profile$sigma_vz[j] = sqrt(mean(grp$vz * grp$vz) - (profile$vz[j] * profile$vz[j])) # radial velocity dispersion, (km/s)
    }
  } # trimming the galaxy to specified rmax and grouping data into bins based on the type of binning

  return(profile)

}
