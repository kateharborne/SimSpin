 # Kate Harborne (last edit - 13/09/2017)
#'Find the 3D effective radius of a simulated galaxy.
#'
#'The purpose of this function is to find the observed effective radius of a simulated galaxy projected at a user supplied inclination. It works
#'by ordering particles by their observed radius and returning the radial coordinate at which half the total number of particles in the galaxy is contained.
#'
#'@param filename The Gadget file containing the particle information of the galaxy to be analysed.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the simulation, 1 - gas, 2 - dark matter, 3 - disc, 4 - bulge, 5 - stars, 6 - boundary.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param fract The fraction of particles to be contained within the radius calculated. Default is 0.5, i.e. Reff.
#'@return Returned is a scaled axis ratio that describes the semimajor and semiminor axes of an ellipse that contains half the total number of particles.
#'@examples
#' \dontrun{
#' find_reff(filename = "path/to/some/snapshot_XXX",
#'           ptype    = 0,
#'           r200     = 10,
#'           inc_deg  = 0)
#' }
#'

find_reff_3d = function(filename, ptype=NA, r200=200, fract=0.5){

  galaxy_data = snapshot::snapread(filename) # reading in the snapshot data into large list
  galaxy_data$part$part_type = rep(0, nrow(galaxy_data$part))
  p = seq(1,6) # all possible particle values
  ppart = cumsum(galaxy_data$head$Npart[which(galaxy_data$head$Nall[p] != 0)]) # present particles

  for (i in 1:length(ppart)){
    if (i == 1){
      galaxy_data$part[1:as.integer(ppart[i]),]$part_type =  which(galaxy_data$head$Nall[p] != 0)[i]
    } else {
      galaxy_data$part[as.integer(ppart[i-1]+1):as.integer(ppart[i]),]$part_type = which(galaxy_data$head$Nall[p] != 0)[i]
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

  galaxy_df  = sim_galaxy(galaxy_data$part, centre=TRUE) # extracting the position and line of sight of galaxy particles for a given inclination
  ntotal     = as.numeric(nrow(galaxy_df)) # total number of particles in the galaxy
  galaxy_cdf = galaxy_df[galaxy_df$r<r200,] # removing particles further than r200 kpc from the centre
  fac = 1
  elli = galaxy_cdf$r[galaxy_cdf$r <= fac] # particles within a sphere of radius 1 kpc
  etotal = as.integer(length(elli)) # number of particles within the ellipse

  if (etotal < ntotal*fract){
    while (etotal < ntotal*fract){
      fac  = fac + 1
      elli = galaxy_cdf$r[galaxy_cdf$r <= fac]
      etotal = length(elli)
    } # if there are fewer than ntotal/2 particles inside, the sphere grows by small amount until it contains more than ntotal/2

    diff = etotal - (ntotal*fract) # checks the difference between the total contained and ntotal/2

    if (diff > (0.005 * ntotal*fract)){
      fac = fac - 1
      elli = galaxy_cdf$r[galaxy_cdf$r <= fac]
      etotal = length(elli)
      while (etotal < ntotal*fract){
        fac = fac + 0.05
        elli = galaxy_cdf$r[galaxy_cdf$r <= fac]
        etotal = length(elli)
      } # if the difference is greater than 0.5% of ntotal/2, the sphere shrinks again by one factor and increases in smaller increments until etotal > ntotal/2 again
    }

    ellipse_axis_ratio = data.frame("reff" = fac)

    return(ellipse_axis_ratio)
  }

  if (etotal > ntotal*fract){
    while (etotal > ntotal*fract){
      fac  = fac - 1
      elli = galaxy_cdf$r[galaxy_cdf$r <= fac]
      etotal = length(elli)
    } # if there are more than ntotal/2 particles inside, the sphere shrinks by small amount until it contains less than ntotal/2

    diff = (ntotal*fract) - etotal # checks the difference between the total contained and ntotal/2

    if (diff > (0.005 * ntotal*fract)){
      fac = fac + 1
      elli = galaxy_cdf$r[galaxy_cdf$r <= fac]
      etotal = length(elli)
      while (etotal > ntotal*fract){
        fac = fac - 0.05
        elli = galaxy_cdf$r[galaxy_cdf$r <= fac]
        etotal = length(elli)
      } # if the difference is greater than 0.5% of ntotal/2, the sphere grows again by one factor and decreases in smaller increments until etotal < ntotal/2 again
    }

    ellipse_axis_ratio = data.frame("reff" = fac)

    return(ellipse_axis_ratio)
  }
}
