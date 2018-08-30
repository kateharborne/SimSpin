# Kate Harborne (last edit - 23/04/18)
#'Find the half mass radius of a simulated galaxy.
#'
#'The purpose of this function is to find the half mass radius of a simulated galaxy. It works by
#' growing a sphere and returning the radial coordinate at which half the total number of particles
#' in the galaxy is contained.
#'
#'@param filename The Gadget file containing the particle information of the galaxy to be analysed.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the
#' simulation, 1 - gas, 2 - dark matter, 3 - disc, 4 - bulge, 5 - stars, 6 - boundary.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param fract The fraction of particles to be contained within the radius calculated. Default is
#' 0.5, i.e. R50.
#'@return Returned is a data.frame that contains the half mass radius, the number of particles
#' contained within that radius, the total number of particles in the simulation and the exact
#' fraction of particles contained within that radius.
#'@examples
#' data      = obs_data_prep(filename = system.file("extdata", 'S0_vignette', package="SimSpin"))
#' ifucube   = ifu_cube(obs_data = data)
#'
#' output = find_r50(filename     = system.file("extdata", 'S0_vignette', package="SimSpin"),
#'                   ptype        = NA,
#'                   r200         = 20)
#'

find_r50 = function(filename, ptype=NA, r200=200, fract=0.5){

  galaxy_data = snapshot::snapread(filename)               # reading the data into large list
  galaxy_data$part$part_type = rep(0, nrow(galaxy_data$part))
  p = seq(1,6)                                             # all possible particle values
  ppart = cumsum(galaxy_data$head$Npart[which(galaxy_data$head$Nall[p] != 0)])
  # present particles

  for (i in 1:length(ppart)){
    if (i == 1){
      galaxy_data$part[1:as.integer(ppart[i]),]$part_type =  which(galaxy_data$head$Nall[p] != 0)[i]
    } else {
      galaxy_data$part[as.integer(ppart[i-1]+1):as.integer(ppart[i]),]$part_type = which(galaxy_data$head$Nall[p] != 0)[i]
    }
  }
  # label the data with particle types

  if (is.na(ptype[1])){ptype = which(galaxy_data$head$Nall[p] != 0)}
  # for all particles, leave ptype = NA
  if (0 %in% galaxy_data$head$Nall[ptype]){
    cat("Particles of ptype =", which(galaxy_data$head$Nall[ptype] %in% 0), " are missing in this model. \n")
    stop("Npart Error")
  }
  # error returned if ptype is not present

  galaxy_data$part = galaxy_data$part[galaxy_data$part$part_type %in% ptype,]
  # leaving only requested ptype
  galaxy_data$head$Npart[p[!p %in% ptype]] = as.integer(0)
  galaxy_data$head$Nall[p[!p %in% ptype]] = as.integer(0)  # removing record in the header

  galaxy_df  = sim_galaxy(galaxy_data$part, centre=TRUE)   # extracting position & LOS velocities
  ntotal     = as.numeric(nrow(galaxy_df))                 # total number of particles in the galaxy
  galaxy_cdf = galaxy_df[galaxy_df$r<r200,]                # removing particles beyond r200
  inside = galaxy_cdf[galaxy_cdf$r < 1,]                   # particles within the sphere
  etotal = as.numeric(dim(inside)[1])                      # number of particles within the ellipse
  fac = 1

  if (etotal < ntotal*fract){
    while (etotal < ntotal*fract){
      fac  = fac + 1
      inside = galaxy_cdf[galaxy_cdf$r < fac,]
      etotal = as.numeric(dim(inside)[1])
    }
    # if there are fewer than ntotal/2 particles inside, the ellipse grows by small amount until
    #  it contains more than ntotal/2

    diff = etotal - (ntotal*fract) # checks the difference between the total contained and ntotal/2

    if (diff > (0.0005 * ntotal*fract)){
      fac = fac - 1
      inside = galaxy_cdf[galaxy_cdf$r < fac,]
      etotal = as.numeric(dim(inside)[1])
      while (etotal < ntotal*fract){
        fac = fac + 0.01
        inside = galaxy_cdf[galaxy_cdf$r < fac,]
        etotal = as.numeric(dim(inside)[1])
      }
      # if the difference is greater than 0.5% of ntotal/2, the ellipse shrinks again by one
      #  factor and increases in smaller increments until etotal > ntotal/2 again
    }

    r50 = data.frame("r50" = fac,
                     "p_inside" = etotal,
                     "ntotal"   = ntotal,
                     "fract"    = etotal/ntotal)

    return(r50)
  }

  if (etotal > ntotal*fract){
    while (etotal > ntotal*fract){
      fac  = fac - 0.01
      inside = galaxy_cdf[galaxy_cdf$r < fac,]
      etotal = as.numeric(dim(inside)[1])
    }
    # if there are more than ntotal/2 particles inside, the ellipse shrinks by small amount until
    #  it contains less than ntotal/2

    diff = (ntotal*fract) - etotal # checks the difference between the total contained and ntotal/2

    if (diff > (0.0005 * ntotal*fract)){
      fac = fac + 0.01
      inside = galaxy_cdf[galaxy_cdf$r < fac,]
      etotal = as.numeric(dim(inside)[1])
      while (etotal > ntotal*fract){
        fac = fac - 0.0005
        inside = galaxy_cdf[galaxy_cdf$r < fac,]
        etotal = as.numeric(dim(inside)[1])
      }
      # if the difference is greater than 0.5% of ntotal/2, the ellipse grows again by one factor
      #  and decreases in smaller increments until etotal < ntotal/2 again
    }

    r50 = data.frame("r50" = fac,
                     "p_inside" = etotal,
                     "ntotal"   = ntotal,
                     "fract"    = etotal/ntotal)

    return(r50)
  }
}
