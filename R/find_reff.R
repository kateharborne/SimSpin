# Kate Harborne (last edit - 23/04/18)
#'Find the effective radius of a simulated galaxy.
#'
#'The purpose of this function is to find the observed effective radius of a simulated galaxy
#' projected at a user supplied inclination. It works by ordering particles by their observed
#' radius and returning the radial coordinate at which half the total number of particles in the
#' galaxy is contained.
#'
#'@param filename The Gadget file containing the particle information of the galaxy to be analysed.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the
#' simulation, 1 - gas, 2 - dark matter, 3 - disc, 4 - bulge, 5 - stars, 6 - boundary.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param inc_deg The observed inclination of the simulated galaxy in degrees.
#'@param fract The fraction of particles to be contained within the radius calculated. Default is
#' 0.5, i.e. Reff.
#'@param axis_ratio A data frame containing the semi-major and semi-minor axes lengths for the
#' observed galaxy, as given by \code{\link{ifu_cube}} or \code{\link{blur_cube}}.
#'@param angular_size The kpc/'' scaling factor output by \code{\link{obs_data_prep}} used to give
#' the axis ratios in both kpc and ''.
#'@return Returned is a scaled axis ratio that describes the semi-major and semi-minor axes of an
#' ellipse that contains half the total number of particles.
#'@examples
#' data      = obs_data_prep(filename = system.file("extdata", 'S0_vignette', package="SimSpin"))
#' ifucube   = ifu_cube(obs_data = data)
#'
#' output = find_reff(filename     = system.file("extdata", 'S0_vignette', package="SimSpin"),
#'                    ptype        = NA,
#'                    r200         = 10,
#'                    inc_deg      = 0,
#'                    axis_ratio   = ifucube$axis_ratio,
#'                    angular_size = data$angular_size)
#'

find_reff = function(filename, ptype=NA, r200=200, inc_deg, fract=0.5, axis_ratio, angular_size){

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

  inc_rad    = inc_deg * (pi / 180)                        # the galaxy inclination in radians
  galaxy_df  = obs_galaxy(galaxy_data$part, centre=TRUE, inc_rad)
                                                           # extracting position & LOS velocities
  ntotal     = as.numeric(nrow(galaxy_df))                 # total number of particles in the galaxy
  galaxy_cdf = galaxy_df[galaxy_df$r_obs<r200,]            # removing particles beyond r200
  elli = galaxy_cdf$ID[((galaxy_cdf$x^2  / axis_ratio$a^2) + (galaxy_cdf$z_obs^2 / axis_ratio$b^2)) <= 1]
                                                           # particles within the ellipse
  etotal = as.integer(length(elli))                        # number of particles within the ellipse
  fac = 1

  if (etotal < ntotal*fract){
    while (etotal < ntotal*fract){
      fac  = fac + 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf$ID[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
      etotal = length(elli)
      }
    # if there are fewer than ntotal/2 particles inside, the ellipse grows by small amount until
    #  it contains more than ntotal/2

    diff = etotal - (ntotal*fract) # checks the difference between the total contained and ntotal/2

    if (diff > (0.005 * ntotal*fract)){
      fac = fac - 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf$ID[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
      etotal = length(elli)
      while (etotal < ntotal*fract){
        fac = fac + 0.0005
        axes = axis_ratio * fac
        elli = galaxy_cdf$ID[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
        etotal = length(elli)
        }
      # if the difference is greater than 0.5% of ntotal/2, the ellipse shrinks again by one
      #  factor and increases in smaller increments until etotal > ntotal/2 again
      }

    ellipse_axis_ratio = data.frame("a_kpc"     = axis_ratio$a * fac,
                                    "b_kpc"     = axis_ratio$b * fac,
                                    "a_arcsec"  = (axis_ratio$a * fac) / angular_size,
                                    "b_arcsec"  = (axis_ratio$b * fac) / angular_size,
                                    "angle"     = 90)

    return(ellipse_axis_ratio)
   }

  else {

    while (etotal > ntotal*fract){
      fac  = fac - 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1,]
      etotal = length(elli)
    }
    # if there are more than ntotal/2 particles inside, the ellipse shrinks by small amount until
    #  it contains less than ntotal/2

    diff = (ntotal*fract) - etotal # checks the difference between the total contained and ntotal/2

    if (diff > (0.005 * ntotal*fract)){
      fac = fac + 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf$ID[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
      etotal = length(elli)
      while (etotal > ntotal*fract){
        fac = fac - 0.0005
        axes = axis_ratio * fac
        elli = galaxy_cdf$ID[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
        etotal = length(elli)
      }
      # if the difference is greater than 0.5% of ntotal/2, the ellipse grows again by one factor
      #  and decreases in smaller increments until etotal < ntotal/2 again
    }

    ellipse_axis_ratio = data.frame("a_kpc"     = axis_ratio$a * fac,
                                    "b_kpc"     = axis_ratio$b * fac,
                                    "a_arcsec"  = (axis_ratio$a * fac) / angular_size,
                                    "b_arcsec"  = (axis_ratio$b * fac) / angular_size,
                                    "angle"     = 90)
    return(ellipse_axis_ratio)
  }
}
