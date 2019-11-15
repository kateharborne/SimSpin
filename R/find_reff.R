# Kate Harborne (last edit - 15/11/19)
#'Find the effective radius of a simulated galaxy.
#'
#'The purpose of this function is to find the observed effective radius of a simulated galaxy
#' projected at a user supplied inclination. It works by ordering particles by their observed
#' radius and returning the radial coordinate at which half the total number of particles in the
#' galaxy is contained.
#'
#'@param simdata The simulation information data.frame output by \code{\link{sim_data}}.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param inc_deg The observed inclination of the simulated galaxy in degrees.
#'@param fract The fraction of particles to be contained within the radius calculated. Default is
#' 0.5, i.e. Reff.
#'@param axis_ratio A data frame containing the semi-major and semi-minor axes lengths for the
#' observed galaxy, as given by \code{\link{obs_imgs}}.
#'@return Returned is a data.frame containing the scaled axis ratio that describes the semi-major
#' and semi-minor axes of an ellipse that contains the specified fraction of the total number of
#' particles.
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' fluxes      = flux_grid(obs_data = data)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#' images      = obs_images(obs_data = data, ifu_datacube = cube)
#'
#' output = find_reff(simdata      = galaxy_data,
#'                    r200         = 10,
#'                    inc_deg      = 0,
#'                    axis_ratio   = images$axis_ratio)
#'

find_reff = function(simdata, r200=200, inc_deg, fract=0.5, axis_ratio){

  result = grepl(paste(c("PartType2", "PartType3", "PartType4"), collapse = "|"), names(simdata))
  # finding the luminous matter within the simulation for imaging

  if (any(result)) {                                       # concatenating the seperate particle
    present = which(result)                                #  data.frames into a single data.frame
    if (length(present) == 1){
      galaxy_data = simdata[[present[1]]]$Part
    } else if (length(present) == 2){
      galaxy_data = rbind(simdata[[present[1]]]$Part, simdata[[present[2]]]$Part)
    } else if (length(present) == 3){
      galaxy_data = rbind(simdata[[present[1]]]$Part, simdata[[present[2]]]$Part, simdata[[present[3]]]$Part)
    }
  } else {
    cat("There are no particles representing luminous matter in this simulation (i.e. no stars, bulge or disc particles). \n")
    stop("LumPart Error") # if no stars, bulge or disc component, stop trying to build IFU cube
  }

  inc_rad    = inc_deg * (pi / 180)                        # the galaxy inclination in radians
  galaxy_data = .reorient_galaxy(galaxy_data)              # reorient galaxy to horizontal
  galaxy_df  = obs_galaxy(galaxy_data, centre=TRUE, inc_rad)
                                                           # extracting position & LOS velocities

  galaxy_cdf = galaxy_df[galaxy_df$r_obs<r200,]            # removing particles beyond r200
  ntotal     = as.numeric(nrow(galaxy_cdf))                # total number of particles in the galaxy within r200
  elli = galaxy_cdf$x[((galaxy_cdf$x^2  / axis_ratio$a^2) + (galaxy_cdf$z_obs^2 / axis_ratio$b^2)) <= 1]
                                                           # particles within the ellipse
  etotal = as.integer(length(elli))                        # number of particles within the ellipse
  fac = 1

  if (etotal < ntotal*fract){
    while (etotal < ntotal*fract){
      fac  = fac + 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf$x[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
      etotal = length(elli)
      }
    # if there are fewer than ntotal/2 particles inside, the ellipse grows by small amount until
    #  it contains more than ntotal/2

    diff = etotal - (ntotal*fract) # checks the difference between the total contained and ntotal/2

    if (diff > (0.005 * ntotal*fract)){
      fac = fac - 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf$x[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
      etotal = length(elli)
      while (etotal < ntotal*fract){
        fac = fac + 0.0005
        axes = axis_ratio * fac
        elli = galaxy_cdf$x[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
        etotal = length(elli)
        }
      # if the difference is greater than 0.5% of ntotal/2, the ellipse shrinks again by one
      #  factor and increases in smaller increments until etotal > ntotal/2 again
      }

    ellipse_axis_ratio = data.frame("a"     = axis_ratio$a * fac,
                                    "b"     = axis_ratio$b * fac)

    return(ellipse_axis_ratio)
   }

  else {

    while (etotal > ntotal*fract){
      fac  = fac - 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf$x[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
      etotal = length(elli)
    }
    # if there are more than ntotal/2 particles inside, the ellipse shrinks by small amount until
    #  it contains less than ntotal/2

    diff = (ntotal*fract) - etotal # checks the difference between the total contained and ntotal/2

    if (diff > (0.005 * ntotal*fract)){
      fac = fac + 0.01
      axes = axis_ratio * fac
      elli = galaxy_cdf$x[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
      etotal = length(elli)
      while (etotal > ntotal*fract){
        fac = fac - 0.0005
        axes = axis_ratio * fac
        elli = galaxy_cdf$x[((galaxy_cdf$x^2  / axes$a^2) + (galaxy_cdf$z_obs^2 / axes$b^2)) <= 1]
        etotal = length(elli)
      }
      # if the difference is greater than 0.5% of ntotal/2, the ellipse grows again by one factor
      #  and decreases in smaller increments until etotal < ntotal/2 again
    }

    ellipse_axis_ratio = data.frame("a"     = axis_ratio$a * fac,
                                    "b"     = axis_ratio$b * fac)

    return(ellipse_axis_ratio)
  }
}
