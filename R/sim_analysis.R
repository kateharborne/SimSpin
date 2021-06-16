# Author: Kate Harborne
# Date: 21/04/2021
# Title: sim_analysis - a function for analysing the particle data from the
# simulation directly
#
#'A function for measuring the particle properties of the input simulation
#'
#'The purpose of this function is to measure the true particle properties of the
#' simulation prior to observation. This function produces a list that includes
#' a summary of the total galaxy, a summary of the half-mass properties, and a
#' data.frame containing the radial trends of a series of properties, computed
#' in a series of spherical shells.
#'
#'@param simspin_file The path to the location of the SimSpin .Rdata file OR
#' output list from \code{make_simspin_file()}.
#'@param type String "stars" (default) or "gas" to specify which set of
#' of particles are used in the property calculations.
#'@param bin_breaks Optional parameter that allows you to specify radial bin
#'  break positions. Default will give bins spaced with varying sized bins:
#'  `seq(0, 9, by=1), seq(12, 51, by=3), seq(61, 101, by=10), seq(151, 501, by=50)`
#'@param half_mass If simulation file contains all particles cutout from a box
#' (rather than just particles from a single galaxy), you can the half-mass
#' value at which the alignment function is run. Numeric length = 1. Default is
#' NA, in which case half the total mass of the suplied simulation data is used.
#'@return Returns a list that contains:
#'\describe{
#'   \item{Properties}{list()}
#'   \item{HalfMassProperties}{list()}
#'   \item{RadialTrends}{data.frame()}}
#'   where \code{Properties} include:
#'   \describe{\item{Type}{Component considered within analysis}
#'             \item{TotalMass}{Total mass (solar)}
#'             \item{MeanAge}{Mean age (Gyr)}
#'             \item{MeanMetallicity}{Mean metallicity (fraction of solar)}
#'             \item{NumberOfParticles}{Total number of particles}}
#'   where \code{HalfMassProperties} include:
#'   \describe{\item{Mass}{Half mass (solar)}
#'             \item{RadiusCircular}{Circularised radius at half-mass (kpc)}
#'             \item{RadiusElliptical}{Elliptical radius at half-mass given shapes p & q (kpc)}
#'             \item{Shape_p}{Axis ratio (b/a) of particles within half-mass}
#'             \item{Shape_q}{Axis ratio (c/a) of particles within half-mass}}
#'   and where \code{RadialTrends} include:
#'   \describe{\item{Radius}{Radial coordinate at the centre of the radial bin (kpc)}
#'             \item{Mass}{Mass (solar) contained within radial bin}
#'             \item{CumulativeMass}{Mass (solar) of all particles contained within radius of this bin}
#'             \item{Density}{Mass density of shell (Msol/kpc^3)}
#'             \item{Age}{Mean age (Gyr) of particles within radial bin}
#'             \item{Metallicity}{Mean metallicity (fraction of solar) of particles within radial bin}
#'             \item{CircularVelocity}{Circular velocity (km/s) due to mass contained within radius of this bin}
#'             \item{RotationalVelocity}{Mean velocity along circular orbits (km/s) in radial bin}
#'             \item{RotationalDispersion}{Standard deviation of velocities in circular orbits (km/s) in radial bin}
#'             \item{VelocityAnisotropy}{Beta parameter of particles within radial bin}
#'             \item{SpinParameter_Bullock}{Bullock et al (2001) spin parameter measured from all particles contained within radius of this bin}
#'             \item{SpecificAngularMomentum}{Mass weighted angular momentum, j (kpc km/s), of particles within radial bin}
#'             \item{Shape_p}{Axis ratio (b/a) of all particles within radius of this bin}
#'             \item{Shape_q}{Axis ratio (c/a) of all particles within radius of this bin},
#'             \item{NumberOfParticles}{Number of particles contained within radial bin}}
#'
#'@examples
#'ss_gadget = system.file("extdata", "SimSpin_example_Gadget_spectra.Rdata",
#'                         package = "SimSpin")
#'props = sim_analysis(simspin_file = ss_gadget)
#'

sim_analysis = function(simspin_file, type = "stars", half_mass = NA, bin_breaks=NA){

  # Reading in SimSpin file data
  if (typeof(simspin_file) == "character"){ # if provided with path to file
    simspin_data = readRDS(simspin_file)
  }
  if (typeof(simspin_file) == "list"){ # if provided with output list
    simspin_data = simspin_file
  }

  # Determining which table to examine
  if (stringr::str_to_lower(type) == "stars"){
    galaxy_data = simspin_data$star_part
  }
  if (stringr::str_to_lower(type) == "gas"){
    galaxy_data = simspin_data$gas_part
    if (is.null(galaxy_data)){
      stop(cat("Error: There are no gas particles contained within this simulation. Try again with type = 'stars'"))
    } else {
      galaxy_data$Age = rep(NA, length(galaxy_data$ID))
    }
  }

  if (is.na(half_mass)){
    half_mass = sum(galaxy_data$Mass)/2
  } else {
    if (half_mass > sum(galaxy_data$Mass)){
      stop(cat("Error: The total mass in the simulation is less than the requested half-mass radius. Something is wrong."))
    }
  }

  if (is.na(bin_breaks[1])){
    lseq = c(seq(0, 9, by=1), seq(12, 51, by=3), seq(61, 101, by=10), seq(151, 501, by=50))
    rbins = length(lseq)
    bin_ends = c(seq(1, 9, by=1), seq(12, 51, by=3), seq(61, 101, by=10), seq(151, 551, by=50))
  } else {
    if (bin_breaks[1]!=0){
      lseq = c(0, bin_breaks)
    } else {
      lseq = bin_breaks
    }
    rbins = length(lseq)
    bin_ends = lseq + c(diff(lseq), diff(lseq)[(rbins-1)])
  }

  analysis_data = list("Properties" = list("Type" = type,
                                           "TotalMass" = sum(galaxy_data$Mass),
                                           "MeanAge" = mean(galaxy_data$Age),
                                           "MeanMetallicity" = mean(galaxy_data$Metallicity),
                                           "NumberOfParticles" = length(galaxy_data$ID)),
                       "HalfMassProperties" = list("Mass" = half_mass,
                                                   "RadiusCircular" = numeric(1),
                                                   "RadiusElliptical" = numeric(1),
                                                   "Shape_p" = numeric(1),
                                                   "Shape_q" = numeric(1)),
                       "RadialTrends" = data.frame("Radius"  = lseq + ((bin_ends - lseq)/2),
                                                   "Mass"    = numeric(rbins),
                                                   "CumulativeMass" = numeric(rbins),
                                                   "Density" = numeric(rbins),
                                                   "Age"     = numeric(rbins),
                                                   "Metallicity" = numeric(rbins),
                                                   "CircularVelocity" = numeric(rbins),
                                                   "RotationalVelocity" = numeric(rbins),
                                                   "RotationalDispersion" = numeric(rbins),
                                                   "VelocityAnisotropy" = numeric(rbins),
                                                   "SpinParameter_Bullock" = numeric(rbins),
                                                   "SpecificAngularMomentum" = numeric(rbins),
                                                   "Shape_p" = numeric(rbins),
                                                   "Shape_q" = numeric(rbins),
                                                   "NumberOfParticles" = numeric(rbins)))

  galaxy_data$r = sqrt((galaxy_data$x^2) + (galaxy_data$y^2) + (galaxy_data$z^2))
  galaxy_data$theta = atan2(y = galaxy_data$y, x = galaxy_data$x)
  galaxy_data$phi = atan2(y = sqrt((galaxy_data$x^2) + (galaxy_data$y^2)), x = galaxy_data$z)
  galaxy_data$v_r = sqrt((galaxy_data$vx^2) + (galaxy_data$vy^2) + (galaxy_data$vz^2))
  galaxy_data$v_theta = atan2(y = galaxy_data$vy, x = galaxy_data$vx)
  galaxy_data$v_phi = atan2(y = sqrt((galaxy_data$vx^2) + (galaxy_data$vy^2)), x = galaxy_data$vz)

  galaxy_data$rbin = cut(galaxy_data$r, breaks=lseq, labels=F)
  particle_ID = NULL # initiallising varible to avoid CRAN error in checks

  galaxy_data_table = data.table::data.table("particle_ID" = seq(1, length(galaxy_data$x)), "rbin_ID"=galaxy_data$rbin)

  # which particles sit in each rbin?
  part_in_rbin = galaxy_data_table[, list(val=list(particle_ID)), by = "rbin_ID"]

  for (bin in 1:rbins){
    binID = part_in_rbin$rbin_ID[bin]
    if (!is.na(binID)){
      sample = galaxy_data[part_in_rbin$val[[bin]],]

      analysis_data$RadialTrends$Mass[binID] = sum(sample$Mass)
      analysis_data$RadialTrends$Age[binID] = mean(sample$Age)
      analysis_data$RadialTrends$Metallicity[binID] = mean(sample$Metallicity)
      analysis_data$RadialTrends$RotationalVelocity[binID] = mean(sample$v_theta)
      analysis_data$RadialTrends$RotationalDispersion[binID] = sd(sample$v_theta)
      analysis_data$RadialTrends$VelocityAnisotropy[binID] = 1 - ((sd(sample$v_theta)^2 + sd(sample$v_phi)^2) / (2 * sd(sample$v_r)^2))
      analysis_data$RadialTrends$SpecificAngularMomentum[binID] = sqrt(sum(angmom_galaxy(sample[,1:8])^2))/sum(sample$Mass)
      analysis_data$RadialTrends$NumberOfParticles[binID] = length(sample$ID)
    } else {
      analysis_data$RadialTrends$Mass[binID] = NA
      analysis_data$RadialTrends$Age[binID] = NA
      analysis_data$RadialTrends$Metallicity[binID] = NA
      analysis_data$RadialTrends$RotationalVelocity[binID] = NA
      analysis_data$RadialTrends$RotationalDispersion[binID] = NA
      analysis_data$RadialTrends$VelocityAnisotropy[binID] = NA
      analysis_data$RadialTrends$SpecificAngularMomentum[binID] = NA
      analysis_data$RadialTrends$NumberOfParticles[binID] = 0
    }

  }

  bin_vol = (4/3)*pi*(bin_ends^3)
  analysis_data$RadialTrends$CumulativeMass = cumsum(analysis_data$RadialTrends$Mass)
  analysis_data$RadialTrends$Density = analysis_data$RadialTrends$CumulativeMass / bin_vol
  analysis_data$RadialTrends$CircularVelocity = sqrt(.g_in_kpcMsolkms2 * analysis_data$RadialTrends$CumulativeMass / bin_ends)

  for (cbin in 1:rbins){
    sample = galaxy_data[galaxy_data$r <= bin_ends[cbin],]
    vc = analysis_data$RadialTrends$CircularVelocity[cbin]
    M  = analysis_data$RadialTrends$CumulativeMass[cbin]
    analysis_data$RadialTrends$SpinParameter_Bullock[cbin] = sqrt(sum(angmom_galaxy(sample[,1:8])^2)) / (sqrt(2)*M*vc*bin_ends[cbin])
    shapes = .measure_pqj(galaxy_data = list("star_part" = galaxy_data), half_mass = M)
    analysis_data$RadialTrends$Shape_p[cbin] = shapes$p
    analysis_data$RadialTrends$Shape_q[cbin] = shapes$q
  }

  shapes = .measure_pqj(galaxy_data = list("star_part" = galaxy_data), half_mass = analysis_data$HalfMassProperties$Mass)
  analysis_data$HalfMassProperties$Shape_p = shapes$p
  analysis_data$HalfMassProperties$Shape_q = shapes$q

  hm_props = .new_half_mass_data(galaxy_data = galaxy_data, p = shapes$p, q = shapes$q, half_mass = analysis_data$HalfMassProperties$Mass)
  analysis_data$HalfMassProperties$RadiusCircular = max(hm_props$r)
  analysis_data$HalfMassProperties$RadiusElliptical = max(sqrt((hm_props$x*hm_props$x) + ((hm_props$y/shapes$p)*(hm_props$y/shapes$p)) + ((hm_props$z/shapes$q)*(hm_props$z/shapes$q))))

  return(analysis_data)
}
