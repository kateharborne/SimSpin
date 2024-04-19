# Author: Kate Harborne
# Co-author: Alice Serene
# Date: 10/01/2023
# Title: Testing the make_simspin_file.R code
#'Reformatting isolated galaxy simulations to contain spectra.
#'
#'The purpose of this function is to construct a SimSpin file containing the
#' mock spectra for each particle contained within the galaxy simulation file.
#' If the snapshot provided is from a cosmological simulation, the SEDs
#' generated will be with respect to the Stellar Formation Time/Age, Metallicity
#' and Initial Mass of each stellar particle. If the system is an N-body model,
#' stellar particles are assumed to have an age and metallicity as provided to
#' the function as \code{disk_age}, \code{bulge_age}, \code{disk_Z} and
#' \code{bulge_Z}. Returned is an .Rdata file in a SimSpin readable format.
#'
#'@param filename The path to the snapshot file.
#'@param cores The number of cores across which to multi-thread the problem.
#'@param disk_age The age of the disk particles in Gyr.
#'@param bulge_age The age of the bulge particles in Gyr.
#'@param disk_Z The metallicity of the disk particles as a mass fraction (mass
#' of all metal elements above He over the total mass).
#'@param bulge_Z The metallicity of the bulge particles as a mass fraction (mass
#' of all metal elements above He over the total mass).
#'@param template The stellar templates from which to derive the SEDs. Options
#' include "BC03lr" (GALEXEV low resolution, Bruzual & Charlot 2003), "BC03hr"
#' (GALEXEV high resolution, Bruzual & Charlot 2003) or "EMILES" (Vazdekis et
#' al, 2016).
#'@param write_to_file Boolean to specify whether the list produced should be
#' written to a ".Rdata" file or output to the environment. Default is TRUE, so
#' that files can be re-observed without having to generate spectra each time.
#'@param output The path at which the output file is written. If not provided,
#' file will be written at the location of the input filename with the addition
#' of "_spectra.Rdata".
#'@param overwrite If true, and the file already exists at the output location,
#' a new file will be written over the old one.
#'@param centre If simulation file contains all particles cutout from a box
#' (rather than just particles from a single galaxy), you can specify the point
#' around which the view should be centred. Numeric length = 3. Default is NA,
#' in which case the system is centred around the median position. Specified in
#' units of kpc.
#'@param half_mass If simulation file contains all particles cutout from a box
#' (rather than just particles from a single galaxy), you can the half-mass
#' value at which the alignment function is run. Numeric length = 1. Default is
#' NA, in which case half the total mass of the supplied simulation data is
#' used. Specified in units of solar mass.
#'@param sph_spawn_n Numeric describing the number of gas particles with which
#' to sample the SPH smoothing length. Default is 1, which will not spawn
#' additional gas particles. Increasing this value increases the number of
#' particles used to model the gas distribution. This value may need to be
#' tested for convergence depending on the resolution of the grid used to image
#' the gas properties at the `build_datacube()` stage.
#'@return Returns an .Rdata file that contains a list of particle positions,
#' velocities, and spectral weights (or a list containing the same information
#' to the environment without writing to file, when `write_to_file = F`).
#'@examples
#'ss_file = make_simspin_file(filename = system.file("extdata",
#'                                                   "SimSpin_example_Gadget",
#'                                                    package = "SimSpin"),
#'                            write_to_file = FALSE)
#'

make_simspin_file = function(filename, cores=1, disk_age=5, bulge_age=10,
                             disk_Z=0.024, bulge_Z=0.001, template="BC03lr",
                             write_to_file=TRUE, output, overwrite = F,
                             centre=NA, half_mass=NA, sph_spawn_n=1){

  header = list("InputFile" = filename,
                "OutputFile" = NULL,
                "Type" = character(1),
                "Template" = character(1),
                "Template_LSF" = numeric(1),
                "Template_waveres" = numeric(1),
                "Centre" = centre,
                "HalfMass" = half_mass,
                "TotalStellarMass" = numeric(1),
                "TotalGasMass" = numeric(1),
                "Alignment" = "Default",
                "Npart" = numeric(1),
                "SmoothingN" = sph_spawn_n,
                "Origin" = paste0("SimSpin_v", packageVersion("SimSpin")),
                "Date" = Sys.time())

  temp_name = stringr::str_to_upper(template)

  if (write_to_file){
    if (missing(output)){
    output = paste(sub('\\..*', '', filename), "_", temp_name, ".Rdata", sep="")
    }
    if (file.exists(output) & !overwrite){
      stop(cat("FileExists Error:: SimSpin file already exists at: ", output, "\n",
               "If you wish to overwrite this file, please specify 'overwrite=T'. \n"))
    }
    header$OutputFile = output
  }

  if(temp_name == "BC03LR" | temp_name == "BC03"){
    temp = SimSpin::BC03lr
    header$Template = "BC03lr"
    header$Template_LSF = 3 # as according to Bruzual & Charlot (2003) MNRAS 344, pg 1000-1028
    header$Template_waveres = min(diff(temp$Wave))
  } else if (temp_name == "BC03HR"){
    temp = SimSpin::BC03hr
    header$Template = "BC03hr"
    header$Template_LSF = 3 # as according to Bruzual & Charlot (2003) MNRAS 344, pg 1000-1028
    header$Template_waveres = min(diff(temp$Wave))
  } else if (temp_name == "EMILES"){
    temp = SimSpin::EMILES
    header$Template = "EMILES"
    header$Template_LSF = 2.51 # as according to Vazdekis et al (2016) MNRAS 463, pg 3409-3436
    header$Template_waveres = min(diff(temp$Wave))
  } else {
    stop(cat("Error: template specified is unavailable.", "\n",
             "Please specify template = 'BC03', 'BC03lr', 'BC03hr' or 'EMILES'"))
  }

  if (sph_spawn_n%%1!=0){
    stop(cat("Error: sph_spawn_n must be a whole number.", "\n",
             "Please specify an integer value for sph_spawn_n."))
  }

  file_type = .get_file_type(filename)

  if (file_type == "hdf5"){
    galaxy_data = .read_hdf5(filename, cores)
  } else if (file_type == "gadget_binary") {
    galaxy_data = .read_gadget(filename)
  } else if (file_type == "tipsy_binary_big") {
    galaxy_data = .read_tipsy(filename, endian = "big", verbose=F)
  } else if (file_type == "tipsy_binary_little"){
    galaxy_data = .read_tipsy(filename, endian = "little", verbose=F)
  }

  header$Type = galaxy_data$head$Type

  Npart_sum = cumsum(galaxy_data$head$Npart) # Particle indices of each type

  header$Npart = sum(galaxy_data$head$Npart, na.rm=T)

  galaxy_data = .centre_galaxy(galaxy_data, centre) # centering the galaxy based on stellar particles

  if (header$Type != "nbody"){
    if (!is.na(half_mass)){ # if a half mass is requested, this alignment is specified by the user
      header$Alignment = "Specified"
    }
    galaxy_data = .align_galaxy(galaxy_data, half_mass) # align 3D shape of galaxy
    gc()
  } else {
    header$Alignment = "None"
    if (!is.na(half_mass[1])){
      warning(c("Input SimSpin file contains an N-body model. \n",
                "These models are not aligned by default so input half-mass will not be used for alignment. \n",
                "If you wish to align the model at the radius containing `half_mass`, modify the header element `Type` within the SimSpin file."))
    }
  }

  if (!is.null(galaxy_data$ssp)){ # if there are stellar particles in the file at all
                                  #  adding info to stellar_particle data

    # If a particle has a metallicity of 0, remove from sample
    Z0_int = which(galaxy_data$ssp$Metallicity == 0)
    if (length(Z0_int)!=0){ # remove if there are any Z = 0
      galaxy_data$star_part = galaxy_data$star_part[-c(Z0_int),]
      galaxy_data$ssp       = galaxy_data$ssp[-c(Z0_int),]
    }

    # ensure that Ages are described relative to the time "now"
    galaxy_data$ssp$Age = galaxy_data$ssp$Age - as.numeric(.SFTtoAge(galaxy_data$head$Time, cores = 1))

    # reassign any age==0 particles to have a small non-zero age
    galaxy_data$ssp$Age[galaxy_data$ssp$Age==0] = 1e-9

    # if there are greater than 1e6 stellar particles, bin stellar particles on age/Z position
    if (length(galaxy_data$star_part$ID) > 1e6){
      age_grid = 10^(seq(log10(min(galaxy_data$ssp$Age))-0.02, log10(max(galaxy_data$ssp$Age))+0.02, by = 0.02))
      Z_grid   = 10^(seq(log10(min(galaxy_data$ssp$Metallicity))-0.1, log10(max(galaxy_data$ssp$Metallicity))+0.1, by = 0.1))

      age_cen  = 10^(seq(log10(min(galaxy_data$ssp$Age))-0.01, log10(max(galaxy_data$ssp$Age))+0.01, by = 0.02))
      Z_cen    = 10^(seq(log10(min(galaxy_data$ssp$Metallicity))-0.05, log10(max(galaxy_data$ssp$Metallicity))+0.05, by = 0.1))

      AZ       = data.frame("ages" = rep(age_cen, length = length(age_cen)*length(Z_cen)),
                            "metallicities" = rep(Z_cen, each = length(age_cen)),
                            "id" = seq(1, length(age_cen)*length(Z_cen)))

      az_pos = cut(galaxy_data$ssp$Age, breaks = age_grid, labels = F) +
        (length(age_cen) * cut(galaxy_data$ssp$Metallicity, breaks = Z_grid, labels = F)) - length(age_cen)

      AZ_bins = AZ[sort(unique(az_pos)),]
      AZ_bins$sed_id = seq(1, length(AZ_bins$id))

      # adding info to stellar_particle data
      galaxy_data$star_part[, sed_id := AZ_bins$sed_id[match(az_pos,AZ_bins$id)], ]
      galaxy_data$star_part[, Metallicity := galaxy_data$ssp$Metallicity, ]
      galaxy_data$star_part[, Age := galaxy_data$ssp$Age, ]
      galaxy_data$star_part[, Initial_Mass := galaxy_data$ssp$Initial_Mass, ]

    } else {

      galaxy_data$star_part[, Metallicity := galaxy_data$ssp$Metallicity, ]
      galaxy_data$star_part[, Age := galaxy_data$ssp$Age, ]
      galaxy_data$star_part[, Initial_Mass := galaxy_data$ssp$Initial_Mass, ]
      galaxy_data$star_part[, sed_id := seq(1, length(galaxy_data$star_part$ID)), ]

    }

    sed  = .spectral_weights(Metallicity = galaxy_data$ssp$Metallicity,
                             Age = galaxy_data$ssp$Age,
                             Template = temp, cores = cores)
    # returns a list of spectra ids from the template set and the associated weights for A and Z

  } else if (header$Type == "nbody"){ # if working with N-body

    n_disk = galaxy_data$head$Npart[3]; n_bulge = galaxy_data$head$Npart[4] # number of disk and bulge particles
    n_stars = n_disk + n_bulge # total number of "stars"

    galaxy_data$star_part$Initial_Mass = galaxy_data$star_part$Mass[Npart_sum[2]+1:Npart_sum[4]]/2 # assuming the initial mass is half of the current mass
    galaxy_data$star_part$Age          = numeric(n_stars)
    galaxy_data$star_part$Metallicity  = numeric(n_stars)
    galaxy_data$star_part$sed_id       = numeric(n_stars)

    if (n_disk > 0 & n_bulge > 0){ # assigning ages and metallicities to disk and bulge particles (if present in snap)
      galaxy_data$star_part$Age[1:n_disk] = disk_age
      galaxy_data$star_part$Age[(n_disk+1):n_stars] = bulge_age
      galaxy_data$star_part$Metallicity[1:n_disk] = disk_Z
      galaxy_data$star_part$Metallicity[(n_disk+1):n_stars] = bulge_Z
      galaxy_data$star_part$sed_id[1:n_disk] = 1
      galaxy_data$star_part$sed_id[(n_disk+1):n_stars] = 2

      sed = .spectral_weights(Metallicity = c(disk_Z, bulge_Z),
                              Age = c(disk_age, bulge_age),
                              Template = temp, cores = cores)

    } else if (n_disk > 0 & n_bulge == 0){
      galaxy_data$star_part$Age = disk_age
      galaxy_data$star_part$Metallicity = disk_Z
      galaxy_data$star_part$sed_id = 1

      sed = .spectral_weights(Metallicity = disk_Z,
                              Age = disk_age,
                              Template = temp, cores = cores)

    } else if (n_disk == 0 & n_bulge > 0){
      galaxy_data$star_part$Age = bulge_age
      galaxy_data$star_part$Metallicity = bulge_Z
      galaxy_data$star_part$sed_id = 1

      sed = .spectral_weights(Metallicity = bulge_Z,
                              Age = bulge_age,
                              Template = temp, cores = cores)

      gc()

    }

  } else {sed = NULL} # if only gas in the file

  if (galaxy_data$head$Type == "EAGLE" |
      galaxy_data$head$Type == "Magneticum" |
      galaxy_data$head$Type == "Horizon-AGN" |
      galaxy_data$head$Type == "Illustris-TNG" |
      galaxy_data$head$Type == "Tipsy"){
    if (length(galaxy_data$gas_part$SmoothingLength)>0 & sph_spawn_n>1){ # if we need to spawn gas particles because we are working with SPH models

      gas_part_names = names(galaxy_data$gas_part)
      new_gas_part = stats::setNames(data.table::data.table(matrix(0,
                                                                   ncol = length(gas_part_names),
                                                                   nrow = (sph_spawn_n*galaxy_data$head$NumPart_Total[1]))),
                                     gas_part_names) # generate a new DF containing original gas_part names

      new_gas_part[, ID := seq(1, galaxy_data$head$NumPart_Total[1]*sph_spawn_n),]
      new_gas_part[, vx := rep(galaxy_data$gas_part$vx, each=sph_spawn_n),]
      new_gas_part[, vy := rep(galaxy_data$gas_part$vy, each=sph_spawn_n),]
      new_gas_part[, vz := rep(galaxy_data$gas_part$vz, each=sph_spawn_n),]
      new_gas_part[, SmoothingLength := rep(0, galaxy_data$head$NumPart_Total[1]*sph_spawn_n),]
      new_gas_part[, Metallicity := rep(galaxy_data$gas_part$Metallicity, each=sph_spawn_n),]
      #new_gas_part[, Carbon := rep(galaxy_data$gas_part$Carbon, each=sph_spawn_n),]
      new_gas_part[, Hydrogen := rep(galaxy_data$gas_part$Hydrogen, each=sph_spawn_n),]
      new_gas_part[, Oxygen := rep(galaxy_data$gas_part$Oxygen, each=sph_spawn_n),]

      kernel = character(1) # choosing the kernel relevant for the
      if (galaxy_data$head$Type == "EAGLE"){
        kernel = "WC2"
      } else if (galaxy_data$head$Type == "Tipsy"){
        kernel = "WC4"
      } else if (galaxy_data$head$Type == "Magneticum"){
        kernel = "WC6"
      } else if (galaxy_data$head$Type == "HorizonAGN" | galaxy_data$head$Type == "Illustris-TNG"){
        kernel = "M4"
      } else {
        kernel = "WC2"
      }

      if (cores > 1){
        galaxy_data$gas_part  = .sph_spawn_mc(galaxy_data$gas_part, new_gas_part, sph_spawn_n, kernel, cores)
      } else {
        galaxy_data$gas_part  = .sph_spawn(galaxy_data$gas_part, new_gas_part, sph_spawn_n, kernel)
      }

    }
  }

  header$TotalStellarMass = sum(galaxy_data$star_part$Mass, na.rm=T)
  header$TotalGasMass = sum(galaxy_data$gas_part$Mass, na.rm=T)

  simspin_file = list("header"    = header,
                      "star_part" = galaxy_data$star_part,
                      "gas_part"  = galaxy_data$gas_part,
                      "spectral_weights" = sed)

  if (write_to_file){
    saveRDS(simspin_file, file = output, compress = "xz")
    return(message("SimSpin file written to: ", output, "\n"))
  } else {
    return(simspin_file)
  }

}
