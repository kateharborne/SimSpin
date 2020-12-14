# Author: Kate Harborne
# Date: 22/10/2020
# Title: Testing the make_simspin_file.R code
#'Reformating isolated galaxy simulations to contain spectra.
#'
#'The purpose of this function is to construct a SimSpin file containing the
#' mock spectra for each particle contained within the galaxy simulation file.
#' If the snapshot provided is from a cosmological simulation, the SEDs
#' generated will be with respect to the Stellar Formation Time/Age, Metallicity
#' and Initial Mass of each stellar particle. If the system is an N-body model,
#' stellar particles are assumed to have an age and metallicity as provided to
#' the function as \code{disk_age}, \code{bulge_age}, \code{disk_Z} and
#' \code{bulge_Z}. Returned is an .fst file in a SimSpin readable format.
#'
#'@param filename The path to the snapshot file.
#'@param cores The number of cores across which to multi-thread the problem.
#'@param disk_age The age of the disk particles in Gyr.
#'@param bulge_age The age of the bulge particles in Gyr.
#'@param disk_Z The metallicity of the disk particles in Gyr.
#'@param bulge_Z The metallicity of the bulge particles in Gyr.
#'@param template The stellar templates from which to derive the SEDs. Options
#' include "BC03lr" (GALEXEV low resolution, Bruzual & Charlot 2003), "BC03hr"
#' (GALEXEV high resolution, Bruzual & Charlot 2003) or "EMILES" (Vazdekis et
#' al, 2016).
#'@param output The path at which the output file is written. If not provided,
#' file will be written at the location of the input filename with the addition
#' of "_spectra.fst".
#'@param overwrite If true, and the file already exists at the output location,
#' a new file will be written over the old one.
#'@param verbose Boolean to switch on code progress updates.
#'@return Returns an .fst file that contains a matrix of particle positions,
#' velocities, and spectra.
#'@examples
#'\dontrun{
#'make_simspin_file(filename = system.file("extdata", "SimSpin_example_EAGLE.hdf5",
#'                                          package = "SimSpin"),
#'                  output=tempfile())
#'}
#'

make_simspin_file = function(filename, cores=1, disk_age=5, bulge_age=10,
                             disk_Z=0.024, bulge_Z=0.001, template="BC03lr",
                             output, overwrite = F, verbose = F){

  temp_name = stringr::str_to_upper(template)

  if(missing(output)){
    output = paste(sub('\\..*', '', filename), "_", temp_name, ".Rdata", sep="")
  }

  if(file.exists(output) & !overwrite){
    stop(cat("FileExists Error:: SimSpin file already exists at: ", output, "\n",
               "If you wish to overwrite this file, please specify 'overwrite=T'. \n"))
  }

  if(temp_name == "BC03LR" | temp_name == "BC03"){
    temp = ProSpect::BC03lr
  } else if (temp_name == "BC03HR"){
    temp = ProSpect::BC03hr
  } else if (temp_name == "EMILES"){
    temp = ProSpect::EMILES
  } else {
    stop("Error: template specified is unavailable. \n Please specify template = 'BC03', 'BC03lr', 'BC03hr' or 'EMILES'")
  }

  galaxy_data = tryCatch(expr = {.read_gadget(filename, verbose)},
                         error = function(e){.read_hdf5(filename, cores, verbose)})

  Npart_sum = cumsum(galaxy_data$head$Npart) # Particle indices of each type

  galaxy_data = .centre_galaxy(galaxy_data) # centering the galaxy based on stellar particles

  galaxy_data = .align_galaxy(galaxy_data) # align 3D shape of galaxy

  if(!"ssp" %in% names(galaxy_data)){ # if the SSP field does not come from the snapshot file, must be working with N-body

    n_disk = galaxy_data$head$Npart[3]; n_bulge = galaxy_data$head$Npart[4] # number of disk and bulge particles
    n_stars = n_disk + n_bulge # total number of "stars"
    galaxy_data$ssp = data.frame("Initial_Mass"=numeric(n_stars), "Age"=numeric(n_stars),
                                 "Metallicity"=numeric(n_stars))
    galaxy_data$ssp$Initial_Mass = galaxy_data$star_part$Mass[Npart_sum[2]+1:Npart_sum[4]]/2 # assuming the initial mass is half of the current mass

    if (n_disk > 0 & n_bulge > 0){ # assigning ages and metallities to disk and bulge particles (if present in snap)
      galaxy_data$ssp$Age[1:n_disk] = disk_age
      galaxy_data$ssp$Age[(n_disk+1):n_stars] = bulge_age
      galaxy_data$ssp$Metallicity[1:n_disk] = disk_Z
      galaxy_data$ssp$Metallicity[(n_disk+1):n_stars] = bulge_Z
      AZ_bins = data.frame("ages" = c(disk_age, bulge_age),
                           "metallicities" = c(disk_Z, bulge_Z),
                           "id" = c(1,2),
                           "sed_id" = c(1,2))
      az_pos = c(rep(1, n_disk), rep(2, n_bulge))
    } else if (n_disk > 0 & n_bulge == 0){
      galaxy_data$ssp$Age = disk_age
      galaxy_data$ssp$Metallicity = disk_Z
      AZ_bins = data.frame("ages" = disk_age,
                           "metallicities" = disk_Z,
                           "id" = 1,
                           "sed_id" = 1)
      az_pos = rep(1, n_disk)
    } else if (n_disk == 0 & n_bulge > 0){
      galaxy_data$ssp$Age = bulge_age
      galaxy_data$ssp$Metallicity = bulge_Z
      AZ_bins = data.frame("ages" = bulge_age,
                           "metallicities" = bulge_Z,
                           "id" = 1,
                           "sed_id" = 1)
      az_pos = rep(1, n_bulge)
    }

  }

  # If a particle has a metallicity of 0, remove from sample?
  Z0_int = which(galaxy_data$ssp$Metallicity == 0)
  if (length(Z0_int)!=0){ # remove if there are any Z = 0
    galaxy_data$star_part = galaxy_data$star_part[-c(Z0_int),]
    galaxy_data$ssp       = galaxy_data$ssp[-c(Z0_int),]
  }

  # now binning stellar particles based on their A/Z position
  if (length(unique(galaxy_data$ssp$Age)) > 2){ # if not an N-body simulation
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

  }

  # adding info to stellar_particle data
  galaxy_data$star_part$sed_id = AZ_bins$sed_id[match(az_pos,AZ_bins$id)]
  galaxy_data$star_part$Initial_Mass = galaxy_data$ssp$Initial_Mass

  sed  = .spectra(Metallicity = AZ_bins$metallicities, Age = AZ_bins$ages, Template = temp, cores = cores) # returns a list

  simspin_file = list("star_part" = galaxy_data$star_part,
                      "spectra"   = sed,
                      "wave"      = temp$Wave)

  saveRDS(simspin_file, file = output)

  return(message("SimSpin file written to: ", output, "\n"))
}


