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

  if(missing(output)){
    output = paste(sub('\\..*', '', filename), "_spectra.fst", sep="")
  }

  if(file.exists(output) & !overwrite){
    stop(cat("FileExists Error:: SimSpin file already exists at: ", output, "\n",
               "If you wish to overwrite this file, please specify 'overwrite=T'. \n"))
  }

  if(stringr::str_to_upper(template) == "BC03LR" | stringr::str_to_upper(template) == "BC03"){
    temp = ProSpect::BC03lr
  } else if (stringr::str_to_upper(template) == "BC03HR"){
    temp = ProSpect::BC03hr
  } else if (stringr::str_to_upper(template) == "EMILES"){
    temp = ProSpect::EMILES
  } else {
    stop("Error: template specified is unavailable. \n Please specify template = 'BC03', 'BC03lr', 'BC03hr' or 'EMILES'")
  }

  galaxy_data = tryCatch(expr = {.read_gadget(filename, verbose)},
                         error = function(e){.read_hdf5(filename, cores, verbose)})

  Npart_sum = cumsum(galaxy_data$head$Npart) # Particle indices of each type

  galaxy_data$part = cen_galaxy(galaxy_data$part) # centering the galaxy

  galaxy_data = .align(galaxy_data) # align angular momentum vector with z-axis

  if(!"ssp" %in% names(galaxy_data)){ # if the SSP field does not come from the snapshot file, must be working with N-body

    n_disk = galaxy_data$head$Npart[3]; n_bulge = galaxy_data$head$Npart[4] # number of disk and bulge particles
    n_stars = n_disk + n_bulge # total number of "stars"
    galaxy_data$ssp = data.frame("Initial_Mass"=numeric(n_stars), "Age"=numeric(n_stars),
                                 "Metallicity"=numeric(n_stars))
    galaxy_data$ssp$Initial_Mass = galaxy_data$part$Mass[Npart_sum[2]+1:Npart_sum[4]]/2 # assuming the initial mass is half of the current mass

    if (n_disk > 0 & n_bulge > 0){ # assigning ages and metallities to disk and bulge particles (if present in snap)
      galaxy_data$ssp$Age[1:n_disk] = disk_age
      galaxy_data$ssp$Age[(n_disk+1):n_stars] = bulge_age
      galaxy_data$ssp$Metallicity[1:n_disk] = disk_Z
      galaxy_data$ssp$Metallicity[(n_disk+1):n_stars] = bulge_Z
    } else if (n_disk > 0 & n_bulge == 0){
      galaxy_data$ssp$Age = disk_age
      galaxy_data$ssp$Metallicity = disk_Z
    } else if (n_disk == 0 & n_bulge > 0){
      galaxy_data$ssp$Age = bulge_age
      galaxy_data$ssp$Metallicity = bulge_Z
    }

  }

  n_stars = length(galaxy_data$ssp$Age) # number of "stellar" particles
  wavelengths = length(temp$Wave)
  id_stars = seq(Npart_sum[2]+1, Npart_sum[5]) # ids of "stellar" particles
  simspin_file = matrix(data=NA, nrow=(7+wavelengths), ncol=n_stars+1)
  simspin_file[1,] = seq(0, n_stars)
  simspin_file[2,] = c(NA, galaxy_data$part$x[id_stars])
  simspin_file[3,] = c(NA, galaxy_data$part$y[id_stars])
  simspin_file[4,] = c(NA, galaxy_data$part$z[id_stars])
  simspin_file[5,] = c(NA, galaxy_data$part$vx[id_stars])
  simspin_file[6,] = c(NA, galaxy_data$part$vy[id_stars])
  simspin_file[7,] = c(NA, galaxy_data$part$vz[id_stars])
  simspin_file[8:(wavelengths+7),1] = temp$Wave
  simspin_file[8:(wavelengths+7),2:(n_stars+1)] = .part_spec(Metallicity = galaxy_data$ssp$Metallicity,
                                                             Age = galaxy_data$ssp$Age,
                                                             Mass = galaxy_data$ssp$Initial_Mass,
                                                             Template = temp, cores = cores)

  fst::write_fst(as.data.frame(simspin_file), path = output, compress = 100)

  message("SimSpin file written to: ", output, "\n")
}


