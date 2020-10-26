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
#'@return Returns an .fst file tat contains a matrix of particle positions,
#' velocities, and spectra.
#'@examples
#'make_simspin_file(filename = system.file("extdata", "SimSpin_example_Gadget",
#'                                          package = "SimSpin"),
#'                  output=tempfile())
#'


make_simspin_file = function(filename, cores=1, disk_age=5, bulge_age=10,
                             disk_Z=0.024, bulge_Z=0.001, template="BC03lr",
                             output, overwrite = F, ...){

  if(missing(output)){
    output = paste(sub('\\..*', '', filename), "_spectra.fst", sep="")
  }

  if(file.exists(output) & !overwrite){
    stop(cat("FileExists Error:: SimSpin file already exists at: ", output, "\n",
               "If you wish to overwrite this file, please specify 'overwrite=T'. \n"))
  }

  galaxy_data = tryCatch(expr = {.read_gadget(filename, ...)},
                         error = function(e){.read_hdf5(filename, cores, ...)})

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

  if(stringr::str_to_upper(template) == "BC03LR" | stringr::str_to_upper(template) == "BC03"){
    temp = ProSpect::BC03lr
  } else if (stringr::str_to_upper(template) == "BC03HR"){
    temp = ProSpect::BC03hr
  } else if (stringr::str_to_upper(template) == "EMILES"){
    temp = ProSpect::EMILES
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

  return(cat("SimSpin file written to: ", output, "\n"))
}

# Function for reading in Gadget binary files
.read_gadget = function(f, verbose = FALSE){
  data = file(f, "rb") # open file for reading in binary mode

  block         = readBin(data, "integer", n=1) #block size field, giving the length of the header
  if(block!=256){close(data); stop("Not a binary file. Trying HDF5 reader...")}
  if(verbose){cat("Reading in header.\n")}
  Npart         = readBin(data, "integer", n=6) # number of particles of each type in this file
  Massarr       = readBin(data, "numeric", n=6, size=8) # mass of each particle type. set to 0 for present particles means read mass in mass block.
  Time          = readBin(data, "numeric", n=1, size=8) # time of output, or expansion factor for cosmological simulations
  Redshift      = readBin(data, "numeric", n=1, size=8) # z = 1/(a-1)
  FlagSfr       = readBin(data, "integer", n=1) # flag for star formation
  FlagFeedback  = readBin(data, "integer", n=1) # flag for feedback
  Nall          = readBin(data, "integer", n=6) # total number of particles of each type in the simulation
  FlagCooling   = readBin(data, "integer", n=1) # flag for cooling
  NumFiles      = readBin(data, "integer", n=1) # number of files in each snapshot
  BoxSize       = readBin(data, "numeric", n=1, size=8) # box size if periodic boundary conditions are used
  Omega0        = readBin(data, "numeric", n=1, size=8) # matter density at z=0 in units of critical density
  OmegaLambda   = readBin(data, "numeric", n=1, size=8) # vacuum energy density at z=0 in units of critical density
  HubbleParam   = readBin(data, "numeric", n=1, size=8) # the Hubble constant in units of 100 kms^-1Mpc^-1
  FlagAge       = readBin(data, "integer", n=1) # Creation time of stars
  FlagMetals    = readBin(data, "integer", n=1) # Metalliticy values
  NallHW        = readBin(data, "integer", n=6) # For simulations with more than 2^32 particles
  flag_entr_ics = readBin(data, "integer", n=1) # flags that ICs contain entropy rather than thermal energy in the U block
  empty         = readBin(data, "integer", n=15) # unused empty bytes at the end of the header
  block         = readBin(data, "integer", n=1)

  # Reading in the positions block
  block         = readBin(data, "integer", n=1)
  if(verbose){cat("Reading in positions for", block/4/3, "particles. \n")}
  pos           = readBin(data, "numeric", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  # Reading in the velocities block
  block         = readBin(data, "integer", n=1)
  if(verbose){cat("Reading in velocities for", block/4/3, "particles. \n")}
  vel           = readBin(data, "numeric", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  # Reading in the ID's block
  block         = readBin(data, "integer", n=1)
  if(verbose){cat("Reading in IDs for", block/4, "particles. \n")}
  id            = readBin(data, "integer", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  # Reading in the mass block
  block         = readBin(data, "integer", n=1)
  if(verbose){cat("Reading in masses for", block/4, "particles. \n")}
  masses        = readBin(data, "numeric", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  close(data)

  extract = ((1:floor(sum(Npart)))*3)-2 # giving integers of x/vx within pos/vel
  part = data.frame("ID" = id,
                    "x" = pos[extract], "y"=pos[extract+1], "z"=pos[extract+2],
                    "vx" = vel[extract], "vy" = vel[extract+1], "vz"=vel[extract+2],
                    "Mass" = masses)
  head = list("Npart" = Npart, "Massarr" = Massarr, "Time" = Time, "Redshift" = Redshift,
              "FlagSfr" = FlagSfr, "FlagFeedback" = FlagFeedback, "Nall" = Nall,
              "FlagCooling" = FlagCooling, "NumFiles" = NumFiles, "BoxSize" = BoxSize,
              "Omega0" = Omega0, "OmegaLambda" = OmegaLambda, "HubbleParam" = HubbleParam,
              "FlagAge" = FlagAge, "FlagMetals" = FlagMetals, "NallHW" = NallHW,
              "flag_entr_ics" = flag_entr_ics)

  if(verbose){cat("Done reading Gadget snapshot file. \n")}
  return(list(part=part, head=head))

}

# Function for reading in Gadget HDF5 files and EAGLE HDF5 files
.read_hdf5   = function(f, cores, verbose = FALSE){

  data          = hdf5r::h5file(f, mode="r")

  if(length(hdf5r::list.attributes(data[["Header"]])) > 17){eagle = T}else{eagle=F} # determining if EAGLE input (based on number of parameters in Header)
  if(verbose & eagle){cat("EAGLE input snapshot detected.\n")}

  if(verbose){cat("Reading in header.\n")}
  Npart         = hdf5r::h5attr(data[["Header"]], "NumPart_ThisFile")
  Massarr       = hdf5r::h5attr(data[["Header"]], "MassTable")
  Time          = hdf5r::h5attr(data[["Header"]], "Time")
  Redshift      = hdf5r::h5attr(data[["Header"]], "Redshift")
  FlagSfr       = hdf5r::h5attr(data[["Header"]], "Flag_Sfr")
  FlagFeedback  = hdf5r::h5attr(data[["Header"]], "Flag_Feedback")
  Nall          = hdf5r::h5attr(data[["Header"]], "NumPart_Total")
  FlagCooling   = hdf5r::h5attr(data[["Header"]], "Flag_Cooling")
  NumFiles      = hdf5r::h5attr(data[["Header"]], "NumFilesPerSnapshot")
  BoxSize       = hdf5r::h5attr(data[["Header"]], "BoxSize")
  Omega0        = hdf5r::h5attr(data[["Header"]], "Omega0")
  OmegaLambda   = hdf5r::h5attr(data[["Header"]], "OmegaLambda")
  HubbleParam   = hdf5r::h5attr(data[["Header"]], "HubbleParam")
  FlagAge       = hdf5r::h5attr(data[["Header"]], "Flag_StellarAge")
  FlagMetals    = hdf5r::h5attr(data[["Header"]], "Flag_Metals")
  NallHW        = hdf5r::h5attr(data[["Header"]], "NumPart_Total_HighWord")
  if (eagle){
    flag_entr_ics = 0L
  } else {
    flag_entr_ics = hdf5r::h5attr(data[["Header"]], "Flag_Entropy_ICs")
  }

  head = list("Npart" = Npart, "Massarr" = Massarr, "Time" = Time, "Redshift" = Redshift,
              "FlagSfr" = FlagSfr, "FlagFeedback" = FlagFeedback, "Nall" = Nall,
              "FlagCooling" = FlagCooling, "NumFiles" = NumFiles, "BoxSize" = BoxSize,
              "Omega0" = Omega0, "OmegaLambda" = OmegaLambda, "HubbleParam" = HubbleParam,
              "FlagAge" = FlagAge, "FlagMetals" = FlagMetals, "NallHW" = NallHW,
              "flag_entr_ics" = flag_entr_ics)

  present_types  = which(Npart > 0) # Which PartType groups will be present in the file?
  Npart_sum      = cumsum(Nall) # Particle indices of each type
  Ntotal         = sum(Nall) # total number of all particle types
  mass_excpt     = which(Massarr > 0) # are any masses listed in the header?

  part = data.frame("ID" = numeric(Ntotal),
                    "x" = numeric(Ntotal), "y" = numeric(Ntotal), "z" = numeric(Ntotal),
                    "vx" = numeric(Ntotal), "vy" = numeric(Ntotal), "vz" = numeric(Ntotal),
                    "Mass" = numeric(Ntotal))

  pos = array(NA, dim=c(3, Ntotal))
  vel = array(NA, dim=c(3, Ntotal))

  if(verbose){cat("Reading particle properties. \n")}

  for (i in 1:length(present_types)){
    if (i == 1){
      i_start = 1; i_end = Npart_sum[present_types[i]]
      } else {
      i_start = Npart_sum[present_types[i-1]]+1; i_end = Npart_sum[present_types[i]]
      } # computing the indices of the particles of each present PartType
    pos[,i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_types[i]-1, "/Coordinates", sep="")]]) # reading x, y, z coordinates
    if (eagle){vel[,i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_types[i]-1, "/Velocity", sep="")]])}else{
      vel[,i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_types[i]-1, "/Velocities", sep="")]])
    } # reading vx, vy, vz velocities (called "Velocity" in EAGLE snapshots and "Velocities" Gadget HDF5 files)
    part$ID[i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_types[i]-1, "/ParticleIDs", sep="")]]) # reading particle IDs
    if (length(mass_excpt)!=0 & present_types[i] %in% mass_excpt){ # if masses have been specified in header, use these
      part$Mass[i_start:i_end] = Massarr[present_types[i]]
    } else { # else, read masses from PartType ("Mass" in EAGLE snapshots and "Masses" in Gadget HDF5 files)
      if (eagle){part$Mass[i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_types[i]-1, "/Mass", sep="")]])}else{
        part$Mass[i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_types[i]-1, "/Masses", sep="")]])
      }
    }
  }

  if (eagle){ # reading the details from EAGLE files for simple stellar population
    ssp = list("Initial_Mass"=numeric(length=Nall[5]), "Age"=numeric(length = Nall[5]),
               "Metallicity"=numeric(length=Nall[5]))
    ssp$Initial_Mass = hdf5r::readDataSet(data[["PartType4/InitialMass"]])
    Stellar_Formation_Time = hdf5r::readDataSet(data[["PartType4/StellarFormationTime"]])
    ssp$Age = as.numeric(.SFTtoAge(a = Stellar_Formation_Time, cores = cores))
    ssp$Metallicity = hdf5r::readDataSet(data[["PartType4/SmoothedMetallicity"]])

    a_coor = hdf5r::h5attr(data[["PartType4/Coordinates"]], "aexp-scale-exponent")
    h_coor = hdf5r::h5attr(data[["PartType4/Coordinates"]], "h-scale-exponent")
    a_velo = hdf5r::h5attr(data[["PartType4/Velocity"]], "aexp-scale-exponent")
    h_velo = hdf5r::h5attr(data[["PartType4/Velocity"]], "h-scale-exponent")
    a_mass = hdf5r::h5attr(data[["PartType4/Mass"]], "aexp-scale-exponent")
    h_mass = hdf5r::h5attr(data[["PartType4/Mass"]], "h-scale-exponent")

    pos = pos * head$Time^(a_coor) * head$HubbleParam^(h_coor) * 1e3
    vel = vel * head$Time^(a_velo) * head$HubbleParam^(h_velo)
    part$Mass = part$Mass * head$Time^(a_mass) * head$HubbleParam^(h_mass)
  }

  part$x = pos[1,]; part$vx = vel[1,]
  part$y = pos[2,]; part$vy = vel[2,]
  part$z = pos[3,]; part$vz = vel[3,]

  hdf5r::h5close(data)

  if(verbose){cat("Done reading HDF5 snapshot file. \n")}
  if(eagle){return(list(part=part, head=head, ssp=ssp))}else{return(list(part=part, head=head))}

}

# Function for computing the stellar age from the formation time in parallel
.SFTtoAge = function(a, cores=1){
  c1 = snow::makeCluster(cores)
  doSNOW::registerDoSNOW(c1)
  output = foreach(i = 1:length(a), .packages = "celestial")%dopar%{celestial::cosdistTravelTime((1 / a[i]) - 1)}
  closeAllConnections()

  return(output)
}

# Function to generate spectra
.part_spec = function(Metallicity, Age, Mass, Template, cores){
  c1 = snow::makeCluster(cores)
  doSNOW::registerDoSNOW(c1)
  part_spec = foreach(i = 1:length(Metallicity), .packages = c("ProSpect", "SimSpin"))%dopar%{
    Z = as.numeric(ProSpect::interp_quick(Metallicity[i], Template$Z, log = TRUE))
    A = as.numeric(ProSpect::interp_quick(Age[i]*1e9, Template$Age, log = TRUE))

    weights = data.frame("hihi" = Z[4] * A[4],   # ID_lo = 1, ID_hi = 2, wt_lo = 3, wt_hi = 4
                         "hilo" = Z[4] * A[3],
                         "lohi" = Z[3] * A[4],
                         "lolo" = Z[3] * A[3])

    part_spec = array(data = NA, dim = c(1, length(Template$Wave)))

    part_spec = ((Template$Zspec[[Z[2]]][A[2],] * weights$hihi) +
                 (Template$Zspec[[Z[2]]][A[1],] * weights$hilo) +
                 (Template$Zspec[[Z[1]]][A[2],] * weights$lohi) +
                 (Template$Zspec[[Z[1]]][A[1],] * weights$lolo)) * Mass[i] * 1e10

    return(part_spec)
  }
  closeAllConnections()
  output = matrix(unlist(part_spec, use.names=FALSE), ncol = length(Metallicity), byrow = FALSE)

  return(output)
}

# Function to flip galaxy if Jz is upside-down
.flip = function(galaxy_data){

    rot_mat = matrix(c(1,0,0,0,-1,0,0,0,-1), nrow = 3, ncol=3)

    new_coor =  rot_mat %*% rbind(galaxy_data$part$x, galaxy_data$part$y, galaxy_data$part$z)
    new_vel =  rot_mat %*% rbind(galaxy_data$part$vx, galaxy_data$part$vy, galaxy_data$part$vz)

    galaxy_data$part$x = new_coor[1,]; galaxy_data$part$y = new_coor[2,]; galaxy_data$part$z = new_coor[3,];
    galaxy_data$part$vx = new_vel[1,]; galaxy_data$part$vy = new_vel[2,]; galaxy_data$part$vz = new_vel[3,];

    return(galaxy_data)
}

# Function to align galaxy
.align = function(galaxy_data){

  Npart_sum = cumsum(galaxy_data$head$Npart) # Particle indices of each type

  # To align the galaxy edge-on, align the J vector of inner 10kpc of disk with
  # the z axis. Inner 10kc of disk is chosen using preferentially gas, then stars,
  # then disk particles (in the case that the former does not exist in the sim).
  no_angmom = galaxy_data$head$Npart[1] # number of gas particles for computing angmom
  angmom_ids = c(1, no_angmom) # ids of gas particles
  if (no_angmom == 0){
    no_angmom = galaxy_data$head$Npart[5] # number of stellar particles
    angmom_ids = c(Npart_sum[4]+1, Npart_sum[5]) # ids of star particles
    if (no_angmom == 0){
      no_angmom = galaxy_data$head$Npart[3] # number of disk particles
      angmom_ids = c(Npart_sum[2]+1,Npart_sum[3]) # ids of disk particles
    }
  }
  # pulling out the disk from which to compute J
  r = r_galaxy(galaxy_data$part[angmom_ids[1]:angmom_ids[2],]) # compute radial coordinates
  J = angmom_galaxy(galaxy_data$part[angmom_ids[1]:angmom_ids[2],][r < 0.33*max(r),]) # compute J
  J_norm = matrix(J/(sqrt(J[1]^2 + J[2]^2 + J[3]^2)), nrow=1, ncol=3)

  v = c(J_norm[2], -J_norm[1], 0) # unit vector normal to J and z-axis, about which we want to rotate
  c = J_norm[3] # giving cos(angle)
  s = sqrt(v[1]^2 + v[2]^2) # giving sin(angle)

  if (J_norm[3] == -1){
    galaxy_data = .flip(galaxy_data)
    return(galaxy_data)
  }
  if (J_norm[3] == 1){
    return(galaxy_data)
  }

  v_x = matrix(data = c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), nrow = 3, ncol = 3) # skew-symmetric cross product
  I = matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3) # identity matrix
  rot_mat = I + v_x + (1/(1+c))*(v_x %*% v_x) # rotation matrix via Rodrigues Rotation Formula: wikipedia.org/wiki/Rodrigues'_rotation_formula

  new_coor =  rot_mat %*% rbind(galaxy_data$part$x, galaxy_data$part$y, galaxy_data$part$z)
  new_vel =  rot_mat %*% rbind(galaxy_data$part$vx, galaxy_data$part$vy, galaxy_data$part$vz)

  galaxy_data$part$x = new_coor[1,]; galaxy_data$part$y = new_coor[2,]; galaxy_data$part$z = new_coor[3,];
  galaxy_data$part$vx = new_vel[1,]; galaxy_data$part$vy = new_vel[2,]; galaxy_data$part$vz = new_vel[3,];

  return(galaxy_data)
}
