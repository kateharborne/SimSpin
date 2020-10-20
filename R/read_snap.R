# Kate Harborne - 16/10/2020
# SimSpin v2.0.0 - read_snapshot() function

library(hdf5r)

read_snap = function(filename, ...){

  galaxy_data = tryCatch(expr = {.read_gadget(filename, ...)},
                         error = function(e){.read_hdf5(filename, ...)})

  return(galaxy_data)
}

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

.read_hdf5   = function(f, verbose = FALSE){

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

  part$x = pos[1,]; part$vx = vel[1,]
  part$y = pos[2,]; part$vy = vel[2,]
  part$z = pos[3,]; part$vz = vel[3,]

  if (eagle){ # reading the details from EAGLE files for simple stellar population
    ssp = list("Initial_Mass"=numeric(length=Nall[5]), "Stellar_Formation_Time"=numeric(length = Nall[5]),
               "Metallicity"=numeric(length=Nall[5]))
    ssp$Initial_Mass = hdf5r::readDataSet(data[["PartType4/InitialMass"]])
    ssp$Stellar_Formation_Time = hdf5r::readDataSet(data[["PartType4/StellarFormationTime"]])
    ssp$Metallicity = hdf5r::readDataSet(data[["PartType4/SmoothedMetallicity"]])
  }

  hdf5r::h5close(data)

  if(verbose){cat("Done reading HDF5 snapshot file. \n")}
  if(eagle){return(list(part=part, head=head, ssp=ssp))}else{return(list(part=part, head=head))}

}
