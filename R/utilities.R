# Author: Kate Harborne
# Date: 27/10/2020
# Title: Utilities functions (i.e. hidden functions from the user)

# Some useful constants:
.lsol_to_erg    = 3.828e33
.mpc_to_cm      = 3.08568e+24
.speed_of_light = 299792.458
.cm_to_kpc = 3.24078e-22
.cms_to_kms = 1e-5
.g_to_msol = 5.02785e-34
.gcm3_to_msolkpc3 = 1.477e+31
.g_constant_cgs = 6.67430e-11
.g_in_kpcMsolkms2 = 4.3009e-6

# Functions for computing weighted means
.meanwt = function(x,wt){
  return(sum(x*wt, na.rm=T)/sum(wt,na.rm=T))
} # weighted mean

.varwt = function(x, wt, xcen){
  if (missing(xcen)){xcen = .meanwt(x,wt)}
  return(sum(wt*(x - xcen)^2, na.rm=T)/sum(wt, na.rm=T))
} # weighted variance

# A function for combining multiple results from a parallel loop
.comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Function for reading in Gadget binary files
.read_gadget = function(f){
  data = file(f, "rb") # open file for reading in binary mode

  block         = readBin(data, "integer", n=1) #block size field, giving the length of the header
  if(block!=256){close(data); stop("Not a binary file. Trying HDF5 reader...")}
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
  pos           = readBin(data, "numeric", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  # Reading in the velocities block
  block         = readBin(data, "integer", n=1)
  vel           = readBin(data, "numeric", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  # Reading in the ID's block
  block         = readBin(data, "integer", n=1)
  id            = readBin(data, "integer", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  # Reading in the mass block
  block         = readBin(data, "integer", n=1)
  masses        = readBin(data, "numeric", n=block/4, size=4)
  block         = readBin(data, "integer", n=1)

  close(data)

  extract = ((1:floor(sum(Npart)))*3)-2 # giving integers of x/vx within pos/vel

  part = data.frame("ID" = id,          # the particle data table
                    "x" = pos[extract], "y"=pos[extract+1], "z"=pos[extract+2],
                    "vx" = vel[extract], "vy" = vel[extract+1], "vz"=vel[extract+2],
                    "Mass" = masses*1e10)

  head = list("Npart" = c(Npart[1], 0, Npart[3], Npart[4], Npart[5], 0), # number of gas and stars
              "Time" = Time, "Redshift" = Redshift, # relevent simulation data
              "Nall" = Nall, "Type"="nbody") # number of particles in the original file


  Npart_sum = cumsum(Npart) # cumulative number of each particle type

  star_part = part[(Npart_sum[2]+1):Npart_sum[5],]

  if (Npart[1] > 0){
    gas_part = part[1:Npart_sum[1],]
  } else {gas_part = NULL}

  return(list(star_part = star_part, gas_part = gas_part, head = head))

}

# Functions for reading in HDF5 files
.read_hdf5   = function(f, cores=1){

  data = hdf5r::h5file(f, mode="r")

  # Read in all attributes listed in the header
  header_attr = hdf5r::list.attributes(data[["Header"]])

  if(length(header_attr) < 23){gadget2 = T}else{gadget2=F}
  if(length(header_attr) == 27){eagle = T}else{eagle=F} # determining if EAGLE input (based on number of parameters in Header)
  if(length(header_attr) == 23){magneticum = T}else{magneticum=F}

  # Create a list to store each variable
  head = vector("list", length(header_attr))
  names(head) = header_attr

  # Read in each variable and store in list
  for (i in 1:length(header_attr)){
    head[[i]] = hdf5r::h5attr(data[["Header"]], paste0(header_attr[i]))
  }

  # Read particle data differently depending on the simulation being read in...
  if (gadget2){output = .gadget2_read_hdf5(data, head)}
  if (eagle){output = .eagle_read_hdf5(data, head, cores)}
  if (magneticum){output = .magneticum_read_hdf5(data, head, cores)}

  hdf5r::h5close(data)

  return(output)
}

.gadget2_read_hdf5 = function(data, head){

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if ("PartType0" %in% groups){ # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])
    n_gas_prop = length(PT0_attr)
    gas = vector("list", n_gas_prop)
    names(gas) = PT0_attr

    for (i in 1:n_gas_prop){
      gas[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType0/",PT0_attr[i])]])
    }

    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.frame("ID" = gas$ParticleIDs,
                          "x"  = if(one_p_flag){gas$Coordinates[1]}else{gas$Coordinates[1,]},# Coordinates in kpc
                          "y"  = if(one_p_flag){gas$Coordinates[2]}else{gas$Coordinates[2,]},
                          "z"  = if(one_p_flag){gas$Coordinates[3]}else{gas$Coordinates[3,]},
                          "vx"  = if(one_p_flag){gas$Velocity[1]}else{gas$Velocity[1,]}, # Velocities in km/s
                          "vy"  = if(one_p_flag){gas$Velocity[2]}else{gas$Velocity[2,]},
                          "vz"  = if(one_p_flag){gas$Velocity[3]}else{gas$Velocity[3,]},
                          "Mass" = gas$Mass*1e10, # Mass in solar masses
                          "SFR" = gas$StarFormationRate,
                          "Density" = gas$Density, # Density in Msol/kpc^3
                          "Temperature" = gas$Temperature,
                          "SmoothingLength" = gas$SmoothingLength) # Smoothing length in kpc

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType2" %in% groups & "PartType3" %in% groups){ # If both bulge and disk stars are present

    PT2_attr = hdf5r::list.datasets(data[["PartType2"]])
    n_star_prop = length(PT2_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT2_attr

    for (i in 1:n_star_prop){
        stars[[i]] = hdf5r::readDataSet(data[[paste0("PartType2/",PT2_attr[i])]])
    }

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    disk_part = data.frame("ID" = stars$ParticleIDs,
                           "x"  = if(one_p_flag){stars$Coordinates[1]}else{stars$Coordinates[1,]}, # Coordinates in kpc
                           "y"  = if(one_p_flag){stars$Coordinates[2]}else{stars$Coordinates[2,]},
                           "z"  = if(one_p_flag){stars$Coordinates[3]}else{stars$Coordinates[3,]},
                           "vx"  = if(one_p_flag){stars$Velocities[1]}else{stars$Velocities[1,]}, # Velocities in km/s
                           "vy"  = if(one_p_flag){stars$Velocities[2]}else{stars$Velocities[2,]},
                           "vz"  = if(one_p_flag){stars$Velocities[3]}else{stars$Velocities[3,]},
                           "Mass" = stars$Masses*1e10) # Mass in solar masses

    remove(stars); remove(PT2_attr)

    PT3_attr = hdf5r::list.datasets(data[["PartType3"]])
    n_star_prop = length(PT3_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT3_attr

    for (i in 1:n_star_prop){
      stars[[i]] = hdf5r::readDataSet(data[[paste0("PartType3/",PT3_attr[i])]])
    }

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.frame("ID" = c(disk_part$ID, stars$ParticleIDs),
                           "x"  = c(disk_part$x, if(one_p_flag){stars$Coordinates[1]}else{stars$Coordinates[1,]}), # Coordinates in kpc
                           "y"  = c(disk_part$y, if(one_p_flag){stars$Coordinates[2]}else{stars$Coordinates[2,]}),
                           "z"  = c(disk_part$z, if(one_p_flag){stars$Coordinates[3]}else{stars$Coordinates[3,]}),
                           "vx"  = c(disk_part$vx, if(one_p_flag){stars$Velocities[1]}else{stars$Velocities[1,]}), # Velocities in km/s
                           "vy"  = c(disk_part$vy, if(one_p_flag){stars$Velocities[2]}else{stars$Velocities[2,]}),
                           "vz"  = c(disk_part$vz, if(one_p_flag){stars$Velocities[3]}else{stars$Velocities[3,]}),
                           "Mass" = c(disk_part$Mass, stars$Masses*1e10)) # Mass in solar masses

    remove(stars); remove(PT3_attr); remove(disk_part)

  } else if ("PartType2" %in% groups){
    PT2_attr = hdf5r::list.datasets(data[["PartType2"]])
    n_star_prop = length(PT2_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT2_attr

    for (i in 1:n_star_prop){
      stars[[i]] = hdf5r::readDataSet(data[[paste0("PartType2/",PT2_attr[i])]])
    }

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.frame("ID" = stars$ParticleIDs,
                           "x"  = if(one_p_flag){stars$Coordinates[1]}else{stars$Coordinates[1,]}, # Coordinates in kpc
                           "y"  = if(one_p_flag){stars$Coordinates[2]}else{stars$Coordinates[2,]},
                           "z"  = if(one_p_flag){stars$Coordinates[3]}else{stars$Coordinates[3,]},
                           "vx"  = if(one_p_flag){stars$Velocities[1]}else{stars$Velocities[1,]}, # Velocities in km/s
                           "vy"  = if(one_p_flag){stars$Velocities[2]}else{stars$Velocities[2,]},
                           "vz"  = if(one_p_flag){stars$Velocities[3]}else{stars$Velocities[3,]},
                           "Mass" = stars$Masses*1e10) # Mass in solar masses

    remove(stars); remove(PT2_attr)
  } else if ("PartType3" %in% groups){

    PT3_attr = hdf5r::list.datasets(data[["PartType3"]])
    n_star_prop = length(PT3_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT3_attr

    for (i in 1:n_star_prop){
      stars[[i]] = hdf5r::readDataSet(data[[paste0("PartType3/",PT3_attr[i])]])
    }

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.frame("ID" = stars$ParticleIDs,
                           "x"  = if(one_p_flag){stars$Coordinates[1]}else{stars$Coordinates[1,]}, # Coordinates in kpc
                           "y"  = if(one_p_flag){stars$Coordinates[2]}else{stars$Coordinates[2,]},
                           "z"  = if(one_p_flag){stars$Coordinates[3]}else{stars$Coordinates[3,]},
                           "vx"  = if(one_p_flag){stars$Velocities[1]}else{stars$Velocities[1,]}, # Velocities in km/s
                           "vy"  = if(one_p_flag){stars$Velocities[2]}else{stars$Velocities[2,]},
                           "vz"  = if(one_p_flag){stars$Velocities[3]}else{stars$Velocities[3,]},
                           "Mass" = stars$Masses*1e10) # Mass in solar masses

    remove(stars); remove(PT3_attr);

  } else {star_part = NULL}

  Npart = head$NumPart_ThisFile
  head = list("Npart" = c(Npart[1], 0, Npart[3], Npart[4], Npart[5], 0), # number of gas and stars
              "Time" = head$Time, "Redshift" = head$Redshift, # relevent simulation data
              "Nall" = head$NumPart_Total, "Type"="nbody") # number of particles in the original file

  return(list(star_part=star_part, gas_part=gas_part, head=head))

}

.eagle_read_hdf5 = function(data, head, cores){

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if ("PartType0" %in% groups){ # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])
    n_gas_prop = length(PT0_attr)
    gas = vector("list", n_gas_prop)
    names(gas) = PT0_attr

    for (i in 1:n_gas_prop){
      aexp = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "aexp-scale-exponent")
      hexp = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "h-scale-exponent")
      cgs  = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "CGSConversionFactor")
      gas[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType0/",PT0_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.frame("ID" = gas$ParticleIDs,
                          "x"  = if(one_p_flag){gas$Coordinates[1]*.cm_to_kpc}else{gas$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                          "y"  = if(one_p_flag){gas$Coordinates[2]*.cm_to_kpc}else{gas$Coordinates[2,]*.cm_to_kpc},
                          "z"  = if(one_p_flag){gas$Coordinates[3]*.cm_to_kpc}else{gas$Coordinates[3,]*.cm_to_kpc},
                          "vx"  = if(one_p_flag){gas$Velocity[1]*.cms_to_kms}else{gas$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                          "vy"  = if(one_p_flag){gas$Velocity[2]*.cms_to_kms}else{gas$Velocity[2,]*.cms_to_kms},
                          "vz"  = if(one_p_flag){gas$Velocity[3]*.cms_to_kms}else{gas$Velocity[3,]*.cms_to_kms},
                          "Mass" = gas$Mass*.g_to_msol, # Mass in solar masses
                          "SFR" = gas$StarFormationRate,
                          "Density" = gas$Density*.gcm3_to_msolkpc3, # Density in Msol/kpc^3
                          "Temperature" = gas$Temperature,
                          "SmoothingLength" = gas$SmoothingLength*.cm_to_kpc, # Smoothing length in kpc
                          "Metallicity" = gas$SmoothedMetallicity,
                          "Carbon" = gas$`SmoothedElementAbundance/Carbon`,
                          "Hydrogen" = gas$`SmoothedElementAbundance/Hydrogen`,
                          "Oxygen" = gas$`SmoothedElementAbundance/Oxygen`)

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType4" %in% groups){
    PT4_attr = hdf5r::list.datasets(data[["PartType4"]])
    n_star_prop = length(PT4_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT4_attr

    for (i in 1:n_star_prop){
      aexp = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "aexp-scale-exponent")
      hexp = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "h-scale-exponent")
      cgs  = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "CGSConversionFactor")
      stars[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType4/",PT4_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.frame("ID" = stars$ParticleIDs,
                           "x"  = if(one_p_flag){stars$Coordinates[1]*.cm_to_kpc}else{stars$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                           "y"  = if(one_p_flag){stars$Coordinates[2]*.cm_to_kpc}else{stars$Coordinates[2,]*.cm_to_kpc},
                           "z"  = if(one_p_flag){stars$Coordinates[3]*.cm_to_kpc}else{stars$Coordinates[3,]*.cm_to_kpc},
                           "vx"  = if(one_p_flag){stars$Velocity[1]*.cms_to_kms}else{stars$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                           "vy"  = if(one_p_flag){stars$Velocity[2]*.cms_to_kms}else{stars$Velocity[2,]*.cms_to_kms},
                           "vz"  = if(one_p_flag){stars$Velocity[3]*.cms_to_kms}else{stars$Velocity[3,]*.cms_to_kms},
                           "Mass" = stars$Mass*.g_to_msol) # Mass in solar masses

    ssp = data.frame("Initial_Mass" = stars$InitialMass*.g_to_msol,
                     "Age" = as.numeric(.SFTtoAge(a = stars$StellarFormationTime, cores = cores)),
                     "Metallicity" = stars$SmoothedMetallicity)

    remove(stars); remove(PT4_attr)

  } else {star_part=NULL; ssp=NULL}

  head$Type = "EAGLE"
  return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))

}

.magneticum_read_hdf5 = function(data, head, cores){

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if ("PartType0" %in% groups){ # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])
    n_gas_prop = length(PT0_attr)
    gas = vector("list", n_gas_prop)
    names(gas) = PT0_attr

    for (i in 1:n_gas_prop){
      aexp = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "aexp-scale-exponent")
      hexp = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "h-scale-exponent")
      cgs  = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "CGSConversionFactor")
      gas[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType0/",PT0_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.frame("ID" = gas$ParticleIDs,
                          "x"  = if(one_p_flag){gas$Coordinates[1]*.cm_to_kpc}else{gas$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                          "y"  = if(one_p_flag){gas$Coordinates[2]*.cm_to_kpc}else{gas$Coordinates[2,]*.cm_to_kpc},
                          "z"  = if(one_p_flag){gas$Coordinates[3]*.cm_to_kpc}else{gas$Coordinates[3,]*.cm_to_kpc},
                          "vx"  = if(one_p_flag){gas$Velocity[1]*.cms_to_kms}else{gas$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                          "vy"  = if(one_p_flag){gas$Velocity[2]*.cms_to_kms}else{gas$Velocity[2,]*.cms_to_kms},
                          "vz"  = if(one_p_flag){gas$Velocity[3]*.cms_to_kms}else{gas$Velocity[3,]*.cms_to_kms},
                          "Mass" = gas$Mass*.g_to_msol, # Mass in solar masses
                          "SFR" = gas$StarFormationRate,
                          "Density" = gas$Density*.gcm3_to_msolkpc3, # Density in Msol/kpc^3
                          "Temperature" = gas$Temperature,
                          "SmoothingLength" = gas$SmoothingLength*.cm_to_kpc, # Smoothing length in kpc
                          "Metallicity" = (colSums(gas$Metallicity[2:11,]))/(gas$Mass),
                          "Carbon" = gas$Metallicity[2,] / gas$Mass,
                          "Hydrogen" = (gas$Mass - colSums(gas$Metallicity)) / gas$Mass,
                          "Oxygen" =  gas$Metallicity[4,] / gas$Mass)

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType4" %in% groups){
    PT4_attr = hdf5r::list.datasets(data[["PartType4"]])
    n_star_prop = length(PT4_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT4_attr

    for (i in 1:n_star_prop){
      aexp = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "aexp-scale-exponent")
      hexp = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "h-scale-exponent")
      cgs  = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "CGSConversionFactor")
      stars[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType4/",PT4_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.frame("ID" = stars$ParticleIDs,
                           "x"  = if(one_p_flag){stars$Coordinates[1]*.cm_to_kpc}else{stars$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                           "y"  = if(one_p_flag){stars$Coordinates[2]*.cm_to_kpc}else{stars$Coordinates[2,]*.cm_to_kpc},
                           "z"  = if(one_p_flag){stars$Coordinates[3]*.cm_to_kpc}else{stars$Coordinates[3,]*.cm_to_kpc},
                           "vx"  = if(one_p_flag){stars$Velocity[1]*.cms_to_kms}else{stars$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                           "vy"  = if(one_p_flag){stars$Velocity[2]*.cms_to_kms}else{stars$Velocity[2,]*.cms_to_kms},
                           "vz"  = if(one_p_flag){stars$Velocity[3]*.cms_to_kms}else{stars$Velocity[3,]*.cms_to_kms},
                           "Mass" = stars$Mass*.g_to_msol) # Mass in solar masses

    ssp = data.frame("Initial_Mass" = stars$InitialMass*.g_to_msol,
                     "Age" = as.numeric(.SFTtoAge(a = stars$StellarFormationTime, cores = cores)),
                     "Metallicity" = (colSums(stars$Metallicity[2:11,]))/(stars$Mass))

    remove(stars); remove(PT4_attr)

  } else {star_part=NULL; ssp=NULL}

  head$Type = "Magneticum"

  return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))

}

# Function for computing the stellar age from the formation time in parallel
.SFTtoAge = function(a, cores=1){
  cosdist = function(x) { return (celestial::cosdistTravelTime((1 / x) - 1)); }
  if (cores > 1) {
    doParallel::registerDoParallel(cores = cores)
    i = integer()
    output = foreach(i = 1:length(a), .packages = "celestial") %dopar% { cosdist(a[i]) }
    closeAllConnections()
  }
  else {
    output = lapply(a, cosdist)
  }
  return(output)
}

# Function to centre all galaxy particles based on stellar particle positions
.centre_galaxy = function(galaxy_data, centre=NA){
  if (!is.na(centre[1])){ # if an external centre is provided, use this to centre positions
    stellar_data = galaxy_data$star_part
    gas_data = galaxy_data$gas_part

    stellar_data$x = stellar_data$x - centre[1]
    stellar_data$y = stellar_data$y - centre[2]
    stellar_data$z = stellar_data$z - centre[3]
    star_r2 = stellar_data$x^2 + stellar_data$y^2 + stellar_data$z^2
    star_vcen = c(median(stellar_data$vx[star_r2 < 100]), # using the median velocities within
                  median(stellar_data$vy[star_r2 < 100]), #  10kpc of the galaxy centre to define
                  median(stellar_data$vz[star_r2 < 100])) #  the central velocity
    stellar_data$vx = stellar_data$vx - star_vcen[1]
    stellar_data$vy = stellar_data$vy - star_vcen[2]
    stellar_data$vz = stellar_data$vz - star_vcen[3]

    gas_data$x = gas_data$x - centre[1]
    gas_data$y = gas_data$y - centre[2]
    gas_data$z = gas_data$z - centre[3]
    gas_data$vx = gas_data$vx - star_vcen[1]
    gas_data$vy = gas_data$vy - star_vcen[2]
    gas_data$vz = gas_data$vz - star_vcen[3]

    galaxy_data$star_part = stellar_data
    galaxy_data$gas_part = gas_data

  } else if (!is.null(galaxy_data$star_part)){
    stellar_data = cen_galaxy(galaxy_data$star_part) # centering and computing medians for stellar particles
    galaxy_data$star_part = stellar_data$part_data
    if (!is.null(galaxy_data$gas_part)){ # if gas is present, centering these particles based on stellar medians
      gas_data = galaxy_data$gas_part
      gas_data$x = gas_data$x - stellar_data$xcen
      gas_data$y = gas_data$y - stellar_data$ycen
      gas_data$z = gas_data$z - stellar_data$zcen
      gas_data$vx = gas_data$vx - stellar_data$vxcen
      gas_data$vy = gas_data$vy - stellar_data$vycen
      gas_data$vz = gas_data$vz - stellar_data$vzcen
      galaxy_data$gas_part = gas_data
    }
  } else {
    gas_data = cen_galaxy(galaxy_data$gas_part)
    galaxy_data$gas_part$x = gas_data$part_data$x
    galaxy_data$gas_part$y = gas_data$part_data$y
    galaxy_data$gas_part$z = gas_data$part_data$z
    galaxy_data$gas_part$vx = gas_data$part_data$vx
    galaxy_data$gas_part$vy = gas_data$part_data$vy
    galaxy_data$gas_part$vz = gas_data$part_data$vz
  }
  return(galaxy_data)
}

# Functions for computing vector properties
.vector_mag = function(v){
  # Returns the magnitude of vector, v
  return(sqrt(sum(v^2)))
}

.vector_angle = function(v1,v2){
  # Returns the angle between vectors v1 and v2 in radians
  return(acos((v1%*%v2) / (.vector_mag(v1) * .vector_mag(v2))))
}

.vector_unit = function(v){
  # Returns a unit vector along the direction of vector v
  return(v/.vector_mag(v))
}

# Functions for rotating galaxies
.rot_mat_ang = function(v, angle){
  # Function for generating a rotation matrix that will rotate vector v by some angle
  x = v[1]; y = v[2]; z = v[3]
  co = cos(angle); si = sin(angle)

  R = rbind(c(co+(x*x*(1.-co)), (x*y*(1.-co))-(z*si), (x*z*(1.-co))+(y*si)),
            c((x*y*(1.-co))+(z*si), co+(y*y*(1.-co)), (y*z*(1.-co))-(x*si)),
            c((z*x*(1.-co))-(y*si), (z*y*(1.-co))+(x*si), co+(z*z*(1.-co))))

  return(R)
}

.rot_mat_vec = function(v1, v2){
  # Function for generation a rotation matrix that rotates vector v1 to match vector v2
  u1 = .vector_unit(v1); u2 = .vector_unit(v2)
  angle = -1 * .vector_angle(u1, u2)
  v = pracma::cross(u2, u1) / .vector_mag(pracma::cross(u2, u1))

  return(.rot_mat_ang(v, angle))
}

# Functions for measuring 3D shape
.new_half_mass_data = function(galaxy_data, p, q, half_mass){
  # function for getting all particles within the half mass radius (ordered by ellipsoid radii)
  x = galaxy_data$x; y = galaxy_data$y; z = galaxy_data$z
  if (is.na(half_mass)){
    half_mass = sum(galaxy_data$Mass) / 2
  }
  ellip_radius = sqrt((x*x) + ((y/p)*(y/p)) + ((z/q)*(z/q)))

  int_order = order(ellip_radius) # get the indicies of the radii in order (low to high)
  ordered_galaxy_data = galaxy_data[int_order,]
  cum_mass  = cumsum(ordered_galaxy_data$Mass) # cumulative sum of mass given this order
  half_mass_ind = which(abs(cum_mass - half_mass) == min(abs(cum_mass - half_mass)))[1] # at what radius does this half-mass now occur?

  return(ordered_galaxy_data[1:half_mass_ind,])
}

.ellipsoid_tensor = function(galaxy_data, p, q){
  # Computing the weighted ellipsoid tensor
  x = galaxy_data$x; y = galaxy_data$y; z = galaxy_data$z

  ellip_radius = sqrt((x*x) + ((y/p)*(y/p)) + ((z/q)*(z/q)))

  M = array(data = 0.0, dim = c(3,3))

  M[1,1] = sum((galaxy_data$Mass * x * x) / ellip_radius)
  M[1,2] = sum((galaxy_data$Mass * x * y) / ellip_radius)
  M[1,3] = sum((galaxy_data$Mass * x * z) / ellip_radius)
  M[2,1] = M[1,2]
  M[2,2] = sum((galaxy_data$Mass * y * y) / ellip_radius)
  M[2,3] = sum((galaxy_data$Mass * y * z) / ellip_radius)
  M[3,1] = M[1,3]
  M[3,2] = M[2,3]
  M[3,3] = sum((galaxy_data$Mass * z * z) / ellip_radius)

  return(M)
}

.ellipsoid_ratios_p_q = function(galaxy_data, p, q){
  # Function for calculating the p & q values from the ellipsoid tensor
  M = .ellipsoid_tensor(galaxy_data, p, q)
  eig = eigen(M)
  p = sqrt(eig$values[2]/eig$values[1])
  q = sqrt(eig$values[3]/eig$values[1])
  yax = eig$vectors[,2]
  zax = eig$vectors[,3]

  return(list("eigenvalues"= eig$values, "p" = p, "q" = q, "y_axis" = yax, "z_axis" = zax, "ellipsoid_tensor" = M))
}

# Function to iteratively find the shape and align at the half-mass stellar radius
.measure_pqj = function(galaxy_data, half_mass, abort_count=50){
  # Set up - we begin by assuming a sphere
  a = 1; b = 1; c = 1
  p = b/a; q = c/a
  aborted = 0; flag = 0
  cnt = 1
  temp_p = numeric(); temp_q = numeric()

  # Select all particles within initial half-mass (spherical) of stellar
  hm_galaxy_data = .new_half_mass_data(galaxy_data$star_part, p, q, half_mass)

  while (flag == 0){
    fit_ellip = .ellipsoid_ratios_p_q(hm_galaxy_data, p, q)
    temp_p[cnt] = fit_ellip$p # recording the axis ratios at this iteration
    temp_q[cnt] = fit_ellip$q

    # Check if current value is close to (or equal to) the last 10
    # iterations. If so return current p and q. The reason is that
    # sometimes the algorithm will end up jumping back and forth
    # between two similar values
    if (cnt > 10){
      last_10p = abs(temp_p[(cnt-9):cnt] - fit_ellip$p)
      last_10q = abs(temp_q[(cnt-9):cnt] - fit_ellip$q)
      if (all(last_10p < 0.01) & all(last_10q < 0.01)){
        flag = 1
      }
    }

    # Abort if iteration limit is reached, output current p and q.
    # The default abort count is 50, usually it doesn't take too long
    # to converge. Sometimes it just doesn't converge... rare, but I threw these out
    if (cnt > abort_count){
      aborted = 1
      flag = 1
    }

    # Check if current z-axis is in the same direction as
    # the unit vector (0,0,1). If not, rotate such that it is
    if (all.equal(.vector_unit(fit_ellip$z_axis), c(0,0,1)) != TRUE){
      Rz = .rot_mat_vec(fit_ellip$z_axis, c(0,0,1)) # Compute first the rotation to z
      v21 = as.numeric(Rz %*% fit_ellip$y_axis) # Then the next required rotation for new angle
      Ry = .rot_mat_vec(v21, c(0,1,0)) # to the y axis too

      new_star_coor_1 =  Rz %*% rbind(galaxy_data$star_part$x, galaxy_data$star_part$y, galaxy_data$star_part$z)
      new_star_vel_1  =  Rz %*% rbind(galaxy_data$star_part$vx, galaxy_data$star_part$vy, galaxy_data$star_part$vz)

      new_star_coor_2 = Ry %*% new_star_coor_1
      new_star_vel_2  = Ry %*% new_star_vel_1

      galaxy_data$star_part$x = new_star_coor_2[1,]; galaxy_data$star_part$y = new_star_coor_2[2,];
      galaxy_data$star_part$z = new_star_coor_2[3,]
      galaxy_data$star_part$vx = new_star_vel_2[1,]; galaxy_data$star_part$vy = new_star_vel_2[2,];
      galaxy_data$star_part$vz = new_star_vel_2[3,]

      if (!is.null(galaxy_data$gas_part)){ # if gas is present, aligning this based on the stellar coordinates
        new_gas_coor_1 =  Rz %*% rbind(galaxy_data$gas_part$x, galaxy_data$gas_part$y, galaxy_data$gas_part$z)
        new_gas_vel_1  =  Rz %*% rbind(galaxy_data$gas_part$vx, galaxy_data$gas_part$vy, galaxy_data$gas_part$vz)
        new_gas_coor_2 = Ry %*% new_gas_coor_1
        new_gas_vel_2  = Ry %*% new_gas_vel_1
        galaxy_data$gas_part$x = new_gas_coor_2[1,]; galaxy_data$gas_part$y = new_gas_coor_2[2,];
        galaxy_data$gas_part$z = new_gas_coor_2[3,]
        galaxy_data$gas_part$vx = new_gas_vel_2[1,]; galaxy_data$gas_part$vy = new_gas_vel_2[2,];
        galaxy_data$gas_part$vz = new_gas_vel_2[3,]
      }

    }

    hm_galaxy_data = .new_half_mass_data(galaxy_data$star_part, fit_ellip$p, fit_ellip$q, half_mass)
    p = fit_ellip$p
    q = fit_ellip$q
    cnt = cnt + 1

  }

  return(list("galaxy_data" = galaxy_data, "p" = mean(tail(temp_p, n=6)), "q" = mean(tail(temp_q, n=6))))

}

# Function to align full galaxy based on the stellar particles
.align_galaxy = function(galaxy_data, half_mass=NA){
  if (is.null(galaxy_data$star_part)){ # if there are no stellar particles (just gas), use these
    dummy_data = list(star_part = galaxy_data$gas_part,
                      gas_part= galaxy_data$star_part,
                      head = galaxy_data$head,
                      ssp = galaxy_data$ssp)
    dummy = .measure_pqj(dummy_data, half_mass)
    data = list(galaxy_data = vector("list"))
    data$galaxy_data = list(star_part = dummy$galaxy_data$gas_part,
                            gas_part  = dummy$galaxy_data$star_part,
                            head      = dummy$galaxy_data$head,
                            ssp       = dummy$galaxy_data$ssp)
  } else {
    data = .measure_pqj(galaxy_data, half_mass)
  }
  return(data$galaxy_data)
}

# Functions for smoothing SPH kernels
.wendland_c2 = function(r){ # SPH smoothing kernel used in EAGLE
  return((21/(2*pi))*((1-r)^4)*((4*r) + 1))
} # input a radial position, r
# returns the corresponding kernel weight at that radius

.wendland_c6 = function(r){ # SPH smoothing kernel used in Magneticum
  return((1365/(64*pi))*((1-r)^8)*(1 + (8*r) + (25*r^2) + (32*r^3)))
} # input a radial position, r
  # returns the corresponding kernel weight at that radius

.generate_uniform_sphere = function(number_of_points, kernel="WC2"){

  # Function for generating random coordinates that
  # uniformly sample the volume of a sphere and computing their corresponding
  # weights

  # input the number of new particles you would like to spawn ("number_of_points")
  # returns a data.frame containing the longitude and latitude of those new
  # points, the radial coordinate "r" as a function of the softening length "h"
  # and the corresponding SPH kernel weight normalised so that the total of the
  # weights sums to 1 (i.e. forcing mass conservation).

  xyz = cbind(stats::rnorm(number_of_points), stats::rnorm(number_of_points), stats::rnorm(number_of_points))
  r = stats::runif(number_of_points, min = 0, max = 1)^(1/3)
  den = sqrt((xyz[,1]^2) + (xyz[,2]^2) + (xyz[,3]^2))
  xyz_norm = (r*xyz)/den
  # method for calculating a randomly distibution of n points uniformly
  # across a spherical volume

  sph = sphereplot::car2sph(xyz_norm) # convert to spherical coordinates
  if (kernel == "WC2"){
    weights = .wendland_c2(sph[,3])
    sph_kernel = data.frame("long" = sph[,1], "lat" = sph[,2],
                            "r/h" = sph[,3], "weight" = weights/sum(weights))
  }
  if (kernel == "WC6"){
    weights = .wendland_c6(sph[,3])
    sph_kernel = data.frame("long" = sph[,1], "lat" = sph[,2],
                            "r/h" = sph[,3], "weight" = weights/sum(weights))

  }

  return(sph_kernel)
}

# Function to generate spectra (w/o mass weighting)
.spectra = function(Metallicity, Age, Template, cores){
  f = function(metallicity, age) {
    Z = as.numeric(ProSpect::interp_quick(metallicity, Template$Z, log = TRUE))
    A = as.numeric(ProSpect::interp_quick(age * 1e9, Template$Age, log = TRUE))

    weights = data.frame("hihi" = Z[4] * A[4],   # ID_lo = 1, ID_hi = 2, wt_lo = 3, wt_hi = 4
                         "hilo" = Z[4] * A[3],
                         "lohi" = Z[3] * A[4],
                         "lolo" = Z[3] * A[3])

    part_spec = array(data = 0.0, dim = c(1, length(Template$Wave)))

    part_spec = ((Template$Zspec[[Z[2]]][A[2],] * weights$hihi) +
                 (Template$Zspec[[Z[2]]][A[1],] * weights$hilo) +
                 (Template$Zspec[[Z[1]]][A[2],] * weights$lohi) +
                 (Template$Zspec[[Z[1]]][A[1],] * weights$lolo))

    return(part_spec)
  }

  if (length(Age) == 1){
    spectra = list(f(Metallicity, Age))
    return(spectra)
  } else {

    if (cores > 1) {
      doParallel::registerDoParallel(cores = cores)
      i = integer()
      part_spec = foreach::foreach(i = 1:length(Metallicity), .packages = c("ProSpect", "SimSpin", "foreach")) %dopar% {
        f(Metallicity[i], Age[i])
      }
      closeAllConnections()
      output = part_spec
    } else {
      part_spec = mapply(f, Metallicity, Age)
      output = lapply(seq_len(ncol(part_spec)), function(i) part_spec[,i])
    }
    return(output)

  }
}

.circular_ap=function(sbin){
  ap_region = matrix(data = 0, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  x = matrix(data = rep(seq(1,sbin), each=sbin), nrow = sbin, ncol = sbin)
  y = matrix(data = rep(seq(sbin,1), sbin), nrow = sbin, ncol = sbin)
  xx = x - xcentre; yy = y - ycentre
  rr = sqrt(xx^2 + yy^2)
  ap_region[rr<= sbin/2] = 1
  return(as.vector(ap_region))
}

.hexagonal_ap=function(sbin){
  ap_region = matrix(data = 0, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  for (x in 1:sbin){
    for (y in 1:sbin){
      xx = x - xcentre
      yy = y - ycentre
      rr = (2 * (sbin / 4) * (sbin * sqrt(3) / 4)) - ((sbin / 4) ) * abs(yy) - ((sbin * sqrt(3) / 4)) * abs(xx)
      if ((rr >= 0) && (abs(xx) < sbin/2) && (abs(yy) < (sbin  * sqrt(3) / 4))){
        ap_region[x,y] = 1
      }
    }
  }
  return(as.vector(ap_region))
}

.interpolate_spectra = function(shifted_wave, spectra, wave_seq){ # function for interpolating the spectra onto a new grid
  shifted_spectra = vector(mode = "list", length = dim(shifted_wave)[1])
  for(j in 1:dim(shifted_wave)[1]){
    shifted_spectra[[j]] = stats::approx(x = shifted_wave[j,], y = spectra[j,], xout = wave_seq, rule=1)[2]
  }
  output = matrix(unlist(shifted_spectra, use.names=FALSE), nrow = dim(shifted_wave)[1], byrow = TRUE)
  spaxel_spectra = colSums(output)

  return(spaxel_spectra)
}

.sum_velocities = function(galaxy_sample, observation){
  vel_diff = function(lum, vy){diff((lum * pnorm(observation$vbin_edges, mean = vy,
                                                     sd = observation$vbin_error)))}

  bins = mapply(vel_diff, galaxy_sample$luminosity, galaxy_sample$vy)

  return(rowSums(bins))

}

.sum_gas_velocities = function(galaxy_sample, observation){
  vel_diff = function(mass, vy){diff((mass * pnorm(observation$vbin_edges, mean = vy,
                                                     sd = observation$vbin_error)))}

  bins = mapply(vel_diff, galaxy_sample$Mass, galaxy_sample$vy)

  return(rowSums(bins))

}

# Function to apply LSF to spectra
.lsf_convolution = function(observation, luminosity, lsf_sigma){

  scaled_sigma = lsf_sigma / observation$wave_res # scaling the size of the gaussian to match the pixel dimensions
  kernel = stats::dnorm(seq(-scaled_sigma*5,scaled_sigma*5,length.out = 25), mean = 0, sd = scaled_sigma)
  lsf_gauss = kernel/sum(kernel)
  lum = stats::convolve(luminosity, lsf_gauss, type="open")
  end = (length(luminosity) + length(lsf_gauss) - 1) - 12

  return(lum[13:end])
}

# Function to add noise
.add_noise = function(luminosity, S2N){
  noise_level = min(luminosity) / S2N
  noise = stats::rpois(length(luminosity), lambda = noise_level)
  noisey_lum = luminosity + (stats::rnorm(length(luminosity), mean = 0, sd=1) * noise)
}

# Function to generate a Gaussian kernel
.gaussian_kernel = function(m, n, sigma){
  dim = pracma::meshgrid(-((m-1)/2):((m-1)/2), -((n-1)/2):((n-1)/2))
  hg = exp(-(dim$X^2 + dim$Y^2)/(2*sigma^2))
  kernel = hg / sum(hg)
  return(kernel)
}

# Functions for computing spaxel properties
# spectral mode -
.spectral_spaxels = function(part_in_spaxel, wavelength, observation, galaxy_data, simspin_data, verbose){

  spectra = matrix(data = 0.0, ncol = observation$wave_bin, nrow = observation$sbin^2)
  vel_los = array(data = 0.0, dim = observation$sbin^2)
  dis_los = array(data = 0.0, dim = observation$sbin^2)
  lum_map = array(data = 0.0, dim = observation$sbin^2)
  part_map = array(data=0, dim = observation$sbin^2)

  for (i in 1:(dim(part_in_spaxel)[1])){ # computing the spectra at each occupied spatial pixel position

    num_part = length(part_in_spaxel$val[[i]]) # number of particles in spaxel
    part_map[part_in_spaxel$spaxel_ID[i]] = num_part

    # if number is greater than the particle limit
    #if (num_part >= observation$particle_limit){
      galaxy_sample = galaxy_data[part_in_spaxel$val[[i]],]
      intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]), nrow = num_part, byrow = T) *
        (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra

      # pulling wavelengths and using doppler formula to compute the shift in
      #   wavelengths caused by LOS velocity
      wave = matrix(data = rep(wavelength, num_part), nrow = num_part, byrow=T)
      wave_shift = ((galaxy_sample$vy / .speed_of_light) * wave) + wave

      # interpolate each shifted wavelength to telescope grid of wavelengths
      #   and sum to one spectra
      luminosity = .interpolate_spectra(shifted_wave = wave_shift, spectra = intrinsic_spectra,
                                        wave_seq = observation$wave_seq)

      if (observation$LSF_conv){ # should the spectra be degraded for telescope LSF?
        luminosity = .lsf_convolution(observation=observation, luminosity=luminosity, lsf_sigma=observation$lsf_sigma)
      }

      if (!is.na(observation$signal_to_noise)){ # should we add noise?
        luminosity = .add_noise(luminosity, observation$signal_to_noise)
      }

      # transform luminosity into flux detected at telescope
      #    flux in units erg/s/cm^2/Ang
      spectra[part_in_spaxel$spaxel_ID[i],] = (luminosity*.lsol_to_erg) /
                                                  (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) / (1 + observation$z)
      lum_map[part_in_spaxel$spaxel_ID[i]] = ProSpect::bandpass(wave = observation$wave_seq,
                                                                flux = spectra[part_in_spaxel$spaxel_ID[i],],
                                                                filter = observation$filter, flux_in = "wave",
                                                                flux_out = "wave")
      vel_los[part_in_spaxel$spaxel_ID[i]] = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
      dis_los[part_in_spaxel$spaxel_ID[i]] = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))

    #} else { # if insufficient particles in spaxel
    #  spectra[part_in_spaxel$spaxel_ID[i],] = rep(NA, observation$wave_bin)
    #  lum_map[part_in_spaxel$spaxel_ID[i]] = NA
    #  vel_los[part_in_spaxel$spaxel_ID[i]] = NA
    #  dis_los[part_in_spaxel$spaxel_ID[i]] = NA
    #}
    if (verbose){cat(i, "... ", sep = "")}
  }
  return(list(spectra, lum_map, vel_los, dis_los, part_map))
}

.spectral_spaxels_mc = function(part_in_spaxel, wavelength, observation, galaxy_data, simspin_data, verbose, cores){

  spectra = matrix(data = 0.0, ncol = observation$wave_bin, nrow = observation$sbin^2)
  vel_los = array(data = 0.0, dim = observation$sbin^2)
  dis_los = array(data = 0.0, dim = observation$sbin^2)
  lum_map = array(data = 0.0, dim = observation$sbin^2)
  part_map = array(data=0, dim = observation$sbin^2)

  doParallel::registerDoParallel(cores)

  i = integer()
  output = foreach(i = 1:(dim(part_in_spaxel)[1]), .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list())) %dopar% {

                     num_part = length(part_in_spaxel$val[[i]])
                     part_map = num_part
                     # if the number of particles in the spaxel is greater than the particle limit
                     #if (num_part >= observation$particle_limit){
                       galaxy_sample = galaxy_data[part_in_spaxel$val[[i]],]
                       intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]),
                                                  nrow = num_part, byrow = T) *
                                                  (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra

                       # pulling wavelengths and using doppler formula to compute the shift in
                       #   wavelengths caused by LOS velocity
                       wave = matrix(data = rep(wavelength, num_part), nrow = num_part, byrow=T)
                       wave_shift = ((galaxy_sample$vy / .speed_of_light) * wave) + wave

                       # interpolate each shifted wavelength to telescope grid of wavelengths
                       #   and sum to one spectra
                       luminosity = .interpolate_spectra(shifted_wave = wave_shift,
                                                         spectra = intrinsic_spectra,
                                                         wave_seq = observation$wave_seq)

                       if (observation$LSF_conv){ # should the spectra be degraded for telescope LSF?
                         luminosity = .lsf_convolution(observation=observation, luminosity=luminosity,
                                                       lsf_sigma=observation$lsf_sigma)
                       }
                       if (!is.na(observation$signal_to_noise)){ # should we add noise?
                         luminosity = .add_noise(luminosity, observation$signal_to_noise)
                       }

                       # transform luminosity into flux detected at telescope
                       #    flux in units erg/s/cm^2/Ang
                       spectra = (luminosity*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) / (1 + observation$z)
                       lum_map = ProSpect::bandpass(wave = observation$wave_seq,
                                                    flux = spectra,
                                                    filter = observation$filter,
                                                    flux_in = "wave", flux_out = "wave")
                       vel_los = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
                       dis_los = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
                     #} else { # if insufficient particles in spaxel
                      # spectra = rep(NA, observation$wave_bin)
                      # lum_map = NA
                      # vel_los = NA
                      # dis_los = NA
                     #}
                     result = list(spectra, lum_map, vel_los, dis_los, part_map)
                     return(result)
                     closeAllConnections()
                   }

  spectra[part_in_spaxel$spaxel_ID,] = matrix(unlist(output[[1]]), ncol=observation$wave_bin, byrow = T)
  lum_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[2]]))
  vel_los[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[3]]))
  dis_los[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[4]]))
  part_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[5]]))

  return(list(spectra, lum_map, vel_los, dis_los, part_map))
}


# stellar velocity mode -
.velocity_spaxels = function(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, mass_flag){

  vel_spec = matrix(data = 0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  lum_map  = array(data = 0.0, dim = observation$sbin^2)
  age_map  = array(data = 0.0, dim = observation$sbin^2)
  met_map  = array(data = 0.0, dim = observation$sbin^2)
  part_map = array(data=0, dim = observation$sbin^2)

  for (i in 1:(dim(part_in_spaxel)[1])){

    num_part = length(part_in_spaxel$val[[i]]) # number of particles in spaxel
    part_map[part_in_spaxel$spaxel_ID[i]] = num_part

    # if number is greater than the particle limit
    #if (num_part >= observation$particle_limit){
      galaxy_sample = galaxy_data[part_in_spaxel$val[[i]],]

      if (mass_flag){

        galaxy_sample$luminosity = galaxy_sample$Mass # replace luminosity weighting with a mass weight

      } else {
        intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]), nrow = num_part, byrow = T) *
          (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra

        # transform luminosity into flux detected at telescope
        #    flux in units erg/s/cm^2/Ang
        spectral_flux = (intrinsic_spectra*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) / (1 + observation$z)

        # computing the r-band luminosity per particle from spectra
        galaxy_sample$luminosity = apply(spectral_flux, 1, ProSpect::bandpass,
                                         wave = simspin_data$wave,
                                         filter = observation$filter, flux_in = "wave", flux_out = "wave")
      }


      # adding the "gaussians" of each particle to the velocity bins
      vel_spec[part_in_spaxel$spaxel_ID[i],] = .sum_velocities(galaxy_sample = galaxy_sample, observation = observation)
      lum_map[part_in_spaxel$spaxel_ID[i]] = sum(galaxy_sample$luminosity)
      vel_los[part_in_spaxel$spaxel_ID[i]] = .meanwt(galaxy_sample$vy, galaxy_sample$luminosity)
      dis_los[part_in_spaxel$spaxel_ID[i]] = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$luminosity))
      age_map[part_in_spaxel$spaxel_ID[i]] = .meanwt(galaxy_sample$Age, galaxy_sample$Initial_Mass)
      met_map[part_in_spaxel$spaxel_ID[i]] = .meanwt(galaxy_sample$Metallicity, galaxy_sample$Initial_Mass)

    #} else { # if insufficient particles in spaxel
    #  vel_spec[part_in_spaxel$spaxel_ID[i],] = rep(NA, observation$vbin)
    #  lum_map[part_in_spaxel$spaxel_ID[i]] = NA
    #  vel_los[part_in_spaxel$spaxel_ID[i]] = NA
    #  dis_los[part_in_spaxel$spaxel_ID[i]] = NA
    #  age_map[part_in_spaxel$spaxel_ID[i]] = NA
    #  met_map[part_in_spaxel$spaxel_ID[i]] = NA

    #}
    if (verbose){cat(i, "... ", sep = "")}

  }

  return(list(vel_spec, lum_map, vel_los, dis_los, age_map, met_map, part_map))
}

.velocity_spaxels_mc = function(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, cores, mass_flag){

  vel_spec = matrix(data = 0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  lum_map  = array(data = 0.0, dim = observation$sbin^2)
  age_map  = array(data = 0.0, dim = observation$sbin^2)
  met_map  = array(data = 0.0, dim = observation$sbin^2)
  part_map = array(data=0, dim = observation$sbin^2)

  doParallel::registerDoParallel(cores)

  i = integer()
  output = foreach(i = 1:(dim(part_in_spaxel)[1]), .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list(), list(), list())) %dopar% {

                     num_part = length(part_in_spaxel$val[[i]])
                     part_map = num_part
                     # if the number of particles in the spaxel is greater than the particle limit
                     #if (num_part >= observation$particle_limit){
                       galaxy_sample = galaxy_data[part_in_spaxel$val[[i]],]

                       if (mass_flag){

                         galaxy_sample$luminosity = galaxy_sample$Mass

                       } else {
                         intrinsic_spectra = matrix(unlist(simspin_data$spectra[galaxy_sample$sed_id]),
                                                  nrow = num_part, byrow = T) *
                         (galaxy_sample$Initial_Mass * 1e10) # reading relavent spectra
                         # transform luminosity into flux detected at telescope
                         #    flux in units erg/s/cm^2/Ang
                         spectral_flux = (intrinsic_spectra*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) / (1 + observation$z)

                         # computing the r-band luminosity per particle from spectra
                         galaxy_sample$luminosity = apply(spectral_flux, 1, ProSpect::bandpass,
                                                          wave = simspin_data$wave,
                                                          filter = observation$filter, flux_in = "wave", flux_out = "wave")
                         }

                       # adding the "gaussians" of each particle to the velocity bins
                       vel_spec = .sum_velocities(galaxy_sample = galaxy_sample, observation = observation)
                       lum_map = sum(galaxy_sample$luminosity)
                       vel_los = .meanwt(galaxy_sample$vy, galaxy_sample$luminosity)
                       dis_los = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$luminosity))
                       age_map = .meanwt(galaxy_sample$Age, galaxy_sample$Initial_Mass)
                       met_map = .meanwt(galaxy_sample$Metallicity, galaxy_sample$Initial_Mass)

                    # } else { # if insufficient particles in spaxel
                      # vel_spec = rep(NA, observation$vbin)
                      # lum_map = NA
                      # vel_los = NA
                      # dis_los = NA
                      # age_map = NA
                      # met_map = NA
                     #}
                     result = list(vel_spec, lum_map, vel_los, dis_los, age_map, met_map, part_map)
                     return(result)
                     closeAllConnections()
                   }

  vel_spec[part_in_spaxel$spaxel_ID,] = matrix(unlist(output[[1]]), ncol=observation$vbin, byrow = T)
  lum_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[2]]))
  vel_los[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[3]]))
  dis_los[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[4]]))
  age_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[5]]))
  met_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[6]]))
  part_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[7]]))

  return(list(vel_spec, lum_map, vel_los, dis_los, age_map, met_map, part_map))

}

# gas velocity mode -
.gas_velocity_spaxels = function(part_in_spaxel, observation, galaxy_data, simspin_data, verbose){

  vel_spec = matrix(data = 0.0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  mass_map = array(data = 0.0, dim = observation$sbin^2)
  SFR_map  = array(data = 0.0, dim = observation$sbin^2)
  Z_map    = array(data = 0.0, dim = observation$sbin^2)
  OH_map   = array(data = 0.0, dim = observation$sbin^2)

  for (i in 1:(dim(part_in_spaxel)[1])){

    num_part = length(part_in_spaxel$val[[i]]) # number of particles in spaxel

    # if number is greater than the particle limit
    #if (num_part >= observation$particle_limit){
      galaxy_sample = galaxy_data[part_in_spaxel$val[[i]],]

      # adding the "gaussians" of each particle to the velocity bins
      vel_spec[part_in_spaxel$spaxel_ID[i],] = .sum_gas_velocities(galaxy_sample = galaxy_sample, observation = observation)
      mass_map[part_in_spaxel$spaxel_ID[i]]  = sum(galaxy_sample$Mass)
      vel_los[part_in_spaxel$spaxel_ID[i]]   = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
      dis_los[part_in_spaxel$spaxel_ID[i]]   = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
      SFR_map[part_in_spaxel$spaxel_ID[i]]   = .meanwt(galaxy_sample$SFR, galaxy_sample$Mass)
      Z_map[part_in_spaxel$spaxel_ID[i]]     = log10(mean(galaxy_sample$Metallicity)/0.0127)
      OH_map[part_in_spaxel$spaxel_ID[i]]    = log10(mean(galaxy_sample$Oxygen/galaxy_sample$Hydrogen))+12

    #} else { # if insufficient particles in spaxel
    #  vel_spec[part_in_spaxel$spaxel_ID[i],] = rep(NA, observation$vbin)
    #  mass_map[part_in_spaxel$spaxel_ID[i]]  = NA
    #  vel_los[part_in_spaxel$spaxel_ID[i]]   = NA
    #  dis_los[part_in_spaxel$spaxel_ID[i]]   = NA
    #  SFR_map[part_in_spaxel$spaxel_ID[i]]   = NA
    #  Z_map[part_in_spaxel$spaxel_ID[i]]     = NA
    #  OH_map[part_in_spaxel$spaxel_ID[i]]    = NA
    #}
    if (verbose){cat(i, "... ", sep = "")}

  }

  return(list(vel_spec, mass_map, vel_los, dis_los, SFR_map, Z_map, OH_map))
}

.gas_velocity_spaxels_mc = function(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, cores){

  vel_spec = matrix(data = 0.0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  mass_map  = array(data = 0.0, dim = observation$sbin^2)
  SFR_map  = array(data = 0.0, dim = observation$sbin^2)
  Z_map    = array(data = 0.0, dim = observation$sbin^2)
  OH_map   = array(data = 0.0, dim = observation$sbin^2)

  doParallel::registerDoParallel(cores)

  i = integer()
  output = foreach(i = 1:(dim(part_in_spaxel)[1]), .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list(), list(), list())) %dopar% {

                     num_part = length(part_in_spaxel$val[[i]])
                     # if the number of particles in the spaxel is greater than the particle limit
                    # if (num_part >= observation$particle_limit){
                       galaxy_sample = galaxy_data[part_in_spaxel$val[[i]],]

                       # adding the "gaussians" of each particle to the velocity bins
                       vel_spec = .sum_gas_velocities(galaxy_sample = galaxy_sample, observation = observation)
                       mass_map = sum(galaxy_sample$Mass)
                       vel_los  = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
                       dis_los  = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
                       SFR_map  = .meanwt(galaxy_sample$SFR, galaxy_sample$Mass)
                       Z_map    = log10(mean(galaxy_sample$Metallicity)/0.0127)
                       OH_map   = log10(mean(galaxy_sample$Oxygen/galaxy_sample$Hydrogen))+12


                     #} else { # if insufficient particles in spaxel
                    #   vel_spec = rep(NA, observation$vbin)
                    #   mass_map = NA
                    #   vel_los  = NA
                    #   dis_los  = NA
                    #   SFR_map  = NA
                    #   Z_map    = NA
                    #   OH_map   = NA
                    # }
                     result = list(vel_spec, mass_map, vel_los, dis_los, SFR_map, Z_map, OH_map)
                     return(result)
                     closeAllConnections()
                   }

  vel_spec[part_in_spaxel$spaxel_ID,] = matrix(unlist(output[[1]]), ncol=observation$vbin, byrow = T)
  mass_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[2]]))
  vel_los[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[3]]))
  dis_los[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[4]]))
  SFR_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[5]]))
  Z_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[6]]))
  OH_map[part_in_spaxel$spaxel_ID] = matrix(unlist(output[[7]]))

  return(list(vel_spec, mass_map, vel_los, dis_los, SFR_map, Z_map, OH_map))

}

# spawn gas particles -
.sph_spawn = function(gas_part, new_gas_part, sph_spawn_n, kernel){
# Function for generating "n" number of gas particles (specified by
# "sph_spawn_n") smoothed across some volume by the relevent
# SPH kernel. This function feeds into the process at the
# "make_simspin_file()" stage.

  no_gas = length(gas_part$ID)

  for (each in 1:no_gas){ # for each particle
    ind1 = ((each*sph_spawn_n)-sph_spawn_n)+1; ind2 = (each*sph_spawn_n)
    part = gas_part[each,]
    # pull the data relevent to that particle from the original data.frame

    rand_pos = .generate_uniform_sphere(sph_spawn_n, kernel = kernel)
    # distribute that particle randomly across the SPH kernel volume
    # as a function of smoothing length
    rand_pos$r = rand_pos$r.h * part$SmoothingLength
    # use the particle's specific smoothing length to scale the radial
    # positions of the particle.
    new_xyz = sphereplot::sph2car(cbind(rand_pos$long, rand_pos$lat, rand_pos$r))
    # convert spherical coordinates back into cartesian coords

    new_gas_part$ID[ind1:ind2] = as.numeric(paste0(part$ID, 1:sph_spawn_n))
    new_gas_part$x[ind1:ind2] = part$x+new_xyz[,1]
    new_gas_part$y[ind1:ind2] = part$y+new_xyz[,2]
    new_gas_part$z[ind1:ind2] = part$z+new_xyz[,3]
    new_gas_part$Mass[ind1:ind2] = part$Mass*rand_pos$weight
    # in the new data.frame of particle properties, assign their
    # new positions and masses scaled by the kernel weight.
  }

  return(new_gas_part)
}

.sph_spawn_mc = function(gas_part, new_gas_part, sph_spawn_n, kernel, cores){

  # Parallel version of the function ".sph_spawn" above.
  doParallel::registerDoParallel(cores)
  no_gas = length(gas_part$ID)
  ID = numeric(no_gas); x = numeric(no_gas); y = numeric(no_gas)
  z = numeric(no_gas); Mass = numeric(no_gas) # initialising

  i = integer()
  output = foreach(i = 1:no_gas, .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list())) %dopar% {

                     part = gas_part[i,]

                     rand_pos = .generate_uniform_sphere(sph_spawn_n, kernel = kernel)
                     rand_pos$r = rand_pos$r.h * part$SmoothingLength
                     new_xyz = sphereplot::sph2car(cbind(rand_pos$long, rand_pos$lat, rand_pos$r))

                     ID = as.numeric(paste0(part$ID, 1:sph_spawn_n))
                     x = part$x+new_xyz[,1]
                     y = part$y+new_xyz[,2]
                     z = part$z+new_xyz[,3]
                     Mass = part$Mass*rand_pos$weight

                     return(list(ID, x, y, z, Mass))
                     closeAllConnections()
                   }

  new_gas_part$ID = as.numeric(unlist(output[[1]]))
  new_gas_part$x = as.numeric(unlist(output[[2]]))
  new_gas_part$y = as.numeric(unlist(output[[3]]))
  new_gas_part$z = as.numeric(unlist(output[[4]]))
  new_gas_part$Mass = as.numeric(unlist(output[[5]]))

  return(new_gas_part)
}
