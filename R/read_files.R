# Author: Kate Harborne
# Co-author: Alice Serene
# Date: 20/02/2025
# Title: read_files.R
# Description: Hidden from user functions for reading various input simulation
#              files, including
#                 - Gadget binaries
#                 - Tipsy binaries
#                 - Gadget HDF5
#                 - EAGLE
#                 - IllustrisTNG
#                 - Magneticum
#                 - HorizonAGN
#                 - Colibre

.get_file_type = function(f){

  # Input files could be in binary or HDF5 format. Investigate format before
  # pointing to the correct read function.

  data = file(f, "rb") # open file for reading in binary mode

  block = readBin(data, "integer", n=1)

  close(data)

  if (block == 256){

    return(output = "gadget_binary")

  } else {

    return(output = "hdf5")

  }

}


# Function for reading in Gadget binary files
.read_gadget = function(f){
  data = file(f, "rb") # open file for reading in binary mode

  block         = readBin(data, "integer", n=1) #block size field, giving the length of the header
  if(block!=256){close(data); stop("Not a Gadget binary file.")}
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

  part = data.table::data.table("ID" = id,          # the particle data table
                                "x" = pos[extract], "y" = pos[extract+1],
                                "z" = pos[extract+2],
                                "vx" = vel[extract], "vy" = vel[extract+1],
                                "vz"=vel[extract+2],
                                "Mass" = masses*1e10) # masses in Msol

  head = list("Npart" = c(Npart[1], 0, Npart[3], Npart[4], Npart[5], 0), # number of gas and stars
              "Time" = Time, "Redshift" = Redshift, # relevent simulation data
              "Nall" = Nall, "Type"="nbody", # number of particles in the original file
              "H0" = HubbleParam*100, "OmegaM" = Omega0, "OmegaL" = OmegaLambda,
              "OmegaR" = NA)

  Npart_sum = cumsum(Npart) # cumulative number of each particle type

  star_part = part[(Npart_sum[2]+1):Npart_sum[5],]

  if (Npart[1] > 0){
    gas_part = part[1:Npart_sum[1],]
  } else {gas_part = NULL}

  return(list(star_part = star_part, gas_part = gas_part, head = head))

}

# Functions for reading in HDF5 files
.read_hdf5 = function(f, cores=1){

  data = hdf5r::h5file(f, mode="r")

  # Read in all attributes listed in the header
  header_attr = hdf5r::list.attributes(data[["Header"]])

  # Create a list to store each variable
  head = vector("list", length(header_attr))
  names(head) = header_attr

  # Read in each variable and store in list
  for (i in 1:length(header_attr)){
    head[[i]] = hdf5r::h5attr(data[["Header"]], paste0(header_attr[i]))
  }

  # check for missing header fields
  required_headers = c("BoxSize", "Redshift", "HubbleParam", "MassTable", "Omega0", "OmegaLambda")

  if (!all(required_headers %in% names(head))){
    stop("Error. Missing a required header field. \n
         One of `BoxSize`, `Redshift`, `HubbleParam`, `MassTable`, `Omega0` or `OmegaLambda` is missing. \n
         See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#header for more details.")
  }

  other_headers = c("NumPart_ThisFile", "NumPart_Total")
  if (!any(other_headers %in% names(head))){
    stop("Error. Missing a required header field. \n
         One of `NumPart_ThisFile` or `NumPart_Total` is missing. \n
         See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#header for more details.")
  }


  # default (if header if blank) is a gadget file.
  if(is.null(head$RunLabel) && is.null(head$SimulationName)){
    gadget2 = T
    eagle = F
    magneticum = F
    horizonagn = F
    illustristng = F
    colibre = F
  } else {

    if ("SimulationName" %in% names(head)){
      gadget2 = F
      eagle = F
      magneticum = F
      horizonagn = F
      colibre = F
      if(stringr::str_detect(stringr::str_to_lower(head$SimulationName), "tng")){illustristng = T}else{illustristng = F}
    } else {
      gadget2 = F
      if(stringr::str_detect(stringr::str_to_lower(head$RunLabel), "tng")){illustristng = T}else{illustristng = F}
      if(stringr::str_detect(stringr::str_to_lower(head$RunLabel), "eagle")){eagle = T}else{eagle=F} # determining if EAGLE input (based on number of parameters in Header)
      if(stringr::str_detect(stringr::str_to_lower(head$RunLabel), "magneticum")){magneticum = T}else{magneticum=F}
      if(stringr::str_detect(stringr::str_to_lower(head$RunLabel), "horizon")){horizonagn = T}else{horizonagn = F}
      if(stringr::str_detect(stringr::str_to_lower(head$RunLabel), "colibre")){colibre = T}else{colibre = F}
    }

  }

  # Read particle data differently depending on the simulation being read in...
  if (gadget2){output = .gadget2_read_hdf5(data, head)}
  if (eagle){output = .eagle_read_hdf5(data, head, cores)}
  if (magneticum){output = .magneticum_read_hdf5(data, head, cores)}
  if (horizonagn){output = .horizonagn_read_hdf5(data, head, cores)}
  if (illustristng){output = .illustristng_read_hdf5(data, head, cores)}
  if (colibre){output = .colibre_read_hdf5(data, head, cores)}

  hdf5r::h5close(data)

  return(output)
}

.gadget2_read_hdf5 = function(data, head){

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if (!("PartType2" %in% groups) &
      !("PartType3" %in% groups) &
      ("PartType4" %in% groups)){
    stop("Error. SimSpin is trying to process the input simulation as an N-body file. \n
          No stars are present in PartType2 or PartType3, but stars are present in PartType4. These stars will be missed from the output. \n
          Is this meant to be a Hydrodyanmical model? See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#header `RunLabel` for more info.")
  }

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

    gas_part = data.table::data.table("ID" = gas$ParticleIDs,
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

    disk_part = data.table::data.table("ID" = stars$ParticleIDs,
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

    star_part = data.table::data.table("ID" = c(disk_part$ID, stars$ParticleIDs),
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

    star_part = data.table::data.table("ID" = stars$ParticleIDs,
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

    star_part = data.table::data.table("ID" = stars$ParticleIDs,
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
              "Time" = head$Time, "Redshift" = head$Redshift, # relevant simulation data
              "Nall" = head$NumPart_Total, "Type"="nbody") # number of particles in the original file

  return(list(star_part=star_part, gas_part=gas_part, head=head))

}

.eagle_read_hdf5 = function(data, head, cores){

  head$Type = "EAGLE"
  head$H0 = head$HubbleParam * 100
  names(head)[names(head) == "Omega0"] = "OmegaM"
  names(head)[names(head) == "OmegaLambda"] = "OmegaL"
  head$OmegaR = 0

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if ("PartType0" %in% groups){ # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])

    expected_names_gas = c("Coordinates", "Density", "Mass", "ParticleIDs",
                           "ElementAbundance/Oxygen", "SmoothedElementAbundance/Oxygen",
                           "ElementAbundance/Hydrogen", "SmoothedElementAbundance/Hydrogen",
                           "SmoothedMetallicity", "Metallicity",
                           "StarFormationRate", "Velocity", "SmoothingLength",
                           "Temperature", "InternalEnergy")
    PT0_attr = PT0_attr[which(PT0_attr %in% expected_names_gas)] # trim list to only read in necessary data sets

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

    gas = .check_names(gas, type="EAGLE")
    eagle_gas_names = c("SmoothingLength", "Temperature", "InternalEnergy")
    if (!all(eagle_gas_names %in% names(gas))){
      stop("Error. Missing a necessary dataset for EAGLE PartType0. \n
           Either `SmoothingLength`, `Temperature`, or `InternalEnergy`. \n
           See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#parttype0 for more info.")
    }

    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.table::data.table("ID" = gas$ParticleIDs,
                                      "x"  = if(one_p_flag){gas$Coordinates[1]*.cm_to_kpc}else{gas$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                      "y"  = if(one_p_flag){gas$Coordinates[2]*.cm_to_kpc}else{gas$Coordinates[2,]*.cm_to_kpc},
                                      "z"  = if(one_p_flag){gas$Coordinates[3]*.cm_to_kpc}else{gas$Coordinates[3,]*.cm_to_kpc},
                                      "vx"  = if(one_p_flag){gas$Velocity[1]*.cms_to_kms}else{gas$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                      "vy"  = if(one_p_flag){gas$Velocity[2]*.cms_to_kms}else{gas$Velocity[2,]*.cms_to_kms},
                                      "vz"  = if(one_p_flag){gas$Velocity[3]*.cms_to_kms}else{gas$Velocity[3,]*.cms_to_kms},
                                      "Mass" = gas$Mass*.g_to_msol, # Mass in solar masses
                                      "SFR" = gas$StarFormationRate*(.g_to_msol/.s_to_yr), #SFR in Msol/yr
                                      "Density" = gas$Density*.gcm3_to_msolkpc3, # Density in Msol/kpc^3
                                      "Temperature" = gas$Temperature,
                                      "SmoothingLength" = gas$SmoothingLength*.cm_to_kpc, # Smoothing length in kpc
                                      "ThermalDispersion" = sqrt((gas$InternalEnergy*.cms_to_kms)*(.adiabatic_index - 1)),
                                      "Metallicity" = gas$Metallicity,
                                      "Hydrogen" = gas$`ElementAbundance/Hydrogen`,
                                      "Oxygen" = gas$`ElementAbundance/Oxygen`)

    gas_part$ThermalDispersion[gas_part$Temperature <= 1e4] = 11

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType4" %in% groups){
    PT4_attr = hdf5r::list.datasets(data[["PartType4"]])

    expected_names_stars = c("Coordinates", "InitialMass", "Mass", "ParticleIDs",
                             "Metallicity", "SmoothedMetallicity",
                             "StellarFormationTime", "Velocity")
    PT4_attr = PT4_attr[which(PT4_attr %in% expected_names_stars)] # trim list to only read in necessary data sets


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

    stars = .check_names(stars, type="EAGLE")

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.table::data.table("ID" = stars$ParticleIDs,
                                       "x"  = if(one_p_flag){stars$Coordinates[1]*.cm_to_kpc}else{stars$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                       "y"  = if(one_p_flag){stars$Coordinates[2]*.cm_to_kpc}else{stars$Coordinates[2,]*.cm_to_kpc},
                                       "z"  = if(one_p_flag){stars$Coordinates[3]*.cm_to_kpc}else{stars$Coordinates[3,]*.cm_to_kpc},
                                       "vx"  = if(one_p_flag){stars$Velocity[1]*.cms_to_kms}else{stars$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                       "vy"  = if(one_p_flag){stars$Velocity[2]*.cms_to_kms}else{stars$Velocity[2,]*.cms_to_kms},
                                       "vz"  = if(one_p_flag){stars$Velocity[3]*.cms_to_kms}else{stars$Velocity[3,]*.cms_to_kms},
                                       "Mass" = stars$Mass*.g_to_msol) # Mass in solar masses

    ssp = data.table::data.table("Initial_Mass" = stars$InitialMass*.g_to_msol,
                                 "Age" = as.numeric(.SFTtoAge(a=stars$StellarFormationTime,
                                                              cores=cores, H0=head$H0,
                                                              OmegaM=head$OmegaM,
                                                              OmegaL=head$OmegaL,
                                                              OmegaR=head$OmegaR)),
                                 "Metallicity" = stars$Metallicity)

    remove(stars); remove(PT4_attr)

  } else {star_part=NULL; ssp=NULL}

  return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))

}

.magneticum_read_hdf5 = function(data, head, cores){

  head$Type = "Magneticum"
  head$H0 = head$HubbleParam * 100
  names(head)[names(head) == "Omega0"] = "OmegaM"
  names(head)[names(head) == "OmegaLambda"] = "OmegaL"
  head$OmegaR = 0

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if ("PartType0" %in% groups){ # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])

    expected_names_gas = c("Coordinates", "Density", "Mass", "ParticleIDs",
                           "Metallicity",  "StarFormationRate", "Velocity",
                           "SmoothingLength", "Temperature", "InternalEnergy")
    PT0_attr = PT0_attr[which(PT0_attr %in% expected_names_gas)] # trim list to only read in necessary data sets

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

    gas = .check_names(gas, type="Magneticum")

    magneticum_gas_names = c("SmoothingLength", "Temperature", "InternalEnergy")
    if (!all(magneticum_gas_names %in% names(gas))){
      stop("Error. Missing a necessary dataset for Magneticum PartType0. \n
           Either `SmoothingLength`, `Temperature` or `InternalEnergy`. \n
           See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#parttype0 for more info.")
    }

    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.table::data.table("ID" = gas$ParticleIDs,
                                      "x"  = if(one_p_flag){gas$Coordinates[1]*.cm_to_kpc}else{gas$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                      "y"  = if(one_p_flag){gas$Coordinates[2]*.cm_to_kpc}else{gas$Coordinates[2,]*.cm_to_kpc},
                                      "z"  = if(one_p_flag){gas$Coordinates[3]*.cm_to_kpc}else{gas$Coordinates[3,]*.cm_to_kpc},
                                      "vx"  = if(one_p_flag){gas$Velocity[1]*.cms_to_kms}else{gas$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                      "vy"  = if(one_p_flag){gas$Velocity[2]*.cms_to_kms}else{gas$Velocity[2,]*.cms_to_kms},
                                      "vz"  = if(one_p_flag){gas$Velocity[3]*.cms_to_kms}else{gas$Velocity[3,]*.cms_to_kms},
                                      "Mass" = gas$Mass*.g_to_msol, # Mass in solar masses
                                      "SFR" = gas$StarFormationRate*(.g_to_msol/.s_to_yr), #SFR in Msol/yr
                                      "Density" = gas$Density*.gcm3_to_msolkpc3, # Density in Msol/kpc^3
                                      "Temperature" = gas$Temperature,
                                      "SmoothingLength" = gas$SmoothingLength*.cm_to_kpc, # Smoothing length in kpc
                                      "ThermalDispersion" = sqrt((gas$InternalEnergy*.cms_to_kms)*(.adiabatic_index - 1)),
                                      "Metallicity" = gas$Metallicity,
                                      "Hydrogen" = gas$`ElementAbundance/Hydrogen`,
                                      "Oxygen" =  gas$`ElementAbundance/Oxygen`)

    gas_part$ThermalDispersion[gas_part$Temperature <= 1e4] = 11

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType4" %in% groups){
    PT4_attr = hdf5r::list.datasets(data[["PartType4"]])

    expected_names_stars = c("Coordinates", "InitialMass", "Mass", "ParticleIDs",
                             "Metallicity", "StellarFormationTime", "Velocity")
    PT4_attr = PT4_attr[which(PT4_attr %in% expected_names_stars)] # trim list to only read in necessary data sets

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

    stars = .check_names(stars, type="Magneticum")

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.table::data.table("ID" = stars$ParticleIDs,
                                       "x"  = if(one_p_flag){stars$Coordinates[1]*.cm_to_kpc}else{stars$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                       "y"  = if(one_p_flag){stars$Coordinates[2]*.cm_to_kpc}else{stars$Coordinates[2,]*.cm_to_kpc},
                                       "z"  = if(one_p_flag){stars$Coordinates[3]*.cm_to_kpc}else{stars$Coordinates[3,]*.cm_to_kpc},
                                       "vx"  = if(one_p_flag){stars$Velocity[1]*.cms_to_kms}else{stars$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                       "vy"  = if(one_p_flag){stars$Velocity[2]*.cms_to_kms}else{stars$Velocity[2,]*.cms_to_kms},
                                       "vz"  = if(one_p_flag){stars$Velocity[3]*.cms_to_kms}else{stars$Velocity[3,]*.cms_to_kms},
                                       "Mass" = stars$Mass*.g_to_msol) # Mass in solar masses

    ssp = data.table::data.table("Initial_Mass" = stars$InitialMass*.g_to_msol,
                                 "Age" = as.numeric(.SFTtoAge(a=stars$StellarFormationTime,
                                                              cores=cores, H0=head$H0,
                                                              OmegaM=head$OmegaM,
                                                              OmegaL=head$OmegaL,
                                                              OmegaR=head$OmegaR)),
                                 "Metallicity" = stars$Metallicity)

    remove(stars); remove(PT4_attr)

  } else {star_part=NULL; ssp=NULL}

  return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))

}

.horizonagn_read_hdf5 = function(data, head, cores){

  head$Type = "Horizon-AGN"
  head$H0 = head$HubbleParam * 100
  names(head)[names(head) == "Omega0"] = "OmegaM"
  names(head)[names(head) == "OmegaLambda"] = "OmegaL"
  head$OmegaR = 0

  head$Time = 1/(1+head$Redshift)

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if ("PartType0" %in% groups){ # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])

    expected_names_gas = c("Coordinates", "Density", "Mass", "ParticleIDs",
                           "ElementAbundance/Oxygen", "SmoothedElementAbundance/Oxygen",
                           "ElementAbundance/Hydrogen", "SmoothedElementAbundance/Hydrogen",
                           "Metallicity", "SmoothedMetallicity",
                           "StarFormationRate", "Velocity", "Temperature",
                           "Pressure")
    PT0_attr = PT0_attr[which(PT0_attr %in% expected_names_gas)] # trim list to only read in necessary data sets

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

    gas = .check_names(gas, type="HAGN")
    horizon_gas_names = c("Temperature", "Pressure")
    if (!all(horizon_gas_names %in% names(gas))){
      stop("Error. Missing a necessary dataset for HorizonAGN PartType0. \n
           Missing `Temperature` or `Pressure`. \n
           See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#parttype0 for more info.")
    }

    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.table::data.table("ID" = seq(1, length(gas$ParticleIDs)),
                                      "x"  = if(one_p_flag){gas$Coordinates[1]*.cm_to_kpc}else{gas$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                      "y"  = if(one_p_flag){gas$Coordinates[2]*.cm_to_kpc}else{gas$Coordinates[2,]*.cm_to_kpc},
                                      "z"  = if(one_p_flag){gas$Coordinates[3]*.cm_to_kpc}else{gas$Coordinates[3,]*.cm_to_kpc},
                                      "vx"  = if(one_p_flag){gas$Velocity[1]*.cms_to_kms}else{gas$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                      "vy"  = if(one_p_flag){gas$Velocity[2]*.cms_to_kms}else{gas$Velocity[2,]*.cms_to_kms},
                                      "vz"  = if(one_p_flag){gas$Velocity[3]*.cms_to_kms}else{gas$Velocity[3,]*.cms_to_kms},
                                      "Mass" = gas$Mass*.g_to_msol, # Mass in solar masses
                                      "SFR" = gas$StarFormationRate*(.g_to_msol/.s_to_yr), #SFR in Msol/yr
                                      "Density" = gas$Density*.gcm3_to_msolkpc3, # Density in Msol/kpc^3
                                      "Temperature" = gas$Temperature,
                                      "ThermalDispersion" = sqrt((gas$Pressure*.gcm1_to_msolkm1)/(gas$Density*.gcm3_to_msolkm3)),
                                      "SmoothingLength" = 2*(((3/(4*pi))*((gas$Mass*.g_to_msol) / (gas$Density*.gcm3_to_msolkpc3)))^(1/3)), # smoothing length based on mass/density in units of kpc
                                      "Metallicity" = gas$Metallicity,
                                      "Hydrogen" = gas$`ElementAbundance/Hydrogen`,
                                      "Oxygen" =  gas$`ElementAbundance/Oxygen`)

    gas_part$ThermalDispersion[gas_part$Temperature <= 1e4] = 11

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType4" %in% groups){
    PT4_attr = hdf5r::list.datasets(data[["PartType4"]])

    expected_names_stars = c("Coordinates", "InitialMass", "Mass", "ParticleIDs",
                             "Metallicity", "SmoothedMetallicity",
                             "StellarFormationTime", "Velocity")
    PT4_attr = PT4_attr[which(PT4_attr %in% expected_names_stars)] # trim list to only read in necessary data sets

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

    stars = .check_names(stars, type="HAGN")

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.table::data.table("ID" = stars$ParticleIDs,
                                       "x"  = if(one_p_flag){stars$Coordinates[1]*.cm_to_kpc}else{stars$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                       "y"  = if(one_p_flag){stars$Coordinates[2]*.cm_to_kpc}else{stars$Coordinates[2,]*.cm_to_kpc},
                                       "z"  = if(one_p_flag){stars$Coordinates[3]*.cm_to_kpc}else{stars$Coordinates[3,]*.cm_to_kpc},
                                       "vx"  = if(one_p_flag){stars$Velocity[1]*.cms_to_kms}else{stars$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                       "vy"  = if(one_p_flag){stars$Velocity[2]*.cms_to_kms}else{stars$Velocity[2,]*.cms_to_kms},
                                       "vz"  = if(one_p_flag){stars$Velocity[3]*.cms_to_kms}else{stars$Velocity[3,]*.cms_to_kms},
                                       "Mass" = stars$Mass*.g_to_msol) # Mass in solar masses

    ssp = data.table::data.table("Initial_Mass" = stars$InitialMass*.g_to_msol,
                                 "Age" = as.numeric(.SFTtoAge(a=stars$StellarFormationTime,
                                                              cores=cores, H0=head$H0,
                                                              OmegaM=head$OmegaM,
                                                              OmegaL=head$OmegaL,
                                                              OmegaR=head$OmegaR)),
                                 "Metallicity" = stars$Metallicity)

    remove(stars); remove(PT4_attr)

  } else {star_part=NULL; ssp=NULL}

  return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))

}

.illustristng_read_hdf5 = function(data, head, cores){

  head$Type = "Illustris-TNG"
  head$H0 = head$HubbleParam * 100
  names(head)[names(head) == "Omega0"] = "OmegaM"
  names(head)[names(head) == "OmegaLambda"] = "OmegaL"
  head$OmegaR = 0

  head$Time = 1/(1+head$Redshift)

  if (!"NumPart_Total" %in% names(head)){
    head$NumPart_Total = head$NumPart_ThisFile
  }

  groups = hdf5r::list.groups(data)                         # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")]  # Pick out PartTypeX groups

  if ("PartType0" %in% groups){                             # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])    # get list of fields

    expected_names_gas = c("Coordinates", "Density", "Masses", "ParticleIDs",
                           "GFM_Metals", "GFM_Metallicity",
                           "StarFormationRate", "Velocities",
                           "ElectronAbundance", "InternalEnergy")
    PT0_attr = PT0_attr[which(PT0_attr %in% expected_names_gas)] # trim list to only read in necessary data sets

    n_gas_prop = length(PT0_attr)                           # how many fields?
    gas = vector("list", n_gas_prop)                        # make a vector
    names(gas) = PT0_attr                                   # assign field names to elements

    # loop through the fields,
    # assign them their conversion factor attributes,
    # and then convert their data
    for (i in 1:n_gas_prop){
      aexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "a_scaling")},
                      error = function(e){NULL})
      hexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "h_scaling")},
                      error = function(e){NULL})
      cgs  = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "to_cgs")},
                      error = function(e){NULL})
      if (any(is.null(aexp), is.null(hexp), is.null(cgs))){
        aexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "aexp-scale-exponent")},
                        error = function(e){NULL})
        hexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "h-scale-exponent")},
                        error = function(e){NULL})
        cgs  = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "CGSConversionFactor")},
                        error = function(e){NULL})
        if (any(is.null(aexp), is.null(hexp), is.null(cgs))){
          stop("Error: Attributes listed incorrectly. \n
               Please check that the scaling factors for conversion between comoving and physical coordinates are included in the input file as attributes for each dataset. \n
               See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#hydrodynamical-simulations for more details.")}
      }
      if (cgs == 0) {cgs = 1}                                                   # converting 0 values to 1
      gas[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType0/",PT0_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    gas = .check_names(gas, type="TNG")
    illustris_gas_names = c("InternalEnergy", "ElectronAbundance")
    if (!all(illustris_gas_names %in% names(gas))){
      stop("Error. Missing a necessary dataset for IllustrisTNG PartType0. \n
           Missing `InternalEnergy` or `ElectronAbundance`. \n
           See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#parttype0 for more info.")
    }

    # check if the Coordinate field has any dimensions
    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.table::data.table("ID" = gas$ParticleIDs,
                                      "x"  = if(one_p_flag){gas$Coordinates[1]*.cm_to_kpc}else{gas$Coordinates[1,]*.cm_to_kpc},   # Coordinates in kpc
                                      "y"  = if(one_p_flag){gas$Coordinates[2]*.cm_to_kpc}else{gas$Coordinates[2,]*.cm_to_kpc},
                                      "z"  = if(one_p_flag){gas$Coordinates[3]*.cm_to_kpc}else{gas$Coordinates[3,]*.cm_to_kpc},
                                      "vx"  = if(one_p_flag){gas$Velocity[1]*.cms_to_kms}else{gas$Velocity[1,]*.cms_to_kms},  # Velocity in km/s
                                      "vy"  = if(one_p_flag){gas$Velocity[2]*.cms_to_kms}else{gas$Velocity[2,]*.cms_to_kms},
                                      "vz"  = if(one_p_flag){gas$Velocity[3]*.cms_to_kms}else{gas$Velocity[3,]*.cms_to_kms},
                                      "Mass" = gas$Mass*.g_to_msol,                                                             # Mass in solar Mass
                                      "SFR" = gas$StarFormationRate*(.g_to_msol/.s_to_yr),                                        # SFR in Msol/yr
                                      "Density" = gas$Density*.gcm3_to_msolkpc3,                                                  # Density in Msol/kpc^3
                                      "Temperature" = (.adiabatic_index-1)*(gas$InternalEnergy/.Boltzmann_constant)*(4*.mass_of_proton/(1 + gas$`ElementAbundance/Hydrogen`*(3 + 4*gas$ElectronAbundance))),
                                      "ThermalDispersion" = sqrt((gas$InternalEnergy*.cms_to_kms)*(.adiabatic_index - 1)),
                                      "SmoothingLength" = 2*(((3/(4*pi))*((gas$Mass*.g_to_msol) / (gas$Density*.gcm3_to_msolkpc3)))^(1/3)), # smoothing length based on mass/density in units of kpc
                                      "Metallicity" = gas$Metallicity,
                                      "Hydrogen" = gas$`ElementAbundance/Hydrogen`,
                                      "Oxygen" = gas$`ElementAbundance/Oxygen`)

    gas_part$ThermalDispersion[gas_part$Temperature <= 1e4] = 11

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType4" %in% groups){

    PT4_attr = hdf5r::list.datasets(data[["PartType4"]])

    expected_names_stars = c("Coordinates", "GFM_InitialMass", "Masses", "ParticleIDs",
                             "GFM_Metallicity", "GFM_StellarFormationTime", "Velocities")
    PT4_attr = PT4_attr[which(PT4_attr %in% expected_names_stars)] # trim list to only read in necessary data sets

    n_star_prop = length(PT4_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT4_attr

    for (i in 1:n_star_prop){
      aexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "a_scaling")},
                      error = function(e){NULL})
      hexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "h_scaling")},
                      error = function(e){NULL})
      cgs  = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "to_cgs")},
                      error = function(e){NULL})
      if (any(is.null(aexp), is.null(hexp), is.null(cgs))){
        aexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "aexp-scale-exponent")},
                        error = function(e){NULL})
        hexp = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "h-scale-exponent")},
                        error = function(e){NULL})
        cgs  = tryCatch(expr = {hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "CGSConversionFactor")},
                        error = function(e){NULL})
        if (any(is.null(aexp), is.null(hexp), is.null(cgs))){
          stop("Error: Attributes listed incorrectly. \n
               Please check that the scaling factors for conversion between comoving and physical coordinates are included in the input file as attributes for each dataset. \n
               See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#hydrodynamical-simulations for more details.")}
      }
      if (cgs == 0) {cgs = 1}                                                   # converting 0 values to 1
      stars[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType4/",PT4_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    stars = .check_names(stars, type="TNG")

    # check if the Coordinate field has any dimensions
    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.table::data.table("ID" = stars$ParticleIDs,
                                       "x"  = if(one_p_flag){stars$Coordinates[1]*.cm_to_kpc}else{stars$Coordinates[1,]*.cm_to_kpc},  # Coordinates in kpc
                                       "y"  = if(one_p_flag){stars$Coordinates[2]*.cm_to_kpc}else{stars$Coordinates[2,]*.cm_to_kpc},
                                       "z"  = if(one_p_flag){stars$Coordinates[3]*.cm_to_kpc}else{stars$Coordinates[3,]*.cm_to_kpc},
                                       "vx"  = if(one_p_flag){stars$Velocity[1]*.cms_to_kms}else{stars$Velocity[1,]*.cms_to_kms}, # Velocity in km/s
                                       "vy"  = if(one_p_flag){stars$Velocity[2]*.cms_to_kms}else{stars$Velocity[2,]*.cms_to_kms},
                                       "vz"  = if(one_p_flag){stars$Velocity[3]*.cms_to_kms}else{stars$Velocity[3,]*.cms_to_kms},
                                       "Mass" = stars$Mass*.g_to_msol,                                                              # Mass in solar Mass
                                       "SFT" = stars$StellarFormationTime)

    ssp = data.table::data.table("Initial_Mass" = stars$InitialMass*.g_to_msol,
                                 "Age" = as.numeric(.SFTtoAge(a=abs(stars$StellarFormationTime),
                                                              cores=cores, H0=head$H0,
                                                              OmegaM=head$OmegaM,
                                                              OmegaL=head$OmegaL,
                                                              OmegaR=head$OmegaR)),                                 "Metallicity" = stars$Metallicity,
                                 "SFT" = stars$StellarFormationTime)

    # remove stellar wind particles and drop unneeded SFT columns
    star_part = star_part[SFT >= 0]
    star_part = star_part[,SFT:=NULL]
    ssp = ssp[SFT >= 0]
    ssp = ssp[,SFT:=NULL]

    remove(stars); remove(PT4_attr)

  } else {star_part=NULL; ssp=NULL}

  return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))
}

.colibre_read_hdf5 = function(data, head, cores){

  head$Type = "Colibre"
  head$H0 = head$HubbleParam * 100
  names(head)[names(head) == "Omega0"] = "OmegaM"
  names(head)[names(head) == "OmegaLambda"] = "OmegaL"
  head$OmegaR = 0

  groups = hdf5r::list.groups(data) # What particle data is present?
  groups = groups[stringr::str_detect(groups, "PartType")] # Pick out PartTypeX groups

  if ("PartType0" %in% groups){ # If gas particles are present in the file

    PT0_attr = hdf5r::list.datasets(data[["PartType0"]])

    expected_names_gas = c("Coordinates", "Densities", "Masses", "ParticleIDs",
                           "ElementMassFractions/Oxygen", "ElementMassFractions/Hydrogen",
                           "MetalMassFractions",
                           "StarFormationRates", "Velocities", "SmoothingLengths",
                           "Temperatures", "InternalEnergies")
    PT0_attr = PT0_attr[which(PT0_attr %in% expected_names_gas)] # trim list to only read in necessary data sets

    n_gas_prop = length(PT0_attr)
    gas = vector("list", n_gas_prop)
    names(gas) = PT0_attr

    for (i in 1:n_gas_prop){
      aexp = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "a-scale exponent")
      hexp = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "h-scale exponent")
      cgs  = hdf5r::h5attr(data[[paste0("PartType0/",PT0_attr[i])]], "Conversion factor to CGS (not including cosmological corrections)")
      gas[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType0/",PT0_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    gas = .check_names(gas, type="Colibre")

    colibre_gas_names = c("SmoothingLength", "Temperature", "InternalEnergy")
    if (!all(colibre_gas_names %in% names(gas))){
      stop("Error. Missing a necessary dataset for COLIBRE PartType0. \n
           Either `SmoothingLength`, `Temperature`, or `InternalEnergy`. \n
           See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#parttype0 for more info.")
    }

    one_p_flag = FALSE
    if (is.null(dim(gas$Coordinates))){one_p_flag = TRUE}

    gas_part = data.table::data.table("ID" = gas$ParticleIDs,
                                      "x"  = if(one_p_flag){gas$Coordinates[1]*.cm_to_kpc}else{gas$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                      "y"  = if(one_p_flag){gas$Coordinates[2]*.cm_to_kpc}else{gas$Coordinates[2,]*.cm_to_kpc},
                                      "z"  = if(one_p_flag){gas$Coordinates[3]*.cm_to_kpc}else{gas$Coordinates[3,]*.cm_to_kpc},
                                      "vx"  = if(one_p_flag){gas$Velocity[1]*.cms_to_kms}else{gas$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                      "vy"  = if(one_p_flag){gas$Velocity[2]*.cms_to_kms}else{gas$Velocity[2,]*.cms_to_kms},
                                      "vz"  = if(one_p_flag){gas$Velocity[3]*.cms_to_kms}else{gas$Velocity[3,]*.cms_to_kms},
                                      "Mass" = gas$Mass*.g_to_msol, # Mass in solar masses
                                      "SFR" = gas$StarFormationRate*(.g_to_msol/.s_to_yr), #SFR in Msol/yr
                                      "Density" = gas$Density*.gcm3_to_msolkpc3, # Density in Msol/kpc^3
                                      "Temperature" = gas$Temperature,
                                      "SmoothingLength" = gas$SmoothingLength*.cm_to_kpc, # Smoothing length in kpc
                                      "ThermalDispersion" = sqrt((gas$InternalEnergy*.cms_to_kms)*(.adiabatic_index - 1)),
                                      "Metallicity" = gas$Metallicity,
                                      "Hydrogen" = gas$`ElementAbundance/Hydrogen`,
                                      "Oxygen" = gas$`ElementAbundance/Oxygen`)

    gas_part$ThermalDispersion[gas_part$Temperature <= 1e4] = 11

    remove(gas); remove(PT0_attr)

  } else {gas_part=NULL}

  if ("PartType4" %in% groups){
    PT4_attr = hdf5r::list.datasets(data[["PartType4"]])

    expected_names_stars = c("Coordinates", "InitialMasses", "Masses", "ParticleIDs",
                             "MetalMassFractions", "BirthScaleFactors", "Velocities")
    PT4_attr = PT4_attr[which(PT4_attr %in% expected_names_stars)] # trim list to only read in necessary data sets

    n_star_prop = length(PT4_attr)
    stars = vector("list", n_star_prop)
    names(stars) = PT4_attr

    for (i in 1:n_star_prop){
      aexp = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "a-scale exponent")
      hexp = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "h-scale exponent")
      cgs  = hdf5r::h5attr(data[[paste0("PartType4/",PT4_attr[i])]], "Conversion factor to CGS (not including cosmological corrections)")
      stars[[i]] =
        hdf5r::readDataSet(data[[paste0("PartType4/",PT4_attr[i])]]) * head$Time^(aexp) * head$HubbleParam^(hexp) * cgs
    }

    stars = .check_names(stars, type="Colibre")

    one_p_flag = FALSE
    if (is.null(dim(stars$Coordinates))){one_p_flag = TRUE}

    star_part = data.table::data.table("ID" = stars$ParticleIDs,
                                       "x"  = if(one_p_flag){stars$Coordinates[1]*.cm_to_kpc}else{stars$Coordinates[1,]*.cm_to_kpc}, # Coordinates in kpc
                                       "y"  = if(one_p_flag){stars$Coordinates[2]*.cm_to_kpc}else{stars$Coordinates[2,]*.cm_to_kpc},
                                       "z"  = if(one_p_flag){stars$Coordinates[3]*.cm_to_kpc}else{stars$Coordinates[3,]*.cm_to_kpc},
                                       "vx"  = if(one_p_flag){stars$Velocity[1]*.cms_to_kms}else{stars$Velocity[1,]*.cms_to_kms}, # Velocities in km/s
                                       "vy"  = if(one_p_flag){stars$Velocity[2]*.cms_to_kms}else{stars$Velocity[2,]*.cms_to_kms},
                                       "vz"  = if(one_p_flag){stars$Velocity[3]*.cms_to_kms}else{stars$Velocity[3,]*.cms_to_kms},
                                       "Mass" = stars$Mass*.g_to_msol) # Mass in solar masses

    ssp = data.table::data.table("Initial_Mass" = stars$InitialMass*.g_to_msol,
                                 "Age" = as.numeric(.SFTtoAge(a=stars$StellarFormationTime,
                                                              cores=cores, H0=head$H0,
                                                              OmegaM=head$OmegaM,
                                                              OmegaL=head$OmegaL,
                                                              OmegaR=head$OmegaR)),
                                 "Metallicity" = stars$Metallicity)

    remove(stars); remove(PT4_attr)

  } else {star_part=NULL; ssp=NULL}

  return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))

}

# Function to check existing names in a data set and convert if necessary
.check_names = function(particle_list, type){

  current_names = names(particle_list)

  if ("BirthScaleFactors" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "BirthScaleFactors")] <- "StellarFormationTime"
    names(particle_list) <- current_names
  }

  if ("SmoothedMetallicity" %in% current_names & type == "EAGLE" |
      "SmoothedMetallicity" %in% current_names & type == "HAGN"){
    current_names[which(current_names == "SmoothedMetallicity")] <- "Metallicity"
    names(particle_list) <- current_names
  }

  if ("MetalMassFractions" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "MetalMassFractions")] <- "Metallicity"
    names(particle_list) <- current_names
  }

  if ("SmoothedElementAbundance/Hydrogen" %in% current_names & type == "EAGLE" |
      "SmoothedElementAbundance/Hydrogen" %in% current_names & type == "HAGN" ){
    current_names[which(current_names == "SmoothedElementAbundance/Hydrogen")] <- "ElementAbundance/Hydrogen"
    names(particle_list) <- current_names
  }

  if ("ElementMassFractions/Hydrogen" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "ElementMassFractions/Hydrogen")] <- "ElementAbundance/Hydrogen"
    names(particle_list) <- current_names
  }

  if ("SmoothedElementAbundance/Oxygen" %in% current_names & type == "EAGLE" |
      "SmoothedElementAbundance/Oxygen" %in% current_names & type == "HAGN"){
    current_names[which(current_names == "SmoothedElementAbundance/Oxygen")] <- "ElementAbundance/Oxygen"
    names(particle_list) <- current_names
  }

  if ("ElementMassFractions/Oxygen" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "ElementMassFractions/Oxygen")] <- "ElementAbundance/Oxygen"
    names(particle_list) <- current_names
  }

  if ("Velocities" %in% current_names & type == "TNG" |
      "Velocities" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "Velocities")] <- "Velocity"
    names(particle_list) <- current_names
  }

  if ("Densities" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "Densities")] <- "Density"
    names(particle_list) <- current_names
  }

  if ("Masses" %in% current_names & type == "TNG" |
      "Masses" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "Masses")] <- "Mass"
    names(particle_list) <- current_names
  }

  if ("GFM_StellarFormationTime" %in% current_names & type == "TNG"){
    current_names[which(current_names == "GFM_StellarFormationTime")] <- "StellarFormationTime"
    names(particle_list) <- current_names
  }

  if ("StarFormationRates" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "StarFormationRates")] <- "StarFormationRate"
    names(particle_list) <- current_names
  }

  if ("SmoothingLengths" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "SmoothingLengths")] <- "SmoothingLength"
    names(particle_list) <- current_names
  }

  if ("GFM_Metallicity" %in% current_names & type == "TNG"){
    current_names[which(current_names == "GFM_Metallicity")] <- "Metallicity"
    names(particle_list) <- current_names
  }

 if ("GFM_InitialMass" %in% current_names & type == "TNG"){
    current_names[which(current_names == "GFM_InitialMass")] <- "InitialMass"
    names(particle_list) <- current_names
  }

  if ("InitialMasses" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "InitialMasses")] <- "InitialMass"
    names(particle_list) <- current_names
  }

  if ("InternalEnergies" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "InternalEnergies")] <- "InternalEnergy"
    names(particle_list) <- current_names
  }

  if ("Temperatures" %in% current_names & type == "Colibre"){
    current_names[which(current_names == "Temperatures")] <- "Temperature"
    names(particle_list) <- current_names
  }

  if ("GFM_Metals" %in% current_names & type == "TNG"){
    id_to_remove = which(current_names == "GFM_Metals")

    one_p_flag = FALSE
    if (is.null(dim(particle_list$Coordinates))){one_p_flag = TRUE}

    particle_list$`ElementAbundance/Oxygen` = if(one_p_flag){particle_list$GFM_Metals[5]}else{particle_list$GFM_Metals[5,]}
    particle_list$`ElementAbundance/Hydrogen` = if(one_p_flag){particle_list$GFM_Metals[1]}else{particle_list$GFM_Metals[1,]}

    particle_list = particle_list[-id_to_remove]
  }

  if (!is.null(nrow(particle_list$Metallicity)) & type == "Magneticum" |
      length(particle_list$Metallicity)[1] == 11 & length(particle_list$ParticleIDs) != 11 & type == "Magneticum"){

    one_p_flag = FALSE
    if (is.null(dim(particle_list$Coordinates))){one_p_flag = TRUE}

    Metallicity = if(one_p_flag){(sum(particle_list$Metallicity[2:11]))/(particle_list$Mass)}else{(colSums(particle_list$Metallicity[2:11,]))/(particle_list$Mass)}
    Hydrogen    = if(one_p_flag){(particle_list$Mass - sum(particle_list$Metallicity)) / particle_list$Mass}else{(particle_list$Mass - colSums(particle_list$Metallicity)) / particle_list$Mass}
    Oxygen      = if(one_p_flag){particle_list$Metallicity[4] / particle_list$Mass}else{particle_list$Metallicity[4,] / particle_list$Mass}

    particle_list$Metallicity = Metallicity
    particle_list$`ElementAbundance/Oxygen` = Hydrogen
    particle_list$`ElementAbundance/Hydrogen` = Oxygen

  }

  expected_names_gas = c("Coordinates", "Density", "Mass", "ParticleIDs",
                         "ElementAbundance/Oxygen",
                         "ElementAbundance/Hydrogen", "Metallicity",
                         "StarFormationRate", "Velocity")

  expected_names_stars = c("Coordinates", "InitialMass", "Mass", "ParticleIDs",
                           "Metallicity", "StellarFormationTime", "Velocity")

  if (!all(expected_names_gas %in% names(particle_list)) &
      !all(expected_names_stars %in% names(particle_list))){
    stop("Error. A key dataset is missing that is necessary for processing. \n
          Please check https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html#hydrodynamical-simulations for an expected list. \n")
  }

  return(particle_list)

}


# Function for computing the stellar age from the formation time in parallel
.SFTtoAge = function(a, H0, OmegaM, OmegaL, OmegaR, cores=1){
  cosdist = function(x) { return (celestial::cosdistTravelTime((1 / x) - 1, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR)); }
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
