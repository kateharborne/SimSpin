# Author: Kate Harborne
# Date: 27/10/2020
# Title: Utilities functions (i.e. hidden functions from the user)

# Some useful constants:
.lsol_to_erg    = 3.828e33
.mpc_to_cm      = 3.08568e+24
.speed_of_light = 299792.458

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
  cosdist = function(x) { return (celestial::cosdistTravelTime((1 / x) - 1)); }
  if (cores > 1) {
    cl = snow::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    i = integer()
    output = foreach(i = 1:length(a), .packages = "celestial") %dopar% { cosdist(a[i]) }
    closeAllConnections()
  }
  else {
    output = lapply(a, cosdist)
  }
  return(output)
}

# Function to calculate spectral weights
.spectra_weights = function(Metallicity, Age, Template){
  temp_metals = as.numeric(Template[[1]])
  temp_ages   = as.numeric(Template[[2]])
  sw = function(metallicity, age){
    Z = as.numeric(ProSpect::interp_quick(metallicity, temp_metals, log = TRUE))
    A = as.numeric(ProSpect::interp_quick(age * 1e9, temp_ages, log = TRUE))
    return(c(Z, A))
    }
  weights = mapply(sw, Metallicity, Age)
  return(weights)
}

# Function to generate spectra from spectral weights
.compute_spectra = function(Weights, Mass, Template){
  cs = function(w, mass){
    spectra = ((Template$Zspec[[w[2]]][w[6],] * w[4]*w[8]) +
               (Template$Zspec[[w[2]]][w[5],] * w[4]*w[7]) +
               (Template$Zspec[[w[1]]][w[6],] * w[3]*w[8]) +
               (Template$Zspec[[w[1]]][w[5],] * w[3]*w[7])) * 1e10 * mass
  return(spectra)
  }

  intrinsic_spectra = matrix(unlist(mapply(cs, Weights, Mass)), nrow=length(Mass), byrow=T)

  return(intrinsic_spectra)
}

# Function to generate spectra
.part_spec = function(Metallicity, Age, Mass, Template, cores){
  f = function(metallicity, age, mass) {
    Z = as.numeric(ProSpect::interp_quick(metallicity, Template$Z, log = TRUE))
    A = as.numeric(ProSpect::interp_quick(age * 1e9, Template$Age, log = TRUE))

    weights = data.frame("hihi" = Z[4] * A[4],   # ID_lo = 1, ID_hi = 2, wt_lo = 3, wt_hi = 4
                         "hilo" = Z[4] * A[3],
                         "lohi" = Z[3] * A[4],
                         "lolo" = Z[3] * A[3])

    part_spec = array(data = NA, dim = c(1, length(Template$Wave)))

    part_spec = ((Template$Zspec[[Z[2]]][A[2],] * weights$hihi) +
                   (Template$Zspec[[Z[2]]][A[1],] * weights$hilo) +
                   (Template$Zspec[[Z[1]]][A[2],] * weights$lohi) +
                   (Template$Zspec[[Z[1]]][A[1],] * weights$lolo)) * mass * 1e10

    return(part_spec)
  }
  if (cores > 1) {
    cl = snow::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    i = integer()
    part_spec = foreach::foreach(i = 1:length(Metallicity), .packages = c("ProSpect", "SimSpin", "foreach")) %dopar% {
      f(Metallicity[i], Age[i], Mass[i])
    }
    closeAllConnections()
  }
  else {
    part_spec = mapply(f, Metallicity, Age, Mass)
  }
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

.circular_ap=function(sbin){
  ap_region = matrix(data = NA, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  x = matrix(data = rep(seq(1,sbin), each=sbin), nrow = sbin, ncol = sbin)
  y = matrix(data = rep(seq(sbin,1), sbin), nrow = sbin, ncol = sbin)
  xx = x - xcentre; yy = y - ycentre
  rr = sqrt(xx^2 + yy^2)
  ap_region[rr<= sbin/2] = 1
  return(as.vector(ap_region))
}

.hexagonal_ap=function(sbin){
  ap_region = matrix(data = NA, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
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

# Function to apply LSF to spectra
.lsf_convolution = function(observation, luminosity, lsf_sigma){

  scaled_sigma = lsf_sigma / observation$wave_res # scaling the size of the gaussian to match the pixel dimensions
  kernel = dnorm(seq(-scaled_sigma*5,scaled_sigma*5,length.out = 25), mean = 0, sd = scaled_sigma)
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
