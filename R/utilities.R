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

  part = data.frame("ID" = id,          # the particle data table
                    "x" = pos[extract], "y"=pos[extract+1], "z"=pos[extract+2],
                    "vx" = vel[extract], "vy" = vel[extract+1], "vz"=vel[extract+2],
                    "Mass" = masses)

  head = list("Npart" = c(Npart[1], 0, Npart[3], Npart[4], Npart[5], 0), # number of gas and stars
              "Time" = Time, "Redshift" = Redshift, # relevent simulation data
              "Nall" = Nall) # number of particles in the original file

  if(verbose){cat("Done reading Gadget snapshot file. \n")}

  Npart_sum = cumsum(Npart) # cumulative number of each particle type

  if(verbose){cat("Writing stellar particles... \n")}
  star_part = part[(Npart_sum[2]+1):Npart_sum[5],]

  if (Npart[1] > 0){
    if(verbose){cat("Writing gas particles... \n")}
    gas_part = part[1:Npart_sum[1],]
  } else {gas_part = NULL}

  if(verbose){cat("Done! \n")}
  return(list(star_part = star_part, gas_part = gas_part, head = head))

}

# Function for reading in Gadget HDF5 files and EAGLE HDF5 files
.read_hdf5   = function(f, cores=1, verbose = FALSE){

  data = hdf5r::h5file(f, mode="r")

  if(length(hdf5r::list.attributes(data[["Header"]])) > 17){eagle = T}else{eagle=F} # determining if EAGLE input (based on number of parameters in Header)
  if(verbose & eagle){cat("EAGLE input snapshot detected.\n")}

  if(verbose){cat("Reading in header.\n")}
  Npart         = hdf5r::h5attr(data[["Header"]], "NumPart_ThisFile")
  Massarr       = hdf5r::h5attr(data[["Header"]], "MassTable")
  Time          = hdf5r::h5attr(data[["Header"]], "Time")
  Redshift      = hdf5r::h5attr(data[["Header"]], "Redshift")
  Nall          = hdf5r::h5attr(data[["Header"]], "NumPart_Total")
  HubbleParam   = hdf5r::h5attr(data[["Header"]], "HubbleParam")

  n_stellar      = c(0,0,Npart[3],Npart[4],Npart[5],0) # sum of all "stellar" particles in the file

  # First considering the stellar properties. If Bulge and Disk particles are present,
  # we need to make sure that are ordered correctly in the output so that they can
  # be assigned the correct SEDs.
  present_stars  = which(Npart > 0)[which(Npart > 0) %in% c(3,4,5)] # which stellar groups are present in the file?

  stellar_sum    = sum(n_stellar)
  Npart_sum      = cumsum(n_stellar) # particle indices of each type
  mass_excpt     = which(Massarr > 0) # are any masses listed in the header?

  star_part = data.frame("ID" = numeric(stellar_sum),
                         "x" = numeric(stellar_sum), "y" = numeric(stellar_sum), "z" = numeric(stellar_sum),
                         "vx" = numeric(stellar_sum), "vy" = numeric(stellar_sum), "vz" = numeric(stellar_sum),
                         "Mass" = numeric(stellar_sum))

  head = list("Npart" = c(Npart[1], 0, Npart[3], Npart[4], Npart[5], 0), # number of gas and star particles
              "Time" = Time, "Redshift" = Redshift, # relevent simulation data
              "Nall" = Nall)  # number of particles in the original file

  star_pos = array(NA, dim=c(3, stellar_sum))
  star_vel = array(NA, dim=c(3, stellar_sum))

  if(verbose){cat("Reading stellar particle properties. \n")}

  for (i in 1:length(present_stars)){

    if (i == 1){
      i_start = 1; i_end = Npart_sum[present_stars[i]]
    } else {
      i_start = Npart_sum[present_stars[i-1]]+1; i_end = Npart_sum[present_stars[i]]
    } # computing the indices of the particles of each present PartType

    # reading x, y, z coordinates
    star_pos[,i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_stars[i]-1, "/Coordinates", sep="")]])

    # reading vx, vy, vz velocities (called "Velocity" in EAGLE snapshots and "Velocities" Gadget HDF5 files)
    if (eagle){
      star_vel[,i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_stars[i]-1, "/Velocity", sep="")]])
      } else {
      star_vel[,i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_stars[i]-1, "/Velocities", sep="")]])
      }

    # reading particle IDs
    star_part$ID[i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_stars[i]-1, "/ParticleIDs", sep="")]])

    # reading masses
    if (length(mass_excpt)!=0 & present_stars[i] %in% mass_excpt){ # if masses have been specified in header, use these
      star_part$Mass[i_start:i_end] = Massarr[present_stars[i]]
    } else { # else, read masses from PartType ("Mass" in EAGLE snapshots and "Masses" in Gadget HDF5 files)
      if (eagle){
        star_part$Mass[i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_stars[i]-1, "/Mass", sep="")]])
        } else {
        star_part$Mass[i_start:i_end] = hdf5r::readDataSet(data[[paste("PartType", present_stars[i]-1, "/Masses", sep="")]])
      }
    }
  }

  if (eagle){ # reading the details from EAGLE files for simple stellar population
    ssp = data.frame("Initial_Mass"=numeric(length=stellar_sum), "Age"=numeric(length=stellar_sum),
                     "Metallicity"=numeric(length=stellar_sum))
    ssp$Initial_Mass = hdf5r::readDataSet(data[["PartType4/InitialMass"]])
    Stellar_Formation_Time = hdf5r::readDataSet(data[["PartType4/StellarFormationTime"]])
    ssp$Age = as.numeric(.SFTtoAge(a = Stellar_Formation_Time, cores = cores))
    ssp$Metallicity = hdf5r::readDataSet(data[["PartType4/SmoothedMetallicity"]])

    # coordinate transform from co-moving to physical coordinates
    a_coor = hdf5r::h5attr(data[["PartType4/Coordinates"]], "aexp-scale-exponent")
    h_coor = hdf5r::h5attr(data[["PartType4/Coordinates"]], "h-scale-exponent")
    a_velo = hdf5r::h5attr(data[["PartType4/Velocity"]], "aexp-scale-exponent")
    h_velo = hdf5r::h5attr(data[["PartType4/Velocity"]], "h-scale-exponent")
    a_mass = hdf5r::h5attr(data[["PartType4/Mass"]], "aexp-scale-exponent")
    h_mass = hdf5r::h5attr(data[["PartType4/Mass"]], "h-scale-exponent")

    star_pos = star_pos * Time^(a_coor) * HubbleParam^(h_coor) * 1e3
    star_vel = star_vel * Time^(a_velo) * HubbleParam^(h_velo)
    star_part$Mass = star_part$Mass * Time^(a_mass) * HubbleParam^(h_mass)
    ssp$Initial_Mass = ssp$Initial_Mass * Time^(a_mass) * HubbleParam^(h_mass)
  }

  # sort stellar coord & vel into x, y, ..., vz
  if (stellar_sum > 1){
    star_part$x = star_pos[1,]; star_part$vx = star_vel[1,]
    star_part$y = star_pos[2,]; star_part$vy = star_vel[2,]
    star_part$z = star_pos[3,]; star_part$vz = star_vel[3,]
  } else {
    star_part$x = star_pos[1]; star_part$vx = star_vel[1]
    star_part$y = star_pos[2]; star_part$vy = star_vel[2]
    star_part$z = star_pos[3]; star_part$vz = star_vel[3]
  }

  if (Npart[1] > 0){ # If gas is present in the simulation
    if(verbose){cat("Reading gas particle properties. \n")}
    gas_part = data.frame("ID" = numeric(Npart[1]),
                          "x" = numeric(Npart[1]), "y" = numeric(Npart[1]), "z" = numeric(Npart[1]),
                          "vx" = numeric(Npart[1]), "vy" = numeric(Npart[1]), "vz" = numeric(Npart[1]),
                          "Mass" = numeric(Npart[1]))
    gas_pos = array(NA, dim=c(3, Npart[1]))
    gas_vel = array(NA, dim=c(3, Npart[1]))

    gas_part$ID = hdf5r::readDataSet(data[["PartType0/ParticleIDs"]]) # reading particle IDs

    # reading in Masses
    if (length(mass_excpt)!=0 & 1 %in% mass_excpt){
      gas_part$Mass = Massarr[1]
    } else {
      if (eagle){gas_part$Mass = hdf5r::readDataSet(data[["PartType0/Mass"]])} else {
        gas_part$Mass = hdf5r::readDataSet(data[["PartType0/Masses"]])
      }
    }

    # reading in particle positions
    gas_pos = hdf5r::readDataSet(data[["PartType0/Coordinates"]])

    # reading vx, vy, vz velocities (called "Velocity" in EAGLE snapshots and "Velocities" Gadget HDF5 files)
    if (eagle){
      gas_vel= hdf5r::readDataSet(data[["PartType0/Velocity"]])
      } else {
      gas_vel = hdf5r::readDataSet(data[["PartType0/Velocities"]])
      }

    # coordinate transform from co-moving to physical coordinates
    if (eagle){
      a_coor = hdf5r::h5attr(data[["PartType0/Coordinates"]], "aexp-scale-exponent")
      h_coor = hdf5r::h5attr(data[["PartType0/Coordinates"]], "h-scale-exponent")
      a_velo = hdf5r::h5attr(data[["PartType0/Velocity"]], "aexp-scale-exponent")
      h_velo = hdf5r::h5attr(data[["PartType0/Velocity"]], "h-scale-exponent")
      a_mass = hdf5r::h5attr(data[["PartType0/Mass"]], "aexp-scale-exponent")
      h_mass = hdf5r::h5attr(data[["PartType0/Mass"]], "h-scale-exponent")

      gas_pos = gas_pos * Time^(a_coor) * HubbleParam^(h_coor) * 1e3
      gas_vel = gas_vel * Time^(a_velo) * HubbleParam^(h_velo)
      gas_part$Mass = gas_part$Mass * Time^(a_mass) * HubbleParam^(h_mass)
    }

    # sort gas coord & vel into x, y, ..., vz
    if (Npart[1] > 1){
      gas_part$x = gas_pos[1,]; gas_part$vx = gas_vel[1,]
      gas_part$y = gas_pos[2,]; gas_part$vy = gas_vel[2,]
      gas_part$z = gas_pos[3,]; gas_part$vz = gas_vel[3,]
    } else {
      gas_part$x = gas_pos[1]; gas_part$vx = gas_vel[1]
      gas_part$y = gas_pos[2]; gas_part$vy = gas_vel[2]
      gas_part$z = gas_pos[3]; gas_part$vz = gas_vel[3]
    }

  } else {gas_part = NULL} # if no gas in sim, return gas_part = NULL

  hdf5r::h5close(data)

  if (verbose){cat("Done reading HDF5 snapshot file. \n")}
  if (eagle){
    return(list(star_part=star_part, gas_part=gas_part, head=head, ssp=ssp))
    } else {
    return(list(star_part=star_part, gas_part=gas_part, head=head))}
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

# Function to centre all galaxy particles based on stellar particle positions
.centre_galaxy = function(galaxy_data){
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
.new_half_mass_data = function(galaxy_data, p, q){
  # function for getting all particles within the half mass radius (ordered by ellipsoid radii)
  x = galaxy_data$x; y = galaxy_data$y; z = galaxy_data$z
  half_mass = sum(galaxy_data$Mass) / 2

  ellip_radius = sqrt((x*x) + ((y/p)*(y/p)) + ((z/q)*(z/q)))

  int_order = order(ellip_radius) # get the indicies of the radii in order (low to high)
  ordered_galaxy_data = galaxy_data[int_order,]
  cum_mass  = cumsum(ordered_galaxy_data$Mass) # cumulative sum of mass given this order
  half_mass_ind = which(abs(cum_mass - half_mass) == min(abs(cum_mass - half_mass))) # at what radius does this half-mass now occur?

  return(ordered_galaxy_data[1:(half_mass_ind-1),])
}

.ellipsoid_tensor = function(galaxy_data, p, q){
  # Computing the weighted ellipsoid tensor
  x = galaxy_data$x; y = galaxy_data$y; z = galaxy_data$z

  ellip_radius = sqrt((x*x) + ((y/p)*(y/p)) + ((z/q)*(z/q)))

  M = array(data = NA, dim = c(3,3))

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

  return(list("eigenvalues"= eig$values, "p" = p, "q" = q, "y_axis" = yax, "z_axis" = zax))
}

# Function to iteratively find the shape and align at the half-mass stellar radius
.measure_pqj = function(galaxy_data, abort_count=50){
  # Set up - we begin by assuming a sphere
  a = 1; b = 1; c = 1
  p = b/a; q = c/a
  aborted = 0; flag = 0
  cnt = 1
  temp_p = numeric(); temp_q = numeric()

  # Select all particles within initial half-mass (spherical) of stellar
  hm_galaxy_data = .new_half_mass_data(galaxy_data$star_part, p, q)

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

    hm_galaxy_data = .new_half_mass_data(galaxy_data$star_part, fit_ellip$p, fit_ellip$q)
    p = fit_ellip$p
    q = fit_ellip$q
    cnt = cnt + 1
  }

  return(list("galaxy_data" = galaxy_data, "p" = mean(tail(temp_p, n=6)), "q" = mean(tail(temp_q, n=6))))

}

# Function to align full galaxy based on the stellar particles
.align_galaxy = function(galaxy_data){
  data = .measure_pqj(galaxy_data)
  return(data$galaxy_data)
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

    part_spec = array(data = NA, dim = c(1, length(Template$Wave)))

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
      cl = snow::makeCluster(cores)
      doSNOW::registerDoSNOW(cl)
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

.sum_velocities = function(galaxy_sample, observation, cores){
  vel_diff = function(lum, vy_obs){diff((lum * pnorm(observation$vbin_edges, mean = vy_obs, sd = observation$vbin_error)))}

  if (cores > 1){
    cl = snow::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    i = integer()
    bins_list = foreach(i = 1:length(galaxy_sample$luminosity)) %dopar% { vel_diff(lum = galaxy_sample$luminosity[i], vy_obs = galaxy_sample$vy_obs[i]) }
    closeAllConnections()
    bins = matrix(unlist(bins_list, use.names=FALSE), nrow = observation$vbin)
  } else {
    bins = mapply(vel_diff, galaxy_sample$luminosity, galaxy_sample$vy_obs)
  }

  return(rowSums(bins))

}

.particles_to_pixels = function(galaxy_data, occupied, cores){
  if (cores > 1){
    cl = snow::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    x = integer()
    particle_IDs = foreach(x = occupied) %dopar% { which(galaxy_data$pixel_pos == x) }
    closeAllConnections()
  } else {
    particle_IDs = lapply(occupied, function(x) which(galaxy_data$pixel_pos == x))
  }
 return(particle_IDs)
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
