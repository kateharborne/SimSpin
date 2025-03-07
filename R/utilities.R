# Author: Kate Harborne
# Co-author: Alice Serene
# Date: 10/01/2023
# Title: Utilities functions (i.e. functions hidden from the user)

# Some useful constants:
.lsol_to_erg        = 3.828e33
.mpc_to_cm          = 3.08568e+24
.speed_of_light     = 299792.458
.cm_to_kpc          = 3.24078e-22
.cms_to_kms         = 1e-5
.g_to_msol          = 5.02785e-34
.gcm1_to_msolkm1    = 5.02785e-29
.gcm3_to_msolkm3    = 5.02785e-28
.gcm3_to_msolkpc3   = 1.477e+31
.g_constant_cgs     = 6.67430e-11
.g_in_kpcMsolkms2   = 4.3009e-6
.s_to_yr            = 3.171e-8
.mass_of_proton     = 1.67262e-24  # grams
.adiabatic_index    = 5/3          # heat is contained
.Boltzmann_constant = 1.38066e-16  # cm^2 g s^-2 K-1


# globalVariable definitions
globalVariables(c(".N", ":=", "Age", "Carbon", "CellSize", "Density", "filter_luminosity",
                  "Hydrogen", "hcl.colors", "ID", "Initial_Mass", "luminosity", "Mass",
                  "Metallicity", "N", "Oxygen", "par", "pixel_pos", "R", "SFR", "SFT", "SmoothingLength",
                  "sed_id", "Temperature", "ThermalDispersion", "text", "vx", "vy", "vz", "x", "y",
                  "z"))

# Functions for computing weighted means
.meanwt = function(x,wt){
  if (sum(wt,na.rm=T) == 0){
    val = 0
  } else {
    val = sum(x*wt, na.rm=T)/sum(wt,na.rm=T)
  }
  return(val)
} # weighted mean

.varwt = function(x, wt, xcen){
  if (sum(wt,na.rm=T) == 0){
    val = 0
  } else {
    if (missing(xcen)){xcen = .meanwt(x,wt)}
    val = sum(wt*(x - xcen)^2, na.rm=T)/sum(wt, na.rm=T)
  }
  return(val)
} # weighted variance

# A function for fitting a Gaussian Hermit distribution to the LOSVD
.losvd_fit_vsig= function(par, x, losvd){

  vel = par[1]
  sig = par[2]
  k  = 1

  w = (x - vel)/sig

  measured_vlos = (k * exp(-0.5*(w^2)))

  return=sum((measured_vlos-losvd)^2)
}

# A function for fitting a Gaussian Hermit distribution to the LOSVD
.losvd_fit_h3h4 = function(par, x, losvd){

  vel = par[1]
  sig = par[2]
  h3 = par[3]
  h4 = par[4]
  k  = 1

  w = (x - vel)/sig
  H3 = (1/sqrt(6))  * (((2*sqrt(2))* w^3) - ((3*sqrt(2)) * w))
  H4 = (1/sqrt(24)) * ((4* w^4) - (12 * w^2) + 3)

  measured_vlos = (k * exp(-0.5*(w^2))) * (1 + (h3*H3) + (h4*H4))

  return=sum((measured_vlos-losvd)^2)
}

.losvd_out_h3h4 = function(x, vel, sig, h3, h4){

  k=1
  w = (x - vel)/sig
  H3 = (1/sqrt(6))  * (((2*sqrt(2))* w^3) - ((3*sqrt(2)) * w))
  H4 = (1/sqrt(24)) * ((4* w^4) - (12 * w^2) + 3)

  measured_vlos = (k * exp(-0.5*(w^2))) * (1 + (h3*H3) + (h4*H4))
  return(measured_vlos)
}

.losvd_out_vsig = function(x, vel, sig){
  k=1
  w = (x - vel)/sig

  measured_vlos = (k * exp(-0.5*(w^2)))
  return(measured_vlos)
}

# A function for combining multiple results from a parallel loop
.comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Function to centre all galaxy particles based on stellar particle positions
.centre_galaxy = function(galaxy_data, centre=NA){

  if (!is.na(centre[1])){ # if an external centre is provided, use this to centre positions
    stellar_data = galaxy_data$star_part
    gas_data = galaxy_data$gas_part

    # check that provided centre falls within the range of values within the input file
    check_bounds = (centre[1] >= min(stellar_data$x) & centre[1] <= max(stellar_data$x) &
                    centre[2] >= min(stellar_data$y) & centre[2] <= max(stellar_data$y) &
                    centre[3] >= min(stellar_data$z) & centre[3] <= max(stellar_data$z))

    # transform coordinates relative to centre point
    if (check_bounds){
      stellar_data$x = stellar_data$x - centre[1]
      stellar_data$y = stellar_data$y - centre[2]
      stellar_data$z = stellar_data$z - centre[3]
      star_r2 = stellar_data$x^2 + stellar_data$y^2 + stellar_data$z^2
      star_vcen = c(median(stellar_data$vx[star_r2 < 100]), # using the median velocities within
                    median(stellar_data$vy[star_r2 < 100]), # 10kpc of the galaxy centre to define
                    median(stellar_data$vz[star_r2 < 100])) # the central velocity
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

    } else {
      stop(paste0("Error: Requested centre is outside the bounds of the region within the input file. \n
                  Please re-submit with centre specified, c(x,y,z): \n ",
                  min(stellar_data$x), " <= x <= ", max(stellar_data$x), "\n ",
                  min(stellar_data$y), " <= y <= ", max(stellar_data$y), "\n ",
                  min(stellar_data$z), " <= z <= ", max(stellar_data$z), "\n "))
    }
  }

  # what to do when no centre given
  # (first check data exists)
  else if (!is.null(galaxy_data$star_part)){
    stellar_data = cen_galaxy(galaxy_data$star_part) # centering and computing medians for stellar particles
    galaxy_data$star_part = data.table::as.data.table(stellar_data$part_data)

    if (!is.null(galaxy_data$gas_part)){ # if gas is also present, centering these particles based on stellar medians
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

  M[1,1] = sum((galaxy_data$Mass * x * x) / ellip_radius, na.rm = T)
  M[1,2] = sum((galaxy_data$Mass * x * y) / ellip_radius, na.rm = T)
  M[1,3] = sum((galaxy_data$Mass * x * z) / ellip_radius, na.rm = T)
  M[2,1] = M[1,2]
  M[2,2] = sum((galaxy_data$Mass * y * y) / ellip_radius, na.rm = T)
  M[2,3] = sum((galaxy_data$Mass * y * z) / ellip_radius, na.rm = T)
  M[3,1] = M[1,3]
  M[3,2] = M[2,3]
  M[3,3] = sum((galaxy_data$Mass * z * z) / ellip_radius, na.rm = T)

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

    if(!is.na(half_mass[1])){
      check_bounds = (sum(galaxy_data$gas_part$Mass) > half_mass) & # check that the total model contains more mass than the requested half-mass
        (min(galaxy_data$gas_part$Mass) <= half_mass) # check that the requested half mass is greater than a single particle

      if (!check_bounds){
        stop(paste0("Error: Requested half-mass is outside the range possible within the input simulation. \n",
                    "Minimum mass of particle = ", min(galaxy_data$gas_part$Mass), "\n",
                    "Maximum mass of simulation = ", sum(galaxy_data$gas_part$Mass), "\n",
                    "Requested half mass ", half_mass, " is outside this range. \n",
                    "Please resubmit your request with an appropriate value OR with NA for self-fit half-mass."))
      }
    }

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

    if(!is.na(half_mass[1])){
      check_bounds = (sum(galaxy_data$star_part$Mass) > half_mass) & # check that the total model contains more mass than the requested half-mass
                     (min(galaxy_data$star_part$Mass) < half_mass) # check that the requested half mass is greater than a single particle

      if (!check_bounds){
        stop(paste0("Error: Requested half-mass is outside the range possible within the input simulation. \n",
                    "Minimum mass of particle = ", min(galaxy_data$star_part$Mass), "\n",
                    "Maximum mass of simulation = ", sum(galaxy_data$star_part$Mass), "\n",
                    "Requested half mass ", half_mass, " is outside this range. \n",
                    "Please resubmit your request with an appropriate value OR with NA for self-fit half-mass."))
      }
    }

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

.cubic_spline_m4 = function(r){
  return((1/pi)*(((1/4) * ((2 - r)^3)) - ((1 - r)^3)))
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
  # method for calculating a randomly distribution of n points uniformly
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
  if (kernel == "M4"){
    weights = .cubic_spline_m4(sph[,3])
    sph_kernel = data.frame("long" = sph[,1], "lat" = sph[,2],
                            "r/h" = sph[,3], "weight" = weights/sum(weights))

  }
  return(sph_kernel)
}

# interp_quick function from https://github.com/asgr/ProSpect/blob/master/R/utility.R
.interp_quick = function(x, params, log=FALSE){
  if(length(x) > 1){stop('x must be scalar!')}
  if(x < min(params)){
    return(c(ID_lo = 1, ID_hi = 1, wt_lo = 1, wt_hi = 0))
  }
  if(x > max(params)){
    return(c(ID_lo = length(params), ID_hi = length(params), wt_lo = 0, wt_hi = 1))
  }
  if(log){
    params = log(params)
    x = log(x)
  }
  interp = approx(params, 1:length(params), xout=x)$y
  IDlo = floor(interp)
  IDhi = ceiling(interp)
  return(c(ID_lo = IDlo, ID_hi = IDhi, wt_lo = 1-(interp-IDlo), wt_hi = interp-IDlo))
}

.spectral_weights = function(Metallicity, Age, Template, cores){

  f = function(metallicity, age) {
      Z = as.numeric(.interp_quick(metallicity, Template$Z, log = TRUE))
      A = as.numeric(.interp_quick(age * 1e9, Template$Age, log = TRUE))
      # ID_lo = 1, ID_hi = 2, wt_lo = 3, wt_hi = 4
      return(c(Z, A))
    }

    if (length(Age) == 1){

      output = f(Metallicity, Age)
      # "Z_ID_lo", "Z_ID_hi", "Z_wt_lo", "Z_wt_hi"
      # "A_ID_lo", "A_ID_hi", "A_wt_lo", "A_wt_hi"
      output = data.table::as.data.table(output)
      return(output)

    } else {

      if (cores > 1) {
        doParallel::registerDoParallel(cores = cores)
        i = integer()
        output = foreach::foreach(i = 1:length(Metallicity), .packages = c("SimSpin", "foreach")) %dopar% {
          f(Metallicity[i], Age[i])
        }
        closeAllConnections()
        output = data.table::as.data.table(output)
      } else {
        output = mapply(f, Metallicity, Age)
        output = data.table::as.data.table(output)
      }

      return(output)

    }
}

# Function to generate spectra (w/o mass weighting)
.spectra = function(SW, Template){
  # s = function(SW){
    weights = data.frame("hihi" = SW[4] * SW[8],   # Z_ID_lo = 1, Z_ID_hi = 2, Z_wt_lo = 3, Z_wt_hi = 4
                         "hilo" = SW[4] * SW[7],   # A_ID_lo = 5, A_ID_hi = 6, A_wt_lo = 7, A_wt_hi = 8
                         "lohi" = SW[3] * SW[8],
                         "lolo" = SW[3] * SW[7])

    part_spec = array(data = 0.0, dim = c(1, length(Template$Wave)))
    part_spec = ((Template$Zspec[[SW[2]]][SW[6],] * weights$hihi) +
                 (Template$Zspec[[SW[2]]][SW[5],] * weights$hilo) +
                 (Template$Zspec[[SW[1]]][SW[6],] * weights$lohi) +
                 (Template$Zspec[[SW[1]]][SW[5],] * weights$lolo))

    return(part_spec)
  # }

  # if (length(spectral_weights) == 1){
  #   spectra = data.table::data.table(s(spectral_weights))
  #   return(spectra)
  # } else {
  #
  #   if (cores > 1) {
  #     doParallel::registerDoParallel(cores = cores)
  #     i = integer()
  #     part_spec = foreach::foreach(i = 1:length(spectral_weights), .packages = c("SimSpin", "foreach")) %dopar% {
  #       s(spectral_weights[[i]])
  #     }
  #     closeAllConnections()
  #     output = data.table::as.data.table(part_spec)
  #   } else {
  #     output = mapply(s, spectral_weights)
  #     output = data.table::as.data.table(output)
  #   }
  #   return(output)
  #
  # }
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

.sum_velocities = function(galaxy_sample, observation){
  vel_diff = function(lum, vy){diff((lum * pnorm(observation$vbin_edges, mean = vy, sd = 0)))}

  bins = mapply(vel_diff, galaxy_sample$luminosity, galaxy_sample$vy)

  return(rowSums(bins))

}

.sum_gas_velocities = function(galaxy_sample, observation){
  vel_diff = function(mass, vy, sigma_t){diff((mass * pnorm(observation$vbin_edges, mean = vy, sd = sigma_t)))}

  bins = mapply(vel_diff, galaxy_sample$Mass, galaxy_sample$vy, (galaxy_sample$ThermalDispersion/sqrt(3)))

  return(rowSums(bins))

}

# Function to apply LSF to spectra
.lsf_convolution = function(observation, luminosity, lsf_sigma){

  #kernel_radius = (4 * lsf_sigma + 0.5)
  #x = seq(-kernel_radius, kernel_radius, by=observation$wave_res)
  x = seq(-(observation$wave_res*12), (observation$wave_res*12), by=observation$wave_res)
  #x = seq(-kernel_radius, kernel_radius, length.out = 25)
  #phi_x = exp((-0.5 / (lsf_sigma^2)) * (x^2)) / (lsf_sigma * sqrt(2*pi))
  phi_x = exp((-0.5 * (x^2)) / (lsf_sigma^2))
  phi_x = phi_x / sum(phi_x)

  lum = stats::convolve(luminosity, phi_x, type="open")
  end = (length(luminosity) + length(phi_x) - 1) - 12

  return(lum[13:end])
}

# Function to add noise
.add_noise = function(cube, S2N){
  S2N[is.infinite(S2N)] = 0 # removing infinite noise where particles per pixel = 0
  noisey_cube = cube
  for (i in 1:nrow(S2N)){
    for (j in 1:ncol(S2N)){
      if (S2N[i,j]!=0){
        noise = (stats::rnorm(length(cube[i,j,]), mean = 0, sd=1))*S2N[i,j]
        noisey_cube[i,j,] = cube[i,j,]*noise # to be added to the cube
      }
    }
  }
  return(noisey_cube)
}

# Function to generate a Gaussian kernel
.gaussian_kernel = function(m, n, sigma){
  dim = pracma::meshgrid(-((m-1)/2):((m-1)/2), -((n-1)/2):((n-1)/2))
  hg = exp(-(dim$X^2 + dim$Y^2)/(2*sigma^2))
  kernel = hg / sum(hg)
  return(kernel)
}

# Functions for computing spaxel properties

# Taken from https://github.com/asgr/ProSpect/blob/d340c64555ba631257513ea4c99b0069cdebf477/R/utility.R#L123
# to avoid ProSpect dependency

.qdiff=function(vec){
    return(c(0,vec[2:length(vec)]-vec[1:(length(vec)-1)]))
 }

# Taken from https://github.com/asgr/ProSpect/blob/d340c64555ba631257513ea4c99b0069cdebf477/R/photom.R#L281
# to avoid ProSpect dependency and trimmed for the purpose of these internal functions

.bandpass=function(wave, flux, filter){

  response = filter(wave)
  response[is.na(response)] = 0

  wave_diff=abs(.qdiff(wave))

  if (is.null(dim(flux))){
    output = response * wave * flux * wave_diff/sum(response * wave * wave_diff, na.rm = TRUE)
    return(sum(output, na.rm=TRUE))
  } else {
    for (j in 1:dim(flux)[2]){
    set(flux, j = j,
        value = response * wave * flux[[j]] * wave_diff/sum(response * wave * wave_diff, na.rm = TRUE))
    }
    return(as.numeric(colSums(flux, na.rm=TRUE)))
  }

}

.compute_flux = function(observation, galaxy_data, simspin_data,
                         template, verbose, spectra_flag){

  # Function to pre-compute the fluxes per particle and add this to the galaxy
  # data frame. This can then be used quickly to compute an un-binned flux map.
  #
  # observation  :: List. The output list from the observation() function.
  # galaxy_data  :: Data.Frame. The output table of particle details.
  # simspin_data :: List. The simspin_file data, importantly containing the
  #                 spectral mapping for each particle.
  # template     :: List. The spectral templates with which the SimSpin file has
  #                 been built.
  # verbose      :: Boolean. Should the progress of the function be output to
  #                 the user?
  # spectra_flag :: Boolean. Is this an old SimSpin file where the full spectra
  #                 is given, rather than just the template indicies?
  #

  # read original wavelengths of the template spectra and then applying a shift
  # to those spectra due to redshift, z
  if(verbose){cat("Using assigned spectra to compute the flux per particle... \n")}

  lum = numeric(length = nrow(galaxy_data))
  band_lum = numeric(length = nrow(galaxy_data))

  wavelength = template$Wave * (observation$z + 1)
  wave_diff_observed  = .qdiff(observation$wave_seq)
  filter = stats::approxfun(x = observation$filter$wave, y = abs(observation$filter$response))

  for (p in 1:nrow(galaxy_data)){

    # reading particle luminosity in units of Lsol/Ang
    if (spectra_flag == 1){
      intrinsic_spectra = simspin_data$spectra[[galaxy_data$sed_id[p]]][1:length(template$Wave)] *
        (galaxy_data$Initial_Mass[p])

    } else {
      intrinsic_spectra = .spectra(simspin_data$spectral_weights[[galaxy_data$sed_id[p]]],
                                   template) * (galaxy_data$Initial_Mass[p])
    }

    # pulling wavelengths and using doppler formula to compute the shift in
    #   wavelengths caused by LOS velocity
    wave_shift = wavelength * exp((galaxy_data$vy[p] / .speed_of_light))

    # pulling out the wavelengths that would fall within the telescope range
    wave_seq_int = which(wave_shift >= min(observation$wave_edges) & wave_shift <= max(observation$wave_edges))
    wave_diff_intrinsic = .qdiff(wave_shift[wave_seq_int])

    tot_lum = sum(intrinsic_spectra[wave_seq_int] * wave_diff_intrinsic)
    # total luminosity within the wavelength range of the telescope

    # interpolate each shifted wavelength to telescope grid of wavelengths
    #   and sum to one spectra
    part_lum = stats::spline(x = wave_shift, y = intrinsic_spectra, xout = observation$wave_seq)[[2]]

    new_lum = sum(part_lum * wave_diff_observed)
    # total luminosity integrated in wavelength bins

    scale_frac = tot_lum / new_lum
    # scaling factor necessary to conserve flux in the new spectrum

    luminosity = part_lum*scale_frac

    # transform luminosity into flux detected at telescope
    #    flux in units erg/s/cm^2/Ang
    spectral_dist = (luminosity*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) /
      (1 + observation$z)


    lum[p] = sum(spectral_dist, na.rm=T)
    band_lum[p] = .bandpass(wave = observation$wave_seq,
                            flux = spectral_dist,
                            filter = filter)

    if(verbose){if(p == 1){cat("Computed flux from spectra 1, ")}else{cat(paste(p), ", ")}}

  }

  galaxy_data[ , luminosity := lum, ]
  galaxy_data[ , filter_luminosity := band_lum, ]

  if (verbose){cat("\n Done!")}

  return(galaxy_data)

}

.compute_flux_mc = function(observation, galaxy_data, simspin_data,
                             template, verbose, spectra_flag, cores){

  # Function to pre-compute the fluxes per particle and add this to the galaxy
  # data frame. This can then be used quickly to compute an un-binned flux map.
  # This iteration of the function is Parallelised using doParallel.
  #
  # observation  :: List. The output list from the observation() function.
  # galaxy_data  :: Data.Frame. The output table of particle details.
  # simspin_data :: List. The simspin_file data, importantly containing the
  #                 spectral mapping for each particle.
  # template     :: List. The spectral templates with which the SimSpin file has
  #                 been built.
  # verbose      :: Boolean. Should the progress of the function be output to
  #                 the user?
  # spectra_flag :: Boolean. Is this an old SimSpin file where the full spectra
  #                 is given, rather than just the template indicies?
  # cores        :: Integer. The number of cores between which to divide the
  #                 task.

  # read original wavelengths of the template spectra and then applying a shift
  # to those spectra due to redshift, z

  lum = numeric(length = nrow(galaxy_data))
  band_lum = numeric(length = nrow(galaxy_data))

  wavelength = template$Wave * (observation$z + 1)
  wave_diff_observed  = .qdiff(observation$wave_seq)
  filter = stats::approxfun(x = observation$filter$wave, y = abs(observation$filter$response))

  doParallel::registerDoParallel(cores)

  p = integer()
  output = foreach(p = 1:nrow(galaxy_data), .combine='.comb', .multicombine=TRUE,
                   .init=list(list(),list())) %dopar% {

    # reading particle luminosity in units of Lsol/Ang
    if (spectra_flag == 1){
      intrinsic_spectra = simspin_data$spectra[[galaxy_data$sed_id[p]]][1:length(template$Wave)] *
        (galaxy_data$Initial_Mass[p])

    } else {
      intrinsic_spectra = .spectra(simspin_data$spectral_weights[[galaxy_data$sed_id[p]]],
                                   template) * (galaxy_data$Initial_Mass[p])
    }

    # pulling wavelengths and using doppler formula to compute the shift in
    #   wavelengths caused by LOS velocity
    wave_shift = wavelength * exp((galaxy_data$vy[p] / .speed_of_light))

    # pulling out the wavelengths that would fall within the telescope range
    wave_seq_int = which(wave_shift >= min(observation$wave_edges) & wave_shift <= max(observation$wave_edges))
    wave_diff_intrinsic = .qdiff(wave_shift[wave_seq_int])

    tot_lum = sum(intrinsic_spectra[wave_seq_int] * wave_diff_intrinsic)
    # total luminosity within the wavelength range of the telescope

    # interpolate each shifted wavelength to telescope grid of wavelengths
    #   and sum to one spectra
    part_lum = stats::spline(x = wave_shift, y = intrinsic_spectra, xout = observation$wave_seq)[[2]]

    new_lum = sum(part_lum * wave_diff_observed)
    # total luminosity integrated in wavelength bins

    scale_frac = tot_lum / new_lum
    # scaling factor necessary to conserve flux in the new spectrum

    luminosity = part_lum*scale_frac

    # transform luminosity into flux detected at telescope
    #    flux in units erg/s/cm^2/Ang
    spectral_dist = (luminosity*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) /
      (1 + observation$z)


    lum = sum(spectral_dist, na.rm=T)
    band_lum = .bandpass(wave = observation$wave_seq,
                            flux = spectral_dist,
                            filter = filter)

    result = list(lum,
                  band_lum)
    return(result)
    closeAllConnections()

  }

  lum = matrix(unlist(output[[1]]))
  band_lum = matrix(unlist(output[[2]]))

  galaxy_data[ , luminosity := lum, ]
  galaxy_data[ , filter_luminosity := band_lum, ]

  return(galaxy_data)

}


# spectral mode -
.spectral_spaxels = function(part_in_spaxel, wavelength, observation,
                             galaxy_data, simspin_data, template, verbose, spectra_flag){

  spectra = matrix(data = 0.0, ncol = observation$wave_bin, nrow = observation$sbin^2)
  vel_los = array(data = 0.0, dim = observation$sbin^2)
  dis_los = array(data = 0.0, dim = observation$sbin^2)
  ageM_map = array(data = 0.0, dim = observation$sbin^2)
  ageL_map = array(data = 0.0, dim = observation$sbin^2)
  met_map = array(data = 0.0, dim = observation$sbin^2)
  vorbin_map = array(data=0, dim = observation$sbin^2)
  filter = stats::approxfun(x = observation$filter$wave, y = abs(observation$filter$response))

  for (i in 1:nrow(part_in_spaxel)){ # computing the spectra at each occupied spatial pixel position

    num_part = part_in_spaxel$N[i] # number of particles in spaxel
    vorbin_map[part_in_spaxel$pixel_pos[[i]]] = i

    galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
    pixel_size = length(unique(galaxy_sample$pixel_pos))
    galaxy_sample$Mass = galaxy_sample$Mass / pixel_size

    luminosity = numeric(length(observation$wave_seq)) # initiallise a spectrum array for this pixel

    wave_diff_observed  = .qdiff(observation$wave_seq) # compute the wavelength channel widths in telescope

    for (p in 1:num_part){

      if (spectra_flag == 1){
        intrinsic_spectra = simspin_data$spectra[[galaxy_sample$sed_id[p]]][1:length(template$Wave)] *
          (galaxy_sample$Initial_Mass[p])

      } else {
        intrinsic_spectra = .spectra(simspin_data$spectral_weights[[galaxy_sample$sed_id[p]]],
                                     template) * (galaxy_sample$Initial_Mass[p])
      }
      # reading particle luminosity in units of Lsol/Ang

      # pulling wavelengths and using doppler formula to compute the shift in
      #   wavelengths caused by LOS velocity
      wave_shift = wavelength * exp((galaxy_sample$vy[p] / .speed_of_light))

      # pulling out the wavelengths that would fall within the telescope range
      wave_seq_int = which(wave_shift >= min(observation$wave_edges) & wave_shift <= max(observation$wave_edges))
      wave_diff_intrinsic = .qdiff(wave_shift[wave_seq_int])

      tot_lum = sum(intrinsic_spectra[wave_seq_int] * wave_diff_intrinsic)
      # total luminosity within the wavelength range of the telescope

      # interpolate each shifted wavelength to telescope grid of wavelengths
      #   and sum to one spectra
      part_lum = stats::spline(x = wave_shift, y = intrinsic_spectra, xout = observation$wave_seq)[[2]]

      new_lum = sum(part_lum * wave_diff_observed)
      # total luminosity integrated in wavelength bins

      scale_frac = tot_lum / new_lum
      # scaling factor necessary to conserve flux in the new spectrum

      luminosity = luminosity + (part_lum*scale_frac)
    }

    # transform luminosity into flux detected at telescope
    #    flux in units erg/s/cm^2/Ang
    spectral_dist = (luminosity*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) /
                    (1 + observation$z) / pixel_size

    for (bin in 1:length(part_in_spaxel$pixel_pos[[i]])){
      spectra[part_in_spaxel$pixel_pos[[i]][bin], ] = spectral_dist
    }

    vel_los[part_in_spaxel$pixel_pos[[i]]] = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
    dis_los[part_in_spaxel$pixel_pos[[i]]] = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
    ageM_map[part_in_spaxel$pixel_pos[[i]]] = .meanwt(galaxy_sample$Age, galaxy_sample$Mass)
    ageL_map[part_in_spaxel$pixel_pos[[i]]] = .meanwt(galaxy_sample$Age, galaxy_sample$Luminosity)
    met_map[part_in_spaxel$pixel_pos[[i]]] = .meanwt(galaxy_sample$Metallicity, galaxy_sample$Mass)

    if (verbose){cat(i, "... ", sep = "")}
  }
  return(list(spectra,
              vel_los, dis_los, ageM_map, ageL_map, met_map,
              vorbin_map))
}

.spectral_spaxels_mc = function(part_in_spaxel, wavelength, observation, galaxy_data, simspin_data, template, verbose, cores, spectra_flag){

  spectra = matrix(data = 0.0, ncol = observation$wave_bin, nrow = observation$sbin^2)
  vel_los = array(data = 0.0, dim = observation$sbin^2)
  dis_los = array(data = 0.0, dim = observation$sbin^2)
  ageM_map = array(data = 0.0, dim = observation$sbin^2)
  ageL_map = array(data = 0.0, dim = observation$sbin^2)
  met_map = array(data = 0.0, dim = observation$sbin^2)
  vorbin_map = array(data=0, dim = observation$sbin^2)
  filter = stats::approxfun(x = observation$filter$wave, y = abs(observation$filter$response))

  doParallel::registerDoParallel(cores)

  i = integer()
  output = foreach(i = 1:(dim(part_in_spaxel)[1]), .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list(), list(), list())) %dopar% {

                     num_part = part_in_spaxel$N[i]
                     vorbin_map = i

                     galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
                     pixel_size = length(unique(galaxy_sample$pixel_pos))
                     galaxy_sample$Mass = galaxy_sample$Mass / pixel_size

                     luminosity = numeric(length(observation$wave_seq)) # initialize a spectrum array for this pixel

                     wave_diff_observed  = .qdiff(observation$wave_seq) # compute the wavelength channel widths in telescope

                     for (p in 1:num_part){

                       if (spectra_flag == 1){
                         intrinsic_spectra = simspin_data$spectra[[galaxy_sample$sed_id[p]]][1:length(template$Wave)] *
                           (galaxy_sample$Initial_Mass[p])

                       } else {
                         intrinsic_spectra = .spectra(simspin_data$spectral_weights[[galaxy_sample$sed_id[p]]],
                                                      template) * (galaxy_sample$Initial_Mass[p])
                       }
                       # reading particle luminosity in units of Lsol/Ang

                       # pulling wavelengths and using doppler formula to compute the shift in
                       #   wavelengths caused by LOS velocity
                       wave_shift = wavelength * exp((galaxy_sample$vy[p] / .speed_of_light))#((galaxy_sample$vy[p] / .speed_of_light) * wavelength) + wavelength

                       # pulling out the wavelengths that would fall within the telescope range
                       wave_seq_int = which(wave_shift >= min(observation$wave_seq) & wave_shift <= max(observation$wave_seq))
                       wave_diff_intrinsic = .qdiff(wave_shift[wave_seq_int])

                       tot_lum = sum(intrinsic_spectra[wave_seq_int] * wave_diff_intrinsic)
                       # total luminosity within the wavelength range of the telescope

                       # interpolate each shifted wavelength to telescope grid of wavelengths
                       #   and sum to one spectra
                       part_lum = stats::spline(x = wave_shift, y = intrinsic_spectra, xout = observation$wave_seq)[[2]]

                       new_lum = sum(part_lum * wave_diff_observed)
                       # total luminosity integrated in wavelength bins

                       scale_frac = tot_lum / new_lum
                       # scaling factor necessary to conserve flux in the new spectrum

                       luminosity = luminosity + (part_lum*scale_frac)

                     }

                     # transform luminosity into flux detected at telescope
                     #    flux in units erg/s/cm^2/Ang
                     spectra = (luminosity*.lsol_to_erg) / (4 * pi * (observation$lum_dist*.mpc_to_cm)^2) /
                               (1 + observation$z) / pixel_size
                     vel_los = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
                     dis_los = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
                     ageM_map = .meanwt(galaxy_sample$Age, galaxy_sample$Mass)
                     ageL_map = .meanwt(galaxy_sample$Age, galaxy_sample$Luminosity)
                     met_map = .meanwt(galaxy_sample$Metallicity, galaxy_sample$Mass)

                     result = list(spectra,
                                   vel_los, dis_los, ageM_map, ageL_map, met_map,
                                   vorbin_map)
                     return(result)
                     closeAllConnections()
                   }

  spectral_dist = matrix(unlist(output[[1]]), ncol=observation$wave_bin, byrow = T)
  vel_dist = matrix(unlist(output[[2]]))
  dis_dist = matrix(unlist(output[[3]]))
  ageM_dist = matrix(unlist(output[[4]]))
  ageL_dist = matrix(unlist(output[[5]]))
  met_dist = matrix(unlist(output[[6]]))
  vorbin_dist = matrix(unlist(output[[7]]))

  for (bin in 1:nrow(part_in_spaxel)){

    for (p in 1:length(unlist(part_in_spaxel$pixel_pos[[bin]]))){
      spectra[unlist(part_in_spaxel$pixel_pos[[bin]][p]),] = spectral_dist[bin,]
    }

    vel_los[unlist(part_in_spaxel$pixel_pos[[bin]])] = vel_dist[bin]
    dis_los[unlist(part_in_spaxel$pixel_pos[[bin]])] = dis_dist[bin]
    ageM_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = ageM_dist[bin]
    ageL_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = ageL_dist[bin]
    met_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = met_dist[bin]
    vorbin_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = vorbin_dist[bin]
  }

  return(list(spectra,
              vel_los, dis_los, ageM_map, ageL_map, met_map,
              vorbin_map))
}

# stellar velocity mode -
.velocity_spaxels = function(part_in_spaxel, observation, galaxy_data, simspin_data, template, verbose, mass_flag, spectra_flag){

  vel_spec = matrix(data = 0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  ageM_map  = array(data = 0.0, dim = observation$sbin^2)
  ageL_map = array(data = 0.0, dim = observation$sbin^2)
  met_map  = array(data = 0.0, dim = observation$sbin^2)
  vorbin_map = array(data=0, dim = observation$sbin^2)

  for (i in 1:(dim(part_in_spaxel)[1])){

    vorbin_map[part_in_spaxel$pixel_pos[[i]]] = i

    galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
    pixel_size = length(unique(galaxy_sample$pixel_pos))
    galaxy_sample$Mass = galaxy_sample$Mass / pixel_size

    if (mass_flag){

      galaxy_sample[, luminosity := Mass, ]
      galaxy_sample[, filter_luminosity := Mass, ]

    }

    # adding the "gaussians" of each particle to the velocity bins
    v_dist = .sum_velocities(galaxy_sample = galaxy_sample, observation = observation)

    for (bin in 1:length(part_in_spaxel$pixel_pos[[i]])){
      vel_spec[part_in_spaxel$pixel_pos[[i]][bin], ] = v_dist
    }

    vel_los[part_in_spaxel$pixel_pos[[i]]]   = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
    dis_los[part_in_spaxel$pixel_pos[[i]]]   = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
    ageM_map[part_in_spaxel$pixel_pos[[i]]]   = .meanwt(galaxy_sample$Age, galaxy_sample$Mass)
    ageL_map[part_in_spaxel$pixel_pos[[i]]]   = .meanwt(galaxy_sample$Age, galaxy_sample$Luminosity)
    met_map[part_in_spaxel$pixel_pos[[i]]]   = .meanwt(galaxy_sample$Metallicity, galaxy_sample$Mass)

    if (verbose){cat(i, "... ", sep = "")}

  }

  return(list(vel_spec, #lum_map, band_map,
              vel_los, dis_los, ageM_map, ageL_map, met_map, #mass_map, part_map,
              vorbin_map))
}

.velocity_spaxels_mc = function(part_in_spaxel, observation, galaxy_data, simspin_data, template, verbose, cores, mass_flag, spectra_flag){

  vel_spec = matrix(data = 0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  ageM_map  = array(data = 0.0, dim = observation$sbin^2)
  ageL_map  = array(data = 0.0, dim = observation$sbin^2)
  met_map  = array(data = 0.0, dim = observation$sbin^2)
  vorbin_map = array(data=0, dim = observation$sbin^2)

  doParallel::registerDoParallel(cores)

  i = integer()
  output = foreach(i = 1:(dim(part_in_spaxel)[1]), .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list(), list(), list())) %dopar% {

                     vorbin_map = i

                     galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
                     pixel_size = length(unique(galaxy_sample$pixel_pos))
                     galaxy_sample$Mass = galaxy_sample$Mass / pixel_size

                     if (mass_flag){

                       galaxy_sample[, luminosity := Mass, ]
                       galaxy_sample[, filter_luminosity := Mass, ]

                     }

                     # adding the "gaussians" of each particle to the velocity bins
                     vel_spec = .sum_velocities(galaxy_sample = galaxy_sample, observation = observation)
                     vel_los = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
                     dis_los = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
                     ageM_map = .meanwt(galaxy_sample$Age, galaxy_sample$Mass)
                     ageL_map = .meanwt(galaxy_sample$Age, galaxy_sample$Luminosity)
                     met_map = .meanwt(galaxy_sample$Metallicity, galaxy_sample$Mass)

                     result = list(vel_spec,
                                   vel_los, dis_los, ageM_map, ageL_map, met_map,
                                   vorbin_map)
                     return(result)
                     closeAllConnections()
                   }


  v_dist = matrix(unlist(output[[1]]), ncol=observation$vbin, byrow = T)
  vel_dist = matrix(unlist(output[[2]]))
  dis_dist = matrix(unlist(output[[3]]))
  ageM_dist = matrix(unlist(output[[4]]))
  ageL_dist = matrix(unlist(output[[5]]))
  met_dist = matrix(unlist(output[[6]]))
  vorbin_dist = matrix(unlist(output[[7]]))

  for (bin in 1:nrow(part_in_spaxel)){
    for (p in 1:length(unlist(part_in_spaxel$pixel_pos[[bin]]))){
      vel_spec[unlist(part_in_spaxel$pixel_pos[[bin]][p]),] = v_dist[bin,]
    }
    vel_los[unlist(part_in_spaxel$pixel_pos[[bin]])] = vel_dist[bin]
    dis_los[unlist(part_in_spaxel$pixel_pos[[bin]])] = dis_dist[bin]
    ageM_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = ageM_dist[bin]
    ageL_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = ageL_dist[bin]
    met_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = met_dist[bin]
    vorbin_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = vorbin_dist[bin]

  }

  return(list(vel_spec,
              vel_los, dis_los, ageM_map, ageL_map, met_map,
              vorbin_map))

}

# gas velocity mode -
.gas_velocity_spaxels = function(part_in_spaxel, observation, galaxy_data, simspin_data, verbose){

  vel_spec = matrix(data = 0.0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  Z_map    = array(data = 0.0, dim = observation$sbin^2)
  OH_map   = array(data = 0.0, dim = observation$sbin^2)
  vorbin_map = array(data=0, dim = observation$sbin^2)

  for (i in 1:(dim(part_in_spaxel)[1])){

    vorbin_map[part_in_spaxel$pixel_pos[[i]]] = i

    galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
    pixel_size = length(unique(galaxy_sample$pixel_pos))
    galaxy_sample$Mass = galaxy_sample$Mass / pixel_size

    # adding the "gaussians" of each particle to the velocity bins
    v_dist = .sum_gas_velocities(galaxy_sample = galaxy_sample, observation = observation)

    for (bin in 1:length(part_in_spaxel$pixel_pos[[i]])){
      vel_spec[part_in_spaxel$pixel_pos[[i]][bin], ] = v_dist
    }

    vel_los[part_in_spaxel$pixel_pos[[i]]]   = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
    dis_los[part_in_spaxel$pixel_pos[[i]]]   = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
    Z_map[part_in_spaxel$pixel_pos[[i]]]     = log10(mean(galaxy_sample$Metallicity)/0.0127)
    OH_map[part_in_spaxel$pixel_pos[[i]]]    = log10(mean(galaxy_sample$Oxygen/galaxy_sample$Hydrogen))+12

    if (verbose){cat(i, "... ", sep = "")}

  }

  return(list(vel_spec,
              vel_los, dis_los,
              Z_map, OH_map,
              vorbin_map))
}

.gas_velocity_spaxels_mc = function(part_in_spaxel, observation, galaxy_data, simspin_data, verbose, cores){

  vel_spec = matrix(data = 0.0, ncol = observation$vbin, nrow = observation$sbin^2)
  vel_los  = array(data = 0.0, dim = observation$sbin^2)
  dis_los  = array(data = 0.0, dim = observation$sbin^2)
  Z_map    = array(data = 0.0, dim = observation$sbin^2)
  OH_map   = array(data = 0.0, dim = observation$sbin^2)
  vorbin_map = array(data=0, dim = observation$sbin^2)

  doParallel::registerDoParallel(cores)

  i = integer()
  output = foreach(i = 1:(dim(part_in_spaxel)[1]), .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list(), list())) %dopar% {

                     vorbin_map = i

                     galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
                     pixel_size = length(unique(galaxy_sample$pixel_pos))
                     galaxy_sample$Mass = galaxy_sample$Mass / pixel_size

                     # adding the "gaussians" of each particle to the velocity bins
                     vel_spec = .sum_gas_velocities(galaxy_sample = galaxy_sample, observation = observation)
                     vel_los  = .meanwt(galaxy_sample$vy, galaxy_sample$Mass)
                     dis_los  = sqrt(.varwt(galaxy_sample$vy, galaxy_sample$Mass))
                     Z_map    = log10(mean(galaxy_sample$Metallicity)/0.0127)
                     OH_map   = log10(mean(galaxy_sample$Oxygen/galaxy_sample$Hydrogen))+12

                     result = list(vel_spec,
                                   vel_los, dis_los,
                                   Z_map, OH_map,
                                   vorbin_map)
                     return(result)
                     closeAllConnections()
                   }

  v_dist = matrix(unlist(output[[1]]), ncol=observation$vbin, byrow = T)
  vel_dist = matrix(unlist(output[[2]]))
  dis_dist = matrix(unlist(output[[3]]))
  Z_dist = matrix(unlist(output[[4]]))
  OH_dist = matrix(unlist(output[[5]]))
  vorbin_dist = matrix(unlist(output[[6]]))

  for (bin in 1:nrow(part_in_spaxel)){
    for (p in 1:length(unlist(part_in_spaxel$pixel_pos[[bin]]))){
      vel_spec[unlist(part_in_spaxel$pixel_pos[[bin]][p]),] = v_dist[bin,]
    }
    vel_los[unlist(part_in_spaxel$pixel_pos[[bin]])] = vel_dist[bin]
    dis_los[unlist(part_in_spaxel$pixel_pos[[bin]])] = dis_dist[bin]
    Z_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = Z_dist[bin]
    OH_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = OH_dist[bin]
    vorbin_map[unlist(part_in_spaxel$pixel_pos[[bin]])] = vorbin_dist[bin]
  }

  return(list(vel_spec,
              vel_los, dis_los,
              Z_map, OH_map,
              vorbin_map))

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

      new_gas_part[ind1:ind2, x := part$x+new_xyz[,1],]
      new_gas_part[ind1:ind2, y := part$y+new_xyz[,2],]
      new_gas_part[ind1:ind2, z := part$z+new_xyz[,3],]
      new_gas_part[ind1:ind2, Mass := part$Mass*rand_pos$weight,]
      new_gas_part[ind1:ind2, SFR := part$SFR*rand_pos$weight,]
      new_gas_part[ind1:ind2, Density := part$Density*rand_pos$weight,]
      new_gas_part[ind1:ind2, Temperature := part$Temperature*rand_pos$weight,]
      new_gas_part[ind1:ind2, ThermalDispersion := sqrt((part$ThermalDispersion^2)*rand_pos$weight),]
      # in the new data.frame of particle properties, assign their
      # new positions and masses scaled by the kernel weight.
    }

  return(new_gas_part)
}

# same as above but with multiple cores
.sph_spawn_mc = function(gas_part, new_gas_part, sph_spawn_n, kernel, cores){

  # Parallel version of the function ".sph_spawn" above.
  doParallel::registerDoParallel(cores)
  no_gas = length(gas_part$ID)
  x = numeric(no_gas); y = numeric(no_gas)
  z = numeric(no_gas); Mass = numeric(no_gas) # initialising

  i = integer()

  output = foreach(i = 1:no_gas, .combine='.comb', .multicombine=TRUE,
                   .init=list(list(), list(), list(), list(), list(), list(), list(), list())) %dopar% {

                       part = gas_part[i,]

                       rand_pos = .generate_uniform_sphere(sph_spawn_n, kernel = kernel)
                       rand_pos$r = rand_pos$r.h * part$SmoothingLength
                       new_xyz = sphereplot::sph2car(cbind(rand_pos$long, rand_pos$lat, rand_pos$r))

                       x = part$x+new_xyz[,1]
                       y = part$y+new_xyz[,2]
                       z = part$z+new_xyz[,3]
                       Mass = part$Mass*rand_pos$weight
                       SFR  = part$SFR*rand_pos$weight
                       Density = part$Density*rand_pos$weight
                       Temperature = part$Temperature*rand_pos$weight
                       ThermalDispersion = sqrt((part$ThermalDispersion^2)*rand_pos$weight)

                       return(list(x, y, z, Mass, SFR, Density, Temperature, ThermalDispersion))
                       closeAllConnections()
                     }


  new_gas_part[, x := as.numeric(unlist(output[[1]])),]
  new_gas_part[, y := as.numeric(unlist(output[[2]])),]
  new_gas_part[, z := as.numeric(unlist(output[[3]])),]
  new_gas_part[, Mass := as.numeric(unlist(output[[4]])),]
  new_gas_part[, SFR := as.numeric(unlist(output[[5]])),]
  new_gas_part[, Density := as.numeric(unlist(output[[6]])),]
  new_gas_part[, Temperature := as.numeric(unlist(output[[7]])),]
  new_gas_part[, ThermalDispersion := as.numeric(unlist(output[[8]])),]

  return(new_gas_part)
}
