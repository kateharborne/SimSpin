# Kate Harborne (last edit - 13/09/2017)
# Utilities functions

.meanwt=function(x, wt){
  return=sum(x*wt, na.rm = T)/sum(wt, na.rm = T)
} # weighted mean

.varwt=function(x, wt, xcen){
  if(missing(xcen)){xcen = .meanwt(x, wt)}
  return=(sum(wt * (x - xcen)^2, na.rm = T) / (sum(wt, na.rm = T)))
} # weighted variance

.covarwt=function(x, y, wt, xcen, ycen){
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  return=(sum(wt * (x - xcen) * (y - ycen), na.rm = T) / (sum(wt, na.rm = T)))
} # weighted covariance

.cov2eigval=function(sx, sy, sxy){
  b=-sx^2-sy^2
  c=sx^2*sy^2-sxy^2
  return=(list(hi=(-b+sqrt(b^2-4*c))/2,lo=(-b-sqrt(b^2-4*c))/2))
} # covariance eigenvalues

.cov2eigvec=function(sx,sy,sxy){
  eigval=.cov2eigval(sx,sy,sxy)$hi
  eigvec=(sx^2-eigval)/sxy
  return=(eigvec)
}

.losvd = function(v_los, m_v, s_v, h3, h4){
  y = (v_los - m_v) / s_v # value used in the LOSVD equation
  L_v =  (((1 / (sqrt(2 * pi))) * exp(-((y^2) / 2))) / s_v) *
    (1 + (h3 * (1 / sqrt(6)) * ((2 * sqrt(2) * y^3) + 3 * sqrt(2) * y)) + (h4 * (1 / sqrt(24)) * ((4 * y^4) - (12 * y^2) + 3)))
  return(L_v)
} # line of sight velocity distribution

.gh_sim=function(par=c(1, 1), v_los, m_v, s_v){
  fit_los=.losvd(v_los, m_v, s_v, h3=par[1], h4=par[2])
  obs_los=approxfun(density(v_los))
  total_los = obs_los(v_los)
  return=sum(((total_los) - (fit_los))^2)
} # Gauss-Hermite coefficients


.circular_ap=function(sbin){
  ap_region = matrix(data = 0, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  x = matrix(data = rep(seq(1,sbin), each=sbin), nrow = sbin, ncol = sbin)
  y = matrix(data = rep(seq(sbin,1), sbin), nrow = sbin, ncol = sbin)
  xx = x - xcentre; yy = y - ycentre
  rr = sqrt(xx^2 + yy^2)
  ap_region[rr<= sbin/2] = 1
  return(ap_region)
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
  return(ap_region)
}

.SFTtoAge = function(a){
  celestial::cosdistTravelTime((1 / a) - 1)
}

.part_spec = function(Metallicity, Age, Mass){
  Z = ProSpect::interp_param(Metallicity, ProSpect::BC03lr$Z, log = TRUE)
  A = ProSpect::interp_param(Age, ProSpect::BC03lr$Age, log = TRUE)

  weights = data.frame("hihi" = Z$wt_hi * A$wt_hi,
                       "hilo" = Z$wt_hi * A$wt_lo,
                       "lohi" = Z$wt_lo * A$wt_hi,
                       "lolo" = Z$wt_lo * A$wt_lo)

  part_spec = array(data = NA, dim = c(1, length(ProSpect::BC03lr$Wave)))

  part_spec = ((ProSpect::BC03lr$Zspec[[Z$ID_hi]][A$ID_hi,] * weights$hihi) +
               (ProSpect::BC03lr$Zspec[[Z$ID_hi]][A$ID_lo,] * weights$hilo) +
               (ProSpect::BC03lr$Zspec[[Z$ID_lo]][A$ID_hi,] * weights$lohi) +
               (ProSpect::BC03lr$Zspec[[Z$ID_lo]][A$ID_lo,] * weights$lolo)) * Mass

  return(part_spec)
}

.reorient_galaxy_OLD = function(galaxy_data){

  # step 1: define the moment of inertia tensor
  inertiaTensor = matrix(nrow=3, ncol=3)
  inertiaTensor[1,1] = sum(galaxy_data$Mass * (galaxy_data$y^2 + galaxy_data$z^2))
  inertiaTensor[2,2] = sum(galaxy_data$Mass * (galaxy_data$x^2 + galaxy_data$z^2))
  inertiaTensor[3,3] = sum(galaxy_data$Mass * (galaxy_data$x^2 + galaxy_data$y^2))
  inertiaTensor[1,2] = -sum(galaxy_data$Mass*galaxy_data$x*galaxy_data$y)
  inertiaTensor[1,3] = -sum(galaxy_data$Mass*galaxy_data$x*galaxy_data$z)
  inertiaTensor[2,3] = -sum(galaxy_data$Mass*galaxy_data$y*galaxy_data$z)
  inertiaTensor[2,1] = inertiaTensor[1,2]
  inertiaTensor[3,1] = inertiaTensor[1,3]
  inertiaTensor[3,2] = inertiaTensor[2,3]


  # step 2: find eigen vectors and reorder such that x is the major axis
  eigen_vec = eigen(inertiaTensor)$vectors
  eigen_vec = cbind(eigen_vec[,3], eigen_vec[,2], eigen_vec[,1])
  rot_mat = t(eigen_vec) # tranpose to get rotation matrix

  # step 5: rotate coordinates using rotation matrix
  new_coor =  rot_mat %*% rbind(galaxy_data$x, galaxy_data$y, galaxy_data$z)
  new_vel =  rot_mat %*% rbind(galaxy_data$vx, galaxy_data$vy, galaxy_data$vz)

  galaxy_data$x = new_coor[1,]; galaxy_data$y = new_coor[2,]; galaxy_data$z = new_coor[3,];
  galaxy_data$vx = new_vel[1,]; galaxy_data$vy = new_vel[2,]; galaxy_data$vz = new_vel[3,];

  return(galaxy_data)
}

.reorient_galaxy = function(galaxy_data){
  r = r_galaxy(galaxy_data)
  J = angmom_galaxy(galaxy_data[r < 0.33*max(r),]) # compute J
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

  new_coor =  rot_mat %*% rbind(galaxy_data$x, galaxy_data$y, galaxy_data$z)
  new_vel =  rot_mat %*% rbind(galaxy_data$vx, galaxy_data$vy, galaxy_data$vz)

  galaxy_data$x = new_coor[1,]; galaxy_data$y = new_coor[2,]; galaxy_data$z = new_coor[3,];
  galaxy_data$vx = new_vel[1,]; galaxy_data$vy = new_vel[2,]; galaxy_data$vz = new_vel[3,];

  return(galaxy_data)
}

.flip = function(galaxy_data){

  rot_mat = matrix(c(1,0,0,0,-1,0,0,0,-1), nrow = 3, ncol=3)

  new_coor =  rot_mat %*% rbind(galaxy_data$x, galaxy_data$y, galaxy_data$z)
  new_vel =  rot_mat %*% rbind(galaxy_data$vx, galaxy_data$vy, galaxy_data$vz)

  galaxy_data$x = new_coor[1,]; galaxy_data$y = new_coor[2,]; galaxy_data$z = new_coor[3,];
  galaxy_data$vx = new_vel[1,]; galaxy_data$vy = new_vel[2,]; galaxy_data$vz = new_vel[3,];

  return(galaxy_data)
}

.masslum2flux = function(X, obs__data, image__grid, lengths__grid, lumdist){
  if (lengths__grid[[X]] > 0){ # if there are particles in that grid cell
    lum = sum(obs__data$galaxy_obs$Lum[image__grid[[X]]]) # sum luminosities in each cell
    flux = lum * lumdist * 1e23 # convert to flux in Janskys
    #flux = lum * ProSpect::Lum2FluxFactor(z = redshift, ref="Planck") # in CGS
  } else {flux=0}
  return(flux)
}

.rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

.image_nan <- function(z, zlim, col, na.color='gray', ...){
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.na <- zlim[2] + zstep # new z for NA

  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na

  zlim[2] <- zlim[2] + zstep # extend top limit to include the new value na

  col <- c(col, na.color) # we construct the new color range by including: na.color

  magicaxis::magimage(z=z, zlim=zlim, col=col, ...) # we finally call image(...)
}

.ellipse = function(a, b, x, y, ang_deg){
  ang_rad = ang_deg * (pi/180)
  r = ((x*cos(ang_rad) + y*sin(ang_rad))^2 /a^2) + ((x*sin(ang_rad) - y*cos(ang_rad))^2 / b^2)
  return(r)
}

.assign_flux = function(X, obs__data, image__grid, lengths__grid, temp_filt, redshift, lumdist){
  if (lengths__grid[[X]] > 0){
    flux_each = numeric(lengths__grid[[X]])
    metals = obs__data$galaxy_obs$Metallicity[image__grid[[X]]]
    ages = obs__data$galaxy_obs$Age[image__grid[[X]]] * 1e9
    masses = obs__data$galaxy_obs$Mass[image__grid[[X]]] * 1e10
    for (Y in 1:lengths__grid[[X]]){
      spectra = .part_spec(Metallicity = metals[Y],
                           Age = ages[Y],
                           Mass = masses[Y])
      flux_each[Y] = ProSpect::photom_lum(wave = ProSpect::BC03lr$Wave, lum = spectra, outtype = "Jansky",
                                          filters = temp_filt, z = redshift, ref = "Planck", LumDist_Mpc = lumdist)
    }
    flux = sum(flux_each)
  } else {flux=0}
  return(flux)
}
