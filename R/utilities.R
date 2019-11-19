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

.circular_ap_cut=function(galaxy_df, ap_size){
  trimdata = galaxy_df[galaxy_df$r_obs < (ap_size / 2),]
  return(trimdata)
}

.square_ap_cut=function(galaxy_df, sbin, sbinsize){
  galaxy_cdf  = galaxy_df[abs(galaxy_df$x) < (sbin / 2) * sbinsize,]
  trimdata    = galaxy_cdf[abs(galaxy_cdf$z_obs) < (sbin / 2) * sbinsize,]
  return(trimdata)
}

.hexagonal_ap_cut=function(galaxy_df, sbin, sbinsize){
  galaxy_cdf  = galaxy_df[abs(galaxy_df$x) < (sbin / 2) * sbinsize,]
  galaxy_cdf  = galaxy_cdf[abs(galaxy_cdf$z_obs) < (sbin  * sqrt(3) / 4) * sbinsize,]
  dotprod     = (2 * (sbin / 4) * sbinsize * (sbin * sqrt(3) / 4) * sbinsize) - ((sbin / 4) * sbinsize) * abs(galaxy_cdf$z_obs) - ((sbin * sqrt(3) / 4) * sbinsize) * abs(galaxy_cdf$x)
  trimdata  = galaxy_cdf[dotprod >= 0,]
  return(trimdata)
}

.SFTtoAge = function(x){
  celestial::cosdistTravelTime((1 / x) - 1)
}

.part_spec = function(gal_0){
  Z = ProSpect::interp_param(gal_0$Metallicity, ProSpect::BC03lr$Z, log = FALSE)
  A = ProSpect::interp_param(gal_0$Age, ProSpect::BC03lr$Age, log = FALSE)

  weights = data.frame("hihi" = Z$weight_hi * A$weight_hi,
                       "hilo" = Z$weight_hi * A$weight_lo,
                       "lohi" = Z$weight_lo * A$weight_hi,
                       "lolo" = Z$weight_lo * A$weight_lo)

  part_spec = array(data = NA, dim = c(nrow(weights), length(ProSpect::BC03lr$Wave)))

  for (i in 1:nrow(weights)){
    part_spec[i,] = ((ProSpect::BC03lr$Zspec[[Z$ID_hi[i]]][A$ID_hi[i],] * weights$hihi[i]) + (ProSpect::BC03lr$Zspec[[Z$ID_hi[i]]][A$ID_lo[i],] * weights$hilo[i]) +
                       (ProSpect::BC03lr$Zspec[[Z$ID_lo[i]]][A$ID_hi[i],] * weights$lohi[i]) + (ProSpect::BC03lr$Zspec[[Z$ID_lo[i]]][A$ID_lo[i],] * weights$lolo[i])) * gal_0$Mass[i] * 1e10

  }

  return(part_spec)

}

.reorient_galaxy = function(galaxy_data){

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

.assign_flux = function(X, obs__data, image__grid, lengths__grid, temp_filt, redshift){
      if (lengths__grid[[X]] > 1){ # if there are particles in the X
        temp = data.frame("Age" = obs__data$galaxy_obs$Age[image__grid[[X]]] * 1e9,
                          "Mass" = obs__data$galaxy_obs$Mass[image__grid[[X]]] * 1e10,
                          "cumMass" = rep(NA, length(obs__data$galaxy_obs$Mass[image__grid[[X]]])),
                          "Metallicity" = obs__data$galaxy_obs$Metallicity[image__grid[[X]]])
        temp = temp[order(temp$Age),] # fill a data frame and order by age

      if (length(unique(temp$Age)) != length(temp$Age)){
        # if there are any particles with the same age
        nu = as.integer(which(table(temp$Age) > 1))
        for (n in 1:length(nu)){
          vals = which(temp$Age == temp$Age[nu[n]])
          temp$Mass[vals[1]] = sum(temp$Mass[vals])
          temp$Metallicity[vals[1]] = sum(temp$Metallicity[vals])
          vals = vals[-1]
          temp = temp[-vals,]
          } # sum together all masses and metallicities such that we have a single particle at that age
      }

      if (length(temp$Age) > 1){ # if two particles are summed above to give more than a single particle
          temp$cumMass = cumsum(temp$Mass) # cumulative sum of ages
          tmp_SFHfunc = approxfun(x = c(temp$Age), y = c(0,diff(temp$cumMass)/diff(temp$Age)), yleft=0, yright=0)
          # describe star formation rate in M_sol/yr as a function of age
          tmp_Zfunc = approxfun(x = c(temp$Age), y = c(temp$Metallicity), yleft=0, yright=0)
          # describe metallicity as a function of age
          tmp_SFH = ProSpect::SFHfunc(massfunc = tmp_SFHfunc, Z = tmp_Zfunc, z=redshift, speclib = ProSpect::BC03lr,
                                      outtype="Jansky",
                                      filters = temp_filt, ref="Planck")
          # generate a spectrum
          flux = tmp_SFH$out$out
          # convert to flux
        }
      } else {flux=0}

  return(flux)
}

.masslum2flux = function(X, obs__data, image__grid, lengths__grid, redshift){
  if (lengths__grid[[X]] > 1){ # if there are particles in that grid cell
    lum = sum(obs__data$galaxy_obs$Lum[image__grid[[X]]]) # sum luminosities in each cell
    flux = ProSpect::CGS2Jansky(lum * ProSpect::Lum2FluxFactor(z = redshift, ref="Planck")) # convert to flux in Janskys
  } else {flux=0}
  return(flux)
}

.rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

.image_nan <- function(z, zlim, col, na.color='gray', ...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.na <- zlim[2] + zstep # new z for NA

  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na

  zlim[2] <- zlim[2] + zstep # extend top limit to include the new value  na

  col <- c(col, na.color) # we construct the new color range by including: na.color

  magicaxis::magimage(z=z, zlim=zlim, col=col, ...) # we finally call image(...)
}
