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

