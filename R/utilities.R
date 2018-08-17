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
