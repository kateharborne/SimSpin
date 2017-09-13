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
