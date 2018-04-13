// Kate Harborne (last edit - 13/09/2017)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Constructing galaxy simulation data from the Gadget output file data
//'
//' The purpose of this function is to produce the extra kinematic features for simulation data in spherical polar coordinates. It accepts
//' the Gadget particle information output from the \code{snapshot::snapread} function and returns several additional galaxy
//' properties that are useful for deriving the galaxy kinematics.
//'
//' @param part_data The data frame output by \code{snapshot::snapread} for galaxy simulation data in Gadget format.
//' @param centre A logical that tells the function to centre the galaxy about its centre of mass or not (i.e. TRUE or FALSE).
//' @return Returns a data frame containing the particle \code{$ID}, \code{$x-}, \code{$y-} and \code{$z-}positions and corresponding
//'  velocities (\code{$vx, $vy } and \code{$vz}), along with the spherical polar coordinates (\code{$r}, \code{$theta} and \code{$phi}) and
//'  associated velocities (\code{$vr}, \code{$vtheta} and \code{$vphi}), the cylindrical radial coordinate and its associated velocity (\code{$cr}
//'  and \code{$vcr})  and the mass of each particle (\code{$Mass}) and their angular momentum magnitude (\code{$J}).
//' @examples
//'  \dontrun{
//'   data = snapshot::snapread()
//'
//'   sim_galaxy(part_data = data$part,
//'              centre    = TRUE)
//'  }
//' @export
// [[Rcpp::export]]
Rcpp::List sim_galaxy(Rcpp::DataFrame part_data, bool centre) {

  Rcpp::NumericVector ID        = part_data["ID"];
  Rcpp::NumericVector x         = part_data["x"];
  Rcpp::NumericVector y         = part_data["y"];
  Rcpp::NumericVector z         = part_data["z"];
  Rcpp::NumericVector vx        = part_data["vx"];
  Rcpp::NumericVector vy        = part_data["vy"];
  Rcpp::NumericVector vz        = part_data["vz"];
  Rcpp::NumericVector Mass      = part_data["Mass"];
  Rcpp::NumericVector part_type = part_data["part_type"];

  int n = x.size();
  if(centre == TRUE){
    double xcen = median(x);
    double ycen = median(y);
    double zcen = median(z);
    x = x-xcen;
    y = y-ycen;
    z = z-zcen;
  }
  Rcpp::NumericVector r(n), cr(n), theta(n), phi(n), vr(n), vt(n), vcr(n), vtheta(n), vphi(n), Jx(n), Jy(n), Jz(n);
  for(int i=0; i<n; i++){
    r[i]      = ::sqrt((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]));
    cr[i]     = ::sqrt((x[i] * x[i]) + (y[i] * y[i]));
    theta[i]  = ::acos(z[i] / r[i]);
    phi[i]    = ::atan(y[i] / x[i]);
    vr[i]     = (x[i] * vx[i] + y[i] * vy[i] + z[i] * vz[i]) /r[i];
    vt[i]     = ::sqrt(((y[i] * vz[i] - z[i] * vy[i]) * (y[i] * vz[i] - z[i] * vy[i])) + ((z[i] * vx[i] - x[i] * vz[i]) * (z[i] * vx[i] - x[i] * vz[i])) + ((x[i] * vy[i] - y[i] * vx[i]) * (x[i] * vy[i] - y[i] * vx[i]))) / r[i];
    vcr[i]    = (x[i] * vx[i] + y[i] * vy[i]) /cr[i];
    vtheta[i] = -(r[i] / ::sqrt(1 - ((z[i] * z[i]) / (r[i] * r[i])))) * ((vz[i] / r[i]) - (vr[i] * z[i] / (r[i] * r[i])));
    vphi[i]   = (r[i] / (((y[i] * y[i]) / (x[i] * x[i])) + 1)) * ::sin(theta[i]) * ((vy[i] / x[i]) - (vx[i] * y[i] / (x[i] * x[i])));
    Jx[i]     = Mass[i] * ((y[i] * vz[i]) - (z[i] * vy[i])) * 3.086e16;
    Jy[i]     = Mass[i] * ((z[i] * vx[i]) - (x[i] * vz[i])) * 3.086e16;
    Jz[i]     = Mass[i] * ((x[i] * vy[i]) - (y[i] * vx[i])) * 3.086e16;
  }

  Rcpp::DataFrame df =
    Rcpp::DataFrame::create(Rcpp::Named("x")         = x,
                            Rcpp::Named("y")         = y,
                            Rcpp::Named("z")         = z,
                            Rcpp::Named("vx")        = vx,
                            Rcpp::Named("vy")        = vy,
                            Rcpp::Named("vz")        = vz,
                            Rcpp::Named("r")         = r,
                            Rcpp::Named("cr")        = cr,
                            Rcpp::Named("theta")     = theta,
                            Rcpp::Named("phi")       = phi,
                            Rcpp::Named("vr")        = vr,
                            Rcpp::Named("vt")        = vt,
                            Rcpp::Named("vcr")       = vcr,
                            Rcpp::Named("vtheta")    = vtheta,
                            Rcpp::Named("vphi")      = vphi,
                            Rcpp::Named("Mass")      = Mass,
                            Rcpp::Named("Jx")        = Jx,
                            Rcpp::Named("Jy")        = Jy,
                            Rcpp::Named("Jz")        = Jz,
                            Rcpp::Named("part_type") = part_type);

  return(df);

}
