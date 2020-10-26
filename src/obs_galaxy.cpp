// Kate Harborne (last edit - 26/10/2020)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Constructing galaxy observation data from the Gadget output file data.
//'
//' The purpose of this function is to produce the observable features of simulation data when
//' taking mock IFU observations. It accepts particle information output from the
//' \code{\link{sim_data}} function and returns the observable galaxy properties projected at a
//' user supplied inclination.
//'
//' @param part_data The concatenated data frames output by \code{\link{sim_data}}.
//' @param inc_rad The observed inclination angle in radians.
//' @return Returns a data frame containing the original particle information plus the observed
//'  z-position (\code{$z_obs}), observed radial position (\code{$r_obs}) and the observed line of
//'  sight velocity (\code{$vy_obs}) at the given inclination.
//' @export
// [[Rcpp::export]]
Rcpp::List obs_galaxy(Rcpp::DataFrame part_data, double inc_rad) {

  Rcpp::NumericVector ID        = part_data["ID"];
  Rcpp::NumericVector x         = part_data["x"];
  Rcpp::NumericVector y         = part_data["y"];
  Rcpp::NumericVector z         = part_data["z"];
  Rcpp::NumericVector vx        = part_data["vx"];
  Rcpp::NumericVector vy        = part_data["vy"];
  Rcpp::NumericVector vz        = part_data["vz"];

  int n = x.size();

  Rcpp::NumericVector r(n), z_obs(n), vy_obs(n), r_obs(n);
  for(int i=0; i<n; i++){
    r[i]      = sqrt((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]));
    z_obs[i]  = sin(inc_rad) * z[i] + cos(inc_rad) * y[i];
    vy_obs[i] = cos(inc_rad) * vz[i] - sin(inc_rad) * vy[i];
    r_obs[i]  = sqrt((x[i] * x[i]) + (z_obs[i] * z_obs[i]));
  }

  Rcpp::DataFrame df =
    Rcpp::DataFrame::create(Rcpp::Named("ID")        = ID,
                            Rcpp::Named("x")         = x,
                            Rcpp::Named("y")         = y,
                            Rcpp::Named("z")         = z,
                            Rcpp::Named("vx")        = vx,
                            Rcpp::Named("vy")        = vy,
                            Rcpp::Named("vz")        = vz,
                            Rcpp::Named("r")         = r,
                            Rcpp::Named("z_obs")     = z_obs,
                            Rcpp::Named("r_obs")     = r_obs,
                            Rcpp::Named("vy_obs")    = vy_obs);

  return(df);

}
