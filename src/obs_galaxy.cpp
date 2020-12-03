// Kate Harborne (last edit - 26/10/2020)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Constructing galaxy observation data from the Gadget output file data.
//'
//' The purpose of this function is to produce the observable features of
//' simulation data when taking mock IFU observations.
//'
//' @param part_data A data.frame describing the particles ID, positions, and
//'  velocities.
//' @param inc_rad The observed inclination angle in radians.
//' @return Returns a data frame containing the original particle information plus the observed
//'  z-position (\code{$z_obs}), observed radial position (\code{$r_obs}) and the observed line of
//'  sight velocity (\code{$vy_obs}) at the given inclination.
//' @examples
//'   galaxy_data = data.frame("ID"=1:100, "x"=stats::rnorm(100),
//'                            "y"=stats::rnorm(100), "z"=stats::rnorm(100),
//'                            "vx"=stats::rnorm(100), "vy"=stats::rnorm(100),
//'                            "vz"=stats::rnorm(100))
//'   observed_data = obs_galaxy(galaxy_data, inc_rad = 1.047)
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
  Rcpp::NumericVector Mass         = part_data["Mass"];
  Rcpp::NumericVector sed_id       = part_data["sed_id"];
  Rcpp::NumericVector Initial_Mass = part_data["Initial_Mass"];

  int n = x.size();

  Rcpp::NumericVector r(n), z_obs(n), y_obs(n), vy_obs(n), r_obs(n);
  for(int i=0; i<n; i++){
    r[i]      = sqrt((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]));
    z_obs[i]  = sin(inc_rad) * z[i] + cos(inc_rad) * y[i];
    y_obs[i]  = cos(inc_rad) * z[i] - sin(inc_rad) * y[i];
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
                            Rcpp::Named("Mass")         = Mass,
                            Rcpp::Named("sed_id")       = sed_id,
                            Rcpp::Named("Initial_Mass") = Initial_Mass,
                            Rcpp::Named("r")         = r,
                            Rcpp::Named("z_obs")     = z_obs,
                            Rcpp::Named("y_obs")     = y_obs,
                            Rcpp::Named("r_obs")     = r_obs,
                            Rcpp::Named("vy_obs")    = vy_obs);

  return(df);

}
