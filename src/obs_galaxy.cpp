// Kate Harborne (last edit - 26/10/2020)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Constructing galaxy observation data from the Gadget output file data.
//'
//' The purpose of this function is to produce the observable features of
//' simulation data when taking mock IFU observations.
//'
//' @param part_data A data.frame describing the particle positions and
//'  velocities.
//' @param inc_rad The observed inclination angle in radians.
//' @return Returns a data frame containing the original particle information plus the observed
//'  z-position (\code{$z_obs}), observed radial position (\code{$r_obs}) and the observed line of
//'  sight velocity (\code{$vy_obs}) at the given inclination.
//' @examples
//'   galaxy_data = data.frame("x"=stats::rnorm(100), "y"=stats::rnorm(100),
//'                            "z"=stats::rnorm(100), "vx"=stats::rnorm(100),
//'                            "vy"=stats::rnorm(100), "vz"=stats::rnorm(100))
//'   observed_data = obs_galaxy(galaxy_data, inc_rad = 1.047)
//' @export
// [[Rcpp::export]]
Rcpp::List obs_galaxy(Rcpp::DataFrame part_data, double inc_rad) {

  Rcpp::NumericVector x         = part_data["x"];
  Rcpp::NumericVector y         = part_data["y"];
  Rcpp::NumericVector z         = part_data["z"];
  Rcpp::NumericVector vx        = part_data["vx"];
  Rcpp::NumericVector vy        = part_data["vy"];
  Rcpp::NumericVector vz        = part_data["vz"];

  int n = x.size();

  Rcpp::NumericVector y_obs(n), z_obs(n), vy_obs(n), vz_obs(n);
  for(int i=0; i<n; i++){
    y_obs[i]  = cos(inc_rad) * z[i] - sin(inc_rad) * y[i];
    z_obs[i]  = sin(inc_rad) * z[i] + cos(inc_rad) * y[i];
    vy_obs[i] = cos(inc_rad) * vz[i] - sin(inc_rad) * vy[i];
    vz_obs[i] = sin(inc_rad) * vz[i] + cos(inc_rad) * vy[i];
  }

  Rcpp::DataFrame df =
    Rcpp::DataFrame::create(Rcpp::Named("x")         = x,
                            Rcpp::Named("y")         = y_obs,
                            Rcpp::Named("z")         = z_obs,
                            Rcpp::Named("vx")        = vx,
                            Rcpp::Named("vy")        = vy_obs,
                            Rcpp::Named("vz")        = vz_obs);

  return(df);

}
