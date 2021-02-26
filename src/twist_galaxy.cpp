// Kate Harborne (last edit - 26/10/2020)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Twisting the galaxy relative to the observer.
//'
//' The purpose of this function is to re-orient the galaxy and rotate the
//' observer around the z-axis.
//'
//' @param part_data A data.frame describing the particles positions and
//'  velocities.
//' @param twist_rad The observed inclination angle in radians.
//' @return Returns the original data frame containing the original particle
//'  information but in the new orientation.
//' @examples
//'   galaxy_data = data.frame("x"=stats::rnorm(100), "y"=stats::rnorm(100),
//'                            "z"=stats::rnorm(100), "vx"=stats::rnorm(100),
//'                            "vy"=stats::rnorm(100), "vz"=stats::rnorm(100))
//'   observed_data = twist_galaxy(galaxy_data, twist_rad = 1.047)
//' @export
// [[Rcpp::export]]
Rcpp::List twist_galaxy(Rcpp::DataFrame part_data, double twist_rad) {

  Rcpp::NumericVector x         = part_data["x"];
  Rcpp::NumericVector y         = part_data["y"];
  Rcpp::NumericVector z         = part_data["z"];
  Rcpp::NumericVector vx        = part_data["vx"];
  Rcpp::NumericVector vy        = part_data["vy"];
  Rcpp::NumericVector vz        = part_data["vz"];

  int n = x.size();

  Rcpp::NumericVector x_obs(n), y_obs(n), vx_obs(n), vy_obs(n);
  for(int i=0; i<n; i++){
    x_obs[i]  = cos(twist_rad) * x[i] - sin(twist_rad) * y[i];
    y_obs[i]  = sin(twist_rad) * x[i] + cos(twist_rad) * y[i];
    vx_obs[i] = cos(twist_rad) * vx[i] - sin(twist_rad) * vy[i];
    vy_obs[i] = sin(twist_rad) * vx[i] + cos(twist_rad) * vy[i];
  }

  Rcpp::DataFrame df =
    Rcpp::DataFrame::create(Rcpp::Named("x")         = x_obs,
                            Rcpp::Named("y")         = y_obs,
                            Rcpp::Named("z")         = z,
                            Rcpp::Named("vx")        = vx_obs,
                            Rcpp::Named("vy")        = vy_obs,
                            Rcpp::Named("vz")        = vz);

  return(df);

}
