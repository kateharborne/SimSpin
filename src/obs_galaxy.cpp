// Kate Harborne (last edit - 13/09/2017)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Constructing galaxy observation data from the Gadget output file data.
//'
//' The purpose of this function is to produce the observable features of simulation data when
//' taking mock IFU observations. It accepts the Gadget particle information output from the
//' \code{snapshot::snapread} function and returns the observable galaxy properties projected at a
//' user supplied inclination.
//'
//' @param part_data The data frame output by \code{snapshot::snapread} for galaxy simulation data
//'  from Gadget.
//' @param centre A logical that tells the function to centre the galaxy about its centre of mass or
//'  not (i.e. TRUE or FALSE).
//' @param inc_rad The observed inclination angle in radians.
//' @return Returns a data frame containing the original particle information plus the observed
//'  z-position (\code{$z_obs}), observed radial position (\code{$r_obs}) and the observed line of
//'  sight velocity (\code{$vy_obs}) at the given inclination.
//' @examples
//'   galaxy_file = h5::h5file(system.file("extdata", 'S0_vignette', package="SimSpin"), mode = "r")
//'   galaxy_data = data.frame("x"         = h5::readDataSet(galaxy_file["x"]),
//'                            "y"         = h5::readDataSet(galaxy_file["y"]),
//'                            "z"         = h5::readDataSet(galaxy_file["z"]),
//'                            "vx"        = h5::readDataSet(galaxy_file["vx"]),
//'                            "vy"        = h5::readDataSet(galaxy_file["vy"]),
//'                            "vz"        = h5::readDataSet(galaxy_file["vz"]),
//'                            "Mass"      = h5::readDataSet(galaxy_file["Mass"]),
//'                            "part_type" = h5::readDataSet(galaxy_file["Part_Type"]))
//'   h5::h5close(galaxy_file)
//'
//'   output = obs_galaxy(part_data = galaxy_data,
//'                       centre    = TRUE,
//'                       inc_rad   = 0)
//' @export
// [[Rcpp::export]]
Rcpp::List obs_galaxy(Rcpp::DataFrame part_data, bool centre, double inc_rad) {

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
  Rcpp::NumericVector r(n), z_obs(n), vy_obs(n), r_obs(n);
  for(int i=0; i<n; i++){
    r[i]      = ::sqrt((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]));
    z_obs[i]  = ::sin(inc_rad) * z[i] + ::cos(inc_rad) * y[i];
    vy_obs[i] = ::cos(inc_rad) * vz[i] - ::sin(inc_rad) * vy[i];
    r_obs[i]  = ::sqrt((x[i] * x[i]) + (z_obs[i] * z_obs[i]));
  }

  Rcpp::DataFrame df =
    Rcpp::DataFrame::create(Rcpp::Named("x")         = x,
                            Rcpp::Named("y")         = y,
                            Rcpp::Named("z")         = z,
                            Rcpp::Named("vx")        = vx,
                            Rcpp::Named("vy")        = vy,
                            Rcpp::Named("vz")        = vz,
                            Rcpp::Named("Mass")      = Mass,
                            Rcpp::Named("r")         = r,
                            Rcpp::Named("z_obs")     = z_obs,
                            Rcpp::Named("r_obs")     = r_obs,
                            Rcpp::Named("vy_obs")    = vy_obs,
                            Rcpp::Named("part_type") = part_type);

  return(df);

}
