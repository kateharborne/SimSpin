// Kate Harborne (last edit - 03/09/19)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Centering the galaxy.
//'
//' The purpose of this function is to centre the galaxy such that the origin of the
//' system lies at (0,0,0).
//'
//' @param part_data The concatenated data frames output by \code{\link{sim_data}}.
//' @examples
//'   data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
//'   galaxy_data = rbind(data$PartType2$Part, data$PartType3$Part)
//'
//'   output = cen_galaxy(part_data = galaxy_data)
//' @export
// [[Rcpp::export]]

Rcpp::List cen_galaxy(Rcpp::DataFrame part_data) {

  Rcpp::NumericVector ID        = part_data["ID"];
  Rcpp::NumericVector x         = part_data["x"];
  Rcpp::NumericVector y         = part_data["y"];
  Rcpp::NumericVector z         = part_data["z"];
  Rcpp::NumericVector vx        = part_data["vx"];
  Rcpp::NumericVector vy        = part_data["vy"];
  Rcpp::NumericVector vz        = part_data["vz"];
  Rcpp::NumericVector Mass      = part_data["Mass"];

  double xcen = median(x);
  double ycen = median(y);
  double zcen = median(z);
  double vxcen = median(vx);
  double vycen = median(vy);
  double vzcen = median(vz);
  x = x-xcen;
  y = y-ycen;
  z = z-zcen;
  vx = vx - vxcen;
  vy = vy - vycen;
  vz = vz - vzcen;

  Rcpp::DataFrame df =
    Rcpp::DataFrame::create(Rcpp::Named("ID")        = ID,
                            Rcpp::Named("x")         = x,
                            Rcpp::Named("y")         = y,
                            Rcpp::Named("z")         = z,
                            Rcpp::Named("vx")        = vx,
                            Rcpp::Named("vy")        = vy,
                            Rcpp::Named("vz")        = vz,
                            Rcpp::Named("Mass")      = Mass);

  return(df);

}
