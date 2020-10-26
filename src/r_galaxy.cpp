// Kate Harborne (last edit - 03/09/19)
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Computing radius.
//'
//' The purpose of this function is to compute the radial coordinates.
//'
//' @param part_data The concatenated data frames output by \code{\link{sim_data}}.
//' @export
// [[Rcpp::export]]

Rcpp::NumericVector r_galaxy(Rcpp::DataFrame part_data) {

  Rcpp::NumericVector ID        = part_data["ID"];
  Rcpp::NumericVector x         = part_data["x"];
  Rcpp::NumericVector y         = part_data["y"];
  Rcpp::NumericVector z         = part_data["z"];
  Rcpp::NumericVector vx        = part_data["vx"];
  Rcpp::NumericVector vy        = part_data["vy"];
  Rcpp::NumericVector vz        = part_data["vz"];
  Rcpp::NumericVector Mass      = part_data["Mass"];

  Rcpp::NumericVector r = sqrt((x * x) + (y * y) + (z * z));

  return(r);

}
