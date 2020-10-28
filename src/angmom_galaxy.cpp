#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Computing the angular momentum of the particles.
//'
//' The purpose of this function is to compute the angular momentum vector.
//'
//' @param part_data The particles from which you would like to calculate the
//' angular momentum vector, J.
//' @return The numeric vector describing the magnitudes of each component of
//' the angular momentum.
//' @examples
//'   galaxy_data = data.frame("ID"=1:100, "x"=stats::rnorm(100),
//'                            "y"=stats::rnorm(100), "z"=stats::rnorm(100),
//'                            "vx"=stats::rnorm(100), "vy"=stats::rnorm(100),
//'                            "vz"=stats::rnorm(100), "Mass"=rep(1,100))
//'   J = angmom_galaxy(galaxy_data)
//' @export
// [[Rcpp::export]]

Rcpp::NumericVector angmom_galaxy(Rcpp::DataFrame part_data) {

  Rcpp::NumericVector ID        = part_data["ID"];
  Rcpp::NumericVector x         = part_data["x"];
  Rcpp::NumericVector y         = part_data["y"];
  Rcpp::NumericVector z         = part_data["z"];
  Rcpp::NumericVector vx        = part_data["vx"];
  Rcpp::NumericVector vy        = part_data["vy"];
  Rcpp::NumericVector vz        = part_data["vz"];
  Rcpp::NumericVector Mass      = part_data["Mass"];

  double Jx = sum(Mass * ((y * vz) - (z * vy)));
  double Jy = sum(Mass * ((x * vz) - (z * vx)));
  double Jz = sum(Mass * ((x * vy) - (y * vx)));

  Rcpp::NumericVector J = Rcpp::NumericVector::create(Jx, Jy, Jz);

  return(J);
}
