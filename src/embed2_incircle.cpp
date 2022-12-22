#include <Rcpp.h>
#include "firstzero_ac.hpp"
using namespace Rcpp;

//' Calculate points inside a given circular boundary in a 2-D embedding space from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @param acfv \code{numeric} vector of autocorrelation values
//' @param boundary \code{numeric} boundary for the circle calculation
//' @return \code{numeric} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' embed2_incircle(x)
//'
// [[Rcpp::export]]
double embed2_incircle(NumericVector x, NumericVector acfv, double boundary = NA_REAL) {

  if (NumericVector::is_na(boundary)) {
    warning("`embed2_incircle()` using `boundary = 1`. Set value with `boundary`.");
    boundary = 1;
  }

  double tau = firstzero_ac(x, acfv);
  NumericVector xt = x[Range(0, x.size() - tau - 1)]; // Part of the time series
  NumericVector xtp = x[Range(tau, x.size() - 1)]; // Time-lagged time series
  int N = x.size() - tau; // Length of each time series sub-segment

  // Circles (points inside a given circular boundary)

  return sum(pow(xtp, 2) + pow(xt, 2) < boundary) / N;
}
