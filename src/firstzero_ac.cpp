#include <Rcpp.h>
#include <algorithm>
#include "firstzero_ac.hpp"
using namespace Rcpp;

// Define a helper function for making a logical vector to see if AC value is < 0
 Rcpp::LogicalVector get_negative(Rcpp::NumericVector x) {
   int n = x.size();
   Rcpp::LogicalVector out(n);

   for (int i = 0; i < n; i++) {
     out[i] = x[i] < 0;
   }
   return out;
 }

// Define a helper function for finding index of first `true` value
int get_first_true_index(Rcpp::LogicalVector x) {
  auto it = std::find(x.begin(), x.end(), true);
  if (it != x.end()) {
    return std::distance(x.begin(), it);  // Rcpp uses 0-based indexing
  } else {
    return -1;  // no TRUE values found
  }
}

//' Calculate the first zero crossing of the autocorrelation function from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @param acfv \code{numeric} vector of autocorrelation values
//' @return \code{integer} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' firstzero_ac(x)
//'
// [[Rcpp::export]]
int firstzero_ac(NumericVector x, NumericVector acfv) {
  int N = x.size();
  LogicalVector tau = get_negative(acfv[Range(1, N - 1)]);

  if (tau.size() == 0) {
    return 0;
  } else if (all(is_na(tau))) {
    return 0;
  } else if (!std::any_of(x.begin(), x.end(), [](bool b) { return b; })) {
    return N;
  } else {
    return get_first_true_index(tau) + 1;
  }
}
