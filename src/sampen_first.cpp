#include <Rcpp.h>
#include "sampenc.hpp"
using namespace Rcpp;

//' Calculate Second Sample Entropy of a time series from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @return \code{numeric} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' sampen_first(x)
//'
// [[Rcpp::export]]
double sampen_first(NumericVector x) {
  int M = 5;
  double r = 0.3;
  double sampEn = sampenc(x, M + 1, r);
  return sampEn;
}
