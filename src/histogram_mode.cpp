#include <Rcpp.h>
using namespace Rcpp;

//' Calculate mode of a data vector from software package \code{hctsa}
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
//' histogram_mode(x)
//'
// [[Rcpp::export]]
double histogram_mode(NumericVector x, int numBins = 10) {
  if (numBins <= 0) stop("Invalid number of bins");

  List histdata = hist(x, plot = false, breaks = numBins);
  NumericVector binCenters = histdata["mids"];
  NumericVector counts = histdata["counts"];
  IntegerVector maxIdx = which_max(counts);
  return mean(binCenters[maxIdx]);
}
