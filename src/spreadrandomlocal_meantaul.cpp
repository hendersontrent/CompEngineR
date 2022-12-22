#include <Rcpp.h>
using namespace Rcpp;

//' Calculate bootstrap-based stationarity measure from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @param l \code{integer} denoting the length of local time-series segments to analyse
//' @return \code{numeric} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' spreadrandomlocal_meantaul(x)
//'
// [[Rcpp::export]]
double spreadrandomlocal_meantaul(NumericVector x, int l = 50) {
  if (l > 0.9 * x.size()) {
    warning("The time series is too short. Specify proper segment length in `l`");
    return NA_REAL;
  }

  int numSegs = 100;
  NumericVector qs(numSegs);

  for (int j = 0; j < numSegs; j++) {
    int ist = sample(x.size() - 1 - l, 1)[0];
    int ifh = ist + l - 1;
    IntegerVector rs = seq(ist, ifh);
    NumericVector xsub = x[rs];
    int taul = firstzero_ac(xsub);
    qs[j] = taul;
  }
  return mean(qs, na_rm = true);
}
