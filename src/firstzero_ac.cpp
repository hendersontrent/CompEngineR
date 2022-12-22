#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the first zero crossing of the autocorrelation function from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @param acfv \code{numeric} vector of autocorrelation values
//' @return \code{integer}
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
  IntegerVector tau = which(acfv[Range(1, acfv.size() - 1)] < 0);

  if (tau.size() == 0) // Nothing to see here
    return 0;
  else if (all(is_na(tau))) // All missing
    return 0;
  else if (!any(tau))  // No negatives, so set output to sample size
    return N;
  else // Return lag of first negative
    return tau[0];
}
