#include <Rcpp.h>
using namespace Rcpp;

//' Calculate a simulation of a hypothetical walker moving through the time domain from software package \code{hctsa}
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
//' walker_propcross(x)
//'
// [[Rcpp::export]]
double walker_propcross(NumericVector y) {
  int N = y.size();
  double p = 0.1;
  NumericVector w(N);
  w[0] = 0; // Start at zero

  for (int i = 1; i < N; i++) {
    w[i] = w[i - 1] + p * (y[i - 1] - w[i - 1]);
  }

  double out_sw_propcross = sum((w[Range(0, N - 2)] - y[Range(0, N - 2)]) * (w[Range(1, N - 1)] - y[Range(1, N - 1)]) < 0) / (N - 1);
  return out_sw_propcross;
}
