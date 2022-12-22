#include <Rcpp.h>
using namespace Rcpp;

//' Calculate Second Sample Entropy from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @param M \code{integer} embedding dimension
//' @param r \code{numeric} threshold
//' @return \code{numeric} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' sampenc(x, 6, 0.3)
//'
// [[Rcpp::export]]
double sampenc(NumericVector x, int M = 6, double r = 0.3) {
  int N = x.size();
  NumericVector lastrun(N);
  NumericVector run(N);
  NumericVector A(M);
  NumericVector B(M);

  // Go through each point in the time series, counting matches

  for (int i = 0; i < N - 1; i++) {
    double x1 = x[i];

    // Compare to points through the rest of the time series
    for (int jj = 0; jj < N - i; jj++) {
      // Compare to future index, j
      int j = i + jj;
      // This future point, j, matches the time-series value at i
      if (std::abs(x[j] - x1) < r) {
        run[jj] = lastrun[jj] + 1;
        int M1 = std::min(M, run[jj]);
        A[Range(0, M1 - 1)] += 1;
        if (j < N) B[Range(0, M1 - 1)] += 1;
      } else {
        run[jj] = 0;
      }
    }
    for (int j = 0; j < N - i; j++) {
      lastrun[j] = run[j];
    }
  }

  double p = A[1] / B[0];
  double e = -std::log(p);
  return e;
}
