#include <Rcpp.h>
using namespace Rcpp;

//' Calculate time of first minimum in the autocorrelation function from software package \code{hctsa}
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
//' firstmin_ac(x)
//'
// [[Rcpp::export]]
int firstmin_ac(NumericVector x, NumericVector acfv) {

  int N = x.size();

  // getting acf for all lags
  // possible delay when sample size is too big

  NumericVector autoCorr(N - 1);
  autoCorr[Range(0, N - 2)] = acfv[1:acfv.size() - 1];

  for (int i = 0; i < autoCorr.size(); ++i) {
    if (NumericVector::is_na(autoCorr[i])) {
      warning("No minimum was found.");
      return NA_INTEGER;
    }
    if (i == 1 && autoCorr[1] > autoCorr[0]) {
      return 1;
    } else if (i > 1 && autoCorr[i - 1] > autoCorr[i] && autoCorr[i] < autoCorr[i + 1]) {
      return i;
    }
  }
  return N - 1;
}
