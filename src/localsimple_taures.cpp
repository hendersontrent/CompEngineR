#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the first zero crossing of the autocorrelation function of the residuals from Simple local time-series forecasting from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @param forecastMeth \code{character}
//' @return \code{integer} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' localsimple_taures(x)
//'
// [[Rcpp::export]]
int localsimple_taures(NumericVector x, String forecastMeth) {

  if (forecastMeth == "mean") {
    lp = 1;
  } else{
    lp = 3;
  }

  int N = x.size();
  IntegerVector evalr = seq(lp + 1, N);

  if (lp >= N) stop("Time series too short for forecasting in `localsimple_taures`");

  NumericVector res(evalr.size());
  if (forecastMeth == "mean") {
    for (int i = 0; i < evalr.size(); i++)
      res[i] = mean(y[Range(evalr[i] - lp, evalr[i] - 1)]) - x[evalr[i]];
  }

  if (forecastMeth == "lfit") {
    for (int i = 0; i < evalr.size(); i++) {
      // Fit linear
      NumericVector a = seq(1, lp);
      NumericVector b = y[Range(evalr[i] - lp, evalr[i] - 1)];
      DataFrame df = DataFrame::create(Named("a")=a, Named("b")=b);
      lm lm_ab = lm(b ~ a, data = df);
      res[i] = predict(lm_ab, newdata = DataFrame::create(Named("a") = lp + 1)) - y[evalr[i]];
    }
  }
  int out_taures = firstzero_ac(res);
  return out_taures;
}
