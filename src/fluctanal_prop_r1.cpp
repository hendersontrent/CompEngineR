#include <Rcpp.h>
using namespace Rcpp;

//' Calculate fluctuation analysis from software package \code{hctsa}
//'
//' @importFrom stats lm residuals
//' @param x \code{numeric} vector
//' @return \code{numeric} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' fluctanal_prop_r1(x)
//'
// [[Rcpp::export]]
double fluctanal_prop_r1(NumericVector x) {
  int q = 2;
  int tauStep = 50;
  int k = 1;

  int N = x.size();
  NumericVector x_NA0 = ifelse(!is_na(x), x, 0);

  NumericVector y = cumsum(x_NA0);
  NumericVector taur = unique(round(exp(seq(log(5), log(floor(N / 2)), length = tauStep))));
  int ntau = taur.size();
  if (ntau < 8) {
    stop("This time series is too short to analyse using this fluctuation analysis");
  }

  NumericVector Fl(ntau);

  for (int i = 0; i < ntau; i++) {
    int tau = taur[i];
    NumericVector y_buff = split(y, ceiling(seq_along(y) / tau));

    if (y_buff.size() > floor(y.size() / tau)) {
      y_buff = y_buff[Range(0, y_buff.size() - 2)];
    }

    int nn = y_buff.size() * tau;
    NumericVector tt = seq(1, tau);

    for (int j = 0; j < y_buff.size(); j++) {
      DataFrame df = DataFrame::create(Named("tt")=tt, Named("lmy")=y_buff[j]);
      SEXP lm_tt = lm(lmy ~ tt, data = df);
      y_buff[j] = residuals(lm_tt);
    }

    NumericMatrix tem = sapply(y_buff, range);
    NumericVector y_dt = tem(_, 1) - tem(_, 0);

    Fl[i] = pow(mean(pow(y_dt, k)), 1.0 / k);
  }
  NumericVector logtt = log(taur);
  NumericVector logFF = log(Fl);
  int ntt = ntau;

  NumericVector sserr(ntt, NA_REAL);
  int minPoints = 6;

  for (int i = minPoints; i < ntt - minPoints; i++) {
    NumericVector r1 = seq(1, i);
    DataFrame df1 = DataFrame::create(Named("x")=logtt[r1-1], Named("y")=logFF[r1-1]);
    SEXP p1 = lm(y ~ x, data = df1);
    NumericVector r2 = seq(i, ntt);
    DataFrame df2 = DataFrame::create(Named("x")=logtt[r2-1], Named("y")=logFF[r2-1]);
    SEXP p2 = lm(y ~ x, data = df2);
    sserr[i] = sum(abs(residuals(p1)) + abs(residuals(p2)));
  }

  int breakPt = which_min(sserr) + 1;
  NumericVector r1 = seq(1, breakPt);
  NumericVector r2 = seq(breakPt, ntt);
  double prop_r1 = r1.size() / (double) ntt;
  return prop_r1;
}
