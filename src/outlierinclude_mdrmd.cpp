#include <Rcpp.h>
using namespace Rcpp;

//' Calculate how the median depends on distributional outliers from software package \code{hctsa}
//'
//' @importFrom stats sd median
//' @param x \code{numeric} vector
//' @param zscored \code{Boolean} whether the data should be z-scored before computing the statistics. Defaults to \code{true}
//' @return \code{numeric} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' outlierinclude_mdrmd(x)
//'
// [[Rcpp::export]]
double outlierinclude_mdrmd(NumericVector x, bool zscored = true) {
  if (unique(x).size() == 1) {
    stop("The time series is a constant!");
  }
  double isd;
  if (zscored) {
    NumericVector tmp = ts(c(scale(x)));
    tsp(tmp) = tsp(x);
    x = tmp;
    isd = 1;
  } else {
    isd = sd(x, na_rm = true);
  }
  int N = x.size();
  double inc = 0.01 * isd;
  NumericVector thr = seq(0, max(abs(x), na_rm = true), inc);
  int tot = N;
  if (thr.size() == 0) stop("Peculiar time series");

  NumericVector msDt(thr.size());
  NumericVector msDtp(thr.size());

  for (int i = 0; i < thr.size(); i++) {
    double th = thr[i];
    IntegerVector r = which(abs(x) >= th);

    NumericVector Dt_exc = diff(r);
    msDt[i] = median(r) / (N / 2) - 1;
    msDtp[i] = Dt_exc.size() / tot * 100;
  }

  int trimthr = 2;
  IntegerVector mj = which(msDtp > trimthr)[msDtp.size() - 1];

  if (mj.size() != 0) {
    msDt = msDt[Range(0, mj[0])];
    msDtp = msDtp[Range(0, mj[0])];
    thr = thr[Range(0, mj[0])];
  } else {
    stop("Statistical power is lacking: less than 2% of data included");
  }

  double out_mdrmd = median(msDt);
  return out_mdrmd;
}
