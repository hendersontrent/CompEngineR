#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// Define a sequence generation helper function
Rcpp::NumericVector get_sequence(double start, double end, int length) {
  Rcpp::NumericVector out(length);
  double step = (end - start) / (length - 1);
  for (int i = 0; i < length; i++) {
    out[i] = start + i * step;
  }
  return out;
}

// Define a rounding helper function
Rcpp::IntegerVector round_vector(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::IntegerVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = std::round(x[i]);
  }
  return out;
}

// Define a unique finder helper function
Rcpp::IntegerVector unique_elements(Rcpp::IntegerVector x) {
  std::sort(x.begin(), x.end());  // sort the input vector
  auto it = std::unique(x.begin(), x.end());  // remove duplicates
  x.erase(it, x.end());  // erase the removed elements
  return x;  // return the modified vector
}

// Define C++ version of `split`
Rcpp::List split(Rcpp::NumericVector x, Rcpp::NumericVector f) {
  int n = x.size();
  int m = f.size();
  Rcpp::List out(m + 1);  // create a list with m + 1 elements
  int start = 0;
  for (int i = 0; i < m; i++) {
    auto it = std::find(x.begin() + start, x.end(), f[i]);  // find the next split point
    if (it == x.end()) {  // if the split point is not found
      out[i] = x[Rcpp::Range(start, n - 1)];  // copy the remaining values to the output list
      break;  // exit the loop
    } else {  // if the split point is found
      int end = std::distance(x.begin(), it);  // index of the split point
      out[i] = x[Rcpp::Range(start, end - 1)];  // copy the values up to the split point
      start = end;  // update the start index
    }
  }
  out[m] = x[Rcpp::Range(start, n - 1)];  // copy the remaining values to the last element of the output list
  return out;
}

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
  IntegerVector taur = unique_elements(round_vector(exp(get_sequence(log(5), log(floor(N / 2)), tauStep))));
  int ntau = taur.size();
  if (ntau < 8) {
    stop("This time series is too short to analyse using this fluctuation analysis");
  }

  NumericVector Fl(ntau);

  for (int i = 0; i < ntau; i++) {
    int tau = taur[i];
    List y_buff = List::create(split(y, ceiling(seq_along(y) / tau)));

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
