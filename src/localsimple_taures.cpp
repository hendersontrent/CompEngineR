#include <Rcpp.h>
using namespace Rcpp;

NumericVector OLS(NumericVector x, NumericVector y) {
  int n = x.size();
  double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;

  for (int i = 0; i < n; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_xx += x[i] * x[i];
    sum_xy += x[i] * y[i];
  }

  double mean_x = sum_x / n;
  double mean_y = sum_y / n;

  double beta1 = (sum_xy - n * mean_x * mean_y) / (sum_xx - n * mean_x * mean_x);
  double beta0 = mean_y - beta1 * mean_x;

  NumericVector coefficients = NumericVector::create(beta0, beta1);
  return coefficients;
}

// Calculate the mean of a numeric vector
double mean_(NumericVector y) {
  int n = y.size();
  double sum = 0.0;

  for (int i = 0; i < n; i++) {
    sum += y[i];
  }

  return sum / n;
}

// Calculate the sum of a numeric vector
double sum_(NumericVector a) {
  double m = 0.0;
  int size = a.size();

  for (int i = 0; i < size; i++){
    m += a[i];
  }

  return m;
}

//' Calculate autocorrelation function up to a given lag
 //'
 //' @param y a numerical input vector
 //' @param lag_max integer denoting the largest lag to compute for. Defaults to 1
 //' @param demean Boolean whether to demean the time series. Defaults to true
 //' @author Trent Henderson
 //' @examples
 //' y <- stats::rnorm(100)
 //' outs <- ac(y, 1)
 NumericVector ac(NumericVector y, bool demean = true) {
   int n = y.size();
   NumericVector x = y;

   if (demean) {
     x = y - mean_(y);
   }

   NumericVector result(n - 1);

   for (int k = 0; k < n - 1; ++k) {
     double sum = 0;
     for (int i = 0; i < n - k; ++i) {
       sum += x[i] * x[i + k];
     }
     result[k] = sum / sum_(x * x);
   }

   return result;
 }

//' Find the first zero in ACF
 //'
 //' @param acfv a numerical input vector
 //' @author Trent Henderson
 //' @export
 //' @examples
 //' y <- stats::rnorm(100)
 //' outs <- firstzero_ac(y)
 int firstzero_ac(NumericVector y) {
   int N = y.size();
   NumericVector acfv = ac(y);
   NumericVector acfValues = acfv[Range(1, N - 1)];

   // Find indices where acfValues are less than 0
   IntegerVector tau;

   for (int i = 0; i < acfValues.size(); ++i) {
     if (acfValues[i] < 0) {
       tau.push_back(i + 1); // Adding 1 to convert from 0-based index to 1-based
     }
   }

   if (tau.size() == 0) {
     return 0;
   } else if (std::all_of(tau.begin(), tau.end(), [](int val) { return val == NA_INTEGER; })) {
     return 0;
   } else if (std::none_of(tau.begin(), tau.end(), [](int val) { return val != 0; })) {
     return N;
   } else {
     return tau[0];
   }
 }

//' Calculate localsimple_taures feature
//'
//' @param y a numerical input vector
//' @param forecastMeth character denoting the method to use. Can be one of "mean" or "lfit"
//' @param trainLength integer denoting the number of time-series values to use to forecast the next value. Defaults to 0 for auto-detection
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- localsimple_taures(y, "lfit")
// [[Rcpp::export]]
int localsimple_taures(NumericVector y, std::string forecastMeth = "mean", int trainLength = 0) {
  int lp = 1;

  if (trainLength == 0) {
    lp = firstzero_ac(y);
  }

  int N = y.size();
  NumericVector evalr(N - lp);

  for (int i = 0; i < evalr.size(); i++) {
    evalr[i] = lp + 1 + i;
  }

  if (lp >= N) {
    Rcpp::warning("Time series too short for forecasting in `localsimple_taures`");
    return NA_REAL;
  }

  NumericVector res(evalr.size());

  if (forecastMeth == "mean") {
    for (int i = 0; i < evalr.size(); i++) {
      double sum = 0.0;
      for (int j = evalr[i] - lp - 1; j < evalr[i] - 1; j++) {
        sum += y[j];
      }
      res[i] = sum / lp - y[evalr[i] - 1];
    }
  }

  if (forecastMeth == "lfit") {
    for (int i = 0; i < evalr.size(); i++) {

      NumericVector a(lp);

      for (int j = 0; j < lp; j++) {
        a[j] = j + 1;
      }

      NumericVector b(lp);

      for (int j = 0; j < lp; j++) {
        b[j] = y[evalr[i] - lp + j - 1];
      }

      NumericVector lm_ab = OLS(a, b);
      double prediction = lm_ab[0] + lm_ab[1] * (lp + 1);
      res[i] = prediction - y[evalr[i] - 1];
    }
  }

  int out = firstzero_ac(res);
  return out;
}
