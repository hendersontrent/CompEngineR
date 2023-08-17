#include <Rcpp.h>
using namespace Rcpp;

//' Calculate std1st_der feature
//'
//' @param y a numerical input vector
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- std1st_der(y)
// [[Rcpp::export]]
double std1st_der(NumericVector y) {
  int n = y.size();
  if (n < 2) {
    stop("Time series is too short to compute differences");
  }

  NumericVector yd(n - 1);
  for (int i = 1; i < n; i++) {
    yd[i - 1] = y[i] - y[i - 1];
  }

  double sum = 0.0;

  for (int i = 0; i < n - 1; i++) {
    sum += yd[i];
  }

  double mean = sum / (n - 1);
  double squared_diff_sum = 0.0;

  for (int i = 0; i < n - 1; i++) {
    squared_diff_sum += (yd[i] - mean) * (yd[i] - mean);
  }

  double std_dev = sqrt(squared_diff_sum / (n - 2));
  return std_dev;
}
