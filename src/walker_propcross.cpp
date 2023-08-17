#include <Rcpp.h>
using namespace Rcpp;

//' Calculate walker_propcross feature
//'
//' @param y a numerical input vector
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- walker_propcross(y)
// [[Rcpp::export]]
double walker_propcross(NumericVector y) {
  int N = y.size();
  double p = 0.1;

  NumericVector w(N);
  w[0] = 0;

  for (int i = 1; i < N; ++i) {
    w[i] = w[i - 1] + p * (y[i - 1] - w[i - 1]);
  }

  double propcross = 0.0;

  for (int i = 0; i < N - 1; ++i) {
    if ((w[i] - y[i]) * (w[i + 1] - y[i + 1]) < 0) {
      propcross += 1.0;
    }
  }

  return propcross / (N - 1);
}
