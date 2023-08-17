#include <Rcpp.h>
using namespace Rcpp;

// Calculate the mean of a numeric vector
double mean_(NumericVector y) {
  int n = y.size();
  double sum = 0.0;

  for (int i = 0; i < n; i++) {
    sum += y[i];
  }

  return sum / n;
}

IntegerVector binarize_mean(NumericVector y) {
  int n = y.size();
  double y_mean = mean(y);
  IntegerVector Y(n);

  for (int i = 0; i < n; ++i) {
    Y[i] = (y[i] - y_mean) > 0 ? 1 : 0;
  }

  return Y;
}

double f_entropy(NumericVector x) {
  int n = x.size();
  double entropy = 0.0;

  for (int i = 0; i < n; ++i) {
    if (x[i] > 0) {
      entropy -= x[i] * log(x[i]);
    }
  }

  return entropy;
}

//' Calculate motiftwo_entro3 feature
//'
//' @param y a numerical input vector
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- motiftwo_entro3(y)
// [[Rcpp::export]]
double motiftwo_entro3(NumericVector y) {
  IntegerVector yBin = binarize_mean(y);
  int N = yBin.size();

  if (N < 5) {
    Rcpp::warning("Time series too short");
  }

  LogicalVector r1 = yBin[yBin == 1];
  LogicalVector r0 = yBin[yBin == 0];

  int r1_size = r1.size();
  int r0_size = r0.size();

  LogicalVector r00 = (r0 * (yBin[Range(1, N - 1)] == 0))[Range(0, r0_size - 2)];
  LogicalVector r01 = (r0 * (yBin[Range(1, N - 1)] == 1))[Range(0, r0_size - 2)];
  LogicalVector r10 = (r1 * (yBin[Range(1, N - 1)] == 0))[Range(0, r1_size - 2)];
  LogicalVector r11 = (r1 * (yBin[Range(1, N - 1)] == 1))[Range(0, r1_size - 2)];

  int r00_size = r00.size();
  int r01_size = r01.size();
  int r10_size = r10.size();
  int r11_size = r11.size();

  LogicalVector r000 = (r00 * (yBin[Range(2, N - 1)] == 0))[Range(0, r00_size - 2)];
  LogicalVector r001 = (r00 * (yBin[Range(2, N - 1)] == 1))[Range(0, r00_size - 2)];
  LogicalVector r010 = (r01 * (yBin[Range(2, N - 1)] == 0))[Range(0, r01_size - 2)];
  LogicalVector r011 = (r01 * (yBin[Range(2, N - 1)] == 1))[Range(0, r01_size - 2)];
  LogicalVector r100 = (r10 * (yBin[Range(2, N - 1)] == 0))[Range(0, r10_size - 2)];
  LogicalVector r101 = (r10 * (yBin[Range(2, N - 1)] == 1))[Range(0, r10_size - 2)];
  LogicalVector r110 = (r11 * (yBin[Range(2, N - 1)] == 0))[Range(0, r11_size - 2)];
  LogicalVector r111 = (r11 * (yBin[Range(2, N - 1)] == 1))[Range(0, r11_size - 2)];

  double out_ddd = mean_(r000);
  double out_ddu = mean_(r001);
  double out_dud = mean_(r010);
  double out_duu = mean_(r011);
  double out_udd = mean_(r100);
  double out_udu = mean_(r101);
  double out_uud = mean_(r110);
  double out_uuu = mean_(r111);

  NumericVector ppp = NumericVector::create(out_ddd, out_ddu, out_dud, out_duu, out_udd, out_udu, out_uud, out_uuu);
  double out_hhh = f_entropy(ppp);
  return out_hhh;
}
