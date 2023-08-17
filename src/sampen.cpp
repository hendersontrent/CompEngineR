#include <Rcpp.h>
using namespace Rcpp;

//' Calculate sampenc feature
//'
//' @param y a numerical input vector
//' @param M an integer scalar denoting the embedding dimension
//' @param r a numeric scalar denoting the threshold
//' @return numeric scalar value for the feature
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- sampenc(y)
//'
// [[Rcpp::export]]
double sampenc(NumericVector y, int M = 6, double r = 0.3) {
  int N = y.size();
  NumericVector lastrun(N);
  NumericVector run(N);
  NumericVector A(M);
  NumericVector B(M);

  for (int i = 0; i < N - 1; ++i) {
    double y1 = y[i];

    for (int jj = 0; jj < N - i; ++jj) {
      int j = i + jj;
      if (std::abs(y[j] - y1) < r) {
        run[jj] = lastrun[jj] + 1;
        int M1 = std::min(M, static_cast<int>(run[jj]));

        for (int k = 0; k < M1; ++k) {
          A[k] += 1;
          if (j < N) B[k] += 1;
        }
      } else {
        run[jj] = 0;
      }
    }

    for (int j = 0; j < N - i; ++j) {
      lastrun[j] = run[j];
    }
  }

  double p = A[1] / B[0];
  double e = -std::log(p);
  return e;
}

//' Calculate sampen_first feature
//'
//' @param y a numerical input vector
//' @return numeric scalar value for the feature
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- sampen_first(y)
//'
// [[Rcpp::export]]
double sampen_first(NumericVector y) {
  int M = 5;
  double r = 0.3;
  double sampEn = sampenc(y, M + 1, r);
  return sampEn;
}
