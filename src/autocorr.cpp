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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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

//' Calculate embed2_incircle feature
//'
//' @param y a numerical input vector
//' @param boundary an integer denoting the given circular boundary, setting to 1 or 2 in CompEngine. Default to 1
//' @param acfv vector of autocorrelation, if exist, used to avoid repeated computation
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- embed2_incircle(y)
// [[Rcpp::export]]
double embed2_incircle(NumericVector y, double boundary = NA_REAL) {
  if(R_IsNA(boundary)) {
    Rcpp::warning("`embed2_incircleRcpp()` using `boundary = 1`. Set value with `boundary`.");
    boundary = 1;
  }

  NumericVector acfv = ac(y);
  int tau = firstzero_ac(acfv);
  NumericVector xt = y[Rcpp::Range(0, y.size() - tau - 1)];
  NumericVector xtp = y[Rcpp::Range(tau, y.size() - 1)];
  int N = y.size() - tau;
  double sumInside = 0.0;

  for (int i = 0; i < N; ++i) {
    if (xt[i] * xt[i] + xtp[i] * xtp[i] < boundary) {
      sumInside += 1.0;
    }
  }

  return sumInside / N;
}

//' Calculate ac_9 feature
//'
//' @param y a numerical input vector
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- ac_9(y)
// [[Rcpp::export]]
double ac_9(NumericVector y) {
  NumericVector acfv = ac(y);
  return acfv[10];
}

//' Calculate firstmin_ac feature
//'
//' @param y a numerical input vector
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- firstmin_ac(y)
// [[Rcpp::export]]
int firstmin_ac(NumericVector y) {
  int N = y.size();
  NumericVector acfv = ac(y);
  NumericVector autoCorr = acfv;

  for (int i = 1; i < N - 2; ++i) {
    if (NumericVector::is_na(autoCorr[i])) {
      Rcpp::warning("No minimum was found.");
      return NA_INTEGER;
    }
    if (i == 1 && autoCorr[1] > autoCorr[0]) {
      return 1;
    } else if (i > 1 && autoCorr[i - 2] > autoCorr[i - 1] && autoCorr[i - 1] < autoCorr[i]) {
      return i - 1;
    }
  }

  return N - 1;
}

//' Calculate trev_num feature
//'
//' @param y a numerical input vector
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- trev_num(y)
// [[Rcpp::export]]
double trev_num(NumericVector y) {
  int n = y.size();
  NumericVector yn = y[Range(0, n - 2)];
  NumericVector yn1 = y[Range(1, n - 1)];

  double sum = 0.0;
  int count = 0;

  for (int i = 0; i < n - 1; ++i) {
    if (!NumericVector::is_na(yn[i]) && !NumericVector::is_na(yn1[i])) {
      sum += pow(yn1[i] - yn[i], 3);
      count++;
    }
  }

  if (count > 0) {
    return sum / count;
  } else {
    return NA_REAL;
  }
}
