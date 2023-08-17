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

// Calculate the mean of an integer vector
double mean_int(IntegerVector y) {
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

//' Calculate spreadrandomlocal_meantaul feature
//'
//' @param y a numerical input vector
//' @param l integer denoting the length of local time-series segments to analyse as a positive integer. Defaults to 50
//' @author Trent Henderson
//' @export
//' @examples
//' y <- stats::rnorm(100)
//' outs <- spreadrandomlocal_meantaul(y)
// [[Rcpp::export]]
double spreadrandomlocal_meantaul(NumericVector y, int l = 50) {
   int numSegs = 100;
   int N = y.size();

   if (l > 0.9 * N) {
     Rcpp::warning("This time series is too short. Specify proper segment length in `l`");
     return NA_REAL;
   }

   IntegerVector qs(numSegs);

   for (int j = 0; j < numSegs; j++) {
     int ist = Rcpp::sample(N - 1 - l, 1)[0];
     int ifh = ist + l - 1;
     IntegerVector rs(ifh - ist);

     for (int k = 0; j < rs.size(); j++) {
       rs[k] = ist + k;
     }

     NumericVector ysub = y[rs];
     int taul = firstzero_ac(ysub);
     qs[j] = taul;
   }

   double mean_qs = mean_int(qs);
   return mean_qs;
 }
