#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::NumericVector fastLmResid(Rcpp::IntegerMatrix X, Rcpp::NumericVector y);
Rcpp::NumericVector fastLmResid2(Rcpp::NumericMatrix X, Rcpp::NumericVector y);
