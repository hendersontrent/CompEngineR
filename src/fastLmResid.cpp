#include <RcppArmadillo.h>
#include "fastLmResid.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate coefficients and their standard error during linear regression and return the residuals as per \code{RcppArmadillo}
 //'
 //' @param x \code{numeric} vector
 //' @return \code{numeric} scalar
 //' @references Dirk Eddelbuettel, Conrad Sanderson (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Computational Statistics and Data Analysis, Volume 71, March 2014, pages 1054-1063. URL http://dx.doi.org/10.1016/j.csda.2013.02.005
 //' @author Trent Henderson
 //'
 // [[Rcpp::export]]
 Rcpp::NumericVector fastLmResid2(Rcpp::NumericMatrix X, Rcpp::NumericVector y) {
   // Transform into RcppArmadillo objects
   arma::mat X2 = as<arma::mat>(wrap(X));
   arma::colvec y2 = as<arma::colvec>(wrap(y));
   // Dimension information
   int n = X2.n_rows, p = X2.n_cols;
   // Fit model y ~ X
   arma::colvec coef = arma::solve(X2, y2);
   // Compute the residuals
   arma::vec res = y - X*coef;
   Rcpp::NumericVector res2 = as<NumericVector>(wrap(res));
   return res2
 }
