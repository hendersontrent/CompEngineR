#include <Rcpp.h>
using namespace Rcpp;

//' Calculate local motifs in a binary symbolization of the time series from software package \code{hctsa}
//'
//' @param x \code{numeric} vector
//' @return \code{numeric} scalar
//' @references Hyndman R, Kang Y, Montero-Manso P, Talagala T, Wang E, Yang Y, O'Hara-Wild M (2022). _tsfeatures: Time Series Feature Extraction_. R package version 1.1, <https://CRAN.R-project.org/package=tsfeatures>.
//' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
//' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
//' @author Trent Henderson
//' @export
//' @examples
//' x <- rnorm(100)
//' motiftwo_entro3(x)
//'
// [[Rcpp::export]]
double motiftwo_entro3(NumericVector x) {

  NumericVector xBin = binarize_mean(x);
  int N = xBin.size();
  if (N < 5) warning("Time series is too short for motiftwo_entro3");

  // Create all possible combinations of three binary values

  LogicalVector r000 = xBin[Range(0, N - 3)] == 0 & xBin[Range(1, N - 2)] == 0 & xBin[Range(2, N - 1)] == 0;
  LogicalVector r001 = xBin[Range(0, N - 3)] == 0 & xBin[Range(1, N - 2)] == 0 & xBin[Range(2, N - 1)] == 1;
  LogicalVector r010 = xBin[Range(0, N - 3)] == 0 & xBin[Range(1, N - 2)] == 1 & xBin[Range(2, N - 1)] == 0;
  LogicalVector r011 = xBin[Range(0, N - 3)] == 0 & xBin[Range(1, N - 2)] == 1 & xBin[Range(2, N - 1)] == 1;
  LogicalVector r100 = xBin[Range(0, N - 3)] == 1 & xBin[Range(1, N - 2)] == 0 & xBin[Range(2, N - 1)] == 0;
  LogicalVector r101 = xBin[Range(0, N - 3)] == 1 & xBin[Range(1, N - 2)] == 0 & xBin[Range(2, N - 1)] == 1;
  LogicalVector r110 = xBin[Range(0, N - 3)] == 1 & xBin[Range(1, N - 2)] == 1 & xBin[Range(2, N - 1)] == 0;
  LogicalVector r111 = xBin[Range(0, N - 3)] == 1 & xBin[Range(1, N - 2)] == 1 & xBin[Range(2, N - 1)] == 1;

  // Calculate the mean of each combination

  double out_ddd = mean(r000);
  double out_ddu = mean(r001);
  double out_dud = mean(r010);
  double out_duu = mean(r011);
  double out_udd = mean(r100);
  double out_udu = mean(r101);
  double out_uud = mean(r110);
  double out_uuu = mean(r111);

  NumericVector ppp = {out_ddd, out_ddu, out_dud, out_duu, out_udd, out_udu, out_uud, out_uuu};
  double out_hhh = f_entropy(ppp);
  return out_hhh;
}
