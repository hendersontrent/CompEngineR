#include <Rcpp.h>
using namespace Rcpp;

//' Calculate mode of a data vector from software package \code{hctsa}
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
//' histogram_mode(x)
//'
// [[Rcpp::export]]
double histogram_mode(NumericVector x, int numBins = 10) {
  // Find the range of the data

  double x_min = x[0];
  double x_max = x[0];
  for (int i = 1; i < x.size(); ++i) {
    if (x[i] < x_min) x_min = x[i];
    if (x[i] > x_max) x_max = x[i];
  }

  // Compute the bin width

  double bin_width = (x_max - x_min) / numBins;

  // Initialize the counts and breaks vectors

  NumericVector counts(numBins);
  NumericVector breaks(numBins + 1);
  NumericVector midpoints(numBins);

  // Compute the bin midpoints

  for (int i = 0; i < numBins; ++i) {
    breaks[i] = x_min + i * bin_width;
    midpoints[i] = (breaks[i] + breaks[i + 1]) / 2;
  }

  breaks[numBins] = x_max;

  // Find the positions of the maximums

  IntegerVector maximum_positions;
  for (int i = 0; i < counts.size(); ++i) {
    if (counts[i] == max(counts)) {
      maximum_positions.push_back(i);
    }
  }

  // Return mean in case of duplicates

  double out = Rcpp::mean(maximum_positions);
  return out;
}
