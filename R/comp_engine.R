#' Calculate every time-series feature in the package
#' @param x \code{numeric} vector
#' @return \code{data.frame} containing feature names and values
#' @author Trent Henderson
#' @export
#' @examples
#' x <- rnorm(100)
#' outs <- comp_engine(x)
#'

comp_engine <- function(x){

  if(anyNA(x) || !is.numeric(x)){
    stop("Input time series vector x should not have any missing or non-numeric values.")
  }

  autocor_fs <- autocorr_features(x)
  pred_fs <- pred_features(x)
  station_fs <- station_features(x)
  dist_fs <- dist_features(x)
  scal_fs <- scal_features(x)
  outs <- rbind(autocor_fs, pred_fs, station_fs, dist_fs, scal_fs)
  return(outs)
}
