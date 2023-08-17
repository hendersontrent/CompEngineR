#' Calculate every time-series feature in the package
#'
#' @param y \code{numeric} vector
#' @return \code{data.frame} containing feature names and values
#' @author Trent Henderson
#' @export
#' @examples
#' y <- rnorm(100)
#' outs <- comp_engine(y)
#'

comp_engine <- function(y){

  if(anyNA(y) || !is.numeric(y)){
    stop("Input time series vector x should not have any missing or non-numeric values.")
  }

  autocor_fs <- autocorr_features(y)
  pred_fs <- pred_features(y)
  station_fs <- station_features(y)
  dist_fs <- dist_features(y)
  scal_fs <- scal_features(y)
  outs <- rbind(autocor_fs, pred_fs, station_fs, dist_fs, scal_fs)
  return(outs)
}
