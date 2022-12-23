#' Wrapper function for just autocorr_features
#' @importFrom stats acf
#' @param x \code{numeric} vector
#' @return \code{data.frame} containing feature names and values
#' @author Trent Henderson
#' @export
#' @examples
#' x <- rnorm(100)
#' outs <- autocorr_features(x)
#'

autocorr_features <- function(x){

  if(anyNA(x) || !is.numeric(x)){
    stop("Input time series vector x should not have any missing or non-numeric values.")
  }

  acfv <- as.vector(stats::acf(x, length(x) - 1, plot = FALSE, na.action = na.pass)$acf)

  outs <- data.frame(names = c("embed2_incircle_1", "embed2_incircle_2", "ac_9", "firstmin_ac",
                               "trev_num", "motiftwo_entro3", "walker_propcross"),
                     values = c(embed2_incircle(x, 1, acfv = acfv),
                                embed2_incircle(x, 2, acfv = acfv),
                                ac_9(x, acfv),
                                firstmin_ac(x, acfv),
                                trev_num(x),
                                motiftwo_entro3(x),
                                walker_propcross(x)))

  return(outs)
}


#' Wrapper function for just pred_features
#' @param x \code{numeric} vector
#' @return \code{data.frame} containing feature names and values
#' @author Trent Henderson
#' @export
#' @examples
#' x <- rnorm(100)
#' outs <- pred_features(x)
#'

pred_features <- function(x){

  if(anyNA(x) || !is.numeric(x)){
    stop("Input time series vector x should not have any missing or non-numeric values.")
  }

  outs <- data.frame(names = c("localsimple_mean1", "localsimple_lfitac", "sampen_first"),
                     values = c(localsimple_taures(x, "mean"),
                                localsimple_taures(x, "lfit"),
                                sampen_first(x)))

  return(outs)
}


#' Wrapper function for just station_features
#' @param x \code{numeric} vector
#' @return \code{data.frame} containing feature names and values
#' @author Trent Henderson
#' @export
#' @examples
#' x <- rnorm(100)
#' outs <- station_features(x)
#'

station_features <- function(x){

  if(anyNA(x) || !is.numeric(x)){
    stop("Input time series vector x should not have any missing or non-numeric values.")
  }

  outs <- data.frame(names = c("std1st_der", "spreadrandomlocal_meantaul_50", "spreadrandomlocal_meantaul_ac2"),
                     values = c(std1st_der(x),
                                spreadrandomlocal_meantaul(x, 50),
                                spreadrandomlocal_meantaul(x, "ac2")))

  return(outs)
}


#' Wrapper function for just dist_features
#' @param x \code{numeric} vector
#' @return \code{data.frame} containing feature names and values
#' @author Trent Henderson
#' @export
#' @examples
#' x <- rnorm(100)
#' outs <- dist_features(x)
#'

dist_features <- function(x){

  if(anyNA(x) || !is.numeric(x)){
    stop("Input time series vector x should not have any missing or non-numeric values.")
  }

  outs <- data.frame(names = c("histogram_mode_10", "outlierinclude_mdrmd"),
                     values = c(histogram_mode(x), outlierinclude_mdrmd(x)))

  return(outs)
}


#' Wrapper function for just scal_features
#' @importFrom Rcatch22 SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1
#' @param x \code{numeric} vector
#' @return \code{data.frame} containing feature names and values
#' @author Trent Henderson
#' @export
#' @examples
#' x <- rnorm(100)
#' outs <- scal_features(x)
#'

scal_features <- function(x){

  if(anyNA(x) || !is.numeric(x)){
    stop("Input time series vector x should not have any missing or non-numeric values.")
  }

  outs <- data.frame(names = c("fluctanal_prop_r1"),
                     values = c(Rcatch22::SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(x)))

  return(outs)
}
