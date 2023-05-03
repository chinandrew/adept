#' @noRd
#'
sliding_cov_fast <- function(short, long) {
  n <- length(short)
  len_diff <- length(long) - n
  return(convolve_cpp(long, rev(short / (n - 1)))[n:(n + len_diff)])
}

#' @noRd
#'
sliding_cov_fft <- function(short, long) {
  n <- length(short)
  len_diff <- length(long) - n
  return(convolve(long, c(short / (n - 1), rep(0,len_diff)))[1:(len_diff+1)])
}


#' @noRd
#'
sliding_cor_store_sd <- function(short, long) {
  return(sliding_cor_store_sd_cpp(short, long))
}

#' @noRd
#'
sliding_cor_sd <- function(short, long, sds) {
  return(sliding_cor_sd_cpp(short, long, sds))
}
