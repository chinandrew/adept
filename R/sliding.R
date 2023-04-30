#' @noRd
#'
sliding_cov_fast <- function(short, long) {
  n <- length(short)
  len_diff <- length(long) - n
  return(convolve_cpp(long, rev(short / (n - 1) - sum(short) / n / (n - 1)))[n:(n + len_diff)])
}

#' @noRd
#'
sliding_cor_fast <- function(short, long) {
  return(sliding_cor_cpp(short, long, sd(short)))
}

