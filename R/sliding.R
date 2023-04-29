#' @noRd
#'
sliding_cov2 <- function(short, long) {
  n = length(short)
  len_diff = length(long) - n
  convolve_cpp(long, rev(short / (n - 1) - sum(short) / n / (n - 1)))[n:(n + len_diff)]
}

#' @noRd
#'
sliding_cor2 <- function(short, long) {
  sliding_cor_cpp(short, long, sd(short))
}

