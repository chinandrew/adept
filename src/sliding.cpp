#include <Rcpp.h>
using namespace Rcpp;

// From https://dirk.eddelbuettel.com/code/rcpp/Rcpp-attributes.pdf
// [[Rcpp::export]]
NumericVector convolve_cpp(const NumericVector a, const NumericVector b) {
  int na = a.size(), nb = b.size();
  int nab = na + nb - 1;
  NumericVector xab(nab);
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++)
      xab[i + j] += a[i] * b[j];
  return xab;
}

// Optimized version of
// https://github.com/cran/dvmisc/blob/master/src/sliding_cor_c.cpp
// [[Rcpp::export]]
NumericVector sliding_cor_cpp(const NumericVector shortvec,
                              const NumericVector longvec, double sd_shortvec) {

  // Get vector lengths and initialize output vector
  int length_longvec = longvec.size();
  int n = shortvec.size();
  int n_minus1 = n - 1;
  int out_length = length_longvec - n_minus1;
  NumericVector out(out_length);

  // Calculate sum of short vector divided by n nsquared
  double term2 = sum(shortvec) / n / n_minus1;

  // Loop through longvec. For each segment calculate the sum of the products
  // with shortvec and record the covariance
  NumericVector longvec_current(n);
  longvec_current = longvec[Range(0, n_minus1)];
  double sum_longvec_current = sum(longvec_current);
  double mean_longvec_current = sum_longvec_current / n;
  double sum_products = 0;
  double ss_longvec_current = 0;
  for (int b = 0; b < n; ++b) {
    double longvec_current_b = longvec[b];
    sum_products += longvec_current_b * shortvec[b];
    ss_longvec_current += pow(longvec_current_b - mean_longvec_current, 2);
  }
  double sd_longvec_current = sqrt(ss_longvec_current / n_minus1);
  out[0] = (sum_products / n_minus1 - sum_longvec_current * term2) /
    sd_shortvec / sd_longvec_current;

  for (int a = 1; a < out_length; ++a) {
    sum_longvec_current -= longvec[a - 1];
    sum_longvec_current += longvec[a + n_minus1];
    mean_longvec_current = sum_longvec_current / n;
    sum_products = 0;
    ss_longvec_current = 0;
    for (int b = 0; b < n; ++b) {
      double longvec_current_b = longvec[a + b];
      sum_products += longvec_current_b * shortvec[b];
      ss_longvec_current += pow(longvec_current_b - mean_longvec_current, 2);
    }
    sd_longvec_current = sqrt(ss_longvec_current / n_minus1);
    out[a] = (sum_products / n_minus1 - sum_longvec_current * term2) /
      sd_shortvec / sd_longvec_current;
  }
  return out;
}
