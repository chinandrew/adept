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
// that requires having the longvec standar deviations passed in.
// Assumes shortvec is mean 0 with sd 1.
// [[Rcpp::export]]
NumericVector sliding_cor_sd_cpp(const NumericVector shortvec,
                                      const NumericVector longvec,
                                      const NumericVector sd_longvec_current) {

  // Get vector lengths and initialize output vector
  int length_longvec = longvec.size();
  int n = shortvec.size();
  int n_minus1 = n - 1;
  int out_length = length_longvec - n_minus1;
  NumericVector out(out_length);

  // Loop through longvec. For each segment calculate the sum of the products
  // with shortvec and record the covariance
  NumericVector longvec_current(n);
  longvec_current = longvec[Range(0, n_minus1)];
  double sum_longvec_current = sum(longvec_current);
  double sum_products = 0;
  for (int b = 0; b < n; ++b) {
    double longvec_current_b = longvec[b];
    sum_products += longvec_current_b * shortvec[b];
  }
  out[0] = (sum_products / n_minus1) / sd_longvec_current[0];

  for (int a = 1; a < out_length; ++a) {
    sum_longvec_current -= longvec[a-1];
    sum_longvec_current += longvec[a+n_minus1];
    sum_products = 0;
    for (int b = 0; b < n; ++b) {
      double longvec_current_b = longvec[a+b];
      sum_products += longvec_current_b * shortvec[b];
    }
    out[a] = (sum_products / n_minus1) / sd_longvec_current[a];
  }
  return out;
}

// Optimized version of
// https://github.com/cran/dvmisc/blob/master/src/sliding_cor_c.cpp
// that also returns the longvec standard deviations for reuse in
// sliding_cor_sd_cpp().
// Assumes shortvec is mean 0 with sd 1.
// [[Rcpp::export]]
List sliding_cor_store_sd_cpp(const NumericVector shortvec, const NumericVector longvec) {

  // Get vector lengths and initialize output vector
  int length_longvec = longvec.size();
  int n = shortvec.size();
  int n_minus1 = n - 1;
  int out_length = length_longvec - n_minus1;
  NumericVector out(out_length);
  NumericVector sds(out_length);

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
  if (sd_longvec_current < 1e-10){
    out[0] =  NA_REAL;
    sds[0] =  NA_REAL; // could also set to 0 to get Infs downstream
  } else{
    out[0] = (sum_products / n_minus1) / sd_longvec_current;
    sds[0] = sd_longvec_current;
  }
  for (int a = 1; a < out_length; ++a) {
    sum_longvec_current -= longvec[a-1];
    sum_longvec_current += longvec[a+n_minus1];
    mean_longvec_current = sum_longvec_current / n;
    sum_products = 0;
    ss_longvec_current = 0;
    for (int b = 0; b < n; ++b) {
      double longvec_current_b = longvec[a+b];
      sum_products += longvec_current_b * shortvec[b];
      ss_longvec_current += pow(longvec_current_b - mean_longvec_current, 2);
    }
    if (sd_longvec_current < 1e-10){
      out[a] =  NA_REAL;
      sds[a] =  NA_REAL;
    } else{
      out[a] = (sum_products / n_minus1) / sd_longvec_current;
      sds[a] = sd_longvec_current;
    }
    sd_longvec_current = sqrt(ss_longvec_current / n_minus1);

  }
  return List::create(_["cor"] = out, _["sds"] = sds);
}
