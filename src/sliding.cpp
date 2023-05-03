#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
List update_pmax_max_cpp(NumericVector new_vec, int new_idx, NumericVector old_vec, NumericVector old_idxs) {
  // int max_i = 0;
  NumericVector out = clone(old_vec);
  NumericVector out_idx = clone(old_idxs);
  for (int i = 0; i < out.length(); ++i) {
    if (new_vec[i] > old_vec[i]){
      out_idx[i] = new_idx;
      out[i] = new_vec[i];
    }
  }
  return List::create(_["pmax"]=out,
                      _["idx"]=out_idx);
}


// Adapted from https://stackoverflow.com/a/66020829,
// Vectors should be of the same length.
// [[Rcpp::export]]
List pmax_max_cpp(List args) {
  // int max_i = 0;
  int max_j = 0;
  NumericVector tmp = args[0];
  NumericVector out = clone(tmp);
  NumericVector out_idx (out.length(), 1);
  double current_max = out[0];
  int n_arg = args.length();
  int n_vec = out.length();
  for (int i = 1; i < n_arg; ++i) {
    NumericVector pa = args[i];
    for (int j = 0; j < n_vec; ++j) {
      if (pa[j] > out[j]) {
        out[j] = pa[j];
        out_idx[j] = i + 1;
      }
    }
  }
  return List::create(_["pmax"]=out,
                      _["idx"]=out_idx);
}


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
  double test = 0;
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
