#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::NumericMatrix fast_matrix_product(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B) {
  int nrow = A.nrow();
  int ncol = A.ncol();
  
  Rcpp::NumericMatrix C(nrow, ncol);
  
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      double sum = 0;
      for (int k = 0; k < ncol; ++k) {
        sum += A(i, k) * B(k, j);
      }
      C(i, j) = sum;
    }
  }
  
  return C;
}

