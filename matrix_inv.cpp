// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// Define the function for matrix inversion using Eigen
// [[Rcpp::export]]
Eigen::MatrixXd inverseMatrix(Eigen::Map<Eigen::MatrixXd> A) {
  Eigen::MatrixXd A_inv = A.inverse();
  return A_inv;
}