#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title invC matrix inversion using RcppArmadillo
//'
//' @description invC uses BLAS to conduct a matrix inversion and
//'     appears to be about 73 percent faster than R's solve when used
//'     in a length-based model context.
//'
//' @param x the matrix to be inverted
//'
//' @export
// [[Rcpp::export("invC")]]
arma::mat invC(const arma::mat& x) {
  arma::mat y;
  y=arma::inv(x);
  return(y);
}

//' @title mvC matrix vector product using RcppArmadillo
//'
//' @description mvC multiple the vector rhs by the matrix lhs, which
//'     is used all the time in a length-based stock assessment model.
//'     mvC appears to be 25 percent faster than R's inmat X invect
//'     when used in a length-based model context.
//'
//' @param inmat the matrix, usualy the growth transition matrix
//' @param invect the vector to be multiplied. Usually the numbers-at-
//'     size in a length based model
//' @export
// [[Rcpp::export("mvC")]]
arma::vec mvC(const arma::mat& inmat, const arma::vec& invect) {
  return inmat * invect;
}

//' @title svvC sum of two multiplied vectors using RcppArmadillo
//'
//' @description svvC multiple the vector rhs by the vector lhs and
//'     sums the resulting vector.
//'
//' @param lhs the first vector, possibly numbers at size
//' @param rhs the second vector, possibly maturity x weight-at-length
//' @export
// [[Rcpp::export("svvC")]]
double svvC(const arma::vec lhs, const arma::vec rhs){
  double ans = arma::as_scalar(lhs.t() * rhs);
  return ans;
}
