#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace RcppEigen;

// [[Rcpp::export]]
NumericVector mvC(NumericMatrix x, NumericVector y) {

   MatrixXd G = RccpEigen::as<MatrixXd>(x);
   VectorXd N = RcppEigen::as<VectorXd>(y);

   Eigen::MatrixXd prod = G*N;
   return(wrap(prod));

}
