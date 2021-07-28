
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector accessLOL(List L) {
  double rsum;
  NumericMatrix l1 = as<NumericMatrix>(L["one"]);
  NumericVector l2 = as<NumericVector>(L["two"]);
  int n = l2.length();
  NumericVector Nt(n);
  for (int i = 0; i < n; i++) {
    rsum = 0.0;
    for (int j = 0; j < n; j++) {
      rsum = rsum + (l1(i,j) * l2(j));
    }
    Nt(i) = rsum;
  }
  return Nt;
}

etest <- cxxfunction(signature(tm="NumericMatrix",tm2="NumericMatrix"),
                               plugin="RcppEigen",
                               body="
         NumericMatrix tm22(tm2);
         NumericMatrix tmm(tm);

                       const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
                       const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd> >(tm22));

                       Eigen::MatrixXd prod = ttm*ttm2;
                       return(wrap(prod));
                       ")

  set.seed(123)
  M1 <- matrix(sample(1e3),ncol=50)
  M2 <- matrix(sample(1e3),nrow=50)


