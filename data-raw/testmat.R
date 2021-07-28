library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(microbenchmark)


#sourceCpp("C:/Users/Malcolm/Dropbox/rcode2/testRcpp/meanC.cpp")


sourceCpp("C:/Users/User/Dropbox/A_Code/aMSE/data-raw/Rcpp_functions.cpp")


one <- matrix(rnorm(25),nrow=5,ncol=5)
two <- rnorm(5)



x <- list(one=one,two=two)



accessLOL(x)


a <- c(0.0,1.5,2.0,1.5,0.0,
       0.8,0,0,0,0,
       0,0.65,0,0,0,
       0,0,0.45,0,0,
       0,0,0,0.2,0)
mat <- matrix(a,nrow=5,ncol=5,byrow=TRUE)
mat
vect <- c(100,0,0,0,0)

dat <- list(one=mat,two=vect)
accessLOL(dat)

microbenchmark(
  dat$one %*% dat$two,
  accessLOL(dat),
  times=100
)

sourceCpp("C:/Users/User/Dropbox/A_Code/aMSE/data-raw/testeigen.cpp")









library(inline)

code <- 'arma::mat x = Rcpp::as<arma::mat>(X);
         int n = as<int>(N);
         for (int i=0; i<n; i++) arma::mat y = x*x;
         return R_NilValue;
         '
code <- 'arma::mat x = Rcpp::as<arma::mat>(X);
         arma::vec y = Rcpp::as<arma::vec>(Y);
         arma::mat z = x*y;
         return R_NilValue;
         '



fun_Rcpp <- cxxfunction(signature(X="numeric", Y="numeric"),
                        body=code, plugin="RcppArmadillo")

fun_R <- function(x, y){
  return( x %*% y )
}  ## Simple setup

n <- 1000
p <- 5
x <- matrix(rnorm(p^2), p,p)
y <- x[1,]
microbenchmark(
  x %*% y,
  fun_R(x,y),
  times=100
)



