#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Element-wise Soft Thresholding Function for Matrix
//' @description Soft-thresholding function for a matrix which is the solution of
//' \deqn{\argmin_Y ||Y-X||_2^2/2+\lambda ||Y||_1.}
//' @param \code{X}, a Matrix
//' @param \code{lambda}, a scalar
//' @return Matrix after threholding
//' @export
// [[Rcpp::export]]
arma::mat soft(arma::mat X, double lambda=0){
  lambda = std::abs(lambda);
  arma::mat Y = (X>=lambda)%(X-lambda)+(X<=(-lambda))%(X+lambda);
  return Y;
}

//' @title standard Lasso estimator of the projection matrix
//' @description ADMM algorithm for high dimensional projection matrix with lasso
//' \deqn{\argmin_W tr(\frac{1}{2}W^{\top}\Sigma_z W-\Sigma_{zx}^{\top}W)+\lambda ||W||_1.}
//' @param \code{Z}, a \eqn{n*p_z} input data matrix
//' @param \code{covzx}, a \eqn{p_z*p_x} covariance matrix between \code{Z} and \code{X}.
//' @param \code{err_abs}, the absolute tolerance precision used to stop the convergence of ADMM.
//' @param \code{maxIter}, Maximum number of iterations. Default is 1000.
//' @param \code{lambda}, a scalar
//' @param \code{rho}, initial step parameter for ADMM.
//' @return A list with components
//' \item{Omega}{a list of sparse \eqn{p_z*p_x} matrices corresponding to lambda.}
//' \item{lambda}{the used lambda for the solution path.}
//' \item{niter}{the number of each iterations for each number of lambda.}
//' @export
// [[Rcpp::export]]
Rcpp::List projlasso(arma::mat Z, arma::mat covzx, arma::vec lambda, double err_abs=10^(-4), int maxIter = 1000, double rho=1){
  int n=Z.n_rows;
  int px = covzx.n_cols;
  int pz = Z.n_cols;
  int m = (pz<n)*pz+(pz>=n)*n;
  int nlambda = lambda.size();
  
  /*Centering int*/
  arma::mat dZ = (arma::eye(n,n)-arma::ones(n,n)/n)*Z/sqrt(n-1);
  arma::mat U;
  arma::mat V;
  arma::vec Tau;
  svd(U, Tau, V, dZ.t());
  U = U.cols(0,m-1);
  Tau = Tau%Tau;
  arma::mat Tau1 = arma::diagmat(Tau/(Tau+rho));
  arma::mat D = U*Tau1*U.t();
  
  Rcpp::List Ome_all(nlambda);
  arma::vec niter = lambda;
  
  /*initialization*/
  arma::mat A = covzx;
  arma::mat B = arma::zeros(pz,px);
  arma::mat Ome;
  arma::mat A1;
  arma::mat C;
  double lamb;
  double error=1;
  
  /*iterations*/
  for(int j=0; j<nlambda; ++j){
    lamb = lambda(j);
    int k=0;
    while(((k<maxIter)&&(error>err_abs))||(k==0)){
      k = k+1;
      A1 = A;
      C = covzx+rho*(A-B);
      Ome = 1/rho*C-1/rho*D*C;
      A = soft(Ome+B, lamb/rho);
      B = Ome-A+B;
      error = mean(mean(abs(A-A1)));
    }
    Ome_all(j) = A;
    niter(j) = k;
  }
  return Rcpp::List::create(Rcpp::Named("Ome_all") = Ome_all,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("lambda") = lambda);
}












