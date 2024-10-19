#' The test procedure for partial parameteres under the smoothed quantile regression
#' @param y response vector with length \eqn{n}.
#' @param X data matrix of dimension \eqn{n \times p_{x}}.
#' @param Z data matrix of dimension \eqn{n \times p_{z}}. 
#' @param tau the quantile level with value between (0,1). The default is 0.5.
#' @param kernel A character string specifying the choice of kernel function., the default is "Gaussian".
#' Choices are "Gaussian", "uniform" and "logistic"
#' @param h The bandwidth parameter for kernel smoothing. The default is \eqn{\max\{0.05,0.5*(\log p_z/n)^{0.3}\}}
#' @return A list with components
#' \item{tau}{The chosen quantile level}
#' \item{kernel}{The chosen kernel}
#' \item{statistic}{The value of the proposed U-statistic}
#' \item{standard.deviation}{The value of the standard deviation of the test statistic}
#' \item{statistic_scale}{The statistic finally used for the test following N(0,1)}
#' \item{pvalue}{p-value of the test statistic}
testsqr = function(y,X,Z,tau=0.5,kernel="Gaussian",h=max(0.05,0.5*(log(ncol(Z))/nrow(X))^{0.3})){
  n=length(y)
  betazh = conquer.cv.reg(Z,y,kernel = kernel,penalty = "lasso",h=h,tau = tau)$coeff.min
  err = as.vector(y-betazh[1]-Z%*%betazh[-1])
  if(kernel=="Gaussian"){
    phi = pnorm(-err/h)-tau}else if(kernel=="uniform"){
      phi = (-err/h+1)/2*(abs(err/h)<=1)+(abs(err/h)>1)
    }else if(kernel=="logistic"){
      phi = exp(-err/h)/(1+exp(-err/h))^2
    }
  #statistic
  Xw = X*phi
  a = apply(Xw, 2, sum)
  T1 = sum((a)^2)-sum(Xw^2)
  T1 = T1/n/(n-1)
  #variance estimation
  G1 = Xw%*%t(Xw)
  trs1 = 2*(sum(G1^2)-sum(diag(G1)^2))/n/(n-1)
  sig1 = sqrt(trs1)/n
  z = T1/sig1
  pvalue = 2*(1-pnorm(abs(z)))
  return(list(statistic = T1, standard.deviation = sig1, statistic_scale = z, pvalue = pvalue, tau=tau, kernel=kernel))
}





#' The projected version of testsqr
#' @param X data matrix of dimension n*px.
#' @param Z data matrix of dimension n*pz.
#' @param tau the quantile level with value between (0,1). The default is 0.5.
#' @param kernel A character string specifying the choice of kernel function., the default is "Gaussian".
#' Choices are "Gaussian", "uniform" and "logistic"
#' @param h The bandwidth parameter for kernel smoothing. The default is \eqn{\max\{0.05,0.5*(\log p_z/n)^{0.3}\}}
#' @param K the number of folds. Default is 5.
#' @param lambda user supplied tuning parameter; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. 
#' Default is sqrt(log(max(pz,px))/n).
#' @param err the precision used to stop the convergence. Default is 1e-5. 
#' Iterations stop when average absolute parameter change is less than \code{err}.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A list with components
#' \item{tau}{The chosen quantile level}
#' \item{kernel}{The chosen kernel}
#' \item{statistic}{The value of the proposed U-statistic}
#' \item{standard.deviation}{The value of the standard deviation of the test statistic}
#' \item{statistic_scale}{The statistic finally used for the test following N(0,1)}
#' \item{pvalue}{p-value of the test statistic}
testsqr.p = function(y,X,Z,tau=0.5,kernel="Gaussian",h=max(0.05,0.5*sqrt(log(ncol(Z)))/nrow(X)^{0.3}),
                     K=5,lambda = NULL, lambda.min=sqrt(log(max(ncol(X),ncol(Z)))/nrow(X)), nlambda=50, rho=1, err=1e-5, maxIter=1e3){
  betazh = conquer.cv.reg(Z,y,kernel = kernel,penalty = "lasso",h=h,tau = tau)$coeff.min
  err = as.vector(y-betazh[1]-Z%*%betazh[-1])
  if(kernel=="Gaussian"){
    phi = pnorm(-err/h)-tau}else if(kernel=="uniform"){
      phi = (-err/h+1)/2*(abs(err/h)<=1)+(abs(err/h)>1)
    }else if(kernel=="logistic"){
      phi = exp(-err/h)/(1+exp(-err/h))^2
    }
  #projection matrix 
  w = cvhiqr(X,Z,K=K,lambda=lambda,lambda.min = lambda.min,nlambda = nlambda,rho = rho,err = err,maxIter = maxIter)
  
  #projection test statistic
  Xwp = (X-Z%*%w)*phi
  ap = apply(Xwp, 2, sum)
  Tp = sum((ap)^2)-sum(Xwp^2)
  Tp = Tp/n/(n-1)
  #variance estimation
  Gp = Xwp%*%t(Xwp)
  trsp = 2*(sum(Gp^2)-sum(diag(Gp)^2))/n/(n-1)
  sigp = sqrt(trsp)/n
  zp = Tp/sigp
  pvalue = 2*(1-pnorm(abs(zp)))
  return(list(statistic = Tp, standard.deviation = sigp, statistic_scale = zp, pvalue = pvalue, tau=tau, kernel=kernel))
}
