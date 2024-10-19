#' Efficient admm algorithm for projection matrix estimation
#' @param X data matrix of dimension n*px.
#' @param Z data matrix of dimension n*pz.
#' @param lambda user supplied tuning parameter; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. 
#' Default is sqrt(log(p)/n).
#' @param err the precision used to stop the convergence. Default is 1e-5. 
#' Iterations stop when average absolute parameter change is less than \code{err}.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A list with components
#' \item{Omega}{a list of sparse pz*px matrices corresponding to lambda.}
#' \item{lambda}{the used lambda for the solution path.}
#' \item{niter}{the number of iterations for each element of lambda.}
hiqr <- function(X,Z,lambda = NULL, lambda.min=sqrt(log(max(ncol(X),ncol(Z)))/nrow(X)),nlambda=50,err=10^(-5),maxIter =1000,rho=1){
  px = ncol(X)
  pz = ncol(Z)
  n = nrow(X)
  
  Sn = cov(Z);
  A<-abs(Sn-diag(diag(Sn)));
  if (is.null(lambda)){lambda=exp(seq(log(lambda.min),0,length.out =nlambda))*max(A)}
  lambda<-sort(lambda,decreasing =TRUE)
  
  covzx<- cov(Z,X)
  
  obj <- projlasso(Z,covzx,lambda,err,maxIter,rho)
  return(obj)
}


#' Cross-validation function for hiqr
#' @param X data matrix of dimension n*px.
#' @param Z data matrix of dimension n*pz.
#' @param K the number of folds. Default is 5.
#' @param lambda user supplied tuning parameter; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. 
#' Default is sqrt(log(p)/n).
#' @param err the precision used to stop the convergence. Default is 1e-5. 
#' Iterations stop when average absolute parameter change is less than \code{err}.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A list with components
#' \item{Ome}{the estimated pz*px precision matrix.}
#' \item{cvlambda}{the chosen lambda by cross-validation.}
#' \item{lambda}{the used lambda list for cross-validation.}
#' \item{cvloss}{the empirical loss of cross-validation related to lambda.}
cvhiqr = function(X,Z,K=5,lambda = NULL, lambda.min=sqrt(log(max(ncol(X),ncol(Z)))/nrow(X)), nlambda=50, rho=1, err=1e-5, maxIter=1e3){
  n = nrow(X)
  Sn = cov(Z);
  A = abs(Sn-diag(diag(Sn)))
  if(is.null(lambda)){lambda=exp(seq(log(lambda.min),0,length.out = nlambda))*max(A)}
  lambda<-sort(lambda,decreasing =T)
  
  fold<-split(1:n, rep(1:K, length =n))
  CV<-matrix(0,nrow=K,ncol=nlambda)
  for (k in 1:K){
    xtrain=X[-fold[[k]],]
    xtest=X[fold[[k]],]
    ztrain=Z[-fold[[k]],]
    ztest=Z[fold[[k]],]
    zs<-scale(ztest,center =TRUE,scale =FALSE)/sqrt(nrow(ztest)-1);
    obj<-hiqr(xtrain,ztrain,lambda =lambda,err=10*err, maxIter = 100)
    for (i in 1:nlambda){
      hOme<-as.matrix(obj$Ome_all[[i]])
      CV[k,i]=sum((zs%*%hOme)^2)/2-sum(hOme*cov(ztest,xtest))
    }
  }
  mcv<-apply(CV, 2, mean)
  lambda<-obj$lambda
  cvlambda<-lambda[which.min(mcv)]
  hOme = hiqr(X,Z,lambda =cvlambda,err=err, maxIter = maxIter,rho=rho)
  return(list(Ome = as.matrix(hOme$Ome_all[[1]]), cvlambda = cvlambda, lambda = lambda, CVloss = mcv))
}






