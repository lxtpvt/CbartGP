
#'
#' @export
arimaSigamInv<-function(n, rous, sigma){
  A = diag(n)
  n.rou = length(rous)
  for (i in 2:n) {
    for (j in 1:n.rou) {
      A[i,i-j] = -rous[j]
    }
  }
  SigmaInv = sigma^(-2)*t(A)%*%A
  return(SigmaInv)
}
#'
#' @export
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}
#'
#' @export
maternCov<-function(v,phi,sigma,coords){
  if(!(v %in% c(1/2,3/2,5/2))){
    return("v= 1/2, 3/2, or 5/2 !")
  }
  D <- as.matrix(dist(coords))
  if(v==1/2){
    R <- exp(-D/phi)
  }else if(v==3/2){
    R <- (1+sqrt(3)*D/phi)*exp(-sqrt(3)*D/phi)
  }else{
    R <- (1+sqrt(5)*D/phi+5*D^2/(3*phi^2))*exp(-sqrt(5)*D/phi)
  }
  Sigma = (sigma^2)*R
  SigmaInv = chol2inv(Sigma)
  return(list(cov=Sigma,covInv=SigmaInv))
}
#'
#' @export
proposedMSE<-function(n,Sigma,tau){
  w <- rmvn(1, rep(0,n), Sigma)
  e <- rnorm(n, 0, tau)
  a = w+e
  return(sum(a^2)/n)
}
#'
#' @export
weightedResiduals<-function(data, w){
  coords<-data$coords
  X<-data$X
  y<-data$y
  
  # the de-mean residuals
  residuals.y.demean = y-mean(y)
  # BART
  bartmd = wbart(X,y)
  residuals.bartmd = y-bartmd$yhat.train.mean
  

  # Remove the mean trend
  zeromeany=NULL
  for (wi in w) {
    temp = wi*residuals.y.demean + (1-wi)*residuals.bartmd
    zeromeany = cbind(zeromeany,temp)
  }
  return(zeromeany)
}

#'
#' @export


