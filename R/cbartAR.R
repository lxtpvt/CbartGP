#' CBART-GP for time series data.
#'
#'
#'
#' @param x.train    Explanatory variables for training (in sample) data. May be a matrix or a data frame, with (as usual) rows corresponding to observations and columns to variables. If a variable is a factor in a data frame, it is replaced with dummies. Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor. %\code{makeind} is used to generate the dummies. \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.train}.
#' @param y.train Continuous dependent variable for training (in sample) data.
#' @param x.test Explanatory variables for test (out of sample) data. Should have same structure as x.train. \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
#' @param weights Weights for searching.
#' @param seed seed.
#' 
#' @details 
#'
#' @return A list that stores the results of CBART-GP for 1D data.
#'
#' @references The paper
#' \emph{Gaussian processes Correlated Bayesian Additive Regression Trees}.
#'
#' @examples
#'
#'
#' @export
cbartAR1 <- function(x.train, y.train, x.test=matrix(0.0,0,0), weights=seq(0, 1, length.out = 6), seed=1)
{
  resBART = wbart(x.train=x.train, y.train=y.train, x.test=x.test, message=F)
  residuals.BART = y.train-resBART$yhat.train.mean
  residuals.demean = y.train-mean(y.train)
  
  arima.SS = c()
  new.residual.SS = c()
  cbart.SS = c()
  cbart.list = list()
  AR1.list = list()
  mle.list = list()
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(weights), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for (i in 1:length(weights)) {
    new.residual = weights[i]*residuals.demean+(1-weights[i])*residuals.BART
    ar.mod = arima(x=new.residual, order = c(1, 0L, 0L), include.mean=F)
    AR1.list = append(AR1.list, list(ar.mod))
    new.rou = as.vector(ar.mod$coef)[1]
    new.sigma2 = ar.mod$sigma2
    mle.list = append(mle.list, list(c(new.rou, new.sigma2)))
    
    arimaSigamInv(length(y.train),new.rou,sqrt(new.sigma2))->SigmaInv.new
    set.seed(seed) #it is MCMC, set the seed!!
    cbart_arc.new = cbart(x.train = x.train, y.train = y.train, SigmaInv=SigmaInv.new, x.test=x.test, message=F)
    cbart.list = append(cbart.list, list(cbart_arc.new))
    new.residual.SS=c(new.residual.SS,sum(new.residual^2))
    cbart.SS=c(cbart.SS,sum((y.train-cbart_arc.new$yhat.train.mean)^2))
    arima.SS = c(arima.SS,sum(ar.mod$residuals^2))
    setTxtProgressBar(pb, i)
  }
  close(pb) # Close the connection
  res = list(cbarts = cbart.list, AR1=AR1.list, SS.delta.vec = abs(cbart.SS-(new.residual.SS-arima.SS)), mles=mle.list)
  return(res)
}

#' Prediction AR1
#' 
#' @param residuals.tr The residuals (y-f_CBART) of training data. 
#' @param predCbart.te The CBART prediction on testing data.
#' @param xtr 1D positions of training data.
#' @param xte 1D positions of testing data.
#' @param rou The maximum likelihood estimations of AR(1)'s parameter rou.
#' 
#' @details 
#'
#' @return A vector of the predicted values of the testing data.
#'
#' @references The paper
#' \emph{Correlated Bayesian Additive Regression Trees with Gaussian Processes for Regression Analysis of Dependent Data}.
#' 
#' @export
predCbartAR1 <- function(residuals.tr, predCbart.te, xtr, xte, rou){
  
  yte = c()
  for (i in 1:length(xte)) {
    a = xte[i]-xtr
    if(length(which(a==0))>0){
      yte = c(yte, residuals.tr[which(a==0)])
    }
    if(length(which(a>0))>0){
      max(which(a>0))->leftId.tr
      yte = c(yte,predCbart.te[i]+rou*residuals.tr[leftId.tr])
    }
    if(length(which(a>=0))==0){
      min(which(a<0))->rightId.tr
      yte = c(yte,predCbart.te[i]+(1/rou)*residuals.tr[rightId.tr])
    }
  }
  return(yte)
}
