#' CBART-GP for spatial data.
#'
#'
#' @param coords 2D spatial coordinates.
#' @param x.train Explanatory variables for training (in sample) data. May be a matrix or a data frame, with (as usual) rows corresponding to observations and columns to variables. If a variable is a factor in a data frame, it is replaced with dummies. Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' @param y.train Continuous dependent variable for training (in sample) data.
#' @param x.test Explanatory variables for test (out of sample) data. Should have same structure as x.train. \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
#' @param method Method for estimating MLE of Matern parameters. The methods "geoR" and "fields" correspond to the related R package names.
#' @param weights Weights for searching.
#' @param seed seed.
#' 
#' @details 
#'
#' @return A list that stores the results of CBART-GP for 2D data.
#'
#' @references The paper
#' \emph{Correlated Bayesian Additive Regression Trees with Gaussian Processes for Regression Analysis of Dependent Data}.
#'
#'
#'
#' @export
cbartMatern <- function(coords, x.train, y.train, x.test=matrix(0.0,0,0), method="geoR", weights=seq(0, 1, length.out = 6), seed=1)
{
  resBART = wbart(x.train=x.train, y.train=y.train, x.test=x.test, message=F)
  residuals.BART = y.train-resBART$yhat.train.mean
  residuals.demean = y.train-mean(y.train)
  n.y.train = length(y.train)
  
  matern.SS = c()
  new.residual.SS = c()
  cbart.2D.SS = c()

  cbart.list = list()
  mle.list = list()
  md.dep.list = list()
  
  D <- as.matrix(dist(coords))

  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(weights), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for (i in 1:length(weights)) {
    new.residual = weights[i]*residuals.demean+(1-weights[i])*residuals.BART
    # The MLE of matern covariance function
    geodata.2D = list(coords=coords, data=new.residual)
    if(method=="geoR"){
      maten.mle = suppressWarnings(suppressMessages(geoR::likfit(geodata=geodata.2D, coords=geodata.2D$coords, data=geodata.2D$data,
                                                           ini.cov.pars = c(0.5, 0.5), kappa = .5, messages=F)))
      matern.residuals = resid(maten.mle)
      md.dep = maten.mle
      ntree = 100
      ndpost=1000L
      nskip=100L
    }else if(method=="fields"){
      fit.uk2 <- suppressWarnings(suppressMessages(fields::Krig(x=coords, Y=new.residual, verbose=F)))
      maten.mle = list(cov.pars=c(fit.uk2$sigma.MLE^2, 1), nugget=fit.uk2$tauHat.MLE^2)
      matern.residuals = fit.uk2$residuals
      md.dep = fit.uk2
      ntree = 50
      ndpost=50
      nskip=10
    }else{
      return('Choose method of MLE estimation: "geoR" or "fields".')
    }

    mle.list = append(mle.list, list(c(maten.mle$cov.pars,maten.mle$nugget)))
    # The inverse covariance matrix
    C <- maten.mle$cov.pars[1]*exp(-D/maten.mle$cov.pars[2])
    sigmadiag = maten.mle$nugget*diag(n.y.train)
    Sigma = sigmadiag+C
    SigmaInv <- chol2inv(chol(Sigma))
    # run CBART
    set.seed(seed)
    cbart.2D = cbart(x.train = x.train, y.train = y.train, SigmaInv=SigmaInv, x.test=x.test,
                     ntree = ntree, ndpost=ndpost, nskip=nskip, message=F)
    
    new.residual.SS=c(new.residual.SS,sum(geodata.2D$data^2))
    cbart.2D.SS=c(cbart.2D.SS,sum((y.train-cbart.2D$yhat.train.mean)^2))
    matern.SS = c(matern.SS,sum(matern.residuals^2))
    
    cbart.list = append(cbart.list,list(cbart.2D))
    md.dep.list = append(md.dep.list,list(md.dep))
    setTxtProgressBar(pb, i)
  }
  close(pb) # Close the connection
  res = list(cbarts = cbart.list, SS.delta.vec = abs(cbart.2D.SS-(new.residual.SS-matern.SS)), mles=mle.list, md.dep=md.dep.list)
  return(res)

}

#' Prediction Matern
#' 
#' @param residuals.tr The residuals (y-f_CBART) of training data. 
#' @param predCbart.te The CBART prediction on testing data.
#' @param coords.tr 2D spatial coordinates of training data.
#' @param coords.te 2D spatial coordinates of testing data.
#' @param matern.cov The maximum likelihood estimations of matern covariance function's parameters.
#' @param method NULL, "geoR", "fields"
#' 
#' @details 
#'
#' @return A vector of the predicted values of the testing data.
#'
#' @references The paper
#' \emph{Correlated Bayesian Additive Regression Trees with Gaussian Processes for Regression Analysis of Dependent Data}.
#' 
#' @export
predCbartMatern <- function(residuals.tr, predCbart.te, coords.tr, coords.te, method=NULL, matern.cov=NULL){
  
  yte = NULL
  if(is.null(method)){
    if(is.null(matern.cov)){
      return("The matern.cov should not be NULL !")
    }
    D.tr <- as.matrix(dist(coords.tr))
    K.tr <- matern.cov[1]*exp(-D.tr/matern.cov[2])
    K.tr.inv <- chol2inv(chol(K.tr))
    
    coord.all = rbind(coords.tr,coords.te)
    D.all = as.matrix(dist(coord.all))
    D.sub = D.all[(dim(coords.tr)[1]+1):dim(coord.all)[1],1:dim(coords.tr)[1]]
    K.sub = matern.cov[1]*exp(-D.sub/matern.cov[2])
    yte = predCbart.te + K.sub%*%K.tr.inv%*%residuals.tr
    
  }else if(method$name == "geoR"){
    
    kc.data = list(coords=coords.tr, data=residuals.tr)
    class(kc.data)<-"geodata"
    kc <- geoR::krige.conv(kc.data, locations = coords.te, krige = geoR::krige.control(obj.model = method$model))
    yte = predCbart.te + kc$predict
    
  }else if (method$name == "fields"){
    
    fields.pred = fields::predict.Krig(object=method$model, x=coords.te, drop.Z=TRUE)
    yte = predCbart.te + as.vector(fields.pred)
    
  }

  return(yte)
}
