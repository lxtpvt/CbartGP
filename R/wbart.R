#' BART for continuous outcomes
#'
#'
#'
#' @param x.train    Explanatory variables for training (in sample) data. May be a matrix or a data frame, with (as usual) rows corresponding to observations and columns to variables. If a variable is a factor in a data frame, it is replaced with dummies. Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor. %\code{makeind} is used to generate the dummies. \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.train}.
#' @param y.train Continuous dependent variable for training (in sample) data.
#' @param x.test Explanatory variables for test (out of sample) data. Should have same structure as x.train. \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
#' @param theta Set \eqn{theta} parameter; zero means random.
#' @param omega Set \eqn{omega} parameter; zero means random.
#' @param a Sparse parameter for \eqn{Beta(a, b)} prior: \eqn{0.5<=a<=1} where lower values inducing more sparsity.
#' @param b Sparse parameter for \eqn{Beta(a, b)} prior; typically, \eqn{b=1}.
#' @param rho Sparse parameter: typically \eqn{rho=p} where \eqn{p} is the number of covariates under consideration.
#' @param augment Whether data augmentation is to be performed in sparse variable selection.
#' @param xinfo ou can provide the cutpoints to BART or let BART choose them for you. To provide them, use the xinfo argument to specify a list (matrix) where the items (rows) are the covariates and the contents of the items (columns) are the cutpoints.
#' @param usequants If \code{usequants=FALSE}, then the cutpoints in \code{xinfo} are generated uniformly; otherwise, if \code{TRUE}, uniform quantiles are used for the cutpoints.
#' @param cont Whether or not to assume all variables are continuous.
#' @param rm.const Whether or not to remove constant variables.
#' @param sigest The prior for the error variance is inverted chi-squared (the standard conditionally conjugate prior). The prior is specified by choosing the degrees of freedom, a rough estimate of the corresponding standard deviation and a quantile to put this rough estimate at. If sigest=NA then the rough estimate will be the usual least squares estimator. Otherwise the supplied value will be used.
#' @param sigdf Degrees of freedom for error variance prior.
#' @param sigquant The quantile of the prior that the rough estimate (see sigest) is placed at. The closer the quantile is to 1, the more aggresive the fit will be as you are putting more prior weight on error standard deviations (\eqn{sigma}) less than the rough estimate.
#' @param k For numeric y, k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5. The response (y.train) is internally scaled to range from -.5 to .5. k is the number of prior standard deviations \eqn{f(x)} is away from +/-3. The bigger k is, the more conservative the fitting will be.
#' @param power Power parameter for tree prior.
#' @param base Base parameter for tree prior.
#' @param sigmaf The SD of f.
#' @param lambda The scale of the prior for the variance.
#' @param fmean BART operates on y.train centered by fmean.
#' @param w Vector of weights which multiply the standard deviation.
#' @param ntree The number of trees in the sum.
#' @param numcut The number of possible values of c (see usequants). If a single number if given, this is used for all variables. Otherwise a vector with length equal to ncol(x.train) is required, where the \eqn{i^{th}}{i^th} element gives the number of c used for the \eqn{i^{th}}{i^th} variable in x.train. If usequants is false, numcut equally spaced cutoffs are used covering the range of values in the corresponding column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the corresponding columns of x.train - 1) c values are used.
#' @param ndpost The number of posterior draws returned.
#' @param nskip Number of MCMC iterations to be treated as burn in.
#' @param nkeeptrain Number of MCMC iterations to be returned for train data.
#' @param nkeeptest Number of MCMC iterations to be returned for test data.
#' @param nkeeptestmean Number of MCMC iterations to be returned for test mean.
#' @param nkeeptreedraws Number of MCMC iterations to be returned for tree draws.
#' @param printevery As the MCMC runs, a message is printed every printevery draws.
#' @param keepevery Every keepevery draw is kept to be returned to the user.
#' @param transposed When running \code{wbart} in parallel, it is more memory-efficient to transpose \code{x.train} and \code{x.test}, if any, prior to calling \code{mc.wbart}.
#' 
#' 
#' @details 
#'
#' @return A list that stores the results of CBART.
#'
#' @references The paper
#' \emph{Gaussian processes Correlated Bayesian Additive Regression Trees}.
#'
#' @examples
#'
#'
#' @export
wbart=function(
x.train, y.train, x.test=matrix(0.0,0,0),
sparse=FALSE, theta=0, omega=1,
a=0.5, b=1, augment=FALSE, rho=NULL,
xinfo=matrix(0.0,0,0), usequants=FALSE,
cont=FALSE, rm.const=TRUE,
sigest=NA, sigdf=3, sigquant=.90,
k=2.0, power=2.0, base=.95,
sigmaf=NA, lambda=NA,
fmean=mean(y.train),
w=rep(1,length(y.train)),
ntree=200L, numcut=100L,
ndpost=1000L, nskip=100L, keepevery=1L,
nkeeptrain=ndpost, nkeeptest=ndpost,
nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
printevery=100L, transposed=FALSE, message=TRUE
)
{
#--------------------------------------------------
#data
n = length(y.train)

if(!transposed) {
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(x.test)>0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
    }
    rm.const <- temp$rm.const
    grp <- temp$grp
    rm(temp)
}
else {
    rm.const <- NULL
    grp <- NULL
}

if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')

p = nrow(x.train)
np = ncol(x.test)
if(length(rho)==0) rho=p
if(length(rm.const)==0) rm.const <- 1:p
if(length(grp)==0) grp <- 1:p

##if(p>1 & length(numcut)==1) numcut=rep(numcut, p)

y.train = y.train-fmean
#--------------------------------------------------
#set nkeeps for thinning
if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
   nkeeptrain=ndpost
   cat('*****nkeeptrain set to ndpost\n')
}
if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
   nkeeptest=ndpost
   cat('*****nkeeptest set to ndpost\n')
}
if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
   nkeeptestmean=ndpost
   cat('*****nkeeptestmean set to ndpost\n')
}
if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
   nkeeptreedraws=ndpost
   cat('*****nkeeptreedraws set to ndpost\n')
}
#--------------------------------------------------
#prior
nu=sigdf
if(is.na(lambda)) {
   if(is.na(sigest)) {
      if(p < n) {
         df = data.frame(t(x.train),y.train)
         lmf = lm(y.train~.,df)
         sigest = summary(lmf)$sigma
      } else {
         sigest = sd(y.train)
      }
   }
   qchi = qchisq(1.0-sigquant,nu)
   lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
} else {
   sigest=sqrt(lambda)
}

if(is.na(sigmaf)) {
   tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
} else {
   tau = sigmaf/sqrt(ntree)
}
#--------------------------------------------------
ptm <- proc.time()
#call
res = wbart_inf(
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #pxn training data x
            y.train,   #pxn training data x
            x.test,   #p*np test data x
            ntree,
            numcut,
            ndpost*keepevery,
            nskip,
            power,
            base,
            tau,
            nu,
            lambda,
            sigest,
            w,
            sparse,
            theta,
            omega,
            grp,
            a,
            b,
            rho,
            augment,
            nkeeptrain,
            nkeeptest,
            nkeeptestmean,
            nkeeptreedraws,
            printevery,
            xinfo,
            message
)
    
res$proc.time <- proc.time()-ptm
    
res$mu = fmean
res$yhat.train.mean = res$yhat.train.mean+fmean
res$yhat.train = res$yhat.train+fmean
res$yhat.test.mean = res$yhat.test.mean+fmean
res$yhat.test = res$yhat.test+fmean
if(nkeeptreedraws>0)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
##res$nkeeptreedraws=nkeeptreedraws
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
attr(res, 'class') <- 'wbart'
return(res)
}
