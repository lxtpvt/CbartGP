% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cbart.R
\name{cbart}
\alias{cbart}
\title{Correlated BART for continuous outcomes}
\usage{
cbart(
  x.train,
  y.train,
  SigmaInv,
  x.test = matrix(0, 0, 0),
  theta = 0,
  omega = 1,
  a = 0.5,
  b = 1,
  augment = FALSE,
  rho = NULL,
  xinfo = matrix(0, 0, 0),
  usequants = FALSE,
  cont = FALSE,
  rm.const = TRUE,
  k = 2,
  power = 2,
  base = 0.95,
  sigmaf = NA,
  lambda = NA,
  ntree = 200L,
  numcut = 100L,
  ndpost = 1000L,
  nskip = 100L,
  keepevery = 1L,
  nkeeptrain = ndpost,
  nkeeptest = ndpost,
  nkeeptestmean = ndpost,
  nkeeptreedraws = ndpost,
  printevery = 100L,
  transposed = FALSE,
  message = TRUE
)
}
\arguments{
\item{x.train}{Explanatory variables for training (in sample) data. May be a matrix or a data frame, with (as usual) rows corresponding to observations and columns to variables. If a variable is a factor in a data frame, it is replaced with dummies. Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor. %\code{makeind} is used to generate the dummies. \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.train}.}

\item{y.train}{Continuous dependent variable for training (in sample) data.}

\item{x.test}{Explanatory variables for test (out of sample) data. Should have same structure as x.train. \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.}

\item{theta}{Set \eqn{theta} parameter; zero means random.}

\item{omega}{Set \eqn{omega} parameter; zero means random.}

\item{a}{Sparse parameter for \eqn{Beta(a, b)} prior: \eqn{0.5<=a<=1} where lower values inducing more sparsity.}

\item{b}{Sparse parameter for \eqn{Beta(a, b)} prior; typically, \eqn{b=1}.}

\item{augment}{Whether data augmentation is to be performed in sparse variable selection.}

\item{rho}{Sparse parameter: typically \eqn{rho=p} where \eqn{p} is the number of covariates under consideration.}

\item{usequants}{If \code{usequants=FALSE}, then the cutpoints in \code{xinfo} are generated uniformly; otherwise, if \code{TRUE}, uniform quantiles are used for the cutpoints.}

\item{cont}{Whether or not to assume all variables are continuous.}

\item{rm.const}{Whether or not to remove constant variables.}

\item{k}{For numeric y, k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5. The response (y.train) is internally scaled to range from -.5 to .5. k is the number of prior standard deviations \eqn{f(x)} is away from +/-3. The bigger k is, the more conservative the fitting will be.}

\item{power}{Power parameter for tree prior.}

\item{base}{Base parameter for tree prior.}

\item{sigmaf}{The SD of f.}

\item{lambda}{The scale of the prior for the variance.}

\item{ntree}{The number of trees in the sum.}

\item{numcut}{The number of possible values of c (see usequants). If a single number if given, this is used for all variables. Otherwise a vector with length equal to ncol(x.train) is required, where the \eqn{i^{th}}{i^th} element gives the number of c used for the \eqn{i^{th}}{i^th} variable in x.train. If usequants is false, numcut equally spaced cutoffs are used covering the range of values in the corresponding column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the corresponding columns of x.train - 1) c values are used.}

\item{ndpost}{The number of posterior draws returned.}

\item{nskip}{Number of MCMC iterations to be treated as burn in.}

\item{keepevery}{Every keepevery draw is kept to be returned to the user.}

\item{nkeeptrain}{Number of MCMC iterations to be returned for train data.}

\item{nkeeptest}{Number of MCMC iterations to be returned for test data.}

\item{nkeeptestmean}{Number of MCMC iterations to be returned for test mean.}

\item{nkeeptreedraws}{Number of MCMC iterations to be returned for tree draws.}

\item{printevery}{As the MCMC runs, a message is printed every printevery draws.}

\item{transposed}{When running \code{wbart} in parallel, it is more memory-efficient to transpose \code{x.train} and \code{x.test}, if any, prior to calling \code{mc.wbart}.}
}
\value{
A list that stores the results of CBART.
}
\description{
Correlated BART for continuous outcomes
}
\references{
The paper
\emph{Gaussian processes Correlated Bayesian Additive Regression Trees}.
}
