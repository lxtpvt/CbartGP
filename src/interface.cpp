
#include "cdbartSigmaFix.h"
#include "wbart.h"
using namespace Rcpp;


// [[Rcpp::export]]
List cdbartSigmaFix_inf(
    size_t n,            //number of observations in training data
    size_t p,		//dimension of x
    size_t np,		//number of observations in test data
    NumericVector ixv,		//x, train,  pxn (transposed so rows are contiguous in memory)
    NumericVector iyv,		//y, train,  nx1
    NumericVector ixpv,		//x, test, pxnp (transposed so rows are contiguous in memory)
    Eigen::MatrixXd& SigmaInv,
    size_t m,		//number of trees
    IntegerVector numcutv,		//number of cut points
    size_t nd,		//number of kept draws (except for thinnning ..)
    size_t burn,		//number of burn-in draws skipped
    double mybeta,
    double alpha,
    double tau,
    bool dart,
    double theta,
    double omega,
    double a,
    double b,
    double rho,
    bool aug,
    size_t nkeeptrain,
    size_t nkeeptest,
    size_t nkeeptestme,
    size_t nkeeptreedraws,
    size_t printevery,
    bool showMessage
)
{

  double* ix = &ixv[0];
  double* iy = &iyv[0];
  double* ixp = &ixpv[0];
  int* numcut = &numcutv[0];
  double* trmean=new double[n];
  double* temean=new double[np];

  cdbartSigmaFix(n, p, np, ix, iy, ixp, m, numcut, nd, burn, mybeta, alpha, tau, dart, theta, omega, a, b,
    rho, aug, nkeeptrain, nkeeptest, nkeeptestme, nkeeptreedraws, printevery, trmean, temean, SigmaInv, showMessage);

  NumericVector trmeanv = NumericVector(trmean,trmean+n);
  NumericVector temeanv = NumericVector(temean,temean+np);

  //--------------------------------------------------
  //use wrap to return computed totdim to R as part of a list
  List ret; //list to return
  ret["yhat.train.mean"] = trmeanv;
  ret["yhat.test.mean"] = temeanv;
  return ret;
}



// [[Rcpp::export]]
List wbart_inf(
    size_t n,            //number of observations in training data
    size_t p,		//dimension of x
    size_t np,		//number of observations in test data
    NumericVector ixv,		//x, train,  pxn (transposed so rows are contiguous in memory)
    NumericVector iyv,		//y, train,  nx1
    NumericVector ixpv,		//x, test, pxnp (transposed so rows are contiguous in memory)
    size_t m,		//number of trees
    IntegerVector numcutv,		//number of cut points
    size_t nd,		//number of kept draws (except for thinnning ..)
    size_t burn,		//number of burn-in draws skipped
    double mybeta,
    double alpha,
    double tau,
    double nu,
    double lambda,
    double sigma,
    NumericVector iwv,
    bool dart,
    double theta,
    double omega,
    IntegerVector grpv,
    double a,
    double b,
    double rho,
    bool aug,
    size_t nkeeptrain,
    size_t nkeeptest,
    size_t nkeeptestme,
    size_t nkeeptreedraws,
    size_t printevery,
    Rcpp::NumericMatrix xifo,
    bool showMessage
)
{
  double* ix = &ixv[0];
  double* iy = &iyv[0];
  double* ixp = &ixpv[0];
  int* numcut = &numcutv[0];
  
  double* iw = &iwv[0];
  int* grp = &grpv[0];
  
  Rcpp::List ret;
  ret = cwbart(n, p, np, ix, iy, ixp, m, numcut, nd,
        burn, mybeta, alpha, tau, nu, lambda, sigma, iw,
        dart, theta, omega, grp, a, b, rho, aug, nkeeptrain,
        nkeeptest, nkeeptestme, nkeeptreedraws, printevery, xifo, showMessage);
  
  return ret;
}
