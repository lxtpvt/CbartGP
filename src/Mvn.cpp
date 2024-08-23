

#include "Mvn.h"

Mvn::Mvn(const Eigen::VectorXd& mu, const Eigen::MatrixXd& s)
{
	this->mean = mu;
	this->Sigma = s;
}

double Mvn::pdf(const Eigen::VectorXd& x)
{
  double n = x.rows();
  double sqrt2pi = std::sqrt(2 * PI);
  double quadform  = (x - mean).transpose() * Sigma.inverse() * (x - mean);
  double norm = std::pow(sqrt2pi, - n) *
                std::pow(Sigma.determinant(), - 0.5);

  return norm * exp(-0.5 * quadform);
}

void Mvn::getTransform(Eigen::MatrixXd& mat, Eigen::MatrixXd& sqrtmat)
{
	Eigen::LLT<Eigen::MatrixXd> cholSolver(mat);

	// We can only use the cholesky decomposition if 
	// the covariance matrix is symmetric, pos-definite.
	// But a covariance matrix might be pos-semi-definite.
	// In that case, we'll go to an EigenSolver
	if (cholSolver.info()==Eigen::Success) {
		// Use cholesky solver
		sqrtmat = cholSolver.matrixL();
	} else {

		// Use eigen solver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(mat);
		sqrtmat = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
	}
}


void Mvn::sample(Eigen::MatrixXd& samples, int n)
{
	int sz = this->Sigma.rows();
	Eigen::MatrixXd normTransform(sz,sz);
	getTransform(this->Sigma,normTransform);
	samples = (normTransform * Eigen::MatrixXd::NullaryExpr(sz,n,randN)).colwise() + mean;
}