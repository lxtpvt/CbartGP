

#ifndef __EIGENMULTIVARIATENORMAL_HPP
#define __EIGENMULTIVARIATENORMAL_HPP

#include <Eigen/Dense>
#include <random>

#define PI 3.141592653589793238462643383280

namespace Eigen {
namespace internal {
template<typename Scalar>
struct scalar_normal_dist_op
{
  static std::mt19937 rng;                        // The uniform pseudo-random algorithm
  mutable std::normal_distribution<Scalar> norm; // gaussian combinator
  
  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)
    
    template<typename Index>
    inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
    inline void seed(const uint64_t &s) { rng.seed(s); }
};

template<typename Scalar>
std::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };

} // end namespace internal
} // end namespace Eigen

class Mvn
{
public:
	Mvn(const Eigen::VectorXd& mu, const Eigen::MatrixXd& s);
	~Mvn(){};
	void getMean(Eigen::VectorXd& mu){mu = mean;};
	void getCovar(Eigen::MatrixXd& covar){covar = Sigma;};
	void setMean(const Eigen::VectorXd& mu) { mean = mu; }
	void setCovar(const Eigen::MatrixXd& covar) { Sigma = covar; }
	void getTransform(Eigen::MatrixXd& mat, Eigen::MatrixXd& sqrtmat);
	double pdf(const Eigen::VectorXd& x);
	void sample(Eigen::MatrixXd& samples, int n);

private:
	Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
	Eigen::VectorXd mean;
	Eigen::MatrixXd Sigma;
};
#endif