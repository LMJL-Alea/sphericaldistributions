#include "animaUniformDistribution.h"

#include <Eigen/Core>

// directive for openMP
#ifdef _OPENMP
#include <omp.h>
#endif

//' The Uniform distribution
//'
//' @param x A numeric matrix of shape \eqn{n \times 3} where \eqn{n} is the
//'   number of samples and the columns represent the x, y, z coordinates
//'   of the samples.
//' @param mu A numeric matrix of shape \eqn{1 \times 3} representing the mean
//'   axis of the Uniform distribution.
//' @param kappa A positive numeric value representing the concentration
//'   parameter of the Uniform distribution.
//' @param log A logical value indicating whether the log density should be
//'   returned. Defaults to `false`.
//' @param n An integer value indicating the number of samples to generate.
//'
//' @return
//' - `duni` returns a numeric vector of length \eqn{n} containing the
//' density values of the Uniform distribution.
//' - `puni` returns a numeric vector of length \eqn{n} containing the
//' cumulative density values of the Uniform distribution.
//' - `quni` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the quantile values of the Uniform distribution.
//' - `runi` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the random samples from the Uniform distribution.
//'
//' @examples
//' mu <- c(1, 0, 0)
//' kappa <- 10
//' n <- 100
//' spl <- runi(n, mu, kappa)
//' duni(spl, mu, kappa)
//' puni(spl, mu, kappa)
//'
//' @name Uniform

//' @export
//' @rdname Uniform
// [[Rcpp::export]]
Eigen::VectorXd duni(const Eigen::MatrixXd &x,
                    const Eigen::RowVectorXd &mu,
                    double kappa,
                    bool log = false) {
 using distr = anima::UniformDistribution;
 distr uniDistr;
 uniDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 uniDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

 for (unsigned int i = 0;i < n;++i)
 {
   if (log)
     res(i) = uniDistr.GetLogDensity(xr.row(i));
   else
     res(i) = uniDistr.GetDensity(xr.row(i));
 }
 return res;
}

//' @export
//' @rdname Uniform
// [[Rcpp::export]]
Eigen::VectorXd puni(const Eigen::MatrixXd &x,
                    const Eigen::RowVectorXd &mu,
                    double kappa) {
 using distr = anima::UniformDistribution;
 distr uniDistr;
 uniDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 uniDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

 for (unsigned int i = 0;i < n;++i)
   res(i) = uniDistr.GetCumulative(xr.row(i));
 return res;
}

//' @export
//' @rdname Uniform
// [[Rcpp::export]]
Eigen::MatrixXd runi(unsigned int n,
                    const Eigen::RowVectorXd &mu,
                    double kappa) {
 using distr = anima::UniformDistribution;
 distr uniDistr;
 distr::SampleType samples(n, 3);
 uniDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 uniDistr.SetConcentrationParameter(kappa);
 distr::GeneratorType generator(std::time(0));
 uniDistr.Random(samples, generator);
 return samples;
}

// [[Rcpp::export]]
Eigen::RowVectorXd mean_uni_impl(const Eigen::MatrixXd &x) {
 using distr = anima::UniformDistribution;
 distr uniDistr;
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 uniDistr.Fit(xr, "");
 return uniDistr.GetMeanDirection();
}
