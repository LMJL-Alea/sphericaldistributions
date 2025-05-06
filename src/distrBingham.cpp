#include "animaBinghamDistribution.h"

#include <Eigen/Core>

// directive for openMP
#ifdef _OPENMP
#include <omp.h>
#endif

//' The Bingham distribution
//'
//' @param x A numeric matrix of shape \eqn{n \times 3} where \eqn{n} is the
//'   number of samples and the columns represent the x, y, z coordinates
//'   of the samples.
//' @param mu A numeric matrix of shape \eqn{1 \times 3} representing the mean
//'   axis of the Bingham distribution.
//' @param kappa A positive numeric value representing the concentration
//'   parameter of the Bingham distribution.
//' @param log A logical value indicating whether the log density should be
//'   returned. Defaults to `false`.
//' @param n An integer value indicating the number of samples to generate.
//'
//' @return
//' - `dbin` returns a numeric vector of length \eqn{n} containing the
//' density values of the Bingham distribution.
//' - `pbin` returns a numeric vector of length \eqn{n} containing the
//' cumulative density values of the Bingham distribution.
//' - `qbin` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the quantile values of the Bingham distribution.
//' - `rbin` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the random samples from the Bingham distribution.
//'
//' @examples
//' mu <- c(1, 0, 0)
//' kappa <- 10
//' n <- 100
//' spl <- rbin(n, mu, kappa)
//' dbin(spl, mu, kappa)
//' pbin(spl, mu, kappa)
//'
//' @name Bingham

//' @export
//' @rdname Bingham
// [[Rcpp::export]]
Eigen::VectorXd dbin(const Eigen::MatrixXd &x,
                    const Eigen::RowVectorXd &mu,
                    double kappa,
                    bool log = false) {
 using distr = anima::BinghamDistribution;
 distr binDistr;
 binDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 binDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

 for (unsigned int i = 0;i < n;++i)
 {
   if (log)
     res(i) = binDistr.GetLogDensity(xr.row(i));
   else
     res(i) = binDistr.GetDensity(xr.row(i));
 }
 return res;
}

//' @export
//' @rdname Bingham
// [[Rcpp::export]]
Eigen::VectorXd pbin(const Eigen::MatrixXd &x,
                    const Eigen::RowVectorXd &mu,
                    double kappa) {
 using distr = anima::BinghamDistribution;
 distr binDistr;
 binDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 binDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

 for (unsigned int i = 0;i < n;++i)
   res(i) = binDistr.GetCumulative(xr.row(i));
 return res;
}

//' @export
//' @rdname Bingham
// [[Rcpp::export]]
Eigen::MatrixXd rbin(unsigned int n,
                    const Eigen::RowVectorXd &mu,
                    double kappa) {
 using distr = anima::BinghamDistribution;
 distr binDistr;
 distr::SampleType samples(n, 3);
 binDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 binDistr.SetConcentrationParameter(kappa);
 distr::GeneratorType generator(std::time(0));
 binDistr.Random(samples, generator);
 return samples;
}

// [[Rcpp::export]]
Eigen::RowVectorXd mean_bin_impl(const Eigen::MatrixXd &x) {
 using distr = anima::BinghamDistribution;
 distr binDistr;
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 binDistr.Fit(xr, "");
 return binDistr.GetMeanDirection();
}
