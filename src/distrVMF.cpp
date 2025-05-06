#include "animaVonMisesFisherDistribution.h"

#include <Eigen/Core>

// directive for openMP
#ifdef _OPENMP
#include <omp.h>
#endif

//' The Von-Mises Fisher distribution
//'
//' @param x A numeric matrix of shape \eqn{n \times 3} where \eqn{n} is the
//'   number of samples and the columns represent the x, y, z coordinates
//'   of the samples.
//' @param mu A numeric matrix of shape \eqn{1 \times 3} representing the mean
//'   axis of the Von Mises-Fisher distribution.
//' @param kappa A positive numeric value representing the concentration
//'   parameter of the Von Mises-Fisher distribution.
//' @param log A logical value indicating whether the log density should be
//'   returned. Defaults to `false`.
//' @param n An integer value indicating the number of samples to generate.
//'
//' @return
//' - `dvmf` returns a numeric vector of length \eqn{n} containing the
//' density values of the Von Mises-Fisher distribution.
//' - `pvmf` returns a numeric vector of length \eqn{n} containing the
//' cumulative density values of the Von Mises-Fisher distribution.
//' - `qvmf` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the quantile values of the Von Mises-Fisher distribution.
//' - `rvmf` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the random samples from the Von Mises-Fisher distribution.
//'
//' @examples
//' mu <- c(1, 0, 0)
//' kappa <- 10
//' n <- 100
//' spl <- rvmf(n, mu, kappa)
//' dvmf(spl, mu, kappa)
//' pvmf(spl, mu, kappa)
//'
//' @name VonMisesFisher

//' @export
//' @rdname VonMisesFisher
// [[Rcpp::export]]
Eigen::VectorXd dvmf(const Eigen::MatrixXd &x,
                    const Eigen::RowVectorXd &mu,
                    double kappa,
                    bool log = false) {
 using distr = anima::VonMisesFisherDistribution;
 distr vmfDistr;
 vmfDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 vmfDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

 for (unsigned int i = 0;i < n;++i)
 {
   if (log)
     res(i) = vmfDistr.GetLogDensity(xr.row(i));
   else
     res(i) = vmfDistr.GetDensity(xr.row(i));
 }
 return res;
}

//' @export
//' @rdname VonMisesFisher
// [[Rcpp::export]]
Eigen::VectorXd pvmf(const Eigen::MatrixXd &x,
                    const Eigen::RowVectorXd &mu,
                    double kappa) {
 using distr = anima::VonMisesFisherDistribution;
 distr vmfDistr;
 vmfDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 vmfDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

 for (unsigned int i = 0;i < n;++i)
   res(i) = vmfDistr.GetCumulative(xr.row(i));
 return res;
}

//' @export
//' @rdname VonMisesFisher
// [[Rcpp::export]]
Eigen::MatrixXd rvmf(unsigned int n,
                    const Eigen::RowVectorXd &mu,
                    double kappa) {
 using distr = anima::VonMisesFisherDistribution;
 distr vmfDistr;
 distr::SampleType samples(n, 3);
 vmfDistr.SetMeanDirection(Eigen::RowVector3d(mu));
 vmfDistr.SetConcentrationParameter(kappa);
 distr::GeneratorType generator(std::time(0));
 vmfDistr.Random(samples, generator);
 return samples;
}

// [[Rcpp::export]]
Eigen::RowVectorXd mean_vmf_impl(const Eigen::MatrixXd &x) {
 using distr = anima::VonMisesFisherDistribution;
 distr vmfDistr;
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 vmfDistr.Fit(xr, "");
 return vmfDistr.GetMeanDirection();
}
