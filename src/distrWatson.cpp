#include "animaWatsonDistribution.h"

#include <Eigen/Core>

#include <boost/math/tools/roots.hpp>

//' The Watson distribution
//'
//' @param x A numeric matrix of shape \eqn{n \times 3} where \eqn{n} is the
//'   number of samples and the columns represent the x, y, z coordinates
//'   of the samples.
//' @param mu A numeric matrix of shape \eqn{1 \times 3} representing the mean
//'   axis of the Watson distribution.
//' @param kappa A positive numeric value representing the concentration
//'   parameter of the Watson distribution.
//' @param log A logical value indicating whether the log density should be
//'   returned. Defaults to `false`.
//' @param n An integer value indicating the number of samples to generate.
//'
//' @return
//' - `dwatson` returns a numeric vector of length \eqn{n} containing the
//' density values of the Watson distribution.
//' - `pwatson` returns a numeric vector of length \eqn{n} containing the
//' cumulative density values of the Watson distribution.
//' - `qwatson` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the quantile values of the Watson distribution.
//' - `rwatson` returns a numeric matrix of shape \eqn{n \times 3} containing
//' the random samples from the Watson distribution.
//'
//' @examples
//' mu <- c(1, 0, 0)
//' kappa <- 10
//' n <- 100
//' spl <- rwatson(n, mu, kappa)
//' dwatson(spl, mu, kappa)
//' pwatson(spl, mu, kappa)
//'
//' @name watson

//' @export
//' @rdname watson
// [[Rcpp::export]]
Eigen::VectorXd dwatson(const Eigen::MatrixXd &x,
                        const Eigen::RowVectorXd &mu,
                        double kappa,
                        bool log = false) {
 using distr = anima::WatsonDistribution;
 distr watsonDistr;
 watsonDistr.SetMeanAxis(Eigen::RowVector3d(mu));
 watsonDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);
 for (unsigned int i = 0;i < n;++i)
 {
   if (log)
     res(i) = watsonDistr.GetLogDensity(xr.row(i));
   else
     res(i) = watsonDistr.GetDensity(xr.row(i));
 }
 return res;
}

//' @export
//' @rdname watson
// [[Rcpp::export]]
Eigen::VectorXd pwatson(const Eigen::MatrixXd &x,
                        const Eigen::RowVectorXd &mu,
                        double kappa) {
 using distr = anima::WatsonDistribution;
 distr watsonDistr;
 watsonDistr.SetMeanAxis(Eigen::RowVector3d(mu));
 watsonDistr.SetConcentrationParameter(kappa);
 Eigen::Ref<const Eigen::MatrixX3d> xr(x);
 unsigned int n = xr.rows();
 Eigen::VectorXd res(n);
 for (unsigned int i = 0;i < n;++i)
   res(i) = watsonDistr.GetCumulative(xr.row(i));
 return res;
}

//' @export
//' @rdname watson
// [[Rcpp::export]]
Eigen::MatrixXd rwatson(unsigned int n,
                        const Eigen::RowVectorXd &mu,
                        double kappa) {
 using distr = anima::WatsonDistribution;
 distr watsonDistr;
 distr::SampleType samples(n, 3);
 watsonDistr.SetMeanAxis(Eigen::RowVector3d(mu));
 watsonDistr.SetConcentrationParameter(kappa);
 distr::GeneratorType generator(std::time(0));
 watsonDistr.Random(samples, generator);
 return samples;
}

// [[Rcpp::export]]
Eigen::RowVectorXd mean_watson_impl(const Eigen::MatrixXd &x) {
  using distr = anima::WatsonDistribution;
  distr watsonDistr;
  Eigen::Ref<const Eigen::MatrixX3d> xr(x);
  watsonDistr.Fit(xr, "");
  return watsonDistr.GetMeanAxis();
}
