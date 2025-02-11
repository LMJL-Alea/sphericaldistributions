#include "animaWatsonDistribution.h"

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
arma::vec dwatson(const arma::mat &x, const arma::rowvec &mu, double kappa, bool log = false) {
 using distr = anima::WatsonDistribution;
 distr watsonDistr;
 watsonDistr.SetMeanAxis(mu);
 watsonDistr.SetConcentrationParameter(kappa);
 unsigned int n = x.n_rows;
 arma::vec res(n);
 for (unsigned int i = 0;i < n;++i)
 {
   if (log)
     res(i) = watsonDistr.GetLogDensity(x.row(i));
   else
     res(i) = watsonDistr.GetDensity(x.row(i));
 }
 return res;
}

//' @export
//' @rdname watson
// [[Rcpp::export]]
arma::vec pwatson(const arma::mat &x, const arma::rowvec &mu, double kappa) {
 using distr = anima::WatsonDistribution;
 distr watsonDistr;
 watsonDistr.SetMeanAxis(mu);
 watsonDistr.SetConcentrationParameter(kappa);
 unsigned int n = x.n_rows;
 arma::vec res(n);
 for (unsigned int i = 0;i < n;++i)
   res(i) = watsonDistr.GetCumulative(x.row(i));
 return res;
}

//' @export
//' @rdname watson
// [[Rcpp::export]]
arma::mat qwatson(const arma::vec &p, const arma::rowvec &mu, double kappa) {
 using distr = anima::WatsonDistribution;
 distr watsonDistr;
 watsonDistr.SetMeanAxis(mu);
 watsonDistr.SetConcentrationParameter(kappa);
 unsigned int n = p.n_elem;
 arma::mat res(n, 3);
 // for (unsigned int i = 0;i < n;++i)
 // {
 //   auto f = [&watsonDistr, &p, i](arma::vec x) { return watsonDistr.GetCumulative(x) - p[i]; };
 //   // use boost to find the root
 //   boost::uintmax_t max_iter = 100;
 //   boost::math::tools::eps_tolerance<double> tol(30);
 //   arma::vec x = arma::randu<arma::vec>(3);
 //   boost::uintmax_t iters = boost::math::tools::newton_raphson_iterate(f, x, x, tol, max_iter);
 //   boost::math::tools::bracket_and_solve_root(f, x, x, tol, max_iter);
 //   res.row(i) = x.t();
 // }
 return res;
}

//' @export
//' @rdname watson
// [[Rcpp::export]]
arma::mat rwatson(unsigned int n, const arma::rowvec &mu, double kappa) {
 using distr = anima::WatsonDistribution;
 distr watsonDistr;
 distr::SampleType samples(n, 3);
 watsonDistr.SetMeanAxis(mu);
 watsonDistr.SetConcentrationParameter(kappa);
 distr::GeneratorType generator(std::time(0));
 watsonDistr.Random(samples, generator);
 return samples;
}

// [[Rcpp::export]]
arma::rowvec mean_watson_impl(const arma::mat &x) {
  using distr = anima::WatsonDistribution;
  distr watsonDistr;
  watsonDistr.Fit(x, "");
  return watsonDistr.GetMeanAxis();
}
