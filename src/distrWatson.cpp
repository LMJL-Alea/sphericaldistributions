#include "animaWatsonDistribution.h"

#include <cpp11eigen.hpp>

#include <boost/math/tools/roots.hpp>

//' @export
//' @rdname watson
[[cpp11::register]]
cpp11::doubles dwatson_impl(const cpp11::doubles_matrix<> &x,
                            const cpp11::doubles &mu,
                            double kappa,
                            bool log = false) {
  Eigen::MatrixX3d xr = as_Matrix(x);
  Eigen::Map<Eigen::RowVector3d> mur(
      reinterpret_cast<double *>(REAL(mu.data())), 1, 3);

  using distr = anima::WatsonDistribution;
  distr watsonDistr;
  watsonDistr.SetMeanAxis(mur);
  watsonDistr.SetConcentrationParameter(kappa);

  unsigned int n = xr.rows();
  Eigen::VectorXd res(n);
  for (unsigned int i = 0;i < n;++i)
  {
    if (log)
      res(i) = watsonDistr.GetLogDensity(xr.row(i));
    else
      res(i) = watsonDistr.GetDensity(xr.row(i));
  }

  return as_doubles(res);
}

//' @export
//' @rdname watson
[[cpp11::register]]
cpp11::doubles pwatson_impl(const cpp11::doubles_matrix<> &x,
                            const cpp11::doubles &mu,
                            double kappa) {
  Eigen::MatrixX3d xr = as_Matrix(x);
  Eigen::Map<Eigen::RowVector3d> mur(
      reinterpret_cast<double *>(REAL(mu.data())), 1, 3);

  using distr = anima::WatsonDistribution;
  distr watsonDistr;
  watsonDistr.SetMeanAxis(mur);
  watsonDistr.SetConcentrationParameter(kappa);

  unsigned int n = xr.rows();
  Eigen::VectorXd res(n);
  for (unsigned int i = 0;i < n;++i)
    res(i) = watsonDistr.GetCumulative(xr.row(i));

  return as_doubles(res);
}

//' @export
//' @rdname watson
[[cpp11::register]]
cpp11::doubles_matrix<> rwatson_impl(unsigned int n,
                                     const cpp11::doubles &mu,
                                     double kappa) {
  Eigen::Map<Eigen::RowVector3d> mur(
      reinterpret_cast<double *>(REAL(mu.data())), 1, 3);

  using distr = anima::WatsonDistribution;
  distr watsonDistr;

  distr::SampleType samples(n, 3);
  watsonDistr.SetMeanAxis(mur);
  watsonDistr.SetConcentrationParameter(kappa);
  distr::GeneratorType generator(std::time(0));
  watsonDistr.Random(samples, generator);

  return as_doubles_matrix(samples);
}

[[cpp11::register]]
cpp11::doubles mean_watson_impl(const cpp11::doubles_matrix<> &x) {
  Eigen::MatrixX3d xr = as_Matrix(x);
  using distr = anima::WatsonDistribution;
  distr watsonDistr;
  watsonDistr.Fit(xr, "");
  cpp11::writable::doubles res(3);
  double* res_data = REAL(res);
  std::memcpy(res_data, watsonDistr.GetMeanAxis().data(), 3 * sizeof(double));
  return res;
}
