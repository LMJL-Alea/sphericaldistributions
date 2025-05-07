#include "animaVonMisesFisherDistribution.h"

#include <cpp11.hpp>
#include <cpp11eigen.hpp>

// directive for openMP
#ifdef _OPENMP
#include <omp.h>
#endif

[[cpp11::register]]
cpp11::doubles dvmf_impl(const cpp11::doubles_matrix<> &x,
                         const cpp11::doubles &mu,
                         double kappa,
                         bool log = false) {
  Eigen::MatrixX3d xr = as_Matrix(x);
  Eigen::Map<Eigen::RowVector3d> mur(
      reinterpret_cast<double *>(REAL(mu.data())), 1, 3);

  using distr = anima::VonMisesFisherDistribution;
  distr vmfDistr;
  vmfDistr.SetMeanDirection(mur);
  vmfDistr.SetConcentrationParameter(kappa);

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

  return as_doubles(res);
}

[[cpp11::register]]
cpp11::doubles pvmf_impl(const cpp11::doubles_matrix<> &x,
                         const cpp11::doubles &mu,
                         double kappa) {
  Eigen::MatrixX3d xr = as_Matrix(x);
  Eigen::Map<Eigen::RowVector3d> mur(
      reinterpret_cast<double *>(REAL(mu.data())), 1, 3);

  using distr = anima::VonMisesFisherDistribution;
  distr vmfDistr;
  vmfDistr.SetMeanDirection(mur);
  vmfDistr.SetConcentrationParameter(kappa);

  unsigned int n = xr.rows();
  Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

  for (unsigned int i = 0;i < n;++i)
    res(i) = vmfDistr.GetCumulative(xr.row(i));

  return as_doubles(res);
}

[[cpp11::register]]
cpp11::doubles_matrix<> rvmf_impl(unsigned int n,
                                  const cpp11::doubles &mu,
                                  double kappa) {
  Eigen::Map<Eigen::RowVector3d> mur(
      reinterpret_cast<double *>(REAL(mu.data())), 1, 3);

  using distr = anima::VonMisesFisherDistribution;
  distr vmfDistr;

  distr::SampleType samples(n, 3);
  vmfDistr.SetMeanDirection(mur);
  vmfDistr.SetConcentrationParameter(kappa);
  distr::GeneratorType generator(std::time(0));
  vmfDistr.Random(samples, generator);

  return as_doubles_matrix(samples);
}

[[cpp11::register]]
cpp11::doubles mean_vmf_impl(const cpp11::doubles_matrix<> &x) {
  Eigen::MatrixX3d xr = as_Matrix(x);
  using distr = anima::VonMisesFisherDistribution;
  distr vmfDistr;
  vmfDistr.Fit(xr, "");
  cpp11::writable::doubles res(3);
  double* res_data = REAL(res);
  std::memcpy(res_data, vmfDistr.GetMeanDirection().data(), 3 * sizeof(double));
  return res;
}
