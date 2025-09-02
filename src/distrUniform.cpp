#include "animaUniformDistribution.h"

#include <cpp11.hpp>
#include <cpp11eigen.hpp>

#include <Eigen/Core>

// directive for openMP
#ifdef _OPENMP
#include <omp.h>
#endif

[[cpp11::register]]
cpp11::doubles dunisph_impl(const cpp11::doubles_matrix<> &x,
                            bool log = false) {
  Eigen::MatrixX3d xr = as_Matrix(x);

  using distr = anima::UniformDistribution;
  distr uniDistr;

  unsigned int n = xr.rows();
  Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

  for (unsigned int i = 0; i < n; ++i) {
    if (log)
      res(i) = uniDistr.GetLogDensity(xr.row(i));
    else
      res(i) = uniDistr.GetDensity(xr.row(i));
  }
  return as_doubles(res);
}

[[cpp11::register]]
cpp11::doubles punisph_impl(const cpp11::doubles_matrix<> &x) {
  Eigen::MatrixX3d xr = as_Matrix(x);

  using distr = anima::UniformDistribution;
  distr uniDistr;

  unsigned int n = xr.rows();
  Eigen::VectorXd res(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

  for (unsigned int i = 0; i < n; ++i)
    res(i) = uniDistr.GetCumulative(xr.row(i));
  return as_doubles(res);
}

[[cpp11::register]]
cpp11::doubles_matrix<> runisph_impl(unsigned int n) {
  using distr = anima::UniformDistribution;
  distr uniDistr;

  distr::SampleType samples(n, 3);
  distr::GeneratorType generator(std::time(0));
  uniDistr.Random(samples, generator);

  return as_doubles_matrix(samples);
}
