#pragma once

#include <ctime>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <cpp11.hpp>
#include <cpp11eigen.hpp>

namespace anima {

class BaseDistribution {
public:
  using Self = BaseDistribution;
  using ValueType = Eigen::RowVector3d;
  using SampleType = Eigen::MatrixX3d;
  using RotationMatrixType = Eigen::Matrix3d;
  using GeneratorType = std::mt19937;
  using RealUniformDistributionType = std::uniform_real_distribution<double>;

  BaseDistribution() {}

  virtual double GetDensity(const ValueType &x) = 0;
  virtual double GetLogDensity(const ValueType &x) = 0;
  virtual double GetCumulative(const ValueType &x) = 0;
  virtual void Fit(const SampleType &sample, const std::string &method) = 0;
  virtual void Random(SampleType &sample, GeneratorType &generator) = 0;
  virtual ValueType GetMean() = 0;
  virtual double GetVariance() = 0;
  virtual double GetDistance(Self *otherDistribution) = 0;

protected:
  double GetEpsilon() {
    return std::sqrt(std::numeric_limits<double>::epsilon());
  }
  bool BelongsToSupport(const ValueType &x) {
    return std::abs(x.norm() - 1.0) < this->GetEpsilon();
  }

  void StandardizeSphericalAngles(const double &theta, const double &phi,
                                  double &theta_std, double &phi_std) {
    // Ensure theta_std is in [0, pi], phi_std is in [0, 2pi]
    theta_std = theta;
    phi_std = phi;

    theta_std = std::fmod(theta_std, 2.0 * M_PI);
    if (theta_std > M_PI) {
      theta_std = 2.0 * M_PI - theta_std;
      phi_std += M_PI;
    }

    phi_std = std::fmod(phi_std, 2.0 * M_PI);
  }
};

} // namespace anima
