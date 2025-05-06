#pragma once

#include "animaBaseDistribution.h"

#include <boost/math/distributions/beta.hpp>

namespace anima
{

class VonMisesFisherDistribution : public BaseDistribution
{
public:
  using BetaDistributionType = boost::math::beta_distribution<double>;

  VonMisesFisherDistribution()
  {
    m_MeanDirection[0] = 0;
    m_MeanDirection[1] = 0;
    m_MeanDirection[2] = 1;
    m_ConcentrationParameter = 1.0;
    m_BesselRatio = 1.0;
    m_RealUniformDistribution = RealUniformDistributionType(0.0, 1.0);
    m_BetaDistribution = BetaDistributionType(1.0, 1.0);
  }

  double GetDensity(const ValueType &x);
  double GetLogDensity(const ValueType &x);
  double GetCumulative(const ValueType &x);
  void Fit(const SampleType &sample, const std::string &method);
  void Random(SampleType &sample, GeneratorType &generator);
  ValueType GetMean() { return m_BesselRatio * m_MeanDirection; }
  double GetVariance() { return this->GetCovarianceMatrix().trace(); }
  double GetDistance(Self *otherDistribution);

  void SetMeanDirection(const ValueType &x);
  ValueType GetMeanDirection() { return m_MeanDirection; }

  void SetConcentrationParameter(const double &x);
  double GetConcentrationParameter() { return m_ConcentrationParameter; }

  double GetBesselRatio() { return m_BesselRatio; }
  RotationMatrixType GetCovarianceMatrix();

private:
  void SampleFromVMFDistribution(ValueType &resVec, GeneratorType &generator);
  void SampleFromVMFDistributionNumericallyStable(ValueType &resVec, GeneratorType &generator);
  ValueType m_MeanDirection;
  double m_ConcentrationParameter;
  double m_BesselRatio;
  RotationMatrixType m_NorthToMeanAxisRotationMatrix;
  RealUniformDistributionType m_RealUniformDistribution;
  BetaDistributionType m_BetaDistribution;
  const unsigned int m_AmbientDimension = 3;
};

} // end of namespace anima
