#include "animaUniformDistribution.h"
#include "animaVectorOperations.h"

namespace anima
{

double UniformDistribution::GetDensity(const ValueType &x)
{
  if (!this->BelongsToSupport(x))
    return 0.0;
  return std::exp(this->GetLogDensity(x));
}

double UniformDistribution::GetLogDensity(const ValueType &x)
{
  if (!this->BelongsToSupport(x))
    cpp11::stop("The log-density is not defined outside the support.");

  return -std::log(4.0 * M_PI);
}

double UniformDistribution::GetCumulative(const ValueType &x)
{
  if (!this->BelongsToSupport(x))
    cpp11::stop("The CDF is not defined outside the support.");

  ValueType sphCoords;
  anima::TransformCartesianToSphericalCoordinates(x, sphCoords);
  double thetaVal = sphCoords[0];
  while (thetaVal > M_PI)
    thetaVal -= (2.0 * M_PI);
  while (thetaVal < 0)
    thetaVal += (2.0 * M_PI);
  double phiVal = sphCoords[1];
  while (phiVal > 2.0 * M_PI)
    phiVal -= (2.0 * M_PI);
  while (phiVal < 0)
    phiVal += 2.0 * M_PI;

  return phiVal / (2.0 * M_PI) * (1.0 - std::cos(thetaVal)) / 2.0;
}

void UniformDistribution::Random(SampleType &sample, GeneratorType &generator)
{
  RealUniformDistributionType unifDistr(0.0, 1.0);
  ValueType sphCoords, carCoords;
  sphCoords[2] = 1.0;
  unsigned int nSamples = sample.rows();
  for (unsigned int i = 0; i < nSamples; ++i)
  {
    sphCoords[0] = 2.0 * std::asin(std::sqrt(unifDistr(generator)));
    sphCoords[1] = 2.0 * M_PI * unifDistr(generator);
    anima::TransformSphericalToCartesianCoordinates(sphCoords, carCoords);
    sample.row(i) = carCoords;
  }
}

UniformDistribution::ValueType UniformDistribution::GetMean()
{
  ValueType meanValue;
  meanValue.fill(0.0);
  return meanValue;
}

double UniformDistribution::GetDistance(Self *otherDistribution)
{
  /**
   * \fn double UniformDistribution::GetDistance(Self *otherDistribution)
   *
   * \author Aymeric Stamm
   * \date November 2023
   *
   * \param otherDistribution A pointer specifying another object of class `UniformDistribution`.
   *
   * \return A numeric value storing the symmetric Kullback-Leibler divergence with the
   * input spherical uniform distribution.
   */

  return 0.0;
}

} // end of namespace anima
