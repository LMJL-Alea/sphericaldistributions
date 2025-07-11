#include "animaVonMisesFisherDistribution.h"
#include "animaBesselFunctions.h"
#include "animaVectorOperations.h"
#include "animaRotationOperations.h"

#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace anima
{

void VonMisesFisherDistribution::SetMeanDirection(const ValueType &val)
{
  if (!this->BelongsToSupport(val))
  {
    std::string msg = "The mean direction parameter of the von Mises Fisher distribution should be a unit vector but its norm is ";
    msg += std::to_string(val.norm());
    msg += ", and the mean direction parameter is ";
    msg += std::to_string(val[0]);
    msg += ", ";
    msg += std::to_string(val[1]);
    msg += ", ";
    msg += std::to_string(val[2]);
    msg += ".";
    cpp11::message(msg.c_str());
    cpp11::stop("The mean direction parameter of the von Mises Fisher distribution should be of unit norm.");
  }


  m_MeanDirection[0] = 0.0;
  m_MeanDirection[1] = 0.0;
  m_MeanDirection[2] = 1.0;
  // Compute rotation matrix to bring [0,0,1] on meanAxis
  m_NorthToMeanAxisRotationMatrix = anima::GetRotationMatrixFromVectors(m_MeanDirection, val);
  m_MeanDirection = val;
}

void VonMisesFisherDistribution::SetConcentrationParameter(const double &val)
{
  if (val < 0)
    cpp11::stop("The concentration parameter of the von Mises Fisher distribution should be non-negative.");
  m_ConcentrationParameter = val;
  m_BesselRatio = anima::bessel_ratio_i_lower_bound(val, static_cast<double>(m_AmbientDimension) / 2.0);
}

double VonMisesFisherDistribution::GetDensity(const ValueType &x)
{
  if (!this->BelongsToSupport(x))
    return 0.0;

  double kappa = this->GetConcentrationParameter();
  ValueType meanDirection = this->GetMeanDirection();

  if (kappa < std::sqrt(this->GetEpsilon()))
    return std::exp(kappa * meanDirection.dot(x)) / (4.0 * M_PI);

  double tmpVal = kappa * (meanDirection.dot(x) - 1.0);
  double resVal = std::exp(tmpVal);
  resVal *= kappa;
  tmpVal = 1.0 - std::exp(-2.0 * kappa);
  resVal /= (2.0 * M_PI * tmpVal);

  return resVal;
}

double VonMisesFisherDistribution::GetLogDensity(const ValueType &x)
{
  if (!this->BelongsToSupport(x))
    cpp11::stop("The log-density of the von Mises Fisher distribution is not defined for arguments outside the 2-sphere.");

  return std::log(this->GetDensity(x));
}

double VonMisesFisherDistribution::GetCumulative(const ValueType &x)
{
  if (!this->BelongsToSupport(x))
    cpp11::stop("The CDF is not defined outside the support.");

  ValueType carCoords, sphCoords;
  carCoords.fill(0.0);
  for (unsigned int i = 0; i < m_AmbientDimension; ++i)
    for (unsigned int j = 0; j < m_AmbientDimension; ++j)
      carCoords[i] += m_NorthToMeanAxisRotationMatrix(j, i) * x[j];
  anima::TransformCartesianToSphericalCoordinates(carCoords, sphCoords);

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

  double phiCumul = phiVal / (2.0 * M_PI);
  double thetaCumul = (1.0 - std::exp(-m_ConcentrationParameter * (1.0 - std::cos(thetaVal))));
  thetaCumul /= (1.0 - std::exp(-2.0 * m_ConcentrationParameter));

  return phiCumul * thetaCumul;
}

void VonMisesFisherDistribution::Fit(const SampleType &sample, const std::string &method)
{
  /**********************************************************************************************
   * \fn      void VonMisesFisherDistribution::Fit(std::vector<arma::vec3>,
   *                                               std::mt19937 &generator)
   *
   * \brief	Closed-form approximations of the maximum likelihood estimators for the mean
   *          direction and concentration parameter of the von Mises Fisher distribution using
   *          the procedure described in Sra, A short note on parameter approximation for von
   *          Mises-Fisher distributions: and a fast implementation of Is(x), Computational
   *          Statistics, 2011.
   *
   * \author	Aymeric Stamm
   * \date	October 2023
   *
   * \param	sample A numeric matrix of shape n x 3 specifying a sample of size n drawn from the
   *          von Mises Fisher distribution on the 2-sphere.
   * \param	method A string specifying the estimation method. Unused here.
   **********************************************************************************************/

  unsigned int numberOfObservations = sample.rows();
  double ambientDimension = static_cast<double>(m_AmbientDimension);

  // Eq. (2)
  ValueType meanDirection;
  meanDirection.fill(0.0);
  for (unsigned int i = 0; i < numberOfObservations; ++i)
    for (unsigned int j = 0; j < m_AmbientDimension; ++j)
      meanDirection[j] += sample(i, j);
  double normValue = meanDirection.norm();
  meanDirection /= normValue;
  double resultantValue = normValue / static_cast<double>(numberOfObservations);

  if (std::abs(resultantValue - 1.0) < this->GetEpsilon())
  {
    // All observations are perfectly aligned
    // Kappa should be close to 0.
    this->SetMeanDirection(meanDirection);
    this->SetConcentrationParameter(this->GetEpsilon());
    return;
  }

  // Eq. (4)
  double concentrationParameter = resultantValue * (ambientDimension - resultantValue * resultantValue);
  concentrationParameter /= (1.0 - resultantValue * resultantValue);

  // Eq. (6)
  bool continueLoop = true;
  while (continueLoop)
  {
    double oldConcentrationParameter = concentrationParameter;
    double besselRatio = anima::bessel_ratio_i_lower_bound(concentrationParameter, ambientDimension / 2.0);
    double tmpValue = besselRatio - resultantValue;
    tmpValue /= (1.0 - besselRatio * besselRatio - (ambientDimension - 1.0) / concentrationParameter * besselRatio);
    concentrationParameter -= tmpValue;
    if (std::abs(concentrationParameter - oldConcentrationParameter) < this->GetEpsilon())
      continueLoop = false;
  }

  this->SetMeanDirection(meanDirection);
  this->SetConcentrationParameter(concentrationParameter);
}

void VonMisesFisherDistribution::Random(SampleType &sample, GeneratorType &generator)
{
  /**********************************************************************************************
   * \fn void Random(std::vector<arma::vec3>, std::mt19937 &generator)
   *
   * \brief	Sample from the Watson distribution using the procedure described in Fisher et al.,
   *          Statistical Analysis of Spherical Data, Cambridge University Press, 1993, pp. 59.
   *
   * \author	Aymeric Stamm
   * \date	October 2013
   *
   * \param	sample    A numeric matrix of shape n x 3 storing a sample of size n drawn from the
   *                    Watson distribution on the 2-sphere.
   * \param	generator A pseudo-random number generator.
   **********************************************************************************************/

  unsigned int sampleSize = sample.rows();
  ValueType sampleValue;

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
#endif

  for (unsigned int i = 0; i < sampleSize; ++i)
  {
    if (m_ConcentrationParameter > 700.0)
      this->SampleFromVMFDistributionNumericallyStable(sampleValue, generator);
    else
      this->SampleFromVMFDistribution(sampleValue, generator);
    sample.row(i) = sampleValue;
  }
}

double VonMisesFisherDistribution::GetDistance(Self *otherDistribution)
{
  /**
   * \fn double VonMisesFisherDistribution::GetDistance(VonMisesFisherDistribution *otherDistribution)
   *
   * \author Aymeric Stamm
   * \date November 2023
   *
   * \param otherDistribution An object of class `VonMisesFisherDistribution`.
   *
   * \return A numeric value storing the symmetric Kullback-Leibler divergence to the input von Mises
   * Fisher distribution. This is achieved following Kitagawa & Rowley (2022), von Mises-Fisher
   * distributions and their statistical divergence, arXiv:2202.05192v1
   * (https://arxiv.org/pdf/2202.05192v1.pdf).
   */

  VonMisesFisherDistribution *vmfDistr = dynamic_cast<VonMisesFisherDistribution *>(otherDistribution);
  ValueType otherMeanDirection = vmfDistr->GetMeanDirection();

  double otherConcentrationParameter = vmfDistr->GetConcentrationParameter();
  double otherBesselRatio = vmfDistr->GetBesselRatio();

  double thisToOtherDist = 0.0;
  for (unsigned int i = 0; i < m_AmbientDimension; ++i)
    thisToOtherDist += (m_ConcentrationParameter * m_MeanDirection[i] - otherConcentrationParameter * otherMeanDirection[i]) * m_MeanDirection[i];
  thisToOtherDist *= m_BesselRatio;

  double otherToThisDist = 0.0;
  for (unsigned int i = 0; i < m_AmbientDimension; ++i)
    otherToThisDist += (otherConcentrationParameter * otherMeanDirection[i] - m_ConcentrationParameter * m_MeanDirection[i]) * otherMeanDirection[i];
  otherToThisDist *= otherBesselRatio;

  return thisToOtherDist + otherToThisDist;
}

RotationMatrixType VonMisesFisherDistribution::GetCovarianceMatrix()
{
  /**
   * \fn arma::mat33 VonMisesFisherDistribution::GetCovarianceMatrix()
   *
   * \author Aymeric Stamm
   * \date November 2023
   *
   * \return A numeric matrix of size `m_AmbientDimension x m_AmbientDimension` storing the
   * covariance matrix of the von Mises Fisher distribution. This is achieved following
   * Kitagawa & Rowley (2022), von Mises-Fisher distributions and their statistical divergence,
   * arXiv:2202.05192v1 (https://arxiv.org/pdf/2202.05192v1.pdf).
   */

  RotationMatrixType covarianceMatrix;
  double diagConstant = m_BesselRatio / m_ConcentrationParameter;
  double offDiagConstant = 1.0 - static_cast<double>(m_AmbientDimension) * m_BesselRatio / m_ConcentrationParameter - m_BesselRatio * m_BesselRatio;

  for (unsigned int i = 0; i < m_AmbientDimension; ++i)
  {
    covarianceMatrix(i, i) = diagConstant + offDiagConstant * m_MeanDirection[i] * m_MeanDirection[i];

    for (unsigned int j = i + 1; j < m_AmbientDimension; ++j)
    {
      double tmpValue = offDiagConstant * m_MeanDirection[i] * m_MeanDirection[j];
      covarianceMatrix(i, j) = tmpValue;
      covarianceMatrix(j, i) = tmpValue;
    }
  }

  return covarianceMatrix;
}

void VonMisesFisherDistribution::SampleFromVMFDistribution(ValueType &resVec, GeneratorType &generator)
{
  /**
   * \fn void VonMisesFisherDistribution::SampleFromVMFDistribution(arma::vec3 &resVec, &generator)
   *
   * \brief Sample from the VMF distribution following Ulrich, G. (1984). Computer generation of
   * distributions on the m‐sphere. Journal of the Royal Statistical Society: Series C (Applied
   * Statistics), 33(2), 158-163.
   *
   * \author Aymeric Stamm
   * \date October 2013
   *
   * \param resVec An object of type arma::vec3 that will store a value sampled from
   * the VMF distribution.
   * \param generator An object of type std::mt19937 specifying which random number generator to
   * use.
   */

  double tmpVal = std::sqrt(m_ConcentrationParameter * m_ConcentrationParameter + 1.0);
  double b = (-2.0 * m_ConcentrationParameter + 2.0 * tmpVal) / 2.0;
  double a = (1.0 + m_ConcentrationParameter + tmpVal) / 2.0;
  double d = 4.0 * a * b / (1.0 + b) - 2.0 * std::log(2.0);

  double T = 1.0;
  double U = std::exp(d);
  double W = 0;

  while (2.0 * std::log(T) - T + d < std::log(U))
  {
    double Z = boost::math::quantile(m_BetaDistribution, m_RealUniformDistribution(generator));
    U = m_RealUniformDistribution(generator);
    tmpVal = 1.0 - (1.0 - b) * Z;
    T = 2.0 * a * b / tmpVal;
    W = (1.0 - (1.0 + b) * Z) / tmpVal;
  }

  double theta = 2.0 * M_PI * m_RealUniformDistribution(generator);
  resVec[0] = std::sqrt(1.0 - W * W) * std::cos(theta);
  resVec[1] = std::sqrt(1.0 - W * W) * std::sin(theta);
  resVec[2] = W;

  // Rotate to bring everything back around meanDirection
  resVec = resVec * m_NorthToMeanAxisRotationMatrix.transpose();

  double resNorm = resVec.norm();
  resVec /= resNorm;

  if (std::abs(resNorm - 1.0) > this->GetEpsilon())
  {
    std::string msg = "VMF distribution sampling failed: resNorm = ";
    msg += std::to_string(resNorm);
    msg += " != 1.0. This is a bug in the VMF sampler.";
    cpp11::message(msg.c_str());
    cpp11::stop("The VMF sampler should generate points on the 2-sphere.");
  }

}

void VonMisesFisherDistribution::SampleFromVMFDistributionNumericallyStable(ValueType &resVec, GeneratorType &generator)
{
  /**
   * \fn void VonMisesFisherDistribution::SampleFromVMFDistributionNumericallyStable(arma::vec3 &resVec, &generator)
   *
   * \brief Sample from the VMF distribution following Jakob, W. (2012). Numerically stable sampling
   * of the von Mises-Fisher distribution on Sˆ2 (and other tricks). Interactive Geometry Lab, ETH
   * Zürich, Tech. Rep, 6.
   *
   * \author Aymeric Stamm
   * \date October 2013
   *
   * \param resVec An object of type arma::vec3 that will store a value sampled from
   * the VMF distribution.
   * \param generator An object of type std::mt19937 specifying which random number generator to
   * use.
   */

  double xi = m_RealUniformDistribution(generator);
  double W = 1.0 + (std::log(xi) + std::log(1.0 - (xi - 1.0) * std::exp(-2.0 * m_ConcentrationParameter) / xi)) / m_ConcentrationParameter;
  double theta = 2.0 * M_PI * m_RealUniformDistribution(generator);

  resVec[0] = std::sqrt(1.0 - W * W) * std::cos(theta);
  resVec[1] = std::sqrt(1.0 - W * W) * std::sin(theta);
  resVec[2] = W;

  // Rotate to bring everthing back around meanDirection
  resVec = resVec * m_NorthToMeanAxisRotationMatrix.transpose();

  double resNorm = resVec.norm();
  resVec /= resNorm;

  if (std::abs(resNorm - 1.0) > this->GetEpsilon())
    cpp11::stop("The VMF sampler should generate points on the 2-sphere.");
}

} // end of namespace anima
