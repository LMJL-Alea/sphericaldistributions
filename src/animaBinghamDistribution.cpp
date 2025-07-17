#include "animaBinghamDistribution.h"
#include "animaErrorFunctions.h"
#include "animaKummerFunctions.h"
#include "animaRotationOperations.h"
#include "animaVectorOperations.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

namespace anima
{

// Using Gauss-Legendre to estimate the CDF and normalization constant.
static void gauss_legendre(int n, double a, double b, std::vector<double>& x, std::vector<double>& w)
    {
        x.resize(n);
        w.resize(n);
        // Hardcoded 16-point rule for [-1,1], then scale to [a,b]
        static const double xs[8] = {
            0.989400934991649932596154173450,
            0.944575023073232576077988415535,
            0.865631202387831743880467897712,
            0.755404408355003033895101194847,
            0.617876244402643748446671764049,
            0.458016777657227386342419442984,
            0.281603550779258913230460501460,
            0.095012509837637440185319335425
        };
        static const double ws[8] = {
            0.027152459411754094851780572456,
            0.062253523938647892862843836994,
            0.095158511682492784809925107602,
            0.124628971255533872052476282192,
            0.149595988816576732081501730547,
            0.169156519395002538189312079030,
            0.182603415044923588866763667969,
            0.189450610455068496285396723208
        };
        for (int i = 0; i < 8; ++i) {
            x[i] = 0.5 * ((b - a) * xs[i] + (b + a));
            x[n-1-i] = 0.5 * (-(b - a) * xs[i] + (b + a));
            w[i] = w[n-1-i] = 0.5 * (b - a) * ws[i];
        }
    }


BinghamDistribution::BinghamDistribution()
    {
        m_MeanDirection = {0.0, 0.0, 1.0};
        m_ConcentrationParameter = {0.0, 0.0, 0.0};
        m_OrientationMatrix.setIdentity();
        UpdateNormalizationConstant();
    }

void BinghamDistribution::SetMeanDirection(const ValueType &val)
{
    if (!this->BelongsToSupport(val))
        Rcpp::Rcerr << "The mean axis parameter of the Bingham distribution should be of unit norm." << std::endl;

    m_MeanDirection[0] = 0.0;
    m_MeanDirection[1] = 0.0;
    m_MeanDirection[2] = 1.0;
    // Compute rotation matrix to bring [0,0,1] on meanAxis
    m_NorthToMeanAxisRotationMatrix = anima::GetRotationMatrixFromVectors(m_MeanAxis, val);
    m_MeanAxis = val;
}

void BinghamDistribution::SetConcentrationParameter(const double& val)
{
    m_ConcentrationParameter = val;
}

BinghamDistribution::ValueType BinghamDistribution::GetMeanDirection() const
{
    return m_MeanDirection;
}

BinghamDistribution::ValueType BinghamDistribution::GetConcentrationParameter() const
{
    return m_ConcentrationParameter;
}

double BinghamDistribution::GetDensity(const ValueType &x)
{
    if (!BelongsToSupport(x))
        return 0.0;
    ValueType y = m_OrientationMatrix.transpose() * x;
    double exponent = m_ConcentrationParameter[0]*y[0]*y[0] + m_ConcentrationParameter[1]*y[1]*y[1] + m_ConcentrationParameter[2]*y[2]*y[2];
    return std::exp(exponent) / m_NormalizationConstant;
}

double BinghamDistribution::GetLogDensity(const ValueType &x)
{
    if (!BelongsToSupport(x))
        throw std::runtime_error("Log-density not defined outside the sphere.");
    ValueType y = m_OrientationMatrix.transpose() * x;
    double exponent = m_ConcentrationParameter[0]*y[0]*y[0] + m_ConcentrationParameter[1]*y[1]*y[1] + m_ConcentrationParameter[2]*y[2]*y[2];
    return exponent - std::log(m_NormalizationConstant);
}

double BinghamDistribution::GetCumulative(const ValueType &x)
{
    if (!BelongsToSupport(x))
        throw std::runtime_error("CDF not defined outside the sphere.");

    // Convert x to spherical coordinates (theta, phi)
    anima::Vector3 sphCoords;
    anima::TransformCartesianToSphericalCoordinates(x, sphCoords);
    double theta_max = sphCoords[0];
    double phi_max = sphCoords[1];
    // Ensure theta in [0, pi], phi in [0, 2pi]
    while (theta_max > M_PI)
        theta_max -= 2.0 * M_PI;
    while (theta_max < 0)
        theta_max += 2.0 * M_PI;
    while (phi_max > 2.0 * M_PI)
        phi_max -= 2.0 * M_PI;
    while (phi_max < 0)
        phi_max += 2.0 * M_PI;

    // Integrate density over (theta in [0, theta_max], phi in [0, phi_max])
    return ComputeCumulativeIntegral(theta_max, phi_max) / m_NormalizationConstant;
}

double BinghamDistribution::ComputeNormalizationConstant(int n_theta, int n_phi) const
{
    std::vector<double> theta, w_theta;
    gauss_legendre(n_theta, 0.0, M_PI, theta, w_theta);
    double dphi = 2.0 * M_PI / n_phi;
    double norm_const = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        double th = theta[i];
        double sth = std::sin(th);
        double cth = std::cos(th);
        for (int j = 0; j < n_phi; ++j) {
            double phi = j * dphi;
            double cphi = std::cos(phi);
            double sphi = std::sin(phi);
            // Rotate to principal axes
            ValueType x;
            x[0] = sth * cphi;
            x[1] = sth * sphi;
            x[2] = cth;
            x = m_OrientationMatrix * x;
            double expo = m_ConcentrationParameter[0]*x[0]*x[0] + m_ConcentrationParameter[1]*x[1]*x[1] + m_ConcentrationParameter[2]*x[2]*x[2];
            norm_const += std::exp(expo) * sth * w_theta[i] * dphi;
        }
    }
    return norm_const;
}

double BinghamDistribution::ComputeCumulativeIntegral(double theta_max, double phi_max, int n_theta, int n_phi) const
{
    std::vector<double> theta, w_theta;
    gauss_legendre(n_theta, 0.0, theta_max, theta, w_theta);
    double dphi = phi_max / n_phi;
    double cum = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        double th = theta[i];
        double sth = std::sin(th);
        double cth = std::cos(th);
        for (int j = 0; j < n_phi; ++j) {
            double phi = j * dphi;
            double cphi = std::cos(phi);
            double sphi = std::sin(phi);
            // Rotate to principal axes
            ValueType x;
            x[0] = sth * cphi;
            x[1] = sth * sphi;
            x[2] = cth;
            x = m_OrientationMatrix * x;
            double expo = m_ConcentrationParameter[0]*x[0]*x[0] + m_ConcentrationParameter[1]*x[1]*x[1] + m_ConcentrationParameter[2]*x[2]*x[2];
            cum += std::exp(expo) * sth * w_theta[i] * dphi;
        }
    }
    return cum;
}

double BinghamDistribution::GetCumulative(const ValueType &x)
{
    if (!BelongsToSupport(x))
        throw std::runtime_error("CDF not defined outside the sphere.");

    // Convert x to spherical coordinates (theta, phi)
    ValueType sphCoords;
    anima::TransformCartesianToSphericalCoordinates(x, sphCoords);
    double theta_max = sphCoords[0];
    double phi_max = sphCoords[1];
    // Ensure theta in [0, pi], phi in [0, 2pi]
    while (theta_max > M_PI)
        theta_max -= (2.0 * M_PI);
    while (theta_max < 0)
        theta_max += (2.0 * M_PI);
    while (phi_max > 2.0 * M_PI)
        phi_max -= (2.0 * M_PI);
    while (phi_max < 0)
        phi_max += 2.0 * M_PI;

    // Integrate density over (theta in [0, theta_max], phi in [0, phi_max])
    return ComputeCumulativeIntegral(theta_max, phi_max) / m_NormalizationConstant;
}

void BinghamDistribution::Fit(const SampleType &sample, const std::string &method)
{
    // Moment-based estimation
    MatrixType S;
    S.setZero();
    for (const auto& x : sample)
        S += anima::OuterProduct(x, x);
    S /= static_cast<double>(sample.size());

    // Eigen decomposition
    ValueType eigval;
    MatrixType eigvec;
    Eigen::Matrix3d S_eigen;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            S_eigen(i, j) = S(i, j); // Copy your matrix to Eigen

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(S_eigen);
    if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Eigen decomposition failed!");
        }

    // Eigenvalues (sorted in increasing order)
    Eigen::Vector3d evals = eigensolver.eigenvalues();
    // Eigenvectors: columns are eigenvectors
    Eigen::Matrix3d evecs = eigensolver.eigenvectors();

    // Copy back to your types
    for (int i = 0; i < 3; ++i) eigval[i] = evals(i);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            eigvec(i, j) = evecs(i, j);

    SetMeanDirection(eigvec.getColumn(2));
    ValueType z;
    z[0] = std::log(eigval[0]/eigval[2]);
    z[1] = std::log(eigval[1]/eigval[2]);
    z[2] = 0.0;
    SetConcentrationParameter(z);
}

void BinghamDistribution::Random(SampleType &sample, GeneratorType &generator)
{
    // Rejection sampling from uniform sphere
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    for (auto& x : sample)
    {
        while (true)
        {
            double u = unif(generator);
            double v = unif(generator);
            double theta = 2.0 * M_PI * u;
            double phi = std::acos(2.0 * v - 1.0);
            ValueType y;
            y[0] = std::sin(phi) * std::cos(theta);
            y[1] = std::sin(phi) * std::sin(theta);
            y[2] = std::cos(phi);
            double d = GetDensity(y);
            double M = 1.0 / m_NormalizationConstant;
            if (unif(generator) < d / M)
            {
                x = y;
                break;
            }
        }
    }
}

BinghamDistribution::ValueType BinghamDistribution::GetMean() const
{
    ValueType meanValue;
    meanValue.fill(0.0);
    return meanValue;
}

double BinghamDistribution::GetDistance(Self *otherDistribution)
{
    const unsigned int numberOfMonteCarloSamples = 10000;
    SampleType thisBinghamSample(numberOfMonteCarloSamples, 3), otherBinghamSample(numberOfMonteCarloSamples, 3);
    GeneratorType generator;
  
    this->Random(thisBinghamSample, generator);
    BinghamDistribution *binghamDistr = dynamic_cast<BinghamDistribution *>(otherDistribution);
    binghamDistr->Random(otherBinghamSample, generator);
  
    double thisKLValue = 0.0, otherKLValue = 0.0;
    for (unsigned int i = 0; i < numberOfMonteCarloSamples; ++i)
    {
      thisKLValue += this->GetLogDensity(thisBinghamSample.row(i));
      thisKLValue -= binghamDistr->GetLogDensity(thisBinghamSample.row(i));
      otherKLValue += binghamDistr->GetLogDensity(otherBinghamSample.row(i));
      otherKLValue -= this->GetLogDensity(otherBinghamSample.row(i));
    }
  
    thisKLValue /= static_cast<double>(numberOfMonteCarloSamples);
    otherKLValue /= static_cast<double>(numberOfMonteCarloSamples);
  
    return thisKLValue + otherKLValue;
  }

void BinghamDistribution::UpdateNormalizationConstant()
{
    // Approximate normalization constant for 3D Bingham with diagonal Z
    // In practice, use numerical integration or lookup tables.
    // Here, use a crude approximation:
    m_NormalizationConstant = 8.0 * M_PI * std::exp(m_ConcentrationParameter.max());
}

bool BinghamDistribution::BelongsToSupport(const ValueType& x) const
{
    return std::abs(x.norm() - 1.0) < 1e-8;
}

} // namespace anima
