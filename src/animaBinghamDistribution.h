#pragma once

#include "animaBaseDistribution.h"
#include "animaVectorOperations.h"
#include "animaRotationOperations.h"
#include <random>
#include <string>
#include <vector>

namespace anima
{

    class BinghamDistribution : public BaseDistribution
    {
    public:
    
      BinghamDistribution();
    
      double GetDensity(const ValueType &x);
      double GetLogDensity(const ValueType &x);
      double GetCumulative(const ValueType &x);
      void Fit(const SampleType &sample, const std::string &method);
      void Random(SampleType &sample, GeneratorType &generator);
      ValueType GetMean();
      double GetDistance(Self *otherDistribution);
    
      void SetMeanAxis(const ValueType &x);
      ValueType GetMeanAxis() { return m_MeanAxis; }
    
      void SetConcentrationParameters(const ValueType &x) { m_ConcentrationParameters = x; }
      ValueType GetConcentrationParameters() { return m_ConcentrationParameters; }
    
    private:
      ValueType m_MeanAxis;
      ValueType m_ConcentrationParameters;
      RotationMatrixType m_OrientationMatrix; // columns are principal axes
      RotationMatrixType m_NorthToMeanAxisRotationMatrix;
      double m_NormalizationConstant;
      
      void UpdateNormalizationConstant();
      bool BelongsToSupport(const ValueType &x) const;

      // Helper for normalization and cumulative
      double ComputeNormalizationConstant(int n_theta = 16, int n_phi = 32) const;
      double ComputeCumulativeIntegral(double theta_max, double phi_max, int n_theta = 16, int n_phi = 32) const;
    };

} // end of namespace anima
    