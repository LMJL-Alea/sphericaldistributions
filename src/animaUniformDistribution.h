#pragma once

#include "animaBaseDistribution.h"

namespace anima
{

class UniformDistribution : public BaseDistribution
{
public:

  UniformDistribution() {}

  bool BelongsToSupport(const ValueType &x);
  double GetDensity(const ValueType &x);
  double GetLogDensity(const ValueType &x);
  double GetCumulative(const ValueType &x);
  void Fit(const SampleType &sample, const std::string &method) { return; }
  void Random(SampleType &sample, GeneratorType &generator);
  ValueType GetMean();
  double GetVariance() { return 0.0; }
  double GetDistance(Self *otherDistribution);
};

} // end of namespace
