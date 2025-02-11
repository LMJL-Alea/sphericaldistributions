#pragma once

#include <cmath>

namespace anima
{

class DawsonIntegrand
{
public:
  void SetXValue(double val) { m_XValue = val; }

  double operator()(const double t)
  {
    double inexpValue = m_XValue * m_XValue * (t * t - 1.0);
    return std::exp(inexpValue);
  }

private:
  double m_XValue;
};

double EvaluateDawsonIntegral(const double x, const bool scaled = false);
double EvaluateDawsonFunctionNR(double x);
double EvaluateDawsonFunction(double x);
double EvaluateWImFunction(double x);
double EvaluateWImY100Function(double y100, double x);

} // end namespace anima
