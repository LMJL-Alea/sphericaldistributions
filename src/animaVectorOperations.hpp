#pragma once

#include "animaVectorOperations.h"

namespace anima {

template <class VectorType>
void TransformCartesianToSphericalCoordinates(const VectorType &v, VectorType &resVec)
{
  double normV = arma::norm(v);
  resVec = v / normV;

  if (resVec[2] >= 1.0)
  {
    resVec[0] = 0;
    resVec[1] = 0;
  }
  else if (resVec[2] <= -1.0)
  {
    resVec[0] = M_PI;
    resVec[1] = 0;
  }
  else
  {
    double theta = std::acos(resVec[2]);
    double phi = std::atan2(resVec[1],resVec[0]);

    if (phi < 0.0)
      phi += (2.0 * M_PI);

    resVec[0] = theta;
    resVec[1] = phi;
  }

  resVec[2] = normV;
}

template <class VectorType>
void TransformSphericalToCartesianCoordinates(const VectorType &v, VectorType &resVec)
{
  resVec[0] = v[2] * std::sin(v[0]) * std::cos(v[1]);
  resVec[1] = v[2] * std::sin(v[0]) * std::sin(v[1]);
  resVec[2] = v[2] * std::cos(v[0]);
}

} // end of namespace anima
