#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace anima {

using RotationMatrixType = Eigen::Matrix3d;

template <class VectorType>
RotationMatrixType GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction)
{
  Eigen::Quaterniond qVector;
  qVector.setFromTwoVectors(first_direction, second_direction);
  return qVector.toRotationMatrix();
}

} // end namespace anima
