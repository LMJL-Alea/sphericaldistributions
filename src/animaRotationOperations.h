#pragma once

#include <Eigen/Core>

namespace anima {

using RotationMatrixType = Eigen::Matrix3d;

template <class VectorType>
RotationMatrixType GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction);

const unsigned int m_AmbientDimension = 3;

} // end namespace anima

#include "animaRotationOperations.hpp"
