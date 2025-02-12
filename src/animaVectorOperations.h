#pragma once

namespace anima {

template <class VectorType>
void TransformCartesianToSphericalCoordinates(const VectorType &v, VectorType &resVec);

template <class VectorType>
void TransformSphericalToCartesianCoordinates(const VectorType &v, VectorType &resVec);

} // end of namespace anima

#include "animaVectorOperations.hpp"
