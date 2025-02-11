#pragma once

#include <RcppArmadillo.h>

namespace anima {

template <class VectorType>
void pairingToQuaternion(const VectorType &inputPoint, const VectorType &inputTransformedPoint, arma::mat44 &outputMatrix);

template <class VectorType>
arma::mat33 computeRotationFromQuaternion(const VectorType &eigenVector);

template <class VectorType>
arma::mat33 GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction);

const unsigned int m_AmbientDimension = 3;

} // end namespace anima

#include "animaRotationOperations.hpp"
