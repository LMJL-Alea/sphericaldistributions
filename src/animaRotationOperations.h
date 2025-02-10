#include <RcppArmadillo.h>

namespace anima {

void pairingToQuaternion(const arma::vec3 &inputPoint, const arma::vec3 &inputTransformedPoint, arma::mat44 &outputMatrix);
arma::mat33 computeRotationFromQuaternion(const arma::vec4 &eigenVector);
arma::mat33 GetRotationMatrixFromVectors(const arma::vec3 &first_direction, const arma::vec3 &second_direction);

const unsigned int m_AmbientDimension = 3;

} // end namespace anima
