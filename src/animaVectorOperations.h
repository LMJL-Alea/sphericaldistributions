#pragma once

#include <RcppArmadillo.h>

namespace anima {

void TransformCartesianToSphericalCoordinates(const arma::vec3 &v, arma::vec3 &resVec);
void TransformSphericalToCartesianCoordinates(const arma::vec3 &v, arma::vec3 &resVec);

} // end of namespace anima
