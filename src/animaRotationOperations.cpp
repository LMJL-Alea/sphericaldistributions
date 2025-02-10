#include "animaRotationOperations.h"

namespace anima {

void pairingToQuaternion(const arma::vec3 &inputPoint, const arma::vec3 &inputTransformedPoint, arma::mat44 &outputMatrix)
{
  outputMatrix.fill(0.0);

  for (unsigned int i = 0;i < m_AmbientDimension;++i)
  {
    switch (i)
    {
    case 0:
      outputMatrix(0,1) = inputPoint[i] - inputTransformedPoint[i];
      outputMatrix(2,3) = - (inputPoint[i] + inputTransformedPoint[i]);
      outputMatrix(1,0) = - outputMatrix(0,1);
      outputMatrix(3,2) = - outputMatrix(2,3);
      break;

    case 1:
      outputMatrix(0,2) = inputPoint[i] - inputTransformedPoint[i];
      outputMatrix(1,3) = inputPoint[i] + inputTransformedPoint[i];
      outputMatrix(2,0) = - outputMatrix(0,2);
      outputMatrix(3,1) = - outputMatrix(1,3);
      break;

    case 2:
      outputMatrix(0,3) = inputPoint[i] - inputTransformedPoint[i];
      outputMatrix(1,2) = - (inputPoint[i] + inputTransformedPoint[i]);
      outputMatrix(2,1) = - outputMatrix(1,2);
      outputMatrix(3,0) = - outputMatrix(0,3);
      break;

    default:
      break;
    }
  }
}

arma::mat33 computeRotationFromQuaternion(const arma::vec4 &eigenVector)
{
  double normVector = 0;

  for (unsigned int i = 0;i < eigenVector.size();++i)
    normVector += eigenVector[i]*eigenVector[i];

  arma::mat33 resVal;
  resVal.fill(0.0);

  resVal(0,0) = (1.0/normVector) * (eigenVector[0]*eigenVector[0] + eigenVector[1]*eigenVector[1] - eigenVector[2]*eigenVector[2] - eigenVector[3]*eigenVector[3]);
  resVal(0,1) = (2.0/normVector) * (eigenVector[1]*eigenVector[2] - eigenVector[0]*eigenVector[3]);
  resVal(0,2) = (2.0/normVector) * (eigenVector[1]*eigenVector[3] + eigenVector[0]*eigenVector[2]);

  resVal(1,0) = (2.0/normVector) * (eigenVector[1]*eigenVector[2] + eigenVector[0]*eigenVector[3]);
  resVal(1,1) = (1.0/normVector) * (eigenVector[0]*eigenVector[0] - eigenVector[1]*eigenVector[1] + eigenVector[2]*eigenVector[2] - eigenVector[3]*eigenVector[3]);
  resVal(1,2) = (2.0/normVector) * (eigenVector[2]*eigenVector[3] - eigenVector[0]*eigenVector[1]);

  resVal(2,0) = (2.0/normVector) * (eigenVector[1]*eigenVector[3] - eigenVector[0]*eigenVector[2]);
  resVal(2,1) = (2.0/normVector) * (eigenVector[2]*eigenVector[3] + eigenVector[0]*eigenVector[1]);
  resVal(2,2) = (1.0/normVector) * (eigenVector[0]*eigenVector[0] - eigenVector[1]*eigenVector[1] - eigenVector[2]*eigenVector[2] + eigenVector[3]*eigenVector[3]);

  return resVal;
}

arma::mat33 GetRotationMatrixFromVectors(const arma::vec3 &first_direction, const arma::vec3 &second_direction)
{
  arma::mat44 tmpMatrix;
  tmpMatrix.fill(0.0);
  anima::pairingToQuaternion(first_direction, second_direction, tmpMatrix);
  arma::mat44 AMatrix = tmpMatrix.t() * tmpMatrix;

  // Needs two points, providing one on normal vector (cross product)
  arma::vec3 tmpVec;
  tmpVec[0] = first_direction[1] * second_direction[2] - first_direction[2] * second_direction[1];
  tmpVec[1] = first_direction[2] * second_direction[0] - first_direction[0] * second_direction[2];
  tmpVec[2] = first_direction[0] * second_direction[1] - first_direction[1] * second_direction[0];

  double normTmpVec = 0;
  for (unsigned int i = 0;i < m_AmbientDimension;++i)
    normTmpVec += tmpVec[i] * tmpVec[i];

  normTmpVec = std::sqrt(normTmpVec);
  if (normTmpVec < 1.0e-8)
  {
    tmpVec[0] = 0;
    tmpVec[1] = 0;
    tmpVec[2] = 1;
  }

  anima::pairingToQuaternion(tmpVec, tmpVec, tmpMatrix);
  AMatrix += tmpMatrix.t() * tmpMatrix;

  arma::mat44 eVecs;
  arma::vec4 eVals;
  arma::eig_sym(eVals, eVecs, AMatrix);

  arma::mat33 rotationMatrix = anima::computeRotationFromQuaternion(eVecs.row(0));

  return rotationMatrix;
}

} // end namespace anima
