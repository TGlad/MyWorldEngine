#include "Matrices.h"

// TODO: this is way more conversion than is comfortable
StateVector MinkowskiMatrix::transformState(const StateVector &localState)
{
  MinkowskiMatrix mat(Event(0.0, localState.position), LorentzMatrix(localState.velocity));
  // how do velocities transform? I guess we really want to apply special relativity velocity addition
  MinkowskiMatrix result = mat * *this;
  return StateVector(result.event().position(), result.lorentzMatrix().fourVector().velocity());
}

Vector3d MinkowskiMatrix::transformPosition(const Vector3d &localPosition)
{
  return (*this * FourVector(0, localPosition)).position();
}

Vector3d MinkowskiMatrix::transformDelta(const Vector3d &localOffset)
{
  return (lorentzMatrix() * FourVector(0, localPosition)).position();
}

Vector3d MinkowskiMatrix::transformAngularDelta(const Vector3d &localOffset)
{
  return transformDelta(localOffset); // TODO: fix. What does boost do to this, and scale shouldn't affect it
}

void MinkowskiMatrix::boost(const Vector3d &vel)
{
  Vector3d stretchDir = vel;
  stretchDir.normalize();
  Vector4d stretchAxis(1.0, stretchDir[0], stretchDir[1], stretchDir[2]);
  stretchAxis.normalize();
  Vector4d squashAxis = -stretchAxis;
  squashAxis[0] = 1.0;

  FourVelocity fourVel(vel);
  double stretch = fourVel.dot(stretchAxis) / stretchAxis[0];
  double squash = 1.0/stretch;

  for (int i = 0; i<4; i++)
  {
    col(i) += stretchAxis * ((stretch-1.0) * col(i).dot(stretchAxis));
    col(i) += squashAxis * ((squash-1.0) * col(i).dot(squashAxis));
  }
}