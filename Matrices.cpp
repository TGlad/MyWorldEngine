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

