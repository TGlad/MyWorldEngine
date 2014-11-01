#include "StandardIncludes.h"
#include "Vectors.h"

class SpatialMatrix : public Matrix<double, 3, 3>
{
  double scale() // return the determinant?
  {
    return (row(0).cross(row(1))).dot(row(2));
  }
};

class LorentzMatrix : public Matrix<double, 4, 4>
{
  LorentzMatrix(const Vector3d &velocity)
  {
    setZero();
    boost(velocity);
  }
  LorentzMatrix(const Quaternion &orientation, const Vector3d &velocity)
  {
    setZero();
    block<3,3,>(1,1) = orientation.toMatrix3();
    boost(velocity);
  }
  LorentzMatrix(const Quaternion &orientation)
  {
    setZero();
    block<3,3,>(1,1) = orientation.toMatrix3();
  }
  void boost(const Vector3d &vel);
  double scale(){ return spatialMatrix().scale() / massEnergy(); }
  double massEnergy(){ return fourVector().norm(); }
  FourVector fourVector()
  {
    return (FourVector)block<4,1>(0,0);
  }
  SpatialMatrix spatial()
  {
    return (SpatialMatrix)block<3,3>(1,1);
  }
};

class MinkowskiMatrix : public Matrix<double, 5, 5>
{
  MinkowskiMatrix(const Event &eventVec, LorentzMatrix &lorentz)
  {
    setZero();
    lorentzMatrix() = lorentz;
    event() = eventVec;
  }
  Event &event()
  {
    return (Event &)block<4,1>(0,4);
  }
  LorentzMatrix &lorentz()
  {
    return (LorentzMatrix &)block<4,4>(0, 0);
  }
  FourVector operator *(const FourVector &fourVec) const
  {
    Vector5d res = *this * Vector5d(fourVec[0], fourVec[1], fourVec[2], fourVec[3], 1.0); 
    return FourVector(res[0], res[1], res[2], res[3]);
  }

  StateVector transformState(const StateVector &localState);
  Vector3d transformPosition(const Vector3d &localPosition);
  Vector3d transformDelta(const Vector3d &localOffset);
  Vector3d transformAngularDelta(const Vector3d &localOffset);
};