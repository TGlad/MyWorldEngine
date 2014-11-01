#include "StandardIncludes.h"
#include "Vectors.h"

typedef Matrix<double, 4, 4> PMatrix;

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
    block<3,3,>(1,1) = orientation.toMatrix3();
    boost(velocity);
  }
  void boost(const Vector3d &vel)
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
  double scale(){ return spatialMatrix().scale() / massEnergy(); }
  double massEnergy(){ return fourVector().norm(); }
  FourVector fourVector()
  {
    return (FourVector)block<4,1>(0,0);
  }
  SpatialMatrix spatialMatrix()
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
  LorentzMatrix &lorentzMatrix()
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