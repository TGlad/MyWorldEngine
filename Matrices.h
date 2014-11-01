#include "StandardIncludes.h"

typedef Matrix<double, 4, 4> PMatrix;

class SpatialMatrix : public Matrix<double, 3, 3>
{
};

class LorentzMatrix : public Matrix<double, 4, 4>
{
  LorentzMatrix(const Quaternion &orientation, const Vector3d &velocity)
  {
    block(0,0,3,3) = orientation.toMatrix3();
    boost(velocity);
  }
  void boost(const Vector3d &vel)
  {
  }
  FourVector fourVector()
  {
    return (FourVector)block(4,0,4,4);
  }
  SpatialMatrix spatialMatrix()
  {
    return (SpatialMatrix)block(0,0,3,3);
  }
};

class MinkowskiMatrix : public Matrix<double, 5, 5>
{
  MinkowskiMatrix(const Event &event, LorentzMatrix &lorentz)
  {
    lorentzMatrix() = lorentz;
  }
  Event &event()
  {
    return (Event &)block(4,0,4,4);
  }
  LorentzMatrix &lorentzMatrix()
  {
    return (LorentzMatrix &)block(0, 0, 4, 4);
  }
};