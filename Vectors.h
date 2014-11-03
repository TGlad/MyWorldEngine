#include "StandardIncludes.h"

// concrete classes will make things clearer, and since they are only interfaces, one can convert to raw matrix types to combine them.
// may multiply by mass or not
class FourVector : public Vector4d
{
  FourVector(double temporalVal, const Vector3d &spatialVec)
  {
    temporal() = temporalVal;
    spatial() = spatialVec;
  }
  double squaredNorm()
  {
    Vector4d temp = *this;
    temp[0] = -temp[0];
    return dot(temp);
  }
  double norm()
  {
    return sqrt(squaredNorm());
  }

  // not necessarily position and time, it could be momentum and mass, or power and force. 
  inline Vector3d &spatial(){ return block<3,1>(1,0); }
  inline double &temporal(){ return (*this)[0]; }

  // however these are reasonable to provide as they are scale invariant ratios, so give some sort of derivative
  Vector3d velocity(){ return spatial() / temporal(); }
  Vector3d properVelocity(){ return spatial() / norm(); }
  Vector3d lorentzFactor(){ return time() / norm(); }
  FourVelocity fourVelocity(){ return FourVelocity(lorentzFactor(), properVelocity()); }
}

// no mass included
class Event : public FourVector
{
  Event(double time)
  {
    set(time, 0,0,0);
  }
  Event(double time, const Vector3d &position)
  {
    set(time, position[0], position[1], position[2]);
  }
  inline Vector3d &position(){ return block<3,1>(1,0); }
  inline double &time(){ return (*this)[0]; }
  double properTime(){ return norm(); }
};

// combines proper velocity and lorentz factor. It's norm should be 1, so it is like a direction.
class FourVelocity : public FourVector
{
  FourVelocity(double lorentzFactor, const Vector3d &properVel)
  {
    set(lorentzFactor, properVel[0], properVel[1], properVel[2]);
  }
  FourVelocity(const Vector3d &velocity)
  {
    FourVector vec(1.0, velocity);
    *this = vec.fourVelocity();
  }
  inline double &lorentzFactor(){ return (*this)[0]; }
  inline Vector3d &properVelocity(){ return block<3,1>(1,0); }
};
