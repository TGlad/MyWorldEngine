#include "StandardIncludes.h"

struct StateVector
{
  StateVector(const Vector3d &pos, const Vector3d &vel) : position(pos), velocity(vel) {}
  vector3d position;
  vector3d velocity;
  StateVector operator +(const StateVector &vec)
  {
    return StateVector(position + vec.position, velocity + vec.velocity);
  }
  StateVector operator -(const StateVector &vec)
  {
    return StateVector(position - vec.position, velocity - vec.velocity);
  }
  StateVector operator /(double d)
  {
    return StateVector(position / d, velocity / d);
  }
};

// concrete classes will make things clearer, and since they are only interfaces, one can convert to raw matrix types to combine them.
// may multiply by mass or not
class FourVector : public Vector4d
{
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
  inline Vector3d &spatial(){ return block(1,0,3,0); }
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
  Event(double time, const Vector3d &position)
  {
    set(time, position[0], position[1], position[2]);
  }
  inline Vector3d &position(){ return block(1,0,3,0); }
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
  inline double &lorentzFactor(){ return (*this)[0]; }
  inline Vector3d &properVelocity(){ return block(1,0,3,0); }
};
