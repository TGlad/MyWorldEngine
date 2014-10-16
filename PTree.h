#include "standardIncludes.h"
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> PMatrix;

struct PVector
{
  Vector3d position;
  Vector3d velocity;

  PVector(const Vector3d &position, const Vector3d &velocity) : position(position), velocity(velocity) {}
  Vector6d vec(){ Vector6d vec; vec << position, velocity; return vec; }
};


// implements a parallelotope tree, i.e. 6d boxes with transforms
struct PNode
{
  PVector extents; // box-like extents
  PVector state; // position and velocity
  PMatrix transform;

  vector<PNode> children;
};