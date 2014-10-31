#include "standardIncludes.h"
typedef Matrix<double, 4, 4> PMatrix;

struct Extents
{
  vector3d position;
  vector3d velocity;
};
struct Contact;

// implements a parallelotope tree, i.e. 6d boxes with transforms. Each node's data is relative to its parent
struct PNode
{
  PNode(const Vector3d &pos, const Quaternion &rot, const Extents &extents, const Vector3d &vel, const Vector3d &angVel) : position(pos), rotation(rot), extents(extents), velocity(vel), angularVelocity(angVel) {}

  Vector3d position;
  Quaternion rotation;
  Extents extents; // box-like extents 
  Vector3d velocity, angularVelocity;

  // converts state data into a single transform matrix, for concatenating
  PMatrix getTransform()
  {
    PMatrix mat;
    mat.identity();
    mat.block(0,3,3,3) = position;
    mat.block(0,0,2,2) = rotation.toMatrix();
    // now boost the matrix here
    return mat;
  }

  vector<PNode> children;
  vector<Contact> contacts;
};

struct Contact
{
  PNode *nodes[2];
  Vector3d normal;
  double depth;
  double depthVelocity;  // together these give a time and speed of impact
};