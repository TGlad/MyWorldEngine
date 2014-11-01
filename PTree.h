#include "standardIncludes.h"
#include "Vectors.h"
#include "Matrices.h"

struct Contact;


// implements a parallelotope tree, i.e. 6d boxes with transforms. Each node's data is relative to its parent
struct PTree
{
  PTree(const Vector3d &pos, const Quaternion &rot, const Extents &extents, const Vector3d &vel, const Vector3d &angVel, double scale) : position(pos), rotation(rot), extents(extents), velocity(vel), angularVelocity(angVel), scale(scale), mass(1.0), momentOfInertia(1.0) {}

  // relative to parent. These are summary data of the contents of this tree
  Vector3d position;
  Vector3d velocity; // since velocities don't sum in special relativity, should we use a rapidity vector here instead? No because they're only additive when parallel
  double scale;
  Quaternion rotation;
  Vector3d angularVelocity;
  double mass;
  double momentOfInertia; // will be replaced with an inertial tensor at some point

  void generateContacts();
  void updateState();

  // relative to scale and rotation above
  StateVector minBound; 
  StateVector maxBound; 
  StateVector boxCentre(){ return (minBound + maxBound) / 2.0; } 
  StateVector boxExtents(){ return (maxBound - minBound) / 2.0; } 
  void calculateBounds(StateVector &minB, StateVector &maxB);

  // converts state data into a single transform matrix, for concatenating
  PMatrix getTransform();
  StateVector transformState(const StateVector &localState);
  Vector3d transformPosition(const Vector3d &localPosition);
  Vector3d transformDelta(const Vector3d &localOffset);
  Vector3d transformAngularDelta(const Vector3d &localOffset);

  vector<PTree> children;
  vector<Contact> contacts;
};

struct Contact
{
  PTree *trees[2];
  Vector3d normal;
  double depth;
  double depthVelocity;  // together these give a time and speed of impact
};