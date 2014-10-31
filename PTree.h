#include "standardIncludes.h"
typedef Matrix<double, 4, 4> PMatrix;

struct Extents
{
  Extents(const Vector3d &pos, const Vector3d &vel) : position(pos), velocity(vel) {}
  vector3d position;
  vector3d velocity;
};
struct Contact;

// implements a parallelotope tree, i.e. 6d boxes with transforms. Each node's data is relative to its parent
struct PTree
{
  PTree(const Vector3d &pos, const Quaternion &rot, const Extents &extents, const Vector3d &vel, const Vector3d &angVel, double scale) : position(pos), rotation(rot), extents(extents), velocity(vel), angularVelocity(angVel), scale(scale), mass(1.0), momentOfInertia(1.0) {}

  // relative to parent. These are summary data of the contents of this tree
  Vector3d position;
  Vector3d velocity;
  double scale;
  Quaternion rotation;
  Vector3d angularVelocity;
  double mass;
  double momentOfInertia; // will be replaced with an inertial tensor at some point

  // relative to scale and rotation above
  Extents extents; 

  // converts state data into a single transform matrix, for concatenating
  PMatrix getTransform()
  {
    PMatrix mat;
    mat.identity();
    mat.block(0,3,3,3) = position;
    mat.block(0,0,2,2) = rotation.toMatrix();
    // now boost the matrix here
    mat *= scale;
    return mat;
  }
  Vector3d transformPosition(const Vector3d &localPosition)
  {
    Vector3d trans = (rotation * localPosition)*scale;
    // now boost the vector
    return trans + position;
  }
  Vector3d transformDelta(const Vector3d &localOffset)
  {
    Vector3d trans = (rotation * localOffset)*scale;
    // now boost the vector
    return trans;
  }
  Vector3d transformAngularDelta(const Vector3d &localOffset)
  {
    Vector3d trans = (rotation * localOffset);  // scale doesn't affect angular measurements
    // now boost the vector. What does a boost do to an angular measurement?
    return trans;
  }


  generateContacts()
  {

  }
  // acts to re-adjust the state to be a good summary of children, for cases where game update of simulation is not completely thorough
  updateState() 
  {
    double timeDelta = 1.0/60.0;
    // 1. mass, position and velocity, rotation and angular velocity
    // note that we don't assume position is derivable from vel or visa versa, they aren't always.
    Vector3d delta(0,0,0), velOffset(0,0,0); 
    mass = 0.0;
    vector<Vector3d> offsets(children.size());
    for (int i = 0; i<children.size(); i++)
    {
      offsets[i] = -children[i].position; // we should record old orientations too, but this is extra cost
      children[i].updateState();
      delta += children[i].position * children[i].mass;
      velOffset += children[i].velocity * children[i].mass;
      mass += children[i].mass;
      offsets[i] += children[i].position;
    }
    delta /= mass;
    velOffset /= mass; 
    angleDelta /= momentOfInertia;

    // move so this tree is centred at centre of mass of children
    position += transformDelta(delta);
    velocity += transformDelta(velOffset);
    for (int i = 0; i<children.size(); i++)
    {
      children[i].position -= delta;
      children[i].velocity -= velOffset;
    }

    // 2. angular velocity
    // we assume angle integrates the angular vel. This will be replaced with an 'absoluteOrientation' algorithm at some point
    // at this point the child positions and velocities are relative to the centre of mass
    double rotationTotalMoment = 0.0;
    momentOfInertia = 0.0;
    Vector3d angVelOffset(0,0,0), angleDelta(0,0,0);
    for (int i = 0; i<children.size(); i++)
    {
      double childMoment = children[i].position.squaredNorm()*children[i].mass;
      angleDelta += children[i].position.cross(offsets[i])*children[i].mass; // plus rotational deltas, but these are excess computation
      angVelOffset += children[i].position.cross(children[i].velocity)*children[i].mass + children[i].angularVelocity*children[i].momentOfInertia;
      rotationTotalMoment += childMoment;
      momentOfInertia += childMoment + children[i].momentOfInertia;
    }
    angleDelta /= rotationTotalMoment;
    angVelOffset /= momentOfInertia;

    // so rotate our state
    angularVelocity += transformAngularDelta(angVelOffset);
    Quaternion rotationDelta(angleDelta);
    rotation = rotation * rotationDelta; // apply rotation over the top of our local angle delta

    // now rotate children back by this amount
    for (int i = 0; i<children.size(); i++)
    {
      children[i].velocity -= angVelOffset.cross(children[i].position);
      children[i].angularVelocity -= angVelOffset;
      children[i].rotation *= ~rotationDelta;
      children[i].position = ~rotationDelta.rotateVector(children[i].position);
    }
  }

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