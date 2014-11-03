#include "standardIncludes.h"
#include "Vectors.h"
#include "Matrices.h"

struct Contact;

// implements a spherical hyperboloid tree, more precisely this is a tree where each element bounds its children by a boosted 4d one-sheeted hyperboloid.
// The bound remains spherical in 3d (the length contraction is ignored in the bound, even though length contraction occurs on the children). 
// The hyperboloid represents linearised angular rotation of the tree, which has the property that angular rotations compose with translations, 
// allowing relative motion to be dealt with exactly between sub-trees. The special case of a cone represents expansion or contraction of the tree. 
// Our world objects are therefore capable of Minkowski transforms, linearised angular velocities plus dilation/contraction; but not squashed.
// A light cone is another special case of this shape, and so is a light-hyperboloid (set of events of equally far away).
struct STree
{
  STree(const Vector3d &pos, const Quaternion &rot, const Extents &extents, const Vector3d &vel, const Vector3d &angVel, double scale) : position(pos), rotation(rot), extents(extents), velocity(vel), angularVelocity(angVel), scale(scale), mass(1.0), momentOfInertia(1.0) {}

  // relative to parent. These are summary data of the contents of this tree
  MinkowskiMatrix matrix;

  // accessors from matrix. TODO: make these efficient!
  Vector3d position(){ return matrix.event().position(); }
  Vector3d velocity(){ return matrix.lorentz().fourVector().velocity(); } 
  FourVelocity fourVelocity(){ return (FourVelocity)matrix.lorentz().fourVector(); } 
  double scale(){ return matrix.lorentz().scale() }; // TODO: fix. Can we have a just spatial scale? or must we scale time too?
  Quaternion calculateRotation(){ return Quaternion(matrix.lorentz().spatial()); }
  
  // Note that we require this in its boosted, parent frame, not in its local frame
  calculateRadiusAndExpansionRate(double &radius, double &expansionRate){ /* TODO: calculate these */ } 
  Vector3d angularVelocity;
  double mass;
  double momentOfInertia; // will be replaced with an inertial tensor at some point

  void generateContacts();
  void updateState();

  // relative to scale above
  double radius; 

  vector<STree> children;
  vector<Contact> contacts;
protected:
  bool sphericalHyperboloidIntersection(const STree &child1, const STree &child2, Event &contactEvent, Vector3d &normal);
};

struct Contact
{
  Contact(const STree &p1, const STree &p2, Event &contactEvent, FourVector &normal)
  {
    trees[0] = &p1;
    trees[1] = &p2;
    event = contactEvent;
    this->normal = normal;
  }
  STree *trees[2];
  Vector3d normal;
  Event event;
};