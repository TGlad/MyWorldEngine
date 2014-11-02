#include "PTree.h"

bool STree::sphericalHyperboloidIntersection(const MinkowskiMatrix &mat1, double radius1, const MinkowskiMatrix &mat2, double radius2)
{
  const MinkowskiMatrix &mat[2] = {mat1, mat2};
  double radius[2] = {radius1, radius2};
  for (int m = 0; m<2; m++) // for each tree
  {

  }
}

void STree::generateContacts(const STree &tree1, const STree &tree2)
{
}

// recursively generate contacts between nodes
void STree::generateContacts()
{
  double depth;
  FourVector normal;
  for (int i = 0; i<children.size(); i++)
  {
    for (int j = i+1; j<children.size(); j++)
    {
      if (sphericalHyperboloidIntersection(children[i].matrix, children[i].radius, children[j].matrix, children[j].radius, depth, normal))
      {
        contacts.push_back(Contact(children[i], children[j], depth, normal));
        // need to generate contacts between these two subtrees now
        generateContacts(children[i], children[j]);
      }
    }
  }
  for (int i = 0; i<children.size(); i++)
    generateContacts();
}

// acts to re-adjust the state to be a good summary of children, for cases where game update of simulation is not completely thorough
void STree::updateState() 
{
  double timeDelta = 1.0/60.0;
  // 1. mass, position and velocity, rotation and angular velocity
  // note that we don't assume position is derivable from vel or visa versa, they aren't always.
  Vector3d delta(0,0,0);
  FourVelocity fourVel(0,0,0,0); 
  mass = 0.0;
  for (int i = 0; i<children.size(); i++)
  {
    children[i].updateState();
    delta += children[i].position() * children[i].mass;
    fourVel += children[i].fourVelocity() * children[i].mass;
    mass += children[i].mass;
  }
  delta /= mass;
  Vector3d velOffset = fourVel.velocity(); // note that this isn't an exact average because SR velocities are non-commutative

  // move so this tree is centred at centre of mass of children
  position() += matrix.transformDelta(delta);
  matrix.lorentz().boost(matrix.transformDelta(velOffset));
  MinkowskiMatrix boost(Event(0,-delta), LorentzMatrix(-velOffset));
  for (int i = 0; i<children.size(); i++)
    children[i].matrix *= boost; // note that boost also changes its position and time value

  // 2. angular velocity
  // we assume angle integrates the angular vel. This will be replaced with an 'absoluteOrientation' algorithm at some point
  // at this point the child positions and velocities are relative to the centre of mass
  momentOfInertia = 0.0;
  Vector3d angVelOffset(0,0,0);
  Matrix3d scatterMatrix;
  scatterMatrix.setZero();
  for (int i = 0; i<children.size(); i++)
  {
    scatterMatrix += children[i].position() * children[i].position().transpose() * children[i].mass;
    angVelOffset += children[i].position().cross(children[i].velocity())*children[i].mass + children[i].angularVelocity*children[i].momentOfInertia;
    momentOfInertia += children[i].position().squaredNorm()*children[i].mass + children[i].momentOfInertia;
  }
  angVelOffset /= momentOfInertia;

  // so rotate our state
  angularVelocity += matrix.transformAngularDelta(angVelOffset);

  // now rotate children back by this amount
  for (int i = 0; i<children.size(); i++)
  {
    children[i].angularVelocity -= angVelOffset;
    MinkowskiMatrix unboost(Event(0), LorentzMatrix(-angVelOffset.cross(children[i].position())));
    children[i].matrix *= unboost;
  }

  // next we need to adjust the bounding radius to match the contents
  radius = 0.0;
  for (int i = 0; i<children.size(); i++)
    radius = max(radius, children[i].position().magnitude() + children[i].radius*children[i].scale);
}

void STree::calculateBounds(Bound &minB, Bound &maxB)
{
  FourVector state = matrix * boxCentre().fourVec;
  FourVector extents = boxExtents().fourVec;
  FourVector vs[4] = {matrix * (FourVector(extents[0], 0, 0, 0)),
                      matrix * (FourVector(0, extents[1], 0, 0)),
                      matrix * (FourVector(0, 0, extents[2], 0)),
                      matrix * (FourVector(0, 0, 0, extents[3]))};
  FourVector parentExtents;
  for (int i = 0; i<4; i++)
    parentExtents[i] = abs(vs[0][i]) + abs(vs[1][i]) + abs(vs[2][i]) + abs(vs[3][i]);
  minB.fourVec = state - parentExtents;
  maxB.fourVec = state + parentExtents;
  /* // how we do this part is uncertain...
  minB.properVelocity = state.velocity - extents.velocity;
  minB.properVelocity = state.velocity - extents.velocity;
  */
}
