#include "PTree.h"

void PTree::calculateBounds(StateVector &minB, StateVector &maxB)
{
  StateVector state = transformState(boxCentre());
  StateVector extents = boxExtents();
  Vector3d vs[3] = {transformDelta(Vector3d(extents.position[0], 0, 0)),
                    transformDelta(Vector3d(0, extents.position[1], 0)),
                    transformDelta(Vector3d(0, 0, extents.position[2]))};
  Vector3d parentExtents;
  for (int i = 0; i<3; i++)
    parentExtents[i] = abs(vs[0][i]) + abs(vs[1][i]) + abs(vs[2][i]);
  minB.position = state.position - parentExtents;
  maxB.position = state.position + parentExtents;
  minB.velocity = state.velocity - extents.velocity;
  minB.velocity = state.velocity - extents.velocity;
}

// converts state data into a single transform matrix, for concatenating
PMatrix PTree::getTransform()
{
  PMatrix mat;
  mat.identity();
  mat.block(0,3,3,3) = position;
  mat.block(0,0,2,2) = rotation.toMatrix();
  // now boost the matrix here
  mat *= scale;
  return mat;
}
StateVector PTree::transformState(const StateVector &localState)
{
  // how do velocities transform? I guess we really want to apply special relativity velocity addition
  return StateVector(transformPosition(localState.position), localState.velocity + velocity);
}
Vector3d PTree::transformPosition(const Vector3d &localPosition)
{
  Vector3d trans = (rotation * localPosition)*scale;
  // now boost the vector
  return trans + position;
}
Vector3d PTree::transformDelta(const Vector3d &localOffset)
{
  Vector3d trans = (rotation * localOffset)*scale;
  // now boost the vector
  return trans;
}
Vector3d PTree::transformAngularDelta(const Vector3d &localOffset)
{
  Vector3d trans = (rotation * localOffset);  // scale doesn't affect angular measurements
  // now boost the vector. What does a boost do to an angular measurement?
  return trans;
}


void PTree::generateContacts()
{

}
// acts to re-adjust the state to be a good summary of children, for cases where game update of simulation is not completely thorough
void PTree::updateState() 
{
  double timeDelta = 1.0/60.0;
  // 1. mass, position and velocity, rotation and angular velocity
  // note that we don't assume position is derivable from vel or visa versa, they aren't always.
  Vector3d delta(0,0,0), velOffset(0,0,0); 
  mass = 0.0;
  for (int i = 0; i<children.size(); i++)
  {
    children[i].updateState();
    delta += children[i].position * children[i].mass;
    velOffset += children[i].velocity * children[i].mass;
    mass += children[i].mass;
  }
  delta /= mass;
  velOffset /= mass; 

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
  Vector3d angVelOffset(0,0,0);
  Matrix3d scatterMatrix;
  scatterMatrix.setZero();
  for (int i = 0; i<children.size(); i++)
  {
    scatterMatrix += children[i].position * children[i].position.transpose() * children[i].mass;
    angVelOffset += children[i].position.cross(children[i].velocity)*children[i].mass + children[i].angularVelocity*children[i].momentOfInertia;
    momentOfInertia += children[i].position.squaredNorm()*children[i].mass + children[i].momentOfInertia;
  }
  angVelOffset /= momentOfInertia;

  // now calculate approximate new biggest eigenvector from scatterMatrix
  Vector3d currentBiggestEigen = rotation * Vector3d(1,0,0); // x axis as biggest
  Vector3d newBiggestEigen = scatterMatrix * currentBiggestEigen;
  
  // we only rotate this node towards the biggest eigenvalued eigenvector. Good for approximating long shapes, but
  // we would need a full eigendecomposition to properly diagonalise an inertial tensor from the children. 
  // this is a bit expensive at the moment, especially as proper physics simulation is not yet implemented.

  // non-unit length rotation from cross and dot of vectors
  Vector3d cross = currentBiggestEigen.cross(newBiggestEigen);
  Quaternion rotationDelta(currentBiggestEigen.dot(newBiggestEigen), cross[0], cross[1], cross[2]);
  rotationDelta[0] += rotationDelta.norm();
  rotationDelta.normalize();

  // so rotate our state
  angularVelocity += transformAngularDelta(angVelOffset);
  rotation = rotation * rotationDelta; // apply rotation over the top of our local angle delta

  // now rotate children back by this amount
  for (int i = 0; i<children.size(); i++)
  {
    children[i].velocity -= angVelOffset.cross(children[i].position);
    children[i].angularVelocity -= angVelOffset;
    children[i].rotation *= ~rotationDelta;
    children[i].position = ~rotationDelta.rotateVector(children[i].position);
  }

  // next we need to adjust the bounding extents to match the contents
  Extents minExtents(Vector3d(-1e10, -1e10, -1e10), Vector3d(-1e10, -1e10, -1e10));
  Extents maxExtents(Vector3d(1e10, 1e10, 1e10), Vector3d(1e10, 1e10, 1e10));

  for (int i = 0; i<children.size(); i++)
  {
    StateVector childMin, childMax;
    children[i].calculateBounds(childMin, childMax);
    minExtents.position = minVector(minExtents.position, childMin.position);
    maxExtents.position = maxVector(maxExtents.position, childMax.position);
    minExtents.velocity = minVector(minExtents.velocity, childMin.velocity);
    maxExtents.velocity = maxVector(maxExtents.velocity, childMax.velocity);
  }
}

vector<PTree> children;
vector<Contact> contacts;
};