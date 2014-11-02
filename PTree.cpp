#include "PTree.h"

// note that mat1 and mat2 are box centres and extent1, extent2 are extents from this
// TODO: decide whether space should be R4 or M3+1 for calculations. 
bool ParallelochoronIntersection(const MinkowskiMatrix &mat1, const Bound &extent1, const MinkowskiMatrix &mat2, const Bound &extent2, double &maxDepth, FourVector &maxDepthNormal)
{
  // we could convert one into a rectangle, but this involves matrix inverse, so easier to work with both in original format

  maxDepth = -1e10;
  const MinkowskiMatrix &mat[2] = {mat1, mat2};
  const Bound &extent[2] = {extent1, extent2}; 
  LorentzMatrix axes[2] = {mat[0].lorentz() * diagonal(extent[0].fourVec), mat[1].lorentz() * diagonal(extent[1].fourVec)};
  for (int m = 0; m<2; m++)
  {
    FourVector toOther = mat[1].event() - mat[0].event();
    // 1. corners versus volumes
    for (int c = 0; c<16; c++) // count the corners
    {
      FourVector corner(0,0,0,0);
      for (int j = 0; j<4; j++)
        corner += axes[m].col(j) * (2.0*(double)((c>>j)&1) - 1.0);
//      if (corner.dot(toOther) < 0.0)
//        continue; // only need to do half the corners (the ones in direction towards other)
      corner += mat[m].lorentz().fourVector();

      // volumes in 4d have a normal, in this case we are using a cartesian space to intersect... not sure but maybe it is possible in Minkowski space
      for (int v = 0; v<4; v++) // count the volumes
      {
        FourVector volumeCentre = axes[1-m].col(v);
        if (volumeCentre.dot(toOther) < 0.0)
          volumeCentre = -volumeCentre;
        FourVector volumeNormal = volumeCentre;
        volumeNormal.normalize();
        volumeCentre += mat[1-m].lorentz().fourVector();

        double depth = (volumeCentre - corner).dot(volumeNormal);
        if (depth > maxDepth)
        {
          maxDepth = depth;
          maxDepthNormal = volumeNormal;
        }
      }
    }

    // edges vs faces
    for (int d = 0; d<4; d++) // count the directions
    {
      FourVector edgeDirection = mat[m].lorentz().col(d);

      for (int f = 0; f<6; f++) // count the face directions
      {
        // axes making up face
        int d1[6] = {0,0,0,1,1,2};
        int d2[6] = {1,2,3,2,3,3};
        // other axes which provide face position
        int d3[6] = {2,1,1,3,2,0}; 
        int d4[6] = {3,3,2,0,0,1};
        FourVector faceDirection[2] = {axes[1-m].col(d1[f]), axes[1-m].col(d1[f])};

        // now get a normal
        FourVector edgeFaceNormal(1,2,3,4);
        edgeFaceNormal -= edgeDirection * edgeFaceNormal.dot(edgeDirection) / edgeDirection.squaredNorm();
        edgeFaceNormal -= faceDirection[0] * edgeFaceNormal.dot(faceDirection[0]) / edgeDirection.squaredNorm();
        edgeFaceNormal -= faceDirection[1] * edgeFaceNormal.dot(faceDirection[1]) / edgeDirection.squaredNorm();
        edgeFaceNormal.normalize();

        for (int e = 0; e<8; e++) // count the edge centres
        {
          FourVector edgeCentre = mat[m].lorentz().fourVector();
          for (int j = 0; j<3; j++)
            edgeCentre += axes[m].col((d+j+1)%4) * (2.0*(double)((e>>j)&1) - 1.0);
        
          FourVector faceCentre(0,0,0,0);
          for (int c = 0; c<4; c++)
          {
            double scale1[4] = {1,1,-1,-1};
            double scale2[4] = {1,-1,1,-1};
            faceCentre += axes[1-m].col(d3[f])*scale1[c];
            faceCentre += axes[1-m].col(d4[f])*scale2[c];
            // can we cull out some face centres here?

            FourVector faceNormal = edgeFaceNormal.dot(faceCentre) > 0.0 ? edgeFaceNormal : -edgeFaceNormal;
            faceCentre += mat[1-m].lorentz().fourVector();

            double depth = (faceCentre - edgeCentre).dot(faceNormal);
            if (depth > maxDepth)
            {
              maxDepth = depth;
              maxDepthNormal = faceNormal;
            }
          }
        }
      }
    }
  }
  return maxDepth > 0.0;
}

void PTree::generateContacts(const PTree &tree1, const PTree &tree2)
{
}

// recursively generate contacts between nodes
void PTree::generateContacts()
{
  vector<MinkowskiMatrix> centredMats(children.size());
  for (int i = 0; i<children.size(); i++)
  {
    centredMats[i] = children[i].matrix + children[i].matrix.lorentz() * children[i].boxCentre();
    centredMats[i] = 
  }
  double depth;
  FourVector normal;
  for (int i = 0; i<children.size(); i++)
  {
    for (int j = i+1; j<children.size(); j++)
    {
      if (ParallelochoronIntersection(centredMats[i], children[i].boxExtents(), centredMats[j], children[j].boxExtents(), depth, normal))
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
void PTree::updateState() 
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

  // now calculate approximate new biggest eigenvector from scatterMatrix
  Vector3d currentBiggestEigen = matrix.transformDelta(Vector3d(1,0,0)); // x axis as biggest
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
  angularVelocity += matrix.transformAngularDelta(angVelOffset);
  matrix *= MinkowskiMatrix(Event(0), LorentzMatrix(rotationDelta));
  MinkowskiMatrix unrotate(Event(0), LorentzMatrix(~rotationDelta));

  // now rotate children back by this amount
  for (int i = 0; i<children.size(); i++)
  {
    children[i].matrix *= unrotate;
    children[i].angularVelocity -= angVelOffset;
    MinkowskiMatrix unboost(Event(0), LorentzMatrix(~rotationDelta, -angVelOffset.cross(children[i].position())));
    children[i].matrix *= unboost;
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

void PTree::calculateBounds(Bound &minB, Bound &maxB)
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
