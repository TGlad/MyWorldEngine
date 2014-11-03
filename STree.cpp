#include "PTree.h"
static double timeDelta = 1.0/60.0;

bool STree::sphericalHyperboloidIntersection(const STree &child1, const STree &child2, Event &contactEvent, Vector3d &normal)
{
  // for now we calculate these each call. 
  double radius[2], expansionRate[2];
  const STree &child = {child1, child2};
  for (int i = 0; i<2; i++)
  {
    // get linearised expansion
    child[i].calculateRadiusAndExpansionRate(child[i].time(), radius[i], expansionRate[i]);
    radius[i] -= expansionRate[i]*child[i].time; // in parent's time
  }
  // formula for radius at time t. But we linearise it
  // double r = sqrt(sqr((t - minRadiusTime)*expansionVelocity) + sqr(minRadius));

  Vector3d relVel = child[1].velocity() - child[0].velocity();
  Vector3d relPos = (child[1].position() - child[1].time*child[1].velocity()) - (child[0].position() - child[0].time*child[0].velocity()); // we shift the position to get into the parent's time frame
  // Solve: at^2 + bt + c = 0
  // coefficients of t^2
  double a = relVel.squaredNorm() - sqr(expansionRate[0]) - sqr(expansionRate[1]);
  // coefficients of t
  double b = 2.0*(relVel.dot(relPos) - expansionRate[0]*radius[0] - expansionRate[1]*radius[1]);
  // constant coefficients
  double c = relPos.squaredNorm() - sqr(radius[0]) - sqr(radius[1]);

  double root = sqrt(b*b - 4*a*c);
  double t[2] = {(-b + root) / (2.0*a), (-b - root) / (2.0*a)};
  for (int i = 0; i<2; i++)
  {
    Vector3d pos0 = child[0].position + child[0].velocity*(t[i] - child[0].time);
    Vector3d pos1 = child[1].position + child[1].velocity*(t[i] - child[1].time);
    normal = pos1 - pos0;
    normal.normalize();
    contactEvent.position = pos0 + normal * child[0].radius;
    contactEvent.time = contactEvent.position.magnitude();
    if (t[i] > contactEvent.time && t[i] < contactEvent.time + timeDelta)
      return true; // if we want, we could check the proximity, and repeat the procedure to get a better time when expansion is not so linear
  }
  return false;
}

void STree::generateContacts(const STree &tree1, const STree &tree2)
{
}

// recursively generate contacts between nodes
void STree::generateContacts()
{
  Event contactEvent;
  Vector3d normal;
  for (int i = 0; i<children.size(); i++)
  {
    for (int j = i+1; j<children.size(); j++)
    {
      if (sphericalHyperboloidIntersection(children[i], children[j], contactEvent, normal))
      {
        contacts.push_back(Contact(children[i], children[j], contactEvent, normal));
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

