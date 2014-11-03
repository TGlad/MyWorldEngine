// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include<Eigen>
using namespace Eigen;
using namespace std;

inline minVector(const Vector3d &u, const Vector3d &v)
{
  return Vector3d(min(u[0], v[0]), min(u[1], v[1]), min(u[2], v[2]));
}

inline maxVector(const Vector3d &u, const Vector3d &v)
{
  return Vector3d(max(u[0], v[0]), max(u[1], v[1]), max(u[2], v[2]));
}
template<class T>
inline T sqr(const T &a)
{
  return a*a;
}


// TODO: reference additional headers your program requires here
