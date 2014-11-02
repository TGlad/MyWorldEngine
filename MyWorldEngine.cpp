#include "standardIncludes.h"
#include "PTree.h"

int _tmain(int argc, _TCHAR* argv[])
{
  // 1. build a tree
  STree tree(Vector3d(0,0,0), Quaternion(1,0,0,0), Extents(1,1,1, 1,1,1), Vector3d(0,0,0), Vector3d(0,0,0), 1.0);
  STree child[4];
  for (int i = 0; i<4; i++)
  {
    child[i] = tree;
    child[i].position = Vector3d(i,0,0);
    child[i].scale *= 0.5;
    tree.children.push_back(child[i]);
  }

	return 0;
}

