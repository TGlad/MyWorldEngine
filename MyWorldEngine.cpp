#include "standardIncludes.h"
#include "PTree.h"

int _tmain(int argc, _TCHAR* argv[])
{
  // 1. build a tree
  PNode node(Vector3d(0,0,0), Vector3d(1,1,1), PIdentity());
  PNode child[4];
  for (int i = 0; i<4; i++)
  {
    child[i] = node;
    child[i].state.position = Vector3d(i,0,0);
    child[i].transform *= 0.5;
    node.children.push_back(child[i]);
  }

	return 0;
}

