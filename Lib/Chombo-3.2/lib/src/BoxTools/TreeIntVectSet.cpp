#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TreeIntVectSet.H"
#include "Vector.H"
#include "Arena.H"
#include "RefCountedPtr.H"
#include "DisjointBoxLayout.H"
#include "ProblemDomain.H"

#include "SPMD.H"
#include "Tuple.H"
#include "NamespaceHeader.H"

using std::ostream;
using std::cout;
using std::endl;

// The implementation of the tree here is a little
// weird for some.  It is a strictly 2-ary tree.
// the terminal nodes of the tree represent a 2-ary
// box of information, not just a point.  it uses
// the bit pattern in the pointer data holder to
// store the true/false information for the points in the
// 2-ary box. This requires 8 bits for the 3D case, a bit
// for each point in the box.  hence the appearance of
// the magic '256' in the code (2^8).  TreeNode pointers
// that have a value less than '256' are not pointing at
// memory, but holding the leaf information.

//
//  Further details of how non-recursive TreeNode design works:
//
//  (for a 2D tree)
//
//    (m_tree)
//           + -- 0
//
//             (a)+ - 0
//                    1
//                1   +
//                    1  <------you are here
//
//                + - + - 0
//                    0   1
//                    0   1
//                    0   0
//
//  for the node indicated, the 'index' vector  would contain
//
//  index=[ 0      1    3  ...............]
//parents=[&m_tree &a   ..................]
//
// or directly refered to as m_tree.nodes[1].nodes[3]
//
//  the tree indicates a covering of an index space in either
//  1 or 0.   1 is stored in the tree by pointing at the static
//  'full' data member, 0 is stored as a 0.
//
//  every 'nodes' member of a TreeNode object can be either
//
//   0, &full, or a pointer .
//
//   The interpretation of the tree depends on what m_spanBox is.
//     nodes[i] indicates whether the i'th quadrant of the parent
//     Box is full, empty, or needs to be parsed deeper.
//
//

// NOTE:  Had to replace the static const integer nodeSize with
// a macro named TIVS_NODESIZE to compile with xlC on AIX.
// Currently, this only affects this file and the .H  (ndk)

// handy macro for tree parsing operations....

#define NEXT(index_local, depth)                      \
      index_local[0] = 0;                             \
      index_local[depth]++;                           \
          while (index_local[depth] == TIVS_NODESIZE)  \
            {                                         \
              index_local[depth] = 0;                 \
              depth--;                                \
              index_local[depth]++;                   \
            }

void TreeIntVectSet::define(const Box& a_box)
{

  clearTree(m_tree);
  m_minBox = a_box;
  m_depth = 1;
  if (a_box.isEmpty())
    {
      m_spanBox = Box();
      return;
    }

  //    int size = 0;
  //    for (int i=0; i<SpaceDim; ++i) size = Max(size, a_box.size(i));
  //    int clippedSize = (size/2)*2;
  //    if (clippedSize < size) size = clippedSize+2;
  //
  //    int minBoxSize = 2;
  //    m_depth = 1;
  //    while (minBoxSize < size)
  //      {
  //        minBoxSize*=2;
  //        m_depth++;
  //      }
  //
  //    int m = index.size();
  //    m = Max(m_depth+1, m);
  //    if (m > index.size())
  //      {
  //        index.resize(m);
  //        parents.resize(m);
  //        boxes.resize(m);
  //      }
  //
  //    m_spanBox = Box(a_box.smallEnd(),
  //                    a_box.smallEnd()
  //                    +(minBoxSize-1)*IntVect::Unit);

  *this |= a_box;
}

void TreeIntVectSet::clear()
{
  clearTree(m_tree);
  m_spanBox = Box();
  m_minBox  = Box();
  m_depth = 1;
}

TreeIntVectSet& TreeIntVectSet::operator=(const TreeIntVectSet& rhs)
{
  define(rhs);
  return *this;
}

void TreeIntVectSet::define(const TreeIntVectSet& a_tivs)
{
  if (this == &a_tivs) return;

  clearTree(m_tree);

  m_spanBox = a_tivs.m_spanBox;
  m_minBox  = a_tivs.m_minBox;
  m_depth  = a_tivs.m_depth;

  if (a_tivs.m_tree.nodes == 0 || a_tivs.m_tree.nodes == &full)
    {
      m_tree = a_tivs.m_tree;
      return;
    }

  expandNode(m_tree);

  // perform direct Tree-to-Tree copy operation.
  // this is a little more involved than just
  // adding the boxes from one to the other, but
  // it should be much more efficient and faster

  cloneNode(a_tivs.m_tree, m_tree);

}

void TreeIntVectSet::cloneNode(const TreeNode& src, TreeNode& dest)
{
  static Vector<const TreeNode*> otherParents;
  if (parents.size() > otherParents.size())
    otherParents.resize(parents.size());

  if (src.nodes ==0 || src.nodes == &full)
    {
      dest=src;
      return;
    }
  if (dest.nodes == 0 || dest.nodes == &full)
    {
      expandNode(dest);
    }

  otherParents[0] = &(src);
  parents[0]  = &(dest);
  index[1] = 0;
  int depth = 1;
  while (depth!=0)
    {
      const TreeNode* otherParent  = otherParents[depth-1];
      TreeNode* thisParent   = (TreeNode*)parents[depth-1];
      int ind = index[depth];
      const TreeNode& otherCurrent = otherParent->nodes[ind];
      TreeNode& thisCurrent  = thisParent->nodes[ind];

      if (otherCurrent.nodes == 0 || otherCurrent.nodes == &full) // terminal
        {
          thisCurrent = otherCurrent;
          nextNode(depth);
        }
      else  // not terminal node, add nodes and traverse deeper
        {
          expandNode(thisCurrent);
          otherParents[depth] = &(otherCurrent);
          parents[depth]  = &(thisCurrent);
          depth++;
          index[depth] = 0;
        }
    }
}

void TreeIntVectSet::clearTree(TreeNode& tree)
{
  // quick check to eliminate the work of tiny tree deletion
  if (tree.nodes == 0 || tree.nodes == &full)
  {
    tree.nodes = 0;
    return;
  }

  static Vector<TreeNode*> parents_local;
  static Vector<int>       index_local;
  parents_local.resize(TreeIntVectSet::parents.size());
  index_local.resize(TreeIntVectSet::index.size());
  // Line below added by petermc, 11 Dec 2001, to avoid crash
  if (index_local.size() == 0)
  {
    tree.nodes = 0;
    return;
  }

  index_local[0] = 0;

  int depth = 1;

  parents_local[0] = &(tree);
  index_local[1] = 0;
  while (depth != 0)
    {
      TreeNode* parent =  parents_local[depth-1];
      int ind = index_local[depth];
      //Vector<int>& indexRef = index;
      //Vector<TreeNode*>& parentRef = parents;
      TreeNode& current = parent->nodes[ind];
      if (current.nodes == 0 || current.nodes == &full)
        {
          index_local[depth]++;
          while (index_local[depth] == TIVS_NODESIZE)
            {
              index_local[depth] = 0;
              //delete[] parents[depth-1]->nodes;
              treeNodePool->returnPtr(parents_local[depth-1]->nodes);
              depth--;
              index_local[depth]++;
            }
        }
      else
        {
          parents_local[depth] = &(current);
          depth++;
          index_local[depth] = 0;
        }
    }
  tree.nodes = 0;

}

void TreeIntVectSet::refine(int iref)
{
  if (iref == 1) return;
  CH_assert(iref >= 1);
  int refinements = 1;
  int r = 2;
  while (r < iref)
     {
       refinements++;
       r*=2;
     }

  CH_assert(r == iref); // check iref for power of 2
  m_spanBox.refine(iref);
  m_minBox.refine(iref);
  m_depth+=refinements;
  int m=index.size();
  m = Max(m_depth+1, m);
  if (m > index.size())
    {
      index.resize(m);
      parents.resize(m);
      boxes.resize(m);
          bufferOffset.resize(m);
    }
}

void TreeIntVectSet::shift(const IntVect& iv)
{
  static Vector<Box> boxs;
  int size;

  createBoxes(boxs, size);
  clear();
  for (int i=0; i < size; ++i)
    {
      boxs[i].shift(iv);
    }
  for (int i=0; i<size; ++i)
    {
      *this |= boxs[i];
    }
}

void TreeIntVectSet::grow(int igrow)
{
  if (igrow == 0 ) return;
  if (igrow < 0) MayDay::Error("TreeIntVectSet::grow(int) called with negative value");

  static Vector<Box> boxs;
  int size;

  createBoxes(boxs, size);
  clear();
  for (int i=0; i < size; ++i)
    {
      boxs[i].grow(igrow);
    }
  for (int i=0; i<size; ++i)
    {
      *this |= boxs[i];
    }
}

void TreeIntVectSet::grow(int idir, int igrow)
{
  if (igrow == 0 ) return;
  if (igrow < 0) MayDay::Error("TreeIntVectSet::grow(int) called with negative value");
  static Vector<Box> boxs;
  int size;

  createBoxes(boxs, size);
  clear();
  for (int i=0; i<size; ++i)
    {
      boxs[i].grow(idir, igrow);
    }
  for (int i=0; i<size; ++i)
    {
      *this |= boxs[i];
    }
}

void TreeIntVectSet::growHi()
{
  static Vector<Box> boxs;
  int size;

  createBoxes(boxs, size);
  clear();
  for (int i=0; i < size; ++i)
    {
      boxs[i].setBig(boxs[i].bigEnd() + 1);
      *this |= boxs[i];
    }
}

void TreeIntVectSet::growHi(const int a_dir)
{
  static Vector<Box> boxs;
  int size;

  createBoxes(boxs, size);
  clear();
  for (int i=0; i < size; ++i)
    {
      boxs[i].growHi(a_dir);
      *this |= boxs[i];
    }
}

/*
TreeIntVectSet TreeIntVectSet::chop(int dir, int chop_pnt)
{
  if (m_minBox.smallEnd(dir) >= chop_pnt)
    {
      TreeIntVectSet rtn(*this);
      clear();
      return rtn;
    }
  if (m_minBox.bigEnd(dir) < chop_pnt)
    {
      return TreeIntVectSet();
    }

  static Vector<Box> boxs;
  int size;
  createBoxes(boxs, size);

  clear();
  TreeIntVectSet rtn;
  for (int i=0; i<size; ++i)
    {
      Box& box = boxs[i];
      if (box.smallEnd(dir) >= chop_pnt)
        rtn |= box;
      else if (box.bigEnd(dir) < chop_pnt)
        *this |= box;
      else
        {
          Box hi = box.chop(dir, chop_pnt);
          rtn |= hi;
          *this |= box;
        }
    }
  return rtn;
}

*/

TreeIntVectSet TreeIntVectSet::chop(int dir, int chop_pnt)
{
  if (m_minBox.smallEnd(dir) >= chop_pnt)
    {
      TreeIntVectSet rtn(*this);
      clear();
      return rtn;
    }
  if (m_minBox.bigEnd(dir) < chop_pnt)
    {
      return TreeIntVectSet();
    }

  Box min  = m_minBox;
  Box chop = m_spanBox;

  TreeIntVectSet rtn(*this);

  Box chopThis = chop.chop(dir, chop_pnt);
  rtn   -= chop;
  *this -= chopThis;

//  rtn.recalcMinBox();
//  recalcMinBox();

  m_minBox = min;
  rtn.m_minBox = m_minBox.chop(dir, chop_pnt);

  return rtn;

}

void TreeIntVectSet::chop(int a_dir, int a_chop_pnt, TreeIntVectSet& a_hi)
{
  if (m_minBox.smallEnd(a_dir) >= a_chop_pnt)
    {
      this->swap(a_hi);
      return;
    }
  if (m_minBox.bigEnd(a_dir) < a_chop_pnt)
    {
      a_hi.clear();
      return;
    }
  a_hi.clear();

  Box min  = m_minBox;
  Box chop = m_spanBox;

  Box hiBox = chop.chop(a_dir, a_chop_pnt);
  remove(hiBox, &a_hi);

}

// void TreeIntVectSet::chop(int a_dir, int a_chop_pnt, TreeIntVectSet& a_hi)
// {
//   if (m_minBox.smallEnd(a_dir) >= a_chop_pnt)
//     {
//       this->swap(a_hi);
//       return;
//     }
//   if (m_minBox.bigEnd(a_dir) < a_chop_pnt)
//     {
//       a_hi.clear();
//       return;
//     }
//   a_hi.clear();
//   Vector<Box> vbox;
//   int size;
//   createBoxes(vbox, size);
//   clear();
//   for (int i=0; i<vbox.size(); i++)
//     {
//       Box& b = vbox[i];
//       if (b.bigEnd(a_dir) < a_chop_pnt)
//      {
//        this->operator|=(b);
//      }
//       else if (b.smallEnd(a_dir) >= a_chop_pnt)
//      {
//        a_hi |= b;
//      }
//       else
//      {
//        Box hi = b.chop(a_dir, a_chop_pnt);
//        a_hi |= hi;
//        this->operator|=(b);
//      }
//     }
// }

void TreeIntVectSet::swap(TreeIntVectSet& a_other)
{
  TreeNode tmp1 = a_other.m_tree;
  a_other.m_tree = m_tree;
  m_tree = tmp1;

  int tmp2 = a_other.m_depth;
  a_other.m_depth = m_depth;
  m_depth = tmp2;

  Box tmp3 = a_other.m_minBox;
  a_other.m_minBox = m_minBox;
  m_minBox = tmp3;

  Box tmp4 = m_spanBox;
  a_other.m_spanBox = m_spanBox;
  m_spanBox = tmp4;
}

/*
  void TreeIntVectSet::coarsen(int icoarse)
  {
    static  Vector<Box> boxs;
   int size;

    createBoxes(boxs, size);
    clear();
    for (int i=0; i<size; ++i)
      {
        Box& box = boxs[i];
        box.coarsen(icoarse);
        this->operator|=(box);
      }
    compact();
  }
*/

void TreeIntVectSet::coarsen(int icoarse)
{
  if (icoarse == 1) return;
  CH_assert(icoarse >= 1);
  int coarsenings = 1;
  int c = 2;
  while (c < icoarse)
  {
    c*=2;
    coarsenings++;
  }
  CH_assert(c == icoarse); // check icoarse for power of 2

  compact();
  if (m_tree.nodes == 0) return;
  if (m_tree.nodes == &full)
  {
    m_spanBox.coarsen(icoarse);
    return;
  }
  if (coarsenings > m_depth)
    {
      clearTree(m_tree);
      m_tree.nodes = &full;
      m_spanBox.coarsen(icoarse);
      m_minBox.coarsen(icoarse);
      return;
    }

  // first, pure tree manipulation, turn all unit sized terminal nodes
  // into gathered nodes.
  //      ie.
  // if depth == m_depth-1  && TreeNode.node != 0 or full
  //     TreeNode.nodes = [full,0,0,full]  -> TreeNode.nodes = full

  parents[0] = &m_tree;

  index[1] = 0;
  index[0] = 0;
  int depth = 1;
  while (depth != 0)
    {
      TreeNode& current = parents[depth-1]->nodes[index[depth]];
      if (depth == m_depth-coarsenings-1)
        {
          if (current.nodes != 0)
            {
              clearTree(current);
              current.nodes = &full;
            }
        }
      else if (!(current.nodes == 0 || current.nodes == &full))
        {
          parents[depth] = &current;
          depth++;
          index[depth] = -1;
        }
      nextNode(depth);
    }

  m_spanBox.coarsen(icoarse);
  m_minBox.coarsen(icoarse);
  m_depth-=coarsenings;

}

void TreeIntVectSet::trimCoarsen(int icoarse)
{
  if (icoarse == 1) return;
  CH_assert(icoarse >= 1);
  int coarsenings = 1;
  int c = 2;
  while (c < icoarse)
  {
    c*=2;
    coarsenings++;
  }
  CH_assert(c == icoarse); // check icoarse for power of 2

  compact();
  if (m_tree.nodes == 0) return;
  if (m_tree.nodes == &full)
  {
    m_spanBox.coarsen(icoarse);
    return;
  }
  if (coarsenings > m_depth)
    {
      clearTree(m_tree);
      m_tree.nodes = &full;
      m_spanBox.coarsen(icoarse);
      m_minBox.coarsen(icoarse);
      return;
    }

  // first, pure tree manipulation, turn all unit sized terminal nodes
  // into gathered nodes.  This prunes tree of all partial leaves at m_depth
  //      ie.
  // if depth == m_depth-1  && TreeNode.node != full
  //     TreeNode.nodes = [full,0,0,full]  -> TreeNode.nodes = 0

  parents[0] = &m_tree;

  index[1] = 0;
  index[0] = 0;
  int depth = 1;
  while (depth != 0)
    {
      TreeNode& current = parents[depth-1]->nodes[index[depth]];
      if (depth == m_depth-coarsenings-1)
        {
          if (current.nodes != &full)
            {
              clearTree(current);
              current.nodes = 0;
            }
        }
      else if (!(current.nodes == 0 || current.nodes == &full))
        {
          parents[depth] = &current;
          depth++;
          index[depth] = -1;
        }
      nextNode(depth);
    }

  m_spanBox.coarsen(icoarse);
  m_minBox.coarsen(icoarse);
  m_depth-=coarsenings;

}

void TreeIntVectSet::expandNode(TreeNode& a_node)
{
  // Vector<int>&  indexRef = index;
  TreeNode* tmp = a_node.nodes;
  CH_assert(tmp == 0 || tmp == &full);
  //a_node.nodes = new TreeNode[TIVS_NODESIZE];
  a_node.nodes = (TreeNode*)(treeNodePool->getPtr());
  //printf("alloc e: %p\n",a_node.nodes);
  for (int i=0; i<TIVS_NODESIZE; ++i)
    {
      a_node.nodes[i].nodes = tmp;
    }

}

TreeIntVectSet& TreeIntVectSet::operator|=(const Box& a_box)
{
  if (a_box.isEmpty()) return *this;

  if (m_tree.nodes == 0 || m_tree.nodes == &full)
    {
      expandNode(m_tree);
    }
  m_minBox.minBox(a_box);

  // one half of the real smarts of this class, the other half being
  // the growTree function

  while (!m_spanBox.contains(a_box)) growTree();

  //Vector<int>& indexRef = index; //used for debugging
  //Vector<Box>& boxRef   = boxes; //  ..............
  index[1] = 0;
  parents[0] = &(m_tree);
  boxes[0] = m_spanBox;
  int depth = 1;

  while (depth != 0)
    {
      TreeNode* parent = parents[depth-1];
      int ind = index[depth];
      TreeNode& current = parent->nodes[ind];
      Box& parentBox  = boxes[depth-1];
      Box& currentBox = boxes[depth];
      quadrantBox(parentBox, index[depth], currentBox);
      // clear up the two most common cases quickly

      // case 0 a_box does not touch this quadrant, or this whole
      //quadrant is already marked 'full'...no change
      if (current.nodes == &full || !a_box.intersectsNotEmpty(currentBox))
        {
          nextNode(depth);
        }
      // case 1 currentBox is completely covered by a_box.
      else if (currentBox.smallEnd() >= a_box.smallEnd() &&
              currentBox.bigEnd()   <= a_box.bigEnd())
        {
          clearTree(current);
          current.nodes = &full;
          nextNode(depth);
        }

      // OK, now things are more tricky, partial intersection of
      // boxes.  Need to parse deeper
      else
        {
          if (current.nodes == 0 || current.nodes == &full)
            expandNode(current);
          parents[depth] = &(current);
          depth++;
          index[depth]= 0;
        }
    }

  return *this;

}

void TreeIntVectSet::transfer(TreeNode& a_node, const Box& a_box)
{
  while (!m_spanBox.contains(a_box)) growTree();
  index[1] = 0;
  parents[0] = &(m_tree);
  boxes[0] = m_spanBox;
  int depth = 1;

  while (depth != 0)
    {
      TreeNode* parent = parents[depth-1];
      int ind = index[depth];
      TreeNode& current = parent->nodes[ind];
      Box& parentBox  = boxes[depth-1];
      Box& currentBox = boxes[depth];
      quadrantBox(parentBox, index[depth], currentBox);
      if (currentBox == a_box)
        {
          clearTree(current);
          current = a_node;
          return;
        }
      if (currentBox.contains(a_box))
        {
          if (current.nodes == 0 || current.nodes == &full)
            expandNode(current);
          parents[depth] = &(current);
          depth++;
          index[depth]= 0;
        }
      else
        {
          nextNode(depth);
        }
    }
}

TreeIntVectSet& TreeIntVectSet::operator|=(const TreeIntVectSet& set)
{
  if (set.m_tree.nodes==0) return *this;
  if (set.m_tree.nodes == &full)
    {
      this->operator|=(set.m_spanBox);
      return *this;
    }
  if (m_tree.nodes == 0)
    {
      *this = set;
      return *this;
    }
  if (m_tree.nodes == &full)
    {
      Box b = m_spanBox;
      *this = set;
      *this |= b;
      return *this;
    }

  const TreeNode* yourRoot;
  TreeIntVectSet tmp;

  if (m_depth == set.m_depth)
    {
      yourRoot = &(set.m_tree);
    }
  else if (m_depth < set.m_depth)
    {
      while (m_depth < set.m_depth) growTree();
      yourRoot = &(set.m_tree);
    }
  else
    {
      tmp = set;
      while (tmp.m_depth < m_depth) tmp.growTree();
      yourRoot = &(tmp.m_tree);
    }

  static Vector<const TreeIntVectSet::TreeNode*> parents_other, parents_this;
  static Vector<int> indexSet;
  if (parents.size() > parents_other.size())
    {
      parents_other.resize(parents.size());
      parents_this.resize(parents.size());
      indexSet.resize(parents.size());
    }

  parents_other[0] = yourRoot;
  parents_this[0] = &m_tree;
  indexSet[1] = 0;
  int depth = 1;
  while (depth!=0)
    {
      const TreeNode* otherParent  = parents_other[depth-1];
      TreeNode* thisParent   = (TreeNode*)parents_this[depth-1];
      const TreeNode& otherCurrent = otherParent->nodes[indexSet[depth]];
      TreeNode& thisCurrent  = thisParent->nodes[indexSet[depth]];

      if (otherCurrent.nodes == 0 || thisCurrent.nodes == &full)
        {
          NEXT(indexSet, depth);
        }
      else if (otherCurrent.nodes == &full)
        {
          clearTree(thisCurrent);
          thisCurrent = otherCurrent;
          // leave the rest of this node alone
          NEXT(indexSet, depth);
        }
      else if (thisCurrent.nodes == 0)
        {
          cloneNode(otherCurrent, thisCurrent);
          NEXT(indexSet, depth);
        }
      else  // not terminal node for either one, traverse deeper
        {
          parents_other[depth] = &(otherCurrent);
          parents_this[depth]  = &(thisCurrent);
          depth++;
          indexSet[depth] = 0;
        }
    }
  m_minBox.minBox(set.m_minBox);

  return *this;
}

bool TreeIntVectSet::operator==(const TreeIntVectSet& set) const
{

  if (m_depth  != set.m_depth) return false;

  if (set.m_tree.nodes == 0 || set.m_tree.nodes == &full)
    {
      return set.m_tree.nodes == m_tree.nodes;
    }

  static Vector<const TreeIntVectSet::TreeNode*> parents_other, parents_this;
  static Vector<int> indexSet;
  if (parents.size() > parents_other.size())
    {
      parents_other.resize(parents.size());
      parents_this.resize(parents.size());
      indexSet.resize(parents.size());
    }

  parents_other[0] = &(set.m_tree);
  parents_this[0] = &m_tree;
  indexSet[1] = 0;
  int depth = 1;
  while (depth!=0)
    {
      const TreeNode* otherParent  = parents_other[depth-1];
      TreeNode* thisParent   = (TreeNode*)parents_this[depth-1];
      const TreeNode& otherCurrent = otherParent->nodes[indexSet[depth]];
      TreeNode& thisCurrent  = thisParent->nodes[indexSet[depth]];

      if (otherCurrent.nodes == 0 || otherCurrent.nodes == &full)
        {
          if (thisCurrent.nodes != otherCurrent.nodes) return false;
          NEXT(indexSet, depth);
        }
      else  // not terminal node for otherIVS, traverse deeper
        {
          if (thisCurrent.nodes == 0 || thisCurrent.nodes == &full) return false;
          parents_other[depth] = &(otherCurrent);
          parents_this[depth]  = &(thisCurrent);
          depth++;
          indexSet[depth] = 0;
        }
    }

  return true;
}

bool TreeIntVectSet::operator<(const TreeIntVectSet& a_tivs) const
{
  //
  // Primary criterion: number of IntVects.
  //
  int n0 = this->numPts();
  int n1 = a_tivs.numPts();
  if (n0 < n1)
  {
    return true;
  } else
  if (n0 > n1)
  {
    return false;
  }

  //
  // Secondary criterion: lexicographic comparison of the IntVects, as traversed
  // by TreeIntVectSetIterator.
  //
  TreeIntVectSetIterator it0(*this);
  TreeIntVectSetIterator it1(a_tivs);
  for (; it0.ok(); ++it0, ++it1 )
  {
    CH_assert( it1.ok() );
    if ( it0().lexLT( it1() ) )
    {
      return true;
    } else
    if ( it1().lexLT( it0() ) )
    {
      return false;
    }
  }
  return false;
}

// old version that didn't rely on common centering techniques.....
//
// TreeIntVectSet& TreeIntVectSet::operator|=(const TreeIntVectSet& set)
// {
//
//   if (set.m_tree.nodes==0) return *this;
//   else if (set.m_tree.nodes == &full)
//     {
//       this->operator|=(set.m_spanBox);
//       return *this;
//     }
//
//   //m_minBox.minBox(set.m_minBox);
//
//   int depth = set.m_depth+1;
//   static Vector<int> index_local;
//   static Vector<const TreeIntVectSet::TreeNode*> parents_local;
//   static Vector<Box> boxes_local;
//   if (depth > index_local.size())
//     {
//       index_local.resize(depth);
//       parents_local.resize(depth);
//       boxes_local.resize(depth);
//     }
//   index_local[0] = 0;
//
//   depth = 1;
//   parents_local[0] = &(set.m_tree);
//   boxes_local[0] = set.m_spanBox;
//   index_local[1] = 0;
//   while (depth != 0)
//     {
//       const TreeNode* parent =  parents_local[depth-1];
//       int ind = index_local[depth];
//       const TreeNode& current = parent->nodes[ind];
//       quadrantBox(boxes_local[depth-1], ind, boxes_local[depth]);
//       if (current.nodes == 0 || current.nodes == &full)
//              {
//                if (current.nodes == 0){;}  // do nothing
//                else if (current.nodes == &TreeIntVectSet::full)
//                      {
//                        this->operator|=(boxes_local[depth]);
//                      }
//                index_local[depth]++;
//                while (index_local[depth] == TIVS_NODESIZE)
//                      {
//                        index_local[depth] = 0;
//                        depth--;
//                        index_local[depth]++;
//                      }
//              }
//       else{
//              parents_local[depth] = &(current);
//              depth++;
//              index_local[depth] = 0;
//       }
//     }
//   return *this;
// }

TreeIntVectSet& TreeIntVectSet::operator-=(const TreeIntVectSet& set)
{

  if (set.m_tree.nodes==0) return *this;
  else if (set.m_tree.nodes == &full)
    {
      this->operator-=(set.m_spanBox);
      return *this;
    }
  if (m_tree.nodes == 0) return *this;

  int depth = set.m_depth+1;
  static Vector<int> index_local;
  static Vector<const TreeIntVectSet::TreeNode*> parents_local;
  static Vector<Box> boxes_local;
  if (depth > index_local.size())
    {
      index_local.resize(depth);
      parents_local.resize(depth);
      boxes_local.resize(depth);
    }
  index_local[0] = 0;
  depth = 1;
  parents_local[0] = &(set.m_tree);
  boxes_local[0] = set.m_spanBox;
  index_local[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents_local[depth-1];
      int ind = index_local[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes_local[depth-1], ind, boxes_local[depth]);
      if (boxes_local[depth].intersectsNotEmpty(m_minBox))
        {
          if (current.nodes == 0 || current.nodes == &full)
            {
              if (current.nodes == 0)
              {
                // do nothing
              }
              else if (current.nodes == &TreeIntVectSet::full)
                {
                  this->operator-=(boxes_local[depth]);
                }
              index_local[depth]++;
              while (index_local[depth] == TIVS_NODESIZE)
                {
                  index_local[depth] = 0;
                  depth--;
                  index_local[depth]++;
                }
            }
          else
          {
            parents_local[depth] = &(current);
            depth++;
            index_local[depth] = 0;
          }
        }
      else
        {
          index_local[depth]++;
          while (index_local[depth] == TIVS_NODESIZE)
            {
              index_local[depth] = 0;
              depth--;
              index_local[depth]++;
            }
        }
    }
  return *this;
}

TreeIntVectSet& TreeIntVectSet::operator&=(const TreeIntVectSet& set)
{
  if (set.m_tree.nodes==0 || m_tree.nodes == 0 || !m_minBox.intersects(set.m_minBox))
    {
      clear();
      return *this;
    }

  if (set.m_tree.nodes==&full)
  {
    return *this&=set.m_spanBox;
  }
  else if (m_tree.nodes == &full)
    {
      expandNode(m_tree);
    }

  const TreeNode* yourRoot;
  TreeIntVectSet tmp;

  if (m_depth == set.m_depth)
    {
      yourRoot = &(set.m_tree);
    }
  else if (m_depth < set.m_depth)
    {
      while (m_depth < set.m_depth) growTree();
      yourRoot = &(set.m_tree);
    }
  else
    {
      tmp = set;
      while (tmp.m_depth < m_depth) tmp.growTree();
      yourRoot = &(tmp.m_tree);
    }

  static Vector<const TreeIntVectSet::TreeNode*> parents_other, parents_this;
  static Vector<int> indexSet;
  if (parents.size() > parents_other.size())
    {
      parents_other.resize(parents.size());
      parents_this.resize(parents.size());
      indexSet.resize(parents.size());
    }

  parents_other[0] = yourRoot;
  parents_this[0] = &m_tree;
  indexSet[1] = 0;
  int depth = 1;
  while (depth!=0)
    {
      const TreeNode* otherParent  = parents_other[depth-1];
      TreeNode* thisParent   = (TreeNode*)parents_this[depth-1];
      const TreeNode& otherCurrent = otherParent->nodes[indexSet[depth]];
      TreeNode& thisCurrent  = thisParent->nodes[indexSet[depth]];

      if (otherCurrent.nodes == 0 || thisCurrent.nodes == 0)
        {
          clearTree(thisCurrent);
          NEXT(indexSet, depth);
        }
      else if (otherCurrent.nodes == &full)
        {
          // leave the rest of this node alone
          NEXT(indexSet, depth);
        }
      else if (thisCurrent.nodes == &full) // swap nodes
        {
          clearTree(thisCurrent);
          cloneNode(otherCurrent, thisCurrent);
          NEXT(indexSet, depth);
        }
      else  // not terminal node for either one, traverse deeper
        {
          parents_other[depth] = &(otherCurrent);
          parents_this[depth]  = &(thisCurrent);
          depth++;
          indexSet[depth] = 0;
        }
    }
  m_minBox &= set.m_minBox;

  return *this;
}

/*
TreeIntVectSet& TreeIntVectSet::operator&=(const TreeIntVectSet& set)
{

  if (set.m_tree.nodes==0 || m_tree.nodes == 0 || !m_minBox.intersects(set.m_minBox))
    {
      clear();
      return *this;
    }
  if (set.m_tree.nodes==&full)
    {
      return *this&=set.m_spanBox;
    }
  else if (m_tree.nodes == &full)
    {
      expandNode(m_tree);
    }

  // step one, trim by be minbox;
  *this &= set.m_minBox;

  // now, subtract the complement of the
  int depth = set.m_depth+1;
  static Vector<int> index_local;
  static Vector<const TreeIntVectSet::TreeNode*> parents_local;
  static Vector<Box> boxes_local;
  if (depth > index_local.size())
    {
      index_local.resize(depth);
      parents_local.resize(depth);
      boxes_local.resize(depth);
    }
  index_local[0] = 0;

  depth = 1;
  parents_local[0] = &(set.m_tree);
  boxes_local[0] = set.m_spanBox;
  index_local[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents_local[depth-1];
      int ind = index_local[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes_local[depth-1], ind, boxes_local[depth]);
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              *this -= boxes_local[depth];
            }
          else if (current.nodes == &TreeIntVectSet::full)
            {
              //do nothing.
            }
          index_local[depth]++;
          while (index_local[depth] == TIVS_NODESIZE)
            {
              index_local[depth] = 0;
              depth--;
              index_local[depth]++;
            }
        }
      else
        {
          parents_local[depth] = &(current);
          depth++;
          index_local[depth] = 0;
        }
    }
  return *this;
}
*/

// first pass to help matters a little, next pass can exploit
// the common centering aspect of the tmp and *this tree.  bvs.
void TreeIntVectSet::nestingRegion(int radius, const Box& a_domain, int granularity)
{
  CH_assert(radius>=0);
  if (radius == 0) return;

  Box domain = a_domain;

  if (granularity != 1)
    {
      int nradius = radius/granularity;
      //CH_assert(radius == nradius*granularity);
      if (radius != nradius*granularity) nradius++;
      radius = nradius;
      this->trimCoarsen(granularity);
      domain.coarsen(granularity);
    }

  TreeIntVectSet tmp(*this);
  {
    IntVect lo = m_minBox.smallEnd();
    IntVect hi = m_minBox.bigEnd();
    for (int i=0; i<SpaceDim; ++i)
      {
        if (lo[i] != domain.smallEnd()[i]) lo[i] += radius;
        if (hi[i] != domain.bigEnd()[i])   hi[i] -= radius;
      }
    if (hi<lo)
      {
        this->clear();
        return;
      }
    else
      {
        Box shrink(lo, hi);
        *this &= shrink;
      }
  }

  IntVect lower(-radius * IntVect::Unit), upper(IntVect::Unit * radius);

  int depth  = m_depth+1;
  static Vector<int> index_tmp;
  static Vector<const TreeIntVectSet::TreeNode*> parents_tmp;
  static Vector<Box> boxes_tmp;
  if (depth > index_tmp.size())
    {
      index_tmp  .resize(depth+1);
      parents_tmp.resize(depth+1);
      boxes_tmp  .resize(depth+1);
    }
  index_tmp[0] = 0;

  depth = 1;
  parents_tmp[0] = &(tmp.m_tree);
  boxes_tmp[0] = m_spanBox;
  index_tmp[1] = 0;
  Box clobberBox;
  IntVect& lo = (IntVect&)(clobberBox.smallEnd());
  IntVect& hi = (IntVect&)(clobberBox.bigEnd());
  while (depth != 0)
    {
      const TreeNode* parent =  parents_tmp[depth-1];
      int ind = index_tmp[depth];
      Box&  currentBox = boxes_tmp[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes_tmp[depth-1], ind, currentBox);
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              lo = currentBox.smallEnd();
              hi = currentBox.bigEnd();
              lo.max(domain.smallEnd());
              hi.min(domain.bigEnd());
              if (hi>=lo)
                {
                  hi+=upper;
                  lo+=lower;
                  clobberBox.computeBoxLenNotEmpty();
                  *this -= clobberBox;
                }
            }
          else if (current.nodes == &TreeIntVectSet::full)
            {
              //do nothing.
            }
          index_tmp[depth]++;
          while (index_tmp[depth] == TIVS_NODESIZE)
            {
              index_tmp[depth] = 0;
              depth--;
              index_tmp[depth]++;
            }
        }
      else
        {
          parents_tmp[depth] = &(current);
          depth++;
          index_tmp[depth] = 0;
        }
    }
  if (granularity != 1)
    {
      this->refine(granularity);
    }
}

void TreeIntVectSet::nestingRegion(int radius, const ProblemDomain& a_domain, int granularity)
{
  CH_assert(radius>=0);
  if (radius == 0) return;

  ProblemDomain domain = a_domain;

  if (granularity != 1)
    {
      int nradius = radius/granularity;
      if (nradius != radius*granularity) nradius++;
      radius = nradius;
      this->trimCoarsen(granularity);
      domain.coarsen(granularity);
    }

  const Box& domainBox = domain.domainBox();

  Box domainImaged = domainBox;

  TreeIntVectSet tmp(*this);
  {
    IntVect lo = m_minBox.smallEnd();
    IntVect hi = m_minBox.bigEnd();
    for (int i=0; i<SpaceDim; ++i)
      {
                if (lo[i] != domainBox.smallEnd()[i]) lo[i] += radius;
                if (hi[i] != domainBox.bigEnd()[i])   hi[i] -= radius;
                if (domain.isPeriodic(i))
                {
                  domainImaged.grow(i, radius);
                }
      }
    if (hi<lo)
      {
        this->clear();
        return;
      }
    else
      {
        Box shrink(lo, hi);
        *this &= shrink;
      }
  }

  // secret of periodic handling.  create phantom IntVects in the
  // tmp TreeIntVectSet for the periodic images up to radius distance
  if (domain.isPeriodic())
    {
      Vector<Box> boxRep = tmp.createBoxes();
      ShiftIterator shiftIt = domain.shiftIterator();
      Box intersectBox(domainBox);
      for (int dir=0; dir<SpaceDim; dir++)
        {
          if (domain.isPeriodic(dir))
            intersectBox.grow(dir,-radius);
        }
      IntVect shiftMult(domainBox.size());
      for (int b=0; b<boxRep.size(); ++b)
        {
          Box image = boxRep[b];
          if (!intersectBox.contains(image))
            {
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect(shiftMult*shiftIt());
                  image.shift(shiftVect);
                  if (image.intersects(domainImaged))
                    {
                      tmp |= (image & domainImaged);
                    }
                  image.shift(-shiftVect);
                } // end loop over periodic shifts
            } // end if box needs to be checked for periodic intersections
        } // end loop over boxes

      // this is designed to catch case where this level
      // abuts the physical domain, but there are no periodically
      // shifted boxes to subtract (dfm 11/13/01)
      Box tmpMinBox = tmp.m_minBox;
      IntVect loVect = tmpMinBox.smallEnd();
      IntVect hiVect = tmpMinBox.bigEnd();

      for (int dir=0; dir<SpaceDim; dir++)
      {
        if (domain.isPeriodic(dir))
          {
            if (tmpMinBox.smallEnd()[dir] == domainBox.smallEnd(dir))
              loVect[dir] += radius;
            if (tmpMinBox.bigEnd()[dir] == domainBox.bigEnd(dir))
              hiVect[dir] -= radius;
          } // end if periodic in this direction
      } // end loop over directions

      Box testBox(loVect,hiVect);
      // if testBox != tmpMinBox then at least one box in this IntVect
      // abuts a periodic boundary, but there are no periodic images
      // which abut it.  So, this edge must be shrunk as well
      if (testBox != tmpMinBox)
        *this &= testBox;

    } // end if is periodic

  IntVect lower(-radius * IntVect::Unit), upper(IntVect::Unit * radius);

  int depth  = m_depth+1;
  Vector<int> index_tmp;
  Vector<const TreeIntVectSet::TreeNode*> parents_tmp;
  Vector<Box> boxes_tmp;
  if (depth > index_tmp.size())
    {
      index_tmp  .resize(depth+1);
      parents_tmp.resize(depth+1);
      boxes_tmp  .resize(depth+1);
    }
  index_tmp[0] = 0;

  depth = 1;
  parents_tmp[0] = &(tmp.m_tree);
  boxes_tmp[0] = tmp.m_spanBox;
  index_tmp[0] = 0;
  index_tmp[1] = 0;
  Box clobberBox;
  IntVect& lo = (IntVect&)(clobberBox.smallEnd());
  IntVect& hi = (IntVect&)(clobberBox.bigEnd());
  while (depth != 0)
    {
      const TreeNode* parent =  parents_tmp[depth-1];
      int ind = index_tmp[depth];
      Box&  currentBox = boxes_tmp[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes_tmp[depth-1], ind, currentBox);
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              lo = currentBox.smallEnd();
              hi = currentBox.bigEnd();
              lo.max(domainImaged.smallEnd());
              hi.min(domainImaged.bigEnd());
              if (hi>=lo)
                {
                  hi+=upper;
                  lo+=lower;
                  clobberBox.computeBoxLenNotEmpty();
                  *this -= clobberBox;
                }
            }
          else if (current.nodes == &TreeIntVectSet::full)
            {
              //do nothing.
            }
          index_tmp[depth]++;
          while (index_tmp[depth] == TIVS_NODESIZE)
            {
              index_tmp[depth] = 0;
              depth--;
              index_tmp[depth]++;
            }
        }
      else
        {
          parents_tmp[depth] = &(current);
          depth++;
          index_tmp[depth] = 0;
        }
    }
  if (granularity != 1)
    {
      this->refine(granularity);
    }
}

bool TreeIntVectSet::contains(const IntVect& iv) const
{
  if (m_tree.nodes == 0 || !m_minBox.contains(iv)) return false;
  if (m_tree.nodes == &full) return true;

  //Vector<int>& indexRef = index; //used for debugging
  //Vector<Box>& boxRef   = boxes; //  ..............

  int depth = 1;
  parents[0] = (TreeNode*)&(m_tree);
  boxes[0] = m_spanBox;
  index[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes[depth-1], ind, boxes[depth]);
      if (!boxes[depth].contains(iv))
        {
          index[depth]++;
          CH_assert(index[depth] != TIVS_NODESIZE);
        }
      else
        {
          if (current.nodes == 0) return false;
          if (current.nodes == &full) return true;
          parents[depth] = (TreeNode*)&(current);
          depth++;
          index[depth] = 0;
        }
    }
  return false;
}

bool TreeIntVectSet::isEmpty() const
{
  if (m_tree.nodes == 0) return true;
  if (m_tree.nodes == &full) return false;

  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;
  index[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeNode& current = parent->nodes[ind];
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
          {
            // do nothing
          }
          else if (current.nodes == &full)
          {
            return false;
          }
          nextNode(depth);
        }
      else
        {
          parents[depth] = (TreeNode*)&(current);
          depth++;
          index[depth] = 0;
        }
    }
  ((TreeIntVectSet*)this)->clear();
  return true;
}

TreeIntVectSet& TreeIntVectSet::operator&=(const ProblemDomain& domain)
{
  if (!domain.isPeriodic())
    {
      this->operator&=(domain.domainBox());
      return *this;
    }
  TreeIntVectSet tmp(*this);
  tmp -= domain.domainBox();
  *this &= domain.domainBox();
  TreeIntVectSetIterator it(tmp);
  for (it.begin(); it.ok(); ++it)
    {
      IntVect iv = it();
      if (domain.image(iv)) *this |= iv;
    }
  return *this;
}

TreeIntVectSet& TreeIntVectSet::operator&=(const Box& a_box)
{
  // one half of the real smarts of this class, the other half being
  // the growTree function

  if (m_tree.nodes == 0) return *this;
  if (a_box.contains(m_minBox)) return *this;
  if (m_tree.nodes == &full) expandNode(m_tree);

  //Vector<int>& indexRef = index; //used for debugging
  //Vector<Box>& boxRef   = boxes; //  ..............
  index[0] = 0;
  index[1] = 0;
  parents[0] = &(m_tree);
  boxes[0] = m_spanBox;
  int depth = 1;
  while (depth != 0)
    {
      TreeNode* parent = parents[depth-1];
      int ind = index[depth];
      TreeNode& current = parent->nodes[ind];
      Box& parentBox  = boxes[depth-1];
      Box& currentBox = boxes[depth];
      quadrantBox(parentBox, index[depth], currentBox);
      if (current.nodes == 0)
      {
        nextNode(depth);
      }
      else if (current.nodes == &full)
        {
          if (a_box.contains(currentBox))
            {
              nextNode(depth);
            }
          else if (!a_box.intersectsNotEmpty(currentBox))
            {
              clearTree(current);
              current.nodes = 0;
              nextNode(depth);
            }
          else
            {
              expandNode(current);
              parents[depth] = &(current);
              depth++;
              index[depth]= 0;
            }
        }
      // case 1 currentBox is completely covered by a_box.
      else
        {
          parents[depth] = &(current);
          depth++;
          index[depth]= 0;
        }
    }

  return *this;
}

TreeIntVectSet& TreeIntVectSet::operator-=(const Box& a_box)
{
  remove(a_box, this);
  return *this;
}

void TreeIntVectSet::remove(const Box& a_box, TreeIntVectSet* a_residual)
{
  if (m_tree.nodes == 0 || (a_box & m_minBox).isEmpty()) return;
  if (a_box.contains(m_minBox))
    {
      clear();
      return;
    }
  if (m_tree.nodes == &full)
    {
      expandNode(m_tree);
    }

  //Vector<int>& indexRef = index; //used for debugging
  //Vector<Box>& boxRef   = boxes; //  ..............
  std::list<Box> residual;
  std::list<TreeNode> tnodes;
  index[0] = 0;
  index[1] = 0;
  parents[0] = &(m_tree);
  boxes[0] = m_spanBox;
  int depth = 1;
  while (depth != 0)
    {
      TreeNode* parent = parents[depth-1];
      int ind = index[depth];
      TreeNode& current = parent->nodes[ind];
      Box& parentBox  = boxes[depth-1];
      Box& currentBox = boxes[depth];
      quadrantBox(parentBox, index[depth], currentBox);
      if (a_box.contains(currentBox))
        {
          if (a_residual != this)
            tnodes.push_back(current);
          else
            clearTree(current);
          current.nodes = 0;
          nextNode(depth);
          if (a_residual != this) residual.push_back(currentBox);
        }
      else if (current.nodes == 0)
        {
          nextNode(depth);
        }
      else if (current.nodes == &full)
        {
          if ((depth == m_depth-1) || !a_box.intersectsNotEmpty(currentBox))
            {
              nextNode(depth);
            }
          else
            {
              expandNode(current);
              parents[depth] = &(current);
              depth++;
              index[depth]= 0;
            }
        }
      else
        {
          parents[depth] = &(current);
          depth++;
          index[depth]= 0;
        }
    }
  compact();
  if (a_residual != this)
    {
      std::list<TreeNode>::iterator nodes = tnodes.begin();
      for (std::list<Box>::iterator it = residual.begin(); it!= residual.end(); ++it,++nodes)
        {
          a_residual->transfer(*nodes, *it);
        }
    }
}

void TreeIntVectSet::growTree()
{
  // do the grow operation

  // two cases. top node is terminal, or it has TIVS_NODESIZE leaves to parse out

  if (m_tree.nodes == 0 || m_tree.nodes == &full)
    {
      expandNode(m_tree);
    }
  else
    {
      for (int i=0; i<TIVS_NODESIZE; ++i)
        {
          TreeNode* tmp = m_tree.nodes[i].nodes;
          if (tmp != 0)
            {
              int ind = oppositeQuadrant(i);
              //m_tree.nodes[i].nodes = new TreeNode[TIVS_NODESIZE];
              TreeNode* newPtr = (TreeNode*)(treeNodePool->getPtr());
              m_tree.nodes[i].nodes = newPtr;
              for (int j=0; j< TIVS_NODESIZE; ++j) newPtr[j].nodes = 0;
              //printf("alloc: %p\n", m_tree.nodes[i].nodes );
              //printf("tmp: %p\n", tmp);
              newPtr[ind].nodes = tmp;
            }
        }

    }

  if (m_spanBox.isEmpty())
    m_spanBox = Box(-1*IntVect::Unit, IntVect::Zero);
  else
    m_spanBox.grow(m_spanBox.size(0)/2);

  m_depth++;
  int m = m_depth+1;
  if (m > index.size())
    {
      index.resize(m);
      parents.resize(m);
      boxes.resize(m);
          bufferOffset.resize(m);
    }
}

// Basically, a BoxIterator function without the crap of
// a BoxIterator.
bool TreeIntVectSet::nextIntVect(const Box& box, IntVect& iv)
{

  iv[0]++;
  if (iv[0] > box.bigEnd(0))
    {
#if (CH_SPACEDIM >= 2)
      iv[0] = box.smallEnd(0);
      iv[1]++;
#else
      return false;
#endif
    }
  else
    {
      return true;
    }

#if (CH_SPACEDIM >= 2)
  if (iv[1] > box.bigEnd(1))
    {
#if (CH_SPACEDIM >= 3)
      iv[1] = box.smallEnd(1);
      iv[2]++;
#else
      return false;
#endif
    }
  else
    {
      return true;
    }
#endif

#if (CH_SPACEDIM == 3)
  if (iv[2] > box.bigEnd(2))
    {
      return false;
    }
  return true;
#endif

  // for DIM> 3, implement the remaining directions as a loop
#if (CH_SPACEDIM > 3)
  for (int dir=2; dir<SpaceDim; dir++)
    {
      if (iv[dir] > box.bigEnd(dir))
        {
          if (dir == SpaceDim -1)
            {
              return false;
            }
          else
            {
              iv[dir] = box.smallEnd(dir);
              iv[dir+1]++;
            } // end if dir != SpaceDim-1
        } // end if we've stepped past limit in direction dir
    } // end loop over directions >= 2

  // these assertions are turned off in "regular" Chombo
  CH_assert(iv<=box.bigEnd());
  CH_assert(iv>=box.smallEnd());

  // if we made it here, than we haven't stepped off the box...
  return true;
#endif
}

void TreeIntVectSet::quadrantBox(const Box& inputBox,
                                 int quadrant,
                                 Box& outputQuadrant)
{
  int halfSize = inputBox.size(0)/2;
  if (halfSize == 0)
    {
      outputQuadrant = inputBox;
      return;
    }
  int q   = TIVS_NODESIZE/2;
  //IntVect low, hi;
  int* low = (int*)(outputQuadrant.loVect());
  int* hi  = (int*)(outputQuadrant.hiVect());
  if (quadrant < q)
    {
      low[0] = inputBox.smallEnd(0);
      hi[0]   = inputBox.smallEnd(0)+halfSize-1;
    }
  else
    {
      low[0] = inputBox.smallEnd(0)+halfSize;
      hi[0]   = inputBox.bigEnd(0);
      quadrant-=q;
    }

#if (CH_SPACEDIM >= 2)
  q/=2;
  if (quadrant < q)
    {
      low[1] = inputBox.smallEnd(1);
      hi[1]   = inputBox.smallEnd(1)+halfSize-1;
    }
  else
    {
      low[1] = inputBox.smallEnd(1)+halfSize;
      hi[1]   = inputBox.bigEnd(1);
      quadrant-=q;
    }
#endif

#if (CH_SPACEDIM >= 3)
  q/=2;
  if (quadrant < q)
    {
      low[2] = inputBox.smallEnd(2);
      hi[2]   = inputBox.smallEnd(2)+halfSize-1;
    }
  else
    {
      low[2] = inputBox.smallEnd(2)+halfSize;
      hi[2]   = inputBox.bigEnd(2);
      quadrant -= q;
    }
#endif

#if (CH_SPACEDIM >= 4)
  q/=2;
  if (quadrant < q)
    {
      low[3] = inputBox.smallEnd(3);
      hi[3]   = inputBox.smallEnd(3)+halfSize-1;
    }
  else
    {
      low[3] = inputBox.smallEnd(3)+halfSize;
      hi[3]   = inputBox.bigEnd(3);
      quadrant-=q;
    }
#endif

#if (CH_SPACEDIM >= 5)
  q/=2;
  if (quadrant < q)
    {
      low[4] = inputBox.smallEnd(4);
      hi[4]   = inputBox.smallEnd(4)+halfSize-1;
    }
  else
    {
      low[4] = inputBox.smallEnd(4)+halfSize;
      hi[4]   = inputBox.bigEnd(4);
      quadrant-=q;
    }
#endif

#if (CH_SPACEDIM >= 6)
  q/=2;
  if (quadrant < q)
    {
      low[5] = inputBox.smallEnd(5);
      hi[5]   = inputBox.smallEnd(5)+halfSize-1;
    }
  else
    {
      low[5] = inputBox.smallEnd(5)+halfSize;
      hi[5]   = inputBox.bigEnd(5);
      quadrant-=q;
    }
#endif

  CH_assert(!outputQuadrant.isEmpty());
  outputQuadrant.computeBoxLenNotEmpty();
}

inline
void tab(ostream& os, int tab)
{
  for (; tab>=0 ;tab--) os << " ";
}

int TreeIntVectSet::numPts() const
{
  if (isEmpty()) return 0;
  if (m_tree.nodes == &full) return m_spanBox.numPts();

  int rtn= 0;
  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;
  boxes[0] = m_spanBox;
  index[1] = 0;
  quadrantBox(boxes[depth-1], 0, boxes[depth]);
  long long int bxLength = boxes[depth].size(0);

  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeNode& current = parent->nodes[ind];
      //quadrantBox(boxes[depth-1], ind, boxes[depth]);
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              // do nothing
            }
          else if (current.nodes == &full)
            {
              //  rtn += boxes[depth].numPts();
              // CH_assert(bxNumPts == boxes[depth].numPts());
              rtn +=  D_TERM6(bxLength,*bxLength,*bxLength,*bxLength,*bxLength,*bxLength);
            }
          int d = depth;
          nextNode(depth);
          if (depth < d)
          {
            for (int i=depth; i< d; i++)
              {
                bxLength*=2;
              }
          }
        }
      else
        {
          parents[depth] = (TreeNode*)&(current);
          depth++;
          index[depth] = 0;
          bxLength/=2;
          CH_assert(bxLength != 0);
        }
    }
  return rtn;
}
Vector<Box> TreeIntVectSet::createBoxes() const
{
  static Vector<Box> buffer(200), rtn;
  int size;
  createBoxes(buffer, size);
  rtn.resize(size, Box());
  for (int i=0; i<size; ++i) rtn[i] = buffer[i];
  return rtn;
}

void TreeIntVectSet::createBoxes(Vector<Box>& rtn, int& size) const
{
  size = 0;
  if (m_tree.nodes == 0) return;
  if (m_tree.nodes == &full)
    {
      if (rtn.size() < 1) rtn.resize(1, m_minBox);
      else               rtn[0] = m_minBox;
      size = 1;
      return;
    }
  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;
  boxes[0] = m_spanBox;
  index[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes[depth-1], ind, boxes[depth]);
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              // do nothing
            }
          else if (current.nodes == &full)
            {
              if (rtn.size() > size)
                rtn[size] = boxes[depth];
              else
                rtn.push_back(boxes[depth]);
              ++size;
            }
          nextNode(depth);
        }
      else
        {
          parents[depth] = (TreeNode*)&(current);
          depth++;
          index[depth] = 0;
        }
    }
}

bool TreeIntVectSet::contains(const Box& box) const
{
  if (m_tree.nodes == 0) return false;
  if (!m_minBox.contains(box)) return false;
  if (m_tree.nodes == &full)
    {
      return true;
    }
  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;
  boxes[0] = m_spanBox;
  index[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes[depth-1], ind, boxes[depth]);
      if (!box.intersects(boxes[depth]))
        nextNode(depth);
      else if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              // OK, boxes[depth] is NOT in IntVectSet, but it does
              // intersect 'box', hence we do not 'contain' this box.
              return false;
            }
          else if (current.nodes == &full)
            {
              //do nothing
            }
          nextNode(depth);
        }
      else
      {
        parents[depth] = (TreeNode*)&(current);
        depth++;
        index[depth] = 0;
      }
    }
  return true;
}
void dumpTree(const TreeIntVectSet* set)
{
  if (set == NULL) return;
  if (set->m_tree.nodes == 0) return;
  if (set->m_tree.nodes == &TreeIntVectSet::full)
    {
      cout << set->m_spanBox.smallEnd()
           <<"..."<<set->m_spanBox.bigEnd()<<"\n";
      return;
    }
  int depth = set->m_depth+1;
  static Vector<int> index;
  static Vector<const TreeIntVectSet::TreeNode*> parents;
  static Vector<Box> boxes;
  if (depth > index.size())
    {
      index.resize(depth);
      parents.resize(depth);
      boxes.resize(depth);
    }
  depth = 1;
  parents[0] = &(set->m_tree);
  boxes[0] = set->m_spanBox;
  index[0] = 0;
  index[1] = 0;
  while (depth != 0)
    {
      const TreeIntVectSet::TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeIntVectSet::TreeNode& current = parent->nodes[ind];
      TreeIntVectSet::quadrantBox(boxes[depth-1], ind, boxes[depth]);
      if (current.nodes == 0 || current.nodes == &TreeIntVectSet::full)
        {
          if (current.nodes == 0)
            {
              cout<<depth;
              tab(cout, depth);
              cout <<"0\n";
            }
          else if (current.nodes == &TreeIntVectSet::full)
            {
              cout<<depth;
              tab(cout, depth);
              cout << boxes[depth].smallEnd()
                   <<"..."<<boxes[depth].bigEnd()<<"\n";
            }

          index[depth]++;
          while (index[depth] == TIVS_NODESIZE)
            {
              index[depth] = 0;
              depth--;
              index[depth]++;
            }
        }
      else
        {
          cout<<depth;
          tab(cout, depth);
          cout <<"+\n";
          parents[depth] = &(current);
          depth++;
          index[depth] = 0;
        }
    }
}

void TreeIntVectSet::recalcMinBox() const
{
  Box min = Box();
  if (m_tree.nodes == 0)
    {
      return;
    }
  if (m_tree.nodes == &full)
    {
      min = m_spanBox;
      return;
    }
  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;
  boxes[0] = m_spanBox;
  index[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeNode& current = parent->nodes[ind];
      quadrantBox(boxes[depth-1], ind, boxes[depth]);
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              // do nothing, this region is not part of IntVectSet
            }
          else if (current.nodes == &full)
            {
              min.minBox(boxes[depth]);
            }
          nextNode(depth);
        }
      else
      {
        parents[depth] = (TreeNode*)&(current);
        depth++;
        index[depth] = 0;
      }
    }
  (Box&)m_minBox = min;
}

int TreeIntVectSet::oppositeQuadrant(int ind)
{
  CH_assert(ind < TIVS_NODESIZE);
  CH_assert(ind >= 0);
  return TIVS_NODESIZE+~ind;
}

void TreeIntVectSet::compact() const
{
  if (m_tree.nodes == 0 || m_tree.nodes == &full)
    {
      (Box&)m_minBox=m_spanBox;
      return;
    }
  index[0] = 0;
  int depth = 1;
#if (CH_SPACEDIM <= 3)
  static Vector<Tuple<TreeNode*, 8> > flags;
#elif (CH_SPACEDIM == 4)
  static Vector<Tuple<TreeNode*, 16> > flags;
#elif (CH_SPACEDIM == 5)
  static Vector<Tuple<TreeNode*, 32> > flags;
#elif (CH_SPACEDIM == 6)
  static Vector<Tuple<TreeNode*, 64> > flags;
#else
  bad spacedim
#endif

  flags.resize(index.size());
  parents[0] = (TreeNode*)&(m_tree);
  index[1] = 0;

  while (depth != 0)
    {
      TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      //Vector<int>& indexRef = index;  //debug aids
      //Vector<TreeNode*>& parentRef = parents; //debug aids
      TreeNode& current = parent->nodes[ind];
      flags[depth][ind] = current.nodes;
      if (current.nodes == 0 || current.nodes == &full)
        {
          index[depth]++;
          while (index[depth] == TIVS_NODESIZE)
            {
              int i=0;
              while (i< TIVS_NODESIZE-1 && flags[depth][i] == flags[depth][i+1]) i++;
              if (i == TIVS_NODESIZE-1)
                {
                  TreeNode* parent = parents[depth-1];
                  clearTree(*parent);
                  parent->nodes = flags[depth][0];
                  flags[depth-1][index[depth-1]] = flags[depth][0] ;
                }
              index[depth] = 0;
              depth--;
              index[depth]++;
            }
        }
      else
        {
          parents[depth] = &(current);
          depth++;
          index[depth] = 0;
        }
    }
  // OK, we have the non-redundant tree, now can optimize a little.
  // Trees are origin centered, so can't re-center, but we can
  // see if reverse of growing makes sense

}

int TreeIntVectSet::linearSize() const
{
  compact();

  if (isEmpty())
  {
    return CH_XD::linearSize<int>(m_depth);
  }

  int size = 2* CH_XD::linearSize<Box>(m_minBox) + CH_XD::linearSize<int>(m_depth);

  if (m_tree.nodes == &full)
    return size;

  size += TIVS_NODESIZE * sizeof(int);

  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;
  index[1] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      const TreeNode& current = parent->nodes[ind];

      if (current.nodes == 0 || current.nodes == &full)
        {
          //size += sizeof(int);
          nextNode(depth);
        }
      else
        {
          size += TIVS_NODESIZE * sizeof(int);
          parents[depth] = (TreeNode*)&(current);
          depth++;
          index[depth] = 0;
        }
    }
  return size;

}

void TreeIntVectSet::linearOut(void* const a_outBuf) const
{
  unsigned char* buffer = (unsigned char*)a_outBuf;

  if (isEmpty())
  {
    int* b = (int*) buffer;
    *b = 0;
    return;
  }
  if (m_tree.nodes == &full)
  {
    CH_XD::linearOut<int>(buffer, -1);
  }
  else
  {
    CH_XD::linearOut<int>(buffer, m_depth);
  }
  buffer+= CH_XD::linearSize<int>(m_depth);

  CH_XD::linearOut<Box>(buffer, m_minBox);
  buffer += CH_XD::linearSize<Box>(m_minBox);

  CH_XD::linearOut<Box>(buffer, m_spanBox);
  buffer += CH_XD::linearSize<Box>(m_spanBox);

  if (m_tree.nodes == &full)
  {
    return;
  }

  int* buf = (int*)buffer;

  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;

  index[1] = 0;
  if (bufferOffset.size() < index.size()) bufferOffset.resize(index.size());
  bufferOffset[0] = 0;
  int nextFree = TIVS_NODESIZE;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int   ind = index[depth];
      int   ID  =  bufferOffset[depth-1] + ind;
      const TreeNode& current = parent->nodes[ind];
      if (current.nodes == 0 || current.nodes == &full)
        {
          if (current.nodes == 0)
            {
              buf[ID] = 0;
            }
          else if (current.nodes == &full)
            {
              buf[ID] = 1;
            }
          nextNode(depth);
        }
      else
      {
        buf[ID] = nextFree;
        bufferOffset[depth] = nextFree;
        nextFree += TIVS_NODESIZE;
        parents[depth] = (TreeNode*)&(current);
        depth++;
        index[depth] = 0;
      }
    }

}

void TreeIntVectSet::linearIn(const void* const inBuf)
{
  clear();
  unsigned char* buffer = (unsigned char*)inBuf;

  CH_XD::linearIn(m_depth, buffer);
  buffer+=CH_XD::linearSize(m_depth);
  if (m_depth == 0)
  {
    m_depth = 1;
    return;
  }
  if (m_depth == -1)
    {
      m_tree.nodes = &full;
    }
  if (m_depth != -1 && m_depth > index.size())
    {
      index.resize(m_depth);
      parents.resize(m_depth);
      boxes.resize(m_depth);
      bufferOffset.resize(m_depth);
    }
  CH_XD::linearIn(m_minBox, buffer);
  buffer+=CH_XD::linearSize(m_minBox);
  CH_XD::linearIn(m_spanBox, buffer);
  buffer+=CH_XD::linearSize(m_spanBox);

  if (m_depth == -1)
  {
    m_depth = 1;
    return;
  }
  int* buf = (int*)buffer;

  expandNode(m_tree);

  int depth = 1;
  parents[0] = (TreeNode*)&m_tree;

  index[1] = 0;
  if (bufferOffset.size() < index.size()) bufferOffset.resize(index.size());
  bufferOffset[0] = 0;
  while (depth != 0)
    {
      const TreeNode* parent =  parents[depth-1];
      int ind = index[depth];
      TreeNode& current = parent->nodes[ind];
      int bi = buf[bufferOffset[depth-1] + ind] ;
      if (bi == 0 || bi == 1)
        {
          if (bi == 0)
            {
              current.nodes = 0;
            }
          else if (bi == 1)
            {
              current.nodes = &full;
            }
          bufferOffset[depth]++;
          nextNode(depth);
        }
      else
      {
        bufferOffset[depth] = bi;
        expandNode(current);
        parents[depth] = (TreeNode*)&(current);
        depth++;
        index[depth] = 0;
      }
    }
}

//======================================================================

void TreeIntVectSetIterator::begin()
{
  if (m_ivs->m_tree.nodes == 0)
    {
      m_depth = -1;
      return;
    }
  boxes[0] = m_ivs->m_spanBox;
  TreeIntVectSet::quadrantBox(boxes[0], 0, boxes[1]);
  nodes[0] = (TreeIntVectSet::TreeNode*)&(m_ivs->m_tree);
  m_depth = 1;
  index[0] = 0;
  index[1] = -1;
  if ( m_ivs->m_tree.nodes == &TreeIntVectSet::full)
        {
          m_current = boxes[0].smallEnd();
          m_depth = 0;
          return;
        }
  else
  {
        findNextNode();
  }
  CH_assert((m_depth == -1) || (m_ivs->m_minBox.contains(m_current)));
}

void TreeIntVectSetIterator::findNext()
{
  if (TreeIntVectSet::nextIntVect(boxes[m_depth], m_current))
        return;
  else
        {
          findNextNode();
        }
}

void TreeIntVectSetIterator::findNextNode()
{
  static TreeIntVectSet::TreeNode* full = &TreeIntVectSet::full;
  index[m_depth]++;
  while (index[m_depth] == TIVS_NODESIZE)
        {
          m_depth--;
          index[m_depth]++;
        }
  while (m_depth > 0)
        {
          const TreeIntVectSet::TreeNode* parent =  nodes[m_depth-1];
          CH_assert(index[m_depth] >= 0);
          const TreeIntVectSet::TreeNode& current = parent->nodes[index[m_depth]];
          if (current.nodes != 0)
            TreeIntVectSet::quadrantBox(boxes[m_depth-1], index[m_depth], boxes[m_depth]);
          if (current.nodes == full)
                {
                  m_current = boxes[m_depth].smallEnd();
                  return;
                }
          else if (current.nodes == 0)
                {
                  index[m_depth]++;
                  while (index[m_depth] == TIVS_NODESIZE)
                        {
                          m_depth--;
                          index[m_depth]++;
                        }
                }
          else
          {
                nodes[m_depth] = &(current);
                m_depth++;
                index[m_depth] = 0;
          }
        }
  m_depth = -1;
}

Vector<Box> TreeIntVectSet::boxes(2);
TreeIntVectSet::TreeNode TreeIntVectSet::full;
Vector<int> TreeIntVectSet::index(2);
Vector<int> TreeIntVectSet::bufferOffset(2);
Vector<TreeIntVectSet::TreeNode*> TreeIntVectSet::parents(2);
Pool TreeIntVectSet::treeNodePoolObject(sizeof(TreeNode[TIVS_NODESIZE]), "TreeIntVectSet Pool");
Pool* TreeIntVectSet::treeNodePool = & TreeIntVectSet::treeNodePoolObject;
#include "NamespaceFooter.H"
