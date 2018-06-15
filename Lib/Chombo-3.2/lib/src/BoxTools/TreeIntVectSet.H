#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TREEINTVECTSET_H_
#define _TREEINTVECTSET_H_

#include "Box.H"
#include "Pool.H"
#include "NamespaceHeader.H"

class IntVectSet;
class ProblemDomain;

// xlC doesn't like the static const int inside this class.
// so instead of a nice variable we have a macro: TIVS_NODESIZE
// Currently, this only affects two files:  this one and the .cpp
#ifdef TIVS_NODESIZE
#undef TIVS_NODESIZE
#endif
#if (CH_SPACEDIM==1)
#define TIVS_NODESIZE 2
#elif (CH_SPACEDIM==2)
#define TIVS_NODESIZE 4
#elif (CH_SPACEDIM==3)
#define TIVS_NODESIZE 8
#elif (CH_SPACEDIM==4)
#define TIVS_NODESIZE 16
#elif (CH_SPACEDIM==5)
#define TIVS_NODESIZE 32
#elif (CH_SPACEDIM==6)
#define TIVS_NODESIZE 64
#else
//error -- only tested for dimensions 2 and 3 (ndk)
#error TIVS_NODESIZE is only defined for 1D, 2D, 3D, 4D, 5D, or 6D
#endif

/// IntVectSet implementation based on tree representation.
/**
   For explanations of these functions please look at IntVectSet class
   when the documentation doesn't appear here.

    @see IntVectSet

<HR>

<PRE>
 Further details of how non-recursive TreeNode design works:

 (for a 2D tree)

   (m_tree)
          + -- 0

            (a)+ - 0
                   1
               1   +
                   1  <------you are here

               + - + - 0
                   0   1
                   0   1
                   0   0

 for the node indicated, the 'index' vector  would contain

 index=[ 0      1    3  ...............]
 parents=[&m_tree &a   ..................]

or directly refered to as m_tree.nodes[1].nodes[3]

 the tree indicates a covering of an index space in either
 1 or 0.   1 is stored in the tree by pointing at the static
 'full' data member, 0 is stored as a 0.

 every 'nodes' member of a TreeNode object can be either

  0, &full, or a pointer .

  The interpretation of the tree depends on what m_spanBox is.
    nodes[i] indicates whether the i'th quadrant of the parent
    Box is full, empty, or needs to be parsed deeper.

</PRE>

 */
class TreeIntVectSet
{
public:
  ///
  inline TreeIntVectSet();

  ///
  TreeIntVectSet(const Box&);

  ///
  TreeIntVectSet(const TreeIntVectSet& a_sivs);

  ///
  ~TreeIntVectSet();

  ///
  void define(const Box&);

  ///
  void define(const TreeIntVectSet& a_sivs);

  /// trade internals of two objects
  void swap(TreeIntVectSet& a_other);

  ///
  TreeIntVectSet& operator=(const TreeIntVectSet& a_sivs);

  ///or
  TreeIntVectSet& operator|=(const TreeIntVectSet& a_sivs);

  ///
  TreeIntVectSet& operator|=(const IntVect& a_iv);

  ///
  TreeIntVectSet& operator|=(const Box& a_box);

  ///and
  TreeIntVectSet& operator&=(const TreeIntVectSet& s_sivs);

  ///and
  TreeIntVectSet& operator&=(const Box& a_box);

  ///and
  TreeIntVectSet& operator&=(const ProblemDomain& a_domain);

  ///not
  TreeIntVectSet& operator-=(const TreeIntVectSet& a_sivs);

  ///not
  TreeIntVectSet& operator-=(const IntVect& a_iv);

  ///not
  TreeIntVectSet& operator-=(const Box& a_box);

  ///returns true if
  bool operator==(const TreeIntVectSet& lhs) const;

  /**
      Primary sorting criterion: numPts().
      Secondary sorting criterion: lexigraphical ordering of the IntVects, taken
      in the order traversed by TreeIntVectSetIterator.
      In a total tie, returns false.
  */
  bool operator<(const TreeIntVectSet& a_sivs) const;

  /// Returns Vector<Box> representation of this IntVectSet.
  /**
     Returns Vector<Box> representation of this IntVectSet.
  */
  Vector<Box> createBoxes() const;

  /// Returns Vector<Box> representation of this IntVectSet.
  /**
     Returns Vector<Box> representation of this IntVectSet.
  */
  void createBoxes(Vector<Box>& boxes, int& size) const;

  ///
  bool contains(const IntVect& iv) const;

  ///
  bool contains(const Box& box) const;

  ///
  /**
     somewhat expensive, but shouldn't be ;-)  @see IntVectSet::chop
  */
  TreeIntVectSet chop(int idir, int chop_pnt);

  ///
  /**
     a proper faster chop function that does it all in-place and saves the giant memory cost.
  */
  void chop(int dir, int chop_pnt, TreeIntVectSet& a_hi);

  ///@see IntVectSet::grow
  /**
     expensive
  */
  void grow(int igrow);

  ///
  /**
     expensive
  */
  void grow(int idir, int igrow);

  ///
  /**
     @see IntVectSet::growHi
     expensive
  */
  void growHi();

  ///
  /**
     @see IntVectSet::growHi(int)
     expensive
  */
  void growHi(int a_dir);

  ///fast if iref is power of 2
  /**
     fast if iref is power of 2
  */
  void refine(int iref = 2);

  ///fast if iref is power of 2
  /**
     fast if iref is power of 2
  */
  void coarsen(int iref = 2);

  /// slow operation
  /**

  */
  void shift(const IntVect& iv);

  ///
  void clear();

  ///
  void nestingRegion(int a_radius, const Box& a_domain, int granularity);

  ///
  void nestingRegion(int a_radius, const ProblemDomain& a_domain, int granularity);

  ///
  inline const Box& minBox() const;

  ///
  bool isEmpty() const;

  ///
  int numPts() const;

  friend void dumpTree(const TreeIntVectSet* set);

  ///
  void compact() const;  // logically const operation;

  ///
  void recalcMinBox() const;

  /** \name Linearization routines */
  /*@{*/
  int linearSize() const;

  void linearIn(const void* const inBuf);

  void linearOut(void* const a_outBuf) const;
  /*@}*/

private:

  friend class TreeIntVectSetIterator;
  friend class MeshRefine;

#ifndef DOXYGEN
  struct TreeNode
  {
    TreeNode* nodes;

    TreeNode()
      :
      nodes(0)
    {
    }
  };
#endif

  TreeNode m_tree;
  Box m_minBox, m_spanBox;
  int m_depth;

  void trimCoarsen(int icoarse);

  //===============
  static Pool treeNodePoolObject;
  static Pool* treeNodePool;

  static void   quadrantBox(const Box& inputBox, int quadrant, Box& outputQuadrant);
  static void clearTree(TreeNode& tree);
  static void expandNode(TreeNode& node);
  void growTree();
  void remove(const Box& box, TreeIntVectSet* resdiual);
  void transfer(TreeNode& node, const Box& a_box);

  static int  oppositeQuadrant(int index);
  static bool nextIntVect(const Box& box, IntVect& iv);
  static void nextNode(int& currentDepth);
  static void cloneNode(const TreeNode& src, TreeNode& dest);
  // data structures that are needed again and again in the routines
  // that only change their size during grow() operations.  Keep them
  // around so that we don't have to create and destroy them on every function
  // call
  static Vector<int> index, bufferOffset;
  static Vector<TreeNode*> parents;
  static Vector<Box> boxes;

  // xlC wasn't like this
  // so instead of a nice const static integer, we have a MACRO TIVS_NODESIZE
  //const static int  nodeSize = D_TERM(2,*2,*2);
  static TreeNode full;
  friend struct Flag;
  friend class IntVectSet;
};

void dumpTree(const TreeIntVectSet* set); //in TreeIntVectSet.cpp

class TreeIntVectSetIterator
{
public:
  TreeIntVectSetIterator();
  TreeIntVectSetIterator(const TreeIntVectSet& ivs);
  void define(const TreeIntVectSet& ivs);
  const IntVect& operator()() const ;
  bool ok() const;
  void operator++();
  void begin();
  void end();
  void clear();
private:
  const TreeIntVectSet* m_ivs;
  Vector<const TreeIntVectSet::TreeNode*> nodes;
  Vector<Box>   boxes;
  Vector<int>   index;
  int           m_depth;
  IntVect       m_current;

  void findNextNode();
  void findNext(); // parse tree, starting at [depth, index], set m_current
  //when you find an entry, might run to end of iterator.
};

inline
TreeIntVectSetIterator::TreeIntVectSetIterator():m_ivs(0), m_depth(-1)
{
}

inline
TreeIntVectSetIterator::TreeIntVectSetIterator(const TreeIntVectSet& ivs)
{
  define(ivs);
}

inline
void TreeIntVectSetIterator::clear()
{
  m_depth = -1;
  m_ivs = 0;
}

inline
void TreeIntVectSetIterator::define(const TreeIntVectSet& ivs)
{
  m_ivs = &ivs;
  int max = ivs.index.size();
  //  if (max==10) printf("max=10 \n");
  if (boxes.size() < max)
    {
      boxes.resize(max);
      index.resize(max);
      nodes.resize(max);
    }
  begin();
}

inline
bool TreeIntVectSetIterator::ok() const
{
  return m_depth >= 0;
}

inline
void TreeIntVectSetIterator::end()
{
  m_depth = -1;
}

inline
TreeIntVectSet& TreeIntVectSet::operator|=(const IntVect& iv)
{
  return *this|=Box(iv, iv);
}

inline
TreeIntVectSet& TreeIntVectSet::operator-=(const IntVect& iv)
{
  return *this-=Box(iv, iv);
}

inline
void TreeIntVectSetIterator::operator++()
{
  findNext();
}

inline
const IntVect& TreeIntVectSetIterator::operator()() const
{
  return m_current;
}

//=======================================================

inline
TreeIntVectSet::~TreeIntVectSet()
{
  clearTree(m_tree);
}
inline
TreeIntVectSet::TreeIntVectSet(const TreeIntVectSet& a_tivs)
{
  define(a_tivs);
}

inline
TreeIntVectSet::TreeIntVectSet(const Box& a_box)
{
  define(a_box);
}

inline
TreeIntVectSet::TreeIntVectSet()
{
  m_tree.nodes = 0;
  m_depth=1;

}

inline
const Box&
TreeIntVectSet::minBox() const
{
  return m_minBox;
}

inline void TreeIntVectSet::nextNode(int& depth)
{
  index[0] = 0;
  index[depth]++;
  while (index[depth] == TIVS_NODESIZE)
    {
      index[depth] = 0;
      depth--;
      index[depth]++;
    }
}

#include "NamespaceFooter.H"
#endif //  TREEINTVECTSET_H
