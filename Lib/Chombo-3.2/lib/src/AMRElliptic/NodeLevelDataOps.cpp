#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// BVS, June 26, 2003

#ifdef CH_MPI
#include <string>
#include <sstream>
#endif
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "SPMD.H"
#include "NodeLevelDataOps.H"
#include "NamespaceHeader.H"

NodeLevelDataOps::NodeLevelDataOps()
  :LevelDataOps<NodeFArrayBox>()
{
}

#define ITER(a) for (DataIterator dit = a.dataIterator(); dit.ok(); ++dit) \
                {                                                          \
                   DataIndex d = dit();

#define ENDFOR(a)                                                          \
                }

void NodeLevelDataOps::incr( LevelData<NodeFArrayBox>&       a_lhs,
                             const LevelData<NodeFArrayBox>& a_rhs,
                             Real a_scale)
{
  int numcomp = a_lhs.nComp();
  int  startcomp = 0;
  ITER(a_lhs)
    //this line is diff because it has to be
    Box subbox = a_lhs[d].getFab().box();
    a_lhs[d].plus(a_rhs[d],  subbox, subbox, a_scale, startcomp, startcomp, numcomp);
  ENDFOR(a_lhs);
}

void NodeLevelDataOps:: axby( LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_x,
                              const LevelData<NodeFArrayBox>& a_y, Real a, Real b)
{
  //  CH_assert(a_lhs.disjointBoxLayout() == a_x.disjointBoxLayout());
  ITER(a_lhs)

     NodeFArrayBox& data = a_lhs[d];
     data.copy(a_x[d]);
     data *= a;
     data += a_y[d];

  ENDFOR(a_lhs);
}

void NodeLevelDataOps:: scale(LevelData<NodeFArrayBox>& a_lhs,
                              const Real& a_scale)
{
  ITER(a_lhs)
    a_lhs[d] *= a_scale;
  ENDFOR(a_lhs);
}

#undef ITER
#undef ENDFOR

#include "NamespaceFooter.H"
