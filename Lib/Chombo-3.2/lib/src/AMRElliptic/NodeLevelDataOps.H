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

// We can assume that template class T has null construction.

#ifndef _NODELEVELDATAOPS_H_
#define _NODELEVELDATAOPS_H_

#ifdef CH_MPI
#include <string>
#include <sstream>
#endif
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "SPMD.H"
#include "NodeFArrayBox.H"
#include "LevelDataOps.H"
#include "NamespaceHeader.H"

// default copy constructor and assign are fine.

class NodeLevelDataOps: public LevelDataOps<NodeFArrayBox>
{
public:
  NodeLevelDataOps();

  virtual ~NodeLevelDataOps()
  {
  }

  virtual void incr( LevelData<NodeFArrayBox>& a_lhs,
                     const LevelData<NodeFArrayBox>& a_x,
                     Real a_scale) ;

  virtual void axby( LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_x,
                     const LevelData<NodeFArrayBox>& a_y, Real a_a, Real a_b) ;

  virtual void scale(LevelData<NodeFArrayBox>& a_lhs, const Real& a_scale) ;

protected:
  RefCountedPtr<DataFactory<NodeFArrayBox> > m_levelFactory;
};

#include "NamespaceFooter.H"

#endif
