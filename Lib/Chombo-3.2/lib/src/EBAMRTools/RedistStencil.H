#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves weds oct 3 2001

#ifndef _REDISTSTENCIL_H_
#define _REDISTSTENCIL_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "BaseIVFAB.H"

#include "NamespaceHeader.H"

// This is temporary to test some things about redistribution.
extern bool g_alwaysRedist;

/// stencil for mass redistribution
/**
   The  RedistStencil class holds the stencil for mass redistribution
   at every irregular VoF in a layout.  The default weights
   that the stencil holds are volume weights.  The class
   does allow the flexibility to redefine these weights.
  */
class RedistStencil
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  RedistStencil();

  ///
  ~RedistStencil();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().
     ebisl needs ghost cells out to 3*redistradius
  */
  RedistStencil(const DisjointBoxLayout& a_dbl,
                const EBISLayout&        a_ebisl,
                const ProblemDomain&     a_domain,
                const int&               a_redistRadius);

  ///
  /**
     Defines this object.  Existing information is overriden.
     ebisl needs ghost cells out to 3*redistradius
  */
  void define(const DisjointBoxLayout& a_dbl,
              const EBISLayout&        a_ebisl,
              const ProblemDomain&     a_domain,
              const int&               a_redistRadius,
              bool                     a_do2DStencil = false);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Modify the weights in the stencil by multiplying by
     the inputs in a normalized way.
     If you want mass
     weighting, send in the density.
  */
  void resetWeights(const LevelData<EBCellFAB>& a_modifier,
                    const int&                  a_ivar);

  ///
  /**
     Returns the redistribution stencil at every irregular
     point in input  Box associated with  this DataIndex.
  */
  const BaseIVFAB<VoFStencil >&
  operator[](const DataIndex& datInd) const;

  ///
  /**
     This copy operation detaches copier from copyee.
   */
  void deepCopy(const RedistStencil& a_stenin);

  ///
  /**
   */
  int getRedistRadius() const;

protected:

  LayoutData<BaseIVFAB<VoFStencil > > m_stencil;
  LayoutData<BaseIVFAB<VoFStencil > > m_volsten;

  DisjointBoxLayout m_grids;

  ProblemDomain m_domain;

  EBISLayout m_ebisl;

  bool m_isDefined;
  bool m_hasDefaultWeights;
  bool m_alwaysRedist;

  int m_redistRadius;

private:
  //internal use
  void computePointStencil(VoFStencil&      a_stencil,
                           const VolIndex&  a_srcVoF,
                           const DataIndex& a_datInd,
                           const bool&      a_do2DStencil);

  //disallowed for all the usual reasons
  RedistStencil(const RedistStencil& ebcin)
  {
    MayDay::Error("rdsi 2 invalid operator");
  }
  void operator=(const RedistStencil& fabin)
  {
    MayDay::Error("rdsi 3 invalid operator");
  }
};

#include "NamespaceFooter.H"

#endif
