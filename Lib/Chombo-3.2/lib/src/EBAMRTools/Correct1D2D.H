#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CORRECT1D2D_H_
#define _CORRECT1D2D_H_

#include "DisjointBoxLayout.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "LevelData.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

///
/**
   Given an EB application that uses one algorithm for fluxes on some
   boxes and another on other boxes, this will make the solution's fluxes
   match at box-box boundaries.   In all this class, the winning algorithm
   at the boundary is denoted by 2D, the losing by 1D but this would work
   for whatever polyalgorithm you choose.
 */
class Correct1D2D
{
public:
  static void
  makeIntMap(LevelData< BaseFab<int> >& a_intmap,
             const LayoutData<bool>&    a_is1D,
             const DisjointBoxLayout&   a_dbl);

  Correct1D2D(const EBLevelGrid& a_eblg,
              const LayoutData<bool>& a_is1D,
              int a_nvar);

  ~Correct1D2D()
  {;}


  ///
  /**
     sets buffers to zero
   */
  void setToZero();

  ///
  /**
     increments the 1D (losing) buffer by -flux*scale*sign(side)
     (side is which side of the changed cell we are talking about)
     typically scale = 1/dx[idir]
   */
  void increment1D(const EBFaceFAB& a_1DFlux,
                   const Real&      a_scale,
                   const DataIndex& a_dit);


  ///
  /**
     increments the 2D (winning) buffer by flux*scale*sign(side)
     (side is which side of the changed cell we are talking about)
     typically scale = 1/dx[idir]
   */
  void increment2D(const EBFaceFAB& a_2DFlux,
                   const Real&      a_scale,
                   const DataIndex& a_dit);


  ///
  /**
     subtracts off change in solution due to losing flux
     and adds in change in solution due to winning flux.
   */
  void correctSolution(LevelData<EBCellFAB>& a_U);

protected:

  LevelData< BaseFab<int> > m_1d2d;
  DisjointBoxLayout m_dbl1d[2*SpaceDim];
  DisjointBoxLayout m_dbl2d[2*SpaceDim];
  int m_nvar;
  EBLevelGrid m_eblg;

  //which VoFs actually get changed
  LayoutData< Vector<VolIndex> > m_sets1d[2*SpaceDim];
  LayoutData< Vector<VolIndex> > m_sets2d[2*SpaceDim];

  //change in solution due to 1d (losing) flux
  LevelData<EBCellFAB>  m_deltaU1d[2*SpaceDim];

  //change in solution due to 2d (winning) flux
  LevelData<EBCellFAB>  m_deltaU2d[2*SpaceDim];

  int index(int a_dir, Side::LoHiSide a_sd)
  {
    int retval = a_dir + a_sd*SpaceDim;
    return retval;
  }
private:
  Correct1D2D()
  {
    MayDay::Error("weak construction disallowed");
  }

  Correct1D2D(const Correct1D2D& a_in)
  {
    MayDay::Error("copy construction disallowed");
  }

  void operator=(const Correct1D2D& a_in)
  {
    MayDay::Error("assignment disallowed");
  }
};
#include "NamespaceFooter.H"
#endif
