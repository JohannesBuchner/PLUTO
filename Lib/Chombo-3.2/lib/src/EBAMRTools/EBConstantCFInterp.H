#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCONSTANTCFINTERP_H_
#define _EBCONSTANTCFINTERP_H_

#include "REAL.H"
#include "Box.H"
#include "RefCountedPtr.H"

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "CornerCopier.H"
#include "NamespaceHeader.H"

///
/**
   Class to do constant extrapolation of data onto ONE  ghost cell
   from neighboring internal data.
 */
class EBConstantCFInterp
{
public:
  ///
  EBConstantCFInterp(const DisjointBoxLayout& a_dbl,
                     const EBISLayout&        a_ebisl,
                     const ProblemDomain&     a_domain,
                     const IntVect&           a_nGhost);

  ///
  ~EBConstantCFInterp();

  ///
  /**
     fills ghost cells over coarse fine interpolation with constant extrapolation
   */
  void interpolate(LevelData<EBCellFAB>&      a_soln);

protected:

  DisjointBoxLayout               m_dbl;
  EBISLayout                      m_ebisl;
  ProblemDomain                   m_domain;
  IntVect                         m_nGhost;
  CornerCopier                    m_cornerCopier;
  ///
  void interpolate(EBCellFAB&            a_phi,
                   const DataIndex&      a_datInd,
                   int                   a_idir,
                   Side::LoHiSide        a_hiorlo);


private:

  //strong construction is forced here
  EBConstantCFInterp()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBConstantCFInterp(const EBConstantCFInterp& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBConstantCFInterp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
