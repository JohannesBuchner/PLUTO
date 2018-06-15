#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves july 2009


#ifndef _EBCFDATA_H_
#define _EBCFDATA_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "QuadCFInterp.H"
#include "BaseIVFAB.H"
#include "CornerCopier.H"
#include "VoFIterator.H"
#include "NamespaceHeader.H"

///
/**
   Stuff that objects that do eb/cf will need.
   Used by EBQuadCFInterp and EBTensorCFInterp.
   Not really all that useful outside their context.
 */
class EBCFData
{
public:

  ///
  virtual ~EBCFData();

  ///
  /**
  */
  EBCFData(const DisjointBoxLayout&       a_dblFine,
           const DisjointBoxLayout&       a_dblCoar,
           const EBISLayout&              a_ebislFine,
           const EBISLayout&              a_ebislCoar,
           const ProblemDomain&           a_domainCoar,
           const int&                     a_nref,
           const LayoutData<IntVectSet>&  a_cfivs,
           const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance(),
           bool a_doEBCFCrossing = true,
           bool a_doCornerEdgeIterators = true);



  static void getExtrapSigns(IntVect& a_signs, const IntVect& a_corner, const Box& a_grid) ;

  //deliberately public to avoid the get set stuff.
  EBISLayout        m_ebislCoar;
  EBISLayout        m_ebislCoarsenedFine;
  EBISLayout        m_ebislFine;
  DisjointBoxLayout m_gridsCoar;
  DisjointBoxLayout m_gridsCoarsenedFine;
  DisjointBoxLayout m_gridsFine;
  ProblemDomain     m_domainFine;
  ProblemDomain     m_domainCoar;
  int               m_refRat;

  bool              m_doEBCFCrossing;

  //cache the EB-CF vofiterator, these are on the fine index space
  LayoutData<VoFIterator> m_vofItEBCFLo[SpaceDim];
  LayoutData<VoFIterator> m_vofItEBCFHi[SpaceDim];
  LayoutData<VoFIterator> m_vofItCorners;
  LayoutData<VoFIterator> m_vofItEdges;
  LayoutData<IntVectSet>  m_cornerIVS;
  LayoutData<IntVectSet>  m_edgeIVS;

  //fine index space fine layout
  LayoutData<IntVectSet> m_ebcfivsLo[SpaceDim];
  LayoutData<IntVectSet> m_ebcfivsHi[SpaceDim];

  ///static for reuse
  static void getEBCFIVSGrid(IntVectSet&                a_ebcfivs,
                             const Box&                 a_grid,
                             const int&                 a_idir,
                             const Side::LoHiSide&      a_side,
                             const IntVect&             a_diagGrow,
                             const ProblemDomain&       a_domain,
                             const IntVectSet&          a_cfivs,
                             const EBISBox&             a_ebisBox);

  ///static for reuse
  static void
  getEdgeAndCornerIVS(IntVectSet& a_edgeIVS, IntVectSet& a_cornerIVS,
                      const Box& a_grid, const ProblemDomain& a_domain,
                      const IntVectSet& a_cfivsGrid);


private:

  void defineEdCoIterators(const LayoutData<IntVectSet>& a_cfivs);
  void defineLoHiIterators(const LayoutData<IntVectSet>& a_cfivs);
  ///
  /**
     Builds the IntVectSets for where we need to do
     EB aware CF interpolation.
     Returns a bool that indicates if we need to do EB aware CF interp
  */
  bool getEBCFIVS(const LayoutData<IntVectSet>& a_cfivs);


  //weak construction is discouraged
  EBCFData()
  {
    MayDay::Error("ebcfd invalid operator");
  }

  //disallowed for all the usual reasons
  EBCFData(const EBCFData& ebcin)
  {
    MayDay::Error("ebcfd invalid operator");
  }
  void operator=(const EBCFData& fabin)
  {
    MayDay::Error("ebcfd invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
