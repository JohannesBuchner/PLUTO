#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri aug 17, 2001
// major revision dtg dec 2009
#ifndef _EBMGAVERAGE_H_
#define _EBMGAVERAGE_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

/// Piecewise constant interpolation.
/**
   Finedata += coarse data.
  */
class EBMGAverage
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBMGAverage();

  ///
  ~EBMGAverage();

  EBMGAverage(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout& a_ebislFine,
              const EBISLayout& a_ebislCoar,
              const ProblemDomain& a_domainCoar,
              const int& a_nref,
              const int& a_nvar,
              const EBIndexSpace* ebisPtr,
              const IntVect& a_ghostCellsRHS,
              const bool& a_layoutChanged = true);

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  void define(const DisjointBoxLayout& dblFine,
              const DisjointBoxLayout& dblCoar,
              const EBISLayout& ebislFine,
              const EBISLayout& ebislCoar,
              const ProblemDomain& domainCoar,
              const int& nref,
              const int& nvar,
              const EBIndexSpace* ebisPtr,
              const IntVect& ghostCellsRHS,
              const bool& a_layoutChanged = true);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Piecewise constant interpolation function
     FineData += coarseData
  */
  void
  average(LevelData<EBCellFAB>&       a_coarData,
          const LevelData<EBCellFAB>& a_fineData,
          const Interval&             a_variables);

  ///
  /**
     Piecewise constant interpolation function for unchanged layouts (no agglomeration)
     there are no copies here
     FineData += coarseData
  */
  void
  averageMG(LevelData<EBCellFAB>&       a_coarData,
            const LevelData<EBCellFAB>& a_fineData,
            const Interval&             a_variables);

protected:
  void
  averageFAB(EBCellFAB&       a_coar,
             const Box&       a_box,
             const EBCellFAB& a_refCoar,
             const DataIndex& a_datInd,
             const Interval&  a_variables) const;

  void  setDefaultValues();
  void  defineStencils();

  bool m_isDefined;
  IntVect m_ghost;
  DisjointBoxLayout m_coarGrids;
  DisjointBoxLayout m_fineGrids;
  DisjointBoxLayout m_refinedCoarseGrids;
  Copier m_copier;//this goes from m_refinedCoarseGrids to m_fineGrids
  ProblemDomain m_coarDomain;
  ProblemDomain m_fineDomain;

  EBISLayout m_coarEBISL;
  EBISLayout m_fineEBISL;

  int m_refRat;
  int m_nComp;

  bool m_layoutChanged;
  bool m_coarsenable;

  //LevelData<EBCellFAB> m_refinedCoarseData;
  //if coarsenable, these are the coarsened fine
  // otherwise, refined coarse.
  //actual buffer is no longer member data
  //not defined if !m_layoutChanged
  DisjointBoxLayout m_buffGrids;
  EBISLayout        m_buffEBISL;

  LayoutData<RefCountedPtr<EBStencil> >  m_averageEBStencil;

private:
    //disallowed for all the usual reasons
  EBMGAverage(const EBMGAverage& ebcin)
  {
    MayDay::Error("ebmgi 2 invalid operator");
  }
  void operator=(const EBMGAverage& fabin)
  {
    MayDay::Error("ebmgi 3 invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
