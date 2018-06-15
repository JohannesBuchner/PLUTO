#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCFCOPY_H_
#define _EBCFCOPY_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "Interval.H"

#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBIndexSpace.H"
#include "EBStencil.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

/// Copy coarse values from coarse cells to the fine cells covering them.
/**
    Copy coarse values from coarse cells to the fine cells covering them.
    This was adapted from EBMGInterp to live in EBTools instead of EBAMRTools
    because it's used by EBIndexSpace.
*/
class EBCFCopy
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBCFCopy();

  ///
  ~EBCFCopy();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  EBCFCopy(const DisjointBoxLayout & a_dblFine,
           const DisjointBoxLayout & a_dblCoar,
           const EBISLayout        & a_ebislFine,
           const EBISLayout        & a_ebislCoar,
           const ProblemDomain     & a_domainCoar,
           const int               & a_nref,
           const int               & a_nvar,
           const EBIndexSpace      * a_ebisPtr,
           const IntVect           & a_ghostCellsPhi,
           const bool              & a_layoutChanged = true);


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
  void define(const DisjointBoxLayout & a_dblFine,
              const DisjointBoxLayout & a_dblCoar,
              const EBISLayout        & a_ebislFine,
              const EBISLayout        & a_ebislCoar,
              const ProblemDomain     & a_domainCoar,
              const int               & a_nref,
              const int               & a_nvar,
              const EBIndexSpace      * a_ebisPtr,
              const IntVect           & a_ghostCellsPhi,
              const bool              & a_layoutChanged = true);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Copy coarse data to the fine data covering them (agglomeration allowed).
     fine data = coarse data
  */
  void copy(LevelData<EBCellFAB>       & a_fineData,
            const LevelData<EBCellFAB> & a_coarseData,
            const Interval             & a_variables);

  ///
  /**
     Copy coarse data to the fine data covering them (no agglomeration allowed).
     There are no copyTo's here.
     fine data = coarse data
  */
  void copySameLayout(LevelData<EBCellFAB>       & a_fineData,
                      const LevelData<EBCellFAB> & a_coarseData,
                      const Interval             & a_variables);

protected:
  void copyFAB(EBCellFAB       & a_refCoar,
               const Box       & a_coarBox,
               const EBCellFAB & a_fine,
               const DataIndex & a_datInd,
               const Interval  & a_variables) const;

  void setDefaultValues();

  void defineStencils();

  bool m_isDefined;

  IntVect m_ghost;

  DisjointBoxLayout m_coarGrids;
  DisjointBoxLayout m_fineGrids;

  ProblemDomain m_coarDomain;
  ProblemDomain m_fineDomain;

  EBISLayout m_coarEBISL;
  EBISLayout m_fineEBISL;

  int m_refRat;
  int m_nComp;

  bool m_layoutChanged;
  bool m_coarsenable;

  // LevelData<EBCellFAB> m_refinedCoarseData;
  // if coarsenable, these are the coarsened fine
  // otherwise, refined coarse.
  // actual buffer is no longer member data
  // not defined if !m_layoutChanged
  DisjointBoxLayout m_buffGrids;
  EBISLayout        m_buffEBISL;

  // the following copiers only get defined where needed.
  Copier m_copierFtoRC; // this goes from m_fineGrids to m_refinedCoarseGrids
  Copier m_copierRCtoF; // this goes from m_refinedCoarseGrids to m_fineGrids

  // stencils for prolongation
  LayoutData<RefCountedPtr<EBStencil> >  m_copyEBStencil;

private:
  // disallowed for all the usual reasons
  EBCFCopy(const EBCFCopy& ebcin)
  {
    MayDay::Error("EBCFCopy 2 invalid operator");
  }
  void operator=(const EBCFCopy& fabin)
  {
    MayDay::Error("EBCFCopy 3 invalid operator");
  }

};

#include "NamespaceFooter.H"

#endif
