#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCOARSEN_H_
#define _EBCOARSEN_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "EBLevelGrid.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "VoFIterator.H"
#include "NamespaceHeader.H"
class EBIndexSpace;
class VoFStencil;

/// replaces coarse level data with a coarsening of fine level data.
/**
   This class replaces data at a coarse level of refinement with a
   coarsening of data at a finer level of refinement, in areas where fine
   data is present.  Coarse level data is not modified where fine
   level data is not present.
  */
class EBCoarsen
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBCoarsen();

  ///
  ~EBCoarsen();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     eblgFine, eblgCoar: The fine and coarse layouts
     of the data.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  EBCoarsen(const EBLevelGrid& eblgFine,
            const EBLevelGrid& eblgCoar,
            const int& nref,
            const int& nvar);

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     eblgFine, eblgCoar: The fine and coarse layouts
     of the data.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  void define(const EBLevelGrid& eblgFine,
              const EBLevelGrid& eblgCoar,
              const int& nref,
              const int& nvar);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Replaces a_coarse_data with the coarsening of a_fine_data within the
     coarsening of the domain of the fine level.  This is different from
     averaging in that we are trying to approximate the coarse value at
     the cell center. (averaging approximates the average over the cell)
     Elsewhere, the
     coarse data is unchanged.  It is an error to call if not
     this->isDefined(). The domain of a_fine_data should be
     the same as the fine domain specified in the most recent call to
     define().  It is expected that the coarse and fine level's
     domains are properly nested.  Both a_coarse_data and a_fine_data
     should have the same number of components specified in the most
     recent call to define().
     a_fine_data needs pre-filled ghost cells of size 1 ONLY IF refRat<4

     {\bf Arguments:}\\
     coarData: The data over the coarse layout.\\
     fineData: The data over the fine layout.  \\
     Fine and coarse data must
     have the same number of variables.\\
     variables:  The variables to coarsen.  Those not
     in this range will be left alone.   This range of variables
     must be in both the coarse and fine data.
  */
  void
  coarsenFine(LevelData<EBCellFAB>&       a_coarData,
              const LevelData<EBCellFAB>& a_fineData,
              const Interval&             a_variables);
protected:
  void defineStencil(const LayoutData<IntVectSet>&  a_cfivs);

  void coarsenIrreg(EBCellFAB&       a_coar,
                    const EBCellFAB& a_fine,
                    const DataIndex& a_dit,
                    const Interval&  a_variables);

  void checkStencil(VoFStencil& a_stencil,
                    const Real& a_sum) const;

  bool getVofsSideDir(Vector<VolIndex>&     a_vofList,
                      const VolIndex&       a_vof,
                      const IntVectSet&     a_cfivs,
                      const EBISBox&        a_ebisBox,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const int&            a_stepSize) const;

  void addWeight(VoFStencil&             a_stencil,
                 const Vector<VolIndex>& a_vofList,
                 const Real&             a_weight,
                 const bool&             a_hasVof) const;

  void getCoarsenVoFStencil(VoFStencil&       a_stencil,
                            const EBISBox&    a_ebisBoxCoar,
                            const EBISBox&    a_ebisBoxFine,
                            const VolIndex&   a_vofCoar,
                            const DataIndex&  a_datInd,
                            const IntVectSet& a_cfivs);

  void coarsenFAB(EBCellFAB&       a_coar,
                  const EBCellFAB& a_fine,
                  const DataIndex& a_datInd,
                  const Interval&  a_variables);

  void
  setDefaultValues();

  bool m_isDefined;

  DisjointBoxLayout m_gridsCoar;
  DisjointBoxLayout m_gridsFine;
  DisjointBoxLayout m_coarsenedFineGrids;
  ProblemDomain     m_domainCoar;
  ProblemDomain     m_domainFine;

  EBISLayout m_ebislCoar;
  EBISLayout m_ebislFine;

  EBISLayout m_coarsenedFineEBISL;
  RefCountedPtr<LayoutData<IntVectSet> > m_cfivsPtr;

  int m_refRat;
  int m_nComp;
  LevelData<EBCellFAB> m_coarsenedFineData;

  LayoutData<BaseIVFAB<VoFStencil> > m_coarsenStencil;
  LayoutData<VoFIterator>           m_vofIt;

private:
    //disallowed for all the usual reasons
  EBCoarsen(const EBCoarsen& ebcin)
  {
    MayDay::Error("ebca 2 invalid operator");
  }
  void operator=(const EBCoarsen& fabin)
  {
    MayDay::Error("ebca 3 invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
