#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves  11-02-2001

#ifndef _EBCOARTOFINEREDIST_H_
#define _EBCOARTOFINEREDIST_H_
#include "REAL.H"
#include "Vector.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "IntVectSet.H"
#include "CFStencil.H"
#include "LoHiSide.H"
#include "LevelData.H"
#include "LayoutData.H"
#include "EBLevelGrid.H"
#include "RedistStencil.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

///   Coarse to fine redistribution class.
/**
   Coarse to fine redistribution class.
 */
class EBCoarToFineRedist
{
public:

  friend class EBFluxRegister;

  ///
  /**
     Default constructor.  Leaves object undefined.
  */
  EBCoarToFineRedist();

  ///
  ~EBCoarToFineRedist();

  ///
  /**
      Initialize values of registers  to zero.
  */
  void setToZero();

  ///
  /**
     Increments the register with data from coarseMass.
     This is the full redistribution mass.  Internally
     the class figures out what actually goes to the
     fine level.
  */
  void increment(const BaseIVFAB<Real>& a_coarseMass,
                 const DataIndex& a_coarseDataIndex,
                 const Interval&  a_variables);

  ///
  /**
     Redistribute the data contained in the internal buffers.
  */
  void redistribute(LevelData<EBCellFAB>& a_fineSolution,
                    const Interval& a_variables);

  ///
  bool isDefined() const;

  /// ugly but general define fcn
  /**
      Full define function.
      Define the stencils with volume weights.
      If you want mass weights or whatever,
      use reset weights.
  */
  void define(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout& a_ebislFine,
              const EBISLayout& a_ebislCoar,
              const Box& a_domainCoar,
              const int& a_nref,
              const int& a_nvar,
              int redistRad,
              const EBIndexSpace* ebisPtr);

  ///
  /**
     Interface that uses EBLevelGrid potentially can be
     much faster when the number of boxes gets large.
   */
  void define(const EBLevelGrid& a_eblgFine,
              const EBLevelGrid& a_eblgCoar,
              const int& a_nref,
              const int& a_nvar,
              const int& a_redistRad);


  ///
  /**
  */
  void
  resetWeights(const LevelData<EBCellFAB>& a_modifierCoar,
               const int& a_ivar);
protected:

  //internal use only
  void setDefaultValues();
  void defineDataHolders();

  bool m_isDefined;
  int m_redistRad;
  int m_nComp;
  int m_refRat;
  Box m_domainCoar;

  //input fine grid
  DisjointBoxLayout    m_gridsFine;
  //input coarse grid
  DisjointBoxLayout    m_gridsCoar;
  //input coarse grid
  //CedFine is short for coarsenedfine
  DisjointBoxLayout    m_gridsCedFine;

  LevelData<BaseIVFAB<Real> > m_regsCoar;
  LevelData<BaseIVFAB<Real> > m_regsCedFine;

  LevelData<EBCellFAB> m_densityCedFine;

  LayoutData<BaseIVFAB<VoFStencil> > m_stenCedFine;
  //need both of these to accomodate reweighting the stencil
  LayoutData<BaseIVFAB<VoFStencil> > m_volumeStenc;
  LayoutData<BaseIVFAB<VoFStencil> > m_standardStenc;

  LayoutData<IntVectSet> m_setsCedFine;
  LayoutData<IntVectSet> m_setsCoar;

  //ebisl of input fine grid
  EBISLayout m_ebislFine;
  //ebisl of input coar grid
  EBISLayout m_ebislCoar;
  //ebisl of coarsened fine layout
  EBISLayout m_ebislCedFine;

private:
  //For all the usual reasons,
  //there is no copy constructor for this class.
  //Neither is there an operator= for this class.
  void operator= (const EBCoarToFineRedist&)
  {
    MayDay::Error("invalid operator");
  }
  EBCoarToFineRedist(const EBCoarToFineRedist&)
  {
    MayDay::Error("invalid operator");
  }
};
#include "NamespaceFooter.H"
#endif
