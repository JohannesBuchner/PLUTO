#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves  11-07-2001

#ifndef _EBFINETOCOARREDIST_H_
#define _EBFINETOCOARREDIST_H_
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

///
/**
   A EBFineToCoarRedist  handles all the data choreography
   for fine to coarse redistribution.
 */
class EBFineToCoarRedist
{
public:

  friend class EBFluxRegister;

  ///
  /**
     Modify the weights in the stencil by multiplying by
     the inputs in a normalized way.
     If you want mass
     weighting, send in the density  on the coarse layout.
  */
  void
  resetWeights(const LevelData<EBCellFAB>& a_modifierCoar,
               const int& a_ivar);

  ///
  /**
     Default constructor.  Leaves object undefined.
  */
  EBFineToCoarRedist();

  ///
  ~EBFineToCoarRedist();

  /// uglier but general define fcn
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
              int a_redistRad = 1,
              const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  ///
  /**
     potentially faster define function, especially with large numbers of
     boxes.
   */
  void define(const EBLevelGrid&  a_eblgFine,
              const EBLevelGrid&  a_eblgCoar,
              const int&          a_nref,
              const int&          a_nvar,
              const int&          a_redistRad);


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
     coarse level.
  */
  void increment(const BaseIVFAB<Real>& a_fineMass,
                 const DataIndex& a_fineDataIndex,
                 const Interval&  a_variables);

  ///
  /**
     Redistribute the data contained in the internal buffers.
  */
  void redistribute(LevelData<EBCellFAB>& a_coarSolution,
                    const Interval& a_variables);

  ///
  bool isDefined() const;

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
  //RefCoar is short for refinedCoarse
  DisjointBoxLayout    m_gridsRefCoar;

  //buffers on fine layout
  LevelData<BaseIVFAB<Real> > m_regsFine;
  //buffers on refineed coarse layout
  LevelData<BaseIVFAB<Real> > m_regsRefCoar;

  //need both of these to facilitate mass-weighted
  //redistribution
  LayoutData<BaseIVFAB<VoFStencil> > m_stenRefCoar;
  LayoutData<BaseIVFAB<VoFStencil> > m_volumeStenc;
  LayoutData<BaseIVFAB<VoFStencil> > m_standardStenc;

  LevelData<EBCellFAB> m_densityCoar;
  //location of mass sources on fine layout
  LayoutData<IntVectSet> m_setsFine;
  //location of mass sources on refined coarse layout
  LayoutData<IntVectSet> m_setsRefCoar;

  //ebisl of input fine grid
  EBISLayout m_ebislFine;
  //ebisl of input coar grid
  EBISLayout m_ebislCoar;
  //ebisl of refined coarse layout
  EBISLayout m_ebislRefCoar;

private:
  //For all the usual reasons,
  //there is no copy constructor for this class.
  //Neither is there an operator= for this class.
  void operator= (const EBFineToCoarRedist&)
  {
    MayDay::Error("invalid operator");
  }
  EBFineToCoarRedist(const EBFineToCoarRedist&)
  {
    MayDay::Error("invalid operator");
  }
};
#include "NamespaceFooter.H"
#endif
