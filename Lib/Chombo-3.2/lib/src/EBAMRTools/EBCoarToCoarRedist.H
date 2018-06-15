#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves mon nov 19 2001

#ifndef _EBCOARTOCOARREDIST_H_
#define _EBCOARTOCOARREDIST_H_

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
   A EBCoarToCoarRedist  handles all the data choreography
   for coarse to coarse redistribution ("re-redistribution").
 */
class EBCoarToCoarRedist
{
public:

  friend class EBFluxRegister;

  ///
  /**
     Default constructor.  Leaves object undefined.
  */
  EBCoarToCoarRedist();

  ///
  ~EBCoarToCoarRedist();

  ///
  /**
     Full define function.   \\
     dblFine, dblCoar: The fine and coarse layouts
     of the data. \\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels.\\
     nvar: The number of variables contained in the data
     at each VoF. \\
     coarStencil: The redistribution stencil for this level.
  */
  void define(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout& a_ebislFine,
              const EBISLayout& a_ebislCoar,
              const Box& a_domainCoar,
              const int& a_nref,
              const int& a_nvar,
              int a_redistRad = 1);

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
     the class figures out what actually gets re-redistributed.
  */
  void increment(const BaseIVFAB<Real>& a_coarMass,
                 const DataIndex& a_fineDataIndex,
                 const Interval&  a_variables);

  ///
  /**
   */
  void
  resetWeights(const LevelData<EBCellFAB>& a_modifier,
               const int& a_ivar);

  ///
  /**
     Reredistribute the data contained in the internal buffers.
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
  int  m_redistRad;
  int  m_nComp;
  int  m_refRat;
  Box  m_domainCoar;

  //input fine grid
  DisjointBoxLayout    m_gridsFine;
  //input coarse grid
  DisjointBoxLayout    m_gridsCoar;

  //buffers on refineed coarse layout
  LevelData<BaseIVFAB<Real> > m_regsCoar;

  //location of mass sources on refined coarse layout
  LayoutData<IntVectSet> m_setsCoar;

  //ebisl of input fine grid
  EBISLayout m_ebislFine;
  //ebisl of input coar grid
  EBISLayout m_ebislCoar;

  LayoutData<BaseIVFAB<VoFStencil> > m_stenCoar;
  //need both of these to facilitate mass-weighted
  //redistribution
  LayoutData<BaseIVFAB<VoFStencil> > m_volumeStenc;
  LayoutData<BaseIVFAB<VoFStencil> > m_standardStenc;
  LevelData<EBCellFAB> m_densityCoar;

private:
  //For all the usual reasons,
  //there is no copy constructor for this class.
  //Neither is there an operator= for this class.
  void operator= (const EBCoarToCoarRedist&)
  {
    MayDay::Error("invalid operator");
  }
  EBCoarToCoarRedist(const EBCoarToCoarRedist&)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
