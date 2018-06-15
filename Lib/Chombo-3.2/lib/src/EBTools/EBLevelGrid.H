#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves Nov 14, 2006

#ifndef _EBLEVELGRID_H_
#define _EBLEVELGRID_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

///
/**
   A structure to hold the various components of an EB calculation grid.
  */
class EBLevelGrid
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBLevelGrid();

  ///
  ~EBLevelGrid();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().
  */
  EBLevelGrid(const DisjointBoxLayout& a_dbl,
              const ProblemDomain&     a_domain,
              const int&               a_numGhostEBISL,
              const EBIndexSpace*      a_ebisPtr);

  ///
  EBLevelGrid(const DisjointBoxLayout& a_dbl,
              const EBISLayout&        a_ebisl,
              const ProblemDomain&     a_domain);

  ///
  /**
     Define function that generates EBISL and CFIVS
  */
  void define(const DisjointBoxLayout& a_dbl,
              const ProblemDomain&     a_domain,
              const int&               a_numGhostEBISL,
              const EBIndexSpace*      a_ebisPtr);


  ///
  /**
     only generates cfivs.  set EBISPtr = chombo_ebis::instance
   */
  void define(const DisjointBoxLayout& a_dbl,
              const EBISLayout&        a_ebisl,
              const ProblemDomain&     a_domain);

  ///
  EBLevelGrid(const EBLevelGrid& ebcin);

  ///
  int getGhost() const
  {
    CH_assert(m_isDefined);
    return m_nghost;
  }

  ///
  DisjointBoxLayout getDBL() const
  {
    CH_assert(m_isDefined);
    return m_grids;
  }

  ///
  EBISLayout getEBISL() const
  {
    CH_assert(m_isDefined);
    return m_ebisl;
  }

  ///
  RefCountedPtr<LayoutData<IntVectSet> > getCFIVS() const
  {
    CH_assert(m_isDefined);
    return m_cfivs;
  }

  ///
  bool isDefined() const
  {
    return m_isDefined;
  }

  ///
  const EBIndexSpace* getEBIS() const
  {
    CH_assert(m_isDefined);
    return m_ebisPtr;
  }

  ///
  const IntVectSet& getCoveringIVS() const
  {
    CH_assert(m_isDefined);
    if (!m_isCoveringIVSDefined)
      {
        EBLevelGrid& getAroundConst = (EBLevelGrid&)(*this);
        getAroundConst.defineCoveringIVS();
      }
    return m_coveringIVS;
  }

  ///
  const ProblemDomain& getDomain() const
  {
    CH_assert(m_isDefined);
    return m_domain;
  }

  ///
  /**
     Sets the maximum level of refinement that the EBISLayout
     will have to perform.  Creates and holds new EBISLayouts
     at intermediate levels of refinement.
     Default is one (no refinement done).
  */
  void setMaxRefinementRatio(const int& a_maxRefine)
  {
    CH_assert(m_isDefined);
    m_ebisl.setMaxRefinementRatio(a_maxRefine,m_ebisPtr);
  }

  ///
  /**
     Sets the maximum level of coarsening that the
     will have to perform.  Creates and holds new EBISLayouts
     at intermediate levels of refinement.
     Default is one (no coarsening done).
  */
  void setMaxCoarseningRatio(const int&                a_maxCoarsen,
                             const EBIndexSpace* const a_ebisPtr)
  {
    CH_assert(m_isDefined);
    m_ebisl.setMaxCoarseningRatio(a_maxCoarsen,a_ebisPtr);
  }

  int getMaxCoarseningRatio() const
  {
    CH_assert(m_isDefined);
    return m_ebisl.getMaxCoarseningRatio();
  }

  int
  getMaxRefinementRatio() const
  {
    CH_assert(m_isDefined);
    return m_ebisl.getMaxRefinementRatio();
  }
  ///
  /*
    returns m_grids.coarsenable(a_ref);
  */
  bool
  coarsenable(const int& a_ref) const;

  ///
  /**
     Coarsen eblgFine by a_ref
   */
  friend void coarsen(EBLevelGrid&       a_eblgCoar,
                      const EBLevelGrid& a_eblgFine,
                      const int&         a_ref);


  friend void refine(EBLevelGrid&       a_eblgCoar,
                     const EBLevelGrid& a_eblgFine,
                     const int&         a_ref);

protected:
  void defineCoveringIVS();
  void setDefaultValues();
  bool m_isDefined, m_isCoveringIVSDefined;

  DisjointBoxLayout                      m_grids;
  EBISLayout                             m_ebisl;
  ProblemDomain                          m_domain;
  RefCountedPtr<LayoutData<IntVectSet> > m_cfivs;
  IntVectSet                             m_coveringIVS;
  const EBIndexSpace*                    m_ebisPtr;
  int                                    m_nghost;

private:
};

#include "NamespaceFooter.H"
#endif
