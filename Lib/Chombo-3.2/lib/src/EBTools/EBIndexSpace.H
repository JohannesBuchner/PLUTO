#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#ifndef _EBINDEXSPACE_H_
#define _EBINDEXSPACE_H_

#include "REAL.H"
#include "IntVect.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "CH_HDF5.H"

#include "VolIndex.H"
#include "FaceIndex.H"
#include "EBISBox.H"
#include "EBISLayout.H"
#include "GeometryService.H"
#include "EBISLevel.H"

#include "NamespaceHeader.H"

///
/**
   EBIndexSpace represents the geometric information
   of the domain.    It should be generated through
   the Chombo_EBIS class.  This follows the singleton
   pattern.
*/
class EBIndexSpace
{
public:

  EBIndexSpace();

  ///
  /**
      If a_ncellMax is set, that is the max width of
      an internal grid.  Otherwise use defaults
      of (16 in 3D, 64 in 2d)
   */
  void
  define(const ProblemDomain   & a_domain,
         const RealVect        & a_origin,
         const Real            & a_dx,
         const GeometryService & a_geoserver,
         int                     a_nCellMax = -1,
         int                     a_maxCoarsenings = -1);

#ifdef CH_USE_HDF5
  /// reads in all levels from finest level down
  void
  define(HDF5Handle & a_handle, ProblemDomain a_finestLevel);

  /// defines from one level and then does coarsening
  void
  define(HDF5Handle & a_handle,
         int          a_maxCoarsenings = -1);
#endif

  void
  define(EBISLevel * a_level0,
         int         a_nCellMax = -1,
         int         a_maxCoarsenings = -1);

  void
  define(const ProblemDomain                        & a_entireDomain,
         const RealVect                             & a_origin,
         const Real                                 & a_dx,
         const Vector<RefCountedPtr<EBIndexSpace> > & a_patches,
         const Vector<IntVect>                      & a_offsets,
         int                                          a_maxCoarsenings = -1);

  EBISLevel* buildFirstLevel(const ProblemDomain& a_domain,
                             const RealVect& a_origin,
                             const Real& a_dx,
                             const GeometryService& a_geoserver,
                             int a_nCellMax,
                             int a_maxCoarsenings,
                             bool a_fixRegularNextToMultiValued = true);

  EBISLevel* buildNextLevel(const GeometryService & a_geoserver,
                            bool                    a_fixRegularNextToMultiValued = true);

  ///
  ~EBIndexSpace();

  ///
  /**
    Set a flag that indicates the data defining the EB is distributed.
    This will signal EBISLevel to hand grid generation off to the
    GeometryService.  Function returns previous value of m_distributedData.
    Call this before define() if this functionality is desired.
  */
  bool setDistributedData()
  {
    bool was = m_distributedData;
    m_distributedData = true;
    return was;
  }

  ///
  int numLevels() const;

  ///
  /**
     Get the number of vofs over the entire domain.
     This is blocking as a broadcast and gather are required.
   */
  long long numVoFs(const ProblemDomain& a_domain) const;

  ///
  /**
     return level index of domain.
     return -1 if a_domain does not
     correspond to any refinement of EBIS.
   */
  int getLevel(const ProblemDomain& a_domain) const;

  RealVect getOrigin() const
  {
    return m_ebisLevel[0]->getOrigin();
  }

  ///
  /**
     returns the problem domain box at a level of refinement
  */
  const ProblemDomain& getBox(int level) const
  {
    return m_domainLevel[level];
  }

  ///
  void fillEBISLayout(EBISLayout&              a_ebis,
                      const DisjointBoxLayout& a_grids,
                      const ProblemDomain&     a_domain,
                      const int&               a_nghost) const;

  ///
  /**
     Return true if the define function has been called.
   */
  bool isDefined() const;

  ///
  void clear();

  int getNCellMax() const;
  ///
  /**
     This makes Chombo_EBIS be the only function
     that can call the private constructor.
   */
  friend class Chombo_EBIS;

#ifdef CH_USE_HDF5
  ///
  /**
     Writes out finest level as a default.   Writes at the level you are domain
     you want if a_outputLevel is defined.
   */
  void
  write(HDF5Handle&   a_handle,
        ProblemDomain a_outputLevel = ProblemDomain()) const;

  void
  writeInfo(HDF5Handle& handle) const ;

  ///
  /**
     Outputs ALL levels of the EBIS.
     this is the one to use if you want to do all definition through
     the file (no coarsening).
   */
  void writeAllLevels(HDF5Handle&  a_handle) const;

  ///
  /**
     Define entirely from the input file.   a_finestLevel can be a coarsening
     of the finest level in the file.   No coarsening algorithm is done.
   */
  void readInAllLevels(HDF5Handle & a_handle,  ProblemDomain a_finestLevel);

#endif

  DisjointBoxLayout getIrregGrids(int level) const
  {
    return m_ebisLevel[level]->getIrregGrids(m_domainLevel[level]);
  }

  DisjointBoxLayout getIrregGrids(const ProblemDomain& a_domain) const
  {
    return m_ebisLevel[getLevel(a_domain)]->getIrregGrids(a_domain);
  }

  DisjointBoxLayout getFlowGrids(int level) const
  {
    return m_ebisLevel[level]->getFlowGrids(m_domainLevel[level]);
  }

  DisjointBoxLayout getFlowGrids(const ProblemDomain& a_domain) const
  {
    return m_ebisLevel[getLevel(a_domain)]->getFlowGrids(a_domain);
  }

  DisjointBoxLayout getCoveredGrids(int level) const
  {
    return m_ebisLevel[level]->getCoveredGrids(m_domainLevel[level]);
  }

  DisjointBoxLayout getCoveredGrids(const ProblemDomain& a_domain) const
  {
    return m_ebisLevel[getLevel(a_domain)]->getCoveredGrids(a_domain);
  }

  DisjointBoxLayout levelGrids(int level) const
  {
    return m_ebisLevel[level]->m_grids;
  }

  Real dx(int level) const
  {
    return m_ebisLevel[level]->m_dx;
  }

  IntVectSet irregCells(int depth=0) const
  {
    return m_ebisLevel[depth]->irregCells();
  }

  void setCellMax(int max)
  {
    m_nCellMax = max;
  }

  void resetLevels(int nLevel);

  Vector<RefCountedPtr<EBIndexSpace> > connectedComponents();

  RefCountedPtr<EBIndexSpace> biggestConnectedComponent(int & a_numComponents);

  static bool s_MFSingleBox;

private:
  Vector<RefCountedPtr<EBIndexSpace> > findConnectedComponents(int        & a_numComponents,
                                                               const bool & a_onlyBiggest);

  bool setAllConnectedVoFs(Real&               a_totalVolFrac,
                           EBCellFAB&          a_curEBCellFAB,
                           const EBGraph&      a_curEBGraph,
                           const EBISBox&      a_curEBISBox,
                           const VolIndex&     a_curVoF,
                           const VolIndex&     a_lastVoF,
                           const unsigned int& a_curNum);

  void resetAllConnectedVoFs(EBCellFAB&          a_curEBCellFAB,
                             const EBGraph&      a_curEBGraph,
                             const EBISBox&      a_curEBISBox,
                             const VolIndex&     a_curVoF,
                             const VolIndex&     a_lastVoF);

  //this is private to force the singleton thing.
  //EBIndexSpace();

  int  m_nCellMax;
  bool m_isDefined;

  bool m_distributedData;

  Vector<EBISLevel*> m_ebisLevel;
  Vector<ProblemDomain> m_domainLevel;

  int m_nlevels;

  static Real s_tolerance;
  static bool s_verbose;

  //disallowed for performance reasons
  void operator=(const EBIndexSpace& ebiin)
  {
    MayDay::Error("no assignment of EBIndexSpace allowed");
  }

  EBIndexSpace(const EBIndexSpace& ebiin)
  {
    MayDay::Error("no copy construction of EBIndexSpace allowed");
  }
};

#include "NamespaceFooter.H"

#endif
