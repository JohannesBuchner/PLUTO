#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELFLUXREGISTEREDGE_H_
#define _LEVELFLUXREGISTEREDGE_H_

#include "REAL.H"
#include "Vector.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "IntVectSet.H"
#include "CFStencil.H"
#include "LoHiSide.H"
#include "LevelData.H"
#include "LayoutData.H"

#include "NamespaceHeader.H"

///LevelFluxRegisterEdge -- Flux register for constrained-transport applications
/**
   A LevelFluxRegisterEdge manages the coarse-fine fixup for a face-centered
   field which is defined as a curl of an edge-centered field, performing
   the "reflux-curl" operation described in Balsara(2001) to preserve the
   divergence-free magnetic field in the presence of coarse-fine interfaces.

   This class performs a reflux-curl of edge-centered fluxes to correct a
   face-centered field. This is in contrast to the regular LevelFluxRegister
   class, which performs a reflux-divergence of face-centered fluxes to
   correct a cell-centered field.
 */
class LevelFluxRegisterEdge
{
public:

  ///
  /**
     Default constructor.  Creates an uninitialized LevelFluxRegisterEdge.
  */
  LevelFluxRegisterEdge();

  ///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.
  */
  LevelFluxRegisterEdge(
                    const DisjointBoxLayout& a_dbl,
                    const DisjointBoxLayout& a_dblCoarse,
                    const Box& a_dProblem,
                    int a_nRefine,
                    int a_nComp);

  ///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.
  */
  LevelFluxRegisterEdge(
                    const DisjointBoxLayout& a_dbl,
                    const DisjointBoxLayout& a_dblCoarse,
                    const ProblemDomain& a_dProblem,
                    int a_nRefine,
                    int a_nComp);

  ///
  ~LevelFluxRegisterEdge();

  ///
  /**
     Full define function.  Creates a levels worth of flux registers.
     The values in the flux registers are still undefined, however.
     To zero the fluxregisters, you must call setToZero().
  */
  void define(
              const DisjointBoxLayout& a_dbl,
              const DisjointBoxLayout& a_dblCoarse,
              const Box& a_dProblem,
              int a_nRefine,
              int a_nComp);


  ///
  /**
     Full define function.  Creates a levels worth of flux registers.
     The values in the flux registers are still undefined, however.
     To zero the fluxregisters, you must call setToZero().
  */
  void define(
              const DisjointBoxLayout& a_dbl,
              const DisjointBoxLayout& a_dblCoarse,
              const ProblemDomain& a_dProblem,
              int a_nRefine,
              int a_nComp);

///
/**
   Modifies this LevelFluxRegisterEdge so that it is returned to the
   uninitialized state.  User must now call the full define() before
   using it.

*/
  void
  undefine();

  ///
  /**
      Initialize values of registers  to zero.
  */
  void setToZero();

  ///
  /**
     increments the register with data from coarseFlux, multiplied by scale.
     coarseFlux must contain the edge-centered (in 3d, node centered in
     2d) coarse fluxes in the dir direction
     for the grid   m_coarseLayout[coarseDataIndex].
     By convention, only the low side flux is used to
     avoid double-counting at coarse-fine interfaces.
     This operation is local.
  */
  void incrementCoarse(
                       FArrayBox& a_coarseFlux,
                       Real a_scale,
                       const DataIndex& a_coarseDataIndex,
                       const Interval& a_srcInterval,
                       const Interval& a_dstInterval);

  ///
  /**
     increments the register with data from fineFlux (which is edge-centered
     in 3d, node-centered in 2d), multiplied by scale, for all coarse-fine
     face directions associated with the grid box m_fineLayout[fineDataIndex]
     This operation is local.
  */
  void incrementFine(
                     FArrayBox& a_fineFlux,
                     Real a_scale,
                     const DataIndex& a_fineDataIndex,
                     const Interval& a_srcInterval,
                     const Interval& a_dstInterval);


  ///
  /**
     increments the register with data from fineFlux (which is edge-centered
     in 3d, node-centered in 2d), multiplied by scale.
     a_dir is the normal of the coarse-fine interface, and a_sd determines
     whether we're looking at the high-side or the low-side
     for the grid box m_fineLayout[fineDataIndex]
     This operation is local.
  */
  void incrementFine(
                     FArrayBox& a_fineFlux,
                     Real a_scale,
                     const DataIndex& a_fineDataIndex,
                     const Interval& a_srcInterval,
                     const Interval& a_dstInterval,
                     int a_dir,
                     Side::LoHiSide a_sd);

  ///
  /**
     increments uCoarse with the reflux "CURL" of the
     contents of the flux register.
     Note that there is no srccomp etc here.
     this is done for all components so uCoarse has
     to have the same number of components as input nComp.
     This operation is global and blocking.
  */
  void refluxCurl(
                  LevelData<FluxBox>& a_uCoarse,
                  Real a_scale);

  ///
  /**
  has full define function been called?  return true if so.
  */
  bool isDefined() const;

  //write the contents of all the registers
  void dump();

  //spew the contents of all the low coar registers
  void dumpLoCoar(int idir);

  //spew the contents of all the high coar registers
  void dumpHiCoar(int idir);

  //spew the contents of all the low fine registers
  void dumpLoFine(int idir);

  //spew the contents of all the high fine registers
  void dumpHiFine(int idir);
protected:

  static int index(int dir, Side::LoHiSide side);

  static int getRegComp(const int& faceDir, const int& edgeDir);

  //internal use only
  void setDefaultValues();

  // both fabCoarse and fabFine have FluxBoxes that are
  // defined in the coarse index space.  fabCoarse
  // can ONLY be indexed into by the DataIndex stored in
  // coarToCoarMap.  fabFine is indexed by the same indexing
  // used in the fine DisjointBoxLayout.
  LevelData<EdgeDataBox> m_regCoarse;
  LevelData<FluxBox> m_fabFine[SpaceDim*2];

  // in an attempt to speed things up in the periodic case especially,
  // precompute copiers used to copy from fabFine to coarse-level-indexed
  // fabs
  Vector<Copier> m_crseCopiers;

  /**: reflux locations -- need SpaceDim of these to account
     for the different face directions, which explains the extra
     Vector */
  LayoutData<Vector<Vector<IntVectSet> > > m_refluxLocations[SpaceDim*2];
  LayoutData<Vector<DataIndex>  > m_coarToCoarMap[SpaceDim*2];

  bool m_isDefined;

  ///number of components in register
  int m_nComp;

  ///refinement ratio between levels
  int m_nRefine;

  /// domain at the coarse grid resolution
  ProblemDomain m_domainCoarse;

private:

  ///there is no operator= for this class
  void operator= (const LevelFluxRegisterEdge&);
  ///there is no copy constructor for this class
  LevelFluxRegisterEdge(const LevelFluxRegisterEdge&);
};

#include "NamespaceFooter.H"

#endif
