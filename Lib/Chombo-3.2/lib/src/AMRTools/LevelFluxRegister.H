#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELFLUXREGISTER_H_
#define _LEVELFLUXREGISTER_H_

#include "REAL.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "IntVectSet.H"
#include "LoHiSide.H"
#include "LevelData.H"
#include "LayoutData.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///LevelFluxRegister-A class to encapsulate a levels worth of flux registers.
/**
   A LevelFluxRegister  handles all the data choreography
   necessary to create a levels worth of flux registers.
*/
class LevelFluxRegister
{
protected:

  /// Define status
  enum FluxRegDefineStatus
  {
    FluxRegUndefined      =        0,
    FluxRegDefined        = (1 << 0),      // Indicates define() was called
    FluxRegFineDefined    = (1 << 1),      // m_fineFlux was defined
    FluxRegCoarseDefined  = (1 << 2),      // m_coarFlux was defined
    FluxRegAllDefined     = (1 << 3) - 1   // Everything
  };

public:
  static bool s_verbose;
  ///
  /**
     Default constructor.  Creates an uninitialized levelfluxregister.
  */
  LevelFluxRegister();

  ///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.  a_dProblem is fine-level
     domain.

     If a_scaleFineFluxes is true (default then all fine fluxes are
     scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
  */
  LevelFluxRegister(const DisjointBoxLayout& a_dbl,
                    const DisjointBoxLayout& a_dblCoarse,
                    const ProblemDomain&     a_dProblem,
                    int                      a_nRefine,
                    int                      a_nComp,
                    bool                     a_scaleFineFluxes = true,
                    int                      a_unrefinedDirection=-1);

  /// This is an old style constructor - do not use if at all possible.
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.  a_dProblem is fine-level
     domain.

     If a_scaleFineFluxes is true (default then all fine fluxes are
     scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
  */
  LevelFluxRegister(const DisjointBoxLayout& a_dbl,
                    const DisjointBoxLayout& a_dblCoarse,
                    const Box&               a_dProblem,
                    int                      a_nRefine,
                    int                      a_nComp,
                    bool                     a_scaleFineFluxes = true,
                    int                      a_unrefinedDirection = -1);

  ///
  virtual ~LevelFluxRegister();

  ///
  /**
     Full define function.  Creates a levels worth of flux registers.
     The values in the flux registers are still undefined, however.
     To zero the fluxregisters, you must call setToZero().  a_dProblem
     is fine-level domain. This method may be overridden by a subclass.

     If a_scaleFineFluxes is true then all fine fluxes are
     scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
  */
  virtual void define(const DisjointBoxLayout& a_dbl,
                      const DisjointBoxLayout& a_dblCoarse,
                      const ProblemDomain&     a_dProblem,
                      int                      a_nRefine,
                      int                      a_nComp,
                      bool                     a_scaleFineFluxes,
                      int                      a_unrefinedDirection = -1);

  ///
  /**
     This define function calls the above define(), setting a_scaleFineFluxes to true.
  */
  virtual void define(const DisjointBoxLayout& a_dbl,
              const DisjointBoxLayout& a_dblCoarse,
              const ProblemDomain&     a_dProblem,
              int                      a_nRefine,
              int                      a_nComp);

  /// This is an old style define - do not use if at all possible.
  /**
     Full define function.  Creates a levels worth of flux registers.
     The values in the flux registers are still undefined, however.
     To zero the fluxregisters, you must call setToZero().  a_dProblem
     is fine-level domain.

     If a_scaleFineFluxes is true (default then all fine fluxes are
     scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
  */
  virtual void define(const DisjointBoxLayout& a_dbl,
              const DisjointBoxLayout& a_dblCoarse,
              const Box&               a_dProblem,
              int                      a_nRefine,
              int                      a_nComp,
              bool                     a_scaleFineFluxes = true);

  ///
  /**
     Modifies this LevelFluxRegister so that it is returned to the
     uninitialized state.  User must now call the full define() before
     using it.

  */
  virtual void
  undefine();

  ///
  /**
      Initialize values of registers to zero.
  */
  virtual void setToZero();

  ///
  /**
     increments the register with data from coarseFlux, multiplied by scale.
     coarseFlux must contain the coarse fluxes in the dir direction
     for the grid   m_coarseLayout[coarseDataIndex].  This operation is local.
     Increment coarse does register += -sign(side)*incrScale*coarseFlux
     for both Lo and Hi "side".
  */
  virtual void incrementCoarse(
                       const FArrayBox& a_coarseFlux,
                       Real a_scale,
                       const DataIndex& a_coarseDataIndex,
                       const Interval& a_srcInterval,
                       const Interval& a_dstInterval,
                       int a_dir);

  ///
  /**
     This method implements the register's coarse increment on a single side.
     It may be overridden by a subclass.
  */
  virtual void incrementCoarse(
                               const FArrayBox& a_coarseFlux,
                               Real a_scale,
                               const DataIndex& a_coarseDataIndex,
                               const Interval& a_srcInterval,
                               const Interval& a_dstInterval,
                               int a_dir,
                               Side::LoHiSide a_sd);

  ///
  /**
     increments the register with data from fineFlux, multiplied by scale.
     fineFlux must contain the coarse fluxes in the dir direction for the
     grid m_fineLayout[fineDataIndex].  This operation is local.
     Increment fine does register +=  sign(side)*incrScale*Ave(fineFlux)
     for both Lo and Hi "side".
  */
  virtual void incrementFine(
                     const FArrayBox& a_fineFlux,
                     Real a_scale,
                     const DataIndex& a_fineDataIndex,
                     const Interval& a_srcInterval,
                     const Interval& a_dstInterval,
                     int a_dir);

  ///
  /**
     This method implements the register's fine increment on one side.
     It may be overridden by a subclass.
  */
  virtual void incrementFine(
                             const FArrayBox& a_fineFlux,
                             Real a_scale,
                             const DataIndex& a_fineDataIndex,
                             const Interval& a_srcInterval,
                             const Interval& a_dstInterval,
                             int a_dir,
                             Side::LoHiSide a_sd);

  ///
  /**
     increments uCoarse with the reflux divergence of the
     contents of the flux register.
     Note that there is no srccomp etc here.
     this is done for all components so uCoarse has
     to have the same number of components as input nComp.
     This operation is global and blocking.

     Reflux does cellValue -= refluxScale*registerContents
  */
  virtual void reflux(
              LevelData<FArrayBox>& a_uCoarse,
              Real a_scale);

  ///
  /**
     Interface for multiblock algorithms where, if an interval
     of the components represents a vector, it must be specially
     handled.  If called from a singleblock code, the vector
     interval is just ignored.
  */
  virtual void reflux(
                      LevelData<FArrayBox>& a_uCoarse,
                      const Interval& a_coarseVectorIntv,
                      Real a_scale);

  ///
  /**
     increments uCoarse with the reflux divergence of the
     contents of the flux register. This can be overridden by
     a subclass. This operation is global and blocking.
     Reflux does cellValue -= refluxScale*registerContents
  */
  virtual void reflux(
                      LevelData<FArrayBox>& a_uCoarse,
                      const Interval& a_coarse_interval,
                      const Interval& a_flux_interval,
                      Real a_scale);

  ///same as above with a variable scale multiplied in
  virtual void reflux(
              LevelData<FArrayBox>& a_uCoarse,
              Real a_scale,
              const Interval& a_coarseInterv,
              const Interval& a_fluxInterval,
              const LevelData<FArrayBox>& a_beta);


  bool hasCF(const DataIndex& a_fineDataIndex, Side::LoHiSide) const;
  bool hasCF(const DataIndex& a_coarseIndex) const;

  ///
  /**
     has full define function been called?  return true if so.
  */
  bool isDefined() const;

  bool isAllDefined() const;

  ///
  void poutCoarseRegisters() const;

  ///
  void poutFineRegisters() const;

  LayoutData< Vector<Box> >& getCoarseLocations(int a_idir, Side::LoHiSide a_sd)
  {
    CH_assert(isDefined());
    int index = a_idir + a_sd*CH_SPACEDIM;
    return m_coarseLocations[index];
  }

  Copier& getReverseCopier()
  {
    return m_reverseCopier;
  }

protected:

  // both of these are defined in the coarse index space,
  // but, m_coarFlux uses the coarse grid layout, and
  // m_fineFlux uses the fine grid layout.
  LevelData<FArrayBox> m_coarFlux;
  LevelData<FArrayBox> m_fineFlux;

  //LayoutData<IntVectSet> m_coarseLocations[CH_SPACEDIM*2];
  LayoutData< Vector<Box> > m_coarseLocations[CH_SPACEDIM*2];

  ProblemDomain m_domain;

  int m_isDefined;

  IntVect m_nRefine;
  
  bool m_scaleFineFluxes;

  int m_unrefinedDirection;

  Copier m_reverseCopier;

  bool m_noRealCoarseFineInterface;
private:

  ///there is no operator= for this class
  void operator= (const LevelFluxRegister&);
  ///there is no copy constructor for this class
  LevelFluxRegister(const LevelFluxRegister&);
};
#include "NamespaceFooter.H"
#endif
