#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SCHLICHTINGVELOCITYEBBC_H_
#define _SCHLICHTINGVELOCITYEBBC_H_

#include "RefCountedPtr.H"

#include "BaseEBBC.H"
#include "BaseBCFuncEval.H"
#include "LayoutData.H"
#include "IntVectSet.H"
#include "EBStencil.H"
#include "EBAMRCNSParams.H"
#include "DirichletViscousTensorEBBC.H"

#include "NamespaceHeader.H"

class VoFCache;
class VoFStencil;

///
/**
 */
class SchlichtingVelocityEBBC: public DirichletViscousTensorEBBC
{
public:

  ///
  SchlichtingVelocityEBBC();

  SchlichtingVelocityEBBC(const SchlichtingParams& a_params,
                          const ProblemDomain& a_domain,
                          const EBISLayout&    a_layout,
                          const RealVect&      a_dx);

  /// 
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor);
  ///
  /**
   */
  virtual ~SchlichtingVelocityEBBC();

  //always homogeneous so can simplify
  virtual void
  applyEBFlux(EBCellFAB&                    a_lphi,
              const EBCellFAB&              a_phi,
              VoFIterator&                  a_vofit,
              const LayoutData<IntVectSet>& a_cfivs,
              const DataIndex&              a_dit,
              const RealVect&               a_probLo,
              const RealVect&               a_dx,
              const Real&                   a_factor,
              const bool&                   a_useHomogeneous,
              const Real&                   a_time);

  //always homogeneous so can simplify
  virtual RealVect
  getInhomogeneousContribution(const VolIndex&  a_vof,
                               const EBCellFAB& a_phi,
                               const EBISBox&   a_ebisBox,
                               const DataIndex& a_dit,
                               const Real&      a_dx);
protected:

  SchlichtingParams m_params;
  //returns stencil of flux*normal extrapolated from interior
  bool getSchlichtingStencil(VoFStencil&       a_stencil, 
                             const VolIndex&   a_vof, 
                             const DataIndex&  a_dit,
                             const RealVect&   a_normal, 
                             int a_ivar);
  void
  getAverageFluxStencil(VoFStencil& a_aveFluxSten, 
                        const VolIndex&   a_vof, 
                        const DataIndex&  a_dit,
                        int a_idir, int   a_ivar);
};

///
/**
 */
class SchlichtingVelocityEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  /**
   */
  SchlichtingVelocityEBBCFactory(const SchlichtingParams& a_params);

  ///
  /**
   */
  virtual ~SchlichtingVelocityEBBCFactory();


  ///
  /**
   */
  virtual SchlichtingVelocityEBBC* create(  const ProblemDomain& a_domain,
                                            const EBISLayout&    a_layout,
                                            const RealVect&      a_dx,
                                            const IntVect*       a_ghostCellsPhi,
                                            const IntVect*       a_ghostCellsRhs);

private:

  SchlichtingParams m_params;
};

#include "NamespaceFooter.H"
#endif
