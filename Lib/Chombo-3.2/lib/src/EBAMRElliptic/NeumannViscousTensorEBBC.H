#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEUMANNVISCOUSTENSOREBBC_H_
#define _NEUMANNVISCOUSTENSOREBBC_H_

#include "RefCountedPtr.H"

#include "BaseEBBC.H"
#include "BaseBCFuncEval.H"
#include "NamespaceHeader.H"

///
/**
 */
class NeumannViscousTensorEBBC: public ViscousBaseEBBC
{
public:
  ///
  /**
   */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor)
  {
  }


  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time);
  ///
  /**
   */
  NeumannViscousTensorEBBC(const ProblemDomain& a_domain,
                           const EBISLayout&    a_layout,
                           const RealVect&      a_dx,
                           const IntVect*       a_ghostCellsPhi,
                           const IntVect*       a_ghostCellsRhs);


  //no homogeneous stencil.
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    return NULL;
  }

  ///
  /**
   */
  virtual ~NeumannViscousTensorEBBC();


  ///
  /**
     deprecated
  */
  virtual void getEBFlux(Real&                         a_flux,
                         const VolIndex&               a_vof,
                         const LevelData<EBCellFAB>&   a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&              a_dit,
                         const RealVect&               a_probLo,
                         const RealVect&               a_dx,
                         const bool&                   a_useHomogeneous,
                         const Real&                   a_time,
                         const pair<int,Real>*         a_cacheHint=0 )
  {
    MayDay::Error();
  }

  NeumannViscousTensorEBBC()
  {
  }

private:
};

///
/**
 */
class NeumannViscousTensorEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  /**
   */
  NeumannViscousTensorEBBCFactory()
  {;}

  ///
  /**
   */
  virtual ~NeumannViscousTensorEBBCFactory()
  {;}

  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_flux);

  ///
  /**
   */
  virtual NeumannViscousTensorEBBC* create(const ProblemDomain& a_domain,
                                           const EBISLayout&    a_layout,
                                           const RealVect&      a_dx,
                                           const IntVect*       a_ghostCellsPhi=0,
                                           const IntVect*       a_ghostCellsRhs=0);

private:
  bool m_isFunction;
  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_func;

};

#include "NamespaceFooter.H"
#endif
