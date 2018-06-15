#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEEBBC_H_
#define _BASEEBBC_H_

#include "RealVect.H"
#include "ProblemDomain.H"

#include "VolIndex.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBFluxFAB.H"
#include "EBStencil.H"
#include "Stencils.H"
#include "EBLevelGrid.H"
#include "BaseBCFuncEval.H"
#include "NamespaceHeader.H"


///
/**
 */
class BaseEBBC
{
public:
  BaseEBBC()
  {
    m_dataBased = false;
  }

  virtual ~BaseEBBC()
  {
  }

  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor) = 0;

  ///deprecated interface.
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
    //if you want to call it, you have to implement it.
    MayDay::Error("how did this get called?");
  }

  ///
  /**
     Return a pointer to the homogenous flux stencil for the boundary condition.
     contribution.   In the case where there in no contribution in the homogeneous case
     (ie. homogeneous Neumann) return NULL.
  */
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)= 0;

  ///
  /**
     add change in lphi due to eb flux
  */
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time) = 0;

  virtual void setData( RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_data)
  {
    m_data = a_data;
    m_dataBased = true;
  }

  bool dataBased() const
  {
    return m_dataBased;
  }
protected:

  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_data;
  bool m_dataBased;

};

///
/**
 */
class ViscousBaseEBBC: public BaseEBBC
{
public:
  ///
  /**
   */
  ViscousBaseEBBC()
  {
    m_coefSet  = false;
    m_value = 12345.6789;
    m_func = RefCountedPtr<BaseBCFuncEval>();
    m_isFunction = false;
  }
  virtual ~ViscousBaseEBBC()
  {
  }

  ///
  /**
   */
  void setCoef(EBLevelGrid                                &  a_eblg,
               Real                                       &  a_beta,
               RefCountedPtr<LevelData<BaseIVFAB<Real> > >&  a_eta,
               RefCountedPtr<LevelData<BaseIVFAB<Real> > >&  a_lambda,
               RefCountedPtr<LevelData< EBFluxFAB > >&  a_etaOpen,
               RefCountedPtr<LevelData< EBFluxFAB > >&  a_lambdaOpen
               )
  {
    m_coefSet = true;
    m_beta    = a_beta;
    m_eblg    = a_eblg;
    m_eta     = a_eta;
    m_lambda  = a_lambda;
    m_etaOpen    = a_etaOpen;
    m_lambdaOpen = a_lambdaOpen;
  }

  ///
  /**
   */
  virtual void setValue(Real a_value)
  {
    m_value = a_value;
    m_func = RefCountedPtr<BaseBCFuncEval>();

    m_isFunction = false;
  }


  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_func)
  {
    m_value = 12345.6789;
    m_func = a_func;

    m_isFunction = true;
  }

  void
  getBoundaryGrad(Real                      a_grad[CH_SPACEDIM][CH_SPACEDIM],
                  const VolIndex&           a_vof,
                  const RealVect&           a_dx,
                  const RealVect&           a_probLo,
                  const EBISBox&            a_ebisBox)
  {
    for (int comp = 0; comp < SpaceDim; comp++)
      {
        for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
          {
            Real value;
            if (m_isFunction)
              {
                // Compute the bndryCentroid location in physical coordinates
                RealVect startPt = a_ebisBox.bndryCentroid(a_vof);
                startPt *= a_dx[0];
                startPt += a_probLo;
                RealVect point = EBArith::getVofLocation(a_vof, a_dx, startPt);
                value = m_func->derivative(point, comp, derivDir);
              }
            else
              {
                value = m_value;
              }
            a_grad[comp][derivDir] = value;
          }
      }
  }

  void
  getFluxFromGrad(Real             a_flux[CH_SPACEDIM][CH_SPACEDIM],
                  const Real       a_grad[CH_SPACEDIM][CH_SPACEDIM],
                  const VolIndex&  a_vof,
                  const DataIndex& a_dit)
  {
    //compute divergence at boundary face
    Real divergence = 0;
    for (int divDir = 0; divDir < SpaceDim; divDir++)
      {
        divergence += a_grad[divDir][divDir];
      }

    Real lambda = (*m_lambda)[a_dit](a_vof, 0);
    Real    eta =    (*m_eta)[a_dit](a_vof, 0);

    for (int comp = 0; comp < SpaceDim; comp++)
      {
        for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
          {
            a_flux[comp][derivDir] = eta*(a_grad[comp][derivDir] + a_grad[derivDir][comp]);
            if (comp == derivDir)
              {
                a_flux[comp][derivDir] += lambda*divergence;
              }
          }
      }
  }

  void
  getChangeInSolution(Real                      a_deltaLph[CH_SPACEDIM],
                      const Real                a_flux[CH_SPACEDIM][CH_SPACEDIM],
                      const RealVect&           a_dx,
                      const VolIndex&           a_vof,
                      const DataIndex&          a_dit,
                      const EBISBox&            a_ebisBox)
  {
    Real   beta =      m_beta;
    RealVect normal =  a_ebisBox.normal(a_vof);
    Real areaFrac = a_ebisBox.bndryArea(a_vof);
    //change in solution from one flux is
    // dphi = beta * div F = beta * bndryArea*(f dot normal)/dx
    // = beta *sum_dir(F[comp][derivDir]*n_derivDir/dx
    for (int comp = 0; comp < SpaceDim; comp++)
      {
        a_deltaLph[comp] = 0;
        for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
          {
            a_deltaLph[comp] -= a_flux[comp][derivDir]*beta*areaFrac*normal[derivDir]/a_dx[0];
          }
      }
  }
protected:
  bool m_isFunction;


  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_func;


  EBLevelGrid                                  m_eblg;
  bool                                         m_coefSet;
  Real                                         m_beta;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_eta;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_lambda;

  RefCountedPtr<LevelData<EBFluxFAB > >  m_etaOpen;
  RefCountedPtr<LevelData<EBFluxFAB > >  m_lambdaOpen;

};

///
/**
 */
class ConductivityBaseEBBC: public BaseEBBC
{
public:
  ///
  /**
   */
  ConductivityBaseEBBC()
  {
    m_coefSet  = false;
    m_dataBased = false;
  }

  virtual ~ConductivityBaseEBBC()
  {
  }


  ///
  /**
   */
  void setCoef(EBLevelGrid                                &  a_eblg,
               Real                                       &  a_beta,
               RefCountedPtr<LevelData<BaseIVFAB<Real> > >&  a_bcoe)
  {
    m_coefSet = true;
    m_beta    = a_beta;
    m_eblg    = a_eblg;
    m_bcoe    = a_bcoe;
  }


protected:

  EBLevelGrid                                  m_eblg;
  bool                                         m_coefSet;
  Real                                         m_beta;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_bcoe;

};
///
/**
 */
class BaseEBBCFactory
{
public:
  BaseEBBCFactory()
  {
  }

  virtual ~BaseEBBCFactory()
  {
  }

  ///
  /**
   */
  virtual BaseEBBC* create(const ProblemDomain& a_domain,
                           const EBISLayout&    a_layout,
                           const RealVect&      a_dx,
                           const IntVect*       a_ghostCellsPhi=0,
                           const IntVect*       a_ghostCellsRhs=0 ) = 0;
};

#include "NamespaceFooter.H"
#endif
