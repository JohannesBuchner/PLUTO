#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEDOMAINBC_H_
#define _BASEDOMAINBC_H_

#include "LoHiSide.H"
#include "RealVect.H"
#include "FArrayBox.H"

#include "VolIndex.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBFluxFAB.H"
#include "EBISLayout.H"
#include "EBLevelGrid.H"
#include "TensorCFInterp.H"
#include "BaseBCFuncEval.H"
#include "NamespaceHeader.H"

///
/**
 */
class BaseDomainBC
{
public:
  ///
  /**
   */
  BaseDomainBC()
  {
  }

  ///
  /**
   */
  virtual ~BaseDomainBC()
  {
  }

  ///
  /**
   */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous) = 0;

  // virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)= 0;

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox)
  {
    a_stencil.clear();
  }

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const FaceIndex&       a_face,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox)
  {
    a_stencil.clear();
  }

  ///
  /**
     This is used by the elliptic solvers to get
     the flux of stuff through domain faces.
   */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous) = 0;


  ///
  /**
     A query of whether a_jvof is Dirichlet Domain boundary to a_ivof. a_jvof is the ghost vof; a_ivof is the valid computational vof.
     For now this is mainly used in form_matrix, which needs to know the type of BC at the boundary.
     Default implementation is MayDay::Error();
   */
  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi) const
   {
     MayDay::Error("BaseDomainBC:: isDirichletDom:: default implementation is invalid");
     return true;
   }


  ///
  /**
     This is used by the elliptic solvers to get
     the flux of stuff through domain faces.
   */
  virtual void getInhomFaceFlux(Real&                 a_faceFlux,
                                const VolIndex&       a_vof,
                                const int&            a_comp,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_dx,
                                const int&            a_idir,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time)
  {
    MayDay::Error();
  }

  ///
  /**
     This is used by the projections to get
     a) grad(phi) at domain faces.
   */
  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous) = 0;


  virtual void enforceFaceVel(LevelData<EBFluxFAB>&    a_velocity,
                              const DisjointBoxLayout& a_grids,
                              const EBISLayout&        a_ebisl,
                              const ProblemDomain&     a_domain,
                              const RealVect&          a_dx,
                              const Real&              a_time,
                              const RealVect&          a_origin,
                              const bool&              a_doDivFreeOutflow);

  virtual void enforceFaceVel(LevelData<EBFluxFAB>&    a_velocity,
                              const DisjointBoxLayout& a_grids,
                              const EBISLayout&        a_ebisl,
                              const ProblemDomain&     a_domain,
                              const RealVect&          a_dx,
                              const Real&              a_time,
                              const RealVect&          a_origin,
                              const bool&              a_doDivFreeOutflow,
                              const int&               a_comp);

  ///
  /**
     This is used by the projections to get
     velocity at domain faces.
   */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_vof,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow) = 0;

};

///
/**
 */
class ViscousBaseDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  ViscousBaseDomainBC()
  {
    m_coefSet  = false;
    m_value = 12345.6789;
    m_func = RefCountedPtr<BaseBCFuncEval>();
    m_isFunction = false;
  }

  virtual ~ViscousBaseDomainBC()
  {
  }
  ///
  /**
   */
  void setCoef(EBLevelGrid                         &  a_eblg,
               Real                                &  a_beta,
               RefCountedPtr<LevelData<EBFluxFAB> >&  a_eta,
               RefCountedPtr<LevelData<EBFluxFAB> >&  a_lambda)
  {
    m_coefSet = true;
    m_eblg    = a_eblg;
    m_beta    = a_beta;
    m_eta     = a_eta;
    m_lambda  = a_lambda;
  }

  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_vof,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow)
  {
    //projection thingy does not apply here
    MayDay::Error();
  }
  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous)
  {
    //projection thingy does not apply here
    MayDay::Error();
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

  void getFluxFromGrad(BaseFab<Real>&   a_flux,
                       const FArrayBox& a_grad,
                       const DataIndex& a_dit,
                       const int&       a_idir);

protected:
  bool m_isFunction;


  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_func;

  Real                                  m_beta;
  bool                                  m_coefSet;
  EBLevelGrid                           m_eblg;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_eta;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_lambda;
};

///
/**
 */
class ConductivityBaseDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  ConductivityBaseDomainBC()
  {
    m_coefSet  = false;
  }

  virtual ~ConductivityBaseDomainBC()
  {
  }

  ///
  /**
   */
  void setCoef(EBLevelGrid                         &  a_eblg,
               Real                                &  a_beta,
               RefCountedPtr<LevelData<EBFluxFAB> >&  a_bcoef)
  {
    m_coefSet = true;
    m_eblg    = a_eblg;
    m_beta    = a_beta;
    m_bcoef   = a_bcoef;
  }


protected:
  bool                                  m_coefSet;
  EBLevelGrid                           m_eblg;
  Real                                  m_beta;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_bcoef;
};
///
/**
 */
class BaseDomainBCFactory
{
public:
  ///
  /**
   */
  BaseDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual ~BaseDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual BaseDomainBC* create(const ProblemDomain& a_domain,
                               const EBISLayout&    a_layout,
                               const RealVect&      a_dx) = 0;
};
#include "NamespaceFooter.H"
#endif
