#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIRICHLETPOISSONDOMAINBC_H_
#define _DIRICHLETPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

///
/**
 */
class DirichletPoissonDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  DirichletPoissonDomainBC();

  ///
  /**
   */
  virtual ~DirichletPoissonDomainBC();

  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);

  ///
  /** use this for order of domain boundary
   */
  virtual void setEBOrder(int a_ebOrder);

  ///
  /**
     This is used by the projections to get
     velocity at domain faces where eb intersects.  Sets velocity to a dirichlet value.
   */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);

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
                           const bool&           a_useHomogeneous);

  ///
  /**
   */
  void getHigherOrderFaceFlux(BaseFab<Real>&        a_faceFlux,
                              const BaseFab<Real>&  a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
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
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
   */
  void getInhomFaceFlux(Real&                 a_faceFlux,
                        const VolIndex&       a_vof,
                        const int&            a_comp,
                        const EBCellFAB&      a_phi,
                        const RealVect&       a_probLo,
                        const RealVect&       a_dx,
                        const int&            a_idir,
                        const Side::LoHiSide& a_side,
                        const DataIndex&      a_dit,
                        const Real&           a_time);
  ///
  /**
     A query of whether a_jvof is Dirichlet Domain boundary to a_ivof. a_jvof is the ghost vof; a_ivof is the valid computational vof.
     If the code runs into here it's got to be Dirichlet BC
   */
  bool isDirichletDom(const VolIndex&   a_ivof,
                      const VolIndex&   a_jvof,
                      const EBCellFAB&  a_phi) const
  {
    return true;
  }

  ///
  /**
     This is called by InflowOutflowHelmholtzDomainBC::getFaceFlux,
     which is called by EBAMRPoissonOp::applyOp when EB x domain.
     Used for higher-order velocity boundary conditions.
     This routine iterates over the possible multiple faces at a domain face
     and calls getHigherOrderFaceFlux(...) below for the individual face extrapolation
   */
  void getHigherOrderFaceFlux(Real&                 a_faceFlux,
                              const VolIndex&       a_vof,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useHomogeneous);

  ///
  /**
     This is called by DirichletPoissonDomain::getHigherOrderFaceflux
     which is called by InflowOutflowHelmholtzDomainBC::getFaceFlux.
     This is used by the projections to get extrapolated
     grad(phi) at domain faces for cell-centered gradient.
     Called by getFaceFlux.
   */
  void getHigherOrderFaceFlux(Real&                 a_faceFlux,
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
                              const bool&           a_useHomogeneous);

  ///
  /**
     This is called by InflowOutflowHelmholtzDomainBC::getFaceFlux,
     which is called by EBAMRPoissonOp::applyOp when EB x domain.
     Used for higher-order velocity boundary conditions.
     This routine iterates over the possible multiple faces at a domain face
     and calls getHigherOrderFaceFlux(...) below for the individual face extrapolation
   */
  void getHigherOrderInhomFaceFlux(Real&                 a_faceFlux,
                                   const VolIndex&       a_vof,
                                   const int&            a_comp,
                                   const EBCellFAB&      a_phi,
                                   const RealVect&       a_probLo,
                                   const RealVect&       a_dx,
                                   const int&            a_idir,
                                   const Side::LoHiSide& a_side,
                                   const DataIndex&      a_dit,
                                   const Real&           a_time);

  ///
  /**
     This is called by DirichletPoissonDomain::getHigherOrderFaceflux
     which is called by InflowOutflowHelmholtzDomainBC::getFaceFlux.
     This is used by the projections to get extrapolated
     grad(phi) at domain faces for cell-centered gradient.
     Called by getFaceFlux.
   */
  void getHigherOrderInhomFaceFlux(Real&                 a_faceFlux,
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
                                   const RealVect&       a_centroid);

  ///
  /**
     This is used by the projections to get grad(phi) at domain faces (usually outflow).
     Called by InflowOutflowPoissonDomainBC::getFaceFlux.
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
                              const bool&           a_useHomogeneous);

  ///
  /**
     This is used by the projections to get grad(phi) at domain faces (usually outflow).
     Called by InflowOutflowPoissonDomainBC::getFaceFlux.
   */
  void getInhomFaceGradPhi(Real&                 a_faceFlux,
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
                           const RealVect&       a_centroid);

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const FaceIndex&       a_face,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

  void getFirstOrderFluxStencil(      VoFStencil&      a_stencil,
                                const FaceIndex&       a_face,
                                const int&             a_comp,
                                const RealVect&        a_dx,
                                const int&             a_idir,
                                const Side::LoHiSide&  a_side,
                                const EBISBox&         a_ebisBox);

  void getSecondOrderFluxStencil(      VoFStencil&      a_stencil,
                                 const FaceIndex&       a_face,
                                 const int&             a_comp,
                                 const RealVect&        a_dx,
                                 const int&             a_idir,
                                 const Side::LoHiSide&  a_side,
                                 const EBISBox&         a_ebisBox);

  bool m_isDefined;

private:
  bool m_onlyHomogeneous;
  bool m_isFunctional;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;

  int m_ebOrder;
};

class DirichletPoissonDomainBCFactory: public BaseDomainBCFactory
{
public:
  ///
  /**
   */
  DirichletPoissonDomainBCFactory();

  ///
  /**
   */
  virtual ~DirichletPoissonDomainBCFactory();

  ///
  /**
   */
  virtual DirichletPoissonDomainBC* create(const ProblemDomain& a_domain,
                                           const EBISLayout&    a_layout,
                                           const RealVect&      a_dx);

  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);

  ///
  /**
   */
  virtual void setEBOrder(int a_ebOrder);

private:
  bool m_onlyHomogeneous;
  bool m_isFunctional;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;

  int m_ebOrder;
};
#include "NamespaceFooter.H"
#endif
