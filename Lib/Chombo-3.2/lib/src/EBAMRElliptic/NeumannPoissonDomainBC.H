#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEUMANNPOISSONDOMAINBC_H_
#define _NEUMANNPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"

#include "BaseDomainBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

class NeumannPoissonDomainBC: public BaseDomainBC
{
public:
  NeumannPoissonDomainBC();

  //no homogeneous stencil.
  // virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  // {
  //   return NULL;
  // }

  virtual ~NeumannPoissonDomainBC();

  virtual void setValue(Real a_value);

  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

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
     A query of whether a_jvof is Dirichlet Domain boundary to a_ivof. a_jvof is the ghost vof; a_ivof is the valid computational vof.
     If the code runs into here it's got to be Neumann BC
   */
  bool isDirichletDom(const VolIndex&   a_ivof,
                      const VolIndex&   a_jvof,
                      const EBCellFAB&  a_phi) const
  {
    return false;
  }
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

  virtual void getInhomFaceFlux(Real&                 a_faceFlux,
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
     This is used by the operator to get
     grad(phi) at domain faces.
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
     This is used by the projections to get extrapolated
     grad(phi) at domain faces for cell-centered gradient.
     Called by getFaceFlux.
   */
  virtual void getFaceFluxGradPhi(Real&                 a_faceFlux,
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
   */
  virtual void getInhomFaceGradPhi(Real&                 a_faceFlux,
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
   */
  virtual void getInhomFaceFluxGradPhi(Real&                 a_faceFlux,
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
     This is used by the projections to get
     velocity at domain faces.  This one just returns
     the value in a_vel(a_face).
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

private:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_flux;
};

class NeumannPoissonDomainBCFactory: public BaseDomainBCFactory
{
public:
  NeumannPoissonDomainBCFactory();

  virtual ~NeumannPoissonDomainBCFactory();

  virtual void setValue(Real a_value);

  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  virtual NeumannPoissonDomainBC* create(const ProblemDomain& a_domain,
                                         const EBISLayout&    a_layout,
                                         const RealVect&      a_dx);
private:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_flux;
};

#include "NamespaceFooter.H"
#endif
