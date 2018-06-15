#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIRICHLETCONDUCTIVITYDOMAINBC_H_
#define _DIRICHLETCONDUCTIVITYDOMAINBC_H_

#include "RefCountedPtr.H"

#include "DirichletPoissonDomainBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

class DirichletConductivityDomainBC: public ConductivityBaseDomainBC
{
public:
  DirichletConductivityDomainBC();

  virtual ~DirichletConductivityDomainBC();

  virtual void setValue(Real a_value);

  virtual int whichBC(int                  a_idir,
                      Side::LoHiSide       a_side);

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


  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi)
  {
    return true;
  }

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
  bool m_isFunctional;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;
};

class DirichletConductivityDomainBCFactory: public BaseDomainBCFactory
{
public:
  DirichletConductivityDomainBCFactory();

  virtual ~DirichletConductivityDomainBCFactory();

  virtual void setValue(Real a_value);

  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  virtual DirichletConductivityDomainBC* create(const ProblemDomain& a_domain,
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
