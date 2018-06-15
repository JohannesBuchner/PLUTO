#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEUMANNVISCOUSTENSORDOMAINBC_H_
#define _NEUMANNVISCOUSTENSORDOMAINBC_H_

#include "RefCountedPtr.H"

#include "BaseDomainBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

class NeumannViscousTensorDomainBC: public ViscousBaseDomainBC
{
public:
  NeumannViscousTensorDomainBC();

  virtual ~NeumannViscousTensorDomainBC();

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

private:
  void getDivergenceAndGradient(Real&                 a_divergence,
                                Real*                 a_gradient,
                                const int&            a_idir,
                                const FaceIndex&      a_bndryFace,
                                const VolIndex&       a_vof,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_dx,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time,
                                const bool&           a_useHomogeneous);
};

class NeumannViscousTensorDomainBCFactory: public BaseDomainBCFactory
{
public:
  NeumannViscousTensorDomainBCFactory();

  virtual ~NeumannViscousTensorDomainBCFactory();

  ///
  /**
     Sets value of NORMAL gradient of normal component.  Everything
     else in valueland is set to zero.   If you want something more
     elaborate, use the functional form.
  */
  virtual void setValue(Real a_value);

  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_flux);

  virtual NeumannViscousTensorDomainBC* create(const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx);
private:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_flux;
};

#include "NamespaceFooter.H"
#endif
