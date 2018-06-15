#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIRICHLETVISCOUSTENSORDOMAINBC_H_
#define _DIRICHLETVISCOUSTENSORDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCFuncEval.H"
#include "NamespaceHeader.H"

///
/**
 */
class DirichletViscousTensorDomainBC: public ViscousBaseDomainBC
{
public:
  ///
  /**
   */
  DirichletViscousTensorDomainBC();

  ///
  /**
   */
  virtual ~DirichletViscousTensorDomainBC();


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

class DirichletViscousTensorDomainBCFactory: public BaseDomainBCFactory
{
public:
  ///
  /**
   */
  DirichletViscousTensorDomainBCFactory();

  ///
  /**
   */
  virtual ~DirichletViscousTensorDomainBCFactory();

  ///
  /**
   */
  virtual DirichletViscousTensorDomainBC* create(const ProblemDomain& a_domain,
                                                 const EBISLayout&    a_layout,
                                                 const RealVect&      a_dx);

  ///
  /**
     Sets value of of normal component.  Other components set to zero.
  */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_func);

private:
  bool m_onlyHomogeneous;
  bool m_isFunctional;

  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_func;
};
#include "NamespaceFooter.H"
#endif
