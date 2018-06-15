#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SLIPWALLVISCOUSTENSORDOMAINBC_H_
#define _SLIPWALLVISCOUSTENSORDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCFuncEval.H"

#include "NamespaceHeader.H"

///
/**
 */
class SlipWallViscousTensorDomainBC: public ViscousBaseDomainBC
{
public:
  ///
  /**
   */
  SlipWallViscousTensorDomainBC();

  ///
  /**
   */
  virtual ~SlipWallViscousTensorDomainBC();


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

class SlipWallViscousTensorDomainBCFactory: public BaseDomainBCFactory
{
public:
  ///
  /**
   */
  SlipWallViscousTensorDomainBCFactory();

  ///
  /**
   */
  virtual ~SlipWallViscousTensorDomainBCFactory();

  ///
  /**
   */
  virtual SlipWallViscousTensorDomainBC* create(const ProblemDomain& a_domain,
                                                const EBISLayout&    a_layout,
                                                const RealVect&      a_dx);

private:

};

#include "NamespaceFooter.H"
#endif
