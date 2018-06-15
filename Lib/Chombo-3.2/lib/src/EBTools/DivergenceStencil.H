#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIVERGENCESTENCIL_H_
#define _DIVERGENCESTENCIL_H_

#include "AggStencil.H"
#include "Stencils.H"
#include "Vector.H"
#include "REAL.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBFluxFAB.H"
#include "BaseIVFAB.H"
#include "CH_Timer.H"
#include "EBISBox.H"
#include "RefCountedPtr.H"
#include "NamespaceHeader.H"

/// stencil to evaluate the conservative divergence of a flux
/**
 */
class DivergenceStencil
{
public:
  ///
  /**
     The data holders sent in can be dummies but they must be the
     right size.
   */
  DivergenceStencil(const EBCellFAB       & a_diveFluxFaB,
                    const EBFluxFAB       & a_faceFluxFAB,
                    const BaseIVFAB<Real> & a_bndyFluxFAB,
                    const Box             & a_grid,
                    const EBISBox         & a_ebisBox,
                    const RealVect        & a_dx,
                    bool a_useEBFlux = true);

  ///
  ~DivergenceStencil()
  {
  }

  ///
  /**
     Take the divergence of the incoming fluxes
   */
  void divergence(EBCellFAB             & a_divF,
                  const EBFluxFAB       & a_faceFlux,
                  const BaseIVFAB<Real> & a_bndryFlux,
                  const int             & a_destVar,
                  bool incrementOnly = false);
protected:
  bool     m_useEBFlux;
  Box      m_grid;
  RealVect m_dx;
  RefCountedPtr< AggStencil<EBFaceFAB,       EBCellFAB> >  m_openStencil[CH_SPACEDIM];
  RefCountedPtr< AggStencil<BaseIVFAB<Real>, EBCellFAB> >  m_bdryStencil;

private:
  /// disallowed operators.   Without code because Jeff says that is better.

  ///weak construction is bad.
  DivergenceStencil();

  ///deep copy for this object would kill performance
  DivergenceStencil& operator=(const DivergenceStencil& stenin);

  ///no reason to allow this one either.
  DivergenceStencil(const DivergenceStencil& stenin);

};

#include "NamespaceFooter.H"
#endif
