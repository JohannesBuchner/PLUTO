#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "BoxIterator.H"

#include "NeumannConductivityDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletConductivityEBBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"
#include "SchlichtingTempDomainBC.H"

#include "NamespaceHeader.H"

////
void 
SchlichtingTempDomainBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  //velocity is neumann only at inflow or at slipwalls

  NeumannConductivityDomainBC neumannBC;
  neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
  neumannBC.setValue(0.0);
  neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
}


void 
SchlichtingTempDomainBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  NeumannConductivityDomainBC neumannBC;
  neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
  neumannBC.setValue(0.0);
  neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                        a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
}

#include "NamespaceFooter.H"
