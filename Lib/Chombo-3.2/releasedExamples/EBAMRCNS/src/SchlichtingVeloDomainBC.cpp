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

#include "SchlichtingVeloDomainBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletViscousTensorEBBC.H"
#include "NeumannViscousTensorEBBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"
#include "EBPlanarShockF_F.H"
#include "SlipWallViscousTensorDomainBC.H"

#include "NamespaceHeader.H"

///
class SchlichtingInflowFunc: public BaseBCFuncEval
{
public:
  SchlichtingInflowFunc()
  {
  }

  SchlichtingInflowFunc(const SchlichtingParams& a_params)
  {
    m_paramsSet = true;
    m_params = a_params;
  }

  void setParams(const SchlichtingParams& a_params)
  {
    m_paramsSet = true;
    m_params = a_params;
  }

  virtual ~SchlichtingInflowFunc()
  {
  }

  virtual Real value(const RealVect& a_point,
                     const int&      a_comp) const
  {
    CH_assert(m_paramsSet);
    return m_params.m_initVel[a_comp];
  }

  virtual Real derivative(const RealVect& a_point,
                          const int&      a_comp,
                          const int&      a_derivDir) const
  {
    return 0.;
  }

  SchlichtingParams m_params;
  bool m_paramsSet;
};

////
void 
SchlichtingVeloDomainBC::
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
    bool atInflow = (a_side == Side::Lo);

  if(atInflow)
    {
      RefCountedPtr<BaseBCFuncEval> func(new SchlichtingInflowFunc(m_params));
      DirichletViscousTensorDomainBC diribc;
      diribc.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diribc.setFunction(func);
      diribc.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
      
    }
  else
    {
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void 
SchlichtingVeloDomainBC::
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
  bool atInflow = (a_side == Side::Lo);

  if(atInflow)
    {
      RefCountedPtr<BaseBCFuncEval> func(new SchlichtingInflowFunc(m_params));
      DirichletViscousTensorDomainBC diribc;
      diribc.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diribc.setFunction(func);
      diribc.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);

    }
  else
    {
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

#include "NamespaceFooter.H"
