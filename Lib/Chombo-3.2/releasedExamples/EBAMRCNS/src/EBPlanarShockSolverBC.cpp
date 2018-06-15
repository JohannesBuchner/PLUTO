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

#include "EBPlanarShockSolverBC.H"
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

class DVTInflowFunc: public  BaseBCFuncEval
{
public:
  DVTInflowFunc(Real a_inflowvel, int a_inflowdir) 
  {
    m_inflowvel = a_inflowvel;
    m_inflowdir = a_inflowdir;
  }


  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp) const
  {
    Real retval = 0;
    if(a_comp == m_inflowdir) retval = m_inflowvel;

    return retval;
  }

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir
                          ) const 
  {
    return 0.;
  }
protected:

  Real  m_inflowvel;
  int   m_inflowdir;

};



void 
EBPlanarShockSolverBC::
whereAMI(bool& a_atInflow, 
         bool& a_atOutflow,
         const int&            a_idir, 
         const Side::LoHiSide& a_side)
{

  if(a_idir == m_shockNorm)
    {
      a_atOutflow =  ((  m_shockBackward && a_side == Side::Lo) ||
                      ( !m_shockBackward && a_side == Side::Hi));
      
      a_atInflow =  (( !m_shockBackward && a_side == Side::Lo) ||
                     (  m_shockBackward && a_side == Side::Hi));
    }
  else
    {
      a_atInflow = false;
      a_atOutflow= false;
    }

}

////
void 
EBPlanarShockSolverBC::
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

  bool atInflow, atOutflow;
  whereAMI(atInflow, atOutflow, a_idir, a_side);
  if(atInflow)
    {
      DirichletViscousTensorDomainBC diriBC;
      Real value = 0;
      FORT_GETPOSTSHOCKVEL(CHF_REAL(value));
      RefCountedPtr<BaseBCFuncEval> funk(new DVTInflowFunc (value, m_shockNorm));
      diriBC.setFunction(funk);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);

     // NeumannViscousTensorDomainBC neumBC;
     // neumBC.setValue(0.);
     //
     // neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
     // neumBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(atOutflow)
    {
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(m_slipWall)
    {
//      SlipWallViscousTensorDomainBC slipbc;
//      slipbc.setCoef(m_eblg, m_beta, m_eta, m_lambda);
//      slipbc.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);

      NeumannViscousTensorDomainBC  neumBC;
      Real value = 0;
      neumBC.setValue(value);

      neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      //no slip wall
      DirichletViscousTensorDomainBC diriBC;
      Real value = 0;
      diriBC.setValue(value);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void 
EBPlanarShockSolverBC::
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

  bool atInflow, atOutflow;
  whereAMI(atInflow, atOutflow, a_idir, a_side);

  if(atInflow)
    {
      DirichletViscousTensorDomainBC diriBC;
      Real inflowvel;
      FORT_GETPOSTSHOCKVEL(CHF_REAL(inflowvel));

      RefCountedPtr<BaseBCFuncEval> funk(new DVTInflowFunc (inflowvel, m_shockNorm));
      diriBC.setFunction(funk);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
//      NeumannViscousTensorDomainBC neumBC;
//      neumBC.setValue(0.);
//      neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
//      neumBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
//                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(atOutflow)
    {
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(m_slipWall)
    {
      SlipWallViscousTensorDomainBC slipbc;
      slipbc.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      slipbc.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
//      NeumannViscousTensorDomainBC neumBC;
//      Real value = 0;
//      neumBC.setValue(value);
//      neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
//      neumBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
//                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      //no slip wall
      DirichletViscousTensorDomainBC diriBC;
      Real value = 0;
      diriBC.setValue(value);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

#include "NamespaceFooter.H"
