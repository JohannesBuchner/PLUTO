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

#include "EBPlanarShockTemperatureBC.H"
#include "NeumannConductivityDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletConductivityEBBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"
#include "EBPlanarShockF_F.H"

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
EBPlanarShockTemperatureBC::
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
EBPlanarShockTemperatureBC::
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

  bool atInflow, atOutflow;
  whereAMI(atInflow, atOutflow, a_idir, a_side);
  //  if(atInflow)
  if(0)
    {
      DirichletConductivityDomainBC diriBC;
      Real value = 0;
      FORT_GETPOSTSHOCKTEMP(CHF_REAL(value));
      diriBC.setValue(value);

      diriBC.setCoef(m_eblg, m_beta, m_bcoef);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      NeumannConductivityDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void 
EBPlanarShockTemperatureBC::
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
  //velocity is neumann only at inflow or at slipwalls
  bool atInflow, atOutflow;
  whereAMI(atInflow, atOutflow, a_idir, a_side);

  //  if(atInflow)
  if(0)
    {
      DirichletConductivityDomainBC diriBC;
      Real value;
      FORT_GETPOSTSHOCKTEMP(CHF_REAL(value));
      diriBC.setValue(value);

      diriBC.setCoef(m_eblg, m_beta, m_bcoef);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else 
    {
      NeumannConductivityDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

#include "NamespaceFooter.H"
