#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RealVect.H"
#include "ParmParse.H"

#include "functionsF_F.H"
#include "ViscousBCValue.H"
/**/
ViscousTrigBCValue::ViscousTrigBCValue()
{
  m_isDefined = false;

  m_trig = RealVect::Unit;
}
/**/

ViscousTrigBCValue::~ViscousTrigBCValue()
{
}

/**/
void ViscousTrigBCValue::define(const RealVect& a_trig)
{
  m_isDefined = true;

  m_trig = a_trig;
}
RealVect getTrigRVV()
{
  CH_TIME("PoissonUtilities::getTrigRV");
  RealVect trig;
  ParmParse pp;
  std::vector<Real> trigvec(SpaceDim);
  pp.getarr("trig",trigvec,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      trig[idir] = trigvec[idir];
    }
//   Real pi = 4.*atan(1.0);
//   trig *= pi;
  return trig;
}
/**/
Real
ViscousTrigBCValue::
value(const RealVect& a_point,
      const int&      a_comp) const
{
  Real value;
  ParmParse pp;
  int whichMag;
  bool constant_coef;
  pp.get("use_constant_coef", constant_coef);
  if (constant_coef)
    {
      RealVect     trig = getTrigRVV();
      Real time = 0;
      FORT_GETPHIPOINT(CHF_REAL(value),
                       CHF_CONST_REALVECT(trig),
                       CHF_CONST_REALVECT(a_point),
                       CHF_CONST_REAL(time));
    }
  else
    {
      pp.get("which_vel", whichMag);

      FORT_GETMAGPOINTRESIST(CHF_REAL(value),
                             CHF_CONST_REALVECT(m_trig),
                             CHF_CONST_REALVECT(a_point),
                             CHF_CONST_INT(a_comp),
                             CHF_CONST_INT(whichMag));
    }
  return value;
}
/**/
Real
ViscousTrigBCValue::
derivative(const RealVect&       a_point,
           const int&            a_comp,
           const int&            a_derivDir) const
{
  ParmParse pp;
  int whichMag;
  pp.get("which_vel", whichMag);
  Real value;
  bool constant_coef;
  pp.get("use_constant_coef", constant_coef);
  if (constant_coef)
    {
      RealVect     trig = getTrigRVV();
      Real time = 0;
      RealVect gradient;
      FORT_GETGRADPHIPOINT(CHF_REALVECT(gradient),
                           CHF_CONST_REALVECT(trig),
                           CHF_CONST_REALVECT(a_point),
                           CHF_CONST_REAL(time));
      value = gradient[a_derivDir];
    }
  else
    {
  FORT_GETDVDXPOINTRESIST(CHF_REAL(value),
                          CHF_CONST_REALVECT(m_trig),
                          CHF_CONST_REALVECT(a_point),
                          CHF_INT(a_comp),
                          CHF_INT(a_derivDir),
                          CHF_INT(whichMag));
    }
  return value;
}
/**/
