#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBAMRCNSParams.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "NamespaceHeader.H"

Real 
SchlichtingParams::
distanceAlongAxis(const RealVect& a_point) const
{
  RealVect diffVect = a_point - m_corner;
  Real dist = PolyGeom::dot(diffVect, m_axis);
  return dist;
}


Real 
SchlichtingParams::
distanceFromPlane(const RealVect& a_point) const
{
  if(SpaceDim==3)
    {
      MayDay::Error("this function is wrong in 3d");
    }
  RealVect normalToPlane[CH_SPACEDIM-1];
  PolyGeom::getTangentVectors(normalToPlane, m_axis);

  //I happen to know which plane I am aiming for
  RealVect planeOrigin = m_corner;

  RealVect diffVect = a_point - planeOrigin;
  Real dist = PolyGeom::dot(diffVect, normalToPlane[0]);
  return dist;
}

bool 
SchlichtingParams::
isThisPointSticky(const RealVect& a_point,const RealVect& a_normal) const
{
  //the wall pointing downward is not sticky.  In 3d, same goes for the slip walls on the sides.
  //but I do not know exactly how to characterize that.
  bool retval; 
  if(a_normal[1] < 0.)
    {
      retval = false;
    }
  else
    {
      Real distance = distanceAlongAxis(a_point);
      retval = ((distance > m_stickyStart) && (distance < m_stickyEnd));
    }
  return retval;
}

void
ParseSchlichtingParams(SchlichtingParams& a_params)
{
  ParmParse pp;
  pp.get("gamma"           ,   a_params.m_gamma);
  pp.get("shock_mach"      ,   a_params.m_ms);
  pp.get("sticky_start"    ,   a_params.m_stickyStart);
  pp.get("sticky_end"      ,   a_params.m_stickyEnd);
  pp.get("initial_density" ,   a_params.m_initDense);
  pp.get("initial_pressure",   a_params.m_initPress);
  pp.get("specific_heat"   ,   a_params.m_specificHeat);
  
  Real R = a_params.m_specificHeat*(a_params.m_gamma - 1.0);
  Real denom = a_params.m_initDense*R;
  a_params.m_initTemp = a_params.m_initPress/(denom);

  vector<Real>  schlichtingOriginVect(SpaceDim);
  vector<Real>    schlichtingAxisVect(SpaceDim);

  pp.getarr("schlichting_axis" ,schlichtingAxisVect,   0, SpaceDim);
  pp.getarr("schlichting_origin",schlichtingOriginVect, 0, SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.m_axis[idir]   = schlichtingAxisVect[idir];
      a_params.m_corner[idir] = schlichtingOriginVect[idir];
    }
  Real sum;
  PolyGeom::unifyVector(a_params.m_axis, sum);
  
  Real initialP = a_params.m_initPress;
  Real initialR = a_params.m_initDense;

  Real sound = sqrt(a_params.m_gamma*initialP/initialR);
  Real velMag = sound*a_params.m_ms;
  a_params.m_velMag = velMag;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.m_initVel[idir] = a_params.m_axis[idir]*velMag;
    }
}

ostream&
operator<< (ostream&              a_os,
            const EBAMRCNSParams& a_p)
{
  a_os << "AMRCNSParams:" << endl;

  a_os << "useAirCoefs         =" << a_p.m_useAirCoefs         << endl;
  a_os << "variableCoeff       =" << a_p.m_variableCoeff       << endl;
  a_os << "doDiffusion         =" << a_p.m_doDiffusion         << endl;
  a_os << "slipBoundaries      =" << a_p.m_slipBoundaries      << endl;
  a_os << "useMassRedist       =" << a_p.m_useMassRedist       << endl;
  a_os << "doSmushing          =" << a_p.m_doSmushing          << endl;
  a_os << "redistRad           =" << a_p.m_redistRad           << endl;
  a_os << "tagBufferSize       =" << a_p.m_tagBufferSize       << endl;
  a_os << "verbosity           =" << a_p.m_verbosity           << endl;
  a_os << "domainLength        =" << a_p.m_domainLength        << endl;
  a_os << "initialDtMultiplier =" << a_p.m_initialDtMultiplier << endl;
  a_os << "cfl                 =" << a_p.m_cfl                 << endl;
  a_os << "refineThresh        =" << a_p.m_refineThresh        << endl;
  if(!a_p.m_useAirCoefs)                                             
    {                                                                
      a_os << "specHeatCv          =" << a_p.m_specHeatCv      << endl;
      a_os << "thermalCond         =" << a_p.m_thermalCond     << endl;
      a_os << "viscosityMu         =" << a_p.m_viscosityMu     << endl;
      a_os << "viscosityLa         =" << a_p.m_viscosityLa     << endl;      
    }

  return a_os;
}

#include "NamespaceFooter.H"
