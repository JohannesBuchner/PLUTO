#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRCNSParams_H_
#define _EBAMRCNSParams_H_

#include "Vector.H"
#include "RealVect.H"
#include "REAL.H"
#include "BaseDomainBC.H"
#include "BaseEBBC.H"

#include "NamespaceHeader.H"

///
class EBAMRCNSParams
{
public:
  bool                                m_doLazyRelax;
  bool                                m_useAirCoefs;
  bool                                m_variableCoeff;
  bool                                m_checkMaxMin;
  bool                                m_doDiffusion;
  bool                                m_slipBoundaries;
  bool                                m_useMassRedist;
  bool                                m_doSmushing;
  bool                                m_prepostSmooth;
  bool                                m_explicitHyperbolicSource;
  bool                                m_backwardEuler;
  int                                 m_redistRad;
  int                                 m_tagBufferSize;
  int                                 m_verbosity;
  Real                                m_domainLength;
  Real                                m_initialDtMultiplier;
  Real                                m_cfl;
  Real                                m_refineThresh;
  Real                                m_specHeatCv;
  Real                                m_thermalCond;
  Real                                m_viscosityMu;
  Real                                m_viscosityLa;
  Real                                m_regridSmoothCoef;
  RefCountedPtr<BaseDomainBCFactory>  m_doBCTemp;
  RefCountedPtr<BaseDomainBCFactory>  m_doBCVelo;
  RefCountedPtr<BaseDomainBCFactory>  m_doBCSmooth;
  RefCountedPtr<BaseEBBCFactory>      m_ebBCTemp;
  RefCountedPtr<BaseEBBCFactory>      m_ebBCVelo;
  RefCountedPtr<BaseEBBCFactory>      m_ebBCSmooth;
  EBAMRCNSParams()
  {
    m_regridSmoothCoef    = 1.0e-4;
    m_initialDtMultiplier = 0.1;
    m_cfl                 = 0.5;
    m_redistRad           = 1;
    m_domainLength        = 1.0;
    m_refineThresh        = -123456789;
    m_tagBufferSize       = 1;
    m_verbosity           = 4;
    m_useMassRedist       = true;
    m_doSmushing          = true;
    m_prepostSmooth       = true;
    m_viscosityMu         = 1.0e-3;
    m_viscosityLa         = 1.0e-3;
    m_specHeatCv          = 1.0e0;
    m_thermalCond         = 1.0e-3;
    m_useAirCoefs         = false;
    m_doLazyRelax         = true;
    m_explicitHyperbolicSource = false;
    m_backwardEuler       = true;
    m_variableCoeff = m_useAirCoefs;
    m_checkMaxMin   = false;
  }
  ~EBAMRCNSParams()
  {;}
};

///
class SchlichtingParams
{
public:

  Real         m_gamma;
  Real         m_ms;
  Real         m_stickyStart;
  Real         m_stickyEnd;  
  Real         m_initPress;
  Real         m_initDense;
  Real         m_specificHeat;
  Real         m_initTemp;
  Real         m_velMag;
  RealVect     m_axis;
  RealVect     m_corner;
  RealVect     m_initVel;

  Real distanceAlongAxis(const RealVect& a_point) const;
  Real distanceFromPlane(const RealVect& a_point) const;
  bool isThisPointSticky(const RealVect& a_point,const RealVect& a_normal) const;
};

///get stuff from input file
extern void ParseSchlichtingParams(SchlichtingParams& a_params);

extern ostream&
operator<< (ostream&              a_os,
            const EBAMRCNSParams& a_p);

#include "NamespaceFooter.H"
#endif
