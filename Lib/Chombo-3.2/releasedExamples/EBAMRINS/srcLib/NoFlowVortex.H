#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _NOFLOWVORTEX_H_
#define _NOFLOWVORTEX_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "EBIBC.H"
#include "EBIBCFactory.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"

#include "NamespaceHeader.H"

///
/**
   This class is meant to be  a server for initial boundary conditions
   for all the various stages of BCG INS.
 */
class NoFlowVortex: public EBIBC
{
public:
  ///
  /**
     inside core radius vtheta = strength* r
     outise core radius vtheta = 0.
   */
  NoFlowVortex(const RealVect& a_center,
               const Real&     a_coreRadius,
               const Real&     a_strength);

  ///
  virtual ~NoFlowVortex()
  {;}

  ///
  /**
     Fill input data holder with velocity at time = 0.
  */
  virtual void initializeVelocity(LevelData<EBCellFAB>&    a_velocity,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  ///
  /**
     Fill input data holder with pressure at time = 0.
  */
  virtual void initializePressure(LevelData<EBCellFAB>&    a_pressure,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  ///
  /**
     Fill input data holder with scalar at time = 0.
  */
  virtual void initializeScalar ( LevelData<EBCellFAB>&    a_scalar,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  void getRadius(Real&           a_radius,
                 RealVect&       a_xdiff,
                 const RealVect& a_xval) const;

  void getVelPt(RealVect& a_vel, const RealVect& a_xval, Real a_pi) const;


  void getScalarPt(Real& a_scal, const RealVect& a_xval) const;

  void getXVal(RealVect&       a_xval,
               const RealVect& a_origin,
               const VolIndex& a_vof,
               const RealVect& a_dx) const;
  ///
  /**
     Return pressure boundary conditions for domain.
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getPressBC() const ;

  ///
  /**
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getMACVelBC() const ;

  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
   */
  virtual  RefCountedPtr<EBPhysIBCFactory> getVelAdvectBC(int a_velComp) const ;


  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
   */
  virtual  RefCountedPtr<EBPhysIBCFactory> getScalarAdvectBC(const int&  a_comp) const ;


  ///
  /**
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getVelBC(int a_icomp) const ;

  ///
  /**
     Return velocity boundary conditions for embedded boundary.
   */
  virtual RefCountedPtr<BaseEBBCFactory> getVelocityEBBC(int a_velComp) const ;

  virtual RefCountedPtr<BaseEBBCFactory> getPressureEBBC() const ;
protected:
  RealVect m_center;
  Real     m_coreRadius;
  Real     m_strength;
  RefCountedPtr<BaseDomainBCFactory>  m_domBCPress;
  RefCountedPtr<BaseDomainBCFactory>  m_domBCMACVel;
private:
  /// weak construction
  NoFlowVortex()
  {
    MayDay::Error("invalid operator");
  }
  //we disallow copy construction and assignment for classes with pointered data
};

///
/**
 */
class NoFlowVortexFactory: public EBIBCFactory
{
public:

  ///
  /**
     inside core radius vtheta = strength* r
     outise core radius vtheta = 0.
   */
  NoFlowVortexFactory(const RealVect& a_center,
                         const Real&     a_coreRadius,
                         const Real&     a_strength)
  {
    m_center = a_center;
    m_coreRadius = a_coreRadius;
    m_strength = a_strength;
  }

  ///
  virtual ~NoFlowVortexFactory()
  {
  }

  ///
  /**
   */
  virtual EBIBC* create()  const
  {
    NoFlowVortex* retDerived = new NoFlowVortex(m_center,m_coreRadius, m_strength);
    EBIBC* retBase =  (EBIBC*)retDerived;
    return retBase;
  }

protected:
  RealVect m_center;
  Real     m_coreRadius;
  Real     m_strength;
private:
  /// weak construction is bad but i will allow copy construction and assignment
  NoFlowVortexFactory()
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
