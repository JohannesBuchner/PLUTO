#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _POISSELLETUBE_H_
#define _POISSELLETUBE_H_

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
#include "PoisselleTubeBCValue.H"
#include "ExtrapAdvectBC.H"

#include "UsingNamespace.H"

///
/**
   This class is meant to be  a server for initial boundary conditions
   for all the various stages of BCG INS.
*/
class PoisselleTube: public EBIBC
{
public:
  ///
  /**
     Create poisselle flow in a tube
  */
  PoisselleTube(const Vector<PoisselleTubeBCValue>& a_bcval)
  {
    m_bcval = a_bcval;
  }

  ///
  virtual ~PoisselleTube()
  {;}

  ///
  virtual void initializeVelocity(LevelData<EBCellFAB>&    a_velocity,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  ///
  virtual void initializePressure(LevelData<EBCellFAB>&    a_pressure,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  ///
  virtual void initializeScalar ( LevelData<EBCellFAB>&    a_scalar,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

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
  Vector<PoisselleTubeBCValue> m_bcval;

private:
  /// weak construction is badn
  PoisselleTube()
  {
    MayDay::Error("invalid operator");
  }
};

///
/**
 */
class PoisselleTubeFactory: public EBIBCFactory
{
public:
  ///
  /**
   */
  PoisselleTubeFactory(const RealVect& a_centerPt,
                       const RealVect& a_tubeAxis,
                       const Real&     a_tubeRadius,
                       const Real&     a_maxVel)
  {
    m_bcval.resize(SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        m_bcval[idir] = PoisselleTubeBCValue(a_centerPt, a_tubeAxis, a_tubeRadius, a_maxVel, idir);
      }
  }

  ///
  virtual ~PoisselleTubeFactory()
  {
  }

  ///
  /**
   */
  virtual EBIBC* create()  const
  {
    PoisselleTube* retDerived = new PoisselleTube(m_bcval);
    EBIBC* retBase =  (EBIBC*)retDerived;
    return retBase;
  }

protected:
  Vector<PoisselleTubeBCValue> m_bcval;

private:
  /// weak construction is bad but i will allow copy construction and assignment
  PoisselleTubeFactory()
  {
    MayDay::Error("invalid operator");
  }
};

#endif
