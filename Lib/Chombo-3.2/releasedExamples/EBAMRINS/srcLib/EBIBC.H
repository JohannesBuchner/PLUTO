#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBIBC_H_
#define _EBIBC_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBCFactory.H"
#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseEBBC.H"

#include "NamespaceHeader.H"

///
/**
   This class is meant to be  a server for initial boundary conditions
   for all the various stages of BCG INS.


   This is meant to be a thumbnail sketch of how boundary conditions are
   handled in EBAMRINS.   All initial and boundary conditions are handled by a server
   class called EBIBC (Chombo/lib/src/EBIBC.H).  This class has the
   following functions that need to be filled.   The ones that concern
   a scalar are not included as they are no longer called.  EBAMRINS used
   to advect a scalar.
   \begin{itemize}

   \item initializeVelocity--- Fill the velocity field at time = 0.
   \item initializePressure--- Fill the pressure at time= 0.  There is an
   iteration to make the pressure and the velocity consistent but this
   fills the first guess at that initialization.
   \item getPressBC---Pressure domain boundary conditions.   This is what
   gets used in the projection solve.
   \item getPressureEBBC---EB pressure boundary conditions for the
   projection solve.
   \item getVelBC----Velocity boundary conditions.  This gets used in the
   viscous TGA solve
   \item getVelocityEBBC---EB boundary conditions for the viscous solve.
   \item getMACVelBC---This sets the velocity at faces for the MAC
   projection.
   \item  getVelAdvectBC---boundary conditions for the Godunov advection.

   \end{itemize}

   A good example to look at is InflowOutflowIBC.H (also in
   Chombo/lib/src/EBAMRTimeDependent).
   dtg
 */
class EBIBC
{
public:
  ///
  EBIBC();

  ///
  virtual ~EBIBC();

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
                                  const RealVect&          a_dx) const = 0;

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
                                  const RealVect&          a_dx) const = 0;
  ///
  /**
     Fill input data holder with stress at time = 0.
  */
  virtual void initializeStress(LevelData<EBCellFAB>&    a_stress,
                                const Real&              a_initstress,
                                const Real&              a_inflowstress,
                                const DisjointBoxLayout& a_grids,
                                const EBISLayout&        a_ebisl,
                                const ProblemDomain&     a_domain,
                                const RealVect&          a_origin,
                                const Real&              a_time,
                                const RealVect&          a_dx,
                                const int&               a_comp) const
  {
    //default implementation for those that do not need it
    MayDay::Error("initialize stress not defined");
  }
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
                                  const RealVect&          a_dx) const = 0;

  ///
  /**
     Return pressure boundary conditions for domain.
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getScalarBC() const
  {
    MayDay::Error("default and invalid implementaion of getScalarEBBC");
    //code to make compilers shut up
    return RefCountedPtr<BaseDomainBCFactory>();
  }

  ///
  /**
     Return pressure boundary conditions for domain.
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getPressBC() const = 0;

  ///
  /**
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getVelBC(int a_icomp) const = 0;

  ///
  /**
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getMACVelBC() const = 0;

  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
   */
  virtual  RefCountedPtr<EBPhysIBCFactory> getVelAdvectBC(int a_velComp) const = 0;


  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
   */
  virtual  RefCountedPtr<EBPhysIBCFactory> getScalarAdvectBC(const int&  a_comp) const  = 0;

  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
   */
  virtual  RefCountedPtr<EBPhysIBCFactory> getStressAdvectBC(const int&  a_comp, Real a_inflowstress)
  {
    return getScalarAdvectBC(a_comp);
  }



  ///
  /**
     Return velocity boundary conditions for embedded boundary.
   */
  virtual RefCountedPtr<BaseEBBCFactory> getVelocityEBBC(int a_velComp) const = 0;


  ///
  /**
     Return pressure boundary conditions for embedded boundary.
     this gets a default implementation because not everyone uses it.
   */
  virtual RefCountedPtr<BaseEBBCFactory> getPressureEBBC() const = 0;


  ///
  /**
     Return scalar boundary conditions for embedded boundary.
     this gets a default implementation because not everyone uses it.
   */
  virtual RefCountedPtr<BaseEBBCFactory> getScalarEBBC() const
  {
    MayDay::Error("default and invalid implementaion of getScalarEBBC");
    //code to make compilers shut up
    return RefCountedPtr<BaseEBBCFactory>();
  }
};

#include "NamespaceFooter.H"
#endif
