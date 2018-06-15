#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _INFLOWOUTFLOWIBC_H_
#define _INFLOWOUTFLOWIBC_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "Vector.H"
#include "EBIBC.H"
#include "EBIBCFactory.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "ExtrapAdvectBC.H"
#include "PoiseuilleInflowBCValue.H"

#include "NamespaceHeader.H"

///
/**
   This class is meant to be a server for initial boundary conditions
   for all the various stages of BCG INS.
*/
class InflowOutflowIBC: public EBIBC
{
public:
  ///
  /**
     Create Poiseuille flow in a tube
  */
  InflowOutflowIBC(int a_flowDir,
                   Real a_inflowVel,
                   int a_orderEBBC,
                   IntVect a_doSlipWallsHi,
                   IntVect a_doSlipWallsLo,
                   bool a_doPoiseInflow,
                   bool a_initPoiseData = false,
                   RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                   bool a_doWomersleyInflow = false)
  {
    m_flowDir       = a_flowDir;
    m_inflowVel     = a_inflowVel;
    m_orderEBBC     = a_orderEBBC;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_doPoiseInflow = a_doPoiseInflow;
    m_initPoiseData = a_initPoiseData;
    if (m_initPoiseData)
    {
      CH_assert(m_doPoiseInflow);
    }
    m_poiseInflowFunc = a_poiseBCValue;
    m_isPoiseDefined = false;
    m_doWomersleyInflow = a_doWomersleyInflow;
  }

  ///
  virtual ~InflowOutflowIBC()
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
  void initializePressureGradient(LevelData<EBCellFAB>&    a_gradient,
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
  virtual RefCountedPtr<EBPhysIBCFactory> getVelAdvectBC(int a_velComp) const ;


  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
  */
  virtual RefCountedPtr<EBPhysIBCFactory> getScalarAdvectBC(const int&  a_comp) const ;


  ///
  /**
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getVelBC(int a_icomp) const ;

  ///
  /**
     Return velocity boundary conditions for embedded boundary.
  */
  virtual RefCountedPtr<BaseEBBCFactory> getVelocityEBBC(int a_velComp) const ;

  virtual RefCountedPtr<BaseEBBCFactory> getScalarEBBC() const
  {
    MayDay::Error("default and invalid implementaion of getScalarEBBC");
    //code to make compilers shut up
    return RefCountedPtr<BaseEBBCFactory>();
  }

  virtual RefCountedPtr<BaseEBBCFactory> getPressureEBBC() const ;

  void poiseuilleDefine(const ProblemDomain&  a_domain,
                        const RealVect&       a_dx) const
  {
    CH_assert(SpaceDim==2);
    RealVect poiseuilleAxis = BASISREALV(m_flowDir);
    int tanDir=1;
    if (m_flowDir==1)
    {
      tanDir=0;
    }
//     Tuple<int,SpaceDim-1> tanDirs = PolyGeom::computeTanDirs(idir);
//     for (int itan = 0; itan < SpaceDim-1; itan++)
//       {
//         int tanDir = tanDirs[itan];
//       }
    Real poiseuilleRadius = (Real(a_domain.size(tanDir)))*a_dx[tanDir]/2.;
    RealVect poiseuilleCenter = BASISREALV(tanDir)*poiseuilleRadius;
    if (m_doSlipWallsHi[tanDir]==1 && m_doSlipWallsLo[tanDir]==0)
      {
        poiseuilleRadius = (Real(a_domain.size(tanDir)))*a_dx[tanDir];
        poiseuilleCenter = BASISREALV(tanDir)*poiseuilleRadius;
      }
    else if (m_doSlipWallsHi[tanDir]==0 && m_doSlipWallsLo[tanDir]==1)
      {
        poiseuilleRadius = (Real(a_domain.size(tanDir)))*a_dx[tanDir];
        poiseuilleCenter = RealVect::Zero;
      }
    else if (m_doSlipWallsHi[tanDir]==m_doSlipWallsLo[tanDir] && m_doSlipWallsLo[tanDir]==1)
      {
        MayDay::Error("InflowOutflowIBC.H -- slip walls on both hi and lo sides for Poiseuille...check inputs");
      }
    Real maxVelFactor = 1.5;//planar
//     Real maxVelFactor = 2.0;//axisymmetric tube
    m_poiseInflowFunc = RefCountedPtr<PoiseuilleInflowBCValue>(new PoiseuilleInflowBCValue(poiseuilleCenter, poiseuilleAxis, poiseuilleRadius, maxVelFactor*m_inflowVel, m_flowDir));

    m_isPoiseDefined = true;
  }

// protected:
  int     m_flowDir;
  Real    m_inflowVel;
  int     m_orderEBBC;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  bool    m_doPoiseInflow;
  bool    m_initPoiseData;
  bool    m_doWomersleyInflow;
  mutable bool    m_isPoiseDefined;
  mutable RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;

  /// weak construction is bad, mkay?
  InflowOutflowIBC()
  {
    MayDay::Error("invalid operator");
  }
};

///
/**
 */
class InflowOutflowIBCFactory: public EBIBCFactory
{
public:
  ///
  /**
   */
  InflowOutflowIBCFactory(int a_flowDir,
                          Real a_inflowVel,
                          int a_orderEBBC,
                          IntVect a_doSlipWallsHi=IntVect::Zero,
                          IntVect a_doSlipWallsLo=IntVect::Zero,
                          bool a_doPoiseInflow=false,
                          bool a_initPoiseData=false,
                          RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                          // const RealVect& a_centerPt,
                          // const RealVect& a_tubeAxis,
                          // const Real&     a_tubeRadius,
                          // const Real&     a_maxVel,
                          // const int&      a_velComp,
                          bool a_doWomersleyInflow=false)
  {
    m_flowDir       = a_flowDir;
    m_inflowVel     = a_inflowVel;
    m_orderEBBC     = a_orderEBBC;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_doPoiseInflow = a_doPoiseInflow;
    m_initPoiseData = a_initPoiseData;
    m_poiseInflowFunc = a_poiseBCValue;
    m_doWomersleyInflow = a_doWomersleyInflow;
  }

  ///
  virtual ~InflowOutflowIBCFactory()
  {
  }

  ///
  /**
   */
  virtual EBIBC* create()  const
  {
    InflowOutflowIBC* retDerived = new InflowOutflowIBC(m_flowDir,
                                                        m_inflowVel,
                                                        m_orderEBBC,
                                                        m_doSlipWallsHi,
                                                        m_doSlipWallsLo,
                                                        m_doPoiseInflow,
                                                        m_initPoiseData,
                                                        m_poiseInflowFunc,
                                                        m_doWomersleyInflow);
    EBIBC* retBase =  (EBIBC*)retDerived;
    return retBase;
  }

protected:
  int     m_flowDir;
  Real    m_inflowVel;
  int     m_orderEBBC;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  bool    m_doPoiseInflow;
  bool    m_initPoiseData;
  bool    m_doWomersleyInflow;
  RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;

  /// weak construction is bad but i will allow copy construction and assignment
  InflowOutflowIBCFactory()
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
