#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INFLOWOUTFLOWADVECTBC_H_
#define _INFLOWOUTFLOWADVECTBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"
#include "PoiseuilleInflowBCValue.H"

#include "NamespaceHeader.H"

///for velocity boundary conditions
/**
   sets = inflow vel at inflow face.
   sets = zero at noflow faces.
   sets = extrapval at outflow face
*/
class InflowOutflowAdvectBC : public EBPhysIBC
{
public:
  virtual ~InflowOutflowAdvectBC()
  {;}

  ///
  InflowOutflowAdvectBC(int a_flowDir,
                        Real a_inflowVel,
                        int a_velComp,
                        bool a_doPoiseInflow,
                        IntVect a_doSlipWallsHi,
                        IntVect a_doSlipWallsLo,
                        RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                        bool a_doWomersleyInflow=false)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doPoiseInflow = a_doPoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_doWomersleyInflow = a_doWomersleyInflow;
    m_velComp   = a_velComp;
    m_poiseInflowFunc = a_poiseBCValue;
    m_isDefined = false;
  }

  ///
  virtual void define(const ProblemDomain&  a_domain,
                      const RealVect&       a_dx)
  {
    m_domain    = a_domain;
    m_dx        = a_dx;
    // if (m_doPoiseInflow)
    //   {
    //     poiseuilleDefine();
    //   }
    m_isDefined = true;
  }

  void poiseuilleDefine()
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
    Real poiseuilleRadius = (Real(m_domain.size(tanDir)))*m_dx[tanDir]/2.;
    RealVect poiseuilleCenter = BASISREALV(tanDir)*poiseuilleRadius;
    if (m_doSlipWallsHi[tanDir]==1 && m_doSlipWallsLo[tanDir]==0)
      {
        poiseuilleRadius = (Real(m_domain.size(tanDir)))*m_dx[tanDir];
        poiseuilleCenter = BASISREALV(tanDir)*poiseuilleRadius;
      }
    else if (m_doSlipWallsHi[tanDir]==0 && m_doSlipWallsLo[tanDir]==1)
      {
        poiseuilleRadius = (Real(m_domain.size(tanDir)))*m_dx[tanDir];
        poiseuilleCenter = RealVect::Zero;
      }
    else if (m_doSlipWallsHi[tanDir]==m_doSlipWallsLo[tanDir] && m_doSlipWallsLo[tanDir]==1)
      {
        MayDay::Error("InflowOutflowAdvectBC.H -- slip walls on both hi and lo sides for Poiseuille...check inputs");
      }
    Real maxVelFactor = 1.5;//planar
//     Real maxVelFactor = 2.0;//axisymmetric tube
    m_poiseInflowFunc = RefCountedPtr<PoiseuilleInflowBCValue>(new PoiseuilleInflowBCValue(poiseuilleCenter, poiseuilleAxis, poiseuilleRadius, maxVelFactor*m_inflowVel, m_flowDir));
  }

  ///  For every box in this level, this function is called.
  virtual void fluxBC(EBFluxFAB&            a_flux,
                      const EBCellFAB&      a_Wcenter,
                      const EBCellFAB&      a_Wextrap,
                      const Side::LoHiSide& a_sd,
                      const Real&           a_time,
                      const EBISBox&        a_ebisBox,
                      const DataIndex&      a_dit,
                      const Box&            a_box,
                      const Box&            a_faceBox,
                      const int&            a_dir);

  /// Initialize---no op here
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const
  {
    MayDay::Error("should not be called");
  }

  /// Set boundary slopes --no op here
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir)
  {;}

protected:
  bool           m_isDefined;
  int            m_velComp;
  ProblemDomain  m_domain;
  RealVect       m_dx;
  int m_flowDir;
  Real m_inflowVel;
  bool m_doPoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  bool m_doWomersleyInflow;
  RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;

  //weak construction is bad but I will allow copy construction and assignment
  //because there are no pointers in the internal data
  InflowOutflowAdvectBC()
  {
    MayDay::Error("invalid operator");
  }

};


///
/**
 */
class InflowOutflowAdvectBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~InflowOutflowAdvectBCFactory()
  {;}

  ///
  InflowOutflowAdvectBCFactory(int a_flowDir,
                               Real a_inflowVel,
                               int a_velComp,
                               bool a_doPoiseInflow,
                               IntVect a_doSlipWallsHi,
                               IntVect a_doSlipWallsLo,
                               RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                               bool a_doWomersleyInflow=false)
  {
    m_velComp = a_velComp;
    m_flowDir = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doPoiseInflow = a_doPoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_poiseInflowFunc = a_poiseBCValue;
    m_doWomersleyInflow = a_doWomersleyInflow;
  }

  ///
  EBPhysIBC* create() const
  {
    InflowOutflowAdvectBC* retval =  new InflowOutflowAdvectBC(m_flowDir,
                                                               m_inflowVel,
                                                               m_velComp,
                                                               m_doPoiseInflow,
                                                               m_doSlipWallsHi,
                                                               m_doSlipWallsLo,
                                                               m_poiseInflowFunc,
                                                               m_doWomersleyInflow);
    return static_cast<EBPhysIBC*>(retval);
  }


protected:

  int  m_velComp;
  int  m_flowDir;
  Real m_inflowVel;
  bool m_doPoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  bool m_doWomersleyInflow;
  RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;

  //weak construction is bad but I will allow copy construction and assignment
  InflowOutflowAdvectBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
