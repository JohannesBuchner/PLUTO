#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRTESTCOMMON_H_
#define _EBAMRTESTCOMMON_H_

#include "IntVect.H"
#include "IntVectSet.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "ParmParse.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "MultiSphereIF.H"
#include "EBTensorCFInterp.H"
#include "EBAMRIO.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"
/**********/
///
/**
 */
class BaseAnalytic
{
public:
  ///
  /**
   */
  virtual void veloVal(Real  a_vel[SpaceDim],           const RealVect& a_point) = 0;

  ///
  /**
   */
  virtual void gradVal(Real a_grad[SpaceDim*SpaceDim],  const RealVect& a_point) = 0;

  BaseAnalytic() 
  {;}
  
  virtual ~BaseAnalytic() 
  {;}
};

///
/**
 */
class LinearAnalytic: public BaseAnalytic
{
public:
  ///
  /**
   */
  virtual   void  veloVal(Real  a_vel[SpaceDim],           const RealVect& a_point)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        a_vel[idir] = 0;
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            a_vel[idir] += m_coeff[jdir]*a_point[jdir];
          }
      }
  }

  ///
  /**
   */
  virtual void gradVal(Real a_grad[SpaceDim*SpaceDim], const RealVect& a_point)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            int iindex = TensorCFInterp::gradIndex(idir, jdir);
            a_grad[iindex] = m_coeff[jdir];
          }
      }

  }
  LinearAnalytic() 
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        m_coeff[idir] =  Real(idir + 1);
      }
  }
  virtual ~LinearAnalytic() 
  {;}

  RealVect m_coeff;
};

///
/**
 */
class ConstAnalytic: public BaseAnalytic
{
public:
  ///
  /**
   */
  virtual  void veloVal(Real  a_vel[SpaceDim],           const RealVect& a_point)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            a_vel[idir] = m_coeff[jdir];
          }
      }
  }

  ///
  /**
   */
  virtual void gradVal(Real a_grad[SpaceDim*SpaceDim], const RealVect& a_point)
  {
    for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            int iindex = TensorCFInterp::gradIndex(ivar, jdir);
            a_grad[iindex] = 0.;
          }
      }

  }
  ConstAnalytic() 
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        m_coeff[idir] =  Real(idir + 1);
      }
  }
  virtual ~ConstAnalytic() 
  {;}

  RealVect m_coeff;
};


///
/**
 */
class QuadAnalytic: public BaseAnalytic
{
public:
  ///
  /**
   */
  virtual void veloVal(Real  a_vel[SpaceDim],           const RealVect& a_point)
  {
    for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        a_vel[ivar] = 0;
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            a_vel[ivar] += m_coeff[jdir]*a_point[jdir]*a_point[jdir];
          }
      }
  }

  ///
  /**
   */
  virtual void gradVal(Real a_grad[SpaceDim*SpaceDim], const RealVect& a_point)
  {
    for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            int iindex = TensorCFInterp::gradIndex(ivar, jdir);
            a_grad[iindex] = 2.*m_coeff[jdir]*a_point[jdir];
          }
      }

  }

  ///
  /**
   */
  QuadAnalytic() 
  {
    for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        m_coeff[ivar] =  Real(ivar + 1);
      }
  }
  virtual ~QuadAnalytic() 
  {;}

  RealVect m_coeff;
};
/**********/
class CubeAnalytic: public BaseAnalytic
{
public:
  ///
  /**
   */
  virtual void veloVal(Real  a_vel[SpaceDim],           const RealVect& a_point)
  {
    for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            a_vel[ivar] += m_coeff[jdir]*a_point[jdir]*a_point[jdir]*a_point[jdir];
          }
      }
  }

  ///
  /**
   */
  virtual void gradVal(Real a_grad[SpaceDim*SpaceDim], const RealVect& a_point)
  {
    for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            int iindex = TensorCFInterp::gradIndex(ivar, jdir);
            a_grad[iindex] = 3.*m_coeff[jdir]*a_point[jdir]*a_point[jdir];
          }
      }

  }

  CubeAnalytic() 
  {
    for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        m_coeff[ivar] =  Real(ivar + 1);
      }
  }
  virtual ~CubeAnalytic() 
  {;}

  RealVect m_coeff;
};

///
/**
 */
class EBAMRTestCommon
{
public:
  ///
  /**
   */
  static int compareError(const LevelData<EBCellFAB>& a_errorFine,
                          const EBISLayout&           a_ebislFine,
                          const DisjointBoxLayout&    a_gridsFine,
                          const Box&                  a_domainFine,
                          const LevelData<EBCellFAB>& a_errorCoar,
                          const EBISLayout&           a_ebislCoar,
                          const DisjointBoxLayout&    a_gridsCoar,
                          const Box&                  a_domainCoar);

  ///
  /**
   */
  static int makeEBISL(EBISLayout& a_ebisl,
                       const DisjointBoxLayout& a_grids,
                       const Box& a_domain,
                       const int& a_nghost);

  ///
  /**
   */
  static int
  makeHalfLayouts(DisjointBoxLayout& a_gridsFine,
                  EBISLayout&        a_ebislFine,
                  DisjointBoxLayout& a_gridsCoar,
                  EBISLayout&        a_ebislCoar,
                  const Box& a_domainFine,
                  const Box& a_domainCoar,
                  int a_refrat, int a_halfDir);

  ///
  /**
   */
  static int makeGeometry(Box& a_domain,
                          RealVect& a_dx);

  ///
  /**
   */
  static int 
  checkForZero(const LevelData<EBCellFAB>& a_errorVelo,
               const DisjointBoxLayout&    a_gridsFine,
               const EBISLayout&           a_ebislFine,
               const Box&                  a_domainFine,
               string a_funcname);


  ///
  /**
   */
  static void 
  outputError(const LevelData<EBCellFAB>&   a_errorFine,
              const LevelData<EBCellFAB>&   a_errorCoar,
              const DisjointBoxLayout&      a_gridsFine,
              const DisjointBoxLayout&      a_gridsCoar,
              const Box&                    a_domainFine,
              const Box&                    a_domainCoar,
              const string& a_fileFine,
              const string& a_fileCoar);
};

#include "NamespaceFooter.H"

#endif
