#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _BASEFUNCTIONHOLDER_H_
#define _BASEFUNCTIONHOLDER_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "EBIBC.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "RefCountedPtr.H"
#include "EBIBC.H"
#include "EBIBCFactory.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"

#include "NamespaceHeader.H"

///
/**
   base class to hold all the various functions on the domain surface
 */
class BaseFunctionHolder
{
public:

  BaseFunctionHolder()
  {;}

  virtual ~BaseFunctionHolder()
  {;}

  ///boundary conditions
  virtual RefCountedPtr<BaseBCValue> getPressureBCFunction(int a_idir, Side::LoHiSide a_side) = 0;

  ///boundary conditions
  virtual RefCountedPtr<BaseBCValue> getVelocityBCFunction(int a_idir, Side::LoHiSide a_side, int a_velComp) = 0;

  ///boundary conditions
  virtual RefCountedPtr<BaseBCValue> getVelocityEBBCFunction(int a_velComp) = 0;

  ///boundary conditions
  virtual RefCountedPtr<BaseBCValue> getPressureEBBCFunction() = 0;

  ///boundary conditions
  virtual RefCountedPtr<BaseBCValue> getScalarEBBCFunction(int a_scaComp) = 0;

  /// boundary conditions
  virtual RefCountedPtr<BaseBCValue> getScalarBCFunction(int a_idir, Side::LoHiSide a_side, int a_scaComp) = 0;

  ///initial conditions
  virtual RefCountedPtr<BaseBCValue> getScalarInteriorFunction(int a_scaComp) = 0;

  ///initial conditions
  virtual RefCountedPtr<BaseBCValue> getVelocityInteriorFunction(int a_velComp) = 0;

  ///initial conditions
  virtual RefCountedPtr<BaseBCValue> getPressureInteriorFunction() = 0;

};

#include "NamespaceFooter.H"
#endif
