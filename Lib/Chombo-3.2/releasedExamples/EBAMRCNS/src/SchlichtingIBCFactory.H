#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _SchlichtingIBCFACTORY_H
#define _SchlichtingIBCFACTORY_H

#include "SchlichtingIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "EBAMRCNSParams.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "SchlichtingIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///
/**
 */
class SchlichtingIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~SchlichtingIBCFactory();

  ///
  SchlichtingIBCFactory(const SchlichtingParams&  a_params);


  ///
  EBPhysIBC* create() const;

protected:
  SchlichtingParams  m_params;

private:
  SchlichtingIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const SchlichtingIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  SchlichtingIBCFactory(const SchlichtingIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
