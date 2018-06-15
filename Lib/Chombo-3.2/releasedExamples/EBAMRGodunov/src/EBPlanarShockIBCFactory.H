#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _EBPLANARSHOCKIBCFACTORY_H_
#define _EBPLANARSHOCKIBCFACTORY_H_

#include "EBPlanarShockIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "EBPlanarShockIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///
/**
 */
class EBPlanarShockIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~EBPlanarShockIBCFactory();

  ///
  EBPlanarShockIBCFactory(const Real&     a_gamma,
                          const Real&     a_ms,
                          const Real&     a_center,
                          const int&      a_shocknorm,
                          const bool&     a_shockbackward,
                          const bool&     a_doRZCoords);

  ///
  EBPhysIBC* create() const;

protected:
  Real m_gamma;
  Real m_ms;
  Real m_center;
  int  m_shocknorm;
  bool m_shockbackward;
  bool m_doRZCoords;
private:
  EBPlanarShockIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const EBPlanarShockIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  EBPlanarShockIBCFactory(const EBPlanarShockIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
