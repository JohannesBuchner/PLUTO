#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _EBEXPLOSIONIBCFACTORY_H_
#define _EBEXPLOSIONIBCFACTORY_H_

#include "EBExplosionIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "EBExplosionIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///
/**
 */
class EBExplosionIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~EBExplosionIBCFactory();

  ///
  EBExplosionIBCFactory(const Real&     a_gamma,
                        const Real&     a_size,
                        const Real&     a_p0,
                        const Real&     a_r0,
                        const Real&     a_p1,
                        const Real&     a_r1,
                        const RealVect& a_center,
                        const bool&     a_doRZCoords);

  ///
  EBPhysIBC* create() const;

protected:
  bool     m_doRZCoords;
  Real     m_gamma;
  Real     m_size;
  Real     m_p0;
  Real     m_r0;
  Real     m_p1;
  Real     m_r1;
  RealVect m_center;
private:
  EBExplosionIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const EBExplosionIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  EBExplosionIBCFactory(const EBExplosionIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
