#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _MODIANOIBCFACTORY_H_
#define _MODIANOIBCFACTORY_H_

#include "ModianoIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "ModianoIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///
/**
 */
class ModianoIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~ModianoIBCFactory();

  ///
  ModianoIBCFactory(const Real&     a_gamma,
                    const Real&     a_waveAmp,
                    const Real&     a_waveWidth,
                    const RealVect& a_center,
                    const RealVect& a_waveDir,
                    bool   a_freeStreamOnly,
                    bool   a_useNegWave = false );

  ///
  EBPhysIBC* create() const;

protected:
  Real     m_gamma;
  Real     m_waveAmp;
  Real     m_waveWidth;
  RealVect m_waveDir;
  RealVect m_center;
  bool     m_freeStreamOnly;
  bool     m_useNegWave;

private:
  ModianoIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const ModianoIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  ModianoIBCFactory(const ModianoIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
