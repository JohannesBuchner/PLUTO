#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBTGASTRATEGY_H__
#define _EBTGASTRATEGY_H__

#include "EBImplicitIntegrationStrategy.H"

#include "NamespaceHeader.H"

//! \class EBTGAStrategy
//! This base class implements the first-order Backward Euler method for
//! embedded boundary algorithms.
class EBTGAStrategy: public EBImplicitIntegrationStrategy
{

  public:

  //! Create a new EBTGAStrategy with the desired degree of
  //! implicitness.
  EBTGAStrategy();

  //! Destructor.
  ~EBTGAStrategy();

  void updateSolution(LevelData<EBCellFAB>& a_phiNew,
                      LevelData<EBCellFAB>& a_phiOld,
                      LevelData<EBCellFAB>& a_src,
                      Real a_oldTime,
                      Real a_dt);

  void computeDiffusion(LevelData<EBCellFAB>& a_diffusiveTerm,
                        LevelData<EBCellFAB>& a_phiOld,
                        LevelData<EBCellFAB>& a_src,
                        Real a_oldTime,
                        Real a_dt);

  private:

  Real m_mu1, m_mu2, m_mu3, m_mu4, m_r1;

};

#include "NamespaceFooter.H"
#endif
