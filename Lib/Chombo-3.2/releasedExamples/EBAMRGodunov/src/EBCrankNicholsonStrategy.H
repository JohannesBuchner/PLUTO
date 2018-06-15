#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCRANKNICHOLSONSTRATEGY_H_
#define _EBCRANKNICHOLSONSTRATEGY_H_

#include "EBImplicitIntegrationStrategy.H"

#include "NamespaceHeader.H"

//! \class EBCrankNicholsonStrategy
//! This base class implements the first-order Backward Euler method for
//! embedded boundary algorithms.
class EBCrankNicholsonStrategy: public EBImplicitIntegrationStrategy
{

  public:

  //! Create a new EBCrankNicholsonStrategy with the desired degree of
  //! implicitness.
  //! \param a_implicitness This determines the implicitness of the scheme.
  //!                       0 corresponds to a Forward Euler method,
  //!                       1 to a Backward Euler method, and 1/2 to the
  //!                       2nd-order Crank-Nicholson method.
  explicit EBCrankNicholsonStrategy(Real a_implicitness = 0.5);

  //! Destructor.
  ~EBCrankNicholsonStrategy();

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

  // Degree of implicitness.
  Real m_implicitness;
};

#include "NamespaceFooter.H"
#endif
