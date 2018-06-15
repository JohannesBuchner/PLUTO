#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRNAVIERSTOKESFACTORY_H_
#define _AMRNAVIERSTOKESFACTORY_H_

#include "AMRNavierStokes.H"
#include "AMRLevelFactory.H"
#include "PhysBCUtil.H"

#include "UsingNamespace.H"

class AMRNavierStokesFactory : public AMRLevelFactory
{
public:

  ///
  AMRNavierStokesFactory();

  ///
  virtual AMRLevel* new_amrlevel() const;

  virtual
  ~AMRNavierStokesFactory();

  ///
  void CFL(const Real a_cfl);

  ///
  void refinementThreshold(const Real a_refine_threshold);

  ///
  void setPhysBC(const PhysBCUtil& a_physBC);

  ///
  void setLimitSolverCoarsening(bool a_limitSolverCoarsening);


protected:
  ///
  void setDefaultValues();

  ///
  Real m_cfl;

  ///
  Real m_refine_thresh;

  ///
  PhysBCUtil* m_physBCPtr;

  ///
  bool m_limitSolverCoarsening;
};

#endif
