#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VELBCHOLDER_H_
#define _VELBCHOLDER_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "Tuple.H"

#include "UsingNamespace.H"

// petermc, 14 nov 2007

/// this is a physical BC class designed to handle velocities
class VelBCHolder
/** this class is handles velocities, which are by their nature
    multicomponent */
{
public:

  ///
  VelBCHolder();

  ///
  VelBCHolder(const Tuple<BCHolder, SpaceDim>& a_componentBC);

  ///
  virtual ~VelBCHolder();

  ///
  virtual void applyBCs(LevelData<FArrayBox>& a_state,
                        const DisjointBoxLayout& a_valid,
                        const ProblemDomain& a_domain,
                        const Real a_dx,
                        bool a_homogeneous);

protected:

  ///
  Tuple<BCHolder, SpaceDim> m_componentBC;

private:

};


/// this is a physical BC class designed to handle velocities on edges
class EdgeVelBCHolder
/** this class is handles velocities, which are by their nature
    multicomponent */
{
public:

  ///
  EdgeVelBCHolder();

  ///
  EdgeVelBCHolder(const Tuple<BCHolder, SpaceDim>& a_componentBC);

  ///
  virtual ~EdgeVelBCHolder();

  ///
  virtual void applyBCs(LevelData<FluxBox>& a_state,
                        const DisjointBoxLayout& a_valid,
                        const ProblemDomain& a_domain,
                        const Real a_dx,
                        bool a_homogeneous);

protected:

  ///
  Tuple<BCHolder, SpaceDim> m_componentBC;

private:

};

#endif
