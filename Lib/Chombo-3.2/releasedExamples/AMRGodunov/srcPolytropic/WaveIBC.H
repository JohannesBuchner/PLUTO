#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _WAVEIBC_H_
#define _WAVEIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"
#include "PhysIBC.H"

#include "SolidBCF_F.H"

#include "NamespaceHeader.H"

/// Example of a sinusoidal plane wave plus a constant density
/**
   Parameters:
     a_smallPressure  - Lower limit for pressure (returned)
     a_gamma          - Gamma for polytropic, gamma-law gas
     a_ambientDensity - Ambient density add to the plane-wave
     a_deltaDensity   - Mean of the plane-wave
     a_pressure       - If 0, use isentropic pressure
                        if 1, use the constant pressure of 1.0
     a_waveNumber     - The wave number of the plane-wave
     a_velocity       - Initial velocity of the gas
     a_center         - Position of a maximum of the plane-wave
     a_artvisc        - Artificial viscosity coefficient
 */
class WaveIBC: public PhysIBC
{
public:
  /// Null constructor
  /**
   */
  WaveIBC();

  /// Constructor which defines parameters used by Fortran routines
  /**
   */
  WaveIBC(Real&           a_smallPressure,
          const Real&     a_gamma,
          const Real&     a_ambientDensity,
          const Real&     a_deltaDensity,
          const int&      a_pressure,
          const IntVect&  a_waveNumber,
          const RealVect& a_center,
          const RealVect& a_velocity,
          const Real&     a_artvisc);

  /// Destructor
  /**
   */
  virtual ~WaveIBC();

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon(Real&           a_smallPressure,
                        const Real&     a_gamma,
                        const Real&     a_ambientDensity,
                        const Real&     a_deltaDensity,
                        const int&      a_pressure,
                        const IntVect&  a_waveNumber,
                        const RealVect& a_center,
                        const RealVect& a_velocity,
                        const Real&     a_artvisc);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  PhysIBC *new_physIBC();

  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_U);

  /// Set boundary primitive values.
  /**
   */
  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);


protected:
  // True if the Fortran common block has been set
  bool m_isFortranCommonSet;
};

#include "NamespaceFooter.H"

#endif
