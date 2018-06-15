#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  OldPhysIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _OLDPHYSIBC_H_
#define _OLDPHYSIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "UsingNamespace.H"

/// Physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for a hyperbolic system of PDEs.
 */
class OldPhysIBC
{
public:
  /// Constructor
  /**
   */
  OldPhysIBC();

  /// Destructor
  /**
   */
  virtual ~OldPhysIBC()
  {
  };

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a point to a new OldPhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual OldPhysIBC* new_physIBC() = 0;

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;

  /// Set boundary fluxes
  /**
   */
  virtual void fluxBC(FArrayBox&            a_F,
                      const FArrayBox&      a_W,
                      const FArrayBox&      a_Wextrap,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time) = 0;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  virtual void setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time) = 0;

  /// Apply artificial viscosity at the boundary
  /**
     The default implementation of this currently does nothing.
   */
  virtual void artViscBC(FArrayBox&       a_F,
                         const FArrayBox& a_U,
                         const FArrayBox& a_divVel,
                         const int&       a_dir,
                         const Real&      a_time);

  /// set current time
  virtual void time(Real a_time);

  ///
  virtual Real time() const;


protected:
  /// define() has been called
  bool          m_isDefined;

  // The current level's problem domain
  ProblemDomain m_domain;

  /// The current level's grid spacing
  Real          m_dx;

  /// current solution time
  Real m_time;


private:
  // Disallowed for all the usual reasons
  void operator=(const OldPhysIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  OldPhysIBC(const OldPhysIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
