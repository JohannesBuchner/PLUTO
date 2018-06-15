#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RAMPIBC_H_
#define _RAMPIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "PhysIBC.H"

#include "NamespaceHeader.H"

/// Example of a wedge.
/**
   Parameters:
     a_smallPressure - Lower limit for pressure (returned)
     a_gamma         - Gamma for polytropic, gamma-law gas
     a_alpha         - The angle from horizontal to the wedge edge
     a_ms            - Mach shock number
     a_xcorner       - Position of the corner/tip of the wedge
     a_artvisc       - Artificial viscosity coefficient

   This entire problem is rotated so that the wedge edge is horizontal
   in the problem domain.
*/
class RampIBC: public PhysIBC
{
public:
  /// Null Constructor
  /**
   */
  RampIBC();

  /// Constructor which defines parameters used by Fortran routines
  /**
   */
  RampIBC(Real&       a_smallPressure,
          const Real& a_gamma,
          const Real& a_alpha,
          const Real& a_ms,
          const Real& a_xcorner,
          const Real& a_artvisc);

  /// Destructor
  /**
   */
  virtual ~RampIBC();

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon(Real&       a_smallPressure,
                        const Real& a_gamma,
                        const Real& a_alpha,
                        const Real& a_ms,
                        const Real& a_xcorner,
                        const Real& a_artvisc);

  /// Set the flag m_isFortranCommonSet to true
  /**
     Set the flag m_isFortranCommonSet to true so that new IBCs made with
     new_physIBC() will have this flag set without calling setFortranCommon()
     (this is a clumsy design and should be improved).
   */
  void setFortranCommonSet();

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
  virtual void initialize(LevelData<FArrayBox>& a_U);

  /// Set boundary primitive values.
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
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
  virtual
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
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
