#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
//  PhysIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _PHYSIBC_H_
#define _PHYSIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for a hyperbolic system of PDEs.
 */
class PhysIBC
{
public:
  /// Constructor
  /**
   */
  PhysIBC();

  /// Destructor
  /**
   */
  virtual ~PhysIBC();

  /// Define the object
  /**
   */
  virtual void define(/// problem domain index space
                      const ProblemDomain& a_domain,
                      /// grid spacing
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual PhysIBC* new_physIBC() = 0;

  /// Set up initial conditions
  /**
   */
  virtual void initialize(/// conserved variables
                          LevelData<FArrayBox>& a_U) = 0;

  /// Set boundary primitive values
  /**
   */
  virtual void primBC(/// primitive variables at face centers; this routine replaces values along boundary faces of domain
                      FArrayBox&            a_WGdnv,
                      /// extrapolated primitive variables to a_side of cells in direction a_dir (cell-centered data)
                      const FArrayBox&      a_Wextrap,
                      /// primitive variables at start of time step (cell-centered data)
                      const FArrayBox&      a_W,
                      /// normal direction of the domain where boundary condition fluxes needed
                      const int&            a_dir,
                      /// side of the domain where boundary condition fluxes needed
                      const Side::LoHiSide& a_side,
                      /// physical time of the problem, for time-varying boundary conditions
                      const Real&           a_time) = 0;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one-sided difference
     approximations.  If this function doesn't change them, they will be
     used for the slopes at the boundaries.
   */
  virtual
  void setBdrySlopes(/// face-centered differences of primitive variables
                     FArrayBox&       a_dW,
                     /// primitive variables at start of time step (cell-centered data)
                     const FArrayBox& a_W,
                     /// normal direction of faces
                     const int&       a_dir,
                     /// physical time of the problem, for time-varying boundary conditions
                     const Real&      a_time) = 0;

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(/// face-centered flux, to be adjusted on boundary
                 FArrayBox&       a_F,
                 /// cell-centered conserved variables
                 const FArrayBox& a_U,
                 /// face-centered divergence of cell-centered velocity
                 const FArrayBox& a_divVel,
                 /// normal direction of faces
                 const int&       a_dir,
                 /// physical time of the problem, for time-varying boundary conditions
                 const Real&      a_time) = 0;

  /// This function is called by primBC() to get boundary faces of a Box.
  virtual
  void getBoundaryFaces(/// face-centered box of boundary faces to fill in; subbox of a_dataFaceBox
                        Box&                   a_boundaryBox,
                        /// entire face-centered box
                        const Box&             a_dataFaceBox,
                        /// normal direction of faces
                        const int&             a_dir,
                        /// side of the domain where boundary faces needed
                        const Side::LoHiSide&  a_side);

protected:

  // Has this object been defined
  bool          m_isDefined;

  // The current level's problem domain
  ProblemDomain m_domain;

  // The current level's grid spacing
  Real          m_dx;

private:

  // Disallowed for all the usual reasons
  void operator=(const PhysIBC&);
  PhysIBC(const PhysIBC&);
};

#include "NamespaceFooter.H"
#endif
