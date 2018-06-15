#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeBC.H
// adapted from GhostBC by DTGraves, Mon, July 19, 1999
// petermc, 13 Feb 2001
// petermc, 21 Oct 2002, removed m_inhomogeneous, which can be in
// derived class of FaceNodeBC.

#ifndef _NODEBC_H_
#define _NODEBC_H_

#include "ProblemDomain.H"
#include "NodeFArrayBox.H"
#include "REAL.H"
#include "SPACE.H"
#include "Tuple.H"
#include "Interval.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "NamespaceHeader.H"

/// A class to encapsulate the operations of boundary conditions on a face
/** FaceNodeBC is a class to encapsulate the
    operations of boundary conditions at a face.

    If the solution is phi and the face normal direction is x,
    the boundary conditions usually used can be expressed as<br>
    A*dphi/dx + B*phi = C.<br>
    The derived class provides a function to fill A and B and C.

    These functions are meant to be called by DomainNodeBC only.
*/
class FaceNodeBC
{
public:
  friend class DomainNodeBC;

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor, sets interval to (-1:-1).
   */
  FaceNodeBC()
    :m_components(-1,-1)
  {
  }

  ///
  /** Constructor for face on side <i>a_sd</i> in direction <i>a_idir</i>.
      This constructor sets the components interval to be (0:0).
   */
  FaceNodeBC(int a_dir, Side::LoHiSide a_sd);

  ///
  /** Constructor for face on side <i>a_sd</i> in direction <i>a_idir</i>,
      for data components <i>a_comps</i>.
   */
  FaceNodeBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps);

  ///
  /** Destructor.
   */
  virtual ~FaceNodeBC()
  {
  }

protected:

  ///
  /** Virtual constructor workaround.
      Need this in derived class.
   */
  virtual FaceNodeBC* new_boxBC() const = 0;

  ///
  /** Defines face on side <i>a_sd</i> in direction <i>a_idir</i>.
      Sets the components interval to be (0:0).
   */
  void
  define(int a_dir, Side::LoHiSide a_sd);

  ///
  /** Defines face on side <i>a_sd</i> in direction <i>a_idir</i>,
      for data components <i>a_comps</i>.
   */
  void
  define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps);

  /*@}*/

  /**
     \name Functions to apply boundary conditions
  */
  /*@{*/

  ///
  /** Apply inhomogeneous boundary conditions on this face.
   */
  virtual void
  applyInhomogeneousBCs(/// NODE-centered data to be modified by boundary condition
                        FArrayBox& a_state,
                        /// CELL-centered physical domain
                        const ProblemDomain& a_domain,
                        /// mesh spacing
                        Real a_dx) const;

  ///
  /** Apply inhomogeneous boundary conditions on this face.
   */
  virtual void
  applyInhomogeneousBCs(/// NODE-centered data to be modified by boundary condition
                        FArrayBox& a_state,
                        /// CELL-centered physical domain
                        const Box& a_domain,
                        /// mesh spacing
                        Real a_dx) const;

  ///
  /** Apply homogeneous boundary conditions on this face.
   */
  virtual void
  applyHomogeneousBCs(/// NODE-centered data to be modified by boundary condition
                      FArrayBox& a_state,
                      /// CELL-centered physical domain
                      const ProblemDomain& domain,
                      /// mesh spacing
                      Real a_dx) const;
  ///
  /** Apply homogeneous boundary conditions on this face.
   */
  virtual void
  applyHomogeneousBCs(/// NODE-centered data to be modified by boundary condition
                      FArrayBox& a_state,
                      /// CELL-centered physical domain
                      const Box& domain,
                      /// mesh spacing
                      Real a_dx) const;
  ///
  /** Apply boundary conditions on this face.
   */
  virtual void
  applyEitherBCs(/// NODE-centered data to be modified by boundary condition
                 FArrayBox& a_state,
                 /// CELL-centered physical domain
                 const ProblemDomain& domain,
                 /// mesh spacing
                 Real a_dx,
                 /// flag for homogeneous boundary condition
                 bool a_homogeneous) const;

  ///
  /** Apply boundary conditions on this face.
   */
  virtual void
  applyEitherBCs(/// NODE-centered data to be modified by boundary condition
                 FArrayBox& a_state,
                 /// CELL-centered physical domain
                 const Box& domain,
                 /// mesh spacing
                 Real a_dx,
                 /// flag for homogeneous boundary condition
                 bool a_homogeneous) const;

  ///
  /** Apply boundary condition A*dphi/dx + B*phi = C, with
      the coefficients having been set in fillBCValues().
   */
  virtual void
  applyBCs(/// NODEs on this boundary face
           const Box& a_bcbox,
           /// NODE-centered data to be modified by boundary condition
           FArrayBox& a_state,
           /// coefficients of dphi/dx
           const FArrayBox& a_neumfac,
           /// coefficients of phi
           const FArrayBox& a_dircfac,
           /// constant coefficients
           const FArrayBox& a_inhmval,
           /// mesh spacing
           Real a_dx) const;

  /*@}*/

  /**
     \name Functions to fill in boundary-condition coefficients
  */
  /*@{*/

  ///
  /** Set A and B and C in the boundary condition
      A*dphi/dx + B*phi = C.<br>
      The FArrayBoxes are all based on the NODEs of this face.

      Need this function in derived class.
   */
  virtual void
  fillBCValues(/// coefficients of dphi/dx
               FArrayBox& a_neumfac,
               /// coefficients of phi
               FArrayBox& a_dircfac,
               /// constant coefficients
               FArrayBox& a_inhmval,
               /// mesh spacing
               Real a_dx,
               /// CELL-centered physical domain
               const ProblemDomain& domain) const = 0;

  ///
  /** Set A and B and C in the boundary condition
      A*dphi/dx + B*phi = C.<br>
      The FArrayBoxes are all based on the NODEs of this face.

      Need this function in derived class.
   */
  virtual void
  fillBCValues(/// coefficients of dphi/dx
               FArrayBox& a_neumfac,
               /// coefficients of phi
               FArrayBox& a_dircfac,
               /// constant coefficients
               FArrayBox& a_inhmval,
               /// mesh spacing
               Real a_dx,
               /// CELL-centered physical domain
               const Box& domain) const = 0;

  /*@}*/

  /** low or high face
   */
  Side::LoHiSide m_side;

  /** dimension of face
   */
  int m_direction;

  /** components of data
   */
  Interval m_components;

private:
  FaceNodeBC(const FaceNodeBC&)
    :m_components(-1,-1)
  {
  }

  virtual void operator=(const FaceNodeBC&)
  {
  }
};


/// Class to enforce boundary conditions
/** This class holds a domain's worth of FaceNodeBCs,
    one for each face of the physical domain.
    The user adds a FaceNodeBC-derived class
    for each face to enforce boundary conditions.

    If the solution is phi and the face normal direction is x,
    the boundary conditions at a face are expressed as<br>
    A*dphi/dx + B*phi = C.

    The class derived from FaceNodeBC provides a function to
    to fill A and B and C.  DomainNodeBC calls this function
    and fills the boundary nodes appropriately.
*/
class DomainNodeBC
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Empty constructor.  Need to assign boundary conditions for each face
      with calls to setFaceNodeBC().
   */
  DomainNodeBC();

  ///
  /** Destructor.
   */
  ~DomainNodeBC();

  ///
  /** Assignment, copies input.
   */
  DomainNodeBC& operator=(const DomainNodeBC&);

  ///
  /** Copy constructor.
   */
  DomainNodeBC(const DomainNodeBC&);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /**
     Retrieve FaceNodeBC for a particular face.
   */
  const FaceNodeBC&
  operator() (int direction, Side::LoHiSide side) const;

  /*@}*/

  /**
     \name Setting functions
  */
  /*@{*/

  ///
  /** Set boundary condition object at a face.
   */
  void
  setFaceNodeBC(const FaceNodeBC& a_bc);

  /*@}*/

  /**
     \name Functions to apply boundary conditions
  */
  /*@{*/

  ///
  /** Apply homogeneous boundary conditions on all faces of physical domain.
  */
  void
  applyHomogeneousBCs(/// array on which boundary conditions will be applied
                      NodeFArrayBox& a_state,
                      /// CELL-centered physical domain
                      const ProblemDomain& a_domain,
                      /// mesh spacing
                      Real a_dx) const;

  ///
  /** Apply homogeneous boundary conditions on all faces of physical domain.
  */
  void
  applyHomogeneousBCs(/// array on which boundary conditions will be applied
                      NodeFArrayBox& a_state,
                      /// CELL-centered physical domain
                      const Box& a_domain,
                      /// mesh spacing
                      Real a_dx) const;

  ///
  /** Apply inhomogeneous boundary conditions on all faces of physical domain.
  */
  void
  applyInhomogeneousBCs(/// array on which boundary conditions will be applied
                        NodeFArrayBox& a_state,
                        /// CELL-centered physical domain
                        const ProblemDomain& a_domain,
                        /// mesh spacing
                        Real a_dx) const;

  ///
  /** Apply inhomogeneous boundary conditions on all faces of physical domain.
  */
  void
  applyInhomogeneousBCs(/// array on which boundary conditions will be applied
                        NodeFArrayBox& a_state,
                        /// CELL-centered physical domain
                        const Box& a_domain,
                        /// mesh spacing
                        Real a_dx) const;

  /*@}*/

protected:

  /** Return <tt>true</tt> if the FaceNodeBC for the specified face
      is defined.
   */
  bool isBCDefined(int a_dir, const Side::LoHiSide a_side) const;

  /** Reset the FaceNodeBC for the specified face to NULL.
   */
  void resetFaceNodeBC(const int a_dir, const Side::LoHiSide a_side);

  Tuple<FaceNodeBC*, SpaceDim> m_loBC;
  Tuple<FaceNodeBC*, SpaceDim> m_hiBC;

private:
};

#include "NamespaceFooter.H"

#endif
