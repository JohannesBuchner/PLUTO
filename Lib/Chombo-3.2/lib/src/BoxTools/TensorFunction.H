#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TENSORFUNCTION_H_
#define _TENSORFUNCTION_H_

#include "REAL.H"
#include "RealVect.H"
#include "RealTensor.H"
#include "NamespaceHeader.H"

//! \class TensorFunction
//! This base class represents a tensor function \f$T: \mathbf{R}^D \rightarrow \mathbf{R}^D\otimes\mathbf{R}^D\f$.
class TensorFunction
{
  public:

  //! Base class constructor. Must be called by subclasses.
  //! \param a_homogeneous This flag indicates whether the tensor function
  //!                      is constant in space.
  //! \param a_constant This flag indicates whether the tensor function
  //!                   is constant in time.
  TensorFunction(bool a_homogeneous,
                 bool a_constant);

  //! Destructor.
  virtual ~TensorFunction();

  //! Override this method to evaluate this function at the given
  //! point in space and time.
  //! \param a_x A point in \f$D\f$-dimensional space.
  //! \param a_t The time at which the function is evaluated.
  virtual RealTensor operator()(const RealVect& a_x,
                            Real a_t) const = 0;

  //! This evaluates the function at time 0.
  //! \param a_x A point in \f$D\f$-dimensional space.
  RealTensor operator()(const RealVect& a_x) const
  {
    return operator()(a_x, 0.0);
  }

  //! Returns true if this function is homogeneous, false otherwise.
  bool isHomogeneous() const
  {
    return m_isHomogeneous;
  }

  //! Returns true if this function is constant, false otherwise.
  bool isConstant() const
  {
    return m_isConstant;
  }

  protected:

  bool m_isHomogeneous, m_isConstant;

  private:

  // Disallowed!
  TensorFunction();
  TensorFunction(const TensorFunction&);
  TensorFunction& operator=(const TensorFunction&);
};

#include "NamespaceFooter.H"
#endif
