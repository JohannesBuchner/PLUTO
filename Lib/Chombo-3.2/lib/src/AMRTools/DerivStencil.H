#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DERIVSTENCIL_H_
#define _DERIVSTENCIL_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "IntVect.H"
#include "Vector.H"
#include "NamespaceHeader.H"

///class to encapsulate the operations to create derivs on irreg stencils
/**
   DerivStencil is meant to be used to encapsulate the information
   necessary to take finite difference derivatives at a point in
   space.  Every point in the stencil has a weight.  You add them
   (the boxarrayindex and the weight) at the same time and you can manipulate
   the weights enmasse by real number operations.  Stencils may not
   interact with each other with the same sort of arithmetic because
   that would bring up issues as to what to do when there is
   incomplete intersection between the stencils.

 */
class DerivStencil
{
public:
  ///
  void define();

  ///
  /** return true if  any define function been called.
   */
  bool isDefined() const
  {
    return isdefined;
  }

  ///
  /**
     default constructor;  creates empty vectors
  */
  DerivStencil();

  ///
  /**
     copy constructor;  sets *this = a_dsin
  */
  DerivStencil(const DerivStencil& a_dsin);

  ///
  ~DerivStencil();

  ///
  /**
     make derivstencil empty
  */
  void clear();

  ///
  /**
     return length of vectors
  */
  int size() const;

  ///
  /**
     get iv at ivec
  */
  const IntVect& getIndex(int a_ivec) const;

  ///
  /**
     get weight at ivec
  */
  const Real& getWeight(int a_ivec) const;

  ///
  /**
     add another set if the IntVect is not in the
     stencil already.  \\
     **Add the weight to the existing weight otherwise**
     */
  void accumulate(const IntVect& a_iv,  Real a_weight);

  ///
  /**
     assignment operator
  */
  const DerivStencil& operator=(const DerivStencil& a_dsin);

  ///
  /**
     Multiply each weight by a_facin
     does nothing if vectors are of zero length
  */
  const DerivStencil& operator*=(Real a_facin);

  ///
  /**
     divide each weight by a_denom
     does nothing if vectors are of zero length.
     Be advised ---
     This function does no checking to see if a_denom
     is close to zero.
  */
  const DerivStencil& operator/=(Real a_denom);

  ///
  /**
     add a_facin to each weight
     does nothing if vectors are of zero length
  */
  const DerivStencil& operator+=(Real a_facin);

  ///
  /**
     subtract a_facin from each weight
     does nothing if vectors are of zero length
  */
  const DerivStencil& operator-=(Real a_facin);

protected:
  ///
  Vector<IntVect> m_vectIV;
  ///
  Vector<Real>          m_vectWgt;
  ///
  bool isdefined;
};

#include "NamespaceFooter.H"
#endif
