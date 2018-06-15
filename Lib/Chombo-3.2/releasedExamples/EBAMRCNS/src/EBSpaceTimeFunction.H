#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBSPACETIMEFUNCTION_H_
#define _EBSPACETIMEFUNCTION_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"

#include "NamespaceHeader.H"

/**
 * \class EBSpaceTimeFunction
   This virtual base class represents a function F(x, t) that can be used to represent 
   initial conditions, sources, or anything else that depends upon space and time.
 */
class EBSpaceTimeFunction
{
  public:

  /** Default constructor. You must set the problem domain for this EBSpaceTimeFunction 
   *  by calling the define method before using this object.
   */ 
  EBSpaceTimeFunction();

  /** Constructs a space-time function on a given problem domain.
   *  \param a_domain The problem domain on which this function is defined. 
   */ 
  explicit EBSpaceTimeFunction(const ProblemDomain& a_domain);

  /** Destructor. */
  virtual ~EBSpaceTimeFunction();

  /** Set the problem domain for the space-time function.
   *  \param a_domain The problem domain on which this function is to be defined.
   */
  void setDomain(const ProblemDomain& a_domain);

  /** Evaluate this function on the given EBISLayout, placing the results into the given 
   *  LevelData object.
   *  \param a_data The LevelData object to which the evaluated function values will be 
   *                written.
   *  \param a_indexSpace The EBISLayout upon which the function is evaluated.
   *  \param a_dx The cell spacing used to calculate the spatial position x in evaluating the 
   *              function.
   *  \param a_t The time at which the function will be evaluated.
   */
  void evaluate(LevelData<EBCellFAB>& a_data,
                const EBISLayout& a_indexSpace,
                RealVect a_dx,
                Real a_t) const;

  /** Evaluate this function at time 0 on the given EBISLayout, placing the results into the 
   * given LevelData object.
   *  \param a_data The LevelData object to which the evaluated function values will be 
   *                written.
   *  \param a_indexSpace The EBISLayout upon which the function is evaluated.
   *  \param a_dx The cell spacing used to calculate the spatial position x in evaluating the 
   *              function.
   */
  void evaluate(LevelData<EBCellFAB>& a_data,
                const EBISLayout& a_indexSpace,
                RealVect a_dx) const;

  /** Evaluate this function on the given EBISBox, placing the results into the given 
   *  EBCellFAB.
   *  \param a_FAB The EBCellFAB object to which the evaluated function values will be 
  *                written.
   *  \param a_indexSpaceBox The EBISBox upon which the function is evaluated.
   *  \param a_dx The cell spacing used to calculate the spatial position x in evaluating the 
   *              function.
   *  \param a_t The time at which the function will be evaluated.
   */
  void evaluate(EBCellFAB& a_FAB,
                const EBISBox& a_indexSpaceBox,
                RealVect a_dx,
                Real a_t) const;

  /** Evaluate this function at time 0 on the given EBISBox, placing the results into the given 
   *  EBCellFAB.
   *  \param a_FAB The EBCellFAB object to which the evaluated function values will be 
  *                written.
   *  \param a_indexSpaceBox The EBISBox upon which the function is evaluated.
   *  \param a_dx The cell spacing used to calculate the spatial position x in evaluating the 
   *              function.
   */
  void evaluate(EBCellFAB& a_FAB,
                const EBISBox& a_indexSpaceBox,
                RealVect a_dx) const;

  protected:

  /** Override this method to evaluate the function at all single-valued cells on the 
   *  given EBCellFAB. 
   *  \param a_FAB A FAB for all single-valued cells in a Box.
   *  \param a_box The intersection of the Box on which \a a_FAB is defined with the 
   *               problem domain of the function.
   *  \param a_dx The cell spacing used to calculate the spatial position x in evaluating the 
   *              function.
   *  \param a_t The time at which the function will be evaluated.
   */
  virtual void evaluateOnSingleValuedCells(BaseFab<Real>& a_FAB,
                                           const Box& a_box,
                                           RealVect a_dx,
                                           Real a_t) const = 0;

  /** Override this method to evaluate the function at all multi-valued cells on the 
   *  given EBCellFAB.
   *  \param a_FAB The EBCellFAB to which function values will be written.
   *  \param a_vofIterator A VoFIterator that can be used to iterate over all multi-valued 
   *                       cells for \a a_FAB.
   *  \param a_dx The cell spacing used to calculate the spatial position x in evaluating the 
   *              function.
   *  \param a_t The time at which the function will be evaluated.
   */
  virtual void evaluateOnMultiValuedCells(EBCellFAB& a_FAB,
                                          VoFIterator& a_vofIterator,
                                          RealVect a_dx,
                                          Real a_t) const = 0;

  private:

  // Problem domain.
  ProblemDomain m_domain;

  // Disallowed operations.
  EBSpaceTimeFunction(const EBSpaceTimeFunction&);
  EBSpaceTimeFunction& operator=(const EBSpaceTimeFunction&);
};

// Here's a class whose identity tells the EBAMRCNS to just use IBCs for 
// initial conditions.
class EBSpaceTimeFunctionIBCAdaptor: public EBSpaceTimeFunction
{
  void evaluateOnSingleValuedCells(BaseFab<Real>& a_FAB,
                                   const Box& a_box,
                                   RealVect a_dx,
                                   Real a_t) const {}

  void evaluateOnMultiValuedCells(EBCellFAB& a_FAB,
                                  VoFIterator& a_vofIterator,
                                  RealVect a_dx,
                                  Real a_t) const {}

};

#include "NamespaceFooter.H"
#endif
