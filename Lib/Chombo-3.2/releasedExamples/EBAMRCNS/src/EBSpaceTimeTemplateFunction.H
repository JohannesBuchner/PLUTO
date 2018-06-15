#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBSPACETIMETEMPLATEFUNCTION_H_
#define _EBSPACETIMETEMPLATEFUNCTION_H_

#include "EBSpaceTimeFunction.H"

#include "NamespaceHeader.H"

/**
 * \class EBSpaceTimeTemplateFunction
   This class template implements a function F(x, t) by using a functor to evaluate data 
   at points in space and time. This class is meant to be used as a quick-n-dirty way to 
   implement initial conditions--since it doesn't use Chombo fortran, its performance is 
   not guaranteed to impress!
 */
template <typename Functor>
class EBSpaceTimeTemplateFunction: public EBSpaceTimeFunction
{
  public:

  /** Construct an EBSpaceTimeTemplateFunction from a functor. The domain must be set by
   *  calling the define() method before the function is used.
   *  \param a_functor A functor defining the following methods:
   *     - interval() const, which returns an Interval describing the function's range.
   *     - void operator()(F, X, t) const, where F is an array of length interval().size(), 
   *                                       X is a RealVect storing the coordinates of a spatial 
   *                                       point, and t is the time at which the function is 
   *                                       evaluated.
   */ 
  EBSpaceTimeTemplateFunction(Functor* a_functor):
    EBSpaceTimeFunction(), m_functor(a_functor) {}

  /** Construct an EBSpaceTimeTemplateFunction from a functor and a domain.
   *  \param a_functor A functor defining the following methods:
   *     - interval() const, which returns an Interval describing the function's range.
   *     - void operator()(F, X, t) const, where F is an array of length interval().size(), 
   *                                       X is a RealVect storing the coordinates of a spatial 
   *                                       point, and t is the time at which the function is 
   *                                       evaluated.
   *  \param a_domain The problem domain on which this function is defined. 
   */ 
  EBSpaceTimeTemplateFunction(Functor* a_functor, 
                              const ProblemDomain& a_domain):
    EBSpaceTimeFunction(a_domain), m_functor(a_functor) {}

  /** Destructor. */
  ~EBSpaceTimeTemplateFunction() { delete m_functor; }

  protected:

  // This method uses the functor to evaluate the function on single-valued cells for the 
  // given FAB.
  void evaluateOnSingleValuedCells(BaseFab<Real>& a_FAB,
                                   const Box& a_box,
                                   RealVect a_dx,
                                   Real a_t) const
  {
    // Get the interval from the functor.
    Interval interval = m_functor->interval();

    // We'll store each function value here.
    Vector<Real> data(interval.size());

    for (BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      // Get the spatial coordinate for this VoF.
      IntVect index = bit();
      RealVect xi(index);
      xi *= a_dx;

      // Call the functor at this point and time.
      (*m_functor)(data, xi, a_t);

      // Assign the data to the FAB.
      for (int i = interval.begin(); i <= interval.end(); ++i)
        a_FAB(index, i) = data[i]; 
    }
  }

  // This method uses the functor to evaluate the function on multi-valued cells for the 
  // given EBCellFAB.
  void evaluateOnMultiValuedCells(EBCellFAB& a_FAB,
                                  VoFIterator& a_vofIterator,
                                  RealVect a_dx,
                                  Real a_t) const
  {
    // Get the interval from the functor.
    Interval interval = m_functor->interval();

    // We'll store each function value here.
    Vector<Real> data(interval.size());

    // Loop over the VoFs.
    for (; a_vofIterator.ok(); ++a_vofIterator)
    {
      // Get the spatial coordinate for this VoF.
      VolIndex vi = a_vofIterator();
      IntVect gi = vi.gridIndex();
      RealVect xi(gi);
      xi *= a_dx;

      // Call the functor at this point and time.
      (*m_functor)(data, xi, a_t);

      // Assign the data to the FAB.
      for (int i = interval.begin(); i <= interval.end(); ++i)
        a_FAB(vi, i, 1) = data[i]; // vi is known to be multi-valued!
    }
  }

  private:

  //! The C++ functor that implements the function.
  Functor* m_functor;

  // Disallowed operations.
  EBSpaceTimeTemplateFunction();
  EBSpaceTimeTemplateFunction(const EBSpaceTimeTemplateFunction&);
  EBSpaceTimeTemplateFunction& operator=(const EBSpaceTimeTemplateFunction&);
};

#include "NamespaceFooter.H"
#endif
