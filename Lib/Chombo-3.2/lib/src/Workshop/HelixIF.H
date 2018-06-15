#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HELIXIF_H_
#define _HELIXIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function takes one implicit function and move it along a
    helix in 3D.  In 2D the implicit function is simply returned.

    In 3D, if a_vertical is false, the given function is restricted to the
    x-y plane (z = 0), and this restricted function is rotated around the
    z-axis at the specified rate as z changes.

    In 3D, if a_vertical is true, the given function is restricted to the
    x-y plane (z = 0), but this x-y cross section is oriented vertically and
    swept in a helix at the specified rate.

    In 2D, the given function is simply returned.
 */
class HelixIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying one implicit function to be rotated, the
      rotation rate, and whether the domain is on the inside (a_inside),
      i.e. where the function is negative.
   */
  HelixIF(const BaseIF& a_impFunc,
          const Real&   a_rate,
          const bool&   a_inside,
          const bool&   a_vertical = true);

  /// Copy constructor
  HelixIF(const HelixIF& a_inputIF);

  /// Destructor
  virtual ~HelixIF();

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  ///
  /**
     Pass this call onto the IFs contained in this IF class.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
    m_impFunc->boxLayoutChanged(a_newBoxLayout,a_dx);
  }

protected:
  BaseIF* m_impFunc; // implicit function to rotate

  Real m_rate;       // the rate of rotate as z changes

  bool m_inside;     // inside flag

  bool m_vertical;   // vertical cross section

  Real m_minz;       // for vertical, the bounds on the cross section
  Real m_maxz;

private:
  HelixIF()
  {
    MayDay::Abort("HelixIF uses strong construction");
  }

  void operator=(const HelixIF& a_inputIF)
  {
    MayDay::Abort("HelixIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
