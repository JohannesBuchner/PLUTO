#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TORUSIF_H_
#define _TORUSIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a torus whose axis is parallel to the
    z-axis
 */
class TorusIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying torus major radius (a_majorRadius), minor radius
      (a_minorRadius), center (a_center), and whether the domain is on the
      inside (a_inside).
   */
  TorusIF(const Real&     a_majorRadius,
          const Real&     a_minorRadius,
          const RealVect& a_center,
          const bool&     a_inside);

  /// Copy constructor
  TorusIF(const TorusIF& a_inputIF);

  /// Destructor
  virtual ~TorusIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real&     a_majorRadius,
                         Real&     a_minorRadius,
                         RealVect& a_center,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real&     a_majorRadius,
                         const Real&     a_minorRadius,
                         const RealVect& a_center,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  Real     m_majorRadius;   // major radius
  Real     m_minorRadius;   // minor radius

  RealVect m_center;        // center

  bool     m_inside;        // inside flag

  Real     m_minorRadius2;  // precomputed minor radius squared

private:
  TorusIF()
  {
    MayDay::Abort("TorusIF uses strong construction");
  }

  void operator=(const TorusIF& a_inputIF)
  {
    MayDay::Abort("TorusIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
