#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _SQUARECYLINDERIF_H_
#define _SQUARECYLINDERIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a square (infinite) cylinder
    around the x axis
 */
class SquareCylinderIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying cylinder radius (a_radius), direction of the
      axis (a_direction), a point on the axis (a_point), and whether the
      domain is on the inside (a_inside).
   */
  SquareCylinderIF(const Real&     a_radius,
                   const bool&     a_inside);

  /// Copy constructor
  SquareCylinderIF(const SquareCylinderIF& a_inputIF);

  /// Destructor
  virtual ~SquareCylinderIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real&     a_radius,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real&     a_radius,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  Real     m_radius;      // radius
  bool     m_inside;      // inside flag


private:
  SquareCylinderIF()
  {
    MayDay::Error("SquareCylinderIF uses strong construction");
  }

  void operator=(const SquareCylinderIF& a_inputIF)
  {
    MayDay::Error("SquareCylinderIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
