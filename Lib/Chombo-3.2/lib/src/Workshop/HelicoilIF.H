#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HELICOILIF_H_
#define _HELICOILIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"
#include "HelixIF.H"

#include "NamespaceHeader.H"

///
/**
    In 3D this is an approximation of a circle swept around a helix which is
    always in the plane normal to the tangent of the helix.  If the slope
    of the helix is less than 1.0, this is achieved by sweeping an ellipse
    oriented vertically along a helix.  If the slope is greater than or equal
    to 1.0, this is achieved by sweeping an ellipse oriented horizontally
    along a helix.  The eccentricity of the ellipse is adjusted based on the
    slope.

    In 2D this is an ellipse offset from the origin.
 */
class HelicoilIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying the helix radius, pitch, radius of the swept
      circle, and whether the domain is on the inside (a_inside),
      i.e. where the function is negative.
   */
  HelicoilIF(const Real& a_helixR,
             const Real& a_helixPitch,
             const Real& a_circleR,
             const bool& a_inside);

  /// Copy constructor
  HelicoilIF(const HelicoilIF& a_inputIF);

  /// Destructor
  virtual ~HelicoilIF();

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
    m_helixIF->boxLayoutChanged(a_newBoxLayout,a_dx);
  }

protected:
  Real m_helixR;     // the radius of the helix
  Real m_helixPitch; // the pitch of the helix
  Real m_circleR;    // the radius of the swept circle

  bool m_inside;     // inside flag

  Real m_ellipseR;
  bool m_vertical;

  BaseIF* m_helixIF;

private:
  HelicoilIF()
  {
    MayDay::Abort("HelicoilIF uses strong construction");
  }

  void operator=(const HelicoilIF& a_inputIF)
  {
    MayDay::Abort("HelicoilIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
