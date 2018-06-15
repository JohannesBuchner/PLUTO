#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TILTEDCYLINDERIF_H_
#define _TILTEDCYLINDERIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a tilted (infinite) cylinder
 */
class TiltedCylinderIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying cylinder radius (a_radius), direction of the
      axis (a_direction), a point on the axis (a_point), and whether the
      domain is on the inside (a_inside).
   */
  TiltedCylinderIF(const Real&     a_radius,
                   const RealVect& a_direction,
                   const RealVect& a_point,
                   const bool&     a_inside);

  /// Copy constructor
  TiltedCylinderIF(const TiltedCylinderIF& a_inputIF);

  /// Destructor
  virtual ~TiltedCylinderIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real&     a_radius,
                         RealVect& a_direction,
                         RealVect& a_point,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real&     a_radius,
                         const RealVect& a_direction,
                         const RealVect& a_point,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual bool fastIntersection(const RealVect& a_low, const RealVect& a_high) const
  {return !m_tiltedAxis;}


  virtual GeometryService::InOut InsideOutside(const RealVect& a_low, const RealVect& a_high) const ;

protected:
  Real     m_radius;      // radius
  RealVect m_direction;   // axis direction
  bool     m_tiltedAxis;  // false if straight cylinder
  int      m_coordDir;    // coordinate direction of axis for a straight cylinder
  RealVect m_point;       // point on the axis
  bool     m_inside;      // inside flag

  Real     m_radius2;     // precomputed radius squared
  Real     m_length2;     // precomputed length squared of direction vector

private:
  TiltedCylinderIF()
  {
    MayDay::Abort("TiltedCylinderIF uses strong construction");
  }

  void operator=(const TiltedCylinderIF& a_inputIF)
  {
    MayDay::Abort("TiltedCylinderIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
