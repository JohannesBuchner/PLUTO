#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPHEREIF_H_
#define _SPHEREIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a sphere.
 */
class SphereIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying sphere radius (a_radius), center (a_center), and
      whether the domain is on the inside (a_inside).
   */
  SphereIF(const Real&     a_radius,
           const RealVect& a_center,
           const bool&     a_inside);

  /// Copy constructor
  SphereIF(const SphereIF& a_inputIF);

  /// Destructor
  virtual ~SphereIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real&     a_radius,
                         RealVect& a_center,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real&     a_radius,
                         const RealVect& a_center,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual bool fastIntersection(const RealVect& a_low,
                                const RealVect& a_high) const
  {return true;}


  virtual GeometryService::InOut InsideOutside(const RealVect& a_low,
                                               const RealVect& a_high) const ;

protected:
  Real     m_radius;    // radius
  RealVect m_center;    // center
  bool     m_inside;    // inside flag

  Real     m_radius2;   // precomputed radius squared

private:
  SphereIF()
  {
    MayDay::Abort("SphereIF uses strong construction");
  }

  void operator=(const SphereIF& a_inputIF)
  {
    MayDay::Abort("SphereIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
