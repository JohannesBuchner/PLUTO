#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPHEREARRAYIF_H_
#define _SPHEREARRAYIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
   Constructor specifying an infinite array of spheres
   starting at first center and going in all directions
   spaced out by spacing
 */
class SphereArrayIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying an infinite array of spheres
      starting at first center and going in all directions
      spaced out by spacing
   */
  SphereArrayIF(const Real&     a_radius,
                const RealVect& a_firstCenter,
                const RealVect& a_spacing);

  /// Copy constructor
  SphereArrayIF(const SphereArrayIF& a_inputIF);

  /// Destructor
  virtual ~SphereArrayIF();

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
                                               const RealVect& a_high) const;

protected:
  Real      m_radius;
  RealVect  m_firstCenter;
  RealVect  m_spacing;

  Real     m_radius2;   // precomputed radius squared

private:
  SphereArrayIF()
  {
    MayDay::Abort("SphereArrayIF uses strong construction");
  }

  void operator=(const SphereArrayIF& a_inputIF)
  {
    MayDay::Abort("SphereArrayIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
