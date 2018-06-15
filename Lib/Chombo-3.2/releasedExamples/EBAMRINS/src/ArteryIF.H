#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ARTERYIF_H_
#define _ARTERYIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

// A specialized geometry class to create a few examples of something which
// looks like a segment of a artery.
//
// If "type" is 1, the geometry is a surface of revolution (revolved about
// the axis defined by (x,a_center[1],a_center[2])).  Thus, the centerline
// is a line in the x direction and radius varies along that centerline.
//
// If "type" is 2, the geometry has a varying centerline and radius.
//
// In both cases, the artery starts at the plane x = 0 and ends at the plane
// x = 1.

class ArteryIF: public BaseIF
{
public:
  // "a_type" it the artery type (see above)
  // "a_radius" is the initial radius (at the plane x = 0)
  // "a_amplitude" is the amplitude of the perturbation (zero -> cylinder)
  // "a_start" is the of start the artery in x
  // "a_end" is the end of the artery in x
  // "a_center" is the center of the inlet circle (the x coordinate is ignored)
  // "a_inside" being true means the domain is inside the artery, otherwise it
  //            is outside
  ArteryIF(const int&      a_type,
           const Real&     a_radius,
           const Real&     a_amplitude,
           const Real&     a_start,
           const Real&     a_end,
           const RealVect& a_center,
           const bool&     a_inside);

  ArteryIF(const ArteryIF& a_inputIF);

  virtual ~ArteryIF();

  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  // The artery type (see above)
  int  m_type;

  // The initial radius (at the plane x = 0)
  Real m_radius;

  // The amplitude of the perturbation (zero -> cylinder)
  Real m_amplitude;

  // The start of the artery in x
  Real m_start;

  // The end of the artery in x
  Real m_end;

  // The center of the inlet circle (the x coordinate is ignored)
  RealVect m_center;

  // True means the domain is inside the artery, otherwise it is outside
  bool m_inside;

private:
  ArteryIF()
  {
    MayDay::Error("ArteryIF uses strong construction");
  }

  void operator=(const ArteryIF& a_inputIF)
  {
    MayDay::Error("ArteryIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
