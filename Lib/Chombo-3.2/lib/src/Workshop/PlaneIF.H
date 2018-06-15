#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PLANEIF_H_
#define _PLANEIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a plane.
 */
class PlaneIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying plane normal (a_normal), a point on the plance
      (a_point), and whether the domain is the same side as the normal
      (a_inside).
   */
  PlaneIF(const RealVect& a_normal,
          const RealVect& a_point,
          const bool&     a_inside);

  /// Copy constructor
  PlaneIF(const PlaneIF& a_inputIF);

  /// Destructor
  virtual ~PlaneIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(RealVect& a_normal,
                         RealVect& a_point,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const RealVect& a_normal,
                         const RealVect& a_point,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual bool fastIntersection(const RealVect& a_low,
                                const RealVect& a_high)const
  {return true;}


  virtual GeometryService::InOut InsideOutside(const RealVect& a_low,
                                               const RealVect& a_high) const ;
protected:
  RealVect m_normal;    // plane normal
  RealVect m_point;     // point on the plane
  bool     m_inside;    // inside flag

private:
  PlaneIF()
  {
    MayDay::Abort("PlaneIF uses strong construction");
  }

  void operator=(const PlaneIF& a_inputIF)
  {
    MayDay::Abort("PlaneIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
