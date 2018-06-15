#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEIF_H_
#define _BASEIF_H_

#include "RealVect.H"
#include "ProblemDomain.H"
#include "IndexTM.H"

#include "Notation.H"
#include "GeometryService.H"

#include "NamespaceHeader.H"

///
/**
   This is the base class for an implicit function specification of geometry.
   All that is needed is a constructor/destructor, a method to give the value
   of the function at any point in space (in 2D or 3D), and a factory method
   (these will probably all be very similar).
 */
class BaseIF
{

public:

  /// Default constructor
  BaseIF()
  {
  }

  /// Default destructor
  virtual ~BaseIF()
  {
  }

  ///
  /**
   Return the value of the function at a_point.  When delineating a domain,
   the level set value=0 represents the boundary and value<0 is inside the
   fluid.
  */
  virtual Real value(const RealVect& a_point) const = 0;

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    Real retval= LARGEREALVAL;
    if (a_partialDerivative == IndexTM<int,GLOBALDIM>::Zero)
      {
        retval = value(a_point);
      }
    else
      {
        MayDay::Abort("Not implemented");
      }
    return retval;
  }

  virtual bool fastIntersection(const Box&           a_region,
                                const ProblemDomain& a_domain,
                                const RealVect&      a_origin,
                                const Real&          a_dx) const
  {
    RealVect low, high;
    corners(a_region, a_origin, a_dx, low, high);
    return fastIntersection(low, high);
  }

  virtual bool fastIntersection(const RealVect& a_low,
                                const RealVect& a_high) const
  {return false;}

  virtual GeometryService::InOut InsideOutside(const Box&           a_region,
                                               const ProblemDomain& a_domain,
                                               const RealVect&      a_origin,
                                               const Real&          a_dx) const
  {
    RealVect low, high;
    corners(a_region, a_origin, a_dx, low, high);
    return InsideOutside(low, high);
  }

  virtual GeometryService::InOut InsideOutside(const RealVect& a_low,
                                               const RealVect& a_high) const

  {
    MayDay::Abort("This class has not implemented a fastIntersection operation");
    return GeometryService::Irregular;
  }

  ///
  /**
     Return the value of the function at a_point (of type IndexTM).
  */
  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    MayDay::Abort("'value' called with 'IndexTM' not implemented in BaseIF");
    Real retval = LARGEREALVAL;
    return retval;
  };

  ///
  /**
     Return a newly allocated derived class.  The responsibility
     for deleting the memory is left to the calling function.
  */
  virtual BaseIF* newImplicitFunction() const = 0;

  virtual void print(ostream& out) const
  {
    MayDay::Abort("Print function not implemented");
  };

  static void corners(const Box& a_region, const RealVect& a_origin, const Real& a_dx,
                      RealVect& a_lo, RealVect& a_hi)
  {
    a_lo = a_origin + RealVect(a_region.smallEnd())*a_dx;
    a_hi = a_origin + RealVect(a_region.bigEnd()+IntVect::Unit)*a_dx;
  }

  ///
  /**
     Inform the implicit function, IF, of a BoxLayout change in the
     application that is using it.  If the implicit function uses
     distributed data, this gives it a chance to react to the spatial
     layout of the application that is using it.

     This is a empty implementation for IFs that don't need this
     functionality.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
  }
};

#include "NamespaceFooter.H"

#endif
