#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COORDINATESYSTEM_H_
#define _COORDINATESYSTEM_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif


#include <map>
using std::map;
#include "Vector.H"
#include "REAL.H"
#include "IndexTM.H"

#include "Notation.H"
#include "IFSlicer.H"

using std::endl;
using std::cerr;
#include "NamespaceHeader.H"

template <int dim> class CoordinateSystem
{
public:
  // data types for vector of int and Real in dimension "dim"
  typedef IndexTM<int,dim>  IvDim;
  typedef IndexTM<Real,dim> RvDim;

  // member data
  RvDim m_origin;
  RvDim m_dx;

  // empty constructor
  CoordinateSystem();

  // copy constructor
  CoordinateSystem(const CoordinateSystem<dim>& a_info);

  // makes an CoordinateSystem from a origin and dx
  CoordinateSystem(const RvDim  & a_origin,
                   const RvDim  & a_dx);

  // makes a reduced CoordinateSystem
  CoordinateSystem(const CoordinateSystem<dim+1> & a_hISystem,
                   const int                     & a_idir);

  // destructor
  ~CoordinateSystem();

  // convert "a_point" in coordinate system "a_system" to a point in the
  // current coordinate system.
  RvDim convert(const RvDim                 & a_point,
                const CoordinateSystem<dim> & a_system) const;

  // convert "a_coord" in direction "a_dir" in coordinate system "a_system"
  // to a value in "a_dir" in the current coordinate system.
  Real convertDir(const Real                  & a_coord,
                  const CoordinateSystem<dim> & a_system,
                  const int                   & a_dir) const;

  // output functions
  void print(ostream& out) const;

  // operators
  void operator=(const CoordinateSystem<dim> & a_coordinateSystem);
};

#include "NamespaceFooter.H"

#include "CoordinateSystemImplem.H"

#endif
