#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIONFUNCTIONS_H_
#define _ADVECTIONFUNCTIONS_H_

#include "IntVect.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

/// Velocity function interface.
/**
*/
///
typedef Real (*AdvectionVelocityFunction)(/// location of the function evaluation
                                          const RealVect& a_point,
                                          /// velocity component
                                          const int&      a_velComp);


/// Constant flow example
/**
   Return 1.
 */
Real constantAdvection(/// location of the function evaluation
                       const RealVect& a_point,
                       /// velocity component
                       const int&      a_velComp);


/// Rotating flow example
/**
   At point (x, y, z),
   set d = (x - 1/2)^2 + (y - 1/2)^2 + (z - 1/2)^2.
   \n
   If d < 1/2, then return velocity ((y - 1/2)/d, -(x - 1/2)/d, 0).
   \n
   Otherwise, return velocity zero.
 */
Real rotatingAdvection(/// location of the function evaluation
                       const RealVect& a_point,
                       /// velocity component
                       const int&      a_velComp);

#include "NamespaceFooter.H"

#endif
