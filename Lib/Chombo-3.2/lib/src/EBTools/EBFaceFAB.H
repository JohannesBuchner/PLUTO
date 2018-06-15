#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _EBFACEFAB_H_
#define _EBFACEFAB_H_

#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "BaseEBFaceFAB.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

///
/**
   EBFaceFAB is a BaseEBFaceFAB<Real>
   which includes arithmetic functions.
 */
class EBFaceFAB: public BaseEBFaceFAB<Real>
{
public:
  ///
  EBFaceFAB();

  ///
  /**
    Box going into this needs to be cell-centered.
    Data will exist over the surrounding nodes of the box.
   */
  EBFaceFAB(const EBISBox& a_ebisBox,
            const Box& a_region,
            int a_iDir, int a_nVar);


void define(const EBISBox&  a_ebisBox,
            const Box& a_region,
            int a_iDir, int a_nComp);
  ///
  virtual ~EBFaceFAB();

  ///
  const FArrayBox& getFArrayBox() const;

  ///
  /**
     Does a dynamic cast of the BaseFab<Real>
     to a  FArrayBox for arithmetic purposes.
   */
  FArrayBox& getFArrayBox();

  ///
  /**
     Both fabs need the same ebisBox and region and number of components.
   */
  EBFaceFAB& operator+=(const EBFaceFAB& a_ebfab);


  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  EBFaceFAB& plus(const EBFaceFAB& a_ebfab,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp);


  ///
  /**
     Both fabs need the same ebisBox and region and number of components.
  */
  EBFaceFAB& operator-=(const EBFaceFAB& a_ebfab);

  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  EBFaceFAB& minus(const EBFaceFAB& a_ebfab,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);


  ///
  /**
     Both fabs need the same ebisBox and region and number of components.
  */
  EBFaceFAB& operator*=(const EBFaceFAB& a_ebfab);

  ///
  /**
     Both fabs need the same ebisBox and region.
  */
  EBFaceFAB& mult(const EBFaceFAB& a_ebfab,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp);


  ///
  /**
     Both fabs need the same ebisBox and region and number of components.
  */
  EBFaceFAB& operator/=(const EBFaceFAB& a_ebfab);


  ///
  /**
     Both fabs need the same ebisBox and region.
  */
  EBFaceFAB& divide(const EBFaceFAB& a_ebfab,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp);


  ///
  /**
   */
  EBFaceFAB& operator+=(const Real& a_scalar);

  ///
  /**
   */
  EBFaceFAB& operator-=(const Real& a_scalar);

  ///
  /**
   */
  EBFaceFAB& operator*=(const Real& a_scalar);

  ///
  /**
   */
  EBFaceFAB& operator/=(const Real& a_scalar);

  //! Returns the maximum value for the given component in the FAB.
  Real max(int a_comp = 0) const;

  //! Returns the minimum value for the given component in the FAB.
  Real min(int a_comp = 0) const;

  //needed for leveldata to compile but invalid
  EBFaceFAB(const Box& a_region, int a_nComp)
  {
    MayDay::Error("bogus constructor");
  }
private:
  //verboten for all the usual reasons
  EBFaceFAB(const EBFaceFAB& ebcin)
  {
    MayDay::Error("bogus constructor");
  }
  void operator=(const EBFaceFAB& ebcin)
  {
    MayDay::Error("bogus constructor");
  }

};

#include "NamespaceFooter.H"
#endif
