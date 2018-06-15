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

#ifndef _EBCELLFAB_H_
#define _EBCELLFAB_H_

#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "BaseEBCellFAB.H"
#include "FArrayBox.H"
#include "EBISBox.H"
#include "Box.H"
#include "NamespaceHeader.H"
///
/**
   EBCellFAB is a BaseEBCellFAB<Real>
   which includes arithmetic functions.
 */
class EBCellFAB: public BaseEBCellFAB<Real>
{
public:
  static bool s_verbose;
  ///
  EBCellFAB();

  ///
  /**
    Box going into this needs to be cell-centered.
    Data will exist over the surrounding nodes of the box.
   */
  EBCellFAB(const EBISBox& a_ebisBox,
            const Box& a_region, int a_nVar);


  ///
  void define(const EBISBox& a_ebisBox,
              const Box& a_region, int a_nVar);

  ///
  void setInvalidData(const Real& a_val,
                      const int& a_comp);

  ///
  virtual ~EBCellFAB();

  ///
  const FArrayBox& getFArrayBox() const;

  ///
  /**
     Does a dynamic cast  of the BaseFab<Real>
     to a  FArrayBox for arithmetic purposes.
   */
  FArrayBox& getFArrayBox();

  ///
  /**
     Negate the values in all locations
   */
  EBCellFAB& negate(void);

  ///
  /**
     done over intersection of src, dst
     Both fabs need the same ncomp
   */
  EBCellFAB& operator+=(const EBCellFAB& a_ebfab);


  ///
  /**
     done over intersection of src, dst
  */
  EBCellFAB& plus(const EBCellFAB& a_ebfab,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp);
  ///
  /**
     done over a_region
  */
  EBCellFAB& plus(const EBCellFAB& a_ebfab,
                  const Box& a_region,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp);


  ///
  /**
     done over intersection of src, dst
     Both fabs need the same ncomp
   */
  EBCellFAB& operator-=(const EBCellFAB& a_ebfab);

  ///
  /**
     done over intersection of src, dst
  */
  EBCellFAB& minus(const EBCellFAB& a_ebfab,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);

  ///
  /**
     done over intersection of src, dst
     Both fabs need the same ncomp
  */
  EBCellFAB& operator*=(const EBCellFAB& a_ebfab);

  ///
  /**
     done over intersection of src, dst
  */
  EBCellFAB& mult(const EBCellFAB& a_ebfab,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp);



  ///
  /**
     done over intersection of src, dst
     Both fabs need the same ncomp
  */
  EBCellFAB& operator/=(const EBCellFAB& a_ebfab);


  ///
  /**
     done over intersection of src, dst
  */
  EBCellFAB& divide(const EBCellFAB& a_ebfab,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp);

  ///
  /**
  */
  EBCellFAB& operator+=(const Real& a_scalar);

  ///
  /**
  */
  EBCellFAB& operator-=(const Real& a_scalar);

  ///
  /**
  */
  EBCellFAB& operator*=(const Real& a_scalar);

  ///
  /**
  */
  EBCellFAB& mult(Real a_scalar);

  ///
  /**
  */
  EBCellFAB& operator/=(const Real& a_scalar);

  ///
  /**
     Current FAB += a_src FAB * a_scalar.  Both fabs need the same ebisBox
     and region.
  */
  EBCellFAB& plus(const EBCellFAB& a_src,
                  Real             a_scalar);

  ///
  /**
     Current FAB = a_X FAB * a_A  +  a_Y FAB * a_B.
  */
  EBCellFAB& axby(const EBCellFAB& a_X, const EBCellFAB& a_Y,
                  const Real& a_A, const Real& a_B);

  ///
  /**
     Current FAB = a_X FAB * a_A  +  a_Y FAB * a_B.
  */
  EBCellFAB& axby(const EBCellFAB& a_X, const EBCellFAB& a_Y,
                  const Real& a_A, const Real& a_B,
                  const int& a_destComp,const int& a_xComp,const int& a_yComp);

  //! Returns the maximum value in the FAB for the given component.
  Real max(int a_comp = 0) const;

  //! Returns the minimum value in the FAB for the given component.
  Real min(int a_comp = 0) const;

  /// (Not implemented) Returns the Lp-norm of this EBCellFAB
  /**
     (Not implemented) Returns the Lp-norm of this EBCellFAB using components
     (a_comp : a_comp + a_numcomp - 1).  a_power < 0 -> ERROR.
     a_power = 0  -> infinity norm (max norm).
     a_power = 1  -> L1-norm
     a_power > 1  -> Lp-norm
  */
  virtual Real norm(int a_power,
                    int a_comp,
                    int a_numComp) const;

  /// (Not implemented) Returns the Lp-norm of this EBCellFAB within a region
  /**
     (Not implemented) Returns the Lp-norm of this EBCellFAB using components
     (a_comp : a_comp + a_numcomp - 1) and within the a_subbox.  a_power < 0
     -> ERROR.
     a_power = 0 -> infinity norm (max norm).
     a_power = 1 -> L1-norm
     a_power > 1 -> Lp-norm
  */
  virtual Real norm(const Box& a_subbox,
                    int        a_power,
                    int        a_comp,
                    int        a_numComp) const;

  /// (Not implemented) Returns a sum of powers of a subset of this EBCellFAB
  /**
     (Not implemented) Returns a sum of powers of a subset of this EBCellFAB,
     specifically components a_comp to a_comp+a_numcomp-1 within a_subbox.
     a_power >= 2 only.

  */
  virtual Real sumPow(const Box& a_subbox,
                      int        a_power,
                      int        a_comp,
                      int        a_numComp) const;

  /// (Not implemented) Return the dot product of this EBCellFAB with another
  /**
     (Not implemented) Return the dot product of this EBCellFAB and "ebfab2"
     over their overlap region and all components.
  */
  Real dotProduct(const EBCellFAB& ebfab2) const;

  //needed for leveldata to compile
  //but very invalid
  EBCellFAB(const Box& a_region, int a_nComp)
  {
    MayDay::Error("ebcelfab 1 invalid operator");
  }

  /// good for making temps.
  void clone(const EBCellFAB& a_arg);

private:
  //disallowed for all the usual reasons
  EBCellFAB(const EBCellFAB& ebcin)
  {
    MayDay::Error("ebcellfab 2 invalid operator");
  }
  void operator=(const EBCellFAB& fabin)
  {
    MayDay::Error("ebcellfab 3 invalid operator");
  }

};

void writeVectorLevelName(const Vector<LevelData<EBCellFAB>*>*, Vector<int>* ref, const char*);

#include "NamespaceFooter.H"
#endif
