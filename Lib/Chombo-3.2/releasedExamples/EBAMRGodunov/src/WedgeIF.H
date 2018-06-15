#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _WEDGEIF_H_
#define _WEDGEIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a square (infinite) cylinder
    around the x axis
 */
class WedgeIF: public BaseIF
{
public:
  ///
  /**
     !                        depVar
     !            \            ^
     !           / \           |
     !          /   \          |--->indVar
     !         /     \
     !        /       \
     !-------/         \------------
     center is the point in indvar
   */
  WedgeIF(const RealVect& a_center,
          const Real&     a_slope,
          const int&      a_depVar,
          const int&      a_indVar);

  /// Copy constructor
  WedgeIF(const WedgeIF& a_inputIF);

  /// Destructor
  virtual ~WedgeIF();


  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  RealVect m_center;
  Real     m_slope;
  int      m_depVar;
  int      m_indVar;

private:
  WedgeIF()
  {
    MayDay::Error("WedgeIF uses strong construction");
  }

  void operator=(const WedgeIF& a_inputIF)
  {
    MayDay::Error("WedgeIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
