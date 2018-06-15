#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POLARIF_H_
#define _POLARIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a cylinder in polar coordinates:
    r = a_primaryRadius + a_perturbation*cos(a_frequency*theta)
 */
class PolarIF: public BaseIF
{
public:
  ///
  PolarIF(const Real& a_primaryRadius,
          const Real& a_perturbation,
          const int & a_frequency,
          const bool& a_inside);

  /// Copy constructor
  PolarIF(const PolarIF& a_inputIF);

  /// Destructor
  virtual ~PolarIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real& a_primaryRadius,
                         Real& a_perturbation,
                         int & a_frequency,
                         bool& a_inside)const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real& a_primaryRadius,
                         const Real& a_perturbation,
                         const int & a_frequency,
                         const bool& a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  Real m_primaryRadius;
  Real m_perturbation;
  int  m_frequency;
  bool m_inside;

private:
  PolarIF()
  {
    MayDay::Abort("PolarIF uses strong construction");
  }

  void operator=(const PolarIF& a_inputIF)
  {
    MayDay::Abort("PolarIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
