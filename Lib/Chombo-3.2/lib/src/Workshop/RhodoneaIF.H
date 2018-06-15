#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RHODONEAIF_H_
#define _RHODONEAIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a rhodonea curve.
 */
class RhodoneaIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying rhodonea inner radius, outer radius,
      number of lobes (frequency), center, and whether the domain
      is on the inside.
   */
  RhodoneaIF(const Real&     a_innerRadius,
             const Real&     a_outerRadius,
             const int&      a_frequency,
             const RealVect& a_center,
             const bool&     a_inside);

  /// Copy constructor
  RhodoneaIF(const RhodoneaIF& a_inputIF);

  /// Destructor
  virtual ~RhodoneaIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real&     a_innerRadius,
                         Real&     a_outerRadius,
                         int&      a_frequency,
                         RealVect& a_center,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real&     a_innerRadius,
                         const Real&     a_outerRadius,
                         const int&      a_frequency,
                         const RealVect& a_center,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  Real m_innerRadius;   // inner radius
  Real m_outerRadius;   // outer radius
  int  m_frequency;     // frequency (i.e. number of lobes)
  RealVect m_center;    // center
  bool     m_inside;    // inside flag

private:
  RhodoneaIF()
  {
    MayDay::Abort("RhodoneaIF uses strong construction");
  }

  void operator=(const RhodoneaIF& a_inputIF)
  {
    MayDay::Abort("RhodoneaIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
