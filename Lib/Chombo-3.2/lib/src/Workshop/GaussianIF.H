#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GAUSSIANIF_H_
#define _GAUSSIANIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a gaussian level set.  The gaussian
    level set's dependent variable is in the "upDir" direction and it is a
    function of the other directions.  The gaussian level set is centered
    at
 */
class GaussianIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying the gaussian level set's center (a_origin),
      its height (a_height), the square of its standard deviation (a_sig2),
      its dependent variable direction (a_upDir), and whether the domain
      is the same side as the normal (a_inside).
   */
  GaussianIF(const RealVect& a_origin,
             const Real&     a_height,
             const RealVect& a_sig2,
             const int&      a_upDir,
             const bool&     a_inside);

  /// Copy constructor
  GaussianIF(const GaussianIF& a_inputIF);

  /// Destructor
  virtual ~GaussianIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(RealVect& a_origin,
                         Real&     a_height,
                         RealVect& a_sig2,
                         int&      a_upDir,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const RealVect& a_origin,
                         const Real&     a_height,
                         const RealVect& a_sig2,
                         const int&      a_upDir,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  RealVect m_origin;    // gaussian origin
  Real     m_height;    // height of the gaussian
  RealVect m_sig2;      // standard deviation of gaussian
  int      m_upDir;     // the up direction for the gaussian (the symmetry direction)
  bool     m_inside;    // inside flag

private:
  GaussianIF()
  {
    MayDay::Abort("GaussianIF uses strong construction");
  }

  void operator=(const GaussianIF& a_inputIF)
  {
    MayDay::Abort("GaussianIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
