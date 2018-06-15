#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TYLERCHANNELIF_H_
#define _TYLERCHANNELIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a sphere.
 */
class TylerChannelIF: public BaseIF
{
public:
  ///
  /**
      A C-infinity contraction/expansion channel.

      The geometry is centered at [x_M,y_M] = [(a_x1+a_x2)/2,yDomainLength/2]
      and is symmetric about y = y_M.  The size of the channel goes from
      2*a_y1 to 2*a_y2 as x goes from -infinity to infinity.

      The boundary is given by:

          y(x)     = y_M +/- yDist(x)
          yDist(x) = a_y1*(1-frac(x)) + a_y2*frac(x)
          frac(x)  = (tanh(theta(x)) + 1) / 2
          theta(x) = 4 * (x - x_M)/(a_x2 - a_x1)
   */
  TylerChannelIF(const Real& a_x1,
                 const Real& a_x2,
                 const Real& a_y1,
                 const Real& a_y2,
                 const Real& a_yDomainLength);

  /// Copy constructor
  TylerChannelIF(const TylerChannelIF& a_inputIF);

  /// Destructor
  virtual ~TylerChannelIF();

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  Real m_x1;
  Real m_x2;

  Real m_y1;
  Real m_y2;

  Real m_yDomainLength;

private:
  TylerChannelIF()
  {
    MayDay::Abort("TylerChannelIF uses strong construction");
  }

  void operator=(const TylerChannelIF& a_inputIF)
  {
    MayDay::Abort("TylerChannelIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
