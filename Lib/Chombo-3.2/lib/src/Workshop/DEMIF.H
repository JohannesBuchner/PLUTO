#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DEMIF_H_
#define _DEMIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
   This implicit function specifies a Digital Elevation Model
*/
class DEMIF: public BaseIF
{
public:
  ///
  /**
     Constructor specifying a digital elevation model
  */
  DEMIF(const IntVect&     a_ncell,
        const int&         a_interpType,
        const RealVect&    a_dx,
        const std::string& a_demFile,
        const Real&        a_bottomBuffer,
        const Real&        a_truncElev,
        const Real&        a_highGround,
        const Real&        a_verticalScale);

  /// Copy constructor
  DEMIF(const DEMIF& a_inputIF);

  /// Destructor
  virtual ~DEMIF();

  ///
  /**
     Return the value of the function at a_point.
  */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  /// Return true if DEM is read
  virtual bool readDEM(const bool&        a_justhead,
                       const std::string& a_demFile);

  /// This caches constants based on m_cellsize
  virtual bool cacheConstants();

  ///
  /*
    Return true if DEM is read
  */
  virtual bool fixDEM();


protected:
  int      m_ncols;
  int      m_nrows;
  int      m_NODATA;
  Real     m_xllcorner;
  Real     m_yllcorner;
  Real     m_cellsize;
  Real     m_cellvalue;
  Real     **m_DEM;
  RealVect m_dx;
  IntVect  m_ncell;
  int      m_interpType;
  Real     m_highGround;
  Real     m_minDEM;
  Real     m_maxDEM;
  Real     m_bottomBuffer;
  Real     m_truncElev;
  Real     m_verticalScale;
  bool     m_doCubic;
  Real     m_hx;
  Real     m_hy;
  Real     m_hx2;
  Real     m_hx3;
  Real     m_hy2;
  Real     m_hy3;
private:
  DEMIF()
  {
    MayDay::Abort("DEMIF uses strong construction");
  }

  void operator=(const DEMIF& a_inputIF)
  {
    MayDay::Abort("DEMIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
