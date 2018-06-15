#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STLIF_H_
#define _STLIF_H_

#include <iostream>
#include <fstream>
using std::istream;
using std::ifstream;

#include "MayDay.H"

#include "BaseIF.H"
#include "STLExplorer.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function reads an STL file and uses the polygonal information
    to provide edge intersections.  As such, calling it's "value" function is
    an error.  It is handled specially in "GeometryShop".
 */
class STLIF: public BaseIF
{
public:
  ///
  /**
      Type of data being read
   */
  enum DataType
  {
    Invalid = -1,
    ASCII   =  0,
    Binary      ,
    NUMDATATYPES
  };

  ///
  /**
      Constructor specifying filename (a_filename), the form of the data
      (a_dataType - ASCII or Binary), level set value (a_value), and whether
      inside the domain is where data is less than the level set value
      (a_inside).  Data is read from the file named and a complete ASCII
      header (see above) is expected.
   */
  STLIF(const char* const      a_filename,
        const STLIF::DataType& a_dataType);

  /// Copy constructor
  STLIF(const STLIF& a_inputIF);

  /// Destructor
  virtual ~STLIF();

  ///
  /**
      For STLIF, calling this method is an error.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual STLExplorer* getExplorer() const;

protected:
  void makeExplorer();

  string          m_filename;
  STLIF::DataType m_dataType;

  STLExplorer* m_explorer;

private:
  STLIF()
  {
    MayDay::Abort("STLIF uses strong construction");
  }

  void operator=(const STLIF& a_inputIF)
  {
    MayDay::Abort("STLIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
