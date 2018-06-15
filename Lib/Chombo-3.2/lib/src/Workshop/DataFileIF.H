#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DATAFILEIF_H_
#define _DATAFILEIF_H_

#include <iostream>
#include <fstream>
using std::istream;
using std::ifstream;

#include "MayDay.H"
#include "RealVect.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "RefCountedPtr.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function reads data from a file and produces function values
    from the data read.  The data lies on a rectilinear grid where the number
    of grid points in each direction, the spacing, and the origin are specified
    (by the user or in the file).  Function values within the data region are
    obtained using trilinear interpolation.  Function values outside the data
    region are the maximum value of the data read.

    The data can be read from standard input or a named file.  The data can
    have no header, an ASCII header containing the number of grid points in each
    direction:
    <pre>
         numX numY [numZ]
    </pre>
    or an ASCII header containing the above information followed by the data
    spacing in each direction and the origin of the data - both in physical
    coordinates:
    <pre>
         deltaX  deltaY  [deltaZ]
         originX originY [originZ]
    </pre>
    Following any header is the data either in ASCII or binary.  The data is
    assumed to vary most rapidly in the x coordinate, then the y coordinate,
    and, finally, in the z coordinate (for 3D data).
 */
class DataFileIF: public BaseIF
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
      Constructor specifying the form of the data (a_dataType - ASCII or
      Binary), a level set value (a_value), and whether inside the domain
      is where data is less than the level set value (a_inside).  Data is
      read from standard input and a complete ASCII header (see above) is
      excepted.
   */
  DataFileIF(const DataFileIF::DataType& a_dataType,
             const Real&                 a_value,
             const bool&                 a_inside,
             const bool&                 a_useCubicInterp = false);

  ///
  /**
      Constructor specifying filename (a_filename), the form of the data
      (a_dataType - ASCII or Binary), level set value (a_value), and whether
      inside the domain is where data is less than the level set value
      (a_inside).  Data is read from the file named and a complete ASCII
      header (see above) is expected.
   */
  DataFileIF(const char* const           a_filename,
             const DataFileIF::DataType& a_dataType,
             const Real&                 a_value,
             const bool&                 a_inside,
             const bool&                 a_useCubicInterp = false);

  ///
  /**
      Constructor specifying the form of the data (a_dataType - ASCII or
      Binary), the spacing (a_spacing), the physical origin (a_origin), a
      level set value (a_value), and whether inside the domain is where
      data is less than the level set value (a_inside).  Data is read from
      standard input and an ASCII header (see above) containing the number
      of grid points in each direction is excepted.
   */
  DataFileIF(const DataFileIF::DataType& a_dataType,
             const RealVect&             a_spacing,
             const RealVect&             a_origin,
             const Real&                 a_value,
             const bool&                 a_inside,
             const bool&                 a_useCubicInterp = false);

  ///
  /**
      Constructor specifying filename (a_filename), the form of the data
      (a_dataType - ASCII or Binary), the spacing (a_spacing), the physical
      origin (a_origin), a level set value (a_value), and whether inside
      the domain is where data is less than the level set value (a_inside).
      Data is read from the file named and an ASCII header (see above)
      containing the number of grid points in each direction is excepted.
   */
  DataFileIF(const char* const           a_filename,
             const DataFileIF::DataType& a_dataType,
             const RealVect&             a_spacing,
             const RealVect&             a_origin,
             const Real&                 a_value,
             const bool&                 a_inside,
             const bool&                 a_useCubicInterp = false);

  ///
  /**
      Constructor specifying the form of the data (a_dataType - ASCII or
      Binary), the spacing (a_spacing), the physical origin (a_origin), the
      number of grid points in each direction (a_num), a level set value
      (a_value), and whether inside the domain is where data is less than
      the level set value (a_inside).  Data is read from standard input
      and no ASCII header (see above) is excepted.
   */
  DataFileIF(const DataFileIF::DataType& a_dataType,
             const IntVect&              a_num,
             const RealVect&             a_spacing,
             const RealVect&             a_origin,
             const Real&                 a_value,
             const bool&                 a_inside,
             const bool&                 a_useCubicInterp = false);

  ///
  /**
      Constructor specifying filename (a_filename), the form of the data
      (a_dataType - ASCII or Binary), the spacing (a_spacing), the physical
      origin (a_origin), the number of grid points in each direction (a_num),
      a level set value (a_value), and whether inside the domain is where
      data is less than the level set value (a_inside).  Data is read from
      the file named and no ASCII header (see above) is excepted.
   */
  DataFileIF(const char* const           a_filename,
             const DataFileIF::DataType& a_dataType,
             const IntVect&              a_num,
             const RealVect&             a_spacing,
             const RealVect&             a_origin,
             const Real&                 a_value,
             const bool&                 a_inside,
             const bool&                 a_useCubicInterp = false);

  /// Copy constructor
  DataFileIF(const DataFileIF& a_inputIF);

  ///
  /**
      This is used by the factory (see below) to create a new object.  All
      objects created in this way share a refcounted pointer to their data.

      Constructor specifying a refcounted pointer to the data (a_ascii_data
      or a_binary_data), the no data value (a_noDataValue), the spacing
      (a_spacing), the physical origin (a_origin), the number of grid points
      in each direction (a_num), a level set value (a_value), and whether
      inside the domain is where data is less than the level set value
      (a_inside).
   */
  DataFileIF(const RefCountedPtr<FArrayBox>               a_ascii_data,
             const RefCountedPtr<BaseFab<unsigned char> > a_binary_data,
             const Real&                                  a_noDataValue,
             const IntVect&                               a_num,
             const RealVect&                              a_spacing,
             const RealVect&                              a_origin,
             const Real&                                  a_value,
             const bool&                                  a_inside,
             const bool&                                  a_useCubicInterp = false);

  /// Destructor
  virtual ~DataFileIF();

  ///
  /**
      Return the header information
   */
  virtual void GetHeader(IntVect&  a_num,
                         RealVect& a_spacing,
                         RealVect& a_origin) const;

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real& a_value,
                         bool& a_inside,
                         bool& a_useCubicInterp) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real& a_value,
                         const bool& a_inside,
                         const bool& a_useCubicInterp = false);

  ///
  /**
     value to use when we are outside
  */
  virtual void SetNoDataValue(const Real& a_value);

  ///
  /**
      Return the value of the function at a_point using trilinear interpolation
      of the data.  If a_point is outside the data then return the maximum data
      value.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  void GetFullHeader(IntVect&  a_num,
                      RealVect& a_spacing,
                      RealVect& a_origin)
  {
    a_num = m_num;
    a_spacing = m_spacing;
    a_origin = m_origin;
  }

protected:
  void OpenFile(ifstream&         a_file,
                const char* const a_filename);

  void CloseFile(ifstream& a_file);

  void ReadMinHeader(IntVect& a_num,
                     istream& a_file);

  void ReadFullHeader(IntVect&  a_num,
                      RealVect& a_spacing,
                      RealVect& a_origin,
                      istream&  a_file);

  void ReadData(Real&                       a_maxValue,
                istream&                    a_file,
                const DataFileIF::DataType& a_dataType,
                const IntVect&              a_num);

  void MakeCorners(void);

  IntVect  m_num;       // number of grid points in each direction
  RealVect m_spacing;   // grid spacing (in physical coordinates)
  RealVect m_origin;    // grid origin (in physical coordinates)

  Real     m_value;     // level set value

  bool     m_inside;    // inside less than flag

  // the data - copies all share this
  RefCountedPtr<FArrayBox>               m_ascii_data;  // ASCII data stored as Real
  RefCountedPtr<BaseFab<unsigned char> > m_binary_data; // binary data stored as unsigned char's

  Real m_noDataValue;   // no data value

  IntVectSet m_cornersLinear; // corners of a 1x1 box for linear interpolation
  IntVectSet m_cornersCubic;  // corners of a 4x4 box for cubic  interpolation

  bool m_useCubicInterp;

private:
  DataFileIF()
  {
    MayDay::Abort("DataFileIF uses strong construction");
  }

  void operator=(const DataFileIF& a_inputIF)
  {
    MayDay::Abort("DataFileIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
