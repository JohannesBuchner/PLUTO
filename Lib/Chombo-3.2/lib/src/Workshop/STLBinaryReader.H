#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STLBINARYREADER_H_
#define _STLBINARYREADER_H_

#include <iostream>
#include <string>
using namespace std;

#include "MayDay.H"

#include "STLReader.H"
#include "STLMesh.H"
//#include "STLExplorer.H" // for comparators

#include "NamespaceHeader.H"

/*
 * Reads binary STL files and generates a mesh
 * see http://www.cplusplus.com/doc/tutorial/files/
 */

class STLBinaryReader: public STLReader
{
public:
  /// Constructor - read from standard input
  STLBinaryReader();

  /// Constructor - read from file name
  STLBinaryReader(const string& a_filename);

  /// Destructor
  ~STLBinaryReader();

  /// Return header information
  string* GetHeader() const;

  /// Return number of elements
  void GetNtri(int& a_ntri) const;

  /// Return whether number of elements from header matches file
  void GetNtriMatch(bool& a_ntriMatch) const;

  /// Return pointer to the mesh
  RefCountedPtr<STLMesh> GetMesh() const;

protected:
  void ReadData(istream&   a_file,
                const int offset);

  string* m_header;    // header info
  int     m_ntri;      // number of triangles read in
  bool    m_ntriMatch; // true if m_ntri equals the number of triangles declared in the header

  // actual data - shared by all copies
  RefCountedPtr<STLMesh> m_stlmesh; // pointer to the mesh

private:
  void operator=(const STLBinaryReader& a_inputReader)
  {
    MayDay::Abort("STLBinaryReader doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
