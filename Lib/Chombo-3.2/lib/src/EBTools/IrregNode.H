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

#ifndef _IRREGNODE_H_
#define _IRREGNODE_H_

#include "REAL.H"
#include "IntVect.H"
#include "IntVectSet.H"
#include "VolIndex.H"
#include "Interval.H"
#include "FaceIndex.H"
#include "IntVectSet.H"
#include "RealVect.H"
#include "BaseFab.H"
#include "SPMD.H"
#include "NamespaceHeader.H"
///
/**
   Node for construction of geometric information.
*/
class IrregNode
{
public:
  ///
  IrregNode();

  ///
  ~IrregNode();

  ///return the index into the arc vectors
  static int index(int a_idir, Side::LoHiSide a_side) ;

  ///
  IntVect     m_cell;
  ///data for irregular nodes

  //volume centered data
  ///
  Real         m_volFrac;

  ///
  /**
     each vof has a unique index in the cell
  */
  int          m_cellIndex;

  ///
  RealVect     m_volCentroid;

  ///
  RealVect     m_bndryCentroid;

  //face centered irregular data

  ///
  /**
     Indicies into a_nodes to show connectivity.
     If the arc is to  an irregular
     cell, the index is the unique index of the vof
     in the cell.
     For arcs to regular cells, the arc = -2
     If the arc is to the domain boundary, arc  = -1.
  */
  Vector<int>      m_arc[2*SpaceDim];

  ///
  Vector<Real>     m_areaFrac[2*SpaceDim];

  ///
  Vector<RealVect> m_faceCentroid[2*SpaceDim];

  ///
  /**
      helper function for construction.  makes irregular cell
      that has connectivitity and values like a regular cell, this
      a person can modify as the irregular cell requires.  saves
      some coding in some cases
  */
  void makeRegular(const IntVect& iv);
  void faceReserve(int location, int size);
};


std::ostream& operator<< (std::ostream&  a_os,
                          const IrregNode& a_iv);

#include "NamespaceFooter.H"
#endif
