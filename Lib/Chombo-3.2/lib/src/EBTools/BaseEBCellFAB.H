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

#ifndef _BASEEBCELLFAB_H_
#define _BASEEBCELLFAB_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISBox.H"
#include "BaseIVFAB.H"
#include "IntVectSet.H"
#include "BaseFab.H"
#include "MiniIVFAB.H"
#include "NamespaceHeader.H"

///  Array defined at the VolIndexs of an Box in an EBIS.
/**
   An BaseEBCellFAB is a templated holder for cell-centered
   data over the intersection of a box and an EbisBox.
   At every uncovered VolIndex in this intersection,
   the BaseEBCellFAB contains a specified number of data values.
   At singly valued cells,
   the data is stored internally in a BaseFab.  At multiply-valued
   cells, the data is stored internally in an  BaseIVFAB.
   BaseEBCellFAB provides indexing by VolIndex.
   This class does not provide a copy constructor or assignment operator
   for all the usual  reasons.
*/
template <class T>
class BaseEBCellFAB
{
public:

  /// for AggStencil
  int numDataTypes() const
  {
    return 2;
  }

  /// for AggStencil
  int   dataType(const BaseIndex& a_baseInd) const
  {
    const VolIndex* vofPtr = dynamic_cast<const VolIndex *>(&a_baseInd);
    if (vofPtr == NULL) MayDay::Error("Trying to index into an EBCellFAB with something that is not a VolIndex");
    int retval = 0;
    if (m_multiCells.contains(vofPtr->gridIndex()))
      {
        retval = 1;
      }
    return retval;
  }

  ///for AggStencil
  long  offset(const BaseIndex& a_baseInd, const int a_ivar) const
  {
    const VolIndex* vofPtr = dynamic_cast<const VolIndex *>(&a_baseInd);
    if (vofPtr == NULL) MayDay::Error("Trying to index into an EBCellFAB with something that is not a VolIndex");
    long retval = 0;
    if (m_multiCells.contains(vofPtr->gridIndex()))
      {
        retval = m_irrFAB.offset(*vofPtr, a_ivar);
      }
    else
      {
        retval = m_regFAB.offset(vofPtr->gridIndex(), a_ivar);
      }
    return retval;
  }

  ///for AggStencil
  T* dataPtr(int a_dataType, int a_ivar)
  {
    if (a_dataType == 0)
      {
        return m_regFAB.dataPtr(a_ivar);
      }
    else if (a_dataType == 1)
      {
        return m_irrFAB.dataPtr(a_ivar);
      }
    else
      {
        MayDay::Error("bogus datatype input into dataPtr");
      }
    return NULL; //just avoiding compiler warnings with this
  }
  ///for AggStencil
  const T* dataPtr(int a_dataType, int a_ivar) const
  {
    if (a_dataType == 0)
      {
        return m_regFAB.dataPtr(a_ivar);
      }
    else if (a_dataType == 1)
      {
        return m_irrFAB.dataPtr(a_ivar);
      }
    else
      {
        MayDay::Error("bogus datatype input into dataPtr");
      }
    return NULL; //just avoiding compiler warnings with this
  }

  ///
  BaseEBCellFAB();

  ///
  BaseEBCellFAB(const EBISBox& a_ebisBox, const Box& a_region, int a_nVar);

  ///
  virtual void  define(const EBISBox& a_ebisBox, const Box& a_region, int a_nVar);

  /// aliasing copy-ish constructor.  Used when the 'alias' function is invoked on LevelData
  /**
     Constructs an 'aliased' BaseEBCellFab of the requested interval of the
     argument BaseEBCellFab.  This BaseEBCellFab does not allocate any memory, but
     sets its data pointer into the memory pointed to by the argument
     BaseEBCellFab.  It is the users responsiblity to ensure this aliased
     BaseEBCellFab is not used after the original BaseEBCellFab has deleted its data ptr
     (resize, define(..) called, or destruction, etc.).

     This aliased BaseEBCellFab will also generate side effects (modifying the values
     of data in one will modify the other's data).

     This aliased BaseFab will have a_comps.size() components, starting at zero.
  */
  BaseEBCellFAB(const Interval&   a_comps,
                BaseEBCellFAB<T>& a_original);

  ///
  virtual ~BaseEBCellFAB();

  ///
  void clear();

  ///
  void  setVal(const T& value);

  ///
  void  setVal(int ivar,const T& value);

  ///
  void  setVal(const T& value,
               const Box& a_box,
               int        a_nstart,
               int        a_numcomp);

  ///
  bool isDefined() const;

  ///
  int nComp() const ;

  ///
  const Box& getRegion() const ;

  ///
  const Box& box() const ;

  ///
  const IntVectSet& getMultiCells() const ;

  ///
  const BaseIVFAB<T>& getMultiValuedFAB() const;

  ///
  BaseIVFAB<T>& getMultiValuedFAB() ;

  ///
  const BaseFab<T>& getSingleValuedFAB() const;

  ///
  BaseFab<T>& getSingleValuedFAB();

  /**
  Arg a_isKnownMultiValued should be set to:
     -1 if arg a_ndin is known to be in a single-valued cell;
     +1 if arg a_ndin is known to be in a multi-valued cell;
      0 if you're not sure (in which case the expensive IntVectSet::contains()
        function gets called.

  Note: BVS says third arg is a wart that makes this operator look less like
        a logical array access.
  */
  const T& operator() (const VolIndex& a_ndin, int a_nVarLoc,
                       int a_isKnownMultiValued=0) const;

  /**
  Arg a_isKnownMultiValued should be set to:
     -1 if arg a_ndin is known to be in a single-valued cell;
     +1 if arg a_ndin is known to be in a multi-valued cell;
      0 if you're not sure (in which case the expensive IntVectSet::contains()
        function gets called.
  */
  T& operator() (const VolIndex& a_ndin,int  a_nVarLoc,
                 int a_isKnownMultiValued=0);

  void fill(T* array, const VolIndex& a_ndin, const Interval& a_comps) const;

  void assign(const T* array, const VolIndex& a_ndin, const Interval& a_comps);
  ///
  void copy(const Box& RegionFrom,
            const Interval& destInt,
            const Box& RegionTo,
            const BaseEBCellFAB<T>& source,
            const Interval& srcInt);

  ///
  BaseEBCellFAB<T>& copy(const BaseEBCellFAB<T>& source);

  ///
  /**
     This stuff required by LevelData in parallel:
  */
  static int preAllocatable()
  {
    return 1; // symmetric messaging
  }

  ///
  int size(const Box& R, const Interval& comps) const ;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const ;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///invalid but necessary for leveldata to compile
  BaseEBCellFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for baseebcellfab");
  }

  ///
  const EBISBox& getEBISBox() const;

  ///
  void
  setCoveredCellVal(const T&    a_val,
                    const int&  a_comp,
                    const bool& a_doMulti=true);

protected:

  ///data at multi-valued cells
  //BaseIVFAB<T> m_irrFAB;
  MiniIVFAB<T> m_irrFAB;

  ///data at single-valued cells
  BaseFab<T> m_regFAB;

  EBISBox m_ebisBox;
  //number of data values at each VolIndex
  int m_nComp;

  ///region over which the BaseEBCellFAB exists
  Box m_region;

  ///whether each given cell is multi-valued
  IntVectSet m_multiCells;

  ///is m_multiCells not empty?
  bool m_hasMultiCells;

  ///has full define function been called?
  bool m_isDefined;

protected:
  void
  setDefaultValues();

private:
  //disallowed for perfomance reasons
  void operator= (const BaseEBCellFAB<T>& ebcin)
  {
    MayDay::Error("operator= for baseebcellfab undefined");
  }
  BaseEBCellFAB (const BaseEBCellFAB<T>& ebcin)
  {
    MayDay::Error("copy constructor for baseebcellfab undefined");
  }
};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "BaseEBCellFABI.H"
#endif // CH_EXPLICIT_TEMPLATES

#endif
