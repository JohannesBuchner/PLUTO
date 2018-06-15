#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ARRAYVIEWDATA_H_
#define _ARRAYVIEWDATA_H_

#include "REAL.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "LayoutData.H"
#include "NamespaceHeader.H"

/// Special LayoutData type for use with ArrayView
/** ArrayViewData is a wrapper around LayoutData for BaseFab<Real> or
 *  FArrayBox for use in ArrayView.  It provides random access to the
 *  individual FABs and the ability to directly set the data in them.  It is
 *  intended to work only in ArrayView.  Any other attept to use this class
 *  is unsupported.
 */

class ArrayViewData
{
public:
  ArrayViewData()
    :
    m_layoutdata_ptr(NULL)
  {}

  ~ArrayViewData()
  {}

  ArrayViewData(LayoutData<BaseFab<Real> >* a_layoutdata_ptr)
    :
    m_layoutdata_ptr(a_layoutdata_ptr)
  {}

  ArrayViewData(LayoutData<FArrayBox>* a_layoutdata_ptr)
    :
    m_layoutdata_ptr((LayoutData<BaseFab<Real> >*) a_layoutdata_ptr)
  {}

  /// Access a single FAB in the LayoutData using its internal index
  BaseFab<Real>& operator[](int a_elem)
  {
    return *(m_layoutdata_ptr->m_vector[a_elem]);
  }

  /// Substitute a different FAB into the LayoutData
  /** The caller is responsible for managing the memory used by the FAB that
   * is being replaced.
   */
  void set(int            a_elem,
           BaseFab<Real>* a_fab_ptr)
  {
    m_layoutdata_ptr->m_vector[a_elem] = a_fab_ptr;
  }


private:
  LayoutData<BaseFab<Real> >* m_layoutdata_ptr;
};

#include "NamespaceFooter.H"
#endif
