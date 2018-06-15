#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTG, ANAG, LBNL, 2009

#ifndef _MINIIFFAB_H_
#define _MINIIFFAB_H_

#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "FaceIndex.H"
#include "BaseFab.H"
#include "EBGraph.H"
#include "Arena.H"
#include "NamespaceHeader.H"


///
/**
   MiniIFFAB is a templated data holder
   defined over the Faces of an irregular domain.
*/
template <class T>
class MiniIFFAB
{
public:
  ///
  /**
     Default constructor.  Constructs an uninitialized IFFAB.  You must
     subsequently call {\bf define} before using this IFFAB.
  */
  MiniIFFAB();

  ///
  /**
     Defining constructor.  Specifies the irregular domain
     and the number of data components per face.  The irregular domain
     must lie completely within the EBGraph.  The
     contents are uninitialized.  Calls full define function.
  */
  MiniIFFAB(const Box& a_box,
            const EBGraph& a_ebgraph,
            const int& a_direction,
            const int& a_nvarin);

  ///
  ~MiniIFFAB();

  ///
  /**
     Full define function.  Specifies the irregular domain and the
     number of data components per face.  The irregular domain
     must lie completely within the EBGraph.   The contents
     are uninitialized.  If it has previously been defined, the old definition
     data is overwritten and lost.
  */
  void define(const Box& a_box,
              const EBGraph& a_ebgraph,
              const int& a_direction,
              const int& a_nvarin);

  ///for AggStencil
  long offset(const FaceIndex& a_vof, const int& a_ivar) const
  {
    const T* conval = &((*this)(a_vof, a_ivar));
    long  roffset = (long)(conval-dataPtr(0));
    return roffset;
  }

  ///
  /**
     Set a value at every data location.
  */
  void setVal(const T& value);

  ///
  /**
     Set a value at every data location.
  */
  void setVal(int ivar, const T& value);

  ///
  /**
     define over surrounding nodes of invectset (ivs is cell-centered).
   */
  void copy(const Box& a_intBox,
            const Interval& a_destInterval,
            const Box& a_toBox,
            const MiniIFFAB<T>& a_src,
            const Interval& a_srcInterval);

  ///
  /**
     Remove all data from this MiniIFFAB.
     You must call {\bf define} again in order
     to use it.
  */
  void clear();

  ///invalid but necessary for leveldata to compile
  MiniIFFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for baseebcellfab");
  }

  ///
  /**

     Tells whether this MiniIFFAB has been defined, either with a constructor
     or with {\bf define}.  It must be initialized in order to access its data.
  */
  bool
  isDefined() const
  {
    return (m_isDefined);
  }


  ///
  /**
     Return the number of data components of this MiniIFFAB.
  */
  int nComp() const
  {
    return m_nComp;
  }

  ///
  /**
     Indexing operator.  Return a reference to the contents of this IFFAB,
     at the specified face and data component.  The first component is
     zero, the last is {\em nvar-1}.  The returned object is a modifiable
     lvalue.
  */
  T& operator() (const FaceIndex& a_face,const int& varlocin);
  const T& operator()(const FaceIndex& a_face,const int& varlocin) const;


  ///
  T* dataPtr(const int& a_comp);

  T* getIndex(const FaceIndex& a_face, const int& a_comp) const;
  ///
  const T* dataPtr(const int& a_comp) const;

  static int preAllocatable()
  {
    return 1; // symmetric preallocatable
  }

  ///
  int size(const Box& R, const Interval& comps) const ;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///
  int numFaces() const
  {
    return m_nFaces;
  }

  ///
  const Vector<FaceIndex>& getFaces() const
  {
    return m_faces;
  }

  static bool s_verbose;

  Box m_box;
  EBGraph m_ebgraph;

private:
  void   setDefaultValues();

  void   getFaceSubset(Vector<FaceIndex>& a_subFaces, const Box& a_subBox) const;
protected:

  T* m_data;
  Vector<FaceIndex> m_faces;
  int m_nFaces;
  int m_nComp;
  int m_direction;

  bool    m_isDefined;
private:
  //we disallow copy construction and assigment for data holders
  void operator= (const MiniIFFAB<T>&)
  {
    MayDay::Error("invalid operator");
  }
  MiniIFFAB (const MiniIFFAB<T>&)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "MiniIFFABI.H"
#endif

#endif
