#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTG, ANAG, LBNL

#ifndef _BASEIFFAB_H_
#define _BASEIFFAB_H_

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
   BaseIFFAB is a templated data holder
   defined over the Faces of an irregular domain.
*/
template <class T>
class BaseIFFAB
{
public:
  ///
  /**
     Default constructor.  Constructs an uninitialized IFFAB.  You must
     subsequently call {\bf define} before using this IFFAB.
  */
  BaseIFFAB();

  ///
  /**
     Defining constructor.  Specifies the irregular domain
     and the number of data components per face.  The irregular domain
     must lie completely within the EBGraph.  The
     contents are uninitialized.  Calls full define function.
  */
  BaseIFFAB(const IntVectSet& a_region,
            const EBGraph& a_ebgraph,
            const int& a_direction,
            const int& a_nvarin);

  ///
  ~BaseIFFAB();

  ///
  /**
     Full define function.  Specifies the irregular domain and the
     number of data components per face.  The irregular domain
     must lie completely within the EBGraph.   The contents
     are uninitialized.  If it has previously been defined, the old definition
     data is overwritten and lost.
  */
  void define(const IntVectSet& a_region,
              const EBGraph& a_ebgraph,
              const int& a_direction,
              const int& a_nvarin);

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
            const BaseIFFAB<T>& a_src,
            const Interval& a_srcInterval);

  ///
  /**
     Remove all data from this BaseIFFAB.
     You must call {\bf define} again in order
     to use it.
  */

  void clear();

  ///invalid but necessary for leveldata to compile
  BaseIFFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for baseebcellfab");
  }

  ///
  /**

     Tells whether this BaseIFFAB has been defined, either with a constructor
     or with {\bf define}.  It must be initialized in order to access its data.
  */
  bool
  isDefined() const;

  ///
  /**
     Return the number of face in this BaseIFFAB.
  */
  int numFaces() const;

  ///
  /**
     Return the number of data components of this BaseIFFAB.
  */
  int nComp() const;

  ///

  ///
  /**
     Return the direction of the faces of this BaseIFFAB.
  */
  int direction() const;

  ///

  /**
     Return  the  irregular domain of the BaseIFFAB.
  */
  const IntVectSet& getIVS() const;

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

  ///for AggStencil
  long offset(const FaceIndex& a_vof, const int& a_ivar) const
  {
    const T& conval = (*this)(a_vof, a_ivar);
    long roffset = (int)(&conval);
    roffset -= (int)(&m_data[0]);
    return roffset;
  }

  ///
  const T* dataPtr(const int& a_comp) const;

  const EBGraph& getEBGraph() const;

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
  static void setVerbose(bool a_verbose);
  static bool s_verbose;
private:
  void
  setDefaultValues();

protected:

  //functions to aid in my ongoing battle with
  //cell-centered semantics for a face-centered
  //data structure.
  bool useThisFace(const Box& a_toBox,
                   const FaceIndex& a_face) const;

  void getBoxAndIVS(Box&        a_intBox,
                    IntVectSet& a_ivsIntersect,
                    const Box&  a_toBox) const;


  static Arena* s_Arena;

  //get index vector<int> at the face location
  //using lexographic ordering of vof's indicies
  int getLocalVecIndex(const FaceIndex& a_face) const;

public:
  static void   setSurroundingNodeSemantic(bool a_useSurr);

protected:
  static bool  s_doSurroundingNodeSemantic;
  //T* m_data;
  Vector<T> m_data;
  int m_truesize;
  int m_nComp;
  int m_nFaces;
  int m_direction;

  //object to store pointers to facilitate fast indexing
  BaseFab<T*> m_fab;
  IntVectSet m_ivs;
  EBGraph m_ebgraph;
  bool m_isDefined;
private:
  //we disallow copy construction and assigment for data holders
  void operator= (const BaseIFFAB<T>&)
  {
    MayDay::Error("invalid operator");
  }
  BaseIFFAB (const BaseIFFAB<T>&)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "BaseIFFABI.H"
#endif

#endif
