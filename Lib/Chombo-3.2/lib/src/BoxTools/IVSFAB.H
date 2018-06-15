#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IVSFAB_H_
#define _IVSFAB_H_

#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "IntVect.H"
#include "BaseFab.H"
#include "LayoutData.H"
#include "BoxLayoutData.H"
#include "NamespaceHeader.H"

///
/**
   IVSFAB is a templated
   data holder defined at the Ivs of an irregular domain.
*/
template <class T>
class IVSFAB
{
public:
  ///
  /**
     Default constructor.  Constructs an uninitialized IVFAB.  You must
     subsequently call {\bf define} before using this IVFAB.
  */
  IVSFAB();

  ///
  /**
  */
  IVSFAB(const IntVectSet& a_region,
            const int&        a_nvarin);

  ///
  ~IVSFAB();

  ///
  /**
  */
  void define(const IntVectSet& a_region,
              const int&        a_nvarin);

  ///
  /**
     Set a value at every data location.
  */
  void setVal(const T& value);

  ///
  void copy(const Box& a_fromBox,
            const Interval& a_destInterval,
            const Box& a_toBox,
            const IVSFAB<T>& a_src,
            const Interval& a_srcInterval);

   ///
  /**
         This stuff required by LevelData in parallel:
  */

  static int preAllocatable()
  {
    return 2; // symmetric allocatable.
  }

  ///
  int size(const Box& R, const Interval& comps) const ;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///
  /**
     Remove all data from this IVSFAB.
     You must call {\bf define} again in order
     to use it.
  */
  void clear();

  ///
  /**

     Tells whether this IVSFAB has been defined, either with a constructor
     or with {\bf define}.  It must be initialized in order to access its data.
  */
  bool
  isDefined() const;

  ///
  /**
     Return the number of Ivs in this IVSFAB.
  */
  int numIvs() const;

  ///
  /**
     Return the number of data components of this IVSFAB.
  */
  int nComp() const;

  ///
  /**
     Return  the  irregular domain of the IVSFAB.
  */
  const IntVectSet& getIVS() const;

  ///
  /**
     Indexing operator.  Return a reference to the contents of this IVFAB,
     at the specified Iv and data component.  The first component is
     zero, the last is {\em nvar-1}.  The returned object is a modifiable
     lvalue.
  */
  T& operator() (const IntVect& a_iv,const int& varlocin);
  const T& operator() (const IntVect& a_iv,const int& varlocin) const;

  ///
  /**
     Return a const pointer to the internal data storage of this IVFAB. \\
  */
  const T* dataPtr(const int& a_comp) const;

  ///
  /**
     Return a pointer to the internal data storage of this IVFAB.
  */
  T* dataPtr(const int& a_comp) ;

  ///invalid but necessary for leveldata to compile
  IVSFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called ");
  }

  ///
  /**
     Return a fake "lower left corner of the domain" for use with Fortran.
  */
  const int* loVect() const;

  ///
  /**
     Return a fake "upper right corner of the domain" for use with Fortran.
  */
  const int* hiVect() const;

  ///
  static void setVerbose(bool a_verbose);
private:
  void
  setDefaultValues();

protected:
  //get index into dataPtr
  int getIndex(const IntVect& a_iv,const int& a_comp) const;
  //has to be a pointer because we are planning to ship
  //the damned  thing to fortran.
  T* m_dataPtr;
  int m_nComp;
  int m_nIvs;
  IntVect m_loVect;
  IntVect m_hiVect;

  //this is how we map from vof to vector of data.
  //grab the vof and its intvect tells you the vector<int>
  //the integer in the vof indexes into the vof.
  //that tells you the offset into m_dataPtr
  //you further have to account for the variable number.
  //the getIndex() function encapsulates the process
  BaseFab<int> m_ivmap;
  IntVectSet m_ivs;
  bool m_isDefined;

  static bool s_verbose;
private:
  //disallowed for all the usual reasons
  IVSFAB<T>& operator= (const IVSFAB<T>&);
  IVSFAB (const IVSFAB<T>&);
};

/// Factory class to produce IVSFABs.
/**
   Factory class to produce IVSFABs.
   This is needed for LevelData construction.
 */
template <class T>
class IVSFABFactory
  : public DataFactory< IVSFAB<T> >
{
public:

  /// factory function.
  /**
      Creates a new baseivfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer
      is passed to the user.
      Intersects a_sets[a_dit] from construction
      with input box because leveldata will expect
      that anyway when it tries to do a copy or a linearization.
  */
  virtual IVSFAB<T>* create(const Box& a_box, int a_ncomps,
                            const DataIndex& a_dit) const
  {
    //intersects with input box because leveldata will expect
    //that anyway when it tries to do a copy or a linearization.
    IntVectSet ivsInter = m_sets[a_dit];
    ivsInter &= a_box;
    return new IVSFAB<T>(ivsInter,  a_ncomps);
  }

  ///
  /**
     fills internal data with inputs
  */
  IVSFABFactory(const LayoutData<IntVectSet>& a_sets)
  {
    m_sets.define(a_sets.boxLayout());
    for (DataIterator dit = m_sets.dataIterator(); dit.ok(); ++dit)
      {
        m_sets[dit()] = a_sets[dit()];
      }
  }

  ///
  virtual ~IVSFABFactory()
  {
    ;
  }

private:

  ///
  LayoutData<IntVectSet> m_sets;

  //disallowed constructors and all that
  IVSFABFactory(const IVSFABFactory<T>& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  IVSFABFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const IVSFABFactory<T>& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#include "NamespaceFooter.H"
#include "IVSFABI.H"

#endif
