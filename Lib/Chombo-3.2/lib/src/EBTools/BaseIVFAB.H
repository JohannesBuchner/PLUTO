#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEIVFAB_H_
#define _BASEIVFAB_H_

#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "VolIndex.H"
#include "BaseIndex.H"
#include "BaseFab.H"
#include "EBGraph.H"
#include "NamespaceHeader.H"


///
/**
   BaseIVFAB is a templated
   data holder defined at the VoFs of an irregular domain.
*/
template <class T>
class BaseIVFAB
{
public:
  ///
  /**
     Default constructor.  Constructs an uninitialized IVFAB.  You must
     subsequently call {\bf define} before using this IVFAB.
  */
  BaseIVFAB();

  ///
  /**
     Defining constructor.  Specifies the irregular domain
     and the number of data components per VoF.  The irregular domain
     must lie completely within the EBGraph.  The
     contents are uninitialized.  Calls full define function.
  */
  BaseIVFAB(const IntVectSet& a_region,
            const EBGraph&    a_ebgraph,
            const int&        a_nvarin);

  /// for AggStencil
  int numDataTypes() const
  {
    return 1;
  }
  /// for AggStencil
  int   dataType(const BaseIndex& a_baseInd) const
  {
    return 0;
  }

  T* dataPtr(int a_dataType, int a_ivar)
  {
    return dataPtr(a_ivar);
  }

  const T* dataPtr(int a_dataType, int a_ivar) const
  {
    return dataPtr(a_ivar);
  }



  ///
  /**
     Constructs an 'aliased' BaseIVFAB of the requested interval of the
     argument BaseIVFAB.  This BaseIVFAB does not allocate any memory, but
     sets its data pointer into the memory pointed to by the argument
     BaseIVFAB.  It is the users responsiblity to ensure this aliased
     BaseIVFAB is not used after the original BaseIVFAB has deleted its data ptr
     (resize, define(..) called, or destruction, etc.).

     This aliased BaseIVFAB will also generate side effects (modifying the values
     of data in one will modify the other's data).

     This aliased BaseIVFAB will have a_comps.size() components, starting at zero.
  */
  BaseIVFAB(const Interval& a_comps,
            BaseIVFAB<T>&   a_original);


  ///
  virtual ~BaseIVFAB();

  ///
  /**
     Full define function.  Specifies the irregular domain and the
     number of data components per VoF.  The irregular domain
     must lie completely within the EBGraph.   The contents
     are uninitialized.  If it has previously been defined, the old definition
     data is overwritten and lost.
  */
  virtual void define(const IntVectSet& a_region,
              const EBGraph&    a_ebgraph,
              const int&        a_nvarin);

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
     Set a value at every data location within a_box
  */
  void setVal(const T&   a_value,
              const Box& a_box,
              int        a_nstart,
              int        a_numcomp);
  ///
  void copy(const Box& a_fromBox,
            const Interval& a_destInterval,
            const Box& a_toBox,
            const BaseIVFAB<T>& a_src,
            const Interval& a_srcInterval);

   ///
  /**
         This stuff required by LevelData in parallel:
  */

  static int preAllocatable()
  {
    return 1; // symmetric allocatable.
  }

  ///
  int size(const Box& R, const Interval& comps) const ;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///
  /**
     Remove all data from this BaseIVFAB.
     You must call {\bf define} again in order
     to use it.
  */
  virtual void clear();

  ///
  /**

     Tells whether this BaseIVFAB has been defined, either with a constructor
     or with {\bf define}.  It must be initialized in order to access its data.
  */
  bool
  isDefined() const;

  ///
  /**
     Return the number of VoFs in this BaseIVFAB.
  */
  int numVoFs() const;

  ///
  /**
     Return the number of data components of this BaseIVFAB.
  */
  int nComp() const;

  ///
  /**
     Return  the  irregular domain of the BaseIVFAB.
  */
  const IntVectSet& getIVS() const;

  ///
  const EBGraph& getEBGraph() const;

  ///
  /**
     Indexing operator.  Return a reference to the contents of this IVFAB,
     at the specified VoF and data component.  The first component is
     zero, the last is {\em nvar-1}.  The returned object is a modifiable
     lvalue.
  */
  T& operator() (const VolIndex& a_vof,const int& varlocin);
  const T& operator() (const VolIndex& a_vof,const int& varlocin) const;


  ///invalid but necessary for leveldata to compile
  BaseIVFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for baseebcellfab");
  }

  ///
  static void setVerbose(bool a_verbose);

  ///
  static void setVerboseDebug(bool a_verboseDebug);

  ///
  T* dataPtr(const int& a_comp);

  ///
  const T* dataPtr(const int& a_comp) const;

  ///get index into vector
  /**
     needs to be public so that bulk stencils can be constructed
   */
  virtual T* getIndex(const VolIndex& a_vof,const int& a_comp) const;

  ///for AggStencil
  long offset(const BaseIndex& a_vof, const int& a_ivar) const
  {
    const VolIndex* vof = dynamic_cast< const VolIndex* >(&a_vof);
    if (vof == NULL) MayDay::Error("cast failed:BaseIVFAB only takes vofs for indexing");

    const T& conval = (*this)(*vof, a_ivar);
    const T* conptr = &conval;
    const T* srtptr = dataPtr(0);
    long roffset = (long)(conptr - srtptr);
    return roffset;
  }
  static bool s_verboseDebug;

private:


protected:

  virtual void
  setDefaultValues();

  //has to be this in case someone does a bool
  //T* m_data;
  Vector<T> m_data;

  int m_nComp;
  int m_nVoFs;

  EBGraph m_ebgraph;

  //object to store pointers to facilitate fast indexing
  BaseFab<T*> m_fab;
  IntVectSet m_ivs;
  bool m_isDefined;

  static bool s_verbose;
private:
  //disallowed for all the usual reasons
  void operator= (const BaseIVFAB<T>& a_input)
  {
    MayDay::Error("BaseIVFAB operator = not defined");
  }
  BaseIVFAB(const BaseIVFAB<T>& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "BaseIVFABI.H"
#endif

#endif
