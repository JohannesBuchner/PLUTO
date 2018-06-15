#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOXLAYOUTDATA_H_
#define _BOXLAYOUTDATA_H_

#include "LayoutData.H"
#include "Interval.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "Copier.H"
#include "SPMD.H"
#include "memtrack.H"
#include "NamespaceHeader.H"

//--Forward declarations

class FluxBox;

extern int LinearizationTest;

/// Factory object to data members of a BoxLayoutData container
template <class T> class DataFactory
{
public:

  virtual ~DataFactory()
  {
  }

  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual T* create(const Box& box, int ncomps, const DataIndex& a_datInd) const=0;

  virtual bool callDelete() const
  {
    return true;
  }
};

/// Factory object to data members of a BoxLayoutData container
/**
   creates new T object by invoking operator new

   \code
   template <class T>
   T* DefaultDataFactory<T>::create(const Box& box,
   int ncomps,
   const DataIndex& a_datInd) const
   {
   return new T(box, ncomps);
   }
   \endcode
*/

template <class T> class DefaultDataFactory : public DataFactory<T>
{
public:
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual T* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

};

class FABAliasDataFactory : public DataFactory<FArrayBox>
{
public:
  virtual ~FABAliasDataFactory()
  {}

  FABAliasDataFactory(const LayoutData<Real*>& aliases);
  void define(const LayoutData<Real*>& aliases);
  /// factory function.  creates a new 'T' object using an aliased dataPtr for T
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual FArrayBox* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

protected:
  LayoutData<Real*> aliasPtrs;
};

template<class T> class BoxLayoutData;

template <class T>
class AliasDataFactory : public DataFactory<T>
{
public:
  virtual ~AliasDataFactory()
  {}

  AliasDataFactory(BoxLayoutData<T>* a_original, const Interval& interval);
  void define(BoxLayoutData<T>* a_original, const Interval& interval);
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual T* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

protected:
  BoxLayoutData<T>* m_origPointer;
  Interval          m_interval;
};

//--For aliasing a LevelData<FArrayBox> to a LevelData<FluxBox>

class FABAliasFlBxDataFactory : public DataFactory<FArrayBox>
{
public:
  virtual ~FABAliasFlBxDataFactory()
  {}

  FABAliasFlBxDataFactory(BoxLayoutData<FluxBox>* a_original,
                          const Interval&         a_interval,
                          const int               a_dir);
  void define(BoxLayoutData<FluxBox>* a_original,
              const Interval&         a_interval,
              const int               a_dir);
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual FArrayBox* create(const Box&       a_box,
                            int              a_ncomps,
                            const DataIndex& a_datInd) const;
protected:
  BoxLayoutData<FluxBox> *m_origPointer;
  Interval                m_interval;
  int                     m_dir;
};

template<class T> class LevelData;

template <class T>
class LDOperator
{
public:
  // to keep the compiler happy
  virtual ~LDOperator()
  {
  }

  virtual int size(const T& arg, const Box& b, const Interval& comps) const
  {
    return arg.size(b, comps);
  }

  virtual void linearOut(const T& arg, void* buf, const Box& R,
                         const Interval& comps) const
  {
    arg.linearOut(buf, R, comps);
  }

  virtual void linearIn(T& arg,  void* buf, const Box& R,
                        const Interval& comps)const = 0;

  virtual void op(T& dest,
                  const Box& RegionFrom,
                  const Interval& Cdest,
                  const Box& RegionTo,
                  const T& src,
                  const Interval& Csrc) const = 0;
};

template <class T>
class CopyLDOperator: public LDOperator<T>
{
public:
  // to keep the compiler happy
  virtual ~CopyLDOperator()
  {
  }

  virtual void linearOut(const T& arg, void* buf, const Box& R,
                         const Interval& comps) const
  {
    arg.linearOut(buf, R, comps);
  }
  virtual void linearIn(T& arg,  void* buf, const Box& R,
                        const Interval& comps)const
  {
    arg.linearIn(buf, R, comps);
  }
  virtual void op(T& dest,
                  const Box& RegionFrom,
                  const Interval& Cdest,
                  const Box& RegionTo,
                  const T& src,
                  const Interval& Csrc) const
  {
    if (LinearizationTest == 0)
      {
        dest.copy(RegionFrom, Cdest,RegionTo, src, Csrc);
      }
    else
      {
        int sizeSource = this->size(src,  RegionFrom, Csrc);
        int sizeDest   = this->size(dest, RegionTo,   Cdest);

        if (T::preAllocatable() == 0 || T::preAllocatable() == 1)
          {
            if (sizeSource != sizeDest)
              {
                MayDay::Abort("LinearizationTest failure: dest and source have different sizes");
              }
          }
        Vector<char> buffer(sizeSource);
        void* b = (void*)&(buffer[0]);
        linearOut(src, b, RegionFrom, Csrc);
        linearIn(dest, b, RegionTo,   Cdest);
      }
  }
};

/// Data on a BoxLayout
/**

   A collection of Box-oriented objects.  The arrangement
   of Boxes is given by the underlying BoxLayout object.
   BoxLayoutData attempt to prevent users from maniupulating
   the template class T to become out of synch with the
   boxes in the BoxLayout.  Caveat emptor.

   All access to the data in a BoxLayoutData is forced to be
   data-parallel, by use of the DataIterator.

   class T must provide the following methods:
   - <PRE> T() </PRE>
   - <PRE> T(const Box& box, int comps) </PRE>
   - <PRE> define(const Box& box, int comps) </PRE>
   - <PRE> void copy(const Box& R, const Interval& Cd, const T& source, const Interval Cs)
   Cs and Cd must have same length, but can be displaced </PRE>

   - <PRE> static int preAllocatable()
   returns 0 if the size(...) function is strictly a function of Box and
   Interval, and does not depend on the current state of the T object.
   return 1 if size(..) is symmetric, in that sender and receiver T object
   can size their message buffers, but a static object cannot.
   returns 2 if the object is truly dynamic.  the message size is subject to
   unique object data. </PRE>
   - <PRE> int size(const Box& R, const Interval& comps) </PRE>
   - <PRE> void linearOut(void* buf, const Box& R, const Interval& comps) </PRE>
   - <PRE> void linearIn(void* buf, const Box& R, const Interval& comps) </PRE>

   In sequential codes, many of these functions will not be called, but in the
   general parallel case, they all will be.

*/

template<class T>
class BoxLayoutData : public LayoutData<T>
{
public:
  ///
  BoxLayoutData();

  virtual ~BoxLayoutData();
  ///
  BoxLayoutData(const BoxLayout& boxes, int comps,
                const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  virtual void define(const BoxLayout& boxes, int comps,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  virtual void define(const BoxLayoutData<T>& da,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  /** Define this BoxLayoutData to be the subset of the
      data in da defined by the Interval comps.
  */
  virtual void define(const BoxLayoutData<T>& da, const Interval& comps,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  /// overridden and throws an error.
  virtual void define(const BoxLayout& boxes);

  ///
  int nComp() const
  {
    return m_comps;
  }

  ///
  Interval interval() const
  {
    Interval outint(0, m_comps-1);
    return(outint);
  }

  /// General data copying operation.
  /**
     @param a_destGrids BoxLayout for the destination data holder
     @param a_dest return argument.  upon completion contains a Vector of T objects
     for each Box in a_destGrids that overlaps this->boxLayout()
     @param a_interval range of components to source
     @param a_domain ProblemDomain that this grid and a_destGrids reside on.
     @param factory optional data factory for template classes that have non-trivial construction
  */
  void generalCopyTo(const BoxLayout& a_destGrids,
                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                     const Interval& a_interval,
                     const ProblemDomain& a_domain,
                     const DataFactory<T>& factory = DefaultDataFactory<T>()) const ;

  ///
  void generalCopyTo(const BoxLayout& a_destGrids,
                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                     const Interval& a_interval,
                     const ProblemDomain& a_domain,
                     const Copier& a_copier,
                     const DataFactory<T>& factory = DefaultDataFactory<T>()) const ;

  ///
  /**
     Special version of generalCopyTo that performs increment-on-intersection of
     the destination instead of replacement
  */
  void addTo(const Interval& a_srcComps,
             BoxLayoutData<T>& a_dest,
             const Interval& a_destComps,
             const ProblemDomain& a_domain) const;
  ///
  /**
     Special version of generalCopyTo that performs increment-on-intersection of
     the destination instead of replacement.  This version lets the user send in a
     pre-built Copier.
  */
  void addTo(const Interval& a_srcComps,
             BoxLayoutData<T>& a_dest,
             const Interval& a_destComps,
             const ProblemDomain& a_domain,
             const Copier& a_copier) const;

  ///
  /* User writes a function with the signature:

     <PRE>
     void myfunction(const Box& box, int comps, T& t)
     {
     your code here;
     }
     </PRE>

     They can then hand this off to LayoutData::apply.  This class
     then cycles through all the T objects and invokes this function.  Function
     must not be inline. (I'm still trying to figure out a nice way to send
     in non-static member functions).
  */
  virtual void apply(void (*a_Function)(const Box& box, int comps, T& t));

  ///
  virtual bool isDefined() const;

  virtual void clear() ;

  static int s_verbosity;

protected:
  int             m_comps;
  bool            m_isdefined;

  friend class LevelData<T>;

  void setVector(const BoxLayoutData<T>& da,
                 const Interval& srcComps,
                 const Interval& destComps);

  void allocateGhostVector(const DataFactory<T>& factory,
                           const IntVect& ghost = IntVect::Zero);

  void makeItSo(const Interval&     a_srcComps,
                const BoxLayoutData<T>& a_src,
                BoxLayoutData<T>&   a_dest,
                const Interval&     a_destComps,
                const Copier&       a_copier,
                const LDOperator<T>& a_op = CopyLDOperator<T>()) const;
  void makeItSoBegin(
                     const Interval&     a_srcComps,
                     const BoxLayoutData<T>& a_src,
                     BoxLayoutData<T>&   a_dest,
                     const Interval&     a_destComps,
                     const Copier&       a_copier,
                     const LDOperator<T>& a_op = CopyLDOperator<T>()) const;
  void makeItSoLocalCopy(
                         const Interval&     a_srcComps,
                         const BoxLayoutData<T>& a_src,
                         BoxLayoutData<T>&   a_dest,
                         const Interval&     a_destComps,
                         const Copier&       a_copier,
                         const LDOperator<T>& a_op = CopyLDOperator<T>()) const;

  void makeItSoLocalCopyNoOMP(
                         const Interval&     a_srcComps,
                         const BoxLayoutData<T>& a_src,
                         BoxLayoutData<T>&   a_dest,
                         const Interval&     a_destComps,
                         const Copier&       a_copier,
                         const LDOperator<T>& a_op = CopyLDOperator<T>()) const;
  void makeItSoLocalCopyThread(
                         const Interval&     a_srcComps,
                         const BoxLayoutData<T>& a_src,
                         BoxLayoutData<T>&   a_dest,
                         const Interval&     a_destComps,
                         const Copier&       a_copier,
                         const LDOperator<T>& a_op = CopyLDOperator<T>()) const;
  void makeItSoEnd(
                   BoxLayoutData<T>&   a_dest,
                   const Interval&     a_destComps,
                   const LDOperator<T>& a_op = CopyLDOperator<T>()) const;

  //========================================================================
  //
  // data structures used by makeItSo when we have some
  // data that needs to be moved (ie. there are entries
  // in the 'FROM' or 'TO' CopyIterators)
  //
  void completePendingSends() const;

  void allocateBuffers(const BoxLayoutData<T>& a_src,
                       const Interval& a_srcComps,
                       const BoxLayoutData<T>& a_dest,
                       const Interval& a_destComps,
                       const Copier&   a_copier,
                       const LDOperator<T>& a_op) const;

  void writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                      const Interval& a_srcComps,
                                      const LDOperator<T>& a_op) const;

  void writeSendDataFromMeIntoBuffersThread(const BoxLayoutData<T>& a_src,
                                            const Interval& a_srcComps,
                                            const LDOperator<T>& a_op) const;

  void writeSendDataFromMeIntoBuffersNoOMP(const BoxLayoutData<T>& a_src,
                                           const Interval& a_srcComps,
                                           const LDOperator<T>& a_op) const;

  void postSendsFromMe() const ;

  void postReceivesToMe() const ;

  void unpackReceivesToMe(BoxLayoutData<T>& a_dest,
                          const Interval&   a_destComps,
                          const LDOperator<T>& a_op) const ;

  void unpackReceivesToMeThread(BoxLayoutData<T>&    a_dest,
                                const Interval&      a_destComps,
                                const LDOperator<T>& a_op) const ;

  void unpackReceivesToMeNoOMP(BoxLayoutData<T>&    a_dest,
                               const Interval&      a_destComps,
                               const LDOperator<T>& a_op) const ;

  void unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                 const Interval&   a_destComps,
                                 int ncomp,
                                 const DataFactory<T>& factory,
                                 const LDOperator<T>& a_op) const;


  void unpackReceivesToMe_appendThread(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                       const Interval&   a_destComps,
                                       int ncomp,
                                       const DataFactory<T>& factory,
                                       const LDOperator<T>& a_op) const;

  void unpackReceivesToMe_appendNoOMP (LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                       const Interval&   a_destComps,
                                       int ncomp,
                                       const DataFactory<T>& factory,
                                       const LDOperator<T>& a_op) const;

  bool  doingJustCopy(const LDOperator<T>& a_op) const;

  /** \name Parallel messaging members */
  /*@{*/
  //   mutable void*  m_sendbuffer; // pointer member OK here,
  //                                // since LevelData<T> has no copy
  //   mutable size_t m_sendcapacity;
  //   mutable void*  m_recbuffer;  // pointer member OK here,
  //                                // since LevelData<T> has no copy
  //   mutable size_t m_reccapacity;
  mutable CopierBuffer*  m_buff;
#ifdef CH_MPI

#ifndef DOXYGEN

#endif
  //   mutable std::vector<bufEntry> m_fromMe;
  //   mutable std::vector<bufEntry> m_toMe;

  mutable Vector<MPI_Request>  m_sendRequests,  m_receiveRequests;
  mutable Vector<MPI_Status>   m_receiveStatus, m_sendStatus;
  mutable int numSends, numReceives;
#endif

};

///
/** not actually L-p norms, since it doesn't take into
    account the dx of the system.  A user can take that
    into account or not.

    For p != 0, returns pth root of sum of pth powers over all
    points in all fabs and all components in the interval:

    ( sum [ |A[i][pt,var]|^p : FArrayBox A[i], point pt in A[i].box(), var in interval ] )^(1/p)

    To turn into an L-p norm, one needs to multiply this by dx^(SpaceDim/p).

    For p == 0, returns global max over all points in all fabs and all
    components in the interval:

    max [ |A[i][pt,var]| : FArrayBox A[i], point pt in A[i].box(), var in interval ]

    Some people don't like that this norm is not normalized based on
    number of points in A.  Normalization is your problem.
*/
Real norm(const BoxLayoutData<FArrayBox>& A,
          const Interval& interval,
          const int& p = 2);

//======================================================================
template < >
BaseFab<int>* DefaultDataFactory<BaseFab<int> >::create(const Box& box,
                                                        int ncomps,
                                                        const DataIndex& a_datInd) const;

template < >
FArrayBox* DefaultDataFactory<FArrayBox>::create(const Box& box,
                                                 int ncomps,
                                                 const DataIndex& a_datInd) const;

#include "NamespaceFooter.H"
#include "BoxLayoutDataI.H"

#endif //BOXLAYOUTDATA
