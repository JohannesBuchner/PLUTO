#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MINIIVFAB_H_
#define _MINIIVFAB_H_

#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "VolIndex.H"
#include "BaseFab.H"
#include "EBGraph.H"
#include "BaseIVFAB.H"
#include "NamespaceHeader.H"

///
/**
   MiniIVFAB is a templated
   data holder defined at the VoFs of an irregular domain.

Implemented as just a raw vector of vofs and data, more optimized
for smaller memory footprint and faster linearIn/linearOut.  will
be more brutal for vof-by-vof indexing.
bvs
*/
template <class T>
class MiniIVFAB : public BaseIVFAB<T>
{
public:
  ///
  /**
     Default constructor.  Constructs an uninitialized IVFAB.  You must
     subsequently call {\bf define} before using this IVFAB.
  */
  MiniIVFAB();

  ///
  /**
     Defining constructor.  Specifies the irregular domain
     and the number of data components per VoF.  The irregular domain
     must lie completely within the EBGraph.  The
     contents are uninitialized.  Calls full define function.
  */
  MiniIVFAB(const IntVectSet& a_region,
             const EBGraph&    a_ebgraph,
             const int&        a_nvarin);

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
  MiniIVFAB(const Interval& a_comps,
            BaseIVFAB<T>&   a_original);


  ///
  virtual ~MiniIVFAB();

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

  const Vector<VolIndex>&  getVoFs() const;

  ///invalid but necessary for leveldata to compile
  MiniIVFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for MiniIVFAB");
  }

  ///get index into vector
  /**
     needs to be public so that bulk stencils can be constructed
   */
  virtual T* getIndex(const VolIndex& a_vof,const int& a_comp) const;

private:
  virtual void
  setDefaultValues();

protected:

  Vector<VolIndex> m_vofs;

};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "MiniIVFABI.H"
#endif

#endif
