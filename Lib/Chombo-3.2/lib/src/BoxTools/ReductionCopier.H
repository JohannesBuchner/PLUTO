#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REDUCTIONCOPIER_H_
#define _REDUCTIONCOPIER_H_

#include "DisjointBoxLayout.H"
#include "Pool.H"
#include "Vector.H"
#include "ProblemDomain.H"
#include "Copier.H"
#include "NamespaceHeader.H"

/// Specialized copier for doing reductions in dimensionality
/**
   This is a specialized Copier designed to copy from a multi-dimensional
   DisjointBoxLayout to a smaller-dimensional one by copying all of the
   data in the transverse direction to the destination boxLayout. It is
   assumed that this will be used with a different sort of operator (like
   a summation operator, for example), since a straight copy operation
   wouldn't make much sense.  Anyway, this just defines the intersection
   regions; what you do with them is your problem.

   Note also that both the Spreading and Reduction Copiers ignore periodicity,
   since it can lead to unintended side-effects and multiple-counting of data.
   The way to handle periodicity is to do an exchange on the source data before
   calling the copyTo which does the reduction.
*/

class ReductionCopier : public Copier
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  ReductionCopier()
  {}

  ///
  ReductionCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
                  int a_transverseDir, bool a_exchange = false);

  /// this constructor contains support for periodic BC's
  ReductionCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
                  const ProblemDomain& a_domain,
                  int a_transverseDir,
                  bool  a_exchange  = false);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level
  ReductionCopier(const DisjointBoxLayout& a_level,
                  const BoxLayout& a_dest,
                  const IntVect& a_destGhost,
                  int a_transverseDir,
                  bool  a_exchange = false);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level with support for periodic BCs
  ReductionCopier(const DisjointBoxLayout& a_level,
                  const BoxLayout& a_dest,
                  const ProblemDomain& a_domain,
                  const IntVect& a_destGhost,
                  int a_transverseDir,
                  bool  a_exchange = false);

  ///
  ReductionCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
                  const Vector<int>& a_transverseDir, bool a_exchange = false);

  /// this constructor contains support for periodic BC's
  ReductionCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
                  const ProblemDomain& a_domain,
                  const Vector<int>& a_transverseDir,
                  bool  a_exchange  = false);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level
  ReductionCopier(const DisjointBoxLayout& a_level,
                  const BoxLayout& a_dest,
                  const IntVect& a_destGhost,
                  const Vector<int>& a_transverseDir,
                  bool  a_exchange = false);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level with support for periodic BCs
  ReductionCopier(const DisjointBoxLayout& a_level,
                  const BoxLayout& a_dest,
                  const ProblemDomain& a_domain,
                  const IntVect& a_destGhost,
                  const Vector<int>& a_transverseDir,
                  bool  a_exchange = false);

  ///
  virtual ~ReductionCopier();

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      int a_transverseDir,
                      bool  a_exchange = false);

  /// contains support for periodic BCs
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      int a_transverseDir,
                      bool  a_exchange = false);

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_destGhost,
                      int a_transverseDir,
                      bool  a_exchange = false);

  /// contains support for periodic BCs
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_destGhost,
                      int a_transverseDir,
                      bool  a_exchange = false);

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const Vector<int>& a_transverseDir,
                      bool  a_exchange = false);

  /// contains support for periodic BCs
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const Vector<int>& a_transverseDir,
                      bool  a_exchange = false);

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_destGhost,
                      const Vector<int>& a_transverseDir,
                      bool  a_exchange = false);

  /// contains support for periodic BCs
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_destGhost,
                      const Vector<int>& a_transverseDir,
                      bool  a_exchange = false);

  /// reverse the direction of copy for an already-defined ReductionCopier
  /**
      This function simply replaces the from data with the to-data for
      the already-defined copy operation.
      For now, this function heads right into a MayDay::Error, since what
      you probably want is to use a SpreadingCopier rather than a reversed
      ReductionCopier
  */
  void reverse();

  /// check that this Copier object corresponds to these two BoxLayouts.
  bool check(const DisjointBoxLayout& from, const BoxLayout& to) const;

  int print() const;


protected:

  friend class CopyIterator;

  friend void dumpmemoryatexit();

private:

  // no really good reason to keep this, except for completeness...
  Vector<int> m_transverseDir;

};

//===========================================================================
// end of public interface for ReductionCopier.
//===========================================================================

#include "NamespaceFooter.H"

#endif
