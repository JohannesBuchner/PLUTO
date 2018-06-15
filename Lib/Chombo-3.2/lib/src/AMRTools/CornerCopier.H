#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CORNERCOPIER_H_
#define _CORNERCOPIER_H_

#include "DisjointBoxLayout.H"
#include "Pool.H"
#include "Vector.H"
#include "ProblemDomain.H"
#include "Copier.H"
#include "NamespaceHeader.H"


/// An even stranger (than Copier) thing to copy from ghost cells to corner ghost cells.
/**
   This is a specialized Copier to fill corner ghost cells from ghost
   cells of adjacent boxes.  Calling exchange with this copier will
   fill _only_ the corner ghost cells.  All other exchange operations
   should be performed by calling a separate exchange with a normal
   Copier.
*/

class CornerCopier: public Copier
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  CornerCopier()
  {
    m_ghost = -IntVect::Unit;
  }

  /// Basic defining constructor
  /**  Note that this constructor doesn't take a ghost cell argument,
       so it creates a useless CornerCopier (no intersections are computed).
       Provided for compatibility with Copier base class.
  */
  CornerCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
               bool a_exchange = false);

  /// this constructor contains support for periodic BC's
  /**  Note that this constructor doesn't take a ghost cell argument,
       so it creates a useless CornerCopier (no intersections are computed).
       Provided for compatibility with Copier base class.
  */
  CornerCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               bool  a_exchange  = false);

  /// Copier to copy into corner ghost cells of a_dest from ghost cells of  a_level
  /** Assumption is that a_dest and a_level have the same number of
      ghost cells.
  */
  CornerCopier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const IntVect& a_ghost,
               bool  a_exchange = false);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level with support for periodic BCs
  /** Assumption is that a_dest and a_level have the same number of
      ghost cells.
  */
  CornerCopier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               const IntVect& a_ghost,
               bool  a_exchange = false);

  ///
  virtual ~CornerCopier();

  /// define to match constructor (also doesn't do anything)
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      bool  a_exchange = false);

  /// define to match constructor (also doesn't do anything)
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      bool  a_exchange = false);

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_ghost,
                      bool  a_exchange = false);

  /// contains support for periodic BCs
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_ghost,
                      bool  a_exchange = false);

  ///
  virtual void clear();

  const IntVect& ghost()
  {
    return m_ghost;
  }

protected:

private:

  IntVect m_ghost;
};

std::ostream& operator<<(std::ostream& os, const CornerCopier& copier);


#include "NamespaceFooter.H"
#endif
