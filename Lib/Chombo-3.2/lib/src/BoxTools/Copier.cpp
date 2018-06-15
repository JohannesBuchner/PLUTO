#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DataIterator.H"
#include "IntVect.H"
#include "Copier.H"
#include "MayDay.H"
#include "LayoutIterator.H"
#include "NeighborIterator.H"
#include "SPMD.H"
#include "CH_Timer.H"
#include "memtrack.H"
#include "parstream.H"

#include <vector>
#include "NamespaceHeader.H"

using std::ostream;

Pool Copier::s_motionItemPool(sizeof(MotionItem), "Copier::MotionItem");

CopierBuffer::~CopierBuffer()
{
  clear();
}

void CopierBuffer::clear()
{
  if (m_sendbuffer != NULL) freeMT(m_sendbuffer);
  if (m_recbuffer  != NULL) freeMT(m_recbuffer);
  m_sendbuffer = NULL;
  m_recbuffer  = NULL;
  m_sendcapacity = 0;
  m_reccapacity = 0;
  m_ncomps = 0;
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               bool a_exchange,
               IntVect a_shift)
{
  const ProblemDomain& domain = a_level.physDomain();
  define(a_level, a_dest, domain, IntVect::Zero, a_exchange, a_shift);
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               bool a_exchange,
               IntVect a_shift)
{
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange, a_shift);
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const IntVect& a_destGhost,
               bool a_exchange,
               IntVect a_shift)
{
  const ProblemDomain domain = a_level.physDomain();
  define(a_level, a_dest, domain, a_destGhost, a_exchange, a_shift);
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               const IntVect& a_destGhost,
               bool a_exchange,
               IntVect a_shift)
{
  define(a_level, a_dest, a_domain, a_destGhost, a_exchange, a_shift);
}

void Copier::clear()
{
  for (unsigned int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      s_motionItemPool.returnPtr( m_localMotionPlan[i] );
    }
  for (unsigned int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      s_motionItemPool.returnPtr( m_fromMotionPlan[i] );
    }
  for (unsigned int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      s_motionItemPool.returnPtr( m_toMotionPlan[i] );
    }

  m_localMotionPlan.resize(0);
  m_fromMotionPlan.resize(0);
  m_toMotionPlan.resize(0);
  m_isDefined = false;
  m_buffers.clear();
}

Copier::Copier(const Copier& a_rhs)
{
  *this = a_rhs;
}

Copier& Copier::operator= (const Copier& b)
{
  clear();

  m_localMotionPlan.resize(b.m_localMotionPlan.size());
  for (int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      m_localMotionPlan[i] = new (s_motionItemPool.getPtr()) MotionItem(*(b.m_localMotionPlan[i]));
    }

  m_fromMotionPlan.resize(b.m_fromMotionPlan.size());
  for (int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      m_fromMotionPlan[i] = new (s_motionItemPool.getPtr()) MotionItem(*(b.m_fromMotionPlan[i]));
    }

  m_toMotionPlan.resize(b.m_toMotionPlan.size());
  for (int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      m_toMotionPlan[i] = new (s_motionItemPool.getPtr()) MotionItem(*(b.m_toMotionPlan[i]));
    }

  m_isDefined = true;
  return *this;
}

void Copier::trimMotion(const DisjointBoxLayout& a_exchangedLayout, const IntVect& a_ghost,
                        const Vector<MotionItem*>& a_oldItems, Vector<MotionItem*>& a_newItems)
{
  for (int i = 0; i < a_oldItems.size(); ++i)
    {
      const MotionItem& item = *(a_oldItems[i]);
      const Box& b = a_exchangedLayout[item.toIndex];
      const Box& c = item.toRegion;
      for (int d = 0; d < CH_SPACEDIM; ++d)
        {
          Box g(b);
          g.grow(d, a_ghost[d]);
          if (g.intersectsNotEmpty(c))
            {
              MotionItem* it = new (s_motionItemPool.getPtr()) MotionItem(item);
              a_newItems.push_back(it);
              break;
            }
        }
    }
}

///
void Copier::trimEdges(const DisjointBoxLayout& a_exchangedLayout, const IntVect& a_ghost)
{
  Copier oldCopier = *this;
  clear();

  trimMotion(a_exchangedLayout, a_ghost, oldCopier.m_localMotionPlan, m_localMotionPlan);
  //   pout() << "old Copy operations:" << oldCopier.m_localMotionPlan.size() << "  "
  //         << "new Copy operations:" << m_localMotionPlan.size() << "\n";
  trimMotion(a_exchangedLayout, a_ghost, oldCopier.m_fromMotionPlan, m_fromMotionPlan);
  trimMotion(a_exchangedLayout, a_ghost, oldCopier.m_toMotionPlan, m_toMotionPlan);
}

void Copier::reverse()
{
  for (int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      m_localMotionPlan[i]->reverse();
    }
  for (int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      m_fromMotionPlan[i]->reverse();
    }
  for (int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      m_toMotionPlan[i]->reverse();
    }
  m_fromMotionPlan.swap(m_toMotionPlan);

}

void Copier::coarsen(int a_refRatio)
{
  for (int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      m_localMotionPlan[i]->fromRegion.coarsen(a_refRatio);
      m_localMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }
  for (int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      m_fromMotionPlan[i]->fromRegion.coarsen(a_refRatio);
      m_fromMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }
  for (int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      m_toMotionPlan[i]->fromRegion.coarsen(a_refRatio);
      m_toMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }

}
void Copier::define(const DisjointBoxLayout& a_level,
                    const BoxLayout& a_dest,
                    bool a_exchange,
                    IntVect a_shift)
{
  const ProblemDomain domain = a_level.physDomain();
  define(a_level, a_dest, domain, IntVect::Zero, a_exchange, a_shift);
}

void Copier::define(const DisjointBoxLayout& a_level,
                    const BoxLayout& a_dest,
                    const ProblemDomain& a_domain,
                    bool a_exchange,
                    IntVect a_shift)
{
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange, a_shift);
}

void Copier::define(const DisjointBoxLayout& a_level,
                    const BoxLayout& a_dest,
                    const IntVect& a_ghost,
                    bool a_exchange,
                    IntVect a_shift)
{
  const ProblemDomain domain = a_level.physDomain();
  define(a_level, a_dest, domain, a_ghost, a_exchange, a_shift);
}

void Copier::define(const BoxLayout& a_level,
                    const BoxLayout& a_dest,
                    const ProblemDomain& a_domain,
                    const IntVect& a_ghost,
                    bool a_exchange,
                    IntVect a_shift)
{
  CH_TIME("Copier::define");
  CH_assert(a_level.isClosed());
  CH_assert(a_dest.isClosed());
  //  CH_assert(a_level.checkPeriodic(a_domain));

  clear();
  m_isDefined = true;
  buffersAllocated = false;
  //bool self = a_dest == a_level;
  const BoxLayout& level = a_level;
  const BoxLayout& dest = a_dest;

  // set up vector of dataIndexes to keep track of which
  // "to" boxes are not completely contained within the primary
  // domain.  these boxes are then candidates for filling by
  // periodic images of the "from" data.
  Vector<DataIndex> periodicallyFilledToVect;

  // in order to cull which "from" data may be needed to
  // fill the "to" data, keep track of the radius around the
  // primary domain in which all these cells lie.
  // do this by incrementally growing the domain box and
  // keeping track of what this radius is.
  // just to make things simpler, start off with a radius of one
  Box grownDomainCheckBox = a_domain.domainBox();
  grownDomainCheckBox.grow(1);
  int periodicCheckRadius = 1;

  // since valid regions of the "from" DBL may also be outside
  // the primary domain, need to keep track of whether any of these
  // need to be checked separately.
  Vector<DataIndex> periodicFromVect;
  // use same domain trick here as well
  Box grownFromDomainCheckBox = a_domain.domainBox();
  int periodicFromCheckRadius = 1;

  Box domainBox(a_domain.domainBox());
  bool isPeriodic = false;
  if (!domainBox.isEmpty())
    isPeriodic = a_domain.isPeriodic();

  // (dfm -- 9/13/05) as currently written, the Copier won't correctly
  // handle periodic cases where the number of ghost cells is greater
  // than the width of the domain.  We _should_ do multiple wraparounds,
  // but we don't. So, put in this assertion. We can revisit this if it
  // becomes an issue
  if (isPeriodic)
    {
      for (int dir = 0; dir < SpaceDim; dir++)
        {
          if (a_domain.isPeriodic(dir))
            {
              CH_assert (a_ghost[dir] <= domainBox.size(dir));
            }
        }
    }

  unsigned int myprocID = procID();

  // The following 4 for loops are the result of a performance optimization.
  // When increasing the size of the problem, we found that the code was
  // looping over every destination box for every source box which was N1*N2
  // loop iterations (essentially an N-squared approach).
  // The following code attempts to simply reduce N1 and N2 by first separating
  // the boxes (or LayoutIndexes to boxes) that reside on the current processor.
  // Then the loop to determine which boxes of the first list intersect with
  // which boxes of the second list can be done in N1' * N2' iterations,
  // where N1' is the reduced N1 and N2' is the reduced N2.
  // We have to break up the assigning of MotionItems into two separate
  // loops and be careful about the local copies.  These 4 loops are
  // significantly faster than the original for loop -- _especially_
  // for large problems.  (ndk)

#ifdef CH_MPI  // don't need to do this in serial
  // make a vector of boxes (or LayoutIndexes to boxes) from destination layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorDestDI;
  vector<DataIndex> vectorDestOnProcDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
  {
    vectorDestDI.push_back(DataIndex(to()));
    if (myprocID == dest.procID(to()))
    {
      vectorDestOnProcDI.push_back(DataIndex(to()));
    }
  }

  // make a vector of boxes (or LayoutIndexes to boxes) from "level"/src layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorLevelDI;
  vector<DataIndex> vectorLevelOnProcDI;
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
  {
    vectorLevelDI.push_back(DataIndex(from()));
    if (myprocID == level.procID(from()))
    {
      vectorLevelOnProcDI.push_back(DataIndex(from()));
    }
  }
#else
  // in serial, it's not very interesting as it's all of them.
  vector<DataIndex> vectorDestOnProcDI;
  vector<DataIndex> vectorLevelDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
  {
    vectorDestOnProcDI.push_back(DataIndex(to()));
  }
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
  {
    vectorLevelDI.push_back(DataIndex(from()));
  }
#endif

  bool isSorted = (a_level.isSorted() && a_dest.isSorted());

  // loop over all dest/to DI's on my processor
  for (vector<DataIndex>::iterator vdi = vectorDestOnProcDI.begin();
      vdi != vectorDestOnProcDI.end(); ++vdi)
  {
    // at this point, i know myprocID == toProcID
    const DataIndex todi(*vdi);

    Box ghost(dest[todi]);
    ghost -= a_shift;

    ghost.grow(a_ghost);

    //bool isSorted = (a_level.isSorted() && a_dest.isSorted());
    // then for each level/from DI, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelDI.begin();
        vli != vectorLevelDI.end(); ++vli)
    {
      const DataIndex fromdi(*vli);
      const unsigned int fromProcID = level.procID(fromdi);
      const Box& fromBox = level[fromdi];
      if ((fromBox.bigEnd(0) < ghost.smallEnd(0)) && isSorted )
      {
        //can skip rest cuz we haven't gotten to something interesting
        continue;
      }

      if (ghost.intersectsNotEmpty(fromBox))
      {
        Box srcBox(ghost); // ??
        srcBox &= fromBox; // ??

        Box destBox = srcBox + a_shift;

        MotionItem* item = new (s_motionItemPool.getPtr())
          MotionItem(fromdi, todi, srcBox, destBox);
        if (item == NULL)
        {
          MayDay::Error("Out of Memory in copier::define");
        }
        if (fromProcID == myprocID)
        { // local move
          if (a_exchange && fromdi == todi)
            s_motionItemPool.returnPtr(item);
          else
            m_localMotionPlan.push_back(item);
        }
        else
        {
          item->procID = fromProcID;
          m_toMotionPlan.push_back(item);
        }
      }
      if ((fromBox.smallEnd(0) > ghost.bigEnd(0)) && isSorted)
      {
        //can break out of loop, since we know that the smallEnd
        // of all the remaining boxes are lexigraphically beyond this ghosted box.
        break;
      }

    }
  }

  // Don't need to worry about this in serial as we already
  // took care of the local copy motion items just above.  skip this.
#ifdef CH_MPI
  // loop over all dest/to DI's
  for (vector<DataIndex>::iterator vdi = vectorDestDI.begin();
      vdi != vectorDestDI.end(); ++vdi)
  {

    const DataIndex todi(*vdi);

    Box ghost(dest[todi]);
    ghost -= a_shift;

    ghost.grow(a_ghost);

    const unsigned int toProcID = dest.procID(todi);

    // then for each level/from DI on this processor, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelOnProcDI.begin();
        vli != vectorLevelOnProcDI.end(); ++vli)
    {

      // at this point, i know myprocID == fromProcID

      const DataIndex fromdi(*vli);
      const Box& fromBox = level[fromdi];

      if ((fromBox.bigEnd(0) < ghost.smallEnd(0)) && isSorted)
      {
        //can skip rest cuz we haven't gotten to something interesting
        continue;
      }

      if (ghost.intersectsNotEmpty(fromBox))
      {
        Box srcBox(ghost); // ??
        srcBox &= fromBox; // ??

        Box destBox = srcBox + a_shift;

        if (toProcID == myprocID)
        { // local move
          // don't push back here!  or you will get two.
          //     we already did it above...
          //m_localMotionPlan.push_back(item);
        }
        else
        {
          MotionItem* item = new (s_motionItemPool.getPtr())
            MotionItem(fromdi, todi, srcBox, destBox);
          if (item == NULL)
          {
            MayDay::Error("Out of Memory in copier::define");
          }

          item->procID = toProcID;
          m_fromMotionPlan.push_back(item);
        }
      }
      if ((fromBox.smallEnd(0) > ghost.bigEnd(0)) && isSorted)
      {
        //can break out of loop, since we know that the smallEnd
        // of all the remaining boxes are lexigraphically beyond this ghosted box.
        break;
      }
    }
  }
#endif

  if (isPeriodic && a_shift != IntVect::Zero)
  {
    MayDay::Error("Copier::define - domain periodic and a non-zero shift in the copy is not implemented");
  }

  // put periodic intersection checking in here for "to" boxes
  if (isPeriodic)
  {
    for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
    {
      Box ghost(dest[to()]);
      ghost.grow(a_ghost);
      //unsigned int toProcID = dest.procID(to());  // unused variable

      // only do this if ghost box hangs over domain edge
      if (!domainBox.contains(ghost))
      {
        // add the dataIndex for this box to the list
        // of boxes which we need to come back to
        periodicallyFilledToVect.push_back(DataIndex(to()));
        // now check to see if we need to grow the
        // periodic check radius
        if (!grownDomainCheckBox.contains(ghost))
        {
          // grow the domainCheckBox until it contains ghost
          while (!grownDomainCheckBox.contains(ghost))
          {
            grownDomainCheckBox.grow(1);
            periodicCheckRadius++;
          }
        } // end if we need to grow radius around domain

      } //end if ghost box is not contained in domain
    } // end if periodic
  }

  // Here ends the so-called N-squared optimizations.  the rest is unchanged. (ndk)

  // now do periodic checking, if necessary
  if (isPeriodic)
    {

      // the only "from" boxes we will need to check
      // will be those within periodicCheckRadius of the
      // domain boundary. so, create a box to screen out
      // those which we will need to check.
      Box shrunkDomainBox = a_domain.domainBox();
      shrunkDomainBox.grow(-periodicCheckRadius);

      ShiftIterator shiftIt = a_domain.shiftIterator();
      IntVect shiftMult(domainBox.size());

      // now loop over "from" boxes
      for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
        {
          // first check to see whether we need to look at this box
          const Box& fromBox = level[from()];

          if (!shrunkDomainBox.contains(fromBox))
            {
              unsigned int fromProcID = level.procID(from());

              // check to see if fromBox is contained in domain,
              // if not, add it to the list of fromBoxes we need to
              // go back and check separately to see if it will
              // fill one of the "to" boxes
              if (!domainBox.contains(fromBox))
                {
                  periodicFromVect.push_back(DataIndex(from()));

                  if (!grownFromDomainCheckBox.contains(fromBox))
                    {
                      while (!grownFromDomainCheckBox.contains(fromBox))
                        {
                          grownFromDomainCheckBox.grow(1);
                          periodicFromCheckRadius++;
                        }
                    } // end if we need to grow domain check box
                } // end if fromBox is outside domain

              // now loop over those "to" boxes which were not contained
              // in the domain
              for (int toRef = 0; toRef < periodicallyFilledToVect.size(); toRef++)
                {
                  DataIndex toIndex = periodicallyFilledToVect[toRef];
                  unsigned int toProcID = dest.procID(toIndex);

                  // don't worry about anything that doesn't involve this proc
                  if (toProcID != myprocID && fromProcID != myprocID)
                    {
                      // do nothing
                    }
                  else
                    {
                      Box ghost(dest[toIndex]);
                      ghost.grow(a_ghost);
                      // now need to loop over shift vectors and look at images
                      for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                        {
                          IntVect shiftVect(shiftIt()*shiftMult);
                          ghost.shift(shiftVect);
                          if (ghost.intersectsNotEmpty(fromBox)) // rarely happens
                            {
                              Box intersectBox(ghost);
                              intersectBox &= fromBox;
                              Box toBox(intersectBox);
                              toBox.shift(-shiftVect);
                              MotionItem* item = new (s_motionItemPool.getPtr())
                                MotionItem(DataIndex(from()), DataIndex(toIndex),
                                           intersectBox, toBox);
                              if (item == NULL)
                                {
                                  MayDay::Error("Out of Memory in copier::define");
                                }
                              if (toProcID == fromProcID) // local move
                                m_localMotionPlan.push_back(item);
                              else if (fromProcID == myprocID)
                                {
                                  item->procID = toProcID;
                                  m_fromMotionPlan.push_back(item);
                                }
                              else
                                {
                                  item->procID = fromProcID;
                                  m_toMotionPlan.push_back(item);
                                }

                            } // end if shifted box intersects

                          ghost.shift(-shiftVect);
                        } // end loop over shift vectors
                    } // end if either from box or to box are on this proc
                } // end loop over destination boxes
            } // end if source box is close to domain boundary
        } // end loop over destination boxes

      // now go back through the "from" boxes which were outside
      // the domain and see if they intersect any toBoxes
      if (periodicFromVect.size() != 0)
        {
          // the only "to" boxes we will need to check
          // will be those within periodicCheckRadius of the
          // domain boundary. so, create a box to screen out
          // those which we will need to check.
          shrunkDomainBox = a_domain.domainBox();
          shrunkDomainBox.grow(-periodicFromCheckRadius);

          // now loop over the "to" boxes
          for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
            {
              // first check to see whether we need to look at this box
              Box ghost(dest[to()]);
              ghost.grow(a_ghost);

              if (!shrunkDomainBox.contains(ghost))
                {
                  unsigned int toProcID = a_dest.procID(to());

                  // now loop over those "from" boxes which are not
                  // contained by the domain
                  for (int fromRef = 0; fromRef < periodicFromVect.size(); fromRef++)
                    {
                      DataIndex fromIndex = periodicFromVect[fromRef];
                      const Box& fromBox = level[fromIndex];
                      unsigned int fromProcID = level.procID(fromIndex);

                      // don't worry about anything which doesn't involve
                      // this proc
                      if (toProcID != myprocID && fromProcID != myprocID)
                        {
                          // do nothing
                        }
                      else
                        {
                          // now need to loop over shift vectors and look at images
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                            {
                              IntVect shiftVect(shiftIt()*shiftMult);
                              ghost.shift(shiftVect);
                              if (ghost.intersectsNotEmpty(fromBox))
                                {
                                  Box intersectBox(ghost);
                                  intersectBox &= fromBox;
                                  Box toBox(intersectBox);
                                  toBox.shift(-shiftVect);
                                  MotionItem* item = new (s_motionItemPool.getPtr())
                                    MotionItem(DataIndex(fromIndex), DataIndex(to()),
                                               intersectBox, toBox);
                                  if (item == NULL)
                                    {
                                      MayDay::Error("Out of Memory in copier::define");
                                    }
                                  if (toProcID == fromProcID) // local move
                                    m_localMotionPlan.push_back(item);
                                  else if (fromProcID == myprocID)
                                    {
                                      item->procID = toProcID;
                                      m_fromMotionPlan.push_back(item);
                                    }
                                  else
                                    {
                                      item->procID = fromProcID;
                                      m_toMotionPlan.push_back(item);
                                    }

                                } // end if shifted box intersects

                              ghost.shift(-shiftVect);
                            } // end loop over shift vectors
                        } // end if either from box or to box are on this proc
                    } // end loop over "from" boxes
                } // end if destination box is close to domain boundary
            } // end loop over destination boxes
        } // end if any of the "From" boxes were outside the domain

    } // end if we need to do anything for periodicity
  sort();
}

void Copier::ghostDefine(const DisjointBoxLayout& a_src,
                         const DisjointBoxLayout& a_dest,
                         const ProblemDomain& a_domain,
                         const IntVect& a_srcGhost)
{
  //first, define a regular copier operation
  define(a_dest, a_src,  a_domain, a_srcGhost);

  //now, reverse the direction of the operation.
  reverse();

  sort();
}

void Copier::ghostDefine(const DisjointBoxLayout& a_src,
                         const DisjointBoxLayout& a_dest,
                         const ProblemDomain& a_domain,
                         const IntVect& a_srcGhost,
                         const IntVect& a_destGhost)
{
  CH_TIME("Copier::ghostDefine");
  clear();
  BoxLayout src ; src.deepCopy(a_src);
  src.grow(a_srcGhost);
  src.close();
  define(src, a_dest, a_domain, a_destGhost);
}

void Copier::exchangeDefine(const DisjointBoxLayout& a_grids,
                            const IntVect& a_ghost)
{
  CH_TIME("Copier::exchangeDefine");
  clear();
  DataIterator dit = a_grids.dataIterator();
  NeighborIterator nit(a_grids);
  int myprocID = procID();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = a_grids[dit];
      Box bghost(b);
      bghost.grow(a_ghost);
      for (nit.begin(dit()); nit.ok(); ++nit)
        {
          Box neighbor = nit.box();
          int fromProcID = a_grids.procID(nit());
          if (neighbor.intersectsNotEmpty(bghost))
            {
              Box box(neighbor & bghost);

              MotionItem* item = new (s_motionItemPool.getPtr())
                MotionItem(DataIndex(nit()), dit(), nit.unshift(box), box);
              if (fromProcID == myprocID)
              { // local move
                m_localMotionPlan.push_back(item);
              }
              else
              {
                item->procID = fromProcID;
                m_toMotionPlan.push_back(item);
              }
            }
          neighbor.grow(a_ghost);
          if (neighbor.intersectsNotEmpty(b) && fromProcID != myprocID)
            {
              Box box(neighbor & b);
              MotionItem* item = new (s_motionItemPool.getPtr())
                MotionItem(dit(), DataIndex(nit()), box, nit.unshift(box) );
              item->procID = fromProcID;
              m_fromMotionPlan.push_back(item);
            }
        }

    }
  sort();
}

class MotionItemSorter
{
public:
  inline bool operator()(MotionItem* const & lhs, MotionItem* const& rhs) const;
};

inline bool  MotionItemSorter::operator()(MotionItem* const & lhs, MotionItem* const& rhs) const
{
  if (lhs->procID == rhs->procID)
    {
      const Box& left = lhs->toRegion;
      const Box& right = rhs->toRegion;
      if (left.smallEnd() == right.smallEnd())
        {
          return left.bigEnd().lexLT(right.bigEnd());
        }
      else
        {
          return left < right;
        }
    }
  return lhs->procID < rhs->procID;
}

void Copier::sort()
{
  //return;
  std::vector<MotionItem*>& vlocal = m_localMotionPlan.stdVector();
  std::sort(vlocal.begin(), vlocal.end(), MotionItemSorter());
  std::vector<MotionItem*>& vfrom  = m_fromMotionPlan.stdVector();
  std::sort(vfrom.begin(), vfrom.end(), MotionItemSorter());
  std::vector<MotionItem*>& vto = m_toMotionPlan.stdVector();
  std::sort(vto.begin(), vto.end(), MotionItemSorter());
}

int Copier::print() const
{
  std::cout << *this;
  return 0;
}

int Copier::numLocalCellsToCopy() const
{
  int sum = 0;
  for (unsigned int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      sum += m_localMotionPlan[i]->fromRegion.numPts();
    }
  return sum;
}

int Copier::numFromCellsToCopy() const
{
  int sum = 0;
  for (unsigned int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      sum += m_fromMotionPlan[i]->fromRegion.numPts();
    }
  return sum;
}

int Copier::numToCellsToCopy() const
{
  int sum = 0;
  for (unsigned int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      sum += m_toMotionPlan[i]->fromRegion.numPts();
    }
  return sum;
}

ostream& operator<< (ostream& os, const Copier& copier)
{
  os << "local(" << procID() << "): ";
  for (CopyIterator it(copier, CopyIterator::LOCAL); it.ok(); ++it)
    {
      os << " from" << it().fromRegion << " to " << it().toRegion << '\n'
         << "          ";
    }
  os << "\nfrom(" << procID() << "): ";
  for (CopyIterator it(copier, CopyIterator::FROM); it.ok(); ++it)
    {
      os << " from " << it().fromRegion << " to " << it().toRegion << "[" << it().procID << "]" << '\n'
         << "         ";
    }
  os << "\nto(" << procID() << "): ";
  for (CopyIterator it(copier, CopyIterator::TO); it.ok(); ++it)
    {
      os << " from " << it().fromRegion << " to " << it().toRegion << "[" <<it().procID << "]" << '\n'
         << "       ";
    }
  os << "\n";
  return os;
}

Copier::~Copier()
{
  CH_TIME("~Copier");
  clear();
}

void Copier::setBufferAllocated(bool arg) const
{
  buffersAllocated  = arg;
}

bool Copier::bufferAllocated() const
{
  return buffersAllocated;
}

const ProblemDomain&
Copier::getPhysDomain(const DisjointBoxLayout& a_level) const
{
  return a_level.physDomain();
}

#include "NamespaceFooter.H"
