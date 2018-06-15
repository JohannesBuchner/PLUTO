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
#include "SpreadingCopier.H"
#include "MayDay.H"
#include "LayoutIterator.H"
#include "SPMD.H"
#include "CH_Timer.H"
#include "parstream.H"

#include <vector>
using std::ostream;

#include "NamespaceHeader.H"

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 int a_transverseDir,
                                 bool a_exchange)
{
  // oops -- can't use this since DisjointBoxLayout::physDomain is protected
  //const ProblemDomain& domain = a_level.physDomain();
  // disable this (hopefully deprecated) function
  MayDay::Error("SpreadingCopier constructor without problemDomain is deprecated");

  ProblemDomain domain;
  define(a_level, a_dest, domain, IntVect::Zero,
         a_transverseDir, a_exchange);
}

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 const ProblemDomain& a_domain,
                                 int a_transverseDir,
                                 bool a_exchange)
{
  define(a_level, a_dest, a_domain, IntVect::Zero,
         a_transverseDir, a_exchange);
}

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 const IntVect& a_destGhost,
                                 int a_transverseDir,
                                 bool a_exchange)
{

  const ProblemDomain& domain = getPhysDomain(a_level);
  define(a_level, a_dest, domain, a_destGhost, a_transverseDir,
         a_exchange);
}

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 const ProblemDomain& a_domain,
                                 const IntVect& a_destGhost,
                                 int a_transverseDir,
                                 bool a_exchange)
{
  define(a_level, a_dest, a_domain, a_destGhost, a_transverseDir,
         a_exchange);
}

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 const Vector<int>& a_transverseDir,
                                 bool a_exchange)
{
  // oops -- can't use this since DisjointBoxLayout::physDomain is protected
  //const ProblemDomain& domain = a_level.physDomain();
  // disable this (hopefully deprecated) function
  MayDay::Error("SpreadingCopier constructor without problemDomain is deprecated");

  ProblemDomain domain;
  define(a_level, a_dest, domain, IntVect::Zero,
         a_transverseDir, a_exchange);
}

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 const ProblemDomain& a_domain,
                                 const Vector<int>& a_transverseDir,
                                 bool a_exchange)
{
  define(a_level, a_dest, a_domain, IntVect::Zero,
         a_transverseDir, a_exchange);
}

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 const IntVect& a_destGhost,
                                 const Vector<int>& a_transverseDir,
                                 bool a_exchange)
{

  const ProblemDomain& domain = getPhysDomain(a_level);
  define(a_level, a_dest, domain, a_destGhost, a_transverseDir,
         a_exchange);
}

SpreadingCopier::SpreadingCopier(const DisjointBoxLayout& a_level,
                                 const BoxLayout& a_dest,
                                 const ProblemDomain& a_domain,
                                 const IntVect& a_destGhost,
                                 const Vector<int>& a_transverseDir,
                                 bool a_exchange)
{
  define(a_level, a_dest, a_domain, a_destGhost, a_transverseDir,
         a_exchange);
}

void SpreadingCopier::define(const DisjointBoxLayout& a_level,
                             const BoxLayout& a_dest,
                             int a_transverseDir,
                             bool a_exchange)
{
  Vector<int> dirVect(1,a_transverseDir);
  define(a_level, a_dest, dirVect, a_exchange);
}

void SpreadingCopier::define(const DisjointBoxLayout& a_level,
                             const BoxLayout& a_dest,
                             const ProblemDomain& a_domain,
                             int a_transverseDir,
                             bool a_exchange)
{
  Vector<int> dirVect(1, a_transverseDir);
  define(a_level, a_dest, a_domain, IntVect::Zero, dirVect,
         a_exchange);
}

void SpreadingCopier::define(const DisjointBoxLayout& a_level,
                             const BoxLayout& a_dest,
                             const IntVect& a_ghost,
                             int a_transverseDir,
                             bool a_exchange)
{
  Vector<int> dirVect(1,a_transverseDir);
  const ProblemDomain& domain = getPhysDomain(a_level);
  define (a_level, a_dest, domain, a_ghost, dirVect,
          a_exchange);
}

void SpreadingCopier::define(const BoxLayout& a_level,
                             const BoxLayout& a_dest,
                             const ProblemDomain& a_domain,
                             const IntVect& a_ghost,
                             int a_transverseDir,
                             bool  a_exchange)
{
  Vector<int> dirVect(1,a_transverseDir);
  define(a_level, a_dest,a_domain, a_ghost, dirVect, a_exchange);
}

void SpreadingCopier::define(const DisjointBoxLayout& a_level,
                             const BoxLayout& a_dest,
                             const Vector<int>& a_transverseDir,
                             bool a_exchange)
{
  const ProblemDomain& domain = getPhysDomain(a_level);
  define(a_level, a_dest, domain, IntVect::Zero, a_transverseDir,
         a_exchange);
}

void SpreadingCopier::define(const DisjointBoxLayout& a_level,
                             const BoxLayout& a_dest,
                             const ProblemDomain& a_domain,
                             const Vector<int>& a_transverseDir,
                             bool a_exchange)
{
  define(a_level, a_dest, a_domain, IntVect::Zero, a_transverseDir,
         a_exchange);
}

void SpreadingCopier::define(const DisjointBoxLayout& a_level,
                             const BoxLayout& a_dest,
                             const IntVect& a_ghost,
                             const Vector<int>& a_transverseDir,
                             bool a_exchange)
{
  const ProblemDomain& domain = getPhysDomain(a_level);
  define (a_level, a_dest, domain, a_ghost, a_transverseDir,
          a_exchange);
}

void SpreadingCopier::define(const BoxLayout& a_level,
                             const BoxLayout& a_dest,
                             const ProblemDomain& a_domain,
                             const IntVect& a_ghost,
                             const Vector<int>& a_transverseDir,
                             bool  a_exchange)
{
#ifdef MULTIDIM_TIMER
  CH_TIME("SpreadingCopier::define")
#endif
  m_isDefined = true;
  m_transverseDir = a_transverseDir;

  CH_assert(a_level.isClosed());
  CH_assert(a_dest.isClosed());
  //  CH_assert(a_level.checkPeriodic(a_domain));

  clear();
  buffersAllocated = false;
  //bool self = a_dest == a_level;
  const BoxLayout&         level= a_level;
  const BoxLayout&         dest = a_dest;

  // note that while much of the implementation for the
  // periodic case remains in this function (as copied, more
  // or less, from the original Copier), at least at the moment
  // periodic support is "turned off" for the SpreadingCopier
  // because it's not entirely clear it will do the right thing
  // for this case.  I've left the vast majority of the
  // periodic implementation in place, however, in case we
  // change our minds... (DFM 3/10/09)

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
  // grab index extents of domain in transverse direction
  Vector<int> transverseLo(m_transverseDir.size());
  Vector<int> transverseHi(m_transverseDir.size());

  for (int n=0; n<m_transverseDir.size(); n++)
    {
      // grab index extents of domain in transverse direction
      transverseLo[n] = domainBox.smallEnd(m_transverseDir[n]);
      transverseHi[n] = domainBox.bigEnd(m_transverseDir[n]);
      // allow for the possibility of transverse-direction
      // ghost cells outside of the domain
      transverseLo[n] -= a_ghost[m_transverseDir[n] ];
      transverseHi[n] += a_ghost[m_transverseDir[n] ];
    }

  bool isPeriodic = false;
  // I think the right thing here is that the SpreadingCopier
  // completely ignores periodicity
  //if (!domainBox.isEmpty())
  //isPeriodic = a_domain.isPeriodic();

  // (dfm -- 9/13/05) as currently written, the Copier won't correctly
  // handle periodic cases where the number of ghost cells is greater
  // than the width of the domain.  We _should_ do multiple wraparounds,
  // but we don't. So, put in this assertion. We can revisit this if it
  // becomes an issue
  if (isPeriodic)
    {
      for (int dir=0; dir<SpaceDim; dir++)
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

  // loop over all dest/to DI's on my processor
  for (vector<DataIndex>::iterator vdi=vectorDestOnProcDI.begin();
      vdi != vectorDestOnProcDI.end(); ++vdi)
  {

    // at this point, i know myprocID == toProcID
    const DataIndex todi(*vdi);
    Box ghost(dest[todi]);
    // handle destination ghost cells normally
    ghost.grow(a_ghost);

    //pout() << "GHOST = " << ghost << endl << endl;
    
    // don't do anything unless domain intersects ghost
    if (a_domain.intersects(ghost))
      {
        
        // then for each level/from DI, see if they intersect
        for (vector<DataIndex>::iterator vli = vectorLevelDI.begin();
             vli != vectorLevelDI.end(); ++vli)
          {
            
            const DataIndex fromdi(*vli);
            const unsigned int fromProcID = level.procID(fromdi);
            Box fromBox(level[fromdi]);
            //pout () << "    fromBox = " << fromBox << endl << endl;
            
            //  grow the fromBox if it's near the domain boundary and
            // if there are ghost cells specified.
            for (int dir=0; dir<SpaceDim; dir++)
              {
                bool isTransverseDir = false;
                for (int n=0; n<m_transverseDir.size(); n++)
                  {
                    if (dir == m_transverseDir[n]) isTransverseDir = true;
                  }
                // only need to do this if we specify a ghost radius
                // and if it's not a transverse direction
                if (!isTransverseDir && a_ghost[dir] != 0)
                  {
                    if (fromBox.bigEnd(dir) == domainBox.bigEnd(dir))
                      {
                        fromBox.growHi(dir, a_ghost[dir]);
                      }
                    
                    if (fromBox.smallEnd(dir) == domainBox.smallEnd(dir))
                      {
                        fromBox.growLo(dir, a_ghost[dir]);
                      }
                  }
              }
            
            Box grownFromBox(fromBox);
            // for this Copier, we want to grow the "src" region to cover the
            // entire domain in the transverse direction. In that way, we
            // will be able to catch all of the source data regions which
            // project onto the dest boxes
            // do this for all transverse directions
            for (int n=0; n<m_transverseDir.size(); n++)
              {
                grownFromBox.setSmall(m_transverseDir[n], transverseLo[n]);
                grownFromBox.setBig(m_transverseDir[n], transverseHi[n]);
              }
            
            if (grownFromBox.bigEnd(0) < ghost.smallEnd(0))
              {
                //can skip rest cuz we haven't gotten to something interesting
                continue;
              }
            
            if (ghost.intersectsNotEmpty(grownFromBox))
              {
                Box toBox(ghost); // ??
                toBox&=grownFromBox;   // ??
                
                Box grownToBox(toBox);
                // to get correct fromBox, grow toBox in the transverse direction(s)
                // then intersect with the fromBox to get the appropriate part
                // of the frombox
                
                // do this for all transverse directions
                for (int n=0; n<m_transverseDir.size(); n++)
                  {
                    grownToBox.setSmall(m_transverseDir[n], transverseLo[n]);
                    grownToBox.setBig(m_transverseDir[n], transverseHi[n]);
                  }
                fromBox &= grownToBox;
                
                //pout() << "motionitem: from = " << fromBox
                //       << ", to = " << toBox << endl;
                MotionItem* item = new (s_motionItemPool.getPtr())
                  MotionItem(fromdi, todi, fromBox, toBox);
                if (item == NULL)
                  {
                    MayDay::Error("Out of Memory in copier::define");
                  }
                if (fromProcID == myprocID)
                  {
                    // local move
                    if (a_exchange && fromdi == todi)
                      s_motionItemPool.returnPtr(item);
                    else
                      //pout() << "(local) " << endl;
                      m_localMotionPlan.push_back(item);
                  }
                else
                  {
                    item->procID = fromProcID;
                    //pout () << "(toItem)" << endl;
                    m_toMotionPlan.push_back(item);
                  }
              }
            if (fromBox.smallEnd(0) > ghost.bigEnd(0))
              {
                //can break out of loop, since we know that the smallEnd
                // of all the remaining boxes are lexigraphically beyond this ghosted box.
                break;
              }

          }
      } // end if original "ghost" intersects domain
  }
  
  // Don't need to worry about this in serial as we already
  // took care of the local copy motion items just above.  skip this.
#ifdef CH_MPI
  // loop over all dest/to DI's
  for (vector<DataIndex>::iterator vdi=vectorDestDI.begin();
      vdi != vectorDestDI.end(); ++vdi)
  {

    const DataIndex todi(*vdi);
    Box ghost(dest[todi]);
    // handle destination ghost cells normally
    ghost.grow(a_ghost);

    //pout() << "GHOST2 = " << ghost << endl << endl;
    
    // don't do anything unless domain intersects ghost
    if (a_domain.intersects(ghost))
      {
        const unsigned int toProcID = dest.procID(todi);
        
        // then for each level/from DI on this processor, see if they intersect
        for (vector<DataIndex>::iterator vli = vectorLevelOnProcDI.begin();
             vli != vectorLevelOnProcDI.end(); ++vli)
          {
            
            // at this point, i know myprocID == fromProcID
            
            const DataIndex fromdi(*vli);
            Box fromBox(level[fromdi]);
            
            //pout () << "    fromBox = " << fromBox << endl << endl;
            
            // grow the fromBox if it's near the domain boundary and
            // if there are ghost cells specified.
            for (int dir=0; dir<SpaceDim; dir++)
              {
                bool isTransverseDir = false;
                for (int n=0; n<m_transverseDir.size(); n++)
                  {
                    if (dir == m_transverseDir[n]) isTransverseDir = true;
                  }
                // only need to do this if we specify a ghost radius
                // and if it's not a transverse direction
                if (!isTransverseDir && a_ghost[dir] != 0)
                  {
                    if (fromBox.bigEnd(dir) == domainBox.bigEnd(dir))
                      {
                        fromBox.growHi(dir, a_ghost[dir]);
                      }
                    
                    if (fromBox.smallEnd(dir) == domainBox.smallEnd(dir))
                      {
                        fromBox.growLo(dir, a_ghost[dir]);
                      }
                  }
              }
            
            Box  grownFromBox(fromBox);
            // for this Copier, we want to grow the "src" region to cover the
            // entire domain in the transverse direction. In that way, we
            // will be able to catch all of the source data regions which
            // project onto the dest boxes
            // do this for all transverse directions
            for (int n=0; n<m_transverseDir.size(); n++)
              {
                grownFromBox.setSmall(m_transverseDir[n], transverseLo[n]);
                grownFromBox.setBig(m_transverseDir[n], transverseHi[n]);
              }
            
            if (grownFromBox.bigEnd(0) < ghost.smallEnd(0))
              {
                //can skip rest cuz we haven't gotten to something interesting
                continue;
              }
            
            if (ghost.intersectsNotEmpty(grownFromBox))
              {
                Box toBox(ghost); // ??
                toBox&=grownFromBox;   // ??
                
                Box grownToBox(toBox);
                // to get correct fromBox, grow toBox in the transverse direction(s)
                // then intersect with the fromBox to get the appropriate part
                // of the frombox
                
                // do this for all transverse directions
                for (int n=0; n<m_transverseDir.size(); n++)
                  {
                    grownToBox.setSmall(m_transverseDir[n], transverseLo[n]);
                    grownToBox.setBig(m_transverseDir[n], transverseHi[n]);
                  }
                fromBox &= grownToBox;
                
                if (toProcID == myprocID)
                  {
                    // local move
                    // don't push back here!  or you will get two.
                    //     we already did it above...
                    //m_localMotionPlan.push_back(item);
                  }
                else
                  {
                    //pout() << "motionitem2: from = " << fromBox
                    //       << ", to = " << toBox << endl;
                    MotionItem* item = new (s_motionItemPool.getPtr())
                      MotionItem(fromdi, todi, fromBox, toBox);
                    if (item == NULL)
                      {
                        MayDay::Error("Out of Memory in copier::define");
                      }
                    
                    //pout () << "(from)" << endl;
                    item->procID = toProcID;
                    m_fromMotionPlan.push_back(item);
                  }
              }
            if (fromBox.smallEnd(0) > ghost.bigEnd(0))
              {
                //can break out of loop, since we know that the smallEnd
                // of all the remaining boxes are lexigraphically beyond this ghosted box.
                break;
              }
          }
      } // end if original "ghost" intersects domain
  }
#endif

  // put periodic intersection checking in here for "to" boxes
  if (isPeriodic)
  {
    for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
    {

      Box ghost(dest[to()]);

    // this is somewhat different from the "normal" Copier. We only want
    // to look at ghost cells if they're outside the domain.
    // note that we don't want to do this for domain ghost cells in the
    // transverse (spreading) directions.
    for (int dir=0; dir<SpaceDim; dir++)
      {
        bool isTransverseDir = false;
        for (int n=0; n<m_transverseDir.size(); n++)
          {
            if (dir == m_transverseDir[n]) isTransverseDir = true;
          }
        // only need to do this if we specify a ghost radius
        if (!isTransverseDir && a_ghost[dir] != 0)
          {
            if (ghost.bigEnd(dir) == domainBox.bigEnd(dir))
              {
                ghost.growHi(dir, a_ghost[dir]);
              }

            if (ghost.smallEnd(dir) == domainBox.smallEnd(dir))
              {
                ghost.growLo(dir, a_ghost[dir]);
              }
          }
      }

    // don't do anything unless domain intersects ghost
    if (a_domain.intersects(ghost))
      {
        
        Box destBox(ghost);
        
        // for this Copier, we want to grow the "dest" region to cover the
        // entire domain in the transverse direction. In that way, we
        // will be able to catch all of the source data regions which
        // project onto the dest boxes
        // do this for each transverse directions
        for (int n=0; n<m_transverseDir.size(); n++)
          {
            ghost.setSmall(m_transverseDir[n], transverseLo[n]);
            ghost.setBig(m_transverseDir[n], transverseHi[n]);
          }
        
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
      }// end if original "ghost" intersects domain
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
              for (int toRef=0; toRef<periodicallyFilledToVect.size(); toRef++)
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
                      
                      // this is somewhat different from the "normal" Copier.
                      // We only want
                      // to look at ghost cells if they're outside the domain.
                      // note that we don't want to do this for domain ghost
                      // cells in the transverse (spreading) directions.
                      for (int dir=0; dir<SpaceDim; dir++)
                        {
                          bool isTransverseDir = false;
                          for (int n=0; n<m_transverseDir.size(); n++)
                            {
                              if (dir == m_transverseDir[n]) isTransverseDir = true;
                            }
                          // only need to do this if we specify a ghost radius
                          if (!isTransverseDir && a_ghost[dir] != 0)
                            {
                              if (ghost.bigEnd(dir) == domainBox.bigEnd(dir))
                                {
                                  ghost.growHi(dir, a_ghost[dir]);
                                }

                              if (ghost.smallEnd(dir) == domainBox.smallEnd(dir))
                                {
                                  ghost.growLo(dir, a_ghost[dir]);
                                }
                            }
                        }

                      // don't do anything unless domain intersects ghost
                      if (a_domain.intersects(ghost))
                        {
                          
                          Box destBox(ghost);
                          
                          // for this Copier, we want to grow the "dest"
                          // region to cover the entire domain in the transverse
                          // direction. In that way, we will be able to catch all
                          // of the source data regions which project onto the
                          // dest boxes
                          // do this for each transverse direction
                          for (int n=0; n<m_transverseDir.size(); n++)
                            {
                              ghost.setSmall(m_transverseDir[n], transverseLo[n]);
                              ghost.setBig(m_transverseDir[n], transverseHi[n]);
                            }
                          
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
                                  toBox &= destBox;
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
                        } // end if original "ghost" intersects domain
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

              // this is somewhat different from the "normal" Copier.
              // We only want
              // to look at ghost cells if they're outside the domain.
              // note that we don't want to do this for domain ghost cells
              // in the transverse (spreading) directions.
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  bool isTransverseDir = false;
                  for (int n=0; n<m_transverseDir.size(); n++)
                    {
                      if (dir == m_transverseDir[n]) isTransverseDir = true;
                    }
                  // only need to do this if we specify a ghost radius
                  if (!isTransverseDir && a_ghost[dir] != 0)
                    {
                      if (ghost.bigEnd(dir) == domainBox.bigEnd(dir))
                        {
                          ghost.growHi(dir, a_ghost[dir]);
                        }

                      if (ghost.smallEnd(dir) == domainBox.smallEnd(dir))
                        {
                          ghost.growLo(dir, a_ghost[dir]);
                        }
                    }
                }
              // don't do anything unless domain intersects ghost
              if (a_domain.intersects(ghost))
                {
                  
                  Box destBox(ghost);
                  
                  // for this Copier, we want to grow the "dest" region to
                  // cover the entire domain in the transverse direction.
                  // In that way, we will be able to catch all of the source
                  // data regions which project onto the dest boxes
                  // do this for each transverse direction
                  for (int n=0; n<m_transverseDir.size(); n++)
                    {
                      ghost.setSmall(m_transverseDir[n], transverseLo[n]);
                      ghost.setBig(m_transverseDir[n], transverseHi[n]);
                    }
                  
                  if (!shrunkDomainBox.contains(ghost))
                    {
                      unsigned int toProcID = a_dest.procID(to());
                      
                      // now loop over those "from" boxes which are not
                      // contained by the domain
                      for (int fromRef = 0; fromRef<periodicFromVect.size(); fromRef++)
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
                                      toBox &= destBox;
                                      
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
                } // end if original "ghost" intersects domain
            }  // end loop over destination boxes
        } // end if any of the "From" boxes were outside the domain
      
    } // end if we need to do anything for periodicity
}

void SpreadingCopier::reverse()
{
  // for now, make this a run-time error
  MayDay::Error("SpreadingCopier::reverse disabled -- try ReductionCopier");

  // reverse the direction of the operation.
  for (int i=0; i< m_localMotionPlan.size(); ++i)
    {
      m_localMotionPlan[i]->reverse();
    }
  for (int i=0; i< m_fromMotionPlan.size(); ++i)
    {
      m_fromMotionPlan[i]->reverse();
    }
  for (int i=0; i< m_toMotionPlan.size(); ++i)
    {
      m_toMotionPlan[i]->reverse();
    }
  m_fromMotionPlan.swap(m_toMotionPlan);
}

int SpreadingCopier::print() const
{
  pout() << *this;
  return 0;
}

SpreadingCopier::~SpreadingCopier()
{
  clear();
}

#include "NamespaceFooter.H"
