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
#include "CornerCopier.H"
#include "MayDay.H"
#include "LayoutIterator.H"
#include "SPMD.H"

#include <vector>
#include "NamespaceHeader.H"

using std::ostream;


CornerCopier::CornerCopier(const DisjointBoxLayout& a_level,
                           const BoxLayout& a_dest,
                           bool a_exchange)
{
  // don't call define here, since we know that it shouldn't do anything
  // since we don't have any ghost cells.
  MayDay::Warning("Called CornerCopier::define with no ghost cells!");
  return;

}

CornerCopier::CornerCopier(const DisjointBoxLayout& a_level,
                           const BoxLayout& a_dest,
                           const ProblemDomain& a_domain,
                           bool a_exchange)
{
  // call define here, even though we know that it shouldn't do anything
  // since we don't have any ghost cells.
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange);
}

CornerCopier::CornerCopier(const DisjointBoxLayout& a_level,
                           const BoxLayout& a_dest,
                           const IntVect& a_destGhost,
                           bool a_exchange)
{
  // can't get at the ProbDomain from here, at the moment,
  // so disallow this function
  MayDay::Error("Can't initialize CornerCopier w/o ProblemDomain arg");

}

CornerCopier::CornerCopier(const DisjointBoxLayout& a_level,
                           const BoxLayout& a_dest,
                           const ProblemDomain& a_domain,
                           const IntVect& a_destGhost,
                           bool a_exchange)
{
  define(a_level, a_dest, a_domain, a_destGhost, a_exchange);
}

void CornerCopier::clear()
{
  Copier::clear();
}

void CornerCopier::define(const DisjointBoxLayout& a_level,
                          const BoxLayout& a_dest,
                          bool a_exchange)
{
  // can't get at the ProbDomain from here, at the moment,
  // so disallow this function
  MayDay::Error("Can't initialize CornerCopier w/o ProblemDomain arg");
}

void CornerCopier::define(const DisjointBoxLayout& a_level,
                          const BoxLayout& a_dest,
                          const ProblemDomain& a_domain,
                          bool a_exchange)
{
  // call define, even though it won't do anything
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange);
}

void CornerCopier::define(const DisjointBoxLayout& a_level,
                          const BoxLayout& a_dest,
                          const IntVect& a_ghost,
                          bool a_exchange)
{
  // can't get at the ProbDomain from here, at the moment,
  // so disallow this function
  MayDay::Error("Can't initialize CornerCopier w/o ProblemDomain arg");
}

void CornerCopier::define(const BoxLayout& a_level,
                          const BoxLayout& a_dest,
                          const ProblemDomain& a_domain,
                          const IntVect& a_ghost,
                          bool  a_exchange)
{
  m_ghost = a_ghost;
  bool isSorted = (a_level.isSorted() && a_dest.isSorted());
  // if no ghost cells, warn and return
  if (a_ghost == IntVect::Zero)
  {
    MayDay::Warning("Called CornerCopier::define with no ghost cells!");
    return;
  }

  CH_assert(a_level.isClosed());
  CH_assert(a_dest.isClosed());
  //  CH_assert(a_level.checkPeriodic(a_domain));

  clear();
  m_isDefined = true;
  buffersAllocated = false;
  //bool self = a_dest == a_level;
  const BoxLayout&         level= a_level;
  const BoxLayout&         dest = a_dest;

  // set up vector of dataIndexes to keep track of which
  // "to" boxes are not completely contained within the primary
  // domain.  these boxes are then candidates for filling by
  // periodic images of the "from" data.
  Vector<DataIndex> periodicallyFilledToVect;
  // also keep track of which corner cells we're talking about
  // in that case
  Vector<Box> periodicallyFilledToBoxes;

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
  Vector<Box> periodicFromBoxes;
  // use same domain trick here as well
  Box grownFromDomainCheckBox = a_domain.domainBox();
  int periodicFromCheckRadius = 1;

  Box domainBox(a_domain.domainBox());
  bool isPeriodic = false;
  if (!domainBox.isEmpty())
    isPeriodic = a_domain.isPeriodic();

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

    // here's where things start to differ from the original copy
    // implementation.   Need to loop over SpaceDim, hiLo face, and then
    // hiLo corner in order to isolate corner cells.  This is 4*SpaceDim
    // times as many boxes as the original copier, so this could be
    // noticeably more expensive than a vanilla copier.

    // loop over directions
    for (int dir = 0; dir<SpaceDim; dir++)
    {
      Box cornerGhost;
      Box grownDestBox = dest[todi];
      grownDestBox.grow(a_ghost);
      // loop over hi-lo faces
      SideIterator faceSit;
      for (faceSit.begin(); faceSit.ok(); ++faceSit)
      {
        Box allGhostDir;
        Box tempBox;
        if (faceSit() == Side::Lo)
        {
          allGhostDir = grownDestBox;
          int chopPoint = dest[todi].smallEnd(dir);
          tempBox = allGhostDir.chop(dir, chopPoint);
        }
        else
        {
          tempBox = grownDestBox;
          int chopPoint = dest[todi].bigEnd(dir) +1;
          allGhostDir = tempBox.chop(dir, chopPoint);
        }
        // now we should have all of the ghost cells on
        // the hi or lo face in direction dir

        // now loop over directions for corners
        for (int cornerDir = 0; cornerDir<SpaceDim; cornerDir++)
        {
          if (cornerDir != dir)
          {
            // now loop over hi-lo corners
            SideIterator cornerSit;
            for (cornerSit.begin(); cornerSit.ok(); ++cornerSit)
            {
              // now isolate corner cells by chopping again
              if (cornerSit() == Side::Lo)
              {
                cornerGhost = allGhostDir;
                int chopPoint = dest[todi].smallEnd(cornerDir);
                tempBox = cornerGhost.chop(cornerDir, chopPoint);
              }
              else
              {
                tempBox = allGhostDir;
                int chopPoint = dest[todi].bigEnd(cornerDir)+1;
                cornerGhost = tempBox.chop(cornerDir, chopPoint);
              }

              // now we should have isolated the corner ghost cells
              // for this box in the given directions and faces.


              // then for each level/from DI, see if they intersect
              for (vector<DataIndex>::iterator vli = vectorLevelDI.begin();
                  vli != vectorLevelDI.end(); ++vli)
              {

                const DataIndex fromdi(*vli);
                const unsigned int fromProcID = level.procID(fromdi);

                Box testBox = level[fromdi];
                // if not periodic, skip intersections with myself
                if (isPeriodic || testBox != dest[todi])
                {
                  testBox.grow(a_ghost);
                  if ((testBox.bigEnd(0) < cornerGhost.smallEnd(0)) && isSorted)
                  {
                    // can skip rest cuz we haven't gotten
                    // to something interesting
                    continue;
                  }

                  if ((testBox.smallEnd(0) > cornerGhost.bigEnd(0)) && isSorted)
                  {
                    //can break out of loop, since we know that
                    // the smallEnd of all the remaining boxes
                    // are lexigraphically beyond this ghosted
                    //box.
                    break;
                  }


                  // only look at ghost cells with same faceDir
                  // and hi-lo as the dest

                  Box fromBox;
                  if (faceSit() == Side::Lo)
                  {
                    fromBox = adjCellLo(level[fromdi],
                                        dir,
                                        a_ghost[dir]);
                  }
                  else
                  {
                    fromBox = adjCellHi(level[fromdi],
                                        dir,
                                        a_ghost[dir]);
                  }


                  // now proceed as normal

                  if (cornerGhost.intersectsNotEmpty(fromBox))
                  {
                    Box box(cornerGhost); // ??
                    box&=fromBox;   // ??
                    MotionItem* item = new (s_motionItemPool.getPtr())
                      MotionItem(fromdi, todi, box);
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
                        m_localMotionPlan.push_back(item);
                    }
                    else
                    {
                      item->procID = fromProcID;
                      m_toMotionPlan.push_back(item);
                    }
                  }
                } // end if src and dest boxes not the same
              } // end loop over from boxes
            } // end loop over corner hi-lo
          } // end if corner dir != dir
        }  // end loop over corner dirs
      }  // end loop over hi-lo faces
    } // end loop over face directions
  } // end loop over dest boxes

  // Don't need to worry about this in serial as we already
  // took care of the local copy motion items just above.  skip this.
#ifdef CH_MPI
  // loop over all dest/to DI's
  for (vector<DataIndex>::iterator vdi=vectorDestDI.begin();
      vdi != vectorDestDI.end(); ++vdi)
  {

    const DataIndex todi(*vdi);
    const unsigned int toProcID = dest.procID(todi);

    // do this just like the on-proc ones above
    /// loop over directions
    for (int dir=0; dir<SpaceDim; dir++)
    {
      Box cornerGhost;
      Box grownDestBox = dest[todi];
      grownDestBox.grow(a_ghost);

      // loop over hi-lo faces
      SideIterator faceSit;
      for (faceSit.begin(); faceSit.ok(); ++faceSit)
      {
        Box allGhostDir;
        Box tempBox;
        if (faceSit() == Side::Lo)
        {
          allGhostDir = grownDestBox;
          int chopPoint = dest[todi].smallEnd(dir);
          tempBox = allGhostDir.chop(dir, chopPoint);
        }
        else
        {
          tempBox = grownDestBox;
          int chopPoint = dest[todi].bigEnd(dir) +1;
          allGhostDir = tempBox.chop(dir, chopPoint);
        }
        // now we should have all of the ghost cells on
        // the hi or lo face in direction dir

        // now loop over directions for corners
        for (int cornerDir = 0; cornerDir<SpaceDim; cornerDir++)
        {
          if (cornerDir != dir)
          {
            // now loop over hi-lo corners
            SideIterator cornerSit;
            for (cornerSit.begin(); cornerSit.ok(); ++cornerSit)
            {
              // now isolate corner cells by chopping again
              if (cornerSit() == Side::Lo)
              {
                cornerGhost = allGhostDir;
                int chopPoint = dest[todi].smallEnd(cornerDir);
                tempBox = cornerGhost.chop(cornerDir, chopPoint);
              }
              else
              {
                tempBox = allGhostDir;
                int chopPoint = dest[todi].bigEnd(cornerDir)+1;
                cornerGhost = tempBox.chop(cornerDir, chopPoint);
              }

              // now we should have isolated the corner ghost cells
              // for this box in the given directions and faces.
              // then for each level/from DI on this processor,
              // see if they intersect
              for (vector<DataIndex>::iterator vli = vectorLevelOnProcDI.begin();
                  vli != vectorLevelOnProcDI.end(); ++vli)
              {

                // at this point, i know myprocID == fromProcID

                const DataIndex fromdi(*vli);
                Box testBox = level[fromdi];

                // skip intersections with mysef
                if (testBox != dest[todi])
                {
                  testBox.grow(a_ghost);

                  if ((testBox.bigEnd(0) < cornerGhost.smallEnd(0)) && isSorted)
                  {
                    // can skip rest cuz we haven't gotten
                    // to something interesting
                    continue;
                  }

                  Box fromBox;
                  if (faceSit() == Side::Lo)
                  {
                    fromBox = adjCellLo(level[fromdi],
                                        dir,
                                        a_ghost[dir]);
                  }
                  else
                  {
                    fromBox = adjCellHi(level[fromdi],
                                        dir,
                                        a_ghost[dir]);
                  }

                  // now proceed as normal

                  if (cornerGhost.intersectsNotEmpty(fromBox))
                  {
                    Box box(cornerGhost); // ??
                    box&=fromBox;   // ??

                    if (toProcID == myprocID)
                    {
                      // local move
                      // don't push back here!  or you will
                      // get two.
                      //     we already did it above...
                      //m_localMotionPlan.push_back(item);
                    }
                    else
                    {
                      MotionItem* item = new (s_motionItemPool.getPtr())
                        MotionItem(fromdi, todi, box);
                      if (item == NULL)
                      {
                        MayDay::Error("Out of Memory in copier::define");
                      }

                      item->procID = toProcID;
                      m_fromMotionPlan.push_back(item);
                    }
                  } // end if there's an intersection
                } // end if dest and from are not the same

                if ((testBox.smallEnd(0) > cornerGhost.bigEnd(0)) && isSorted)
                {
                  // can break out of loop, since we know that
                  // the smallEnd
                  // of all the remaining boxes are lexigraphically
                  // beyond this ghosted box.
                  break;
                }
              } // end loop over from boxes
            } // end loop over conrner hi-lo
          } // end if corner dir != dir
        } // end loop over corner directions
      } // end loop over dest faces
    } // end loop over dirs
  } // end loop over dest boxes
#endif

  // put periodic intersection checking in here for "to" boxes
  // (looking for "to" boxes which hang over periodic domain edge)
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

  // Here ends the so-called N-squared optimizations.  the rest is
  // unchanged. (ndk)

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
            // here's where we figure out the corner cells
            Box cornerGhost;
            Box grownDestBox(dest[toIndex]);
            grownDestBox.grow(a_ghost);

            // loop over dir
            for (int dir=0; dir<SpaceDim; dir++)
            {
              SideIterator faceSit;
              for (faceSit.begin(); faceSit.ok(); ++faceSit)
              {
                Box allGhostDir;
                Box tempBox;
                if (faceSit() == Side::Lo)
                {
                  allGhostDir = grownDestBox;
                  int chopPoint = dest[toIndex].smallEnd(dir);
                  tempBox = allGhostDir.chop(dir, chopPoint);
                }
                else
                {
                  tempBox = grownDestBox;
                  int chopPoint = dest[toIndex].bigEnd(dir) +1;
                  allGhostDir = tempBox.chop(dir, chopPoint);
                }
                // now we should have all of the ghost cells on
                // the hi or lo face in direction dir

                // now loop over directions for corners
                for (int cornerDir = 0; cornerDir<SpaceDim; cornerDir++)
                {
                  if (cornerDir != dir)
                  {
                    // now loop over hi-lo corners
                    SideIterator cornerSit;
                    for (cornerSit.begin(); cornerSit.ok(); ++cornerSit)
                    {
                      // now isolate corner cells by chopping again
                      if (cornerSit() == Side::Lo)
                      {
                        cornerGhost = allGhostDir;
                        int chopPoint = dest[toIndex].smallEnd(cornerDir);
                        tempBox = cornerGhost.chop(cornerDir, chopPoint);
                      }
                      else
                      {
                        tempBox = allGhostDir;
                        int chopPoint = dest[toIndex].bigEnd(cornerDir)+1;
                        cornerGhost = tempBox.chop(cornerDir, chopPoint);
                      }

                      // now we should have isolated the corner
                      // ghost cells for this box in the given
                      // directions and faces.  then for each
                      // level/from DI on this processor, see
                      //if they intersect

                      // don't bother with corners that are
                      // not near the boundary
                      if (!shrunkDomainBox.contains(cornerGhost))
                      {

                        // only look at ghost cells in the from
                        // box with the same faceDir and hi-lo
                        Box fromGhost;
                        if (faceSit() == Side::Lo)
                        {
                          fromGhost = adjCellLo(fromBox,
                                                dir,
                                                a_ghost[dir]);
                        }
                        else
                        {
                          fromGhost = adjCellHi(fromBox,
                                                dir,
                                                a_ghost[dir]);
                        }

                        // now need to loop over shift vectors
                        //and look at images -- only do
                        // this if fromGhost is near
                        // boundary as well
                        if (!shrunkDomainBox.contains(fromGhost))
                        {
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                          {
                            IntVect shiftVect(shiftIt()*shiftMult);
                            cornerGhost.shift(shiftVect);
                            if (cornerGhost.intersectsNotEmpty(fromGhost)) // rarely happens
                            {
                              Box intersectBox(cornerGhost);
                              intersectBox &= fromGhost;
                              Box toBox(intersectBox);
                              toBox.shift(-shiftVect);
                              MotionItem* item = new (s_motionItemPool.getPtr())
                                MotionItem(DataIndex(from()),
                                           DataIndex(toIndex),
                                           intersectBox, toBox);
                              if (item == NULL)
                              {
                                MayDay::Error("Out of Memory in copier::define");
                              }
                              if (toProcID == fromProcID) // local move
                              {
                                m_localMotionPlan.push_back(item);
                              }
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

                            cornerGhost.shift(-shiftVect);
                          } // end loop over shift vectors
                        } // end if fromGhost near boundary
                      } // end if corner ghost near bndry
                    } // end loop over corner hi-lo dir
                  } // end if corner dir != faceDir
                } // end loop over corner dirs
              } // end loop over face hi-lo
            } // end loop over face directions
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

          // now do the whole "isolate the corner" dance again.
          // there's probably a better way to do this...

          // loop over directions
          for (int dir = 0; dir<SpaceDim; dir++)
          {
            Box cornerGhost;
            Box grownDestBox = dest[to()];
            grownDestBox.grow(a_ghost);
            // loop over hi-lo faces
            SideIterator faceSit;
            for (faceSit.begin(); faceSit.ok(); ++faceSit)
            {
              Box allGhostDir;
              Box tempBox;
              if (faceSit() == Side::Lo)
              {
                allGhostDir = grownDestBox;
                int chopPoint = dest[to()].smallEnd(dir);
                tempBox = allGhostDir.chop(dir, chopPoint);
              }
              else
              {
                tempBox = grownDestBox;
                int chopPoint = dest[to()].bigEnd(dir) +1;
                allGhostDir = tempBox.chop(dir, chopPoint);
              }
              // now we should have all of the ghost cells on
              // the hi or lo face in direction dir

              // now loop over directions for corners
              for (int cornerDir = 0; cornerDir<SpaceDim; cornerDir++)
              {
                if (cornerDir != dir)
                {
                  // now loop over hi-lo corners
                  SideIterator cornerSit;
                  for (cornerSit.begin(); cornerSit.ok(); ++cornerSit)
                  {
                    // now isolate corner cells by chopping again
                    if (cornerSit() == Side::Lo)
                    {
                      cornerGhost = allGhostDir;
                      int chopPoint = dest[to()].smallEnd(cornerDir);
                      tempBox = cornerGhost.chop(cornerDir, chopPoint);
                    }
                    else
                    {
                      tempBox = allGhostDir;
                      int chopPoint = dest[to()].bigEnd(cornerDir)+1;
                      cornerGhost = tempBox.chop(cornerDir, chopPoint);
                    }

                    // now we should have isolated the corner
                    // ghost cells for this box in the given
                    // directions and faces.
                    // now loop over those "from" boxes which are not
                    // contained by the domain
                    for (int fromRef = 0; fromRef<periodicFromVect.size(); fromRef++)
                    {
                      DataIndex fromIndex = periodicFromVect[fromRef];
                      unsigned int fromProcID = level.procID(fromIndex);

                      // don't worry about anything which doesn't involve
                      // this proc
                      if (toProcID != myprocID && fromProcID != myprocID)
                      {
                        // do nothing
                      }
                      else
                      {
                        // grab ghost cells
                        Box fromBox;
                        if (faceSit() == Side::Lo)
                        {
                          fromBox = adjCellLo(level[fromIndex],
                                              dir,
                                              a_ghost[dir]);
                        }
                        else
                        {
                          fromBox = adjCellHi(level[fromIndex],
                                              dir,
                                              a_ghost[dir]);
                        }


                        // only do the shifting if these
                        // ghost cells are near the periodic
                        // boundary
                        if (!shrunkDomainBox.contains(fromBox))
                        {
                          // now need to loop over shift
                          // vectors and look at images
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                          {
                            IntVect shiftVect(shiftIt()*shiftMult);
                            cornerGhost.shift(shiftVect);
                            if (cornerGhost.intersectsNotEmpty(fromBox))
                            {
                              Box intersectBox(cornerGhost);
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

                            cornerGhost.shift(-shiftVect);
                          } // end loop over shift vectors
                        } // end if from box near boundary
                      } // end if either from box or to box are on this proc
                    } // end loop over "from" boxes
                  } // end loop over corner hi-lo
                } // end if corner dir != faceDir
              } // end loop over cornerDir
            } // end loop over face hi-lo
          } // end loop over face directions
        } // end if destination box is close to domain boundary
      } // end loop over destination boxes
    } // end if any of the "From" boxes were outside the domain

  } // end if we need to do anything for periodicity
}


CornerCopier::~CornerCopier()
{
  clear();
}

#include "NamespaceFooter.H"
