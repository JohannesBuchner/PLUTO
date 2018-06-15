#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "MayDay.H"
#include "DataIterator.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "SliceSpec.H"
#include <list>
#include "CH_Timer.H"
#include "NamespaceHeader.H"

DisjointBoxLayout::DisjointBoxLayout()
  :BoxLayout()
{
}

DisjointBoxLayout::DisjointBoxLayout(const Vector<Box>& a_boxes,
                                     const Vector<int>& a_procIDs)
  :BoxLayout(a_boxes,a_procIDs)
{
  CH_assert(isDisjoint());
  computeNeighbors();

}

DisjointBoxLayout::DisjointBoxLayout(const Vector<Box>& a_boxes,
                                     const Vector<int>& a_procIDs,
                                     const ProblemDomain& a_physDomain)
  :BoxLayout(a_boxes,a_procIDs), m_physDomain(a_physDomain)
{
  CH_assert(isDisjoint());
  computeNeighbors(); // even though BoxLayout::close is virtual, virtual dispatch does not
                      // happen in constructors.
}

void
DisjointBoxLayout::define(const Vector<Box>& a_boxes,
                          const Vector<int>& a_procIDs)
{
  BoxLayout::define(a_boxes,a_procIDs);
  CH_assert(isDisjoint());

}

void
DisjointBoxLayout::define(const Vector<Box>& a_boxes,
                          const Vector<int>& a_procIDs,
                          const ProblemDomain& a_physDomain)
{
  m_physDomain = a_physDomain;
  BoxLayout::define(a_boxes,a_procIDs);
  CH_assert(isDisjoint());

}

void
DisjointBoxLayout::define(const BoxLayout& a_layout)
{
  ((BoxLayout*)this)->operator=(a_layout);
  if (*m_closed)
    {
      CH_assert(isDisjoint());
    }

}

void
DisjointBoxLayout::define(const BoxLayout& a_layout,
                          const ProblemDomain& a_physDomain)
{
  m_physDomain = a_physDomain;

  ((BoxLayout*)this)->operator=(a_layout);
  if (*m_closed)
    {
      CH_assert(isDisjoint());
      if (*m_sorted)
        {
          computeNeighbors();
        }
    }
}

void
DisjointBoxLayout::close()
{
  if (!*m_closed)  //do nothing if already closed
    {
      sort();
      CH_assert(isDisjoint());
      *m_closed = true;
      buildDataIndex();
      m_dataIterator = RefCountedPtr<DataIterator>(
                        new DataIterator(*this, m_layout));
      computeNeighbors();
    }
}
void
DisjointBoxLayout::closeNO()
{
  if (!*m_closed)  //do nothing if already closed
    {
      sort();
      CH_assert(isDisjoint());
      *m_closed = true;
      m_dataIterator = RefCountedPtr<DataIterator>(
                        new DataIterator(*this, m_layout));
      //computeNeighbors(); don't build neighbors
    }
}

void
DisjointBoxLayout::closeNoSort()
{
  if (!*m_closed)  //do nothing if already closed
    {
      //sort(); // don't sort...
      CH_assert(isDisjoint());
      *m_sorted = false;
      *m_closed = true;
      m_dataIterator = RefCountedPtr<DataIterator>(
                        new DataIterator(*this, m_layout));
      //computeNeighbors(); don't build neighbors
    }
}

void
DisjointBoxLayout::closeN(RefCountedPtr<Vector<Vector<std::pair<int, LayoutIndex> > > > neighbors)

{
  if (!*m_closed)  //do nothing if already closed
    {
      sort();
      CH_assert(isDisjoint());
      *m_closed = true;
      m_dataIterator = RefCountedPtr<DataIterator>(
                        new DataIterator(*this, m_layout));
      m_neighbors = neighbors;
    }
}

void DisjointBoxLayout::computeNeighbors()
{
  CH_TIME("DisjointBoxLayout::computeNeighbors");
  const Vector<Entry>& boxes = *m_boxes;
  int n = boxes.size();
  int maxI=0;
  for (int i=0; i<n; ++i)
    {
      maxI = Max(maxI, boxes[i].box.size(0));
    }
  std::list<std::pair<int, LayoutIndex> > periodicImages;
  if (!m_physDomain.isEmpty() && m_physDomain.isPeriodic())
    {
      Box periodicTestBox(m_physDomain.domainBox());
      for (int idir=0; idir<CH_SPACEDIM; idir++)
        {
          if (m_physDomain.isPeriodic(idir))
            {
              periodicTestBox.grow(idir,-1);
            }
        }

      LayoutIterator lit = layoutIterator();
      unsigned int index = 0;
      ShiftIterator shiftIt = m_physDomain.shiftIterator();
      for (lit.reset(); lit.ok(); ++lit,++index)
        {
          const Box& box = this->operator[](lit);
          // if periodic, also need to add periodic images
          // only do this IF we're periodic and box
          // adjacent to  the domain box boundary somewhere
          if (!periodicTestBox.contains(box))
            {
              //IntVect shiftMult(m_physDomain.domainBox().size());
              //Box shiftedBox(box);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  //IntVect shiftVect = shiftMult*shiftIt();
                  //shiftedBox.shift(shiftVect);

                  periodicImages.push_back(std::pair<int, LayoutIndex>(shiftIt.index(), lit()));
                  //shiftedBox.shift(-shiftVect);
                } // end loop over periodic shift directions
            } // end if periodic
        }
    }
  //  OK, for this work, we don't use a standard LayoutIterator, since I want to
  // be pushing the start of the window forward, which LayoutIterator is not built to do.
  unsigned int id = 0;
  unsigned int start = 0;
  unsigned int end   = size();
  m_neighbors = RefCountedPtr<Vector<Vector<std::pair<int, LayoutIndex > > > >(
            new Vector<Vector<std::pair<int, LayoutIndex> > >());
  DataIterator dit = dataIterator();
  m_neighbors->resize(size());
  LayoutIterator lit = layoutIterator();
  const Vector<LayoutIndex>& vecLayoutIndex = *(lit.m_indicies);

  for (DataIterator dit=dataIterator(); dit.ok(); ++dit)
    {
      Box gbox = get(dit());
      gbox.grow(1);
      Vector<std::pair<int, LayoutIndex> >& neighbors = (*m_neighbors)[dit().intCode()];
      int low =  gbox.smallEnd()[0] - maxI;
      int high = gbox.smallEnd()[0] + maxI + 1;
      for (id=start; id<end; ++id)
        {
          if (id != index(dit()))
          {
            //don't include yourself as neighbor
            const Box& b = boxes[id].box;
            const IntVect& s = b.smallEnd();
            if (s[0]< low) start = id+1;
            else if (s[0] > high) id=end;
            if (gbox.intersectsNotEmpty(b))
              {
                neighbors.push_back(std::pair<int, LayoutIndex>(-1, vecLayoutIndex[id]));
              }
          }
        }
      //now run through periodic boxes.
      if(!m_physDomain.isEmpty() && !m_physDomain.domainBox().contains(gbox))
	{
	  std::list<std::pair<int, LayoutIndex> >::iterator it;
	  for (it=periodicImages.begin(); it!=periodicImages.end(); ++it)
	    {
	      Box b = (*m_boxes)[(*it).second.intCode()].box;
	      m_physDomain.shiftIt(b, (*it).first);
	      if (gbox.intersectsNotEmpty(b))
		{
		  pair<int, LayoutIndex> entry(*it);
		  //HERE
		  entry.second = (*it).second;
		  neighbors.push_back(entry);
		}
	    }
	}
    }
}

void
DisjointBoxLayout::deepCopy(const BoxLayout& a_layout)
{
  BoxLayout::deepCopy(a_layout);
  if (!a_layout.isClosed())
    sort();
  CH_assert(isDisjoint());

}

void
DisjointBoxLayout::deepCopy(const DisjointBoxLayout& a_layout)
{
  BoxLayout::deepCopy(a_layout);
  m_physDomain = a_layout.m_physDomain;
  if (!a_layout.isClosed())
    sort();
  CH_assert(isDisjoint());
}

void
DisjointBoxLayout::deepCopy(const BoxLayout& a_layout,
                            const ProblemDomain& a_physDomain)
{
  m_physDomain = a_physDomain;
  BoxLayout::deepCopy(a_layout);
  if (!a_layout.isClosed())
    sort();

  CH_assert(isDisjoint());
}

void
DisjointBoxLayout::degenerate( DisjointBoxLayout& a_to,
                               const SliceSpec& a_sliceSpec ) const
{
  Vector<Box> boxes;
  boxes.reserve( this->size() );
  bool outofbounds;
  for (LayoutIterator it(this->layoutIterator()); it.ok(); ++it)
    {
      Box degenBox;
      (*this)[it()].degenerate( degenBox, a_sliceSpec, &outofbounds );
      if ( !outofbounds )
        {
          boxes.push_back( degenBox );
        }
    }

    // What happens now if boxes.size()==0?
    a_to.defineAndLoadBalance( boxes, 0 );
    a_to.close(); // Do we really need this?
}

void
DisjointBoxLayout::defineAndLoadBalance(const Vector<Box>& a_boxes,
                                        Vector<int> * a_procIDs)
{
  ProblemDomain bogusProbDomain;
  defineAndLoadBalance(a_boxes, a_procIDs, bogusProbDomain);

}

void
DisjointBoxLayout::defineAndLoadBalance(const Vector<Box>& a_boxes,
                                        Vector<int> * a_procIDs,
                                        const ProblemDomain& a_physDomain)
{
    CH_assert( (!a_procIDs) || (a_procIDs->size() == 0) );

    Vector<int> procIDs;
    procIDs.reserve( a_boxes.size() );
    LoadBalance( procIDs, a_boxes );
    if ( a_procIDs )
      {
        *a_procIDs = procIDs; // Could this ever be a performance bottleneck?
      }
    this->define( a_boxes, procIDs, a_physDomain );
}

bool
DisjointBoxLayout::isDisjoint() const
{
  // this is a holder for boxes in the periodic case -- if a box
  // needs to be checked for its periodic images, then save it
  // and come back to it later (to preserve O(N) algorithm)
  Vector<Box> periodicCheckBoxes;
  for (unsigned int i=0; i<size(); ++i)
  {
    const Box& a = m_boxes->operator[](i).box;
    // in the periodic case, no box may be larger than the domain
    if (m_physDomain.isPeriodic())
      {
        if (m_physDomain.domainBox().size() < a.size())
          return false;
        // if a extends out of basic domain box, save it for
        // periodic checking later
        if (!m_physDomain.domainBox().contains(a))
          periodicCheckBoxes.push_back(a);
      }
    for (unsigned int j=i+1; j<size(); ++j)
    {
      bool disjoint = false;
      const Box& b = m_boxes->operator[](j).box;
      D_EXPR6(disjoint = disjoint || a.bigEnd(0) < b.smallEnd(0) || a.smallEnd(0) > b.bigEnd(0),
              disjoint = disjoint || a.bigEnd(1) < b.smallEnd(1) || a.smallEnd(1) > b.bigEnd(1),
              disjoint = disjoint || a.bigEnd(2) < b.smallEnd(2) || a.smallEnd(2) > b.bigEnd(2),
              disjoint = disjoint || a.bigEnd(3) < b.smallEnd(3) || a.smallEnd(3) > b.bigEnd(3),
              disjoint = disjoint || a.bigEnd(4) < b.smallEnd(4) || a.smallEnd(4) > b.bigEnd(4),
              disjoint = disjoint || a.bigEnd(5) < b.smallEnd(5) || a.smallEnd(5) > b.bigEnd(5));
      if (!disjoint)
        {
          return false;
        }
      else if (b.smallEnd(0) > a.bigEnd(0))
        {
          //can skip the rest, since we know that the smallEnd of all
          // the remaining boxes are going to also return true for this test.
          j=size();
        }
    }
  }

  // now come back to boxes we need to check for periodic images
  if (periodicCheckBoxes.size() > 0)
    {
      ShiftIterator shiftIt = m_physDomain.shiftIterator();
      IntVect shiftMult = m_physDomain.domainBox().size();
      for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
        {
          IntVect shiftVect(shiftMult*shiftIt());
          for (unsigned int i=0; i<periodicCheckBoxes.size(); i++)
            {
              Box& a = periodicCheckBoxes[i];
              a.shift(shiftVect);
              // only check for intersections with other boxes if
              // the shifted box intersects the computational domain.
              if (a.intersects(m_physDomain.domainBox()))
                {
                  // now repeat disjointness test with shifted box
                  for (unsigned int j=0; j<size(); ++j)
                    {
                      bool disjoint = false;
                      const Box& b = m_boxes->operator[](j).box;
                      D_EXPR6(disjoint = disjoint || a.bigEnd(0) < b.smallEnd(0) || a.smallEnd(0) > b.bigEnd(0),
                              disjoint = disjoint || a.bigEnd(1) < b.smallEnd(1) || a.smallEnd(1) > b.bigEnd(1),
                              disjoint = disjoint || a.bigEnd(2) < b.smallEnd(2) || a.smallEnd(2) > b.bigEnd(2),
                              disjoint = disjoint || a.bigEnd(3) < b.smallEnd(3) || a.smallEnd(3) > b.bigEnd(3),
                              disjoint = disjoint || a.bigEnd(4) < b.smallEnd(4) || a.smallEnd(4) > b.bigEnd(4),
                              disjoint = disjoint || a.bigEnd(5) < b.smallEnd(5) || a.smallEnd(5) > b.bigEnd(5));
                      if (!disjoint)
                        {
                          return false;
                        }
                      else if (b.smallEnd(0) > a.bigEnd(0))
                        {
                          //can skip the rest, since we know that the smallEnd
                          // of all the remaining boxes are going to also return
                          // true for this test.
                          j=size();
                        }
                    } // end loop over j
                } // end if shifted box intersects computational domain
              // shift box back
              a.shift(-shiftVect);
            } //end loop over i
        } // end loop over shift directions
    } // end if we need to check periodic images

  return true;
}

bool
DisjointBoxLayout::checkPeriodic(const ProblemDomain& a_domain) const
{
  bool goodMatch = true;
  if (a_domain.isPeriodic())
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          if (a_domain.isPeriodic(dir))
          {
            const Box& thisDomBox = m_physDomain.domainBox();
            const Box& inputDomBox = a_domain.domainBox();

            // make sure both are periodic in direction dir
            if (m_physDomain.isPeriodic(dir))
              {
                // both domains are periodic in direction dir -- check
                // that period is the same as well
                if ((thisDomBox.smallEnd(dir) != inputDomBox.smallEnd(dir))
                    ||(thisDomBox.bigEnd(dir) != inputDomBox.bigEnd(dir)))
                  {
                    // periods don't agree
                    goodMatch = false;
                  } // end if periods don't agree
              } // end if both are periodic in direction dir
            else
              {
                // a_domain is periodic, but m_physDomain is not!
                goodMatch = false;
              }
          } else if (m_physDomain.isPeriodic(dir))
          {

            // in this case, a_domain is not periodic, but m_physDomain is!
            goodMatch = false;
          }
        } // end loop over directions

    } // end if a_domain is periodic in any direction
  else if (m_physDomain.isPeriodic())
    {
      // if a_domain is not periodic, then m_physdomain better not be either
      goodMatch = false;
    }

  return goodMatch;
}

bool
DisjointBoxLayout::checkDomains(const DisjointBoxLayout& a_dbl) const
{
  const ProblemDomain& dblDomain = a_dbl.physDomain();
  return checkPeriodic(dblDomain);
}

// Global functions
// ================

void
coarsen(DisjointBoxLayout& a_output, const DisjointBoxLayout& a_input,
        int a_refinement)
{
  if (a_input.size() == 0)
    {
      a_output = a_input;
      return;
    }
  CH_assert(a_input.coarsenable(a_refinement));
  if (!a_input.isClosed())
    {
      MayDay::Error("input to coarsen must be called with closed BoxLayout");
    }
  if (a_output.isClosed())
    {
      MayDay::Error("output of coarsen must be called on open BoxLayout");
    }

  // copy first, then coarsen everything
  // a_output.deepCopy(a_input);
  a_output.m_boxes      = RefCountedPtr<Vector<Entry> >(new Vector<Entry>(*(a_input.m_boxes)));
  a_output.m_layout     = a_input.m_layout;
#ifdef CH_MPI
  a_output.m_dataIndex  = a_input.m_dataIndex;
#endif
  // now coarsen the physDomain
  a_output.m_physDomain = coarsen(a_input.m_physDomain, a_refinement);

  Vector<Entry>& boxes = *(a_output.m_boxes);
  int j=0;
  for (int i=0 ; i<=(int)boxes.size()-4; i+=4)
    {
      boxes[i].box.coarsen(a_refinement);
      boxes[i+1].box.coarsen(a_refinement);
      boxes[i+2].box.coarsen(a_refinement);
      boxes[i+3].box.coarsen(a_refinement);
      j+=4;
    }
  for (; j<boxes.size(); j++) boxes[j].box.coarsen(a_refinement);

  //    for (LayoutIterator it(a_input.layoutIterator()); it.ok(); ++it)
  //      {
  //        a_output.ref(it()).coarsen(a_refinement);
  //     }
  a_output.closeN(a_input.m_neighbors);
}

// we have an easier time with refine, since we know that refinement will
// not change the state of a sort, but, we will play it safe for now
// until this function shows up in the profiler.

void refine(DisjointBoxLayout& a_output,
            const DisjointBoxLayout& a_input,
            int a_refinement)
{
  if (!a_input.isClosed())
    {
      MayDay::Error("input to refine must be called with closed BoxLayout");
    }
  if (a_output.isClosed())
    {
      MayDay::Error("output of refine must be called on open BoxLayout");
    }

  // first copy, then refine everything
  a_output.deepCopy(a_input);

  // start by refining the physDomain
  a_output.m_physDomain = refine(a_input.m_physDomain,a_refinement);

  for (int ivec = 0; ivec < a_output.m_boxes->size(); ivec++)
    {
      (*a_output.m_boxes)[ivec].box.refine(a_refinement);
    }
  a_output.closeN(a_input.m_neighbors);
}

void
adjCellLo(DisjointBoxLayout& a_output,
          const DisjointBoxLayout& a_input,
          int a_dir, int a_len)
{
  if (!a_input.isClosed())
    {
      MayDay::Error("input to adjCellLo must be called with closed BoxLayout");
    }
  if (a_output.isClosed())
    {
      MayDay::Error("output of adjCellLo must be called on open BoxLayout");
    }

  // copy first, then adjCellLo everything
  a_output.deepCopy(a_input);

  for (int ivec = 0; ivec < a_output.m_boxes->size(); ivec++)
    {
      Box startBox = (*a_output.m_boxes)[ivec].box;
      (*a_output.m_boxes)[ivec].box = adjCellLo(startBox ,a_dir, a_len);
    }
  // shouldn't need to resort anything
  a_output.closeNoSort();
}

void
adjCellHi(DisjointBoxLayout& a_output,
          const DisjointBoxLayout& a_input,
          int a_dir, int a_len)
{
  if (!a_input.isClosed())
    {
      MayDay::Error("input to adjCellHi must be called with closed BoxLayout");
    }
  if (a_output.isClosed())
    {
      MayDay::Error("output of adjCellHi must be called on open BoxLayout");
    }

  // copy first, then adjCellHi everything
  a_output.deepCopy(a_input);
  for (int ivec = 0; ivec < a_output.m_boxes->size(); ivec++)
    {
      Box startBox = (*a_output.m_boxes)[ivec].box;
      (*a_output.m_boxes)[ivec].box = adjCellHi(startBox ,a_dir, a_len);
    }
  // don't want to sort here
  a_output.closeNoSort();
}

const
ProblemDomain&
DisjointBoxLayout::physDomain() const
{
  return m_physDomain;
}
#include "NamespaceFooter.H"
