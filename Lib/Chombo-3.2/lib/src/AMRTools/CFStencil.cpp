#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  -*- Mode: C++; Modified: "Tue 02 Feb 1999 14:29:43 by graves"; -*-
#include "LayoutIterator.H"
#include "CFStencil.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

bool
CFStencil::isDefined() const
{
    return m_isDefined;
}

/**
   are there any interpolation points?  \\
   returns false if so.*/
bool
CFStencil::isEmpty() const
{
    CH_assert(m_isDefined);
    return m_fineIVS.isEmpty();
}

/** get fine intvects which need to be interpolated  \\
    This will be empty if isEmpty() returns true*/
const IntVectSet&
CFStencil::getFineIVS() const
{
    CH_assert(m_isDefined);
    return m_fineIVS;
}

/** get coarse intvects that underly fiIVS.
    This will be empty if isEmpty() returns true*/
const IntVectSet&
CFStencil::getCoarIVS() const
{
    CH_assert(m_isDefined);
    return m_coarIVS;
}

void
CFStencil::clear()
{
    m_fineIVS.define();
    m_coarIVS.define();
    setDefaultValues();
}
void
CFStencil::setDefaultValues()
{
    m_isDefined = false;
    m_direction = -777;
    m_isPacked  = false;
    m_packedBox = Box();
}
CFStencil::CFStencil()
{
    setDefaultValues();
}
CFStencil::CFStencil(const Box& a_b, int a_nComp)
{
  setDefaultValues();
  define(a_b, a_nComp);

}
void
CFStencil::define(const Box& a_b, int a_nComp)
{
  clear();
}

CFStencil::~CFStencil()
{
    clear();
}

CFStencil::CFStencil(
    const Box& a_fineDomain,
    const Box& a_grid,
    const DisjointBoxLayout& a_fineBoxes,
    const DisjointBoxLayout& a_coarBoxes,
    int a_refRatio,
    int a_direction,
    Side::LoHiSide a_hiorlo)
{
  ProblemDomain physDomain(a_fineDomain);
  setDefaultValues();
  define(physDomain,  a_grid, a_fineBoxes, a_coarBoxes,
         a_refRatio, a_direction, a_hiorlo);
}

CFStencil::CFStencil(
                     const ProblemDomain& a_fineDomain,
                     const Box& a_grid,
                     const DisjointBoxLayout& a_fineBoxes,
                     const DisjointBoxLayout& a_coarBoxes,
                     int a_refRatio,
                     int a_direction,
                     Side::LoHiSide a_hiorlo)
{
  setDefaultValues();
  define(a_fineDomain, a_grid, a_fineBoxes, a_coarBoxes,
         a_refRatio, a_direction, a_hiorlo);
}

void
CFStencil::define(
                  const Box& a_fineDomain,
                  const Box& a_grid,
                  const DisjointBoxLayout& a_fineBoxes,
                  const DisjointBoxLayout& a_coarBoxes,
                  int a_refRatio,
                  int a_direction,
                  Side::LoHiSide a_hiorlo)
{
  ProblemDomain physDomain(a_fineDomain);
  define(physDomain, a_grid, a_fineBoxes, a_coarBoxes, a_refRatio,
         a_direction, a_hiorlo);
}

void
CFStencil::define(
                  const ProblemDomain& a_fineDomain,
                  const Box& a_grid,
                  const DisjointBoxLayout& a_fineBoxes,
                  const DisjointBoxLayout& a_coarBoxes,
                  int a_refRatio,
                  int a_direction,
                  Side::LoHiSide a_hiorlo)
{
  m_isDefined = true;
  CH_assert(a_refRatio >= 1);
  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert((a_hiorlo == Side::Lo) ||
         (a_hiorlo == Side::Hi));
  CH_assert(!a_fineDomain.isEmpty());

  //set internal vars.  most of these are kept around
  //just to keep the class from having an identity crisis.
  m_direction = a_direction;
  m_hiorlo =  a_hiorlo;

  Box finebox = a_grid;


  //compute intvectset of all points on fine grid that
  //need to be interpolated

  //shift direction
  int hilo = sign(a_hiorlo);

  //create fine stencil
  Box edgebox;
  CH_assert((hilo ==1) || (hilo == -1));
  if (hilo == -1)
    {
      edgebox = adjCellLo(finebox,m_direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,m_direction,1);
    }
  edgebox = a_fineDomain & edgebox;
  if (!edgebox.isEmpty())
    {
      Box periodicTestBox(a_fineDomain.domainBox());
      if (a_fineDomain.isPeriodic())
        {
          for (int idir=0; idir<SpaceDim; idir++)
            {
              if (a_fineDomain.isPeriodic(idir))
                {
                  periodicTestBox.grow(idir,-1);
                }
            }
        }

      m_fineIVS.define(edgebox);
      LayoutIterator lit = a_fineBoxes.layoutIterator();
      for (lit.reset(); lit.ok(); ++lit)
        {
          m_fineIVS -= a_fineBoxes[lit()];
          // if periodic, also need to subtract periodic images
          // only do this IF we're periodic _and_ both boxes
          // adjoin the domain box boundary somewhere
          if (a_fineDomain.isPeriodic() && !periodicTestBox.contains(edgebox)
              && !periodicTestBox.contains(a_fineBoxes[lit()]))
            {
              ShiftIterator shiftIt = a_fineDomain.shiftIterator();
              IntVect shiftMult(a_fineDomain.domainBox().size());
              Box shiftedBox(a_fineBoxes[lit()]);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  m_fineIVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                } // end loop over periodic shift directions
            } // end if periodic
        }
    }

  //ivs where all coarse slopes are defined
  //== coarsened fine ivs
  m_coarIVS.define(m_fineIVS);
  m_coarIVS.coarsen(a_refRatio);
  // this is a trick to get around the lack of a IntVectSet intersection
  // operator which works with a ProblemDomain
  ProblemDomain coardom= coarsen(a_fineDomain, a_refRatio);
  Box domainIntersectBox = m_coarIVS.minBox();
  domainIntersectBox = coardom & domainIntersectBox;
  m_coarIVS &= domainIntersectBox;

  m_packedBox = m_fineIVS.minBox();
  if (m_fineIVS.numPts() == m_packedBox.numPts())
    {
      m_isPacked = true;
    }
  else
    {
      m_isPacked = false;
      m_packedBox = Box();
    }
}

void
CFStencil::define(
                  const ProblemDomain& a_fineDomain,
                  const Box& a_grid,
                  const Vector<Box>& a_periodicVector,
                  int a_refRatio,
                  int a_direction,
                  Side::LoHiSide a_hiorlo)
{
  m_isDefined = true;
  CH_assert(a_refRatio >= 1);
  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert((a_hiorlo == Side::Lo) ||
         (a_hiorlo == Side::Hi));
  CH_assert(!a_fineDomain.isEmpty());

  //set internal vars.  most of these are kept around
  //just to keep the class from having an identity crisis.
  m_direction = a_direction;
  m_hiorlo =  a_hiorlo;

  Box finebox = a_grid;


  //compute intvectset of all points on fine grid that
  //need to be interpolated

  //shift direction
  int hilo = sign(a_hiorlo);

  //create fine stencil
  Box edgebox;
  CH_assert((hilo ==1) || (hilo == -1));
  if (hilo == -1)
    {
      edgebox = adjCellLo(finebox,m_direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,m_direction,1);
    }
  edgebox = a_fineDomain & edgebox;

  if (edgebox.isEmpty()) return;


  int w1 = edgebox.smallEnd()[0];
  int w2 = edgebox.bigEnd()[0];
  m_fineIVS.define(edgebox);
  // moving window loop in i-direction (bvs)
  for (int i=0; i<a_periodicVector.size(); ++i)
    {
      const Box& b = a_periodicVector[i];
      if (b.bigEnd()[0] >= w1)
        {
          m_fineIVS -= b;
          if (b.smallEnd()[0] > w2)
          {
            i=a_periodicVector.size();
          }
        }

    }


  //ivs where all coarse slopes are defined
  //== coarsened fine ivs
  m_coarIVS.define(m_fineIVS);
  m_coarIVS.coarsen(a_refRatio);
  // this is a trick to get around the lack of a IntVectSet intersection
  // operator which works with a ProblemDomain
  ProblemDomain coardom= coarsen(a_fineDomain, a_refRatio);
  Box domainIntersectBox = m_coarIVS.minBox();
  domainIntersectBox = coardom & domainIntersectBox;
  m_coarIVS &= domainIntersectBox;

  m_packedBox = m_fineIVS.minBox();
  if (m_fineIVS.numPts() == m_packedBox.numPts())
    {
      m_isPacked = true;
    }
  else
    {
      m_isPacked = false;
      m_packedBox = Box();
    }
}

void CFStencil::buildPeriodicVector(Vector<Box>& a_periodicVector,
                                    const ProblemDomain& a_fineDomain,
                                    const DisjointBoxLayout& a_fineBoxes)
{
  Box periodicTestBox(a_fineDomain.domainBox());
  if (a_fineDomain.isPeriodic())
    {
      for (int idir=0; idir<SpaceDim; idir++)
        {
          if (a_fineDomain.isPeriodic(idir))
            {
              periodicTestBox.grow(idir,-1);
            }
        }
    }
  a_periodicVector.clear();
  a_periodicVector.reserve(a_fineBoxes.size());

  LayoutIterator lit = a_fineBoxes.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
    {
      const Box& box = a_fineBoxes[lit()];
      a_periodicVector.push_back(box);
      // if periodic, also need to add periodic images
      // only do this IF we're periodic and box
      // adjacent to  the domain box boundary somewhere
      if (a_fineDomain.isPeriodic()
          && !periodicTestBox.contains(box))
        {
          ShiftIterator shiftIt = a_fineDomain.shiftIterator();
          IntVect shiftMult(a_fineDomain.domainBox().size());
          Box shiftedBox(box);
          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
            {
              IntVect shiftVect = shiftMult*shiftIt();
              shiftedBox.shift(shiftVect);
              a_periodicVector.push_back(shiftedBox);
              shiftedBox.shift(-shiftVect);
            } // end loop over periodic shift directions
        } // end if periodic
    }
  a_periodicVector.sort();
}


CFStencil&
CFStencil::operator=(const CFStencil& a_ecfsIn)
{
    if (&a_ecfsIn != this)
    {
        m_direction= a_ecfsIn.m_direction;
        m_hiorlo = a_ecfsIn.m_hiorlo;
        m_dataIndex = a_ecfsIn.m_dataIndex;
        m_fineIVS = a_ecfsIn.m_fineIVS;
        m_coarIVS= a_ecfsIn.m_coarIVS;
    }
    return *this;
}
#include "NamespaceFooter.H"
