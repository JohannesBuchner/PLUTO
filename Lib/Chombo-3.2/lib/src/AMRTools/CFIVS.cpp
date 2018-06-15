#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LayoutIterator.H"
#include "DataIterator.H"

#include "CFIVS.H"
#include "CH_Timer.H"
#include "LoHiSide.H"
#include "ProblemDomain.H"
#include "NeighborIterator.H"
#include "DisjointBoxLayout.H"
#include "NamespaceHeader.H"

//--Definitions for static member data

long long CFIVS::s_packCount = 0;
long long CFIVS::s_sparseCount = 0;


/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/

CFIVS::CFIVS()
{
  setDefaultValues();
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/

CFIVS::~CFIVS()
{
  setDefaultValues();   // Why?
}

/*--------------------------------------------------------------------*/
//  Explicit define based on an IntVectSet
/** Using this means the cells to store in the IVS have been
 *  identified elsewhere!
 *  \param[in]  a_IVS   IVS to use
 *//*-----------------------------------------------------------------*/

void
CFIVS::define(const IntVectSet & a_IVS)
{
  decrementCounts();
  if (a_IVS.isEmpty())  // Avoid compacting
    {
      m_IVS.define();
      m_empty = true;
      m_packed = false;
    }
  else
    {
      m_IVS.define(a_IVS);
      packIVS();
    }
  m_defined = true;
}

/*--------------------------------------------------------------------*/
//  General define for any box of ghost cells using NeighborIterators
/** This routine works on the ghost cells of the box given by
 *  a_grids[a_dataIndex].  The set of ghost cells to consider is given
 *  by a_ghostBox.  That a_ghostBox has any relation to
 *  a_grids[a_dataIndex] is up to you to ensure.
 *  \param[in]  a_dataIndex
 *                      DataIndex for the box that the ghost cells are
 *                      associated with.
 *  \param[in]  a_grids Layout of boxes
 *  \param[in]  a_ghostBox
 *                      Box containing all cells to consider as
 *                      possible ghosts that do not overlap real cells
 *                      on the level.
 *//*-----------------------------------------------------------------*/

void
CFIVS::define(const DataIndex&         a_dataIndex,
              const DisjointBoxLayout& a_grids,
              const Box&               a_ghostBox)
{
  CH_assert(a_grids.isSorted());  //**Not sufficient.  Need neighbours defined!
  decrementCounts();

  m_IVS.define(a_ghostBox);
  m_IVS -= a_grids[a_dataIndex];
  NeighborIterator nit(a_grids);
  for (nit.begin(a_dataIndex); nit.ok(); ++nit)
    {
      m_IVS -= nit.box();
    }

  packIVS();
  m_defined = true;
}

/*--------------------------------------------------------------------*/
//  Coarsen the stored cells
/*--------------------------------------------------------------------*/

void CFIVS::coarsen(int a_ref)
{
  CH_assert(m_defined);
  if (!m_empty)
    {
      m_IVS.coarsen(a_ref);
      m_packedBox.coarsen(a_ref);
    }
}

/*--------------------------------------------------------------------*/
//  For internal use
/*--------------------------------------------------------------------*/

void
CFIVS::setDefaultValues()
{
  m_IVS.define();
  m_packed = false;
  m_empty = true;
  m_defined = false;
}

/*--------------------------------------------------------------------*/
//  Decrement counts (during redefine)
/*--------------------------------------------------------------------*/

void
CFIVS::decrementCounts()
{
  if (m_defined && !m_empty)
    {
      if (m_packed)
        {
          --s_packCount;
        }
      else
        {
          --s_sparseCount;
        }
    }
  CH_assert(s_packCount >= 0);
  CH_assert(s_sparseCount >= 0);
}

/*--------------------------------------------------------------------*/
//  Pack the IVS if possible
/*--------------------------------------------------------------------*/

void
CFIVS::packIVS()
{
  m_IVS.compact();
  if (m_IVS.isEmpty())
    {
      m_empty = true;
      m_packed = false;
    }
  else
    {
      m_empty = false;
      m_packedBox = m_IVS.minBox();
      if (m_packedBox.numPts() == m_IVS.numPts())
        {
          ++s_packCount;
          m_packed = true;
        }
      else
        {
          ++s_sparseCount;
          m_packed = false;
        }
    }
}


/*==============================================================================
 * Legacy members
 *
 * New uses of CFIVS should avoid these routines
 *============================================================================*/


CFIVS::CFIVS(const Box&               a_domain,
             const Box&               a_boxIn,
             const DisjointBoxLayout& a_fineBoxes,
             int                      a_direction,
             Side::LoHiSide           a_hiorlo)
{
  setDefaultValues();
  ProblemDomain probdomain(a_domain);

  define(probdomain, a_boxIn, a_fineBoxes, a_direction, a_hiorlo);

}

CFIVS::CFIVS(const ProblemDomain&     a_domain,
             const Box&               a_boxIn,
             const DisjointBoxLayout& a_fineBoxes,
             int                      a_direction,
             Side::LoHiSide           a_hiorlo)
{
  setDefaultValues();

  define(a_domain, a_boxIn, a_fineBoxes, a_direction, a_hiorlo);

}

void CFIVS::define(const Box&               a_domain,
                   const Box&               a_boxIn,
                   const DisjointBoxLayout& a_fineBoxes,
                   int                      a_direction,
                   Side::LoHiSide           a_hiorlo)
{
  ProblemDomain probdomain(a_domain);

  define(probdomain, a_boxIn, a_fineBoxes, a_direction, a_hiorlo);
}

void CFIVS::define(const ProblemDomain&     a_domain,
                   const Box&               a_boxIn,
                   const DisjointBoxLayout& a_fineBoxes,
                   int                      a_direction,
                   Side::LoHiSide           a_hiorlo)
{
  CH_TIME("CFIVS::define(slow)");
  m_defined = true;

  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_domain.contains(a_boxIn));
  CH_assert(a_fineBoxes.checkPeriodic(a_domain));
  CH_assert((a_hiorlo == Side::Lo) || (a_hiorlo == Side::Hi));

  // create fine stencil
  Box finebox = a_boxIn;
  Box edgebox;

  if (a_hiorlo == Side::Lo)
    {
      edgebox = adjCellLo(finebox,a_direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,a_direction,1);
    }

  edgebox &= a_domain;

  if (!edgebox.isEmpty())
    {
      m_IVS.define(edgebox);

      LayoutIterator lit = a_fineBoxes.layoutIterator();
      Box periodicTestBox(a_domain.domainBox());

      if (a_domain.isPeriodic())
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (a_domain.isPeriodic(idir))
                periodicTestBox.grow(idir,-1);
            }
        }

      for (lit.reset(); lit.ok(); ++lit)
        {
          m_IVS -= a_fineBoxes[lit()];

          // only do this IF we're periodic _and_ both boxes
          // adjoin the domain box boundary somewhere
          if (a_domain.isPeriodic() && !periodicTestBox.contains(edgebox)
              && !periodicTestBox.contains(a_fineBoxes[lit()]))
            {
              ShiftIterator shiftIt = a_domain.shiftIterator();
              IntVect shiftMult = a_domain.domainBox().size();
              Box shiftedBox = a_fineBoxes[lit()];

              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect(shiftMult*shiftIt());

                  shiftedBox.shift(shiftVect);
                  m_IVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }
    }

  packIVS();
}


//-----------------------------------------------------------------------
void CFIVS::define(const ProblemDomain&     a_domain,
                   const Box&               a_boxIn,
                   const Vector<Box>&       a_periodicfineBoxes,
                   int                      a_direction,
                   Side::LoHiSide           a_hiorlo)
{
  CH_TIME("CFIVS::define(fast)");
  m_defined = true;

  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_domain.contains(a_boxIn));
  CH_assert((a_hiorlo == Side::Lo) || (a_hiorlo == Side::Hi));

  // create fine stencil
  Box finebox = a_boxIn;
  Box edgebox;

  if (a_hiorlo == Side::Lo)
    {
      edgebox = adjCellLo(finebox,a_direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,a_direction,1);
    }

  edgebox &= a_domain;

  if (edgebox.isEmpty())
    {
      m_packed = false;
      return;
    }

  m_IVS.define(edgebox);
  int w1 = edgebox.smallEnd()[0];
  int w2 = edgebox.bigEnd()[0];
  for (int i=0; i<a_periodicfineBoxes.size(); ++i)
    {
      const Box& b = a_periodicfineBoxes[i];
      if (b.bigEnd()[0] >= w1)
        {
          m_IVS -= b;
          if (b.smallEnd()[0] > w2)
          {
            i=a_periodicfineBoxes.size();
          }
        }
    }

  packIVS();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void CFIVS::define(const ProblemDomain&     a_domain,
                   const DataIndex&         a_dataIndex,
                   const DisjointBoxLayout& a_grids,
                   int                      a_direction,
                   Side::LoHiSide           a_hiorlo)
{
  CH_TIME("CFIVS::define(NeighborIterator-version)");
  m_defined = true;

  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert(!a_domain.isEmpty());
//  CH_assert(a_domain.contains(a_boxIn));
  CH_assert((a_hiorlo == Side::Lo) || (a_hiorlo == Side::Hi));

  // Create fine stencil
  Box finebox = a_grids[a_dataIndex];
  Box edgebox;

  if (a_hiorlo == Side::Lo)
    {
      edgebox = adjCellLo(finebox,a_direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,a_direction,1);
    }

  edgebox &= a_domain;

  if (edgebox.isEmpty())
    {
      m_packed = false;
      return;
    }

  m_IVS.define(edgebox);
  int w1 = edgebox.smallEnd()[0];
  int w2 = edgebox.bigEnd()[0];
  NeighborIterator nit(a_grids);
  bool isSorted = a_grids.isSorted();
  for (nit.begin(a_dataIndex); nit.ok(); ++nit)
    {
      const Box& b = nit.box();
      if (b.bigEnd()[0] >= w1)
        {
          m_IVS -= b;
          if ((b.smallEnd()[0] > w2) && isSorted)
          {
            // Needed?
            break;
          }
        }
    }

  packIVS();
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
