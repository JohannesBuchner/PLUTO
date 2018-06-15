#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeCFIVS.cpp
// petermc, 1 May 2001
#include "NodeCFIVS.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "NamespaceHeader.H"
using std::cerr;
using std::endl;

bool
NodeCFIVS::isDefined() const
{
  return(m_isdefined);
}

const IntVectSet&
NodeCFIVS::getFineIVS() const
{
  CH_assert(m_isdefined);
  return m_fineInterpIVS;
}

void
NodeCFIVS::setDefaultValues()
{
  m_isdefined = false;
  m_fineInterpIVS.define();
  m_packed = false;
}

NodeCFIVS::NodeCFIVS()
{
  setDefaultValues();
}

NodeCFIVS::~NodeCFIVS()
{
  setDefaultValues();
}

NodeCFIVS::NodeCFIVS(
                     const Box& a_domain,
                     const Box& a_box,
                     const DisjointBoxLayout& a_levelBoxes,
                     int a_idir,
                     Side::LoHiSide a_hiorlo)
{
  setDefaultValues();
  ProblemDomain probdomain(a_domain);
  define(probdomain, a_box, a_levelBoxes, a_idir, a_hiorlo);
}

NodeCFIVS::NodeCFIVS(
                     const ProblemDomain& a_domain,
                     const Box& a_box,
                     const DisjointBoxLayout& a_levelBoxes,
                     int a_idir,
                     Side::LoHiSide a_hiorlo)
{
  setDefaultValues();
  define(a_domain, a_box, a_levelBoxes, a_idir, a_hiorlo);
}

void
NodeCFIVS::define(
                  const Box& a_domain,
                  const Box& a_box,
                  const DisjointBoxLayout& a_levelBoxes,
                  int a_idir,
                  Side::LoHiSide a_hiorlo)
{
  ProblemDomain probdomain(a_domain);
  define(probdomain, a_box, a_levelBoxes, a_idir, a_hiorlo);
}

void
NodeCFIVS::define(
                  const ProblemDomain& a_domain,
                  const Box& a_box,
                  const DisjointBoxLayout& a_levelBoxes,
                  int a_idir,
                  Side::LoHiSide a_hiorlo)
{
  m_isdefined = true;
  CH_assert(a_idir >= 0);
  CH_assert(a_idir < SpaceDim);
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_levelBoxes.checkPeriodic(a_domain));

  // First set faceCells to the CELLS adjacent to a_box,
  // on side a_hiorlo in direction a_idir.
  Box faceCells;
  switch (a_hiorlo)
    {
    case Side::Lo:
      {
        faceCells = adjCellLo(a_box, a_idir, 1);
        break;
      }
    case Side::Hi:
      {
        faceCells = adjCellHi(a_box, a_idir, 1);
        break;
      }
    default:
      {
        cerr << "NodeCFIVS::define(): bogus side" << endl;
        abort();
      }
    }
  // If this face of the box is along physical boundary:  faceboxIn = Empty.
  // Otherwise:  faceboxIn = faceCells.
  // If a_domain is periodic in direction a_idir:  faceboxIn = faceCells.
  Box faceboxIn = faceCells & a_domain;
  // Initialize m_fineInterpIVS with CELLs of faceCells that are in the domain.
  m_fineInterpIVS.define(faceboxIn);
  // Empty faceboxIn means that this face is on physical boundary.
  if (! faceboxIn.isEmpty())
    {
      // Expand m_fineInterpIVS by one CELL in directions other than a_idir.
      for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
        if (idirOther != a_idir) m_fineInterpIVS.grow(idirOther, 1);

      Box periodicTestBox(a_domain.domainBox());
      if (a_domain.isPeriodic())
        for (int idirPeriodic = 0; idirPeriodic < SpaceDim; idirPeriodic++)
          if (a_domain.isPeriodic(idirPeriodic))
            periodicTestBox.grow(idirPeriodic, -1);

      IntVect shiftMult = a_domain.domainBox().size();

      // Remove all CELLs from m_fineInterpIVS that are also in other boxes
      // at this level.
      LayoutIterator lit = a_levelBoxes.layoutIterator();
      for (lit.reset(); lit.ok(); ++lit)
        {
          Box shiftedBox = a_levelBoxes[lit()];
          m_fineInterpIVS -= shiftedBox;
          // added by petermc, 3 Feb 2003
          // only do this IF we're periodic _and_ both boxes
          // adjoin the domain box boundary somewhere
          if (a_domain.isPeriodic() &&
              !periodicTestBox.contains(faceboxIn) &&
              !periodicTestBox.contains(shiftedBox))
            {
              ShiftIterator shiftIt = a_domain.shiftIterator();
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect(shiftMult*shiftIt());
                  shiftedBox.shift(shiftVect);
                  m_fineInterpIVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }

      // Now figure out what the appropriate NODEs are.
      // Since IntVectSet is based on CELLs, not NODEs, we don't use
      // NODEs explicitly but we do figure out their indices.

      if (! m_fineInterpIVS.isEmpty())
        {
          // Get all NODEs on this face.
          // faceNodes is a CELL-centered box but actually contains
          // indices of the NODEs we want.
          Box faceNodes = faceCells;
          // To get indices of NODEs from CELLs, we need to:
          // (1) shift up by unit vector in direction a_idir, if on low side
          //     (not necessary if on high side)
          // (2) expand by one on high side in all non-a_idir directions.

          // (1) If on low side, shift up by unit vector in direction a_idir,
          // in order to get the correct indices for NODEs.
          if (a_hiorlo == Side::Lo)
            faceNodes.shift(a_idir, 1);

          // (2) Expand faceNodes by one on high side in all directions other
          // than a_idir, in order to get all NODEs.
          for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
            if (idirOther != a_idir)
              faceNodes.growHi(idirOther, 1);

          // Set domainInteriorNodes to be the set of all NODEs in interior
          // of a_domain; i.e., not on physical boundary.
          ProblemDomain domainInteriorNodes(a_domain);
          for (int idirAll = 0; idirAll < SpaceDim; idirAll++)
            domainInteriorNodes.growLo(idirAll, -1);

          // Remove from faceNodes any NODEs that are not interior NODEs
          // of a_domain.  (The removed nodes are on the physical boundary.)
          // In periodic directions, faceNodes is not changed.
          faceNodes &= domainInteriorNodes;

          // Now we are done with faceNodes.
          // faceNodes contains the indices of all NODEs that are
          // on this face of this box, and not on the physical boundary.
          // Example:  box [32:39, 24:39, 16:31] face 1 low,
          // faceNodes = [32:40, 24, 16:32].

          // At this point, m_fineInterpIVS contains the indices of
          // all CELLs that are just outside this face of this box,
          // and not contained in any other box.

          if (faceNodes.isEmpty())
            {
              // set to empty IntVectSet
              m_fineInterpIVS.define();
              // Something is wrong, because m_fineInterpIVS is not empty.
//                cerr << "NodeCFIVS::define(): problem" << endl;
//                cerr << "Side " << a_hiorlo << " in direction " << a_idir
//                     << " of box " << a_box << " in domain " << a_domain << endl;
//                cerr << "Empty face contains nodes " << m_fineInterpIVS << endl;
//                abort();
            }
          else
            {
              Box faceNodesExtra = faceNodes;
              for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                if (idirOther != a_idir)
                  faceNodesExtra.grow(idirOther, 1);

              // Now m_fineInterpIVS contains indices of CELLs.
              // We need to get the corresponding indices of NODEs.

              // First get correct indices for NODEs in dimension a_idir,
              // by shifting up by 1 in case of low side,
              // or doing nothing in case of high side.
              if (a_hiorlo == Side::Lo)
                m_fineInterpIVS.shift(BASISV(a_idir));

              // Get correct indices for NODEs in the other dimensions,
              // by taking the union with the whole set shifted up by 1 in
              // each dimension.

              // The complicated code below simply does the following:
              // for all idirOther != a_idir,
              //    m_fineInterpIVS |= (m_fineInterpIVS + BASISV(idirOther));
              for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                if (idirOther != a_idir)
                  {
                    IntVectSet interpIVSLo(m_fineInterpIVS);
                    // The only IntVectSet::shift() function is one that
                    // modifies.  You cannot do
                    // interpIVSHi = shift(m_fineInterpIVS).
                    IntVectSet interpIVSHi(m_fineInterpIVS);
                    interpIVSHi.shift(BASISV(idirOther));

                    // m_fineInterpIVS |= interpIVSHi;
                    // IntVectSet |= may be broken.  In any case, it's slow.
                    m_fineInterpIVS.define(faceNodesExtra);
                    // Now figure out what to eliminate.
                    IntVectSet elimStuff(faceNodesExtra);
                    elimStuff -= interpIVSLo;
                    elimStuff -= interpIVSHi;
                    m_fineInterpIVS -= elimStuff;

                    // Summarizing:
                    // m_fineInterpIVS == faceNodes -
                    //   (faceNodes - m_fineInterpIVS - shift(m_fineInterpIVS))
                    // == m_fineInterpIVS | shift(m_fineInterpIVS).
                  }

              m_fineInterpIVS &= faceNodes;
            }
        }
    }
  m_fineInterpIVS.compact();
  if (m_fineInterpIVS.isEmpty())
    {
      m_empty = true;
      m_packed = false;
    }
  else
    {
      m_empty = false;
      Box container(m_fineInterpIVS.minBox());
      if (container.numPts() == m_fineInterpIVS.numPts())
        {
          m_packCount++;
          m_packed = true;
          m_packedBox = container;
        }
      else
        {
          m_sparseCount++;
          m_packed = false;
        }
    }
}

long long NodeCFIVS::m_packCount = 0;

long long NodeCFIVS::m_sparseCount = 0;

#include "NamespaceFooter.H"
