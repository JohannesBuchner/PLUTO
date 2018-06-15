#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBCFData.H"
#include "EBInterpolateF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "LayoutIterator.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBAlias.H"
#include <algorithm>

#include "NamespaceHeader.H"

/************************************/
EBCFData::~EBCFData()
{
}
/*********************/
void
EBCFData::
getExtrapSigns(IntVect& a_signs, const IntVect& a_corner, const Box& a_grid)
{
  const IntVect& loEnd = a_grid.smallEnd();
  const IntVect& hiEnd = a_grid.bigEnd();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_corner[idir] > hiEnd[idir])
        {
          a_signs[idir] = -1; //estrapolate from the low side if on the high end
        }
      else if (a_corner[idir] < loEnd[idir])
        {
          a_signs[idir] =  1; //estrapolate from the high side if on the low end
        }
      else
        {
#if CH_SPACEDIM==2
          //in 2d, one should not get here by construction
          pout() << "getExtrapSigns a_corner = " << a_corner << ", a_grid = " << a_grid << endl;
          MayDay::Error("getExtrapSigns: corner vof apparently not in the corner");

#elif CH_SPACEDIM==3
          //in 3d you can be along an edge.  only extrapolate in the other directions.
          a_signs[idir] = 0;
#else
          bogus_spacedim_compilation_error();
#endif
        }

    }
}
/************************************/
EBCFData::
EBCFData(const DisjointBoxLayout&       a_gridsFine,
         const DisjointBoxLayout&       a_gridsCoar,
         const EBISLayout&              a_ebislFine,
         const EBISLayout&              a_ebislCoar,
         const ProblemDomain&           a_domainCoar,
         const int&                     a_nref,
         const LayoutData<IntVectSet>&  a_cfivs,
         const EBIndexSpace* const      a_ebisPtr,
         bool a_doEBCFCrossing,
         bool a_doCornerEdgeIterators)
{
  m_refRat    = a_nref;
  m_gridsFine = a_gridsFine;
  m_gridsCoar = a_gridsCoar;
  m_ebislFine = a_ebislFine;
  m_ebislCoar = a_ebislCoar;

  m_domainCoar = a_domainCoar;
  m_domainFine = refine(m_domainCoar, m_refRat);


  if (m_ebislFine.getMaxCoarseningRatio() < m_refRat)
    {
      m_ebislFine.setMaxCoarseningRatio(m_refRat,a_ebisPtr);
    }

  //do EB-Crossing-CF interface tricks
  int nghostEBISL = 4;
  m_doEBCFCrossing = a_doEBCFCrossing;
  //might need this stuff for corners
  m_gridsCoarsenedFine = DisjointBoxLayout();
  coarsen(m_gridsCoarsenedFine, m_gridsFine, m_refRat);
  a_ebisPtr->fillEBISLayout(m_ebislCoarsenedFine,
                            m_gridsCoarsenedFine,
                            m_domainCoar, nghostEBISL);

  if (m_doEBCFCrossing)
    {
      //figure out where we need to do our EB-crossing-tricks (and if we need them at all)
      m_doEBCFCrossing = getEBCFIVS(a_cfivs);

      if (m_doEBCFCrossing)
        {
          defineLoHiIterators(a_cfivs);
        }

    }
  if (a_doCornerEdgeIterators)
    {
      defineEdCoIterators(a_cfivs);
    }

}
/************************************/
void
EBCFData::
getEdgeAndCornerIVS(IntVectSet& a_edgeIVS, IntVectSet& a_cornerIVS,
                    const Box& a_grid, const ProblemDomain& a_domain,
                    const IntVectSet& a_cfivsGrid)
{
  Box grownBox = grow(a_grid, 1);
  grownBox &= a_domain;
  IntVectSet ivsAll(grownBox);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box growBoxSide = a_grid;
      growBoxSide.grow(idir, 1);
      ivsAll -= growBoxSide;
    }

  a_cornerIVS = ivsAll;
  if (SpaceDim == 2)
    {
      //in 2d there are no edges
      a_edgeIVS = IntVectSet();
    }
  else
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit1; sit1.ok(); ++sit1)
            {
              Box faceBox = adjCellBox(a_grid, idir, sit1(), 1);
              for (int jdir = 0; jdir < SpaceDim; jdir++)
                {
                  if (jdir != idir)
                    {
                      for (SideIterator sit2; sit2.ok(); ++sit2)
                        {
                          Box edgeBox = adjCellBox(faceBox, jdir, sit2(), 1);
                          a_cornerIVS -= edgeBox;
                        }
                    }
                }
            }
        }
      a_edgeIVS = ivsAll - a_cornerIVS;
    }

  IntVectSet grownCFIVS= grow(a_cfivsGrid, 1);
  a_cornerIVS &= grownCFIVS;
  a_edgeIVS   &= grownCFIVS;
}
/************************************/
void EBCFData::
getEBCFIVSGrid(IntVectSet&                a_ebcfivs,
               const Box&                 a_grid,
               const int&                 a_idir,
               const Side::LoHiSide&      a_side,
               const IntVect&             a_diagGrow,
               const ProblemDomain&       a_domain,
               const IntVectSet&          a_cfivs,
               const EBISBox&             a_ebisBox)
{

  Box gridSide = adjCellBox(a_grid, a_idir, a_side, 1);
  Box grownBoxSide = gridSide;
  grownBoxSide &= a_domain;

  for (int jdir = 0; jdir < SpaceDim; jdir++)
    {
      if (jdir != a_idir)
        {
          grownBoxSide.grow(jdir, a_diagGrow[jdir]);
        }
    }
  grownBoxSide &= a_domain;

  a_ebcfivs = a_cfivs;
  for (int jdir = 0; jdir < SpaceDim; jdir++)
    {
      if (jdir != a_idir)
        {
          a_ebcfivs.grow(jdir, a_diagGrow[jdir]);
        }
    }

  IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(grownBoxSide);
  ivsIrreg.grow(3);
  ivsIrreg  &= a_domain;
  a_ebcfivs &= ivsIrreg;
  a_ebcfivs &= gridSide;

}
bool
EBCFData::
getEBCFIVS(const LayoutData<IntVectSet>&  a_cfivs)
{
  bool doEBCrossing = false;

  //now for the side ivs in each direction.   If nothing is found to be crossing,
  // then there is no EB crossing
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_ebcfivsLo[idir].define(m_gridsFine);
      m_ebcfivsHi[idir].define(m_gridsFine);
      IntVect ivgrow = IntVect::Zero;

      int ibox = 0;
      for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
        {
          Box grid  = m_gridsFine.get(dit());

          getEBCFIVSGrid(m_ebcfivsLo[idir][dit()], grid, idir, Side::Lo, ivgrow, m_domainFine, a_cfivs[dit()], m_ebislFine[dit()]);
          getEBCFIVSGrid(m_ebcfivsHi[idir][dit()], grid, idir, Side::Hi, ivgrow, m_domainFine, a_cfivs[dit()], m_ebislFine[dit()]);

          //keep track if this  crossing ever happens
          if ((!m_ebcfivsLo[idir][dit()].isEmpty()) ||
             (!m_ebcfivsHi[idir][dit()].isEmpty()))
            {
              doEBCrossing = true;
            }
          ibox++;
        }
    }
  //in the case of parallel, need to check if ANY of the procs
  //have ebcrossing
#ifdef CH_MPI

  int gatherint = 0;
  if (doEBCrossing) gatherint = 1;
  int idoebcf;
  MPI_Allreduce(&gatherint, &idoebcf, 1, MPI_INT,
                MPI_MAX, Chombo_MPI::comm);
  doEBCrossing = (idoebcf==1);

#endif

  return doEBCrossing;
}
/****/
void
EBCFData::
defineEdCoIterators(const LayoutData<IntVectSet>&  a_cfivs)
{
  m_vofItCorners.define(m_gridsFine);
  m_vofItEdges.define(  m_gridsFine);
  m_cornerIVS.define(   m_gridsFine);
  m_edgeIVS.define(     m_gridsFine);
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      getEdgeAndCornerIVS(m_edgeIVS[dit()], m_cornerIVS[dit()], m_gridsFine.get(dit()), m_domainFine, a_cfivs[dit()]);

      m_vofItCorners[dit()].define(m_cornerIVS[dit()], m_ebislFine[dit()].getEBGraph());
      m_vofItEdges[dit()].define(    m_edgeIVS[dit()], m_ebislFine[dit()].getEBGraph());

    }
}

void
EBCFData::
defineLoHiIterators(const LayoutData<IntVectSet>&  a_cfivs)
{
 for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_vofItEBCFLo[idir].define(m_gridsFine);
      m_vofItEBCFHi[idir].define(m_gridsFine);
      
      int ibox = 0;
      for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
        {
          //define all the holders
          const IntVectSet& loEBCFIVS = m_ebcfivsLo[idir][dit()];
          const IntVectSet& hiEBCFIVS = m_ebcfivsHi[idir][dit()];

          m_vofItEBCFLo[idir][dit()].define(loEBCFIVS, m_ebislFine[dit()].getEBGraph());
          m_vofItEBCFHi[idir][dit()].define(hiEBCFIVS, m_ebislFine[dit()].getEBGraph());
          ibox++;
        }

    }
}
/****/

#include "NamespaceFooter.H"
