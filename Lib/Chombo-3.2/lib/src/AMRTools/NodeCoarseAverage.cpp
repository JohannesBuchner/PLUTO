#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeCoarseAverage.cpp
// added from CoarseAverage
// petermc, 28 Nov 2000
// petermc, 11 Apr 2002, added using std::cout and using std::endl

#include "NodeCoarseAverage.H"
#include "DataIterator.H"
#include "IntVectSet.H"
#include "BoxIterator.H"
#include "NodeSetOperations.H"
#include "parstream.H"
#include "NodeAverageF_F.H"
#include "NamespaceHeader.H"
using std::cout;
using std::endl;

// ---------------------------------------------------------
NodeCoarseAverage::NodeCoarseAverage()
  :
  is_defined(false),m_verbose(false)
{
}


// ---------------------------------------------------------
NodeCoarseAverage::~NodeCoarseAverage()
{
}


// ---------------------------------------------------------
NodeCoarseAverage::NodeCoarseAverage(const DisjointBoxLayout& a_gridsFine,
                                     const DisjointBoxLayout& a_gridsCoarse,
                                     int a_numcomps,
                                     int a_refRatio,
                                     const Box& a_domainFine)
  :
  is_defined(false),m_verbose(false)
{
  ProblemDomain probdomain(a_domainFine);
  define(a_gridsFine, a_gridsCoarse, a_numcomps, a_refRatio, probdomain);
}


// ---------------------------------------------------------
NodeCoarseAverage::NodeCoarseAverage(const DisjointBoxLayout& a_gridsFine,
                                     const DisjointBoxLayout& a_gridsCoarse,
                                     int a_numcomps,
                                     int a_refRatio,
                                     const ProblemDomain& a_domainFine)
  :
  is_defined(false),m_verbose(false)
{
  define(a_gridsFine, a_gridsCoarse, a_numcomps, a_refRatio, a_domainFine);
}


// ---------------------------------------------------------
NodeCoarseAverage::NodeCoarseAverage(const DisjointBoxLayout& a_gridsCoarse,
                                     int a_numcomps,
                                     int a_refRatio,
                                     const Box& a_domainFine)
  :
  is_defined(false),m_verbose(false)
{
  ProblemDomain probdomain(a_domainFine);
  define(a_gridsCoarse, a_numcomps, a_refRatio, probdomain);
}


// ---------------------------------------------------------
NodeCoarseAverage::NodeCoarseAverage(const DisjointBoxLayout& a_gridsCoarse,
                                     int a_numcomps,
                                     int a_refRatio,
                                     const ProblemDomain& a_domainFine)
  :
  is_defined(false),m_verbose(false)
{
  define(a_gridsCoarse, a_numcomps, a_refRatio, a_domainFine);
}


// ---------------------------------------------------------
void
NodeCoarseAverage::define(const DisjointBoxLayout& a_gridsFine,
                          const DisjointBoxLayout& a_gridsCoarse,
                          int a_numcomps,
                          int a_refRatio,
                          const Box& a_domainFine)
{
  ProblemDomain probdomain(a_domainFine);
  define(a_gridsFine, a_gridsCoarse, a_numcomps, a_refRatio, probdomain);
}



// ---------------------------------------------------------
void
NodeCoarseAverage::define(const DisjointBoxLayout& a_gridsFine,
                          const DisjointBoxLayout& a_gridsCoarse,
                          int a_numcomps,
                          int a_refRatio,
                          const ProblemDomain& a_domainFine)
{
  m_refRatio = a_refRatio;
  m_numcomps = a_numcomps;
  m_domainCoarse = coarsen(a_domainFine, m_refRatio);

  // DisjointBoxLayout::coarsen()
  coarsen(m_coarsenedGrids, a_gridsFine, m_refRatio);
  m_coarsenedFine.define(m_coarsenedGrids, m_numcomps);

  IntVect extent = (a_refRatio/2) * IntVect::Unit;
  m_refbox = Box(-extent, extent);
  m_weights.define(m_refbox, 1);
  FORT_NODEAVERAGE_GETWEIGHTS(CHF_FRA(m_weights),
                              CHF_CONST_INT(m_refRatio));

  if (m_verbose)
    cout << "IBN NodeCoarseAverage on "
         << m_domainCoarse.domainBox().size()
         << ", from "
         << m_coarsenedGrids.size() << " coarsened grids to "
         << a_gridsCoarse.size() << " coarser-level grids on "
         << m_domainCoarse.domainBox().size()
         << endl;

  interiorBoundaryNodes(m_IVSV, a_gridsCoarse, m_coarsenedGrids, m_domainCoarse);

  if (m_verbose)
    cout << "IBN NodeCoarseAverage on "
         << m_domainCoarse.domainBox().size()
         << ", from "
         << m_coarsenedGrids.size() << " coarsened grids to selves on "
         << m_domainCoarse.domainBox().size()
         << endl;

  interiorBoundaryNodes(m_IVSVsame, m_coarsenedGrids, m_domainCoarse);

  // m_IVSVfull added by petermc, 21 Jul 2003
  fullIntVectSets(m_IVSVfull, m_IVSVsame);

  m_sameGrids = false;
  is_defined = true;
}



// ---------------------------------------------------------
void
NodeCoarseAverage::define(const DisjointBoxLayout& a_gridsCoarse,
                          int a_numcomps,
                          int a_refRatio,
                          const Box& a_domainFine)
{
  ProblemDomain probdomain(a_domainFine);
  define(a_gridsCoarse, a_numcomps, a_refRatio, probdomain);
}


// ---------------------------------------------------------
void
NodeCoarseAverage::define(const DisjointBoxLayout& a_gridsCoarse,
                          int a_numcomps,
                          int a_refRatio,
                          const ProblemDomain& a_domainFine)
{
  m_refRatio = a_refRatio;
  m_numcomps = a_numcomps;
  m_domainCoarse = coarsen(a_domainFine, m_refRatio);

  // DisjointBoxLayout::coarsen()
  // coarsen(m_coarsenedGrids, a_gridsFine, m_refRatio);
  // m_coarsenedFine.define(m_coarsenedGrids, m_numcomps);
  // 11 Apr 2003:  m_coarsenedGrids == a_gridsCoarse when there
  // is no fine-grid argument.
  m_coarsenedGrids.deepCopy(a_gridsCoarse);
  m_coarsenedGrids.close();
  m_coarsenedFine.define(m_coarsenedGrids, m_numcomps);

  IntVect extent = (a_refRatio/2) * IntVect::Unit;
  m_refbox = Box(-extent, extent);
  m_weights.define(m_refbox, 1);
  FORT_NODEAVERAGE_GETWEIGHTS(CHF_FRA(m_weights),
                              CHF_CONST_INT(m_refRatio));

  if (m_verbose)
    cout << "IBN NodeCoarseAverage on "
         << m_domainCoarse.domainBox().size()
         << ", from "
         << m_coarsenedGrids.size() << " coarsened grids to selves on "
         << m_domainCoarse.domainBox().size()
         << endl;

  // interiorBoundaryNodes(m_IVSV, a_gridsCoarse, m_coarsenedGrids, m_domainCoarse);
  interiorBoundaryNodes(m_IVSVsame, m_coarsenedGrids, m_domainCoarse);

  // m_IVSVfull added by petermc, 21 Jul 2003
  fullIntVectSets(m_IVSVfull, m_IVSVsame);

  // Don't need to use m_IVSV.

  m_sameGrids = true;
  is_defined = true;
}


// ---------------------------------------------------------
bool
NodeCoarseAverage::isDefined() const
{
  return ( is_defined );
}


// ---------------------------------------------------------
void
NodeCoarseAverage::setVerbose( bool a_verbose )
{
  m_verbose = a_verbose ;
}

// ---------------------------------------------------------
void
NodeCoarseAverage::averageToCoarse(LevelData<NodeFArrayBox>& a_coarse,
                                   LevelData<NodeFArrayBox>& a_fine)
{
  // this is the method that is called externally.
  // a_fine must have ghost nodes to depth m_refRatio/2.
  // ghost data is required when averaging onto a fine-fine interface.

  CH_assert(is_defined);
  CH_assert(a_fine.ghostVect() >= (m_refRatio / 2) * IntVect::Unit);

  //CH_assert(a_domainCoarse == m_domainCoarse);

  a_fine.exchange(a_fine.interval());

  DataIterator ditc = m_coarsenedGrids.dataIterator();

  // initialize m_coarsenedFine with zero (on boxes in this proc)
  // if the grids are different.
  if (! m_sameGrids)
    for (ditc.begin(); ditc.ok(); ++ditc)
      m_coarsenedFine[ditc()].getFab().setVal(0.0);

  // set m_coarsenedFine on interior nodes of m_coarsenedGrids to
  // weighted average of fine nodes in neighborhood.
  for (ditc.begin(); ditc.ok(); ++ditc)
    {
      // coarsened_fine, fine are cell-centered but represent data on nodes.
      BaseFab<Real>& coarseFab = (m_sameGrids) ?
        a_coarse[ditc()].getFab() : m_coarsenedFine[ditc()].getFab();

      const BaseFab<Real>& fineFab = a_fine[ditc()].getFab();

      // idea for setting m_coarsenedFine on interior nodes:
      // (1) set m_coarsenedFine in interior of this box
      //     (interior of this box is another box).
      // (2) then set m_coarsenedFine at other interior points.
      // ---

      // (1)
      // set m_coarsenedFine at inner nodes of box
      // (note coarsenedFineFab.box(), not m_coarsenedFine.box(),
      // because we use nodes)
      const Box thisBox = m_coarsenedGrids.get(ditc());

      // inner:  cell-centered box, indices correspond to inner nodes
      Box inner(thisBox);
      for (int idir = 0; idir < CH_SPACEDIM; idir++)
        inner.growLo(idir, -1);

      if (!inner.isEmpty())
        FORT_NODEAVERAGE(CHF_FRA(coarseFab),
                         CHF_CONST_FRA(fineFab),
                         CHF_BOX(inner),
                         CHF_CONST_INT(m_refRatio),
                         CHF_CONST_FRA(m_weights));

      // (2)
      // remove inner nodes from interior, leaving only grid-boundary nodes.
      // cell indices of innerCells are node indices of interior nodes.

      // now do averaging for each interior boundary node.
      Vector<IntVectSet>& IVSvec = m_IVSVsame[ditc()];
      BitSet& fullvec = m_IVSVfull[ditc()];

      for (int lcomp = 0; lcomp < IVSvec.size(); lcomp++)
        {
          IntVectSet& IVS = IVSvec[lcomp];
          if (fullvec[lcomp])
            {
              Box container(IVS.minBox());
              FORT_NODEAVERAGE(CHF_FRA(coarseFab),
                               CHF_CONST_FRA(fineFab),
                               CHF_BOX(container),
                               CHF_CONST_INT(m_refRatio),
                               CHF_CONST_FRA(m_weights));
            }
          else
            for (IVSIterator it(IVS); it.ok(); ++it)
              {
                IntVect iv = it();
                // Set coarseFab(iv, :) to average of fineFab(iv+m_refbox, :)
                // Use either Fortran or C++.
                // petermc, 23 Jul 2003:  Fortran method is found to be faster,
                // by a factor of 2.15 for total time on averageToCoarse
                // function.
                // begin Fortran method
                FORT_NODEAVERAGEPOINT(CHF_FRA(coarseFab),
                                      CHF_CONST_FRA(fineFab),
                                      CHF_CONST_INTVECT(iv),
                                      CHF_CONST_INT(m_refRatio),
                                      CHF_CONST_FRA(m_weights));
                // end Fortran method
                // begin C++ method
//                  for (int ivar = 0; ivar < m_numcomps; ivar++)
//                    {
//                      Real csum = 0.;
//                      for (BoxIterator offb(m_refbox); offb.ok(); ++offb)
//                        {
//                          IntVect ivOff = offb();
//                          IntVect ivFine = m_refRatio*iv + ivOff;
//                          csum += m_weights(ivOff, 0) * fineFab(ivFine, ivar);
//                        }
//                      coarseFab(iv, ivar) = csum;
                // end C++ method
              }
        }
    }

  // copy data at interior nodes of m_coarsenedFine to a_coarse data:
  // see NodeFArrayBox for code.

  // IF the grids are not the same, then
  // copy data at interior nodes of m_coarsenedFine to a_coarse:
  // see NodeFArrayBox for code.
  if (! m_sameGrids)
    {
      copyInteriorNodes(a_coarse, m_coarsenedFine, m_IVSV);
    }

  // copyInteriorNodes(a_coarse, m_coarsenedFine, m_domainCoarse);
}

#include "NamespaceFooter.H"
