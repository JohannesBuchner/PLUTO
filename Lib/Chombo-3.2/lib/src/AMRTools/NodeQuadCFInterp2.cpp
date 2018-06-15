#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeQuadCFInterp2.cpp
// original NodeQuadCFInterp by petermc, Mon, Apr 23, 2001
// petermc, 31 Oct 2001
#include "NodeQuadCFInterp2.H"
#include "NamespaceHeader.H"
using std::cerr;
using std::endl;

// ---------------------------------------------------------
void
NodeQuadCFInterp2::setDefaultValues()
{
  m_isDefined = false;
  m_verbose = false;
}


// ---------------------------------------------------------
NodeQuadCFInterp2::NodeQuadCFInterp2()
{
  setDefaultValues();
}


// ---------------------------------------------------------
NodeQuadCFInterp2::NodeQuadCFInterp2(const DisjointBoxLayout& a_grids,
                                     const Box& a_domain,
                                     const LayoutData<NodeCFIVS>* const a_loCFIVS,
                                     const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                                     bool a_interfaceOnly,
                                     int a_interpolationDegree,
                                     int a_ncomp,
                                     bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  setDefaultValues();
  define(a_grids, probdomain, a_loCFIVS, a_hiCFIVS,
         a_interfaceOnly, a_interpolationDegree, a_ncomp, a_verbose);
}


// ---------------------------------------------------------
NodeQuadCFInterp2::NodeQuadCFInterp2(const DisjointBoxLayout& a_grids,
                                     const ProblemDomain& a_domain,
                                     const LayoutData<NodeCFIVS>* const a_loCFIVS,
                                     const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                                     bool a_interfaceOnly,
                                     int a_interpolationDegree,
                                     int a_ncomp,
                                     bool a_verbose)
{
  setDefaultValues();
  define(a_grids, a_domain, a_loCFIVS, a_hiCFIVS,
         a_interfaceOnly, a_interpolationDegree, a_ncomp, a_verbose);
}


// ---------------------------------------------------------
NodeQuadCFInterp2::~NodeQuadCFInterp2()
{
  clearMemory();
}



// ---------------------------------------------------------
void
NodeQuadCFInterp2::define(const DisjointBoxLayout& a_grids,
                          const Box& a_domain,
                          const LayoutData<NodeCFIVS>* const a_loCFIVS,
                          const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                          bool a_interfaceOnly,
                          int a_interpolationDegree,
                          int a_ncomp,
                          bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, probdomain, a_loCFIVS, a_hiCFIVS,
         a_interfaceOnly, a_interpolationDegree, a_ncomp, a_verbose);
}


// ---------------------------------------------------------
void
NodeQuadCFInterp2::define(const DisjointBoxLayout& a_grids,
                          const ProblemDomain& a_domain,
                          const LayoutData<NodeCFIVS>* const a_loCFIVS,
                          const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                          bool a_interfaceOnly,
                          int a_interpolationDegree,
                          int a_ncomp,
                          bool a_verbose)
{
  m_ncomp = a_ncomp;
  m_grids = a_grids;
  m_domain = a_domain;
  m_domainCoarseNodes = surroundingNodes(coarsen(m_domain.domainBox(), 2));
  m_interfaceOnly = a_interfaceOnly;
  m_verbose = a_verbose;

  // Define objects containing the nodes of the coarse/fine interfaces.
  m_loCFIVS = a_loCFIVS;
  m_hiCFIVS = a_hiCFIVS;

  DisjointBoxLayout coarsenedGrids;
  coarsen(coarsenedGrids, m_grids, 2);

  // ghost layer in m_coarseCopy so that interpolation can also use
  // nodes that are just off the interface.
  m_coarseCopy.define(coarsenedGrids, m_ncomp, IntVect::Unit);

  CH_assert((a_interpolationDegree == 1) || (a_interpolationDegree == 2));
  m_interpolationDegree = a_interpolationDegree;

  m_isDefined = true;
}


// ---------------------------------------------------------
bool
NodeQuadCFInterp2::isDefined() const
{
  return(m_isDefined);
}


// ---------------------------------------------------------
void
NodeQuadCFInterp2::clearMemory()
{
  m_isDefined = false;
}



// ---------------------------------------------------------
void
NodeQuadCFInterp2::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
}



// ---------------------------------------------------------
void
NodeQuadCFInterp2::coarseFineInterp(LevelData<NodeFArrayBox>& a_phiFine,
                                    const LevelData<NodeFArrayBox>& a_phiCoarse)
{
  CH_assert(isDefined());

  // Copy coarse phi to m_coarseCopy.
  // petermc, 15 Nov 2002:  m_coarseCopy is built on the coarsened fine
  // grids, so its layout has no cells that are outside the layout of
  // a_phiCoarse.  Therefore, copyTo does not need Copier.
  // By proper nesting, all ghost nodes in the domain will contain data.
  // (Not completely sure that copyTo is OK with the ghost layer.)
  a_phiCoarse.copyTo(a_phiCoarse.interval(), m_coarseCopy, m_coarseCopy.interval());

  // Now interpolate a_phiFine from m_coarseCopy.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      // a_phiFine lives on m_grids.
      FArrayBox& fineFab = a_phiFine[dit()].getFab();
      // m_coarseCopy lives on m_grids coarsened by 2.
      FArrayBox& psiFab = m_coarseCopy[dit()].getFab();

      // petermc, 7 Jun 2001, added "& m_domainCoarseNodes" so you're not
      // interpolating from nodes outside the domain.
      // (psiFabBox is used to determine getLo, getHi)
      Box psiFabBox = psiFab.box() & m_domainCoarseNodes;

      // psiFabProperBox is the node-centered box of psiFab, without ghosts.
      Box psiFabProperBox(m_grids.get(dit()));
      psiFabProperBox.coarsen(2);
      psiFabProperBox.surroundingNodes();

      // petermc, 30 Jan 2002, added psiFabNodes,
      // to contain the indices of nodes for which we have data in psiFab.
      // (This ensures that we are interpolating only from nodes
      // on the interface IF we are in the m_interfaceOnly case.)
      IntVectSet psiFabNodes;

      // psiFabBox is node-centered, but recall that an
      // IntVectSet is supposed to contain only cells.
      // So we start by shifting psiFabBox to cells with
      // the same indices.
      Box psiFabBoxIndices(psiFabBox);
      psiFabBoxIndices.shiftHalf(IntVect::Unit);
      if (! m_interfaceOnly)
        {
          // Set psiFabNodes to contain the indices of all nodes
          // of psiFabBox, because they all contain data from which
          // we can interpolate.  (Not true if m_interfaceOnly == true.)
          psiFabNodes = IntVectSet(psiFabBoxIndices);
        }

      for (int idir = 0; idir < SpaceDim; idir++)
        {
#if (CH_SPACEDIM ==2)
          int idirOther = 1 - idir;
          // IntVect eOther = BASISV(idirOther);
#else
          int i0 = 0, i1 = 2;
          if (idir == 0) i0 = 1;
          if (idir == 2) i1 = 1;
          IntVect e0 = BASISV(i0);
          IntVect e1 = BASISV(i1);
#endif

          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide mySide = sit();

              const NodeCFIVS* cfivs_ptr = NULL;
              switch (mySide)
                {
                case Side::Lo:
                  {
                    cfivs_ptr = &m_loCFIVS[idir][dit()];
                    break;
                  }
                case Side::Hi:
                  {
                    cfivs_ptr = &m_hiCFIVS[idir][dit()];
                    break;
                  }
                default:
                  cerr << "NodeQuadCFInterp2::coarseFineInterp(): bogus side" << endl;
                  abort();
                }

              // interpIVS is the set of nodes on which we set a_phiFine
              // by interpolating psi.
              const IntVectSet& interpIVS = cfivs_ptr->getFineIVS();

              if (! interpIVS.isEmpty() )
                {
                  // petermc, 30 Jan 2002, added psiFabNodes, to contain
                  // the indices of nodes for which we have data in psiFab.
                  // (This ensures that we are interpolating only from nodes
                  // on the interface IF we are in the m_interfaceOnly case.)
                  if (m_interfaceOnly)
                    {
                      psiFabNodes = IntVectSet(interpIVS);
                      psiFabNodes.coarsen(2);

                      // Add nodes on physical boundary faces.
                      // It is not sufficient to do this only for the
                      // face (idir, mySide).
                      // Example:  on domain [0:31]^3,
                      // psiFabBox = [19:32, 7:13, 7:13]
                      // and interpIVS = [40:63, 16, 16:24]
                      // but psiFabNodes = [20:31, 8, 8:12] by coarsening
                      // omits the nodes [32, 8, 8:12].
                      // Add these nodes here.
                      SideIterator sitFace;
                      for (int idirFace = 0; idirFace < SpaceDim; idirFace++)
                        if (! m_domain.isPeriodic(idirFace))
                          for (sitFace.begin(); sitFace.ok(); sitFace.next())
                            {
                              Side::LoHiSide faceSide = sitFace();
                              // domainCoarseBoundary is
                              // physical boundary on this face
                              Box domainCoarseBoundary(m_domainCoarseNodes);
                              switch (faceSide)
                                {
                                case Side::Lo:
                                  {
                                    int extrm = m_domainCoarseNodes.smallEnd(idirFace);
                                    domainCoarseBoundary.setBig(idirFace, extrm);
                                    break;
                                  }
                                case Side::Hi:
                                  {
                                    int extrm = m_domainCoarseNodes.bigEnd(idirFace);
                                    domainCoarseBoundary.setSmall(idirFace, extrm);
                                    break;
                                  }
                                default:
                                  cerr << "NodeQuadCFInterp2::coarseFineInterp(): bogus side" << endl;
                                  abort();
                                }
                              // psiFabBoxBoundary is the part of
                              // the physical boundary
                              // that coincides with this face of psiFabBox.
                              Box psiFabBoxBoundary(psiFabProperBox);
                              psiFabBoxBoundary &= domainCoarseBoundary;

                              // Put these physical boundary nodes into
                              // psiFabNodes, because we can interpolate
                              // from these nodes.
                              if (! psiFabBoxBoundary.isEmpty())
                                {
                                  Box psiFabBoxBoundaryIndices(psiFabBoxBoundary);
                                  psiFabBoxBoundaryIndices.shiftHalf(IntVect::Unit);
                                  // All we do here is:
                                  // psiFabNodes |= psiFabBoxBoundaryIndices;
                                  // But remember that |= is inefficient.
                                  // Implement A |= B as
                                  // A = U - ((U - A) & (U - B))
                                  // where U will include both A and B.
                                  IntVectSet psiFabNodesComp(psiFabBoxIndices);
                                  psiFabNodesComp -= psiFabNodes;
                                  IntVectSet psiFabBoundaryComp(psiFabBoxIndices);
                                  psiFabBoundaryComp -= psiFabBoxBoundaryIndices;
                                  psiFabNodesComp &= psiFabBoundaryComp;
                                  psiFabNodes = IntVectSet(psiFabBoxIndices);
                                  psiFabNodes -= psiFabNodesComp;
                                }
                            } // loop over faces
                    } // if (m_interfaceOnly)

                  IVSIterator ivsit(interpIVS);
                  for (ivsit.begin(); ivsit.ok(); ++ivsit)
                    {
                      IntVect iv = ivsit();
#if (CH_SPACEDIM ==2)
                      interpLine(fineFab, psiFab, psiFabNodes, iv, idirOther);
#else
                      IntVect ivcBase = iv / 2;

                      int ii[2];
                      ii[0] = i0; ii[1] = i1;
                      int parity[2]; // elements are 0 or 1
                      for (int iway = 0; iway < 2; iway++)
                        {
                          int iiway = ii[iway];
                          parity[iway] = iv[iiway] - 2 * ivcBase[iiway];
                        }

                      if (parity[0] == 0 && parity[1] == 0)
                        {
                          // Merely project, iv == 2 * ivcBase.
                          for (int ivar = 0; ivar < m_ncomp; ivar++)
                            fineFab(iv, ivar) = psiFab(ivcBase, ivar);
                        }
                      else if (parity[0] == 0)
                        {
                          interpLine(fineFab, psiFab, psiFabNodes, iv, i1);
                        }
                      else if (parity[1] == 0)
                        {
                          interpLine(fineFab, psiFab, psiFabNodes, iv, i0);
                        }
                      else // parity[0] == 1 == parity[1]
                        {
                          // Fine node is in centre of a coarse 2-D cell.

                          // Initialize with mean of nearest coarse nodes.
                          for (int ivar = 0; ivar < m_ncomp; ivar++)
                            fineFab(iv, ivar) = 0.25 *
                              (psiFab(ivcBase, ivar) +
                               psiFab(ivcBase + e0, ivar) +
                               psiFab(ivcBase + e1, ivar) +
                               psiFab(ivcBase + e0 + e1, ivar));

                          if (m_interpolationDegree >= 2)
                            {
                              // Subtract 1/8 * sum of second derivatives i0, i1.

                              // Since the box should have length at least 2 in
                              // each dimension, you should be able to find
                              // second derivatives in each dimension.
                              for (int irun = 0; irun < 2; irun++)
                                {
                                  // We take derivatives in ii[irun] direction.
                                  int iThis = ii[irun];
                                  int iOther = ii[1-irun];
                                  IntVect e = BASISV(iThis); // e0 or e1
                                  IntVect eOther = BASISV(iOther); // e1 or e0

                                  // Example with iThis == 1, iOther == 0
                                  //
                                  //    O   O   ivcBasePerp + 2*e  (ipar == 1)
                                  //    |   |
                                  //    O---O
                                  //    | x |
                                  //    O---O   ivcBasePerp
                                  //    |   |
                                  //    O   O   ivcBasePerp - e    (ipar == 0)
                                  //
                                  //    0   1  <- iperp

                                  // petermc, 26 June 2002:
                                  // Note that psiFabNodes may contain
                                  // ivcBase+2*e but not ivcBase+2*e+eOther.
                                  // I didn't realize this before, and put the
                                  // gotNext check _before_ the loop over iperp.
                                  // Now fixed.
                                  for (int iperp = 0; iperp < 2; iperp++)
                                    {
                                      IntVect ivcBasePerp = ivcBase + iperp * eOther;

                                      bool gotNext[2];
                                      for (int ipar = 0; ipar < 2; ipar++)
                                        gotNext[ipar] =
                                          psiFabNodes.contains(ivcBasePerp +
                                                               (3*ipar-1) * e);

                                      if (gotNext[0] && gotNext[1])
                                        {
                                          // O---O-x-O---O
                                          for (int ivar = 0; ivar < m_ncomp; ivar++)
                                            fineFab(iv, ivar) -= 1./32. *
                                              (psiFab(ivcBasePerp - e, ivar) -
                                               psiFab(ivcBasePerp, ivar) -
                                               psiFab(ivcBasePerp + e, ivar) +
                                               psiFab(ivcBasePerp + e + e, ivar));
                                        }
                                      else
                                        for (int ipar = 0; ipar < 2; ipar++)
                                          if (gotNext[ipar])
                                            {
                                              IntVect ivcBaseSide =
                                                ivcBasePerp + ipar * e;
                                              // find derivative on one side:
                                              //
                                              // O---O-x-O   O    ipar == 0
                                              //
                                              // O   O-x-O---O    ipar == 1
                                              for (int ivar = 0; ivar < m_ncomp; ivar++)
                                                fineFab(iv, ivar) -= 1./16. *
                                                  (psiFab(ivcBaseSide - e, ivar) -
                                                   2. * psiFab(ivcBaseSide, ivar) +
                                                   psiFab(ivcBaseSide + e, ivar));
                                            }
                                    } // loop over iperp
                                } // loop over irun (derivative of which variable)
                            } // if interpolation degree >= 2
                        } // parity[0] == 1 == parity[1]
#endif
                    } // loop over interpIVS points
                } // if (! interpIVS.isEmpty() )
            } // loop over sides
        } // loop over directions
    } // loop over grids
}



// ---------------------------------------------------------
void
NodeQuadCFInterp2::interpLine(FArrayBox& a_fineFab,
                              const FArrayBox& a_psiFab,
                              const IntVectSet& a_psiFabNodes,
                              const IntVect& a_iv,
                              int a_idirOther)
{
  IntVect eOther = BASISV(a_idirOther);

  IntVect ivc[2];
  for (int iside = 0; iside < 2; iside++)
    ivc[iside] = (a_iv + iside * eOther) / 2;

  int parity = a_iv[a_idirOther] - 2 * ivc[0][a_idirOther];
  if (parity == 0)
    {
      // merely project, ivc[0] == ivc[1] == iv / 2.
      for (int ivar = 0; ivar < m_ncomp; ivar++)
        a_fineFab(a_iv, ivar) = a_psiFab(ivc[0], ivar);
    }
  else // parity == 1, ivc[0] != ivc[1]
    {
      // Initialize with mean of nearest coarse nodes.
      for (int ivar = 0; ivar < m_ncomp; ivar++)
        a_fineFab(a_iv, ivar) = 0.5 *
          (a_psiFab(ivc[0], ivar) + a_psiFab(ivc[1], ivar));

      if (m_interpolationDegree >= 2)
        {
          IntVect ivc2[2];
          bool gotc2[2];
          for (int iside = 0; iside < 2; iside++)
            {
              // ivc2[0] == ivc[0] - eOther
              // ivc2[1] == ivc[1] + eOther
              ivc2[iside] = ivc[iside] + (2*iside-1) * eOther;
              gotc2[iside] = a_psiFabNodes.contains(ivc2[iside]);
            }

          // Now subtract 1/8 * second derivative.
          if (gotc2[0] && gotc2[1])
            {
              // O---O-x-O---O
              for (int ivar = 0; ivar < m_ncomp; ivar++)
                a_fineFab(a_iv, ivar) -= 1./16. *
                  (a_psiFab(ivc2[0], ivar) -
                   a_psiFab(ivc[0], ivar) -
                   a_psiFab(ivc[1], ivar) +
                   a_psiFab(ivc2[1], ivar));
            }
          else
            {
              // check each side
              for (int iside = 0; iside < 2; iside++)
                if (gotc2[iside])
                  {
                    // gotc2[1-iside] is false because if we are here
                    // then gotc2[0], gotc2[1] are not both true.
                    // iside == 0:   O---O-x-O   O
                    // iside == 1:   O   O-x-O---O
                    // second derivative from one side
                    for (int ivar = 0; ivar < m_ncomp; ivar++)
                      a_fineFab(a_iv, ivar) -= 1./8. *
                        (a_psiFab(ivc2[iside], ivar) -
                         2. * a_psiFab(ivc[iside], ivar) +
                         a_psiFab(ivc[1-iside], ivar));
                  }
            }
        } // if interpolation degree >= 2
    }
}

#include "NamespaceFooter.H"
