#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRNodeOp.H"
#include "FORT_PROTO.H"
#include "AMRPoissonOpF_F.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "NodeSetOperations.H"
#include "NodeNorms.H"
#include "NodeSetOperations.H"
#include "NamespaceHeader.H"
#include "NodeLevelMGF_F.H"

void AMRNodeOp::
projectFineInterior(LevelData<NodeFArrayBox>&       a_phi,
                    const LevelData<NodeFArrayBox>& a_phiFine)
{
  int refToFine =  m_refToFiner;
  CH_assert(a_phi.nComp() == 1);

  // coarsen a_phiFine, call it phiFineCoarsened
  LevelData<NodeFArrayBox> phiFineCoarsened(m_coarsenedFineGrids, 1, IntVect::Zero);

  for (DataIterator dit = m_coarsenedFineGrids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phiFineCoarsenedFab = phiFineCoarsened[dit()].getFab();
      const FArrayBox& phiFineFab = a_phiFine[dit()].getFab();

      const Box coarsenedNodes = surroundingNodes(m_coarsenedFineGrids[dit()]);

      for (BoxIterator bit(coarsenedNodes); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          phiFineCoarsenedFab(iv, 0) = phiFineFab(refToFine*iv, 0);
        }
    }

  // Copy phiFineCoarsened to a_phi on interior nodes of coarsenedFineGrids.
  // Code is in NodeFArrayBox.

  copyInteriorNodes(a_phi, phiFineCoarsened, m_IVSVcoarsenedFine);
}
void AMRNodeOp::
createCoarsened(LevelData<NodeFArrayBox>&       a_lhs,
                const LevelData<NodeFArrayBox>& a_rhs,
                const int &                 a_refRat)
{
  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  DisjointBoxLayout dbl = a_rhs.disjointBoxLayout();
  CH_assert(dbl.coarsenable(a_refRat));

  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, dbl, a_refRat);

  a_lhs.define(dblCoarsenedFine, ncomp, a_rhs.ghostVect());
}
  ///
AMRNodeOp::AMRNodeOp()
{
  m_hasFiner = false;
}

/** full define function for AMRLevelOp with both coarser and finer levels */
void AMRNodeOp::define(const DisjointBoxLayout& a_grids,
                       const DisjointBoxLayout& a_gridsFiner,
                       const DisjointBoxLayout& a_gridsCoarser,
                       const Real&              a_dxLevel,
                       int                      a_refRatio,
                       int                      a_refRatioFiner,
                       const ProblemDomain&     a_domain,
                       NodeBCFunc               a_bc)
{
  m_hasFiner = true;
  m_refToFiner = a_refRatioFiner;
  coarsen(m_coarsenedFineGrids, a_gridsFiner, a_refRatioFiner);

  this->define(a_grids, a_gridsCoarser, a_dxLevel, a_refRatio, a_domain, a_bc);
  //have to do this again because the above function overwrites
  m_refToFiner = a_refRatioFiner;
}

/** full define function for AMRLevelOp<LevelData<NodeFArrayBox> > with finer levels, but no coarser */
void AMRNodeOp::define(const DisjointBoxLayout& a_grids,
                       const DisjointBoxLayout& a_gridsFiner,
                       const Real&              a_dxLevel,
                       int                      a_refRatio, // dummy argument
                       int                      a_refRatioFiner,
                       const ProblemDomain&     a_domain,
                       NodeBCFunc               a_bc)
{
  CH_assert(a_refRatio == 1);
  m_hasFiner = true;
  coarsen(m_coarsenedFineGrids, a_gridsFiner, a_refRatioFiner);
  this->define(a_grids, a_dxLevel, a_domain, a_bc); //calls the MG version of define
  m_refToFiner = a_refRatioFiner;
}

void AMRNodeOp::define(const DisjointBoxLayout& a_grids,
                       const DisjointBoxLayout& a_coarse,
                       const Real&              a_dxLevel,
                       int                      a_refRatio,
                       const ProblemDomain&     a_domain,
                       NodeBCFunc               a_bc)
{
  m_refToCoarser = a_refRatio;
  m_dxCrse = a_refRatio*a_dxLevel;
  m_refToFiner = 1;
  m_bc     = a_bc;
  m_domain = a_domain;
  m_domainInteriorNodes = surroundingNodes(m_domain.domainBox());
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_domain.isPeriodic(idir))
        {
          m_domainInteriorNodes.setPeriodic(idir, true);
        }
      else
        {
          m_domainInteriorNodes.setPeriodic(idir, false);
          m_domainInteriorNodes.grow(idir, -1);
        }
    }
  m_dx     = a_dxLevel;
  m_dxCrse = a_refRatio*a_dxLevel;

  //these get set again after define is called
  m_alpha = 0.0;
  m_beta  = 1.0;

  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  setCFIVS(a_grids);
  m_interpWithCoarser.define(a_grids, a_dxLevel, m_domain, m_loCFIVS, m_hiCFIVS, m_refToCoarser, m_bc);
  interiorBoundaryNodes(m_IVSV, a_grids, m_domain);
  exteriorBoundaryNodes(m_IVSVext, m_IVSV, a_grids);
  ProblemDomain coarseDomain = coarsen(m_domain, a_refRatio);
  interiorBoundaryNodes(m_IVSVcoarsened, a_grids, a_coarse, coarseDomain);
  if (m_hasFiner)
    {
      interiorBoundaryNodes(m_IVSVcoarsenedFine, a_grids, m_coarsenedFineGrids, m_domain);
    }

  if (a_grids.coarsenable(2))
    {
      DisjointBoxLayout mgGrids;
      coarsen(mgGrids, a_grids, 2);
      //m_averageOpMG.define(a_grids, mgGrids,  1, 2, m_domain);
      // m_coarsenedGrids is just coarsened version of m_grids, so use
      // the alternative define
      int ref = 2;
      int ncomp = 1;
      m_averageOpMG.define(mgGrids,  ncomp, ref, m_domain);
    }
}

void
AMRNodeOp::setCFIVS(const DisjointBoxLayout& a_grids)
{

  DataIterator dit = a_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(a_grids);
      m_hiCFIVS[idir].define(a_grids);
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box bx(a_grids.get(dit()));
          m_loCFIVS[idir][dit()].define(m_domain, bx, a_grids, idir, Side::Lo);
          m_hiCFIVS[idir][dit()].define(m_domain, bx, a_grids, idir, Side::Hi);
        }
    }

  interiorBoundaryNodes(m_IVSV, a_grids, m_domain);
  fullIntVectSets(m_IVSVfull, m_IVSV);
  exteriorBoundaryNodes(m_IVSVext, m_IVSV, a_grids);
}

void AMRNodeOp::define(const DisjointBoxLayout& a_grids,
                       const Real&              a_dx,
                       const ProblemDomain&     a_domain,
                       NodeBCFunc               a_bc)
{
  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
  m_dxCrse = 2*a_dx;
  m_refToCoarser = 2; // redefined in AMRLevelOp<LevelData<NodeFArrayBox> >::define virtual function.
  m_refToFiner   = 2;

  //these get set again after define is called
  m_alpha = 0.0;
  m_beta  = 1.0;

  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      LayoutData<NodeCFIVS>& lo =  m_loCFIVS[i];
      LayoutData<NodeCFIVS>& hi =  m_hiCFIVS[i];
      lo.define(a_grids);
      hi.define(a_grids);
      for (DataIterator dit(a_grids); dit.ok(); ++dit)
        {
          lo[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Lo);
          hi[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Hi);
        }
    }
  m_domainInteriorNodes = surroundingNodes(m_domain.domainBox());
  for (int idir = 0; idir < SpaceDim; idir++)
    if (m_domain.isPeriodic(idir))
      m_domainInteriorNodes.setPeriodic(idir, true);
    else
      {
        m_domainInteriorNodes.setPeriodic(idir, false);
        m_domainInteriorNodes.grow(idir, -1);
      }
  interiorBoundaryNodes(m_IVSV, a_grids, m_domain);
  fullIntVectSets(m_IVSVfull, m_IVSV);
  exteriorBoundaryNodes(m_IVSVext, m_IVSV, a_grids);

  if (m_hasFiner)
    {
      interiorBoundaryNodes(m_IVSVcoarsenedFine, a_grids, m_coarsenedFineGrids, m_domain);
    }

  if (a_grids.coarsenable(2))
    {
      DisjointBoxLayout mgGrids;
      coarsen(mgGrids, a_grids, 2);
      //m_averageOpMG.define(a_grids, mgGrids,  1, 2, m_domain);
      int ref = 2;
      int ncomp = 1;
      m_averageOpMG.define(mgGrids,  ncomp, ref, m_domain);
    }
}

void AMRNodeOp::residual(LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_phi,
                         const LevelData<NodeFArrayBox>& a_rhs, bool a_homogeneous)
{
  m_levelOps.setToZero(a_lhs);
  applyOp(a_lhs, a_phi, a_homogeneous);
  incr(a_lhs, a_rhs, -1);
  scale(a_lhs, -1.0);

  zeroBoundaryNodes(a_lhs, m_IVSVext);
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void AMRNodeOp::preCond(LevelData<NodeFArrayBox>&       a_phi,
                        const LevelData<NodeFArrayBox>& a_rhs)
{
  Real mult = - (m_dx*m_dx) / (2.0*SpaceDim);
  Interval comps = a_phi.interval();

  // Use copy() now, more efficient than copyTo().
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phiFab = a_phi[dit()].getFab();
      const FArrayBox& rhsFab = a_rhs[dit()].getFab();
      Box bx(rhsFab.box());

      phiFab.copy(bx, comps, bx, rhsFab, comps);

      phiFab *= mult;
    }
  int nrelax = 4;
  relax(a_phi, a_rhs, nrelax);
}

void AMRNodeOp::applyOp(LevelData<NodeFArrayBox>&       a_LofPhi,
                        const LevelData<NodeFArrayBox>& a_phi,
                        bool a_homogeneous )
{
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  LevelData<NodeFArrayBox>& phi = (LevelData<NodeFArrayBox> &) a_phi;
  // homogeneous physical boundary conditions
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      m_bc(phi[dit()],  dbl[dit()], m_domain,  m_dx, a_homogeneous);
    }
  phi.exchange(a_phi.interval(), m_exchangeCopier);

  applyOpOnly(a_LofPhi, phi);
}

void AMRNodeOp::applyOpOnly(LevelData<NodeFArrayBox>&       a_LofPhi,
                            const LevelData<NodeFArrayBox>& a_phi)
{
  // Apply operator to a_phi _at interior nodes_.
  // See NodePoissonOp.ChF for NODEOPLAP().

  CH_assert(a_phi.ghostVect() >= IntVect::Unit); // need this for stencil
  Box alldomainbox = m_domain.domainBox();
  DisjointBoxLayout grids = a_phi.disjointBoxLayout();
  // added by petermc, 21 Jul 2003, used in C++ method below.
  // Real dxinv2 = 1. / (m_dx * m_dx);
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      const Box& thisBox = grids.get(dit());

      const BaseFab<Real>& thisPhiFab = a_phi[dit()].getFab();

      BaseFab<Real>& thisLPhiFab = a_LofPhi[dit()].getFab();
      // initialize with 0, because otherwise it can contain garbage
      thisLPhiFab.setVal(0.);

      // (1)
      // set thisLPhi at inner NODES of box
      // (thisLPhiFab.box() is NODE-centered)

      // const Box& inner = grow(thisLPhiFab.box(), -1);
      // ... but remember thisLPhiFab has ghost cells now.
      Box boxNodes = surroundingNodes(thisBox);
      // get the enclosed nodes
      boxNodes.grow(-1);

      // Box inner(thisBox);
      // for (int dir = 0; dir < CH_SPACEDIM; dir++)
      // inner.growLo(dir, -1);

      if (! boxNodes.isEmpty() )
        {
          // In the periodic case, you must have copied the edge data
          // to the ghost nodes on the other side.
          FORT_NODEOPLAP(CHF_FRA(thisLPhiFab),
                         CHF_CONST_FRA(thisPhiFab),
                         CHF_BOX(boxNodes),
                         CHF_CONST_REAL(m_dx));
        }

      // (2)
      // remove inner nodes from interior, leaving only grid-boundary nodes.
      // cell indices of innerCells are node indices of interior nodes.

      // new:  actually innerCells are nodes.
      // const Box& innerCells = grow(surroundingNodes(thisBox), -1);

      // now do averaging for each interior boundary node.
      Vector<IntVectSet>& IVSvec = m_IVSV[dit()];
      BitSet& fullvec = m_IVSVfull[dit()];

      for (int lcomp = 0; lcomp < IVSvec.size(); lcomp++)
        {
          IntVectSet& IVS = IVSvec[lcomp];
          if (fullvec[lcomp])
            {
              Box container(IVS.minBox());
              FORT_NODEOPLAP(CHF_FRA(thisLPhiFab),
                             CHF_CONST_FRA(thisPhiFab),
                             CHF_BOX(container),
                             CHF_CONST_REAL(m_dx));
            }
          else
            for (IVSIterator it(IVS); it.ok(); ++it)
              {
                IntVect iv = it();
                // Set thisLPhiFab(iv, 0) to operator on thisPhiFab.
                // Use either Fortran or C++.
                FORT_NODEOPLAPPOINT(CHF_FRA(thisLPhiFab),
                                    CHF_CONST_FRA(thisPhiFab),
                                    CHF_CONST_INTVECT(iv),
                                    CHF_CONST_REAL(m_dx));
              }
        }
    }
}

///compute norm over all cells on coarse not covered by finer
Real
AMRNodeOp::AMRNorm(const LevelData<NodeFArrayBox>& a_coarResid,
                   const LevelData<NodeFArrayBox>& a_fineResid,
                   const int& a_refRat,
                   const int& a_ord)

{
  const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();
  const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();

  //create temp and zero out under finer grids
  LevelData<NodeFArrayBox> coarTemp;
  m_levelOps.create(coarTemp, a_coarResid);
  m_levelOps.setToZero(coarTemp);
  m_levelOps.assign(coarTemp, a_coarResid);
  int ncomp = coarTemp.nComp();
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      NodeFArrayBox& coarTempFAB = coarTemp[dit()];
      LayoutIterator litFine = fineGrids.layoutIterator();
      for (litFine.reset(); litFine.ok(); ++litFine)
        {
          Box overlayBox = coarTempFAB.box();
          Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat);

          overlayBox &= coarsenedGrid;
          Box nodeBox = surroundingNodes(overlayBox);
          if (!overlayBox.isEmpty())
            {
              coarTempFAB.getFab().setVal(0.0,nodeBox,0, ncomp);
            }
        }
    }
  //return norm of temp
  return norm(coarTemp, a_ord);
}
void AMRNodeOp::create(    LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_rhs)
{
  m_levelOps.create(a_lhs, a_rhs);
  m_levelOps.setToZero(a_lhs);
}

void AMRNodeOp::assign(    LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_rhs)
{
  m_levelOps.assign(a_lhs, a_rhs);
}

Real AMRNodeOp::dotProduct(const LevelData<NodeFArrayBox>& a_1, const LevelData<NodeFArrayBox>& a_2)
{
  return m_levelOps.dotProduct(a_1, a_2);
}
void AMRNodeOp::incr( LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_x, Real a_scale)
{
  m_levelOps.incr(a_lhs, a_x, a_scale);
}
void AMRNodeOp::axby( LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_x,
                      const LevelData<NodeFArrayBox>& a_y, Real a_a, Real a_b)
{
  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}
void AMRNodeOp::scale(LevelData<NodeFArrayBox>& a_lhs, const Real& a_scale)
{
  m_levelOps.scale(a_lhs, a_scale);
}
void AMRNodeOp::setToZero(LevelData<NodeFArrayBox>& a_lhs)
{
  m_levelOps.setToZero(a_lhs);
}

void AMRNodeOp::relax(LevelData<NodeFArrayBox>& a_e,
                      const LevelData<NodeFArrayBox>& a_residual,
                      int a_iterations)
{
  for (int i=0; i<a_iterations; i++)
    {
      levelGSRB(a_e, a_residual);
      //levelJacobi(a_e, a_residual);
    }
}

void AMRNodeOp::createCoarser(LevelData<NodeFArrayBox>& a_coarse,
                              const LevelData<NodeFArrayBox>& a_fine,
                              bool ghosted)
{
  // CH_assert(!ghosted);
  IntVect ghost = a_fine.ghostVect();
  DisjointBoxLayout dbl;
  CH_assert(dbl.coarsenable(2));
  coarsen(dbl, a_fine.disjointBoxLayout(), 2); //multigrid, so coarsen by 2
  a_coarse.define(dbl, a_fine.nComp(), ghost);
  m_levelOps.setToZero(a_coarse);
}

void AMRNodeOp::restrictResidual(LevelData<NodeFArrayBox>& a_resCoarse,
                                 LevelData<NodeFArrayBox>& a_phiFine,
                                 const LevelData<NodeFArrayBox>& a_rhsFine)
{
  homogeneousCFInterp(a_phiFine);
  LevelData<NodeFArrayBox> residFine;
  m_levelOps.create(residFine, a_phiFine);
  m_levelOps.setToZero(residFine);
  residual(residFine, a_phiFine, a_rhsFine, true);

  m_averageOpMG.averageToCoarse(a_resCoarse, residFine);
}

Real AMRNodeOp::norm(const LevelData<NodeFArrayBox>& a_x, int a_ord)
{
  Interval comps(0,0);
  Real normVal = CH_XD::norm(a_x, m_domain, m_IVSVext, m_dx, comps, a_ord);
  return normVal;
}
/***/

void AMRNodeOp::prolongIncrement(LevelData<NodeFArrayBox>&       a_phiThisLevel,
                                 const LevelData<NodeFArrayBox>& a_correctCoarse,
                                 int a_refRat)
{
  //assumes (as did the old operator) that a_correctCoarse is over the same layout
  //(just coarsened) as a_phiThisLevel
  const DisjointBoxLayout& gridsFine = a_phiThisLevel.disjointBoxLayout();
  const DisjointBoxLayout& gridsCoar = a_correctCoarse.disjointBoxLayout();
  Box boxRef(IntVect::Zero, (a_refRat-1)*IntVect::Unit);
  Box corners(IntVect::Zero, IntVect::Unit);
  FArrayBox weights(corners, boxRef.numPts());
  FORT_NODEINTERPMG_GETWEIGHTS(CHF_CONST_INT(a_refRat),
                               CHF_BOX(boxRef),
                               CHF_FRA(weights));

  for (DataIterator dit = gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      Box crseBox(gridsCoar.get(dit()));

      const FArrayBox& crseFab = a_correctCoarse[dit()].getFab();
      FArrayBox& fineFab = a_phiThisLevel[dit()].getFab();

      FORT_NODEINTERPMG(CHF_FRA(fineFab),
                        CHF_CONST_FRA(crseFab),
                        CHF_BOX(crseBox),
                        CHF_CONST_INT(a_refRat),
                        CHF_BOX(boxRef),
                        CHF_FRA(weights));
    }
}
void AMRNodeOp::prolongIncrement(LevelData<NodeFArrayBox>&       a_phiThisLevel,
                                 const LevelData<NodeFArrayBox>& a_correctCoarse)
{
  int refRat = 2; //this is the multigrid function so always 2
  prolongIncrement(a_phiThisLevel, a_correctCoarse, refRat);

}
/***/
void AMRNodeOp::
levelGSRB(LevelData<NodeFArrayBox>&       a_phi,
          const LevelData<NodeFArrayBox>& a_rhs)
{

  DataIterator dit = a_phi.dataIterator();

  // do red pass, then black pass
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {
      // do coarse/fine and copy bc's.
      // should be done patch by patch.

      homogeneousCFInterp(a_phi);

      //fill in intersection of ghostcells and a_phi's boxes
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      const DisjointBoxLayout& dbl =a_phi.disjointBoxLayout();
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_bc(a_phi[dit()],  dbl[dit()], m_domain,  m_dx, true);
          Box thisBox(surroundingNodes(dbl.get(dit())));
          // removed by petermc, 10 Mar 2003.  ghosts will be used.
          // brought back by petermc, 13 Mar 2003, with new definition
          // of m_domainInteriorNodes depending on boundary conditions.
          thisBox &= m_domainInteriorNodes;
          if (! thisBox.isEmpty() )
            {
              FArrayBox& phiFab = a_phi[dit()].getFab();
              const FArrayBox& rhsFab = a_rhs[dit()].getFab();
              FORT_NODEGSRBLEVELLAP(CHF_FRA(phiFab),
                                    CHF_CONST_FRA(rhsFab),
                                    CHF_BOX(thisBox),
                                    CHF_CONST_REAL(m_dx),
                                    CHF_CONST_INT(whichPass));
            }
        } // end loop through grids

      // deleted by petermc, 29 Jul 2003:  don't need this now
      // a_phi.exchange(a_phi.interval(), m_exchangeCopier);

    } // end loop through red-black
}

void
AMRNodeOp::homogeneousCFInterp(LevelData<NodeFArrayBox>& a_phif)
{

  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

  DataIterator dit = a_phif.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterp(a_phif,datInd,idir,sit());
            }
        }
    }
}

/*************************/
void
AMRNodeOp::homogeneousCFInterp(LevelData<NodeFArrayBox>& a_phif,
                               const DataIndex& a_datInd,
                               int a_idir,
                               Side::LoHiSide a_hiorlo)
{

  const NodeCFIVS* cfivs_ptr = NULL;
  switch (a_hiorlo)
    {
    case Side::Lo:
      {
        cfivs_ptr = &m_loCFIVS[a_idir][a_datInd];
        break;
      }
    case Side::Hi:
      {
        cfivs_ptr = &m_hiCFIVS[a_idir][a_datInd];
        break;
      }
    default:
      {
        cerr << "NodePoissonOp::homogeneousCFInterp(): bogus side" << endl;
        abort();
      }
    }
  FArrayBox& phiFab = a_phif[a_datInd].getFab();
  // set a_phif to 0 on coarse-fine boundary nodes.
  if (! cfivs_ptr->isEmpty() )
  {
    if (cfivs_ptr->isPacked())
      {
        Box cfivsBox(cfivs_ptr->packedBox());
        // shift from CELL centering to NODE centering
        // because the indices of bx actually represent NODEs.
        cfivsBox.shiftHalf(-IntVect::Unit);
        cfivsBox &= phiFab.box();
        // zero out all m_ncomp components starting with 0
        phiFab.setVal(0., cfivsBox, 0, 1);
      }
    else
      {
        IntVectSet interp_ivs = cfivs_ptr->getFineIVS();
        IVSIterator ivsit(interp_ivs);
        for (ivsit.begin(); ivsit.ok(); ++ivsit)
          {
            phiFab(ivsit(), 0) = 0.;
          }
      }
  }
}

void AMRNodeOp::AMRResidualNF(LevelData<NodeFArrayBox>& a_residual,
                              const LevelData<NodeFArrayBox>& a_phi,
                              const LevelData<NodeFArrayBox>& a_phiCoarse,
                              const LevelData<NodeFArrayBox>& a_rhs,
                              bool a_homogeneousPhysBC)
{
  LevelData<NodeFArrayBox>& phi = (LevelData<NodeFArrayBox>&)a_phi;

  m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse, false);
  this->residual(a_residual, a_phi, a_rhs, a_homogeneousPhysBC ); //apply boundary conditions
}

void AMRNodeOp::AMRResidual(LevelData<NodeFArrayBox>&       a_residual,
                            const LevelData<NodeFArrayBox>& a_phiFine,
                            const LevelData<NodeFArrayBox>& a_phi,
                            const LevelData<NodeFArrayBox>& a_phiCoarse,
                            const LevelData<NodeFArrayBox>& a_rhs,
                            bool                            a_homogeneousPhysBC,
                            AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp)
{
  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse,
              a_homogeneousPhysBC, a_finerOp);

  incr(a_residual, a_rhs, -1.0);
  // residual is rhs - L(phi)
  scale(a_residual, -1.0);
  zeroBoundaryNodes(a_residual, m_IVSVext);
}

void AMRNodeOp::AMRResidualNC(LevelData<NodeFArrayBox>& a_residual,
                              const LevelData<NodeFArrayBox>& a_phiFine,
                              const LevelData<NodeFArrayBox>& a_phi,
                              const LevelData<NodeFArrayBox>& a_rhs,
                              bool a_homogeneousPhysBC,
                              AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp)
{

  AMROperatorNC(a_residual, a_phiFine, a_phi,
                a_homogeneousPhysBC,
                a_finerOp);

  axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
  zeroBoundaryNodes(a_residual, m_IVSVext);
}

void AMRNodeOp::AMROperatorNF(LevelData<NodeFArrayBox>& a_LofPhi,
                              const LevelData<NodeFArrayBox>& a_phi,
                              const LevelData<NodeFArrayBox>& a_phiCoarse,
                              bool a_homogeneousPhysBC)
{
  LevelData<NodeFArrayBox>& phi = (LevelData<NodeFArrayBox>&)a_phi;

  m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse, false);
  //apply physical boundary conditions in applyOp
  this->applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC );
}

void AMRNodeOp::AMROperator(LevelData<NodeFArrayBox>&       a_LofPhi,
                            const LevelData<NodeFArrayBox>& a_phiFine,
                            const LevelData<NodeFArrayBox>& a_phi,
                            const LevelData<NodeFArrayBox>& a_phiCoarse,
                            bool                            a_homogeneousPhysBC,
                            AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp)
{
  LevelData<NodeFArrayBox>& phi = (LevelData<NodeFArrayBox>&)a_phi;

  projectFineInterior(phi, a_phiFine);
  m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse, false);
  m_levelOps.setToZero(a_LofPhi);
  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);
}

void AMRNodeOp::AMROperatorNC(LevelData<NodeFArrayBox>& a_LofPhi,
                              const LevelData<NodeFArrayBox>& a_phiFine,
                              const LevelData<NodeFArrayBox>& a_phi,
                              bool a_homogeneousPhysBC,
                              AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp)
{
  LevelData<NodeFArrayBox>& phi = (LevelData<NodeFArrayBox>&)a_phi;

  //AMRNodeOp* finerAMRNodeOp = (AMRNodeOp*) a_finerOp;
  projectFineInterior(phi, a_phiFine);
  //no coarse-fine interpolation here
  m_levelOps.setToZero(a_LofPhi);
  applyOp(a_LofPhi, phi, a_homogeneousPhysBC);

}

void AMRNodeOp::AMRRestrict(LevelData<NodeFArrayBox>& a_resCoarse,
                            const LevelData<NodeFArrayBox>& a_residual,
                            const LevelData<NodeFArrayBox>& a_correction,
                            const LevelData<NodeFArrayBox>& a_coarseCorrection,
                            bool a_skip_res )
{
  CH_assert(!a_skip_res);
  LevelData<NodeFArrayBox> r;
  create(r, a_residual);
  m_levelOps.setToZero(r);

  AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);

  NodeCoarseAverage  averageOp(a_residual.disjointBoxLayout(),
                               a_resCoarse.disjointBoxLayout(),
                               1, m_refToCoarser, m_domain);

  m_levelOps.setToZero(a_resCoarse);
  averageOp.averageToCoarse(a_resCoarse, r);
}

/** a_correction += I[2h->h](a_coarseCorrection) */
void AMRNodeOp::AMRProlong(LevelData<NodeFArrayBox>& a_correction,
                           const LevelData<NodeFArrayBox>& a_coarseCorrection)
{
  DisjointBoxLayout c;
  coarsen(c,  a_correction.disjointBoxLayout(), m_refToCoarser);
  LevelData<NodeFArrayBox> e(c, a_correction.nComp(),a_coarseCorrection.ghostVect());
  a_coarseCorrection.copyTo(e.interval(), e, e.interval());
  prolongIncrement(a_correction, e, m_refToCoarser);
}

void AMRNodeOp::AMRUpdateResidual(LevelData<NodeFArrayBox>& a_residual,
                                  const LevelData<NodeFArrayBox>& a_correction,
                                  const LevelData<NodeFArrayBox>& a_coarseCorrection)
{
  LevelData<NodeFArrayBox> r;
  this->create(r, a_residual);
  m_levelOps.setToZero(r);
  this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  this->assign(a_residual, r);
  zeroBoundaryNodes(a_residual, m_IVSVext);
}

//==========================================================
//
//==========================================================

// MultiGrid define function

void AMRNodeOpFactory::define(const ProblemDomain& a_domain,
                              const DisjointBoxLayout& a_grid,
                              const Real& a_dx,
                              NodeBCFunc a_bc,
                              int maxDepth,
                              Real a_alpha,
                              Real a_beta)
{

  Vector<DisjointBoxLayout> grids(1, a_grid);
  Vector<int> refRatio(1, 2);
  define(a_domain, grids, refRatio, a_dx, a_bc, a_alpha, a_beta);
}

AMRNodeOp*
AMRNodeOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                          int depth,
                          bool homoOnly)
{
  // CH_assert(m_boxes.size()>depth);

  int ref=0;
  for (;ref< m_domains.size(); ref++)
    {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) break;
    }
  CH_assert(ref !=  m_domains.size()); // didn't find domain

  DisjointBoxLayout layout(m_boxes[ref]);
  ProblemDomain domain(m_domains[ref]);
  Real dx = m_dx[ref];

  for (int i=0; i<depth; i++)
    {
      if (!layout.coarsenable(2)) return NULL;
      DisjointBoxLayout dbl;
      coarsen_dbl(dbl, layout, 2);
      layout = dbl;
      dx*=2;
      domain.coarsen(2);
    }

  AMRNodeOp* newOp = new AMRNodeOp;
  newOp->define(layout, dx, domain, m_bc);
  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;
  return newOp;
}

AMRNodeOp* AMRNodeOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{

  AMRNodeOp* newOp = new AMRNodeOp;
  int ref = 0;
  for (;ref< m_domains.size(); ref++)
    {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) break;
    }

  if (ref == 0)
    {
      // coarsest AMR level
      if (m_domains.size() == 1)
        {
          //no finer level
          newOp->define(m_boxes[0],  m_dx[0],
                        a_indexSpace, m_bc);
        }
      else
        {
          //finer level exists but no coarser
          int dummyRat = 1;  //argument so compiler can find right function
          int refToFiner = m_refRatios[0]; //actual refinement ratio
          newOp->define(m_boxes[0],  m_boxes[1], m_dx[0],
                        dummyRat, refToFiner,
                        a_indexSpace, m_bc);
        }
    }
  else if (ref ==  m_domains.size()-1)
    {
      // finest AMR level
      newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1],
                    a_indexSpace, m_bc);
    }
  else if ( ref == m_domains.size())
    {
      MayDay::Error("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");

    }
  else
    {
      // intermediate AMR level, full define
      newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1], m_refRatios[ref], a_indexSpace, m_bc);
    }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;
  return newOp;
}

//  AMR Factory define function
void AMRNodeOpFactory::define(const ProblemDomain& a_coarseDomain,
                              const Vector<DisjointBoxLayout>& a_grids,
                              const Vector<int> a_refRatios,
                              const Real& a_coarsedx,
                              NodeBCFunc a_bc,
                              Real a_alpha,
                              Real a_beta)
{
  m_alpha = a_alpha;
  m_beta = a_beta;
  m_domains.resize(a_grids.size());
  m_boxes=a_grids;
  m_refRatios=a_refRatios;
  m_dx.resize(a_grids.size());
  m_bc = a_bc;
  m_domains[0] = a_coarseDomain;
  m_dx[0] = a_coarsedx;
  for (int i=1; i<a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1]/m_refRatios[i] ;
      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);
    }
}

int AMRNodeOpFactory::refToFiner(const ProblemDomain& a_domain) const

{

  int retval = -1;
  bool found = false;
  for (int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if (m_domains[ilev].domainBox() == a_domain.domainBox())
        {
          retval = m_refRatios[ilev];
          found = true;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}
#include "NamespaceFooter.H"
