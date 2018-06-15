#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ResistivityOp.H"
#include "FORT_PROTO.H"
#include "ResistivityOpF_F.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverageFace.H"
#include "QuadCFInterp.H"
#include "LevelFluxRegister.H"
#include "NamespaceHeader.H"
const int ResistivityOp::s_nComp     = 3;
const int ResistivityOp::s_nGradComp = 3*SpaceDim;
void
rogetMultiColors(Vector<IntVect>& a_colors)
{

#if CH_SPACEDIM==2
  a_colors.resize(4);
  a_colors[0] = IntVect::Zero;//(0,0)
  a_colors[1] = IntVect::Unit;//(1,1)
  a_colors[2] = IntVect::Zero + BASISV(1);//(0,1)
  a_colors[3] = IntVect::Zero + BASISV(0);//(1,0)
#elif CH_SPACEDIM==3
  a_colors.resize(8);
  a_colors[0] = IntVect::Zero;//(0,0,0)
  a_colors[1] = IntVect::Zero + BASISV(0) + BASISV(1);//(1,1,0)
  a_colors[2] = IntVect::Zero + BASISV(1) + BASISV(2);//(0,1,1)
  a_colors[3] = IntVect::Zero + BASISV(0) + BASISV(2);//(1,0,1)
  a_colors[4] = IntVect::Zero + BASISV(1);//(0,1,0)
  a_colors[5] = IntVect::Zero + BASISV(0);//(1,0,0)
  a_colors[6] = IntVect::Zero + BASISV(2);//(0,0,1)
  a_colors[7] = IntVect::Unit;//(1,1,1)
#endif
}
/***/
void
coarsenEta(LevelData<FluxBox> &       a_etaCoar,
           const LevelData<FluxBox> & a_etaFine,
           const int &                a_refToDepth)
{
  CH_assert(a_etaCoar.nComp() == 1);
  CH_assert(a_etaFine.nComp() == 1);
  CoarseAverageFace averageOp(a_etaFine.disjointBoxLayout(), 1, a_refToDepth);

  averageOp.averageToCoarseHarmonic(a_etaCoar, a_etaFine);
}
/***/
void
ResistivityOp::
createCoarser(LevelData<FArrayBox>&       a_coarse,
              const LevelData<FArrayBox>& a_fine,
              bool ghosted)
{
  // CH_assert(!ghosted);
  IntVect ghost = a_fine.ghostVect();
  DisjointBoxLayout dbl;
  coarsen(dbl, a_fine.disjointBoxLayout(), 2); //multigrid, so coarsen by 2
  a_coarse.define(dbl, s_nComp, ghost);
}
/***/
void
ResistivityOp::
computeOperatorNoBCs(LevelData<FArrayBox>& a_lhs,
                     const LevelData<FArrayBox>& a_phi)
{

  //  Real dx = m_dx;
  const DisjointBoxLayout dbl= a_phi.getBoxes();
  DataIterator dit = a_phi.dataIterator();
  //this makes the lhs = alpha*phi
  m_levelOps.setToZero(a_lhs);
  incr(a_lhs, a_phi, m_alpha);
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box gridBox = dbl.get(dit());
      FluxBox flux(gridBox, s_nComp);
      FluxBox& eta = (*m_eta)[dit];

      //operator is alpha I + (divF) = I + beta*div(eta(grad B - grad B^T) + eta I div B )
      //beta and eta are incorporated into the flux F
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Box faceBox = surroundingNodes(gridBox, idir);
          getFlux(flux[idir], a_phi[dit()], m_grad[dit()], eta[idir], faceBox, idir);
          FArrayBox& lhsFAB = a_lhs[dit()];
          const FArrayBox& fluxFAB = flux[idir];

          //beta and eta are part of the flux.
          FORT_ADDDIVFLUXDIRROP(CHF_FRA(lhsFAB),
                                CHF_CONST_FRA(fluxFAB),
                                CHF_BOX(gridBox),
                                CHF_REAL(m_dx),
                                CHF_CONST_INT(idir));
        }
    } // end loop over boxes
}
/***/
void
ResistivityOp::
restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                 LevelData<FArrayBox>&       a_phiFine,
                 const LevelData<FArrayBox>& a_rhsFine )
{
  // temp storage
  LevelData<FArrayBox> resFine;
  create(resFine, a_rhsFine);
  homogeneousCFInterp(a_phiFine);
  //bcs and exchange done within applyOp
  residual(resFine, a_phiFine, a_rhsFine, true);
  const DisjointBoxLayout dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      Box region = dblFine.get(dit());
      a_resCoarse[dit()].setVal(0.0);
      FORT_RESTRICTRESROP(CHF_FRA(a_resCoarse[dit()]),
                          CHF_CONST_FRA(resFine[dit()]),
                          CHF_BOX(region));
    }
}
/***/
void
ResistivityOp::
prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                 const LevelData<FArrayBox>& a_correctCoarse)
{
  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  int mgref = 2; //this is a multigrid func
  for (DataIterator dit = a_phiThisLevel.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_phiThisLevel[dit];
      const FArrayBox& coarse = a_correctCoarse[dit];
      Box region = dbl.get(dit());
      Box cBox = coarsen(region, mgref);

      FORT_PROLONGROP(CHF_FRA(phi),
                      CHF_CONST_FRA(coarse),
                      CHF_BOX(region),
                      CHF_CONST_INT(mgref));

    }
}
/***/
void
ResistivityOp::
AMRResidualNF(LevelData<FArrayBox>&       a_residual,
              const LevelData<FArrayBox>& a_phi,
              const LevelData<FArrayBox>& a_phiCoarse,
              const LevelData<FArrayBox>& a_rhs,
              bool a_homogeneousPhysBC)
{
  this->cfinterp(a_phi, a_phiCoarse);
  this->residual(a_residual, a_phi, a_rhs, a_homogeneousPhysBC ); //apply boundary conditions
}

/***/
void
ResistivityOp::
AMRResidual(LevelData<FArrayBox>&       a_residual,
            const LevelData<FArrayBox>& a_phiFine,
            const LevelData<FArrayBox>& a_phi,
            const LevelData<FArrayBox>& a_phiCoarse,
            const LevelData<FArrayBox>& a_rhs,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)

{
  //fillgrad is called in applyop
  this->cfinterp(a_phi, a_phiCoarse);

  applyOp(a_residual, a_phi, a_homogeneousPhysBC);

  if (a_finerOp != NULL)
    {
      reflux(a_phiFine, a_phi,  a_residual, a_finerOp);
    }
  incr(a_residual, a_rhs, -1.0);
  // residual is rhs - L(phi)
  scale(a_residual, -1.0);
}

/***/
void
ResistivityOp::
AMRResidualNC(LevelData<FArrayBox>&       a_residual,
              const LevelData<FArrayBox>& a_phiFine,
              const LevelData<FArrayBox>& a_phi,
              const LevelData<FArrayBox>& a_rhs,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  //no coarse-fine interpolation here--fillgrad is called in applyop
  m_levelOps.setToZero(m_grad);
  applyOp(a_residual, a_phi, a_homogeneousPhysBC);
  reflux(a_phiFine, a_phi,  a_residual, a_finerOp);
  axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}
/***/
void
ResistivityOp::
AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
            const LevelData<FArrayBox>& a_residual,
            const LevelData<FArrayBox>& a_correction,
            const LevelData<FArrayBox>& a_coarseCorrection, 
            bool a_skip_res )
{  
  CH_assert(!a_skip_res);
  LevelData<FArrayBox> r;
  create(r, a_residual);
  AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  DisjointBoxLayout dblCoar = a_resCoarse.disjointBoxLayout();
  DataIterator dit = a_residual.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& coarse = a_resCoarse[dit];
      const FArrayBox& fine = r[dit];
      const Box& b = dblCoar.get(dit());
      Box refbox(IntVect::Zero,
                 (m_refToCoar-1)*IntVect::Unit);
      FORT_AVERAGE( CHF_FRA(coarse),
                    CHF_CONST_FRA(fine),
                    CHF_BOX(b),
                    CHF_CONST_INT(m_refToCoar),
                    CHF_BOX(refbox)
                    );
    }
}
/***/
void
ResistivityOp::
AMRProlong(LevelData<FArrayBox>& a_correction,
           const LevelData<FArrayBox>& a_coarseCorrection)
{
  DisjointBoxLayout c;
  coarsen(c,  a_correction.disjointBoxLayout(), m_refToCoar);
  LevelData<FArrayBox> eCoar(c, s_nComp, a_coarseCorrection.ghostVect());
  a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  for (DataIterator dit = a_correction.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_correction[dit];
      const FArrayBox& coarse = eCoar[dit];
      Box region = dbl.get(dit());

      FORT_PROLONGROP(CHF_FRA(phi),
                      CHF_CONST_FRA(coarse),
                      CHF_BOX(region),
                      CHF_CONST_INT(m_refToCoar));

    }
}

/***/
void
ResistivityOp::
AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
                  const LevelData<FArrayBox>& a_correction,
                  const LevelData<FArrayBox>& a_coarseCorrection)
{
  LevelData<FArrayBox> r;
  this->create(r, a_residual);
  this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  this->assign(a_residual, r);
}
/***/
void
ResistivityOp::
createCoarsened(LevelData<FArrayBox>&       a_lhs,
                const LevelData<FArrayBox>& a_rhs,
                const int &                 a_refRat)
{
  IntVect ghostVect = a_rhs.ghostVect();

  DisjointBoxLayout dbl = a_rhs.disjointBoxLayout();
  CH_assert(dbl.coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, dbl, a_refRat);

  a_lhs.define(dblCoarsenedFine, s_nComp, a_rhs.ghostVect());
}

/***/
void ResistivityOp::preCond(LevelData<FArrayBox>&       a_phi,
                            const LevelData<FArrayBox>& a_rhs)
{
  relax(a_phi, a_rhs, 10);
}

/***/
void
ResistivityOp::
applyOp(LevelData<FArrayBox>&       a_lhs,
        const LevelData<FArrayBox>& a_phi,
        bool                        a_homogeneous )
{
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
    }

  //contains an exchange
  this->fillGrad(a_phi);

  computeOperatorNoBCs(a_lhs, phi);
}

void ResistivityOp::applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                      const LevelData<FArrayBox>& a_phi)
{
  this->fillGrad(a_phi);
  computeOperatorNoBCs(a_lhs, a_phi);
}

/***/
void
ResistivityOp::
divergenceCC(LevelData<FArrayBox>&       a_div,
             const LevelData<FArrayBox>& a_phi,
             const LevelData<FArrayBox>* a_phiC)
{
  //fill ghost cells of phi
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  //if necessary, apply coarse-fine boundary conditions
  if (a_phiC != NULL)
    {
      cfinterp(a_phi, *a_phiC);
    }
  //then apply boundary conditions at the domain
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, m_dx, false);
    }
  //exchange with neighboring boxes
  phi.exchange(phi.interval(), m_exchangeCopier);
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box  grid = dbl.get(dit());
      a_div[dit()].setVal(0.);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FORT_CELLDIVINCRROP(CHF_FRA1     (a_div[dit()], 0),
                              CHF_CONST_FRA(phi  [dit()]   ),
                              CHF_CONST_REAL(m_dx),
                              CHF_CONST_INT(idir),
                              CHF_BOX(grid));
        }
    }
}
void ResistivityOp::reflux(const LevelData<FArrayBox>& a_phiFine,
                           const LevelData<FArrayBox>& a_phi,
                           LevelData<FArrayBox>& residual,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  ProblemDomain fineDomain = refine(m_domain, m_refToFine);
  LevelFluxRegister levfluxreg(a_phiFine.disjointBoxLayout(),
                               a_phi.disjointBoxLayout(),
                               fineDomain,
                               m_refToFine,
                               s_nComp);

  levfluxreg.setToZero();
  Interval interv(0,s_nComp-1);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab = a_phi[dit];
      const FluxBox& coareta = (*m_eta)[dit];
      const Box& gridBox = a_phi.getBoxes()[dit];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox coarflux;
          Box faceBox = surroundingNodes(gridBox, idir);
          //fillgrad was called before this when applyOp was called
          getFlux(coarflux, coarfab, m_grad[dit()], coareta[idir], faceBox, idir);

          Real scale = 1.0;
          levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                     interv,interv,idir);
        }
    }
  LevelData<FArrayBox>& p = ( LevelData<FArrayBox>&)a_phiFine;

  //interpolate finer phi and compute gradients.
  ResistivityOp* finerAMRPOp = (ResistivityOp*) a_finerOp;
  finerAMRPOp->cfinterp(p, a_phi);
  finerAMRPOp->fillGrad(p);

  IntVect phiGhost = p.ghostVect();

  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = a_phiFine[ditf];
      const FArrayBox& fineTanGrad = finerAMRPOp->m_grad[ditf];
      const FluxBox& fineeta = (*(finerAMRPOp->m_eta))[ditf];
      const Box& gridbox = dblFine.get(ditf());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          int normalGhost = phiGhost[idir];
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide hiorlo = sit();
              Box fabbox;
              Box facebox;

              // assumption here that the stencil required
              // to compute the flux in the normal direction
              // is 2* the number of ghost cells for phi
              // (which is a reasonable assumption, and probably
              // better than just assuming you need one cell on
              // either side of the interface
              // (dfm 8-4-06)
              if (sit() == Side::Lo)
                {
                  fabbox = adjCellLo(gridbox,idir, 2*normalGhost);
                  fabbox.shift(idir, 1);
                  facebox = bdryLo(gridbox, idir,1);
                }
              else
                {
                  fabbox = adjCellHi(gridbox,idir, 2*normalGhost);
                  fabbox.shift(idir, -1);
                  facebox = bdryHi(gridbox, idir, 1);
                }

              // just in case we need ghost cells in the transverse direction
              // (dfm 8-4-06)
              for (int otherDir=0; otherDir<SpaceDim; ++otherDir)
                {
                  if (otherDir != idir)
                    {
                      fabbox.grow(otherDir, phiGhost[otherDir]);
                    }
                }
              CH_assert(!fabbox.isEmpty());

              FArrayBox phifab(fabbox, s_nComp);
              phifab.copy(phifFab);

              FArrayBox fineflux;
              getFlux(fineflux, phifab, fineTanGrad, fineeta[idir], facebox, idir,
                      m_refToFine);

              Real scale = 1.0;
              levfluxreg.incrementFine(fineflux, scale, ditf(),
                                       interv, interv, idir, hiorlo);
            }
        }
    }

  Real scale =  1.0/m_dx;
  levfluxreg.reflux(residual, scale);
}
/****/

/****/
void
ResistivityOp::
loHiCenterFace(Box&                 a_loBox,
               int&                 a_hasLo,
               Box&                 a_hiBox,
               int&                 a_hasHi,
               Box&                 a_centerBox,
               const ProblemDomain& a_eblg,
               const Box&           a_inBox,
               const int&           a_dir)
{
  if (a_eblg.isPeriodic(a_dir))
    {
      Box domainFaceBox = a_eblg.domainBox();
      domainFaceBox.surroundingNodes(a_dir);
      a_hasLo = 0;
      a_hasHi = 0;
      Box inBox = a_inBox;
      inBox &= domainFaceBox;
      a_centerBox = inBox;
    }
  else
    {
      Box domainFaceBox = a_eblg.domainBox();
      domainFaceBox.surroundingNodes(a_dir);

      // Make a copy of the input box which can be modified
      Box inBox = a_inBox;

      inBox &= domainFaceBox;

      a_centerBox = inBox;
      a_centerBox.grow(a_dir, 1);
      a_centerBox &= domainFaceBox;
      a_centerBox.grow(a_dir,-1);

      // See if this chops off the high side of the input box
      Box tmp = inBox;
      tmp.shift(a_dir,1);
      tmp &= domainFaceBox;
      tmp.shift(a_dir,-1);

      // If so, set up the high, one-sided difference box, a_hiBox
      if (tmp != inBox)
        {
          a_hasHi = 1;
          a_hiBox = bdryHi(inBox,a_dir,1);
        }
      else
        {
          a_hasHi = 0;
        }

      // See if this chops off the low side of the input box
      tmp = inBox;
      tmp.shift(a_dir,-1);
      tmp &= domainFaceBox;
      tmp.shift(a_dir,1);

      // If so, set up the low, one-sided difference box, a_loBox
      if (tmp != inBox)
        {
          a_hasLo = 1;
          a_loBox = bdryLo(inBox,a_dir,1);
        }
      else
        {
          a_hasLo = 0;
        }
    }

}
/***/
void
ResistivityOp::
getFaceDivAndGrad(FArrayBox&             a_faceDiv,
                  FArrayBox&             a_faceGrad,
                  const FArrayBox&       a_data,
                  const FArrayBox&       a_gradData,
                  const ProblemDomain&   a_domain,
                  const Box&             a_faceBox,
                  const int&             a_faceDir,
                  const Real             a_dx)
{
  //set divergence to zero everywhere
  a_faceDiv.setVal(0.);

  Box interiorFaceBox ;
  Box loBoxFace, hiBoxFace;
  int hasLo, hasHi;

  loHiCenterFace(loBoxFace, hasLo,
                 hiBoxFace, hasHi,
                 interiorFaceBox,
                 a_domain,
                 a_faceBox,
                 a_faceDir);

  for (int divDir = 0; divDir < SpaceDim; divDir++)
    {
      int gradComp = TensorCFInterp::gradIndex(divDir,divDir);
      FORT_FACEDIVINCRROP(CHF_FRA1(a_faceDiv, 0),
                          CHF_FRA(a_data),
                          CHF_FRA(a_gradData),
                          CHF_BOX(a_faceBox),
                          CHF_BOX(interiorFaceBox),
                          CHF_BOX(loBoxFace),
                          CHF_INT(hasLo),
                          CHF_BOX(hiBoxFace),
                          CHF_INT(hasHi),
                          CHF_REAL(a_dx),
                          CHF_INT(a_faceDir),
                          CHF_INT(divDir),
                          CHF_INT(gradComp));

      //now average cell-centered gradients to the face centers
      //use diffs in data if divDir == faceDir
    }

  for (int divDir = 0; divDir < SpaceDim; divDir++)
    {
      //always 3 components  in mhd land
      for (int velDir = 0; velDir < s_nComp; velDir++)
        {
          int gradcomp = TensorCFInterp::gradIndex(velDir,divDir);
          FORT_GETFACEGRADROP(CHF_FRA1(a_faceGrad, gradcomp),
                              CHF_FRA1(a_gradData, gradcomp),
                              CHF_FRA1(a_data, velDir),
                              CHF_BOX(a_faceBox),
                              CHF_BOX(interiorFaceBox),
                              CHF_BOX(loBoxFace),
                              CHF_INT(hasLo),
                              CHF_BOX(hiBoxFace),
                              CHF_INT(hasHi),
                              CHF_REAL(a_dx),
                              CHF_INT(a_faceDir),
                              CHF_INT(divDir));

        }

    }
}
/***/
void
ResistivityOp::
getFlux(FluxBox&         a_flux,
        const LevelData<FArrayBox>& a_data,
        const Box&       a_grid,
        const DataIndex& a_dit,
        Real             a_scale)
{
  const FArrayBox& data = a_data[a_dit];
  a_flux.define(a_grid, s_nComp);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const FArrayBox& etaFace  = (*m_eta)[a_dit][idir];
      const FArrayBox& gradData =   m_grad[a_dit];
      Box faceBox = a_flux[idir].box();
      Box domFaceBox = surroundingNodes(m_domain.domainBox(), idir);
      faceBox &= domFaceBox;
      getFlux(a_flux[idir], data, gradData, etaFace, faceBox, idir, 1);
      a_flux[idir] *= a_scale;
    }
}

void ResistivityOp::getFlux(FArrayBox&       a_flux,
                            const FArrayBox& a_data,
                            const FArrayBox& a_gradData,
                            const FArrayBox& a_etaFace,
                            const Box&       a_faceBox,
                            int              a_faceDir,
                            int a_ref)
{
  ProblemDomain domain(m_domain);
  domain.refine(a_ref);
  Real dx(m_dx);

  dx /= a_ref;
  CH_assert(a_data.nComp() == s_nComp);
  a_flux.resize(a_faceBox, s_nComp);
  a_flux.setVal(0.);
  FArrayBox  faceDiv(a_faceBox, 1);
  FArrayBox  faceGrad(a_faceBox, s_nGradComp);
  getFaceDivAndGrad(faceDiv, faceGrad, a_data, a_gradData, domain, a_faceBox, a_faceDir, dx);

  //copy the divergence into the diagonal component of the flux
  //diagonal component of the flux is the a_faceDir comp
  int dstcomp = a_faceDir;
  int srccomp = 0;
  int numcomp = 1;

  a_flux.copy(faceDiv, srccomp, dstcomp, numcomp);
  //the flux now contains div v * I
  //now add in grad v - (grad v)^T
  //There are always 3 components.  There are SpaceDim directions.
  for (int component = 0; component < s_nComp; component++)
    {
      if (component != a_faceDir)
        {
          int fluxComp = component;
          {
            int gradComp = TensorCFInterp::gradIndex(component, a_faceDir);
            Real gradsign =  1.0;
            FORT_ADDGRADTOFLUXROP(CHF_FRA(a_flux),
                                    CHF_INT(fluxComp),
                                    CHF_FRA(faceGrad),
                                    CHF_INT( gradComp),
                                    CHF_REAL(gradsign),
                                    CHF_BOX(a_faceBox));
          }
          //don't 3d derivs in 2d
          if (component < SpaceDim)
            {
              int tranComp = TensorCFInterp::gradIndex(a_faceDir, component);
              Real transign = -1.0;
              FORT_ADDGRADTOFLUXROP(CHF_FRA(a_flux),
                                    CHF_INT(fluxComp),
                                    CHF_FRA(faceGrad),
                                    CHF_INT( tranComp),
                                    CHF_REAL(transign),
                                    CHF_BOX(a_faceBox));
            }
        }
    }
  //now flux gets multiplied by eta*beta
  //so that div F = beta*div(eta F)
  for (int icomp = 0; icomp < s_nComp; icomp++)
    {
      int srccomp = 0;
      int dstcomp = icomp;
      a_flux.mult(a_etaFace, srccomp, dstcomp, 1);
    }
  a_flux *= m_beta;
}
/***/
void
ResistivityOp::
residual(LevelData<FArrayBox>&       a_lhs,
         const LevelData<FArrayBox>& a_phi,
         const LevelData<FArrayBox>& a_rhs,
         bool a_homogeneous)
{
  applyOp(a_lhs, a_phi, a_homogeneous);
  incr(a_lhs, a_rhs, -1);
  scale(a_lhs, -1.0);
}
/***/
Real
ResistivityOp::
AMRNorm(const LevelData<FArrayBox>& a_coarResid,
        const LevelData<FArrayBox>& a_fineResid,
        const int& a_refRat,
        const int& a_ord)

{
  const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();
  const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();

  //create temp and zero out under finer grids
  LevelData<FArrayBox> coarTemp;
  m_levelOps.create(coarTemp, a_coarResid);
  m_levelOps.assign(coarTemp, a_coarResid);
  int ncomp = coarTemp.nComp();
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& coarTempFAB = coarTemp[dit()];
      LayoutIterator litFine = fineGrids.layoutIterator();
      for (litFine.reset(); litFine.ok(); ++litFine)
        {
          Box overlayBox = coarTempFAB.box();
          Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat);

          overlayBox &= coarsenedGrid;

          if (!overlayBox.isEmpty())
            {
              coarTempFAB.setVal(0.0,overlayBox,0, ncomp);
            }
        }
    }
  //return norm of temp
  return norm(coarTemp, a_ord);
}
/***/
void
ResistivityOp::
create(LevelData<FArrayBox>&       a_lhs,
       const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.create(a_lhs, a_rhs);
}
/***/
void
ResistivityOp::
assign(LevelData<FArrayBox>&       a_lhs,
       const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.assign(a_lhs, a_rhs);
}
/***/
Real
ResistivityOp::
dotProduct(const LevelData<FArrayBox>& a_1,
           const LevelData<FArrayBox>& a_2)
{
  return m_levelOps.dotProduct(a_1, a_2);
}
/***/
void
ResistivityOp::
incr( LevelData<FArrayBox>&       a_lhs,
      const LevelData<FArrayBox>& a_x,
      Real a_scale)
{
  m_levelOps.incr(a_lhs, a_x, a_scale);
}
/***/
void
ResistivityOp::
axby( LevelData<FArrayBox>&       a_lhs,
      const LevelData<FArrayBox>& a_x,
      const LevelData<FArrayBox>& a_y,
      Real a_a, Real a_b)
{
  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}
/***/
void
ResistivityOp::
scale(LevelData<FArrayBox>& a_lhs,
      const Real& a_scale)
{
  m_levelOps.scale(a_lhs, a_scale);
}
/***/
void
ResistivityOp::
setToZero(LevelData<FArrayBox>& a_lhs)
{
  m_levelOps.setToZero(a_lhs);
}

/***/
void
ResistivityOp::
relax(LevelData<FArrayBox>&       a_phi,
      const LevelData<FArrayBox>& a_rhs,
      int a_iterations)
{
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  LevelData<FArrayBox> lphi;
  m_levelOps.create(lphi, a_rhs);
  // do first red, then black passes
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for (int icolor = 0; icolor < m_colors.size(); icolor++)
        {
          const IntVect& color= m_colors[icolor];
          homogeneousCFInterp(a_phi);

          //after this lphi = L(phi)
          //this call contains bcs and exchange
          applyOp(lphi, a_phi, true);

          for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
            {
              Box dblBox  = dbl.get(dit());

              IntVect loIV = dblBox.smallEnd();
              IntVect hiIV = dblBox.bigEnd();

              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (loIV[idir] % 2 != color[idir])
                    {
                      loIV[idir]++;
                    }
                }

              if (loIV <= hiIV)
                {
                  Box coloredBox(loIV, hiIV);
                  FORT_GSRBROP(CHF_FRA(a_phi[dit]),
                               CHF_CONST_FRA( lphi[dit]),
                               CHF_CONST_FRA(a_rhs[dit]),
                               CHF_CONST_FRA(m_lambda[dit]),
                               CHF_BOX(coloredBox));
                }
            }
        } // end loop through red-black
    }
}
/***/
Real
ResistivityOp::
norm(const LevelData<FArrayBox>& a_x, int a_ord)
{
  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}
/***/
ResistivityOp::
ResistivityOp(const DisjointBoxLayout&                    a_grids,
              const DisjointBoxLayout&                    a_gridsFine,
              const DisjointBoxLayout&                    a_gridsCoar,
              const RefCountedPtr<LevelData<FluxBox> >&   a_eta,
              Real                                        a_alpha,
              Real                                        a_beta,
              int                                         a_refToFine,
              int                                         a_refToCoar,
              const ProblemDomain&                        a_domain,
              const Real&                                 a_dx,
              const Real&                                 a_dxCrse,
              BCFunc                                      a_bc)
{
  rogetMultiColors(m_colors);
  m_grids = a_grids;
  m_eta           =  a_eta;
  m_alpha         =  a_alpha;
  m_beta          =  a_beta;
  m_refToCoar     =  a_refToCoar;
  m_refToFine     =  a_refToFine;
  m_domain        =  a_domain;
  m_dx            =  a_dx;
  m_dxCrse        =  a_dxCrse;
  m_bc            =  a_bc;

  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  if (a_gridsCoar.isClosed())
    {
      m_interpWithCoarser.define(a_grids, &a_gridsCoar, a_dx,
                                 m_refToCoar, s_nComp, m_domain);
    }

  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      LayoutData<CFIVS>& lo =  m_loCFIVS[i];
      LayoutData<CFIVS>& hi =  m_hiCFIVS[i];
      lo.define(a_grids);
      hi.define(a_grids);
      for (DataIterator dit(a_grids); dit.ok(); ++dit)
        {
          lo[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Lo);
          hi[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Hi);
        }
    }

  //define lambda, the relaxation coef
  m_lambda.define(a_grids, s_nComp,     IntVect::Zero);
  m_grad.define(  a_grids, s_nGradComp, IntVect::Unit);
  DataIterator lit = a_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(a_grids);
      m_hiCFIVS[idir].define(a_grids);
      m_loTanStencilSets[idir].define(a_grids);
      m_hiTanStencilSets[idir].define(a_grids);

      for (lit.begin(); lit.ok(); ++lit)
        {
          m_loCFIVS[idir][lit()].define(m_domain, a_grids.get(lit()),
                                        a_grids, idir,Side::Lo);
          m_hiCFIVS[idir][lit()].define(m_domain, a_grids.get(lit()),
                                        a_grids, idir,Side::Hi);

          const IntVectSet& fineIVSlo = m_loCFIVS[idir][lit()].getFineIVS();
          m_loTanStencilSets[idir][lit()].define(fineIVSlo,
                                                 m_domain, idir);

          const IntVectSet& fineIVShi = m_hiCFIVS[idir][lit()].getFineIVS();
          m_hiTanStencilSets[idir][lit()].define(fineIVShi,
                                                 m_domain, idir);

        }
    }
  setLambda();
  m_levelOps.setToZero(m_grad);
}
/***/
void
ResistivityOp::
setLambda()
{
  Real safety = 1.0;
  for (DataIterator dit(m_grids); dit.ok(); ++dit)
    {
      const Box& grid = m_grids.get(dit());
      //first make lambda = alpha - beta*sum_d(etahi + eta lo)/h^2
      m_lambda[dit()].setVal(m_alpha);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FORT_DECRINVLAMBDAROP(CHF_FRA(m_lambda[dit()]),
                                CHF_FRA((*m_eta)[dit()][idir]),
                                CHF_BOX(grid),
                                CHF_REAL(m_beta),
                                CHF_REAL(m_dx),
                                CHF_INT(idir));
        }

      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)
      FORT_INVERTLAMBDAROP(CHF_FRA(m_lambda[dit()]),
                           CHF_REAL(safety),
                           CHF_BOX(grid));
    }
}
/***/
ResistivityOpFactory::
ResistivityOpFactory(const Vector<DisjointBoxLayout>&                     a_grids,
                     const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_eta,
                     Real                                                 a_alpha,
                     Real                                                 a_beta,
                     const Vector<int>&                                   a_refRatios,
                     const ProblemDomain&                                 a_coarseDomain,
                     const Real&                                          a_coarseDx,
                     BCFunc                                               a_bc)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;
  m_eta   = a_eta;
  m_domains.resize(a_grids.size());
  m_boxes=a_grids;
  m_refRatios=a_refRatios;
  m_dx.resize(a_grids.size());
  m_bc = a_bc;
  m_domains[0] = a_coarseDomain;
  m_dx[0] = a_coarseDx;
  for (int i=1; i<a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1]/m_refRatios[i] ;
      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);
    }
}
/***/
MGLevelOp<LevelData<FArrayBox> >*
ResistivityOpFactory::
MGnewOp(const ProblemDomain& a_indexSpace,
        int                  a_depth,
        bool                 a_homoOnly)
{
  int ref = 0;
  bool found = false;
  for (int ivec = 0; ivec < m_domains.size(); ivec++)
    {
      if (a_indexSpace.domainBox() == m_domains[ivec].domainBox())
        {
          found = true;
          ref = ivec;
          break;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }

  DisjointBoxLayout layout(m_boxes[ref]);
  ProblemDomain domain(m_domains[ref]);
  Real dx = m_dx[ref];

  int refToDepth = 1;
  for (int i=0; i< a_depth; i++)
    {
      if (!layout.coarsenable(4)) return NULL;
      Box refCoarDom = domain.domainBox();
      coarsen(refCoarDom, 4);
      refine(refCoarDom, 4);
      if (refCoarDom != domain.domainBox())
        {
          return NULL;
        }
      DisjointBoxLayout dbl;
      coarsen_dbl(dbl, layout, 2);
      layout = dbl;
      dx*=2;
      refToDepth *= 2;
      domain.coarsen(2);
    }

  RefCountedPtr<LevelData<FluxBox> > eta( new LevelData<FluxBox>(layout, 1, IntVect::Zero) );
  coarsenEta(*eta, *m_eta[ref], refToDepth);
  //no coarser or finer grids
  //no reftocoar, reftofine
  Real dxCrse = 2*m_dx[ref];
  if (ref > 0)
    {
      dxCrse = m_dx[ref-1];
    }
  ResistivityOp* newOp = new ResistivityOp(layout, DisjointBoxLayout(), DisjointBoxLayout(),
                                           eta, m_alpha, m_beta, -1, -1,
                                           domain,  dx, dxCrse, m_bc);
  MGLevelOp<LevelData<FArrayBox> >* retval = (MGLevelOp<LevelData<FArrayBox> >*) newOp;
  return retval;
}
void
ResistivityOp::
AMROperatorNC(LevelData<FArrayBox>&       a_LofPhi,
              const LevelData<FArrayBox>& a_phiFine,
              const LevelData<FArrayBox>& a_phi,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  //no coarse-fine interpolation here
  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);
  reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
}
void
ResistivityOp::
AMROperatorNF(LevelData<FArrayBox>& a_LofPhi,
              const LevelData<FArrayBox>& a_phi,
              const LevelData<FArrayBox>& a_phiCoarse,
              bool a_homogeneousPhysBC)
{
  //fillgrad is called in applyop
  cfinterp(a_phi, a_phiCoarse);

  //apply boundary conditions in applyOp
  this->applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC );
}

void
ResistivityOp::
AMROperator(      LevelData<FArrayBox>& a_LofPhi,
                  const LevelData<FArrayBox>& a_phiFine,
                  const LevelData<FArrayBox>& a_phi,
                  const LevelData<FArrayBox>& a_phiCoarse,
                  bool a_homogeneousPhysBC,
                  AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  //fillgrad is called in applyop
  cfinterp(a_phi, a_phiCoarse);

  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);
  if (a_finerOp != NULL)
    {
      reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
    }
}
/***/
AMRLevelOp<LevelData<FArrayBox> >*
ResistivityOpFactory::
AMRnewOp(const ProblemDomain& a_indexSpace)
{
  int ref = 0;
  bool found = false;
  for (int ivec = 0; ivec < m_domains.size(); ivec++)
    {
      if (a_indexSpace.domainBox() == m_domains[ivec].domainBox())
        {
          found = true;
          ref = ivec;
          break;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }

  int refToCoar = 2;
  DisjointBoxLayout dblFine, dblCoar;
  Real dxCrse = 2*m_dx[ref];
  if (ref > 0) //not at coarsest level
    {
      dblCoar   = m_boxes[ref-1];
      refToCoar = m_refRatios[ref-1];
      dxCrse    = m_dx[ref-1];
    }

  if (ref < (m_domains.size()-1)) //not at finest level
    {
      dblFine = m_boxes[ref+1];
    }
  ResistivityOp* newOp = new ResistivityOp(m_boxes[ref], dblFine, dblCoar,
                                           m_eta[ref],   m_alpha, m_beta,
                                           m_refRatios[ref], refToCoar,
                                           m_domains[ref],  m_dx[ref],
                                           dxCrse, m_bc);

  AMRLevelOp<LevelData< FArrayBox> >* retval = (AMRLevelOp<LevelData<FArrayBox> >*) newOp;
  return retval;
}
/***/
int
ResistivityOpFactory::
refToFiner(const ProblemDomain& a_domain) const
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
/**/
void
ResistivityOp::
cellGrad(FArrayBox&             a_gradPhi,
         const  FArrayBox&      a_phi,
         const Box&             a_grid)
{
  //always three components
  CH_assert(a_gradPhi.nComp() == s_nGradComp);
  CH_assert(a_phi.nComp()     == s_nComp);

  for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      for (int phiDir = 0; phiDir < s_nComp;  phiDir++)
        {
          int gradcomp = TensorCFInterp::gradIndex(phiDir,derivDir);
          FORT_CELLGRADROP(CHF_FRA1(a_gradPhi, gradcomp),
                           CHF_CONST_FRA1(a_phi, phiDir),
                           CHF_BOX(a_grid),
                           CHF_CONST_REAL(m_dx),
                           CHF_CONST_INT(derivDir));

        }
    }
}
/**/
void
ResistivityOp::
cfinterp(const LevelData<FArrayBox>&       a_phi,
         const LevelData<FArrayBox>&       a_phiCoarse)
{
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  m_levelOps.setToZero(m_grad);
  if (a_phiCoarse.isDefined())
    {
      m_interpWithCoarser.coarseFineInterp(phi, m_grad, a_phiCoarse);
    }
}
/**/
void
ResistivityOp::
fillGrad(const LevelData<FArrayBox>&       a_phi)
{
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  phi.exchange(phi.interval(), m_exchangeCopier);
  const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
  //compute gradient of phi for parts NOT in ghost cells
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      cellGrad(m_grad[dit()],
               a_phi [dit()],
               grids.get(dit()));
    }
  m_grad.exchange();
}
/**/
void
ResistivityOp::
homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{

  QuadCFInterp::homogeneousCFInterp(a_phif, m_grad,
                                    m_loCFIVS, m_hiCFIVS,
                                    m_dx, m_dxCrse, s_nComp,
                                    m_loTanStencilSets,m_hiTanStencilSets);
}

/***/
#include "NamespaceFooter.H"
