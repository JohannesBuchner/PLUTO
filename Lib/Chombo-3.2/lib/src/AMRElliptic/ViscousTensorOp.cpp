#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ViscousTensorOp.H"
#include "FORT_PROTO.H"
#include "ViscousTensorOpF_F.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverageFace.H"
#include "CoarseAverage.H"
#include "QuadCFInterp.H"
//#include "InterpF_F.H"
#include "NamespaceHeader.H"

// initialize static members here.
// 0 = arithmetic, 1 = harmonic
int ViscousTensorOpFactory::s_coefficientAverageType = 1;
//int ViscousTensorOp::s_prolongType = piecewiseConstant;
int ViscousTensorOp::s_prolongType = linearInterp;

void
vtogetMultiColors(Vector<IntVect>& a_colors)
{

#if CH_SPACEDIM==1
  a_colors.resize(2);
  a_colors[0] = IntVect::Zero;//(0,0)
  a_colors[1] = IntVect::Unit;//(1,1)
#elif CH_SPACEDIM==2
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
coarsenStuff(LevelData<FluxBox> &          a_etaCoar,
             LevelData<FluxBox> &          a_lambdaCoar,
             LevelData<FArrayBox> &        a_acoefCoar,
             const LevelData<FluxBox> &    a_etaFine,
             const LevelData<FluxBox> &    a_lambdaFine,
             const LevelData<FArrayBox>&   a_acoefFine,
             const int &                   a_refToDepth,
             const int &                   a_coefficientAverageType)
{
  CH_assert(a_etaCoar.nComp() == 1);
  CH_assert(a_etaFine.nComp() == 1);
  CoarseAverageFace averageOpFace(a_etaFine.disjointBoxLayout(), 1, a_refToDepth);
  CoarseAverage         averageOp(a_etaFine.disjointBoxLayout(), 1, a_refToDepth);

  if (a_coefficientAverageType == CoarseAverageFace::harmonic)
    {
      averageOpFace.averageToCoarseHarmonic(a_etaCoar,    a_etaFine);
      averageOpFace.averageToCoarseHarmonic(a_lambdaCoar, a_lambdaFine);
    }
  else
    {
      averageOpFace.averageToCoarse(a_etaCoar,    a_etaFine);
      averageOpFace.averageToCoarse(a_lambdaCoar, a_lambdaFine);
    }

  averageOp.averageToCoarse(a_acoefCoar, a_acoefFine);
}
/***/
void
ViscousTensorOp::
createCoarser(LevelData<FArrayBox>&       a_coarse,
              const LevelData<FArrayBox>& a_fine,
              bool ghosted)
{
  // CH_assert(!ghosted);
  IntVect ghost = a_fine.ghostVect();
  DisjointBoxLayout dbl;
  CH_assert(dbl.coarsenable(2));
  coarsen(dbl, a_fine.disjointBoxLayout(), 2); //multigrid, so coarsen by 2
  a_coarse.define(dbl, a_fine.nComp(), ghost);
}
/***/
void
ViscousTensorOp::
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
      for (int ivar = 0; ivar < SpaceDim; ivar++)
        {
          int src = 0; int dst = ivar; int ncomp = 1;
          a_lhs[dit()].mult((*m_acoef)[dit()],src,dst,ncomp);
        }
      Box gridBox = dbl.get(dit());
      FluxBox flux(gridBox, m_ncomp);
      FluxBox& eta    = (*m_eta   )[dit];
      FluxBox& lambda = (*m_lambda)[dit];

      //operator is alpha I + (divF) = alpha I + beta div(eta(grad B + grad B^T) + lambda I div B )
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Box faceBox = surroundingNodes(gridBox, idir);
          getFlux(flux[idir], a_phi[dit()], m_grad[dit()], eta[idir], lambda[idir],faceBox, idir);
          FArrayBox& lhsFAB = a_lhs[dit()];
          const FArrayBox& fluxFAB = flux[idir];

          //beta and eta are part of the flux.
          FORT_ADDDIVFLUXDIRVTOP(CHF_FRA(lhsFAB),
                                 CHF_CONST_FRA(fluxFAB),
                                 CHF_BOX(gridBox),
                                 CHF_REAL(m_dx),
                                 CHF_CONST_INT(m_ncomp),
                                 CHF_CONST_INT(idir));
        }
    } // end loop over boxes
}
/***/
void
ViscousTensorOp::
restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                 LevelData<FArrayBox>&       a_phiFine,
                 const LevelData<FArrayBox>& a_rhsFine)
{
  // temp storage
  LevelData<FArrayBox> resFine;
  create(resFine, a_rhsFine);
  homogeneousCFInterp(a_phiFine);
  //bcs and exchange done within applyOp
  residual(resFine, a_phiFine, a_rhsFine, true);
  int ncomp = SpaceDim;
  const DisjointBoxLayout dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      Box region = dblFine.get(dit());
      a_resCoarse[dit()].setVal(0.0);
      FORT_RESTRICTRESVTOP(CHF_FRA(a_resCoarse[dit()]),
                           CHF_CONST_FRA(resFine[dit()]),
                           CHF_BOX(region),
                           CHF_CONST_INT(ncomp));
    }
}
/***/
void
ViscousTensorOp::
prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                 const LevelData<FArrayBox>& a_correctCoarse)
{
  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  int mgref = 2; //this is a multigrid func

  if (s_prolongType > piecewiseConstant)
    {
      // need to set ghost cells for interpolation
      // note that all BC's are homogeneous, since we're
      // working with the correction.

      // as a first cut, do linear extrapolation everywhere,
      // followed by an exchange
      const DisjointBoxLayout& crseGrids = a_correctCoarse.getBoxes();

      // need to cast away const-ness in order to make this happen
      LevelData<FArrayBox>& crseCorr = *(const_cast<LevelData<FArrayBox>*>(&a_correctCoarse));

      DataIterator crseDit = a_correctCoarse.dataIterator();
      for (crseDit.begin(); crseDit.ok(); ++crseDit)
        {

          const Box crseBox = crseGrids[crseDit];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // don't bother if we're only one cell wide
              if (crseBox.size(dir) > 1)
                {
                  Box loBox = adjCellLo(crseBox, dir, 1);
                  int hiLo = -1;
                  FORT_LINEAREXTRAPVTOP(CHF_FRA(crseCorr[crseDit]),
                                        CHF_BOX(loBox),
                                        CHF_INT(dir),
                                        CHF_INT(hiLo));

                  Box hiBox = adjCellHi(crseBox,dir,1);
                  hiLo = 1;
                  FORT_LINEAREXTRAPVTOP(CHF_FRA(crseCorr[crseDit]),
                                        CHF_BOX(hiBox),
                                        CHF_INT(dir),
                                        CHF_INT(hiLo));
                }
            }
        }
      crseCorr.exchange();
    }

  for (DataIterator dit = a_phiThisLevel.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_phiThisLevel[dit];
      const FArrayBox& coarse = a_correctCoarse[dit];
      Box region = dbl.get(dit());
      Box cBox = coarsen(region, mgref);

      if (s_prolongType == piecewiseConstant)
        {
          FORT_PROLONGVTOP(CHF_FRA(phi),
                           CHF_CONST_FRA(coarse),
                           CHF_BOX(region),
                           CHF_CONST_INT(mgref));
        }
      else if (s_prolongType == linearInterp)
        {
          // first do piecewise-constant
          FORT_PROLONGVTOP(CHF_FRA(phi),
                           CHF_CONST_FRA(coarse),
                           CHF_BOX(region),
                           CHF_CONST_INT(mgref));

          // now compute slopes and increment
          // do we need to add limiting?

          int ncomp = coarse.nComp();
          FArrayBox dirSlopes(cBox, ncomp);
          for (int dir=0; dir<SpaceDim; dir++)
            {

              // don't bother if crse grid is only one cell wide
              if (cBox.size(dir) > 1)
                {

                  FORT_SLOPESVTOP(CHF_FRA(dirSlopes),
                                  CHF_CONST_FRA(coarse),
                                  CHF_BOX(cBox),
                                  CHF_INT(dir));

                  // copy code from PiecewiseLinearFillPatch for now

                  BoxIterator bit(region);
                  for (bit.begin(); bit.ok(); ++bit)
                    {
                      const IntVect& fine_iv = bit();
                      const IntVect crse_iv = coarsen(fine_iv,mgref);
                      const int offset = fine_iv[dir] - mgref*crse_iv[dir];
                      Real interp_coef = -.5 + (offset +.5) / mgref;
                      for (int comp=0; comp < ncomp; ++comp)
                        {
                          phi(fine_iv,comp) +=interp_coef*dirSlopes(crse_iv,comp);
                        }
                    }

#if 0
                  // fortran implementation, using fortran in AMRTools
                  // (DFM 8/26/11) -- in my timing tests, the BoxIterator
                  // version was actually faster; leaving the fortran here
                  // (but commented out) in case somebody later finds that
                  // not to be the case
                  Box refBox(IntVect::Zero,
                             (mgref-1)*IntVect::Unit);

                  // call fortran from AMRTools...
                  FORT_INTERPLINEAR(CHF_FRA(phi),
                                    CHF_CONST_FRA(dirSlopes),
                                    CHF_BOX(cBox),
                                    CHF_CONST_INT(dir),
                                    CHF_CONST_INT(mgref),
                                    CHF_BOX(refBox));
#endif
                } // end if crse box is more than one cell wide

            } // end loop over directions

        }
      else
        {
          MayDay::Error("ViscousTensorOp -- bad prolongation type");
        }
    }
}
/***/
void
ViscousTensorOp::
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
ViscousTensorOp::
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
      // set BCs for fillGrad
      LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phiFine;
      ViscousTensorOp *fineOp = (ViscousTensorOp*)a_finerOp;
      const DisjointBoxLayout& dbl = phi.disjointBoxLayout();
      DataIterator dit = phi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_bc(phi[dit], dbl[dit()], fineOp->m_domain, fineOp->m_dx, a_homogeneousPhysBC, dit());
        }
      
      // fillGrad called on a_phiFine
      reflux(a_phiFine, a_phi,  a_residual, a_finerOp);
    }
  incr(a_residual, a_rhs, -1.0);
  // residual is rhs - L(phi)
  scale (a_residual, -1.0);
}

/***/
void
ViscousTensorOp::
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
  // set BCs for fillGrad
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phiFine;
  ViscousTensorOp *fineOp = (ViscousTensorOp*)a_finerOp;
  const DisjointBoxLayout& dbl = phi.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()], fineOp->m_domain, fineOp->m_dx, a_homogeneousPhysBC, dit());
    }
  // fillGrad called on a_phiFine
  reflux(a_phiFine, a_phi,  a_residual, a_finerOp);
  axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}
/***/
void
ViscousTensorOp::
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
ViscousTensorOp::
AMRProlong(LevelData<FArrayBox>& a_correction,
           const LevelData<FArrayBox>& a_coarseCorrection)
{
  DisjointBoxLayout c;
  coarsen(c,  a_correction.disjointBoxLayout(), m_refToCoar);
  LevelData<FArrayBox> eCoar(c, a_correction.nComp(),a_coarseCorrection.ghostVect());
  a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

  Real dxCrse = m_dx*m_refToCoar;
  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  for (DataIterator dit = a_correction.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_correction[dit];
      FArrayBox& coarse = eCoar[dit];
      Box region = dbl.get(dit());

      // need this if we need to set physical BC's
      // on coarse-level correction.
      ProblemDomain crseDomain(m_domain);
      crseDomain.coarsen(m_refToCoar);

      if (s_prolongType == piecewiseConstant)
        {
          FORT_PROLONGVTOP(CHF_FRA(phi),
                           CHF_CONST_FRA(coarse),
                           CHF_BOX(region),
                           CHF_CONST_INT(m_refToCoar));
        }
      else if (s_prolongType == linearInterp)
        {
          // set homogeneous Physical BC's
          // (all other types are handled by coarse-level copyTo)
          Box cBox = coarsen(region, m_refToCoar);
          bool homogeneousBC=true;
          m_bc(coarse, cBox, crseDomain, dxCrse, homogeneousBC, dit());

          // first do piecewise-constant
          FORT_PROLONGVTOP(CHF_FRA(phi),
                           CHF_CONST_FRA(coarse),
                           CHF_BOX(region),
                           CHF_CONST_INT(m_refToCoar));

          // now compute slopes and increment
          // do we need to add limiting?

          int ncomp = coarse.nComp();
          FArrayBox dirSlopes(cBox, ncomp);
          for (int dir=0; dir<SpaceDim; dir++)
            {

              // don't bother if crse grid is only one cell wide
              if (cBox.size(dir) > 1)
                {

                  FORT_SLOPESVTOP(CHF_FRA(dirSlopes),
                                  CHF_CONST_FRA(coarse),
                                  CHF_BOX(cBox),
                                  CHF_INT(dir));

                  // copy code from PiecewiseLinearFillPatch for now

                  BoxIterator bit(region);
                  for (bit.begin(); bit.ok(); ++bit)
                    {
                      const IntVect& fine_iv = bit();
                      const IntVect crse_iv = coarsen(fine_iv,m_refToCoar);
                      const int offset = fine_iv[dir] - m_refToCoar*crse_iv[dir];
                      Real interp_coef = -.5 + (offset +.5) / m_refToCoar;
                      for (int comp=0; comp < ncomp; ++comp)
                        {
                          phi(fine_iv,comp) +=interp_coef*dirSlopes(crse_iv,comp);
                        }
                    }

#if 0
                  // fortran implementation, using fortran in AMRTools
                  // (DFM 8/26/11) -- in my timing tests, the BoxIterator
                  // version was actually faster; leaving the fortran here
                  // (but commented out) in case somebody later finds that
                  // not to be the case
                  Box refBox(IntVect::Zero,
                             (m_refToCoar-1)*IntVect::Unit);

                  // call fortran from AMRTools...
                  FORT_INTERPLINEAR(CHF_FRA(phi),
                                    CHF_CONST_FRA(dirSlopes),
                                    CHF_BOX(cBox),
                                    CHF_CONST_INT(dir),
                                    CHF_CONST_INT(m_refToCoar),
                                    CHF_BOX(refBox));
#endif
                } // end if crse box is more than one cell wide

            } // end loop over directions

        }
      else
        {
          MayDay::Error("ViscousTensorOp -- bad prolongation type");
        }

    }
}

/***/
void
ViscousTensorOp::
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
ViscousTensorOp::
createCoarsened(LevelData<FArrayBox>&       a_lhs,
                const LevelData<FArrayBox>& a_rhs,
                const int &                 a_refRat)
{
  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  DisjointBoxLayout dbl = a_rhs.disjointBoxLayout();
  CH_assert(dbl.coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, dbl, a_refRat);

  a_lhs.define(dblCoarsenedFine, ncomp, a_rhs.ghostVect());
}

/***/
void ViscousTensorOp::preCond(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs)
{
  //slc : down here, we don't want to change from the default behaviour
  m_relaxTolerance = 0.0;
  m_relaxMinIter = 40;
  relax(a_phi, a_rhs, 1);
}

/***/
void
ViscousTensorOp::
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
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous, dit());
    }

  //contains an exchange
  this->fillGrad(a_phi);

  computeOperatorNoBCs(a_lhs, phi);
}

/***/
void
ViscousTensorOp::
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
      m_bc(phi[dit], dbl[dit()],m_domain, m_dx, false, dit());
    }
  //exchange with neighboring boxes
  phi.exchange(phi.interval(), m_exchangeCopier);
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box  grid = dbl.get(dit());
      a_div[dit()].setVal(0.);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FORT_CELLDIVINCRVTOP(CHF_FRA1     (a_div[dit()], 0),
                               CHF_CONST_FRA(phi  [dit()]   ),
                               CHF_CONST_REAL(m_dx),
                               CHF_CONST_INT(idir),
                               CHF_BOX(grid));
        }
    }
}
void ViscousTensorOp::reflux(const LevelData<FArrayBox>& a_phiFine,
                             const LevelData<FArrayBox>& a_phi,
                             LevelData<FArrayBox>& residual,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  int ncomp = SpaceDim;
  ProblemDomain fineDomain = refine(m_domain, m_refToFine);
  LevelFluxRegister levfluxreg(a_phiFine.disjointBoxLayout(),
                               a_phi.disjointBoxLayout(),
                               fineDomain,
                               m_refToFine,
                               ncomp);

  levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab    =   a_phi    [dit];
      const FluxBox&   coareta    = (*m_eta)   [dit];
      const FluxBox&   coarlambda = (*m_lambda)[dit];
      const Box&       gridBox    =   a_phi.getBoxes()[dit];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox coarflux;
          Box faceBox = surroundingNodes(gridBox, idir);
          //fillgrad was called before this when applyOp was called
          getFlux(coarflux, coarfab, m_grad[dit()], coareta[idir], coarlambda[idir], faceBox, idir);

          Real scale = 1.0;
          levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                     interv,interv,idir);
        }
    }
  LevelData<FArrayBox>& p = ( LevelData<FArrayBox>&)a_phiFine;

  //interpolate finer phi and compute gradients.
  ViscousTensorOp* finerAMRPOp = (ViscousTensorOp*) a_finerOp;
  finerAMRPOp->cfinterp(p, a_phi);
  finerAMRPOp->fillGrad(p);

  IntVect phiGhost = p.ghostVect();

  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = a_phiFine[ditf];
      const FArrayBox& fineTanGrad = finerAMRPOp->m_grad[ditf];
      const FluxBox& fineeta    = (*(finerAMRPOp->m_eta))   [ditf];
      const FluxBox& finelambda = (*(finerAMRPOp->m_lambda))[ditf];
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

              FArrayBox phifab(fabbox, a_phi.nComp());
              phifab.copy(phifFab);

              FArrayBox fineflux;
              getFlux(fineflux, phifab, fineTanGrad,
                      fineeta[idir], finelambda[idir],
                      facebox, idir, m_refToFine);

              Real scale = 1.0;
              levfluxreg.incrementFine(fineflux, scale, ditf(),
                                       interv, interv, idir, hiorlo);
            }
        }
    }

  Real scale =  m_beta/m_dx;
  levfluxreg.reflux(residual, scale);
}
/****/

/****/
void
ViscousTensorOp::
loHiCenterFace(Box&                 a_loBox,
               int&                 a_hasLo,
               Box&                 a_hiBox,
               int&                 a_hasHi,
               Box&                 a_centerBox,
               const ProblemDomain& a_eblg,
               const Box&           a_inBox,
               const int&           a_dir)
{

  // if we're periodic in a_dir, then all boxes are
  // contained in the domain (DFM 2/5/10)
  if (a_eblg.isPeriodic(a_dir) )
    {
      a_centerBox = a_inBox;
      a_hasLo = 0;
      a_hasHi = 0;
    } // end if periodic in a_dir
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
    } // end if not periodic in a_dir
}

/***/
void
ViscousTensorOp::
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
      FORT_FACEDIVINCRVTOP(CHF_FRA1(a_faceDiv, 0),
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
      for (int velDir = 0; velDir < SpaceDim; velDir++)
        {
          int gradcomp = TensorCFInterp::gradIndex(velDir,divDir);
          FORT_GETFACEGRADVTOP(CHF_FRA1(a_faceGrad, gradcomp),
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
void ViscousTensorOp::getFlux(FArrayBox&       a_flux,
                              const FArrayBox& a_data,
                              const FArrayBox& a_gradData,
                              const FArrayBox& a_etaFace,
                              const FArrayBox& a_lambdaFace,
                              const Box&       a_faceBox,
                              int              a_dir,
                              int a_ref)
{
  ProblemDomain domain(m_domain);
  domain.refine(a_ref);
  Real dx(m_dx);

  dx /= a_ref;
  CH_assert(a_data.nComp() == m_ncomp);
  a_flux.resize(a_faceBox, a_data.nComp());
  FArrayBox  faceDiv(a_faceBox, 1);
  FArrayBox  faceGrad(a_faceBox, a_gradData.nComp());
  getFaceDivAndGrad(faceDiv, faceGrad, a_data, a_gradData, domain, a_faceBox, a_dir, dx);

  getFluxFromDivAndGrad(a_flux, faceDiv, faceGrad, a_etaFace, a_lambdaFace, a_faceBox, a_dir);
  a_flux *= m_beta;
}
void ViscousTensorOp::getFluxFromDivAndGrad(FArrayBox&       a_flux,
                                            const FArrayBox& a_faceDiv,
                                            const FArrayBox& a_faceGrad,
                                            const FArrayBox& a_etaFace,
                                            const FArrayBox& a_lambdaFace,
                                            const Box&       a_faceBox,
                                            int              a_dir)

{
  //copy the divergence into the diagonal component of the flux
  //diagonal component of the flux is the a_dir comp
  int dstcomp = a_dir;
  int srccomp = 0;
  int numcomp = 1;

  //this is important because we increment later
  a_flux.setVal(0.);

  //now make diag component of f = div v
  a_flux.copy(a_faceDiv, srccomp, dstcomp, numcomp);
  //the flux now contains div v * I

  //now flux gets multiplied by lambda
  a_flux.mult(a_lambdaFace, srccomp, dstcomp, 1);

  //the flux now contains div v * lambda * I
  //now add in eta*(grad v + (grad v)^T)
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int fluxComp = idir;
      int gradComp = TensorCFInterp::gradIndex(idir,a_dir);
      int tranComp = TensorCFInterp::gradIndex(a_dir, idir);
      FORT_ADDGRADTOFLUXVTOP(CHF_FRA(a_flux),
                             CHF_CONST_FRA1(a_etaFace,0),
                             CHF_CONST_INT(fluxComp),
                             CHF_CONST_FRA(a_faceGrad),
                             CHF_CONST_INT(gradComp),
                             CHF_CONST_INT(tranComp),
                             CHF_BOX(a_faceBox));
    }
}
/***/
void
ViscousTensorOp::
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
ViscousTensorOp::
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
ViscousTensorOp::
create(LevelData<FArrayBox>&       a_lhs,
       const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.create(a_lhs, a_rhs);
}
/***/
void
ViscousTensorOp::
assign(LevelData<FArrayBox>&       a_lhs,
       const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.assign(a_lhs, a_rhs);
}
/***/
Real
ViscousTensorOp::
dotProduct(const LevelData<FArrayBox>& a_1,
           const LevelData<FArrayBox>& a_2)
{
  return m_levelOps.dotProduct(a_1, a_2);
}
/***/
void
ViscousTensorOp::
incr( LevelData<FArrayBox>&       a_lhs,
      const LevelData<FArrayBox>& a_x,
      Real a_scale)
{
  m_levelOps.incr(a_lhs, a_x, a_scale);
}
/***/
void
ViscousTensorOp::
axby( LevelData<FArrayBox>&       a_lhs,
      const LevelData<FArrayBox>& a_x,
      const LevelData<FArrayBox>& a_y,
      Real a_a, Real a_b)
{
  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}
/***/
void
ViscousTensorOp::
scale(LevelData<FArrayBox>& a_lhs,
      const Real& a_scale)
{
  m_levelOps.scale(a_lhs, a_scale);
}
/***/
void
ViscousTensorOp::
setToZero(LevelData<FArrayBox>& a_lhs)
{
  m_levelOps.setToZero(a_lhs);
}

/***/
void
ViscousTensorOp::
relax(LevelData<FArrayBox>&       a_phi,
      const LevelData<FArrayBox>& a_rhs,
      int a_iterations)
{
  CH_TIME("ViscousTensorOp::relax");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

#define HUGE_NORM 1.0e+20
#define TINY_NORM 1.0e-20
  LevelData<FArrayBox> lphi;
  m_levelOps.create(lphi, a_rhs);
  LevelData<FArrayBox> prevPhi;

  if (m_relaxTolerance > TINY_NORM)
    m_levelOps.create(prevPhi, a_rhs);

  // do first red, then black passes
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  int whichIter = 0;
  bool done = false;

  while (whichIter < a_iterations && !done)
    {
      if (whichIter > m_relaxMinIter && m_relaxTolerance >  TINY_NORM)
        assignLocal(prevPhi,a_phi); // no point doing this if we aren't testing

      // Loop over all possibilities (in all dimensions)
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
                  FORT_GSRBVTOP(CHF_FRA(a_phi[dit]),
                                CHF_CONST_FRA( lphi[dit]),
                                CHF_CONST_FRA(a_rhs[dit]),
                                CHF_CONST_FRA(m_relaxCoef[dit]),
                                CHF_BOX(coloredBox),
                                CHF_CONST_INT(m_ncomp));
                }
            }
        } // end loop through red-black

      //slc : if we have done at least some smooths, give up if they
      //      are no longer changing the solution much. This is not
      //      fatal : it usually leads to multigrid happiness.

      if (whichIter > m_relaxMinIter && m_relaxTolerance >  TINY_NORM)
      {
        Real maxPhi =  norm(a_phi,0);
        CH_assert(maxPhi < HUGE_NORM);
        incr(prevPhi,a_phi,-1.0);
        Real maxDPhi = norm(prevPhi, 0);
        CH_assert(maxDPhi < HUGE_NORM);
        done = maxDPhi < m_relaxTolerance * (maxPhi + TINY_NORM);
      }
      whichIter++;
    }
}
/***/
Real
ViscousTensorOp::
norm(const LevelData<FArrayBox>& a_x, int a_ord)
{
  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}

/***/
ViscousTensorOp::
ViscousTensorOp(const DisjointBoxLayout&                    a_grids,
                const DisjointBoxLayout&                    a_gridsFine,
                const DisjointBoxLayout&                    a_gridsCoar,
                const RefCountedPtr<LevelData<FluxBox> >&   a_eta,
                const RefCountedPtr<LevelData<FluxBox> >&   a_lambda,
                const RefCountedPtr<LevelData<FArrayBox> >& a_acoef,
                Real                                        a_alpha,
                Real                                        a_beta,
                int                                         a_refToFine,
                int                                         a_refToCoar,
                const ProblemDomain&                        a_domain,
                const Real&                                 a_dx,
                const Real&                                 a_dxCrse,
                BCFunc                                      a_bc,
                Real                                        a_safety,
                Real                                        a_relaxTolerance,
                int                                         a_relaxMinIter)
{
  BCHolder bc(a_bc);
  define(a_grids, a_gridsFine, a_gridsCoar, a_eta, a_lambda,
         a_acoef, a_alpha, a_beta, a_refToFine, a_refToCoar,
         a_domain, a_dx, a_dxCrse, bc, a_safety, a_relaxTolerance,a_relaxMinIter);

}

/***/
ViscousTensorOp::
ViscousTensorOp(const DisjointBoxLayout&                    a_grids,
                const DisjointBoxLayout&                    a_gridsFine,
                const DisjointBoxLayout&                    a_gridsCoar,
                const RefCountedPtr<LevelData<FluxBox> >&   a_eta,
                const RefCountedPtr<LevelData<FluxBox> >&   a_lambda,
                const RefCountedPtr<LevelData<FArrayBox> >& a_acoef,
                Real                                        a_alpha,
                Real                                        a_beta,
                int                                         a_refToFine,
                int                                         a_refToCoar,
                const ProblemDomain&                        a_domain,
                const Real&                                 a_dx,
                const Real&                                 a_dxCrse,
                BCHolder&                                   a_bc,
                Real                                        a_safety,
                Real                                        a_relaxTolerance,
                int                                         a_relaxMinIter)
{
  define(a_grids, a_gridsFine, a_gridsCoar, a_eta, a_lambda,
         a_acoef, a_alpha, a_beta, a_refToFine, a_refToCoar,
         a_domain, a_dx, a_dxCrse, a_bc, a_safety, a_relaxTolerance,a_relaxMinIter);
}

void
ViscousTensorOp::
define(const DisjointBoxLayout&                    a_grids,
       const DisjointBoxLayout&                    a_gridsFine,
       const DisjointBoxLayout&                    a_gridsCoar,
       const RefCountedPtr<LevelData<FluxBox> >&   a_eta,
       const RefCountedPtr<LevelData<FluxBox> >&   a_lambda,
       const RefCountedPtr<LevelData<FArrayBox> >& a_acoef,
       Real                                        a_alpha,
       Real                                        a_beta,
       int                                         a_refToFine,
       int                                         a_refToCoar,
       const ProblemDomain&                        a_domain,
       const Real&                                 a_dx,
       const Real&                                 a_dxCrse,
       BCHolder&                                   a_bc,
       Real                                        a_safety,
       Real                                        a_relaxTolerance,
       int                                         a_relaxMinIter)

{
  vtogetMultiColors(m_colors);
  m_acoef         =  a_acoef;
  m_lambda        =  a_lambda;
  m_eta           =  a_eta;
  m_alpha         =  a_alpha;
  m_beta          =  a_beta;
  m_refToCoar     =  a_refToCoar;
  m_refToFine     =  a_refToFine;
  m_domain        =  a_domain;
  m_dx            =  a_dx;
  m_dxCrse        =  a_dxCrse;
  m_bc            =  a_bc;
  m_safety        =  a_safety;
  m_relaxTolerance =  a_relaxTolerance;
  m_relaxMinIter = a_relaxMinIter;
  m_ncomp = SpaceDim;

  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  if (a_gridsCoar.isClosed())
    {
      m_interpWithCoarser.define(a_grids, &a_gridsCoar, a_dx,
                                 m_refToCoar, SpaceDim, m_domain);
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
  m_relaxCoef.define(a_grids, SpaceDim,          IntVect::Zero);
  m_grad.define(     a_grids, SpaceDim*SpaceDim, IntVect::Unit);
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
  defineRelCoef();
  m_levelOps.setToZero(m_grad);
}
void
ViscousTensorOp::
defineRelCoef()
{

  DisjointBoxLayout grids = m_relaxCoef.disjointBoxLayout();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      const Box& grid = grids.get(dit());
      for (int ivar = 0; ivar < SpaceDim; ivar++)
        {
          int src = 0; int dst = ivar; int ncomp = 1;
          m_relaxCoef[dit()].copy((*m_acoef)[dit()], src,dst,ncomp);
        }
      m_relaxCoef[dit()] *= m_alpha;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FORT_DECRINVRELCOEFVTOP(CHF_FRA(m_relaxCoef[dit()]),
                                  CHF_FRA((*m_eta)[dit()][idir]),
                                  CHF_FRA((*m_lambda)[dit()][idir]),
                                  CHF_CONST_REAL(m_beta),
                                  CHF_BOX(grid),
                                  CHF_REAL(m_dx),
                                  CHF_INT(idir),
                                  CHF_INT(m_ncomp));
#if 0
          CH_assert(abs(m_relaxCoef[dit()].min(0)) > 1.0e-15
                    && abs(m_relaxCoef[dit()].max(0)) > 1.0e-15);
#endif

        }

      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)

      FORT_INVERTLAMBDAVTOP(CHF_FRA(m_relaxCoef[dit()]),
                            CHF_REAL(m_safety),
                            CHF_BOX(grid),
                            CHF_INT(m_ncomp));

      CH_assert(m_relaxCoef[dit()].norm(0) < 1.0e+15);

    }
}

/***/
ViscousTensorOpFactory::
ViscousTensorOpFactory(const Vector<DisjointBoxLayout>&                     a_grids,
                       const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_eta,
                       const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_lambda,
                       const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_acoef,
                       Real                                                 a_alpha,
                       Real                                                 a_beta,
                       const Vector<int>&                                   a_refRatios,
                       const ProblemDomain&                                 a_coarseDomain,
                       const Real&                                          a_coarseDx,
                       BCFunc                                               a_bc,
                       Real                                                 a_safety,
                       Real                                                 a_relaxTolerance,
                       int                                                  a_relaxMinIter)
{
  setDefaults();
  BCHolder bc(a_bc);
  define(a_grids,a_eta,a_lambda,a_acoef,a_alpha,a_beta,
         a_refRatios,a_coarseDomain,a_coarseDx,bc,a_safety,a_relaxTolerance, a_relaxMinIter);
}

ViscousTensorOpFactory::
ViscousTensorOpFactory(const Vector<DisjointBoxLayout>&                     a_grids,
                       const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_eta,
                       const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_lambda,
                       const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_acoef,
                       Real                                                 a_alpha,
                       Real                                                 a_beta,
                       const Vector<int>&                                   a_refRatios,
                       const ProblemDomain&                                 a_coarseDomain,
                       const Real&                                          a_coarseDx,
                       BCHolder                                            a_bc,
                       Real                                                 a_safety,
                       Real                                                 a_relaxTolerance,
                       int                                                  a_relaxMinIter)

{
  setDefaults();
  define(a_grids,a_eta,a_lambda,a_acoef,a_alpha,a_beta,
         a_refRatios,a_coarseDomain,a_coarseDx,a_bc,a_safety,a_relaxTolerance, a_relaxMinIter);
}

void
ViscousTensorOpFactory::
define(const Vector<DisjointBoxLayout>&                     a_grids,
       const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_eta,
       const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_lambda,
       const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_acoef,
       Real                                                 a_alpha,
       Real                                                 a_beta,
       const Vector<int>&                                   a_refRatios,
       const ProblemDomain&                                 a_coarseDomain,
       const Real&                                          a_coarseDx,
       BCHolder                                             a_bc,
       Real                                                 a_safety,
       Real                                                 a_relaxTolerance,
       int                                                  a_relaxMinIter)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;
  m_eta   = a_eta;
  m_acoef   = a_acoef;
  m_lambda = a_lambda;
  m_domains.resize(a_grids.size());
  m_boxes=a_grids;
  m_refRatios=a_refRatios;
  m_dx.resize(a_grids.size());
  m_bc = a_bc;
  m_safety = a_safety;
  m_relaxTolerance = a_relaxTolerance;
  m_relaxMinIter = a_relaxMinIter;
  m_domains[0] = a_coarseDomain;
  m_dx[0] = a_coarseDx;
  for (int i=1; i<a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1]/m_refRatios[i-1] ;
      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);
    }
}
/***/
ViscousTensorOp*
ViscousTensorOpFactory::
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
      DisjointBoxLayout dbl;
      coarsen_dbl(dbl, layout, 2);
      layout = dbl;
      dx*=2;
      refToDepth *= 2;
      domain.coarsen(2);
    }
  
  RefCountedPtr<LevelData<FluxBox> >     eta;
  RefCountedPtr<LevelData<FluxBox> >  lambda;
  RefCountedPtr<LevelData<FArrayBox> > acoef;

  if (a_depth == 0)
    {      
      // (DFM 5/11/12) copy refCountedPtr to existing storage 
      // instead of allocating new storage
      eta = m_eta[ref];
      lambda = m_lambda[ref];
      acoef = m_acoef[ref];
    }
  else
    {
      // allocated (coarsened) storage and do coarsening
      RefCountedPtr<LevelData<FluxBox> >     etaNew( new LevelData<FluxBox>  (layout, 1, IntVect::Zero) );
      RefCountedPtr<LevelData<FluxBox> >  lambdaNew( new LevelData<FluxBox>  (layout, 1, IntVect::Zero) );
      RefCountedPtr<LevelData<FArrayBox> > acoefNew( new LevelData<FArrayBox>(layout, 1, IntVect::Zero) );
      coarsenStuff(*etaNew, *lambdaNew,  *acoefNew, *m_eta[ref], *m_lambda[ref],  *m_acoef[ref], refToDepth, s_coefficientAverageType);

      eta = etaNew;
      lambda = lambdaNew;
      acoef = acoefNew;
    }
      //no coarser or finer grids
      //no reftocoar, reftofine
  Real dxCrse = 2*m_dx[ref];
  if (ref > 0)
    {
      dxCrse = m_dx[ref-1];
    }
  ViscousTensorOp* newOp = new ViscousTensorOp(layout, DisjointBoxLayout(), DisjointBoxLayout(),
                                               eta, lambda, acoef, m_alpha, m_beta, -1, -1,
                                               domain,  dx, dxCrse, m_bc, m_safety, m_relaxTolerance, m_relaxMinIter);

  return newOp;
}
void
ViscousTensorOp::
AMROperatorNC(LevelData<FArrayBox>&       a_LofPhi,
              const LevelData<FArrayBox>& a_phiFine,
              const LevelData<FArrayBox>& a_phi,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{

  //no coarse-fine interpolation here
  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);

  // set BCs for fillGrad
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phiFine;
  ViscousTensorOp *fineOp = (ViscousTensorOp*)a_finerOp;
  const DisjointBoxLayout& dbl = phi.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()], fineOp->m_domain, fineOp->m_dx, a_homogeneousPhysBC, dit());
    }
  
  // fillGrad called on a_phiFine
  reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
}
void
ViscousTensorOp::
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
ViscousTensorOp::
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
      // set BCs for fillGrad
      LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phiFine;
      ViscousTensorOp *fineOp = (ViscousTensorOp*)a_finerOp;
      const DisjointBoxLayout& dbl = phi.disjointBoxLayout();
      DataIterator dit = phi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_bc(phi[dit], dbl[dit()], fineOp->m_domain, fineOp->m_dx, a_homogeneousPhysBC, dit());
        }
      
      // fillGrad called on a_phiFine
      reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
    }
}
/***/
ViscousTensorOp*
ViscousTensorOpFactory::
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
  // (DFM 3/22/10) -- set this to a bogus value to make sure it isn't
  // used anywhere unless we actually set it in the case where a finer
  // level exists
  int refToFine = -1;
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
      refToFine = m_refRatios[ref];
    }

  ViscousTensorOp* newOp = new ViscousTensorOp(m_boxes[ref], dblFine, dblCoar,
                                               m_eta[ref],   m_lambda[ref], m_acoef[ref],
                                               m_alpha, m_beta,
                                               refToFine, refToCoar,
                                               m_domains[ref],  m_dx[ref],
                                               dxCrse, m_bc, m_safety, m_relaxTolerance, m_relaxMinIter);

  return newOp;
}
/***/
int
ViscousTensorOpFactory::
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

void
ViscousTensorOpFactory::setDefaults()
{
  //s_coefficientAverageType = CoarseAverageFace::harmonic;
  //s_coefficientAverageType = CoarseAverageFace::arithmetic;
}

/**/
void
ViscousTensorOp::
cellGrad(FArrayBox&             a_gradPhi,
         const  FArrayBox&      a_phi,
         const Box&             a_grid)
{
  CH_assert(a_gradPhi.nComp() == SpaceDim*SpaceDim);
  CH_assert(a_phi.nComp() == SpaceDim);

  for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      for (int phiDir = 0; phiDir < SpaceDim; phiDir++)
        {
          int gradcomp = TensorCFInterp::gradIndex(phiDir,derivDir);
          FORT_CELLGRADVTOP(CHF_FRA1(a_gradPhi, gradcomp),
                            CHF_CONST_FRA1(a_phi, phiDir),
                            CHF_BOX(a_grid),
                            CHF_CONST_REAL(m_dx),
                            CHF_CONST_INT(derivDir));

        }
    }
}
/**/
void
ViscousTensorOp::
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
ViscousTensorOp::
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
ViscousTensorOp::
homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{

  QuadCFInterp::homogeneousCFInterp(a_phif, m_grad,
                                    m_loCFIVS, m_hiCFIVS,
                                    m_dx, m_dxCrse, m_ncomp,
                                    m_loTanStencilSets,m_hiTanStencilSets);
}

/***/
#include "NamespaceFooter.H"
