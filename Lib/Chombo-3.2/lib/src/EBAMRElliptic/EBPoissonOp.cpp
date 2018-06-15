#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "LoadBalance.H"
#include "EBArith.H"

#include "EBPoissonOp.H"

#include "EBAMRPoissonOpF_F.H"
#include "InterpF_F.H"
#include "AMRPoissonOpF_F.H"
#include "CH_Timer.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

//////////////
EBPoissonOp::~EBPoissonOp()
{
}
EBPoissonOp::EBPoissonOp()
{
}

//////////////
EBPoissonOp::
EBPoissonOp(   const EBLevelGrid &                  a_eblg,
               const EBLevelGrid &                  a_eblgCoarMG,
               const RefCountedPtr<BaseDomainBC>&   a_domainBC,
               const RefCountedPtr<BaseEBBC>&       a_ebBC,
               const RealVect&                      a_dx,
               const RealVect&                      a_origin,
               const bool&                          a_hasMGObjects,
               const int&                           a_numPreCondIters,
               const int&                           a_relaxType,
               const int&                           a_orderEB,
               const Real&                          a_alpha,
               const Real&                          a_beta,
               const IntVect&                       a_ghostCellsPhi,
               const IntVect&                       a_ghostCellsRHS)
  : m_ghostCellsPhi( a_ghostCellsPhi ),
    m_ghostCellsRHS( a_ghostCellsRHS )
{
  EBArith::getMultiColors(m_colors);
  m_numPreCondIters= a_numPreCondIters;
  m_relaxType      = a_relaxType;
  m_orderEB        = a_orderEB;
  m_eblg            = a_eblg;
  m_domainBC       = a_domainBC;
  m_ebBC           = a_ebBC;
  m_dx             = a_dx;
  m_origin         = a_origin;
  m_alpha          = a_alpha;
  m_beta           = a_beta;
  m_time           = 0.0;

  //pre-compute 1/dx and 1/(dx^2)
  m_invDx  = 1.0/m_dx;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_invDx2[idir] = m_invDx[idir]*m_invDx[idir];
    }

  //dxScale is used for unscaling the boundary area (see GeometryShop)
  // (it would be better if we could just use m_dx[0])
  Real maxDx = 0.0;
  Real dv    = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      dv *= m_dx[idir];
      if ( m_dx[idir] > maxDx )
        {
          maxDx = m_dx[idir];
        }
    }
  m_dxScale = pow(maxDx,SpaceDim-1)/dv;

  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects
  //in this case, we have no coarser levels because this is a single
  //level operator
  m_hasMGObjects = a_hasMGObjects;
  if (m_hasMGObjects)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      int ncomp = 1;
      m_ebInterpMG.define( m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(),   m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain(), mgRef, ncomp, m_eblg.getEBIS(),
                           m_ghostCellsPhi);

      m_ebAverageMG.define(m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain() , mgRef, ncomp, m_eblg.getEBIS(),
                           m_ghostCellsRHS);

    }

  //define stencils for the operator
  defineStencils();
}

/******/
void
EBPoissonOp::
defineStencils()
{
  CH_TIMERS("EBPoissonOp::defineStencils");
  CH_TIMER("opStencil_define", t1);
  CH_TIMER("colorStencil_define", t2);

  // define ebstencil for irregular applyOp
  m_opEBStencil.define(m_eblg.getDBL());
  // create vofstencils for applyOp
  LayoutData<BaseIVFAB<VoFStencil> > opStencil;
  opStencil.define(m_eblg.getDBL());
  //define bc stencils and create flux stencil
  m_ebBC->define((*m_eblg.getCFIVS()), m_dxScale*m_beta);
  LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);
  //create and define colored stencils (2 parts)
  LayoutData<BaseIVFAB<VoFStencil> > colorStencil;
  colorStencil.define(m_eblg.getDBL());
  LayoutData<BaseIVFAB<VoFStencil> > rhsColorStencil;
  rhsColorStencil.define(m_eblg.getDBL());

  m_alphaDiagWeight.define(   m_eblg.getDBL());
  m_betaDiagWeight.define(   m_eblg.getDBL());
  m_vofItIrreg.define( m_eblg.getDBL()); // vofiterator cache

  Box domainBox = m_eblg.getDomain().domainBox();
  Box sideBoxLo[SpaceDim];
  Box sideBoxHi[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_vofItIrregDomLo[idir].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofItIrregDomHi[idir].define( m_eblg.getDBL()); // vofiterator cache for domain hi
    }

  CH_START(t1);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& curBox = m_eblg.getDBL().get(dit());
      const EBISBox& curEBISBox = m_eblg.getEBISL()[dit()];
      const EBGraph& curEBGraph = curEBISBox.getEBGraph();

      IntVectSet notRegular = curEBISBox.getIrregIVS  (curBox);

      BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = opStencil[dit()];
      BaseIVFAB<Real>&       alphaWeight  = m_alphaDiagWeight[dit()];
      BaseIVFAB<Real>&        betaWeight  = m_betaDiagWeight[dit()];
      curStencilBaseIVFAB.define(notRegular,curEBGraph, 1);
      alphaWeight.define(        notRegular,curEBGraph, 1);
      betaWeight.define(         notRegular,curEBGraph, 1);

      //cache the vofIterators
      m_vofItIrreg[dit()].define(notRegular,curEBISBox.getEBGraph());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet loIrreg = notRegular;
          IntVectSet hiIrreg = notRegular;
          loIrreg &= sideBoxLo[idir];
          hiIrreg &= sideBoxHi[idir];
          m_vofItIrregDomLo[idir][dit()].define(loIrreg,curEBISBox.getEBGraph());
          m_vofItIrregDomHi[idir][dit()].define(hiIrreg,curEBISBox.getEBGraph());
        }

      //Operator vofstencil
      VoFIterator& vofit = m_vofItIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();

          VoFStencil& curStencil = curStencilBaseIVFAB(VoF,0);
          getOpVoFStencil(curStencil,curEBISBox,VoF);

          Real& curAlphaWeight  = alphaWeight(VoF,0);
          Real& curBetaWeight  =   betaWeight(VoF,0);

          const Real kappa = curEBISBox.volFrac(VoF);

          curAlphaWeight = kappa;
          curBetaWeight = 0.0;
          for (int i = 0; i < curStencil.size(); i++)
            {
              if (curStencil.vof(i) == VoF)
                {
                  curBetaWeight += curStencil.weight(i);
                  break;
                }
            }

          curStencil *= m_beta;
          if (m_alpha != 0)
            {
              curStencil.add(VoF, kappa*m_alpha);
            }

          const IntVect& iv = VoF.gridIndex();

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Box loSide = bdryLo(m_eblg.getDomain(),idir);
              loSide.shiftHalf(idir,1);
              if (loSide.contains(iv))
                {
                  Real faceAreaFrac = 0.0;
                  Vector<FaceIndex> faces = curEBISBox.getFaces(VoF,idir,Side::Lo);
                  for (int i = 0; i < faces.size(); i++)
                    {
                      faceAreaFrac += curEBISBox.areaFrac(faces[i]);
                    }
                  curBetaWeight += -faceAreaFrac * m_invDx2[idir];
                }

              Box hiSide = bdryHi(m_eblg.getDomain(),idir);
              hiSide.shiftHalf(idir,-1);
              if (hiSide.contains(iv))
                {
                  Real faceAreaFrac = 0.0;
                  Vector<FaceIndex> faces = curEBISBox.getFaces(VoF,idir,Side::Hi);
                  for (int i = 0; i < faces.size(); i++)
                    {
                      faceAreaFrac += curEBISBox.areaFrac(faces[i]);
                    }
                  curBetaWeight += -faceAreaFrac * m_invDx2[idir];
                }
            }

          if (curBetaWeight == 0.0)
            {
              curBetaWeight = -1.0;
            }
          if (fluxStencil != NULL)
            {
              BaseIVFAB<VoFStencil>& fluxStencilBaseIVFAB = (*fluxStencil)[dit()];
              VoFStencil& fluxStencilPt = fluxStencilBaseIVFAB(VoF,0);
              curStencil += fluxStencilPt;
            }
        }//vofit

      //Operator ebstencil
      m_opEBStencil[dit()] = RefCountedPtr<EBStencil>(new EBStencil(m_vofItIrreg[dit()].getVector(), opStencil[dit()], m_eblg.getDBL().get(dit()), m_eblg.getEBISL()[dit()], m_ghostCellsPhi, m_ghostCellsRHS));

    }//dit
  CH_STOP(t1);

  //color vofstencils and ebstencils
  IntVect color = IntVect::Zero;
  IntVect limit = IntVect::Unit;
  color[0]=-1;
  // Loop over all possibilities (in all dimensions)
  CH_START(t2);
  for (int icolor = 0; icolor < m_colors.size(); icolor++)
    {
      m_colorEBStencil[icolor].define(m_eblg.getDBL());
      m_rhsColorEBStencil[icolor].define(m_eblg.getDBL());
      m_vofItIrregColor[icolor].define( m_eblg.getDBL());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_vofItIrregColorDomLo[icolor][idir].define( m_eblg.getDBL());
          m_vofItIrregColorDomHi[icolor][idir].define( m_eblg.getDBL());
        }
      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          IntVectSet ivsColor;
          const EBISBox& curEBISBox = m_eblg.getEBISL()[dit()];
          const EBGraph& curEBGraph = curEBISBox.getEBGraph();
          Box dblBox( m_eblg.getDBL().get(dit()) );

          BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = opStencil[dit()];
          BaseIVFAB<VoFStencil>& colorStencilBaseIVFAB = colorStencil[dit()];
          BaseIVFAB<VoFStencil>& rhsColorStencilBaseIVFAB = rhsColorStencil[dit()];

          const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[dit()];
          const BaseIVFAB<Real>& curBetaWeight  = m_betaDiagWeight[dit()];

          VoFIterator& vofit = m_vofItIrreg[dit()];

          int vofOrdinal = 0;
          for (vofit.reset(); vofit.ok(); ++vofit, ++vofOrdinal)
            {
              const VolIndex& vof = vofit();
              const IntVect& iv = vof.gridIndex();

              bool doThisVoF = true;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (iv[idir] % 2 != color[idir])
                    {
                      doThisVoF = false;
                      break;
                    }
                }

              if (doThisVoF)
                {
                  ivsColor |= iv;
                }
            }

          m_vofItIrregColor[icolor][dit()].define(ivsColor, curEBGraph);
          colorStencilBaseIVFAB.define(ivsColor, curEBGraph, 1);
          rhsColorStencilBaseIVFAB.define(ivsColor, curEBGraph, 1);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              IntVectSet loIrregColor = ivsColor;
              IntVectSet hiIrregColor = ivsColor;
              loIrregColor &= sideBoxLo[idir];
              hiIrregColor &= sideBoxHi[idir];
              m_vofItIrregColorDomLo[icolor][idir][dit()].define(loIrregColor,curEBISBox.getEBGraph());
              m_vofItIrregColorDomHi[icolor][idir][dit()].define(hiIrregColor,curEBISBox.getEBGraph());
            }

          VoFIterator& vofitcolor = m_vofItIrregColor[icolor][dit()];
          for (vofitcolor.reset(); vofitcolor.ok(); ++vofitcolor)
            {
              const VolIndex& vof = vofitcolor();

              VoFStencil& curStencil = curStencilBaseIVFAB(vof,0);
              VoFStencil& colorStencil = colorStencilBaseIVFAB(vof,0);
              VoFStencil& rhsColorStencil = rhsColorStencilBaseIVFAB(vof,0);
              Real weightIrreg = m_alpha*curAlphaWeight(vof,0) + m_beta*curBetaWeight(vof,0);
              colorStencil = curStencil;
              colorStencil *= (-1.0/weightIrreg);
              colorStencil.add(vof, 1.0);
              rhsColorStencil.add(vof, 1.0/weightIrreg);
            }

          Vector<VolIndex> srcVofs = m_vofItIrregColor[icolor][dit()].getVector();

          //color ebstencils
          m_colorEBStencil[icolor][dit()]    = RefCountedPtr<EBStencil>(new EBStencil(srcVofs, colorStencil[dit()]   ,  m_eblg.getDBL().get(dit()), m_eblg.getEBISL()[dit()], m_ghostCellsPhi, m_ghostCellsPhi));
          m_rhsColorEBStencil[icolor][dit()] = RefCountedPtr<EBStencil>(new EBStencil(srcVofs, rhsColorStencil[dit()]     , m_eblg.getDBL().get(dit()) , m_eblg.getEBISL()[dit()], m_ghostCellsRHS, m_ghostCellsPhi));

        }//dit
    }//color
  CH_STOP(t2);
}

void EBPoissonOp::
residual(LevelData<EBCellFAB>&       a_residual,
         const LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         bool                        a_homogeneousPhysBC)
{
  CH_TIME("EBPoissonOp::residual");
  //this is a multigrid operator so only homogeneous CF BC
  //and null coar level
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  DataIterator dit = a_residual.dataIterator();
  applyOp(a_residual,a_phi, a_homogeneousPhysBC, dit );
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}



void EBPoissonOp::
preCond(LevelData<EBCellFAB>&       a_lhs,
        const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBPoissonOp::preCond");

  // Recall that the operator is: alpha*phi + beta*lap(phi)
  // For isotropic-dx Poisson: alpha=0,beta=1 and the
  //    diagonal term of this operator is: 4/h/h in 2D, 6/h/h in 3D,
  // For anisotropic-dx Helmholtz: alpha*phi + beta*lap(phi) and the
  //    diagonal term of this operator is: alpha - beta*(2/dx/dx + 2/dy/dy [+ 2/dz/dz])
  // Inverse of the diagonal term is our initial multiplier.

  getInvDiagRHS(a_lhs,a_rhs);
  //relax(a_lhs, a_rhs, m_numPreCondIters);
  for (int iter = 0; iter < m_numPreCondIters; iter++)
    {
      levelJacobi(a_lhs, a_rhs);
    }

}

void EBPoissonOp::
getInvDiagRHS(LevelData<EBCellFAB>&       a_lhs,
              const LevelData<EBCellFAB>& a_rhs)
{
  //this function computes: a_lhs = (1/diagonal)*a_rhs
  // use this to initialize the preconditioner

  Real mult = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      mult += -2.0 * m_beta * m_invDx2[idir];
    }
  mult = 1.0 / mult;

  //do all cells (assuming mult works for every cell)
  assign(a_lhs,a_rhs);
  scale(a_lhs, mult);

  //now do/fix the irreg cells
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&                   lhs = a_lhs[dit()];
      const EBCellFAB&             rhs = a_rhs[dit()];
      const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[dit()];
      const BaseIVFAB<Real>& curBetaWeight =  m_betaDiagWeight[dit()];

      VoFIterator& vofit = m_vofItIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          Real weightIrreg = m_alpha*curAlphaWeight(VoF,0) + m_beta*curBetaWeight(VoF,0);
          for (int comp = 0; comp < lhs.nComp(); comp++)
            {
              lhs(VoF,comp) = (1.0/weightIrreg)*rhs(VoF,comp);
            }
        }
    }
}

void EBPoissonOp::
applyOp(LevelData<EBCellFAB>&             a_opPhi,
        const LevelData<EBCellFAB>&       a_phi,
        bool                              a_homogeneousPhysBC,
        DataIterator& dit,
        bool a_do_exchange)
{
  CH_TIMERS("EBPoissonOp::applyOp");
  CH_TIMER("regular_apply", t1);
  CH_TIMER("irregular_apply", t2);
  CH_TIMER("eb_bcs_apply", t3);
  CH_TIMER("dom_bcs_apply", t4);
  CH_TIMER("alpha_apply", t5);
  CH_assert(a_opPhi.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);


  CH_assert(a_phi.nComp() == a_opPhi.nComp());

  LevelData<EBCellFAB>& phi = const_cast<LevelData<EBCellFAB>&>(a_phi);
  if (a_do_exchange)
    {
      phi.exchange(phi.interval());
    }

  int nComps = a_phi.nComp();

  for (dit.reset(); dit.ok(); ++dit)
    {

      Box dblBox( m_eblg.getDBL().get(dit()) );
      const EBCellFAB& phifab = phi[dit()];
      Box curPhiBox = phifab.box();
      const BaseFab<Real>& curPhiFAB = phifab.getSingleValuedFAB();

      EBCellFAB& opphifab = a_opPhi[dit()];
      BaseFab<Real>& curOpPhiFAB = opphifab.getSingleValuedFAB();

      //       Interval interv(0, nComps-1);
      CH_START(t5);
      if (m_alpha == 0)
        {
          opphifab.setVal(0.0);
        }
      else
        {
          opphifab.copy(phifab);
          opphifab.mult(m_alpha);
        }
      CH_STOP(t5);

      Box loBox[SpaceDim],hiBox[SpaceDim];
      int hasLo[SpaceDim],hasHi[SpaceDim];
      CH_START(t1);
      applyOpRegularAllDirs( loBox, hiBox, hasLo, hasHi,
                             dblBox, curPhiBox, nComps,
                             curOpPhiFAB,
                             curPhiFAB,
                             a_homogeneousPhysBC,
                             dit(),
                             m_beta);

      CH_STOP(t1);

      CH_START(t2);
      //apply stencil
      m_opEBStencil[dit()]->apply(opphifab, phifab, false);
      CH_STOP(t2);
      //apply inhom boundary conditions
      CH_START(t3);
      if (!a_homogeneousPhysBC)
        {
          const Real factor = m_dxScale*m_beta;
          m_ebBC->applyEBFlux(opphifab, phifab, m_vofItIrreg[dit()], (*m_eblg.getCFIVS()),
                              dit(), m_origin, m_dx, factor,
                              a_homogeneousPhysBC, m_time);
        }
      CH_STOP(t3);

      CH_START(t4);
      int comp = 0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (m_vofItIrregDomLo[idir][dit()].reset(); m_vofItIrregDomLo[idir][dit()].ok();  ++m_vofItIrregDomLo[idir][dit()])
            {
              Real flux;
              const VolIndex& vof = m_vofItIrregDomLo[idir][dit()]();
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi[dit()],
                                      m_origin,m_dx,idir,Side::Lo, dit(), m_time,
                                      a_homogeneousPhysBC);

              opphifab(vof,comp) -= flux * m_beta*m_invDx[idir];
            }
          for (m_vofItIrregDomHi[idir][dit()].reset(); m_vofItIrregDomHi[idir][dit()].ok();  ++m_vofItIrregDomHi[idir][dit()])
            {
              Real flux;
              const VolIndex& vof = m_vofItIrregDomHi[idir][dit()]();
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi[dit()],
                                      m_origin,m_dx,idir,Side::Hi,dit(), m_time,
                                      a_homogeneousPhysBC);

              opphifab(vof,comp) += flux * m_beta*m_invDx[idir];
            }

        }
      CH_STOP(t4);
    }
}

void EBPoissonOp::
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
}

void EBPoissonOp::
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  EBLevelDataOps::assign(a_lhs,a_rhs);
}

Real EBPoissonOp::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
}

void EBPoissonOp::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}

void EBPoissonOp::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}

void EBPoissonOp::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  EBLevelDataOps::scale(a_lhs,a_scale);
}

Real EBPoissonOp::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaNorm(volume,a_rhs,EBLEVELDATAOPS_ALLVOFS,domain,a_ord);
}

void EBPoissonOp::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  EBLevelDataOps::setToZero(a_lhs);
}

void EBPoissonOp::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  EBLevelDataOps::setVal(a_lhs, a_value);
}

// MGLevelOp functions

void EBPoissonOp::
createCoarser(LevelData<EBCellFAB>&       a_coar,
              const LevelData<EBCellFAB>& a_fine,
              bool                        a_ghosted)
{
  CH_assert(a_fine.nComp() == 1);
  const DisjointBoxLayout& dbl = m_eblgCoarMG.getDBL();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), 2);

  int nghost = a_fine.ghostVect()[0];
  EBISLayout coarEBISL;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(coarEBISL,
                          dbl, coarDom, nghost);

  EBCellFactory ebcellfact(coarEBISL);
  a_coar.define(dbl, 1,a_fine.ghostVect(),ebcellfact);
}

void EBPoissonOp::
relax(LevelData<EBCellFAB>&       a_e,
      const LevelData<EBCellFAB>& a_residual,
      int                         a_iterations)
{
  CH_TIME("EBPoissonOp::relax");

  CH_assert(a_e.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);

  CH_assert(a_e.nComp() == 1);
  CH_assert(a_residual.nComp() == 1);

  if (m_relaxType == 0)
    {
      for (int i = 0; i < a_iterations; i++)
        {
          levelJacobi(a_e,a_residual);
        }
    }
  else  if (m_relaxType == 1)
    {
      for (int i = 0; i < a_iterations; i++)
        {
          levelMulticolorGS(a_e,a_residual);
        }
    }
  else
    {
      MayDay::Error("EBPoissonOp::relax - invalid relaxation type");
    }
}
/****/
void EBPoissonOp::
restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                 LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_rhsThisLevel)
{
  CH_TIME("EBPoissonOp::restrictResidual");

  CH_assert(a_resCoar.nComp() == 1);
  CH_assert(a_phiThisLevel.nComp() == 1);
  CH_assert(a_rhsThisLevel.nComp() == 1);

  LevelData<EBCellFAB> resThisLevel;
  bool homogeneous = true;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_rhsThisLevel.ghostVect();

  resThisLevel.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

  // Get the residual on the fine grid
  residual(resThisLevel,a_phiThisLevel,a_rhsThisLevel,homogeneous);

  // now use our nifty averaging operator
  Interval variables(0, 0);
  CH_assert(m_hasMGObjects);
  m_ebAverageMG.average(a_resCoar, resThisLevel, variables);

#ifdef DO_EB_RHS_CORRECTION
  // Apply error-correction modification to restricted RHS
  // Right now this only works with Dirichlet BC's, and only makes sense for the EB

  int correctionType = 0;    // Change this to activate RHS correction

  if (correctionType != 0)
    {
      for (DataIterator dit = a_resCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB&       resFAB = a_resCoar[dit()];                      // Extract the FAB for this box
          const EBISBox&   ebis   = resFAB.getEBISBox();                   // Get the set of all IntVect indices

          // Iterate over the parts of the RHS corresponding to the irregular cells,
          // correcting the RHS in each cell
          VoFIterator ebvofit(ebis.getIrregIVS(ebis.getRegion()), ebis.getEBGraph());
          for (ebvofit.reset(); ebvofit.ok(); ++ebvofit)
            {
              for (int icomp = 0; icomp < resFAB.nComp(); ++icomp)    // For each component of the residual on this VoF
                {
                  if (correctionType == 1)
                    {
                      // Setting the residual to zero on the irregular cells gives convergence,
                      // though the rate isn't as good as the full correction
                      resFAB(ebvofit(), icomp) = 0.0;
                    }
                  else if (correctionType == 2)
                    {
                      // Kludge valid only for pseudo-1D test case.  May not work for you.
                      Real kappa = ebis.volFrac(ebvofit());
                      kappa = (kappa > 0.5 ? kappa : 0.5);    // Floor on kappa to prevent dividing by a tiny number
                      Real rhoCoeff = (kappa + 0.5)/(2.0*kappa);
                      resFAB(ebvofit(), icomp) *= (1.0 - rhoCoeff);
                    }
                  // else silently do nothing
                }
            }
        }
    }
#endif
}

/****/
void EBPoissonOp::
prolongIncrement(LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_correctCoar)
{
  CH_TIME("EBPoissonOp::prolongIncrement");
  Interval vars(0, 0);
  CH_assert(m_hasMGObjects);
  m_ebInterpMG.pwcInterp(a_phiThisLevel, a_correctCoar, vars);
}
///////////
void EBPoissonOp::
getOpVoFStencil(VoFStencil&     a_stencil,
                const EBISBox&  a_ebisbox,
                const VolIndex& a_VoF)
{

  a_stencil.clear();

  Vector<VolIndex> allMonotoneVoFs;
  int radius = 1;
  IntVect timesMoved = IntVect::Zero;
  IntVect pathSign   = IntVect::Zero;

  EBArith::getAllVoFsInMonotonePath(allMonotoneVoFs,timesMoved,pathSign,
                                    a_VoF,a_ebisbox,radius);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      getOpVoFStencil(a_stencil,
                      idir,
                      allMonotoneVoFs,
                      a_ebisbox,
                      a_VoF,
                      false);
    }
}

void EBPoissonOp::
getOpVoFStencil(VoFStencil&             a_stencil,
                const int&              a_idir,
                const Vector<VolIndex>& a_allMonotoneVoFs,
                const EBISBox&          a_ebisbox,
                const VolIndex&         a_VoF,
                const bool&             a_lowOrder)
{

  for (SideIterator sit; sit.ok(); ++sit)
    {
      Side::LoHiSide side = sit();
      Vector<FaceIndex> faces;

      faces = a_ebisbox.getFaces(a_VoF,a_idir,side);

      for (int iface = 0; iface < faces.size(); iface++)
        {
          FaceIndex face = faces[iface];
          VoFStencil faceStencil;

          getOpFaceStencil(faceStencil,a_allMonotoneVoFs,a_ebisbox,a_VoF,
                           a_idir,side,face,a_lowOrder);

          a_stencil += faceStencil;
        }
    }
}

void EBPoissonOp::
getOpFaceStencil(VoFStencil&             a_stencil,
                 const Vector<VolIndex>& a_allMonotoneVoFs,
                 const EBISBox&          a_ebisbox,
                 const VolIndex&         a_VoF,
                 int                     a_dir,
                 const Side::LoHiSide&   a_side,
                 const FaceIndex&        a_face,
                 const bool&             a_lowOrder)
{
  a_stencil.clear();

  const RealVect& faceCentroid = a_ebisbox.centroid(a_face);
  const IntVect& origin = a_VoF.gridIndex();

  IntVect dirs = IntVect::Zero;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != a_dir)
        {
          IntVect ivPlus = origin;
          ivPlus[idir] += 1;
          if (faceCentroid[idir] < 0.0)
            {
              dirs[idir] = -1;
            }
          else if (faceCentroid[idir] > 0.0)
            {
              dirs[idir] = 1;
            }
          else if (m_eblg.getDomain().contains(ivPlus))
            {
              dirs[idir] = 1;
            }
          else
            {
              dirs[idir] = -1;
            }
        }
    }

  bool orderOne = true;

  if (m_orderEB == 0 || a_lowOrder)
    {
      orderOne = false;
    }
  else
    {
      IntVect loVect = dirs;
      loVect.min(IntVect::Zero);

      IntVect hiVect = dirs;
      hiVect.max(IntVect::Zero);

      Box fluxBox(loVect,hiVect);
      IntVectSet fluxIVS(fluxBox);

      IVSIterator ivsit(fluxIVS);
      for (ivsit.begin(); ivsit.ok(); ++ivsit)
        {
          bool curOkay;

          IntVect ivDelta = ivsit();

          IntVect iv1 = ivDelta;
          iv1 += origin;

          VolIndex VoF1;

          curOkay = EBArith::isVoFHere(VoF1,a_allMonotoneVoFs,iv1);

          IntVect iv2 = iv1;
          iv2[a_dir] += sign(a_side);

          VolIndex VoF2;

          curOkay = curOkay && EBArith::isVoFHere(VoF2,a_allMonotoneVoFs,iv2);

          orderOne = orderOne && curOkay;

          if (curOkay)
            {
              VoFStencil curPair;
              curPair.add(VoF1,-m_invDx2[a_dir]);
              curPair.add(VoF2, m_invDx2[a_dir]);

              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (idir != a_dir)
                    {
                      if (ivDelta[idir] == 0)
                        {
                          curPair *= 1 - abs(faceCentroid[idir]);
                        }
                      else
                        {
                          curPair *= abs(faceCentroid[idir]);
                        }
                    }
                }
              a_stencil += curPair;
            }
        }
    }

  if (!orderOne)
    {
      VolIndex VoF1 = a_VoF;
      VolIndex VoF2 = a_face.getVoF(a_side);

      if (a_ebisbox.getRegion().contains(VoF2.gridIndex()))
        {
          a_stencil.clear();

          a_stencil.add(VoF1,-m_invDx2[a_dir]);
          a_stencil.add(VoF2, m_invDx2[a_dir]);
        }
    }
  if (!a_lowOrder)
    {
      a_stencil *= a_ebisbox.areaFrac(a_face);
    }
}

void EBPoissonOp::
levelMulticolorGS(LevelData<EBCellFAB>&       a_phi,
                  const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBPoissonOp::levelMulticolorGS");

  // Loop over all possibilities (in all dimensions)
  for (int icolor = 0; icolor < m_colors.size(); icolor++)
    {
      colorGS(a_phi, a_rhs, m_colors[icolor], icolor);
    }
}

void EBPoissonOp::
colorGS(LevelData<EBCellFAB>&       a_phi,
        const LevelData<EBCellFAB>& a_rhs,
        const IntVect& a_color,
        int icolor)
{
  CH_TIME("EBPoissonOp::colorGS");

  //this is a multigrid operator so only homogeneous CF BC and null coar level
  CH_assert(a_rhs.ghostVect()    == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect()    == m_ghostCellsPhi);

  if (!m_exchangeCopier.isDefined())
    {
      DisjointBoxLayout grids = a_phi.disjointBoxLayout();

      m_exchangeCopier.define(grids, grids, a_phi.ghostVect(),  true);
    }
  a_phi.exchange(m_exchangeCopier);


  Real weight = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0 * m_beta * m_invDx2[idir];
    }
  weight = 1.0 / weight;


  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& phifab = a_phi[dit()];
      m_colorEBStencil[icolor][dit()]->cache(phifab);
    }

  GSColorAllRegular(  a_phi, a_rhs, a_color, weight, true);

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& phifab = a_phi[dit()];
      m_colorEBStencil[icolor][dit()]->uncache(phifab);
    }

  GSColorAllIrregular(a_phi, a_rhs, a_color, true, icolor);
}

void EBPoissonOp::
GSColorAllRegular(LevelData<EBCellFAB>&        a_phi,
                  const LevelData<EBCellFAB>&  a_rhs,
                  const IntVect&               a_color,
                  const Real&                  a_weight,
                  const bool&                  a_homogeneousPhysBC)
{
  CH_TIME("EBPoissonOp::GSColorAllRegular");
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);

  int nComps = a_phi.nComp();
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox(m_eblg.getDBL().get(dit()));
      BaseFab<Real>& phiFAB       = (a_phi[dit()] ).getSingleValuedFAB();
      const BaseFab<Real>& rhsFAB = (a_rhs[dit()] ).getSingleValuedFAB();

      Box loBox[SpaceDim],hiBox[SpaceDim];
      int hasLo[SpaceDim],hasHi[SpaceDim];

      applyDomainFlux(loBox, hiBox, hasLo, hasHi,
                      dblBox, nComps, phiFAB,
                      a_homogeneousPhysBC, dit(),m_beta);

      IntVect loIV = dblBox.smallEnd();
      IntVect hiIV = dblBox.bigEnd();

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (loIV[idir] % 2 != a_color[idir])
            {
              loIV[idir]++;
            }
        }

      if (loIV <= hiIV)
        {
          Box coloredBox(loIV, hiIV);

          for (int comp=0; comp<a_phi.nComp(); comp++)
            {
              FORT_DOALLREGULARMULTICOLOR(CHF_FRA1(phiFAB,comp),
                                          CHF_CONST_FRA1(rhsFAB,comp),
                                          CHF_CONST_REAL(a_weight),
                                          CHF_CONST_REAL(m_alpha),
                                          CHF_CONST_REAL(m_beta),
                                          CHF_CONST_REALVECT(m_dx),
                                          CHF_BOX(coloredBox));
            }
        }
    }
}

void EBPoissonOp::
GSColorAllIrregular(LevelData<EBCellFAB>&        a_phi,
                    const LevelData<EBCellFAB>&  a_rhs,
                    const IntVect&               a_color,
                    const bool&                  a_homogeneousPhysBC,
                    int icolor)
{
  CH_TIME("EBPoissonOp::GSColorAllIrregular");


  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& phifab = a_phi[dit()];
      const EBCellFAB& rhsfab = a_rhs[dit()];

      //phi = (I-lambda*L)phiOld
      m_colorEBStencil[icolor][dit()]->apply(phifab, phifab, false);

      //apply EB bcs to (I-lambda*L)phi (already done in colorStencil += fluxStencil, and hom only here))
      int comp = 0;
      const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[dit()];
      const BaseIVFAB<Real>& curBetaWeight  = m_betaDiagWeight[dit()];
      //apply domain bcs to (I-lambda*L)phi
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (m_vofItIrregColorDomLo[icolor][idir][dit()].reset(); m_vofItIrregColorDomLo[icolor][idir][dit()].ok();  ++m_vofItIrregColorDomLo[icolor][idir][dit()])
            {
              Real flux;
              const VolIndex& vof = m_vofItIrregColorDomLo[icolor][idir][dit()]();
              Real weightIrreg = m_alpha*curAlphaWeight(vof,0) + m_beta*curBetaWeight(vof,0);
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi[dit()],
                                      m_origin,m_dx,idir,Side::Lo, dit(), m_time,
                                      a_homogeneousPhysBC);

              phifab(vof,comp) -= -(1./weightIrreg) * flux * m_beta*m_invDx[idir];
            }
          for (m_vofItIrregColorDomHi[icolor][idir][dit()].reset(); m_vofItIrregColorDomHi[icolor][idir][dit()].ok();  ++m_vofItIrregColorDomHi[icolor][idir][dit()])
            {
              Real flux;
              const VolIndex& vof = m_vofItIrregColorDomHi[icolor][idir][dit()]();
              Real weightIrreg = m_alpha*curAlphaWeight(vof,0) + m_beta*curBetaWeight(vof,0);
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi[dit()],
                                      m_origin,m_dx,idir,Side::Hi,dit(), m_time,
                                      a_homogeneousPhysBC);

              phifab(vof,comp) += -(1./weightIrreg) * flux * m_beta*m_invDx[idir];
            }
        }

      //phi += lambda*rhsfab
      m_rhsColorEBStencil[icolor][dit()]->apply(phifab, rhsfab, true);
    }
}


//////
void EBPoissonOp::
getJacobiRelaxCoeff(LevelData<EBCellFAB>& a_relaxCoeff)
{

  // Diagonal weight for Jacobi on regular grid
  Real weightBase = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weightBase += -2.0 * m_beta * m_invDx2[idir];
    }
  weightBase = 1.0 / weightBase;       // Weight for plain (non-weighted) Jacobi
  Real weightRegular = EBPOISSONOP_JACOBI_OMEGA * weightBase; // ... for weighted Jacobi

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      // Most of the grid is regular, so set the regular weighting by default
      a_relaxCoeff[dit()].setVal(weightRegular);

      // Walk the EB cells, and correct their weighting
      const BaseIVFAB<Real>&      curAlphaWeight = m_alphaDiagWeight[dit()];
      const BaseIVFAB<Real>&      curBetaWeight  = m_betaDiagWeight[dit()];

      VoFIterator& ebvofit = m_vofItIrreg[dit()];
      for (ebvofit.reset(); ebvofit.ok(); ++ebvofit)
        {
          const VolIndex& vof = ebvofit();
          Real weightIrreg = EBPOISSONOP_JACOBI_OMEGA / (m_alpha*curAlphaWeight(vof, 0) + m_beta*curBetaWeight(vof, 0));
          a_relaxCoeff[dit()](vof, 0) = weightIrreg;
        }
    }
}


//////
void EBPoissonOp::
levelJacobi(LevelData<EBCellFAB>&       a_phi,
            const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBPoissonOp::levelJacobi");

  // Overhauled from previous in-place Gauss-Seidel-like method
  // This implementation is very inefficient in terms of memory usage,
  // and could be greatly improved.  It has the advantage, though,
  // of being very simple and easy to verify as correct.

  EBCellFactory factory(m_eblg.getEBISL());
  // Note: this is hardcoded to a single variable (component),
  // like some other code in this class
  LevelData<EBCellFAB> resid(m_eblg.getDBL(), 1, m_ghostCellsRHS, factory);

  residual(resid, a_phi, a_rhs, true);

  LevelData<EBCellFAB> relaxationCoeff(m_eblg.getDBL(), 1, m_ghostCellsRHS, factory);
  getJacobiRelaxCoeff(relaxationCoeff);

  EBLevelDataOps::scale(resid, relaxationCoeff);

  EBLevelDataOps::incr(a_phi, resid, 1.0);
}



/***/
void
EBPoissonOp::
applyOpRegularAllDirs(Box * a_loBox,
                      Box * a_hiBox,
                      int * a_hasLo,
                      int * a_hasHi,
                      Box & a_curDblBox,
                      Box & a_curPhiBox,
                      int a_nComps,
                      BaseFab<Real> & a_curOpPhiFAB,
                      const BaseFab<Real> & a_curPhiFAB,
                      bool a_homogeneousPhysBC,
                      const DataIndex& a_dit,
                      const Real& a_beta)
{
  CH_TIME("EBPoissonOp::applyOpRegularAllDirs");
  CH_assert(m_domainBC != NULL);

  //need to monkey with the ghost cells to account for boundary conditions
  BaseFab<Real>& phiFAB = (BaseFab<Real>&) a_curPhiFAB;
  applyDomainFlux(a_loBox, a_hiBox, a_hasLo, a_hasHi,
                  a_curDblBox, a_nComps, phiFAB,
                  a_homogeneousPhysBC, a_dit,m_beta);

  for (int comp = 0; comp<a_nComps; comp++)
    {

      FORT_REGGET1DLAPLACIAN_INPLACE(CHF_FRA1(a_curOpPhiFAB,comp),
                                     CHF_CONST_FRA1(a_curPhiFAB,comp),
                                     CHF_CONST_REAL(a_beta),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_BOX(a_curDblBox));
    }
}

/***/
void
EBPoissonOp::
applyDomainFlux(Box * a_loBox,
                Box * a_hiBox,
                int * a_hasLo,
                int * a_hasHi,
                Box & a_dblBox,
                int a_nComps,
                BaseFab<Real> & a_phiFAB,
                bool a_homogeneousPhysBC,
                const DataIndex& a_dit,
                const Real& a_beta)
{
  CH_TIME("EBPoissonOp::applyDomainFlux");
  CH_assert(m_domainBC != NULL);

  for (int idir=0; idir<SpaceDim; idir++)
    {

      EBArith::loHi(a_loBox[idir], a_hasLo[idir],
                             a_hiBox[idir], a_hasHi[idir],
                             m_eblg.getDomain(), a_dblBox, idir);

      for (int comp = 0; comp<a_nComps; comp++)
        {

          if (a_hasLo[idir] == 1)
            {
              Box lbox=a_loBox[idir];
              lbox.shift(idir,-1);
              BaseFab<Real> loFaceFlux(a_loBox[idir],a_nComps);
              int side = -1;

              m_domainBC->getFaceFlux(loFaceFlux,a_phiFAB,m_origin,m_dx,idir,Side::Lo,a_dit,m_time,a_homogeneousPhysBC);

              FORT_REGAPPLYDOMAINFLUX_INPLACE(CHF_FRA1(a_phiFAB,comp),
                                              CHF_CONST_FRA1(loFaceFlux,comp),
                                              CHF_CONST_REAL(m_dx[idir]),
                                              CHF_CONST_INT(side),
                                              CHF_CONST_INT(idir),
                                              CHF_BOX(lbox));
            }

          if (a_hasHi[idir] == 1)
            {
              Box hbox=a_hiBox[idir];
              hbox.shift(idir,1);
              BaseFab<Real> hiFaceFlux(a_hiBox[idir],a_nComps);
              int side = 1;

              m_domainBC->getFaceFlux(hiFaceFlux,a_phiFAB,m_origin,m_dx,idir,Side::Hi,a_dit,m_time,a_homogeneousPhysBC);

              FORT_REGAPPLYDOMAINFLUX_INPLACE(CHF_FRA1(a_phiFAB,comp),
                                              CHF_CONST_FRA1(hiFaceFlux,comp),
                                              CHF_CONST_REAL(m_dx[idir]),
                                              CHF_CONST_INT(side),
                                              CHF_CONST_INT(idir),
                                              CHF_BOX(hbox));
            }

        }
    }
}
#include "NamespaceFooter.H"
