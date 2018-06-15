#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoadBalance.H"
#include "EBArith.H"
#include "EBAMRPoissonOp.H"

#include "slowEBCO.H"
#include "EBQuadCFInterp.H"

#include "EBConductivityOpF_F.H"
#include "InterpF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "BCFunc.H"
#include "CH_Timer.H"
#include "BCFunc.H"
#include "EBLevelGrid.H"
#include "EBAlias.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"
bool slowEBCO::s_turnOffBCs = false; //REALLY needs to default to false
Real slowEBCO::getSafety()
{
  Real safety = 0.9;
  return safety;
}
void
slowEBCO::
getDivFStencil(VoFStencil&      a_vofStencil,
               const VolIndex&  a_vof,
               const DataIndex& a_dit)
{
  CH_TIME("slowEBCO::getDivFStencil");
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_dit];
  a_vofStencil.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());
          for (int iface = 0; iface < faces.size(); iface++)
            {
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, faces[iface], a_dit);
              Real areaFrac = ebisBox.areaFrac(faces[iface]);
              fluxStencil *= Real(isign)*areaFrac/m_dx;
              a_vofStencil += fluxStencil;
            }
        }
    }
}
/***/
/// stencil for flux computation.   the truly ugly part of this computation
/// beta and eta are multiplied in here
/****/
void
slowEBCO::
getFluxStencil(VoFStencil&      a_fluxStencil,
               const FaceIndex& a_face,
               const DataIndex& a_dit)
{
  CH_TIME("slowEBCO::getFluxStencil");
  //need to do this by interpolating to centroids
  //so get the stencil at each face center and add with
  //interpolation weights
  FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                     (*m_eblg.getCFIVS())[a_dit],
                                                     m_eblg.getEBISL()[a_dit],
                                                     m_eblg.getDomain());

  a_fluxStencil.clear();
  for (int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      getFaceCenteredFluxStencil(faceCentSten, face, a_dit);
      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
    }
}
/***/
void
slowEBCO::
getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                           const FaceIndex& a_face,
                           const DataIndex& a_dit)
{
  CH_TIME("slowEBCO::getFaceCenteredFluxStencil");
  //face centered gradient is just a centered diff
  int faceDir= a_face.direction();
  a_fluxStencil.clear();

  if (!a_face.isBoundary())
    {
      a_fluxStencil.add(a_face.getVoF(Side::Hi),  1.0/m_dx, 0);
      a_fluxStencil.add(a_face.getVoF(Side::Lo), -1.0/m_dx, 0);
      a_fluxStencil *= (*m_bcoef)[a_dit][faceDir](a_face,0);
    }
  else
    {
      //the boundary condition handles this one.
    }
}
/***/
void
slowEBCO::
setAlphaAndBeta(const Real& a_alpha,
                const Real& a_beta)
{
  CH_TIME("slowEBCO::setAlphaAndBeta");
  m_alpha = a_alpha;
  m_beta  = a_beta;
  defineStencils();
}
/***/
void
slowEBCO::
diagonalScale(LevelData<EBCellFAB> & a_rhs,
              bool a_kappaWeighted)
{
  CH_TIME("slowEBCO::kappaWeight");
  if (a_kappaWeighted)
    EBLevelDataOps::kappaWeight(a_rhs);

  //also have to weight by the coefficient
  for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
    {
      a_rhs[dit()] *= (*m_acoef)[dit()];
    }
}

void
slowEBCO::
divideByIdentityCoef(LevelData<EBCellFAB> & a_rhs)
{
  CH_TIME("slowEBCO::kappaWeight");

  //also have to weight by the coefficient
  for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
    {
      a_rhs[dit()] /= (*m_acoef)[dit()];
    }
}
/***/
slowEBCO::
slowEBCO(const EBLevelGrid &                                  a_eblgFine,
                 const EBLevelGrid &                                  a_eblg,
                 const EBLevelGrid &                                  a_eblgCoar,
                 const EBLevelGrid &                                  a_eblgCoarMG,
                 const RefCountedPtr<EBQuadCFInterp>&                 a_quadCFI,
                 const RefCountedPtr<ConductivityBaseDomainBC>&       a_domainBC,
                 const RefCountedPtr<ConductivityBaseEBBC>&           a_ebBC,
                 const Real    &                                      a_dx,
                 const Real    &                                      a_dxCoar,
                 const int&                                           a_refToFine,
                 const int&                                           a_refToCoar,
                 const bool&                                          a_hasFine,
                 const bool&                                          a_hasCoar,
                 const bool&                                          a_hasMGObjects,
                 const bool&                                          a_layoutChanged,
                 const Real&                                          a_alpha,
                 const Real&                                          a_beta,
                 const RefCountedPtr<LevelData<EBCellFAB> >&          a_acoef,
                 const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
                 const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
                 const IntVect&                                       a_ghostCellsPhi,
                 const IntVect&                                       a_ghostCellsRHS,
                 const int&                                           a_relaxType)
  : m_ghostCellsPhi( a_ghostCellsPhi ),
    m_ghostCellsRHS( a_ghostCellsRHS )
{

  CH_TIME("slowEBCO::ConductivityOp");
  int ncomp = 1;
  m_relaxType = a_relaxType;
  m_acoef    = a_acoef;
  m_bcoef    = a_bcoef;
  m_bcoIrreg = a_bcoIrreg;

  m_quadCFIWithCoar = a_quadCFI;
  m_hasFine        = a_hasFine;
  m_hasCoar        = a_hasCoar;
  m_eblg            = a_eblg;
  m_domainBC       = a_domainBC;
  m_ebBC           = a_ebBC;
  m_dx             = a_dx;
  m_alpha          = a_alpha;
  m_beta           = a_beta;
  m_hasInterpAve = false;
  m_dxCoar         = a_dxCoar;


  if (m_hasFine)
    {
      m_eblgFine       = a_eblgFine;
      m_dxFine         = m_dx/a_refToFine;
      m_refToFine      = a_refToFine;
    }
  else
    {
      m_refToFine    = 1;
    }

  EBCellFactory fact(m_eblg.getEBISL());
  m_relCoef.define(m_eblg.getDBL(), 1, IntVect::Zero, fact);
   if (m_hasCoar)
    {
      m_eblgCoar       = a_eblgCoar;
      m_refToCoar      = a_refToCoar   ;


      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_loCFIVS[idir].define(m_eblg.getDBL());
          m_hiCFIVS[idir].define(m_eblg.getDBL());

          for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
            {
              m_loCFIVS[idir][dit()].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit()),
                                            m_eblg.getDBL(), idir,Side::Lo);
              m_hiCFIVS[idir][dit()].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit()),
                                            m_eblg.getDBL(), idir,Side::Hi);
            }
        }

      //if this fails, then the AMR grids violate proper nesting.
      ProblemDomain domainCoarsenedFine;
      DisjointBoxLayout dblCoarsenedFine;

      int maxBoxSize = 32;
      bool dumbool;
      bool hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarsenedFine,
                                                          domainCoarsenedFine,
                                                          m_eblg.getDBL(),
                                                          m_eblg.getEBISL(),
                                                          m_eblg.getDomain(),
                                                          m_refToCoar,
                                                          m_eblg.getEBIS(),
                                                          maxBoxSize, dumbool);

      //should follow from coarsenable
      if (hasCoarser)
        {
          m_eblgCoarsenedFine = EBLevelGrid(dblCoarsenedFine, domainCoarsenedFine, 4, Chombo_EBIS::instance());
          m_hasInterpAve = true;
          m_ebInterp.define( m_eblg.getDBL(),     m_eblgCoar.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoar.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             a_ghostCellsPhi);
          m_ebAverage.define(m_eblg.getDBL(),     m_eblgCoarsenedFine.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoarsenedFine.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             a_ghostCellsRHS);

        }
      else
        {
          m_hasInterpAve = false;
        }
    }
  else
    {
      m_refToCoar    = 1;
    }

  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects
  m_hasMGObjects = a_hasMGObjects;
  m_layoutChanged = a_layoutChanged;
  if (m_hasMGObjects)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain(), mgRef, ncomp, m_eblg.getEBIS(),
                           a_ghostCellsPhi);
      m_ebAverageMG.define(m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain() , mgRef, ncomp, m_eblg.getEBIS(),
                           a_ghostCellsRHS);

    }

  //define stencils for the operator
  defineStencils();

}
/***/
void
slowEBCO::
calculateRelaxationCoefficient()
{
  //define regular relaxation coefficent
  Real safety = getSafety();
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = m_eblg.getDBL().get(dit());
      const EBCellFAB& acofab = (*m_acoef)[dit()];
      //initialize lambda = alpha*acoef
      m_relCoef[dit()].setVal(0.);
      m_relCoef[dit()].plus(acofab, 0, 0, 1);
      m_relCoef[dit()]*= m_alpha;

      BaseFab<Real>& regRel =   m_relCoef[dit()].getSingleValuedFAB();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          BaseFab<Real>& regBCo = (*m_bcoef)[dit()][idir].getSingleValuedFAB();
          FORT_DECRINVRELCOEFEBCO(CHF_FRA1(regRel,0),
                                  CHF_FRA1(regBCo,0),
                                  CHF_CONST_REAL(m_beta),
                                  CHF_BOX(grid),
                                  CHF_REAL(m_dx),
                                  CHF_INT(idir));
        }

      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)
      FORT_INVERTLAMBDAEBCO(CHF_FRA1(regRel,0),
                            CHF_REAL(safety),
                            CHF_BOX(grid));

      VoFIterator& vofit = m_vofIterIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();

          //add in identity term
          Real alphaWeight = m_alphaDiagWeight[dit()](VoF, 0);
          Real  betaWeight =  m_betaDiagWeight[dit()](VoF, 0);
          alphaWeight *= m_alpha;
          betaWeight  *= m_beta;

          Real diagWeight  = alphaWeight + betaWeight;
          //set irregular relaxation coef
          m_relCoef[dit()](VoF, 0) = safety/diagWeight;
        }
    }
}
/******/
void
slowEBCO::
defineStencils()
{
  CH_TIME("slowEBCO::defineStencils");
  // create ebstencil for irregular applyOp
  m_opEBStencil.define(m_eblg.getDBL());
  // create vofstencils for applyOp and
  LayoutData<BaseIVFAB<VoFStencil> >  opStencil;
  opStencil.define(m_eblg.getDBL());

  Real fakeBeta = 1;
  m_domainBC->setCoef(m_eblg,   fakeBeta ,      m_bcoef   );
  m_ebBC->setCoef(    m_eblg,   fakeBeta ,      m_bcoIrreg);

  Real dxScale = 1.0/m_dx;
  m_ebBC->define((*m_eblg.getCFIVS()), dxScale); //has to happen AFTER coefs are set
  LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);
  //create and define colored stencils (2 parts)
  //define regular relaxation coefficent
  Real safety = getSafety();

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = m_eblg.getDBL().get(dit());
      const EBCellFAB& acofab = (*m_acoef)[dit()];
      //initialize lambda = alpha*acoef
      m_relCoef[dit()].setVal(0.);
      m_relCoef[dit()].plus(acofab, 0, 0, 1);
      m_relCoef[dit()]*= m_alpha;

      BaseFab<Real>& regRel =   m_relCoef[dit()].getSingleValuedFAB();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          BaseFab<Real>& regBCo = (*m_bcoef)[dit()][idir].getSingleValuedFAB();
          FORT_DECRINVRELCOEFEBCO(CHF_FRA1(regRel,0),
                                  CHF_FRA1(regBCo,0),
                                  CHF_CONST_REAL(m_beta),
                                  CHF_BOX(grid),
                                  CHF_REAL(m_dx),
                                  CHF_INT(idir));
        }

      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)
      FORT_INVERTLAMBDAEBCO(CHF_FRA1(regRel,0),
                            CHF_REAL(safety),
                            CHF_BOX(grid));
    }

  m_vofIterIrreg.define(     m_eblg.getDBL()); // vofiterator cache
  m_vofIterMulti.define(     m_eblg.getDBL()); // vofiterator cache
  m_alphaDiagWeight.define(  m_eblg.getDBL());
  m_betaDiagWeight.define(   m_eblg.getDBL());
  Box sideBoxLo[SpaceDim];
  Box sideBoxHi[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box domainBox = m_eblg.getDomain().domainBox();
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_vofIterDomLo[idir].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofIterDomHi[idir].define( m_eblg.getDBL()); // vofiterator cache for domain hi
    }
  EBArith::getMultiColors(m_colors);

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& curBox = m_eblg.getDBL().get(dit());
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const EBGraph& ebgraph = ebisBox.getEBGraph();

      IntVectSet irregIVS = ebisBox.getIrregIVS(curBox);
      IntVectSet multiIVS = ebisBox.getMultiCells(curBox);

      opStencil[dit()].define(irregIVS,ebgraph, 1);

      //cache the vofIterators
      m_alphaDiagWeight[dit()].define(irregIVS,ebisBox.getEBGraph(), 1);
      m_betaDiagWeight [dit()].define(irregIVS,ebisBox.getEBGraph(), 1);
      m_vofIterIrreg   [dit()].define(irregIVS,ebisBox.getEBGraph());
      m_vofIterMulti   [dit()].define(multiIVS,ebisBox.getEBGraph());

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet loIrreg = irregIVS;
          IntVectSet hiIrreg = irregIVS;
          loIrreg &= sideBoxLo[idir];
          hiIrreg &= sideBoxHi[idir];
          m_vofIterDomLo[idir][dit()].define(loIrreg,ebisBox.getEBGraph());
          m_vofIterDomHi[idir][dit()].define(hiIrreg,ebisBox.getEBGraph());
        }

      VoFIterator& vofit = m_vofIterIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();

          VoFStencil& curStencil = opStencil[dit()](VoF,0);

          //bcoef is included here in the flux consistent
          //with the regular
          getDivFStencil(curStencil,VoF, dit());
          if (fluxStencil != NULL)
            {
              BaseIVFAB<VoFStencil>& fluxStencilBaseIVFAB = (*fluxStencil)[dit()];
              //this fills the stencil with the gradient*beta*bcoef
              VoFStencil  fluxStencilPt = fluxStencilBaseIVFAB(VoF,0);
              curStencil += fluxStencilPt;
            }
          Real betaWeight = EBArith::getDiagWeight(curStencil, VoF);
          const IntVect& iv = VoF.gridIndex();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Box loSide = bdryLo(m_eblg.getDomain(),idir);
              loSide.shiftHalf(idir,1);
              Real adjust = 0;
              if (loSide.contains(iv))
                {
                  Real faceAreaFrac = 0.0;
                  Real weightedAreaFrac = 0;
                  Vector<FaceIndex> faces = ebisBox.getFaces(VoF,idir,Side::Lo);
                  for (int i = 0; i < faces.size(); i++)
                    {
                      weightedAreaFrac = ebisBox.areaFrac(faces[i]);
                      weightedAreaFrac *= (*m_bcoef)[dit()][idir](faces[i],0);
                      faceAreaFrac +=  weightedAreaFrac;
                    }
                  adjust += -weightedAreaFrac /(m_dx*m_dx);
                }
              Box hiSide = bdryHi(m_eblg.getDomain(),idir);
              hiSide.shiftHalf(idir,-1);
              if (hiSide.contains(iv))
                {
                  Real faceAreaFrac = 0.0;
                  Real weightedAreaFrac = 0;
                  Vector<FaceIndex> faces = ebisBox.getFaces(VoF,idir,Side::Hi);
                  for (int i = 0; i < faces.size(); i++)
                    {
                      weightedAreaFrac = ebisBox.areaFrac(faces[i]);
                      weightedAreaFrac *= (*m_bcoef)[dit()][idir](faces[i],0);
                      faceAreaFrac +=  weightedAreaFrac;
                    }
                  adjust += -weightedAreaFrac /(m_dx*m_dx);
                }
              betaWeight += adjust;
            }

          //add in identity term
          Real volFrac = ebisBox.volFrac(VoF);
          Real alphaWeight = (*m_acoef)[dit()](VoF, 0);
          alphaWeight *= volFrac;

          m_alphaDiagWeight[dit()](VoF, 0) = alphaWeight;
          m_betaDiagWeight[dit()](VoF, 0)  = betaWeight;
        }

      //Operator ebstencil
      m_opEBStencil[dit()] = RefCountedPtr<EBStencil>
        (new EBStencil(m_vofIterIrreg[dit()].getVector(), opStencil[dit()], m_eblg.getDBL().get(dit()),
                       m_eblg.getEBISL()[dit()], m_ghostCellsPhi, m_ghostCellsRHS, 0, true));
    }//dit

  calculateRelaxationCoefficient();
}
/***/
void
slowEBCO::
residual(LevelData<EBCellFAB>&       a_residual,
         const LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         bool                        a_homogeneousPhysBC)
{
  CH_TIME("slowEBCO::residual");
  //this is a multigrid operator so only homogeneous CF BC
  //and null coar level
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  applyOp(a_residual,a_phi,NULL, a_homogeneousPhysBC, true);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}
/***/
void
slowEBCO::
preCond(LevelData<EBCellFAB>&       a_lhs,
        const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("slowEBCO::preCond");
  EBLevelDataOps::assign(a_lhs, a_rhs);
  EBLevelDataOps::scale(a_lhs, m_relCoef);

  relax(a_lhs, a_rhs, 40);
}
/***/
void
slowEBCO::
applyOp(LevelData<EBCellFAB>&             a_opPhi,
        const LevelData<EBCellFAB>&       a_phi,
        bool                              a_homogeneousPhysBC)
{
  //homogeneous CFBCs because that is all we can do.
  applyOp(a_opPhi, a_phi, NULL, a_homogeneousPhysBC, true);
}
/**/
void
slowEBCO::
applyOp(LevelData<EBCellFAB>&                    a_lhs,
        const LevelData<EBCellFAB>&              a_phi,
        const LevelData<EBCellFAB>* const        a_phiCoar,
        const bool&                              a_homogeneousPhysBC,
        const bool&                              a_homogeneousCFBC)
{
  LevelData<EBCellFAB>& phi = const_cast<LevelData<EBCellFAB>&>(a_phi);
  if (m_hasCoar && (!s_turnOffBCs))
    {
      applyCFBCs(phi, a_phiCoar, a_homogeneousCFBC);
    }
  phi.exchange(phi.interval());

  EBLevelDataOps::setToZero(a_lhs);
  incr( a_lhs, a_phi, m_alpha);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_lhs[dit()].mult((*m_acoef)[dit()], 0, 0, 1);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          incrOpRegularDir(a_lhs[dit()], a_phi[dit()], a_homogeneousPhysBC, idir, dit());
        }
      applyOpIrregular(a_lhs[dit()], a_phi[dit()], a_homogeneousPhysBC, dit());
    }
}
/***/
void
slowEBCO::
incrOpRegularDir(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const int&             a_dir,
                 const DataIndex&       a_datInd)
{
  const Box& grid = m_eblg.getDBL()[a_datInd];
  Box domainFaces = m_eblg.getDomain().domainBox();
  domainFaces.surroundingNodes(a_dir);
  Box interiorFaces = grid;
  interiorFaces.surroundingNodes(a_dir);
  interiorFaces.grow(a_dir, 1);
  interiorFaces &=  domainFaces;
  interiorFaces.grow( a_dir, -1);

  //do flux difference for interior points
  FArrayBox interiorFlux(interiorFaces, 1);
  const FArrayBox& phi  = (FArrayBox&)(a_phi.getSingleValuedFAB());
  getFlux(interiorFlux, phi,  interiorFaces, a_dir, m_dx, a_datInd);

  Box loBox, hiBox, centerBox;
  int hasLo, hasHi;
  EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_eblg.getDomain(),grid, a_dir);

  //do the low high center thing
  BaseFab<Real>& reglhs        = a_lhs.getSingleValuedFAB();
  Box dummyBox(IntVect::Zero, IntVect::Unit);
  FArrayBox domainFluxLo(dummyBox, 1);
  FArrayBox domainFluxHi(dummyBox, 1);

  RealVect origin = RealVect::Zero;
  Real time = 0.0;
  RealVect dxVect = m_dx*RealVect::Unit;
  if (hasLo==1)
    {
      Box loBoxFace = loBox;
      loBoxFace.shiftHalf(a_dir, -1);
      domainFluxLo.resize(loBoxFace, 1);
      if (!s_turnOffBCs)
        {
          //using EBAMRPoissonOp BCs which require a cell centered box.
          domainFluxLo.shiftHalf(a_dir, 1);
          m_domainBC->getFaceFlux(domainFluxLo,phi,origin,dxVect,a_dir,Side::Lo,a_datInd,time,a_homogeneous);
          domainFluxLo *= m_beta;
          domainFluxLo.shiftHalf(a_dir,-1);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          for (BoxIterator boxit(loBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]++;
              domainFluxLo(iv, 0) = interiorFlux(ivn, 0);
            }
        }
    }
  if (hasHi==1)
    {
      Box hiBoxFace = hiBox;
      hiBoxFace.shiftHalf(a_dir, 1);
      domainFluxHi.resize(hiBoxFace, 1);
      if (!s_turnOffBCs)
        {
          //using EBAMRPoissonOp BCs which require a cell centered box.
          domainFluxHi.shiftHalf(a_dir, -1);
          m_domainBC->getFaceFlux(domainFluxHi,phi,origin,dxVect,a_dir,Side::Hi,a_datInd,time,a_homogeneous);
          domainFluxHi *= m_beta;
          domainFluxHi.shiftHalf(a_dir,  1);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          for (BoxIterator boxit(hiBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]--;
              domainFluxHi(iv, 0) = interiorFlux(ivn, 0);
            }
        }
    }
  Real unity = 1.0;  //beta already in the flux
  FORT_INCRAPPLYEBCO(CHF_FRA1(reglhs,0),
                     CHF_CONST_FRA1(interiorFlux, 0),
                     CHF_CONST_FRA1(domainFluxLo, 0),
                     CHF_CONST_FRA1(domainFluxHi, 0),
                     CHF_CONST_REAL(unity),
                     CHF_CONST_REAL(m_dx),
                     CHF_BOX(loBox),
                     CHF_BOX(hiBox),
                     CHF_BOX(centerBox),
                     CHF_CONST_INT(hasLo),
                     CHF_CONST_INT(hasHi),
                     CHF_CONST_INT(a_dir));
}
/*****/
void
slowEBCO::
applyOpIrregular(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const DataIndex&       a_datInd)
{
  RealVect vectDx = m_dx*RealVect::Unit;
  m_opEBStencil[a_datInd]->apply(a_lhs, a_phi, m_alphaDiagWeight[a_datInd], m_alpha, m_beta, false);

  if (!a_homogeneous)
    {
      const Real factor = m_beta/m_dx; //beta and bcoef handled within applyEBFlux
      m_ebBC->applyEBFlux(a_lhs, a_phi, m_vofIterIrreg[a_datInd], (*m_eblg.getCFIVS()),
                          a_datInd, RealVect::Zero, vectDx, factor,
                          a_homogeneous, 0.0);
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int comp = 0;
      for (m_vofIterDomLo[idir][a_datInd].reset(); m_vofIterDomLo[idir][a_datInd].ok();  ++m_vofIterDomLo[idir][a_datInd])
        {
          Real flux;
          const VolIndex& vof = m_vofIterDomLo[idir][a_datInd]();
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero,vectDx,idir,Side::Lo, a_datInd, 0.0,
                                  a_homogeneous);
          //area gets multiplied in by bc operator
          a_lhs(vof,comp) -= flux*m_beta/m_dx;
        }
      for (m_vofIterDomHi[idir][a_datInd].reset(); m_vofIterDomHi[idir][a_datInd].ok();  ++m_vofIterDomHi[idir][a_datInd])
        {
          Real flux;
          const VolIndex& vof = m_vofIterDomHi[idir][a_datInd]();
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero,vectDx,idir,Side::Hi,a_datInd,0.0,
                                  a_homogeneous);
          //area gets multiplied in by bc operator
          a_lhs(vof,comp) += flux*m_beta/m_dx;
        }
    }
}

void
slowEBCO::
applyOpNoBoundary(LevelData<EBCellFAB>&        a_opPhi,
                  const LevelData<EBCellFAB>&  a_phi)
{
  s_turnOffBCs = true;
  applyOp(a_opPhi, a_phi, true);
  s_turnOffBCs = false;
}
/***/
void
slowEBCO::
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
}
/***/
void
slowEBCO::
createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                const LevelData<EBCellFAB>& a_rhs,
                const int &                 a_refRat)
{
  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  CH_assert(m_eblg.getDBL().coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, m_eblg.getDBL(), a_refRat);

  EBISLayout ebislCoarsenedFine;
  IntVect ghostVec = a_rhs.ghostVect();
  //const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), a_refRat);
  m_eblg.getEBIS()->fillEBISLayout(ebislCoarsenedFine, dblCoarsenedFine, coarDom , ghostVec[0]);
  if (m_refToCoar > 1)
    {
      ebislCoarsenedFine.setMaxRefinementRatio(m_refToCoar, m_eblg.getEBIS());
    }

  //create coarsened data
  EBCellFactory ebcellfactCoarsenedFine(ebislCoarsenedFine);
  a_lhs.define(dblCoarsenedFine, ncomp,ghostVec, ebcellfactCoarsenedFine);
}
/***/
void
slowEBCO::
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  EBLevelDataOps::assign(a_lhs,a_rhs);
}
/***/
Real
slowEBCO::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
}
/***/
void
slowEBCO::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}
/***/
void
slowEBCO::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}
/***/
void
slowEBCO::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  EBLevelDataOps::scale(a_lhs,a_scale);
}
/***/
Real
slowEBCO::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
  //  return  EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);
  return EBArith::norm(a_rhs, a_rhs.disjointBoxLayout(), m_eblg.getEBISL(), 0, 1);
}
/***/
void
slowEBCO::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  EBLevelDataOps::setToZero(a_lhs);
}
/***/
void
slowEBCO::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  EBLevelDataOps::setVal(a_lhs, a_value);
}
/***/
void
slowEBCO::
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
/***/
void
slowEBCO::
relax(LevelData<EBCellFAB>&       a_phi,
      const LevelData<EBCellFAB>& a_rhs,
      int                         a_iterations)
{
  if (m_relaxType == 0)
    {
      relaxPoiJac(a_phi, a_rhs, a_iterations);
    }
  else if (m_relaxType == 1)
    {
      relaxGauSai(a_phi, a_rhs, a_iterations);
    }
  else
    {
      MayDay::Error("slowEBCO::bogus relaxtype");
    }
}

void
slowEBCO::
relaxGauSai(LevelData<EBCellFAB>&       a_phi,
            const LevelData<EBCellFAB>& a_rhs,
            int                         a_iterations)
{
  CH_TIME("slowEBCO::relax");

  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);

  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  LevelData<EBCellFAB> lphi;
  create(lphi, a_rhs);
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for (int icolor = 0; icolor < m_colors.size(); icolor++)
        {
          applyHomogeneousCFBCs(a_phi);

          //after this lphi = L(phi)
          //this call contains bcs and exchange
          applyOp(  lphi,  a_phi, true);
          gsrbColor(a_phi, lphi, a_rhs, m_colors[icolor]);
        }
    }
}

void
slowEBCO::
relaxPoiJac(LevelData<EBCellFAB>&       a_phi,
            const LevelData<EBCellFAB>& a_rhs,
            int                         a_iterations)
{
  CH_TIME("slowEBCO::relax");

  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);

  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  LevelData<EBCellFAB> lphi;
  create(lphi, a_rhs);
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      applyHomogeneousCFBCs(a_phi);

      //after this lphi = L(phi)
      //this call contains bcs and exchange
      applyOp(  lphi,  a_phi, true);

      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          lphi[dit()] -=     a_rhs[dit()];
          lphi[dit()] *= m_relCoef[dit()];
          //this is a safety factor because pt jacobi needs a smaller
          //relaxation param
          lphi[dit()] *= -0.5;
          a_phi[dit()] += lphi[dit()];
        }
    }
}
/****/
void
slowEBCO::
gsrbColor(LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>& a_lph,
          const LevelData<EBCellFAB>& a_rhs,
          const IntVect&              a_color)
{

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox  = dbl.get(dit());
      BaseFab<Real>&       regPhi =     a_phi[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regLph =     a_lph[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regRhs =     a_rhs[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regRel = m_relCoef[dit()].getSingleValuedFAB();
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
          FORT_GSRBEBCO(CHF_FRA1(regPhi,0),
                        CHF_CONST_FRA1(regLph,0),
                        CHF_CONST_FRA1(regRhs,0),
                        CHF_CONST_FRA1(regRel,0),
                        CHF_BOX(coloredBox));
        }

      for (m_vofIterMulti[dit()].reset(); m_vofIterMulti[dit()].ok(); ++m_vofIterMulti[dit()])
        {
          const VolIndex& vof = m_vofIterMulti[dit()]();
          const IntVect& iv = vof.gridIndex();

          bool doThisVoF = true;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (iv[idir] % 2 != a_color[idir])
                {
                  doThisVoF = false;
                  break;
                }
            }

          if (doThisVoF)
            {
              Real lph    = a_lph[dit()](vof, 0);
              Real rhs    = a_rhs[dit()](vof, 0);
              Real resid  = rhs - lph;
              Real lambda = m_relCoef[dit()](vof, 0);
              a_phi[dit()](vof, 0) += lambda*resid;
            }
        }
    }
}
/***/
void slowEBCO::
restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                 LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_rhsThisLevel)
{
  CH_TIME("slowEBCO::restrictResidual");

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
  if (m_layoutChanged)
    {
      m_ebAverageMG.average(a_resCoar, resThisLevel, variables);
    }
  else
    {
      m_ebAverageMG.averageMG(a_resCoar, resThisLevel, variables);
    }
}
/****/
void slowEBCO::
prolongIncrement(LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_correctCoar)
{
  CH_TIME("slowEBCO::prolongIncrement");
  Interval vars(0, 0);
  if (m_layoutChanged)
    {
      m_ebInterpMG.pwcInterp(a_phiThisLevel, a_correctCoar, vars);
    }
  else
    {
      m_ebInterpMG.pwcInterpMG(a_phiThisLevel, a_correctCoar, vars);
    }
}
/***/
void
slowEBCO::
applyCFBCs(LevelData<EBCellFAB>&             a_phi,
           const LevelData<EBCellFAB>* const a_phiCoar,
           bool a_homogeneousCFBC)
{
  CH_TIMERS("slowEBCO::applyCFBCs");
  CH_TIMER("inhomogeneous_cfbcs_define",t1);
  CH_TIMER("inhomogeneous_cfbcs_execute",t3);
  CH_TIMER("homogeneous_cfbs",t2);
  CH_assert(a_phi.nComp() == 1);

  if (m_hasCoar)
    {
      if (!a_homogeneousCFBC)
        {
          CH_START(t1);
          if (a_phiCoar==NULL)
            {
              MayDay::Error("cannot enforce inhomogeneous CFBCs with NULL coar");
            }
          //define coarse fine interpolation object on the fly
          //because most operators do not need it
          CH_assert(a_phiCoar->nComp() == 1);
          CH_STOP(t1);

          CH_START(t3);
          Interval interv(0,0);
          m_quadCFIWithCoar->interpolate(a_phi, *a_phiCoar, interv);

          CH_STOP(t3);

        }
      else
        {
          CH_START(t2);
          applyHomogeneousCFBCs(a_phi);
          CH_STOP(t2);
        }
    }
}
/***/
void
slowEBCO::
applyHomogeneousCFBCs(LevelData<EBCellFAB>&   a_phi)
{
  CH_TIME("slowEBCO::applyHomogeneousCFBCs");
  CH_assert(a_phi.nComp() == 1);
  CH_assert( a_phi.ghostVect() >= IntVect::Unit);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); sit.next())
            {
              applyHomogeneousCFBCs(a_phi[dit()],dit(),idir,sit());
            }
        }
    }
}
/***/
void
slowEBCO::
applyHomogeneousCFBCs(EBCellFAB&            a_phi,
                      const DataIndex&      a_datInd,
                      int                   a_idir,
                      Side::LoHiSide        a_hiorlo)
{
  if (m_hasCoar)
    {
      CH_TIMERS("slowEBCO::applyHomogeneousCFBCs2");
      CH_TIMER("packed_applyHomogeneousCFBCs",t1);
      CH_TIMER("unpacked_applyHomogeneousCFBCs",t2);
      CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
      CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
      CH_assert(a_phi.nComp() == 1);
      int ivar = 0;

      const CFIVS* cfivsPtr = NULL;

      if (a_hiorlo == Side::Lo)
        {
          cfivsPtr = &m_loCFIVS[a_idir][a_datInd];
        }
      else
        {
          cfivsPtr = &m_hiCFIVS[a_idir][a_datInd];
        }

      const IntVectSet& interpIVS = cfivsPtr->getFineIVS();
      if (cfivsPtr->isPacked() )
        {
          CH_START(t1);
          const int ihiorlo = sign(a_hiorlo);
          FORT_INTERPHOMO(CHF_FRA(a_phi.getSingleValuedFAB()),
                          CHF_BOX(cfivsPtr->packedBox()),
                          CHF_CONST_REAL(m_dx),
                          CHF_CONST_REAL(m_dxCoar),
                          CHF_CONST_INT(a_idir),
                          CHF_CONST_INT(ihiorlo));

          CH_STOP(t1);
        }
      else
        {
          if (!interpIVS.isEmpty())
            {
              CH_START(t2);
              Real halfdxcoar = m_dxCoar/2.0;
              Real halfdxfine = m_dx/2.0;
              Real xg = halfdxcoar -   halfdxfine;
              Real xc = halfdxcoar +   halfdxfine;
              Real xf = halfdxcoar + 3*halfdxfine;
              Real hf = m_dx;
              Real denom = xf*xc*hf;

              const EBISBox&  ebisBox = m_eblg.getEBISL()[a_datInd];
              const EBGraph&  ebgraph = m_eblg.getEBISL()[a_datInd].getEBGraph();
              for (VoFIterator vofit(interpIVS, ebgraph); vofit.ok(); ++vofit)
                {
                  const VolIndex& VoFGhost = vofit();

                  IntVect ivGhost  = VoFGhost.gridIndex();
                  IntVect ivClose =  ivGhost;
                  IntVect ivFar   =  ivGhost;

                  Vector<VolIndex> farVoFs;
                  Vector<VolIndex> closeVoFs = ebisBox.getVoFs(VoFGhost,
                                                               a_idir,
                                                               flip(a_hiorlo),
                                                               1);
                  bool hasClose = (closeVoFs.size() > 0);
                  bool hasFar = false;
                  Real phic = 0.0;
                  Real phif = 0.0;
                  if (hasClose)
                    {
                      const int& numClose = closeVoFs.size();
                      for (int iVof=0;iVof<numClose;iVof++)
                        {
                          const VolIndex& vofClose = closeVoFs[iVof];
                          phic += a_phi(vofClose,0);
                        }
                      phic /= Real(numClose);

                      farVoFs = ebisBox.getVoFs(VoFGhost,
                                                a_idir,
                                                flip(a_hiorlo),
                                                2);
                      hasFar   = (farVoFs.size()   > 0);
                      if (hasFar)
                        {
                          const int& numFar = farVoFs.size();
                          for (int iVof=0;iVof<numFar;iVof++)
                            {
                              const VolIndex& vofFar = farVoFs[iVof];
                              phif += a_phi(vofFar,0);
                            }
                          phif /= Real(numFar);
                        }
                    }

                  Real phiGhost;
                  if (hasClose && hasFar)
                    {
                      // quadratic interpolation  phi = ax^2 + bx + c
                      Real A = (phif*xc - phic*xf)/denom;
                      Real B = (phic*hf*xf - phif*xc*xc + phic*xf*xc)/denom;

                      phiGhost = A*xg*xg + B*xg;
                    }
                  else if (hasClose)
                    {
                      //linear interpolation
                      Real slope =  phic/xc;
                      phiGhost   =  slope*xg;
                    }
                  else
                    {
                      phiGhost = 0.0; //nothing to interpolate from
                    }
                  a_phi(VoFGhost, ivar) = phiGhost;
                }
              CH_STOP(t2);
            }
        }
    }
}
/***/
int slowEBCO::
refToCoarser()
{
  return m_refToCoar;
}
/***/
int slowEBCO::
refToFiner()
{
  return m_refToFine;
}
/***/
void slowEBCO::
AMRResidual(LevelData<EBCellFAB>&       a_residual,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            const LevelData<EBCellFAB>& a_rhs,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("slowEBCO::AMRResidual");
  CH_TIMER("AMROperator", t1);
  CH_TIMER("axby", t2);
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_residual.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  CH_START(t1);
  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar,
              a_homogeneousPhysBC, a_finerOp);
  CH_STOP(t1);

  //multiply by -1 so a_residual now holds -L(phi)
  //add in rhs so a_residual = rhs - L(phi)
  CH_START(t2);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  CH_STOP(t2);
}
/***/
void slowEBCO::
AMROperator(LevelData<EBCellFAB>&       a_LofPhi,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("slowEBCO::AMROperator");
  CH_TIMER("applyOp", t1);
  CH_TIMER("reflux", t2);
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_LofPhi.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  //apply the operator between this and the next coarser level.
  CH_START(t1);
  applyOp(a_LofPhi, a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);
  CH_STOP(t1);

  //now reflux to enforce flux-matching from finer levels
  if (m_hasFine)
    {
      CH_assert(a_finerOp != NULL);
      CH_START(t2);

      reflux(a_LofPhi, a_phiFine, a_phi, a_finerOp);

      CH_STOP(t2);
    }
}
/***/
void
slowEBCO::
reflux(LevelData<EBCellFAB>& a_residual,
       const LevelData<EBCellFAB>& a_phiFine,
       const LevelData<EBCellFAB>& a_phi,
       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("slowEBCO::fastReflux");
  CH_TIMER("define",t1);
  CH_TIMER("setToZero",t2);
  CH_TIMER("incrementCoar",t3);
  CH_TIMER("incrementFine",t4);
  CH_TIMER("reflux_from_reg",t5);
  Interval interv(0,0);
  CH_START(t1);
  CH_assert(a_phiFine.nComp() == 1);
  int ncomp = 1;
  if (!m_fastFR.isDefined())
    {
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp);
    }
  CH_STOP(t1);
  CH_START(t2);
  m_fastFR.setToZero();
  CH_STOP(t2);
  CH_START(t3);
  incrementFRCoar(m_fastFR, a_phiFine, a_phi);
  CH_STOP(t3);

  CH_START(t4);
  incrementFRFine(m_fastFR, a_phiFine, a_phi, a_finerOp);
  CH_STOP(t4);
  CH_START(t5);

  Real scale = 1.0/m_dx;
  m_fastFR.reflux(a_residual, interv, scale);

  CH_STOP(t5);
}
/***/
void
slowEBCO::
incrementFRCoar(EBFastFR&             a_fluxReg,
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("slowEBCO::incrementFRCoar");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  int ncomp = 1;
  Interval interv(0,0);

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& coarfab = a_phi[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const Box&  box = m_eblg.getDBL().get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //no boundary faces here.

          Box ghostedBox = box;
          ghostedBox.grow(1);
          ghostedBox.grow(idir,-1);
          ghostedBox &= m_eblg.getDomain();

          EBFaceFAB coarflux(ebisBox, ghostedBox, idir, ncomp);

          getFlux(coarflux, coarfab, ghostedBox, box, m_eblg.getDomain(), ebisBox, m_dx, dit(), idir);
          Real scale = 1.0; //beta and bcoef already in flux
          for (SideIterator sit; sit.ok(); ++sit)
            {
              a_fluxReg.incrementCoarseBoth(coarflux, scale, dit(), interv, idir, sit());
            }
        }
    }
}
/***/
void
slowEBCO::
getFlux(EBFluxFAB&                    a_flux,
        const LevelData<EBCellFAB>&   a_data,
        const Box&                    a_grid,
        const DataIndex&              a_dit,
        Real                          a_scale)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box ghostedBox = a_grid;
      ghostedBox.grow(1);
      ghostedBox.grow(idir,-1);
      ghostedBox &= m_eblg.getDomain();


      getFlux(a_flux[idir], a_data[a_dit], ghostedBox, a_grid,
              m_eblg.getDomain(), m_eblg.getEBISL()[a_dit], m_dx, a_dit, idir);

    }
}
/***/
void
slowEBCO::
getFlux(EBFaceFAB&                    a_fluxCentroid,
        const EBCellFAB&              a_phi,
        const Box&                    a_ghostedBox,
        const Box&                    a_fabBox,
        const ProblemDomain&          a_domain,
        const EBISBox&                a_ebisBox,
        const Real&                   a_dx,
        const DataIndex&              a_datInd,
        const int&                    a_idir)
{
  CH_TIME("EBAMRPoissonOp::getFlux");
  //has some extra cells so...
  a_fluxCentroid.setVal(0.);
  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_fluxCentroid.nComp());
  Box cellBox = a_ghostedBox;
  //want only interior faces
  cellBox.grow(a_idir, 1);
  cellBox &= a_domain;
  cellBox.grow(a_idir,-1);

  Box faceBox = surroundingNodes(cellBox, a_idir);

  //make a EBFaceFAB (including ghost cells) that will hold centered gradients
  BaseFab<Real>& regFlux  = a_fluxCentroid.getSingleValuedFAB();
  const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
  const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());

  FORT_GETFLUXEBCO(CHF_FRA1(regFlux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(regPhi, 0),
                   CHF_BOX(faceBox),
                   CHF_CONST_REAL(a_dx),
                   CHF_CONST_INT(a_idir));

  EBFaceFAB fluxCenter(a_ebisBox, a_ghostedBox, a_idir,1);
  fluxCenter.copy(a_fluxCentroid);

  IntVectSet ivsCell = a_ebisBox.getIrregIVS(a_fabBox);
  if (!ivsCell.isEmpty())
    {
      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;

      for (FaceIterator faceit(ivsCell, a_ebisBox.getEBGraph(), a_idir,stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          Real phiHi = a_phi(face.getVoF(Side::Hi), 0);
          Real phiLo = a_phi(face.getVoF(Side::Lo), 0);
          Real fluxFace = bcoebff(face, 0)*(phiHi - phiLo)/a_dx;

          fluxCenter(face, 0) = fluxFace;
        }
      //interpolate from face centers to face centroids
      Box cellBox = a_fluxCentroid.getCellRegion();
      EBArith::interpolateFluxToCentroids(a_fluxCentroid,
                                          fluxCenter,
                                          a_fabBox,
                                          a_ebisBox,
                                          a_domain,
                                          a_idir);
    }

  a_fluxCentroid *= m_beta;
}
/***/
void
slowEBCO::
incrementFRFine(EBFastFR&             a_fluxReg,
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi,
                AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("slowEBCO::incrementFRFine");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(m_hasFine);
  int ncomp = 1;
  Interval interv(0,0);
  slowEBCO& finerEBAMROp = (slowEBCO& )(*a_finerOp);

  //ghost cells of phiFine need to be filled
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&) a_phiFine;
  finerEBAMROp.m_quadCFIWithCoar->interpolate(phiFine, a_phi, interv);
  phiFine.exchange(interv);

  DataIterator ditf = a_phiFine.dataIterator();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const Box&     boxFine = m_eblgFine.getDBL().get(ditf());
      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[ditf()];
      const EBCellFAB& phiFine = a_phiFine[ditf()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); sit.next())
            {
              Box fabBox = adjCellBox(boxFine, idir, sit(), 1);
              fabBox.shift(idir, -sign(sit()));

              Box ghostedBox = fabBox;
              ghostedBox.grow(1);
              ghostedBox.grow(idir,-1);
              ghostedBox &= m_eblgFine.getDomain();

              EBFaceFAB fluxFine(ebisBoxFine, ghostedBox, idir, ncomp);
              finerEBAMROp.getFlux(fluxFine, phiFine, ghostedBox, fabBox, m_eblgFine.getDomain(),
                                   ebisBoxFine, m_dxFine, ditf(), idir);

              Real scale = 1.0; //beta and bcoef already in flux

              a_fluxReg.incrementFineBoth(fluxFine, scale, ditf(), interv, idir, sit());
            }
        }
    }
}
/****/
void
slowEBCO::
getFlux(FArrayBox&                    a_flux,
        const FArrayBox&              a_phi,
        const Box&                    a_faceBox,
        const int&                    a_idir,
        const Real&                   a_dx,
        const DataIndex&              a_datInd)
{
  const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());
  FORT_GETFLUXEBCO(CHF_FRA1(a_flux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(a_phi, 0),
                   CHF_BOX(a_faceBox),
                   CHF_CONST_REAL(a_dx),
                   CHF_CONST_INT(a_idir));

  a_flux *= m_beta;
}

/****/
void
slowEBCO::
AMRResidualNC(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //dummy. there is no coarse when this is called
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);
  LevelData<EBCellFAB> phiC;
  AMRResidual(a_residual, a_phiFine, a_phi, phiC, a_rhs, a_homogeneousPhysBC, a_finerOp);
}
//
void
slowEBCO::
AMRResidualNF(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoar,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC)
{
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);

  AMROperatorNF(a_residual, a_phi, a_phiCoar,
                a_homogeneousPhysBC);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}
/****/
void
slowEBCO::
AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //dummy. there is no coarse when this is called
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);
  LevelData<EBCellFAB> phiC;
  AMROperator(a_LofPhi, a_phiFine, a_phi, phiC,
              a_homogeneousPhysBC, a_finerOp);
}
/***/
void
slowEBCO::
AMROperatorNF(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoar,
              bool a_homogeneousPhysBC)
{
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);

  applyOp(a_LofPhi,a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);

}
/****/
void
slowEBCO::
AMRRestrict(LevelData<EBCellFAB>&       a_resCoar,
            const LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_correction,
            const LevelData<EBCellFAB>& a_coarCorrection, 
            bool a_skip_res)
{
  CH_TIME("slowEBCO::AMRRestrict");
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_correction.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostCellsPhi);

  CH_assert(a_residual.nComp() == 1);
  CH_assert(a_resCoar.nComp() == 1);
  CH_assert(a_correction.nComp() == 1);

  LevelData<EBCellFAB> resThisLevel;
  bool homogeneousPhys = true;
  bool homogeneousCF =   false;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  resThisLevel.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);
  EBLevelDataOps::setVal(resThisLevel, 0.0);

  //API says that we must average(a_residual - L(correction, coarCorrection))
  applyOp(resThisLevel, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);
  incr(resThisLevel, a_residual, -1.0);
  scale(resThisLevel,-1.0);

  //use our nifty averaging operator
  Interval variables(0, 0);
  CH_assert(m_hasInterpAve);
  m_ebAverage.average(a_resCoar, resThisLevel, variables);

}
/**
   compute norm over all cells on coarse not covered by finer
**/
Real
slowEBCO::
AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
        const LevelData<EBCellFAB>& a_fineResid,
        const int& a_refRat,
        const int& a_ord)

{
  CH_TIME("slowEBCO::AMRNorm");
  int ncomp = a_coarResid.nComp();
  Interval interv(0, ncomp-1);
  const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();
  const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();
  CH_assert(coarGrids == m_eblg.getDBL());

  //create temp and zero out under finer grids
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  LevelData<EBCellFAB> coarTemp(coarGrids, ncomp, IntVect::Zero, ebcellfact);
  a_coarResid.copyTo(interv, coarTemp, interv);

  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& coarTempFAB = coarTemp[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];

      LayoutIterator litFine = fineGrids.layoutIterator();
      for (litFine.reset(); litFine.ok(); ++litFine)
        {
          Box overlayBox = coarTempFAB.box();
          Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat);

          overlayBox &= coarsenedGrid;
          if (!overlayBox.isEmpty())
            {
              BaseFab<Real>& regToZeroFAB = coarTempFAB.getSingleValuedFAB();
              FORT_AMRPZEROSUB(CHF_FRA(regToZeroFAB),
                               CHF_BOX(overlayBox),
                               CHF_INT(ncomp));

              IntVectSet ivsZero = ebisBox.getMultiCells(overlayBox);

              for (VoFIterator vofit(ivsZero, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                {
                  for (int ivar =0; ivar < ncomp; ivar++)
                    {
                      coarTempFAB(vofit(), ivar) = 0.0;
                    }
                }
            }
        }
    }
  //return norm of temp
  return norm(coarTemp, a_ord);
}
/***/
void
slowEBCO::
AMRProlong(LevelData<EBCellFAB>&       a_correction,
           const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("slowEBCO::AMRProlong");
  //use cached interpolation object
  Interval variables(0, 0);
  CH_assert(m_hasInterpAve);
  m_ebInterp.pwcInterp(a_correction, a_coarCorrection, variables);
}
/***/
void
slowEBCO::
AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                  const LevelData<EBCellFAB>& a_correction,
                  const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("slowEBCO::AMRUpdateResidual");
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_correction.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostCellsPhi);

  LevelData<EBCellFAB> lcorr;
  bool homogeneousPhys = true;
  bool homogeneousCF   = false;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  lcorr.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

  applyOp(lcorr, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);

  incr(a_residual, lcorr, -1);
}

#include "NamespaceFooter.H"
