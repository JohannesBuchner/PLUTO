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
#include "EBEllipticLoadBalance.H"
#include "EBArith.H"
#include "BRMeshRefine.H"
#include "EBAMRIO.H"

#include "EBAMRPoissonOp.H"
#include "EBQuadCFInterp.H"

#include "EBAMRPoissonOpF_F.H"
#include "InterpF_F.H"
#include "BCFunc.H"
#include "AMRPoissonOpF_F.H"
#include "CH_Timer.H"
#include "BCFunc.H"
#include "EBLevelGrid.H"
#include "EBAlias.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"


#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "parstream.H"
#endif

Real EBAMRPoissonOp::s_time = 0.;
bool EBAMRPoissonOp::s_doLazyRelax = false;
bool EBAMRPoissonOp::s_doInconsistentRelax = false;
bool EBAMRPoissonOp::s_doTrimEdges = false;
bool EBAMRPoissonOp::s_turnOffBCs = false; //REALLY needs to default to false
bool EBAMRPoissonOp::s_doEBEllipticLoadBalance = true; //false for MF
bool EBAMRPoissonOp::s_areaFracWeighted = false;
bool EBAMRPoissonOp::s_doSetListValueResid = false;
int EBAMRPoissonOp::s_numComps = 1;
int EBAMRPoissonOp::s_whichComp = 0;

void
EBAMRPoissonOp::doLazyRelax(bool a_doLazyRelax)
{
  s_doLazyRelax = a_doLazyRelax;
}
void
EBAMRPoissonOp::doEBEllipticLoadBalance(bool a_doEBEllipticLoadBalance)
{
  s_doEBEllipticLoadBalance = a_doEBEllipticLoadBalance;
}
void
EBAMRPoissonOp::areaFracWeighted(bool a_areaFracWeighted)
{
  s_areaFracWeighted = a_areaFracWeighted;
}
//////////////
void
EBAMRPoissonOp::setAlphaAndBeta(const Real& a_alpha,
                                const Real& a_beta)
{
  CH_TIME("EBAMRPoissonOp::setAlphaAndBeta");
  m_alpha = a_alpha*m_aCoef;
  m_beta  = a_beta*m_bCoef;
}
void
EBAMRPoissonOp::diagonalScale(LevelData<EBCellFAB>& a_rhs,
                              bool                  a_kappaWeighted)
{
  CH_TIME("EBAMRPoissonOp::kappaWeight");
  if (a_kappaWeighted)
    {
      EBLevelDataOps::kappaWeight(a_rhs);
    }
  if (s_areaFracWeighted)
    {
      EBLevelDataOps::areaFracScalingWeight(a_rhs);
    }
}

EBAMRPoissonOp::~EBAMRPoissonOp()
{
}

EBAMRPoissonOp::EBAMRPoissonOp()
{
}

//////////////
EBAMRPoissonOp::
EBAMRPoissonOp(const EBLevelGrid &                  a_eblgFine,
               const EBLevelGrid &                  a_eblg,
               const EBLevelGrid &                  a_eblgCoar,
               const EBLevelGrid &                  a_eblgCoarMG,
               const RefCountedPtr<EBQuadCFInterp>& a_quadCFI,
               const RefCountedPtr<BaseDomainBC>&   a_domainBC,
               const RefCountedPtr<BaseEBBC>&       a_ebBC,
               const RealVect&                      a_dx,
               const RealVect&                      a_dxCoar,
               const RealVect&                      a_origin,
               const int&                           a_refToFine,
               const int&                           a_refToCoar,
               const bool&                          a_hasFine,
               const bool&                          a_hasCoar,
               const bool&                          a_hasMGObjects,
               const bool&                          a_layoutChanged,
               const int&                           a_numPreCondIters,
               const int&                           a_relaxType,
               const Real&                          a_alpha,
               const Real&                          a_beta,
               const IntVect&                       a_ghostCellsPhi,
               const IntVect&                       a_ghostCellsRHS,
               int                                  a_testRef)
  : m_testRef( a_testRef ),
    m_ghostCellsPhi( a_ghostCellsPhi ),
    m_ghostCellsRHS( a_ghostCellsRHS )
{
  CH_TIME("EBAMRPoissonOp::EBAMRPoissonOp");
  m_quadCFIWithCoar = a_quadCFI;
  m_hasFine        = a_hasFine;
  m_hasCoar        = a_hasCoar;
  m_numPreCondIters= a_numPreCondIters;
  m_relaxType      = a_relaxType;
  m_eblg            = a_eblg;
  m_domainBC       = a_domainBC;
  m_ebBC           = a_ebBC;
  m_dx             = a_dx;
  m_origin         = a_origin;
  m_alpha          = a_alpha;
  m_aCoef      = a_alpha;
  m_beta           = a_beta;
  m_bCoef       = a_beta;
  m_hasInterpAve = false;
  m_dxCoar         = a_dxCoar;
  m_hasMGObjects = a_hasMGObjects;
  m_layoutChanged = a_layoutChanged;

  //pre-compute 1/dx and 1/(dx^2)
  m_invDx  = 1.0/m_dx;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_invDx2[idir] = m_invDx[idir]*m_invDx[idir];
    }

  defineWithCoarser(a_eblgCoar, a_refToCoar);
  defineWithFiner(  a_eblgFine, a_refToFine);
  defineMGObjects(a_eblgCoarMG);

  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects

  m_exchangeCopier.define(m_eblg.getDBL(), m_eblg.getDBL(), m_ghostCellsPhi,  true);
  if (s_doTrimEdges)
  {
    m_exchangeCopier.trimEdges(m_eblg.getDBL(), m_ghostCellsPhi);
  }

  //define stencils for the operator
  defineStencils();
  /*
    Uncomment this in order to dump a matrix representation of the
    elliptic operator.
  */
  // dumpStencilMatrix();
}
/******/
void
EBAMRPoissonOp::
defineMGObjects(const EBLevelGrid& a_eblgCoarMG)
{
  if (m_hasMGObjects)
    {
      int mgRef = 2;
      int ncomp =1 ;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain(), mgRef, ncomp, m_eblg.getEBIS(),
                           m_ghostCellsPhi, m_layoutChanged);

      m_ebAverageMG.define(m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain() , mgRef, ncomp, m_eblg.getEBIS(),
                           m_ghostCellsRHS, m_layoutChanged);

    }
}
/********/
void
EBAMRPoissonOp::
defineWithCoarser(const EBLevelGrid& a_eblgCoar, const int& a_refToCoar)
{
  if (m_hasCoar)
    {
      int ncomp = 1;
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
                                                          maxBoxSize, dumbool,
                                                          m_testRef);

      //should follow from coarsenable
      if (hasCoarser)
        {
          m_eblgCoarsenedFine = EBLevelGrid(dblCoarsenedFine, domainCoarsenedFine, 4, m_eblg.getEBIS());
          m_hasInterpAve = true;
          m_ebInterp.define( m_eblg.getDBL(),     m_eblgCoar.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoar.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             m_ghostCellsPhi);
          m_ebAverage.define(m_eblg.getDBL(),     m_eblgCoarsenedFine.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoarsenedFine.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             m_ghostCellsRHS);

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
}
/******/
void
EBAMRPoissonOp::
defineWithFiner(const EBLevelGrid& a_eblgFine,
                const int& a_refToFine)
{
  if (m_hasFine)
    {
      int ncomp = 1;
      m_dxFine         = m_dx/a_refToFine;
      m_refToFine      = a_refToFine;
      m_eblgFine       = a_eblgFine;
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp);
      m_hasEBCF = m_fastFR.hasEBCF();

      if (m_hasEBCF)
        {
          defineEBCFStencils();
        }
    }
  else
    {
      m_refToFine    = 1;
    }

}
/******/
void
EBAMRPoissonOp::
defineEBCFStencils()
{
  ///this routine is ugly and complicated.
  //I will attempt to comment it but I fear it is a lost cause
  //because the algorithm is so arcane.
  //We are attempting to only do stuff at the very specific
  //points where there is an embedded boundary crossing a coarse-fine
  //interface.   We happen to know that EBFastFR already has done this
  //choreography and we want to leverage it.

  //EBFastFR has data structures in it that serve as buffers and so on
  //that we will (thankfully) be able to leave alone.   We are only
  //going to access the data structures wherein it identifies which
  //coarse cells are associated with the coarse-fine interface
  // where the EB crosses and use this list to build up  which faces
  // need to be cal

  //important factioid: beta gets multiplied in at the last minute
  //(on evaluation) because it can change as diffusion solvers progress.
  if (m_hasFine && m_hasEBCF)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              //coarse fine stuff is between me and next finer level
              //fine stuff lives over m_eblgfine
              //coar stuff lives over m_eblg
              int index = m_fastFR.index(idir, sit());
              m_stencilCoar[index].define(m_eblg.getDBL());
              m_faceitCoar [index].define(m_eblg.getDBL());

              for (DataIterator dit =      m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
                {
                  Vector<FaceIndex>& facesEBCFCoar =  m_faceitCoar[index][dit()];
                  Vector<VoFStencil>& stencEBCFCoar= m_stencilCoar[index][dit()];
                  Vector<VoFIterator>& vofitlist = m_fastFR.getVoFItCoar(dit(), idir, sit());
                  //first build up the list of the faces
                  for (int ivofit = 0; ivofit < vofitlist.size(); ivofit++)
                    {
                      VoFIterator& vofit = vofitlist[ivofit];
                      for (vofit.reset(); vofit.ok(); ++vofit)
                        {
                          //on the coarse side of the CF interface we are
                          //looking in the flip direction because we look
                          //back at the interface
                          Vector<FaceIndex> facespt = m_eblg.getEBISL()[dit()].getFaces(vofit(), idir, flip(sit()));
                          facesEBCFCoar.append(facespt);
                        }
                    }

                  stencEBCFCoar.resize(facesEBCFCoar.size());
                  const EBISBox& ebisBoxCoar = m_eblg.getEBISL()[dit()];
                  for (int iface = 0; iface < stencEBCFCoar.size(); iface++)
                    {
                      IntVectSet cfivs; //does not apply here
                      EBAMRPoissonOp::getFluxStencil(stencEBCFCoar[iface], facesEBCFCoar[iface],
                                                     ebisBoxCoar, cfivs, m_dx, true);

                    }
                }
            }
        }
    }
}
/******/
void
EBAMRPoissonOp::
defineStencils()
{
  CH_TIME("EBAMRPoissonOp::defineStencils");
  // create ebstencil for irregular applyOp
  m_opEBStencil.define(m_eblg.getDBL());
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_opEBStencilInhomDomLo[idir].define(m_eblg.getDBL());
      m_opEBStencilInhomDomHi[idir].define(m_eblg.getDBL());
    }
  m_invDiagEBStencil.define(m_eblg.getDBL());
 // create vofstencils for applyOp and
  LayoutData<BaseIVFAB<VoFStencil> >&  opStencil = m_opStencil;
  m_rhsSetList.define(m_eblg.getDBL());
  LayoutData<BaseIVFAB<VoFStencil> >  relStencil;
  LayoutData<BaseIVFAB<VoFStencil> >  invDiagStencil;
  opStencil.define(m_eblg.getDBL());
  relStencil.define(m_eblg.getDBL());
  invDiagStencil.define(m_eblg.getDBL());
  //define bc stencils and create eb flux stencil, beta stays out
  //and gets multiplied in later, so does boundary area and dx
  m_ebBC->define((*m_eblg.getCFIVS()), 1./m_dx[0]);
  LayoutData<BaseIVFAB<VoFStencil> >* ebFluxStencil = m_ebBC->getFluxStencil(0);

  m_alphaDiagWeight.define(  m_eblg.getDBL());
  m_betaDiagWeight.define(   m_eblg.getDBL());
  m_one.define(   m_eblg.getDBL());
  m_vofItIrreg.define( m_eblg.getDBL()); // vofiterator cache

  Box domainBox = m_eblg.getDomain().domainBox();
  Box sideBoxLo[SpaceDim];
  Box sideBoxHi[SpaceDim];

  m_cacheInhomDomBCLo.resize(s_numComps);
  m_cacheInhomDomBCHi.resize(s_numComps);
  for (int icomp = 0; icomp < s_numComps; icomp++)
    {
      m_cacheInhomDomBCLo[icomp].resize(SpaceDim);
      m_cacheInhomDomBCHi[icomp].resize(SpaceDim);
    }
  
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_vofItIrregDomLo[idir].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofItIrregDomHi[idir].define( m_eblg.getDBL()); // vofiterator cache for domain hi
      for (int icomp = 0; icomp < s_numComps; icomp++)
        {
          m_cacheInhomDomBCLo[icomp][idir] = RefCountedPtr<LayoutData<EBCellFAB> >(new LayoutData<EBCellFAB>(m_eblg.getDBL()));
          m_cacheInhomDomBCHi[icomp][idir] = RefCountedPtr<LayoutData<EBCellFAB> >(new LayoutData<EBCellFAB>(m_eblg.getDBL()));
        }
    }

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& curBox = m_eblg.getDBL().get(dit());
      const EBISBox& curEBISBox = m_eblg.getEBISL()[dit()];
      const EBGraph& curEBGraph = curEBISBox.getEBGraph();
      Vector<VolIndex>& rhsSetList = m_rhsSetList[dit()];

      IntVectSet notRegular = curEBISBox.getIrregIVS(curBox);

      BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = opStencil[dit()];
      BaseIVFAB<VoFStencil>& relStencilBaseIVFAB = relStencil[dit()];
      BaseIVFAB<VoFStencil>& invDiagStencilBaseIVFAB = invDiagStencil[dit()];
      BaseIVFAB<VoFStencil> opStencilDomLoBaseIVFAB[SpaceDim];
      BaseIVFAB<VoFStencil> opStencilDomHiBaseIVFAB[SpaceDim];
      BaseIVFAB<VoFStencil> opStencilInhomDomLoBaseIVFAB[SpaceDim];
      BaseIVFAB<VoFStencil> opStencilInhomDomHiBaseIVFAB[SpaceDim];
      BaseIVFAB<Real>&       alphaWeight  = m_alphaDiagWeight[dit()];
      BaseIVFAB<Real>&        betaWeight  = m_betaDiagWeight[dit()];
      BaseIVFAB<Real>&        one  = m_one[dit()];
      curStencilBaseIVFAB.define(notRegular,curEBGraph, 1);
      relStencilBaseIVFAB.define(notRegular,curEBGraph, 1);
      invDiagStencilBaseIVFAB.define(notRegular,curEBGraph, 1);
      alphaWeight.define(        notRegular,curEBGraph, 1);
      betaWeight.define(         notRegular,curEBGraph, 1);
      one.define(         notRegular,curEBGraph, 1);

      //cache the vofIterators
      m_vofItIrreg[dit()].define(notRegular,curEBISBox.getEBGraph());

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet loIrreg = notRegular;
          IntVectSet hiIrreg = notRegular;
          loIrreg &= sideBoxLo[idir];
          hiIrreg &= sideBoxHi[idir];
          m_vofItIrregDomLo[idir][dit()].define(loIrreg, curEBGraph);
          m_vofItIrregDomHi[idir][dit()].define(hiIrreg, curEBGraph);
          opStencilDomLoBaseIVFAB[idir].define(loIrreg, curEBGraph, 1);
          opStencilDomHiBaseIVFAB[idir].define(hiIrreg, curEBGraph, 1);
          opStencilInhomDomLoBaseIVFAB[idir].define(loIrreg, curEBGraph, 1);
          opStencilInhomDomHiBaseIVFAB[idir].define(hiIrreg, curEBGraph, 1);

          for (int icomp = 0; icomp < s_numComps; icomp++)
            {
              (*m_cacheInhomDomBCLo[icomp][idir])[dit()].define(curEBISBox,curBox,1);
              (*m_cacheInhomDomBCHi[icomp][idir])[dit()].define(curEBISBox,curBox,1);
            }
         }

      VoFIterator& vofit = m_vofItIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();

          VoFStencil& curStencil = curStencilBaseIVFAB(VoF,0);
          VoFStencil& relStencil = relStencilBaseIVFAB(VoF,0);
          VoFStencil& invDiagStencil = invDiagStencilBaseIVFAB(VoF,0);
          invDiagStencil.clear();

          getDivFStencil(curStencil,VoF, dit(), true);
          getDivFStencil(relStencil,VoF, dit(), !s_doInconsistentRelax);

          Real& curAlphaWeight  = alphaWeight(VoF,0);
          Real& curBetaWeight   =  betaWeight(VoF,0);
          Real& curOne   =  one(VoF,0);

          const Real kappa = curEBISBox.volFrac(VoF);

          curOne = 1.;

          curAlphaWeight = kappa;
          if (s_areaFracWeighted)
            {
              curAlphaWeight *= curEBISBox.areaFracScaling(VoF);
            }

          curBetaWeight = EBArith::getDiagWeight(relStencil, VoF);
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

          //debug take out ebstencil for diagnostic
          //ebFluxStencil = NULL;
          //end debug
          if (ebFluxStencil != NULL)
            {
              BaseIVFAB<VoFStencil>& ebFluxStencilBaseIVFAB = (*ebFluxStencil)[dit()];
              //this fills the stencil with the gradient
              const VoFStencil&  ebFluxStencilPt = ebFluxStencilBaseIVFAB(VoF,0);
              const Real boundaryArea = curEBISBox.bndryArea(VoF);
              //if the stencil returns empty, this means that our
              //geometry is underresolved and we just set the
              //stencil to zero.   This might not work.
              if (ebFluxStencilPt.size() == 0)
                {
                  if (s_doSetListValueResid && boundaryArea > 1.e-12)
                    {
                      curStencil *= 0.0;
                      relStencil *= 0.0;
                      rhsSetList.push_back(VoF);
                    }
                }
              else
                {
                  curStencil += ebFluxStencilPt;
                  relStencil += ebFluxStencilPt;
                }
            }
        }//vofitIrreg

      //add in the homogeneous part of stencil when EB x domain and cache the inhomogeneous part
      BaseIVFAB<VoFStencil>& opStencilBaseIVFAB = opStencil[dit()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          VoFIterator& vofitlo = m_vofItIrregDomLo[idir][dit()];
          for (vofitlo.reset(); vofitlo.ok(); ++vofitlo)
            {
              const VolIndex& vof = vofitlo();

              VoFStencil& opStencilDomLo = opStencilDomLoBaseIVFAB[idir](vof,0);
              opStencilDomLo.clear();
              m_domainBC->getFluxStencil(opStencilDomLo,
                                         vof,
                                         0,//comp
                                         m_dx,
                                         idir,
                                         Side::Lo,
                                         curEBISBox);

              VoFStencil& opStencil = opStencilBaseIVFAB(vof,0);
              opStencil += opStencilDomLo;

              VoFStencil& opStencilInhomDomLo = opStencilInhomDomLoBaseIVFAB[idir](vof,0);
              opStencilInhomDomLo.clear();

              for (int icomp = 0; icomp < s_numComps; icomp++)
                {
                  Real flux;
                  m_domainBC->getInhomFaceFlux(flux,
                                               vof,
                                               icomp,//comp
                                               (*m_cacheInhomDomBCLo[icomp][idir])[dit()],
                                               m_origin,
                                               m_dx,
                                               idir,
                                               Side::Lo,
                                               dit(),
                                               s_time);

                  (*m_cacheInhomDomBCLo[icomp][idir])[dit()](vof, 0) = flux;
                }
            }//vofitlo

          VoFIterator& vofithi = m_vofItIrregDomHi[idir][dit()];
          for (vofithi.reset(); vofithi.ok(); ++vofithi)
            {
              const VolIndex& vof = vofithi();

              VoFStencil& opStencilDomHi = opStencilDomHiBaseIVFAB[idir](vof,0);
              opStencilDomHi.clear();
              m_domainBC->getFluxStencil(opStencilDomHi,
                                         vof,
                                         0,//comp
                                         m_dx,
                                         idir,
                                         Side::Hi,
                                         curEBISBox);

              VoFStencil& opStencil = opStencilBaseIVFAB(vof,0);
              opStencil += opStencilDomHi;

              VoFStencil& opStencilInhomDomHi = opStencilInhomDomHiBaseIVFAB[idir](vof,0);
              opStencilInhomDomHi.clear();

              for (int icomp = 0; icomp < s_numComps; icomp++)
                {
                  Real flux;
                  m_domainBC->getInhomFaceFlux(flux,
                                               vof,
                                               icomp,//comp
                                               (*m_cacheInhomDomBCHi[icomp][idir])[dit()],
                                               m_origin,
                                               m_dx,
                                               idir,
                                               Side::Hi,
                                               dit(),
                                               s_time);

                  (*m_cacheInhomDomBCHi[icomp][idir])[dit()](vof, 0) = flux;
                }
            }//vofithi
        }//idir

      //Operator ebstencil
      m_opEBStencil[dit()] =
        RefCountedPtr<EBStencil>(new EBStencil(m_vofItIrreg[dit()].getVector(),
                                               opStencil[dit()],
                                               m_eblg.getDBL().get(dit()),
                                               m_eblg.getEBISL()[dit()],
                                               m_ghostCellsPhi,
                                               m_ghostCellsRHS,
                                               0,
                                               true));

      //stencil for inverse diagonal weight precond scaling
      m_invDiagEBStencil[dit()] =
        RefCountedPtr<EBStencil>(new EBStencil(m_vofItIrreg[dit()].getVector(),
                                               invDiagStencil[dit()],
                                               m_eblg.getDBL().get(dit()),
                                               m_eblg.getEBISL()[dit()],
                                               m_ghostCellsRHS,//constant rhs
                                               m_ghostCellsPhi,//changing lhs
                                               0,
                                               true));

      //stencils for inhomogeneous domain bcs in irreg cells
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_opEBStencilInhomDomLo[idir][dit()] =
            RefCountedPtr<EBStencil>(new EBStencil(m_vofItIrregDomLo[idir][dit()].getVector(),
                                                   opStencilInhomDomLoBaseIVFAB[idir],
                                                   m_eblg.getDBL().get(dit()),
                                                   m_eblg.getEBISL()[dit()],
                                                   IntVect::Zero,
                                                   m_ghostCellsRHS,
                                                   0,
                                                   true));

          m_opEBStencilInhomDomHi[idir][dit()] =
            RefCountedPtr<EBStencil>(new EBStencil(m_vofItIrregDomHi[idir][dit()].getVector(),
                                                   opStencilInhomDomHiBaseIVFAB[idir],
                                                   m_eblg.getDBL().get(dit()),
                                                   m_eblg.getEBISL()[dit()],
                                                   IntVect::Zero,
                                                   m_ghostCellsRHS,
                                                   0,
                                                   true));
        }
   }//dit

  // create ebstencils for each color
  // Loop over all possibilities (in all dimensions)
  // create and define colored stencils (2 parts)
  LayoutData<BaseIVFAB<VoFStencil> > colorStencil;
  colorStencil.define(m_eblg.getDBL());
  if (m_relaxType == 1 || m_relaxType == 2 || m_relaxType == 3 || m_relaxType == 4 || m_relaxType == 999)
    {
      EBArith::getMultiColors(m_colors);
    }
  else if (m_relaxType == 0)
    {
      m_colors.resize(0);
    }
  else
    {
      MayDay::Error("EBAMRPoissonOp::defineStencils -- invalid relaxation type");
    }
  for (int icolor=0; icolor < m_colors.size(); ++icolor)
    {
      m_colorEBStencil[icolor].define(m_eblg.getDBL());
      m_vofItIrregColor[icolor].define( m_eblg.getDBL());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_colorEBStencilDomLo[icolor][idir].define(m_eblg.getDBL());
          m_colorEBStencilDomHi[icolor][idir].define(m_eblg.getDBL());
          m_vofItIrregColorDomLo[icolor][idir].define( m_eblg.getDBL());
          m_vofItIrregColorDomHi[icolor][idir].define( m_eblg.getDBL());
        }
      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& curEBISBox = m_eblg.getEBISL()[dit()];
          const EBGraph& curEBGraph = curEBISBox.getEBGraph();
          Box dblBox( m_eblg.getDBL().get(dit()) );

          IntVectSet ivsColor(DenseIntVectSet(dblBox, false));

          BaseIVFAB<VoFStencil>& colorStencilBaseIVFAB = colorStencil[dit()];
          BaseIVFAB<VoFStencil>& relStencilBaseIVFAB = relStencil[dit()];
          BaseIVFAB<VoFStencil> colorStencilDomLoBaseIVFAB[SpaceDim];
          BaseIVFAB<VoFStencil> colorStencilDomHiBaseIVFAB[SpaceDim];

          VoFIterator& vofit = m_vofItIrreg[dit()];
          for (vofit.reset(); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const IntVect& iv = vof.gridIndex();

              bool doThisVoF = true;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (iv[idir] % 2 != m_colors[icolor][idir])
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

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              IntVectSet loIrregColor = ivsColor;
              IntVectSet hiIrregColor = ivsColor;
              loIrregColor &= sideBoxLo[idir];
              hiIrregColor &= sideBoxHi[idir];
              colorStencilDomLoBaseIVFAB[idir].define(loIrregColor, curEBGraph, 1);
              colorStencilDomHiBaseIVFAB[idir].define(hiIrregColor, curEBGraph, 1);
              m_vofItIrregColorDomLo[icolor][idir][dit()].define(loIrregColor,curEBGraph);
              m_vofItIrregColorDomHi[icolor][idir][dit()].define(hiIrregColor,curEBGraph);
            }

          VoFIterator& vofitcolor = m_vofItIrregColor[icolor][dit()];
          for (vofitcolor.reset(); vofitcolor.ok(); ++vofitcolor)
            {
              const VolIndex& vof = vofitcolor();

              VoFStencil& colorStencil =   colorStencilBaseIVFAB(vof,0);
              const VoFStencil& relStencil = relStencilBaseIVFAB(vof,0);
              colorStencil = relStencil;
            }

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              VoFIterator& vofitcolorlo = m_vofItIrregColorDomLo[icolor][idir][dit()];
              for (vofitcolorlo.reset(); vofitcolorlo.ok(); ++vofitcolorlo)
                {
                  const VolIndex& vof = vofitcolorlo();

                  VoFStencil& colorStencilDomLo = colorStencilDomLoBaseIVFAB[idir](vof,0);
                  colorStencilDomLo.clear();
                  m_domainBC->getFluxStencil(colorStencilDomLo,
                                             vof,
                                             0,//comp
                                             m_dx,
                                             idir,
                                             Side::Lo,
                                             curEBISBox);

                  VoFStencil& colorStencil =   colorStencilBaseIVFAB(vof,0);
                  colorStencil += colorStencilDomLo;
                }

              VoFIterator& vofitcolorhi = m_vofItIrregColorDomHi[icolor][idir][dit()];
              for (vofitcolorhi.reset(); vofitcolorhi.ok(); ++vofitcolorhi)
                {
                  const VolIndex& vof = vofitcolorhi();

                  VoFStencil& colorStencilDomHi = colorStencilDomHiBaseIVFAB[idir](vof,0);
                  colorStencilDomHi.clear();
                  m_domainBC->getFluxStencil(colorStencilDomHi,
                                             vof,
                                             0,//comp
                                             m_dx,
                                             idir,
                                             Side::Hi,
                                             curEBISBox);

                  VoFStencil& colorStencil =   colorStencilBaseIVFAB(vof,0);
                  colorStencil += colorStencilDomHi;
                }
            }

          Vector<VolIndex> srcVofs = m_vofItIrregColor[icolor][dit()].getVector();

          m_colorEBStencil[icolor][dit()]  =
            RefCountedPtr<EBStencil>(new EBStencil(srcVofs,
                                                   colorStencil[dit()],
                                                   m_eblg.getDBL().get(dit()),
                                                   m_eblg.getEBISL()[dit()],
                                                   m_ghostCellsPhi,
                                                   m_ghostCellsRHS,
                                                   0,
                                                   true));

        }//dit
    }//color
}

//version that does not fill ebislCoar
bool EBAMRPoissonOp::
getCoarserLayouts(DisjointBoxLayout&       a_dblCoar,
                  ProblemDomain&           a_domainCoar,
                  const DisjointBoxLayout& a_dblFine,
                  const EBISLayout&        a_ebislFine,
                  const ProblemDomain&     a_domainFine,
                  int                      a_refToCoar,
                  const EBIndexSpace*      a_ebisPtr,
                  int                      a_maxBoxSize,
                  bool&                    a_layoutChanged,
                  int                      a_testRef)
{
  CH_TIME("EBArith::getCoarserLayouts(with ebisl)");
  bool turnOffAgglomeration = false;
  ParmParse pp;
  if (pp.contains("turn_off_agglomeration"))
    {
      pp.get("turn_off_agglomeration", turnOffAgglomeration);
    }
  // int mgProcFactor = 1;
  // if (pp.contains("mg_proc_factor"))
  //   {
  //     pp.get("mg_proc_factor", mgProcFactor);
  //   }
  // check to see if domain is coarsenable by this amount
  //never want to coarsen down to 1x1
  int testRef = a_testRef*a_refToCoar;
  ProblemDomain testBox = coarsen(a_domainFine, testRef);
  testBox.refine(testRef);
  a_layoutChanged = true;
  if (testBox != a_domainFine)
    {
      //not coarsenable.
      return false;
    }

  a_domainCoar = coarsen(a_domainFine, a_refToCoar);

  int fac = 2;

  if ((a_dblFine.coarsenable(fac*a_refToCoar) && a_dblFine.isClosed() && (a_dblFine.size() > 0)))
    {
      //if we can, just coarsen the grids
      coarsen(a_dblCoar, a_dblFine, a_refToCoar);
      a_layoutChanged = false;
    }
  else if ((!turnOffAgglomeration) && (a_dblFine.size() != 1))
    {//check to see if we are covering all uncovered IntVects with multiple boxes
      unsigned long long numPtsLeft = a_domainFine.domainBox().numPts();
      for (LayoutIterator lit = a_dblFine.layoutIterator(); lit.ok(); ++lit)
        {
          unsigned long long  ptsGrid = a_dblFine[lit()].numPts();
          numPtsLeft -= ptsGrid;
        }
      if (numPtsLeft == 0)
         {
          //if we are covering domain with more than one box, make one big split-up box
          Vector<Box> boxes;
          domainSplit(a_domainCoar, boxes, a_maxBoxSize);
          mortonOrdering(boxes);

          Vector<int> procs;
          if (s_doEBEllipticLoadBalance)
            {
              EBEllipticLoadBalance(procs, boxes, a_domainCoar, false, a_ebisPtr);
            }
          else
            {
              LoadBalance(procs, boxes);
            }
          // for (int iproc; iproc < procs.size(); ++iproc)
          //   {
          //     if ((mgProcFactor*iproc)+1 <= procs.size()+1)
          //       procs[iproc] = mgProcFactor*iproc;
          //   }
          a_dblCoar.define(boxes, procs, a_domainCoar);
        }
      else
        {
          //we are out of ideas
          return false;
        }
    }
  else
    {
      //we are out of ideas
      return false;
    }

  //if we got here, then we have coarser stuff
  return true;
}
/***/
void
EBAMRPoissonOp::
getDivFStencil(VoFStencil&      a_vofStencil,
               const VolIndex&  a_vof,
               const DataIndex& a_dit,
               bool             a_doFaceInterp)
{
  CH_TIME("EBAMRPoissonOp::getDivFStencil");
  const EBISBox& ebisBox  =    m_eblg.getEBISL()[a_dit];
  const IntVectSet& cfivs = (*m_eblg.getCFIVS())[a_dit];

  EBAMRPoissonOp::getDivFStencil(a_vofStencil, a_vof, ebisBox, cfivs, m_dx, a_doFaceInterp);
}
/***/
void
EBAMRPoissonOp::
getDivFStencil(VoFStencil&       a_vofStencil,
               const VolIndex&   a_vof,
               const EBISBox&    a_ebisBox,
               const IntVectSet& a_cfivs,
               const RealVect&   a_dx,
               bool a_doFaceInterp)
{
  a_vofStencil.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real dx =   a_dx[idir];
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = a_ebisBox.getFaces(a_vof, idir, sit());
          for (int iface = 0; iface < faces.size(); iface++)
            {
              VoFStencil fluxStencil;
              EBAMRPoissonOp::getFluxStencil(fluxStencil, faces[iface],
                                             a_ebisBox, a_cfivs, a_dx,
                                             a_doFaceInterp);

              Real areaFrac = a_ebisBox.areaFrac(faces[iface]);
              fluxStencil *= Real(isign)*areaFrac/dx;
              a_vofStencil += fluxStencil;
            }
        }
    }
  if (s_areaFracWeighted)
    {
      a_vofStencil *= a_ebisBox.areaFracScaling(a_vof);
    }
}
/***/
void
EBAMRPoissonOp::
getFluxStencil(VoFStencil&      a_fluxStencil,
               const FaceIndex& a_face,
               const EBISBox&   a_ebisBox,
               const IntVectSet& a_cfivs,
               const RealVect&   a_dx,
               bool              a_doFaceInterp)
{

  if (a_doFaceInterp)
    {
      //need to do this by interpolating to centroids
      //so get the stencil at each face center and add with
      //interpolation weights
      FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                         a_cfivs,
                                                         a_ebisBox,
                                                         a_ebisBox.getDomain());
      a_fluxStencil.clear();
      for (int isten = 0; isten < interpSten.size(); isten++)
        {
          const FaceIndex& face = interpSten.face(isten);
          const Real&    weight = interpSten.weight(isten);
          VoFStencil faceCentSten;
          EBAMRPoissonOp::getFaceCenteredFluxStencil(faceCentSten, face, a_dx);
          faceCentSten *= weight;
          a_fluxStencil += faceCentSten;
        }
    }
  else
    {
      a_fluxStencil.clear();
      EBAMRPoissonOp::getFaceCenteredFluxStencil(a_fluxStencil, a_face, a_dx);
    }
}
/***/
void
EBAMRPoissonOp::
getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                           const FaceIndex& a_face,
                           const RealVect& a_dx)
{
  //CH_TIME("EBAMRPoissonOp::getFaceCenteredFluxStencil");
  //face centered gradient is just a centered diff
  int faceDir= a_face.direction();
  a_fluxStencil.clear();

  Real dx = a_dx[faceDir];
  if (!a_face.isBoundary())
    {
      a_fluxStencil.add(a_face.getVoF(Side::Hi),  1.0/dx, 0);
      a_fluxStencil.add(a_face.getVoF(Side::Lo), -1.0/dx, 0);
    }
  else
    {
      //the boundary condition handles this one.
    }
}

void EBAMRPoissonOp::
residual(LevelData<EBCellFAB>&       a_residual,
         const LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         bool                        a_homogeneousPhysBC)
{
  CH_TIME("EBAMRPoissonOp::residual");
  //this is a multigrid operator so only homogeneous CF BC
  //and null coar level
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);

  if (s_doSetListValueResid && m_alpha == 1)
    {
      setListValue(a_rhs, 0);
      applyOp(a_residual,a_phi, NULL, a_homogeneousPhysBC, true);
      axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
      setListValue(a_residual, 0);
    }
  else
    {
      applyOp(a_residual,a_phi, NULL, a_homogeneousPhysBC, true);
      axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
    }
}

void EBAMRPoissonOp::
setListValue(const LevelData<EBCellFAB>& a_data, Real a_value)
{    
  LevelData<EBCellFAB>& data = const_cast<LevelData<EBCellFAB>&>(a_data);
  for (DataIterator dit = data.dataIterator() ; dit.ok() ; ++dit )
    {
      EBCellFAB& fab = data[dit()];
      const Vector<VolIndex>& setList = m_rhsSetList[dit()];
      setListValue(fab, setList, a_value);
    }
}

void EBAMRPoissonOp::
setListValue(EBCellFAB& a_fab, const Vector<VolIndex>& a_setList, Real a_value)
{
  for (int i = 0; i < a_setList.size(); i++ )
    {
      const VolIndex& vol = a_setList[i];
      a_fab(vol,0) = a_value;
    }
}

void EBAMRPoissonOp::
setRhsSetList(const LayoutData<Vector<VolIndex> >& a_list)
{   
  for (DataIterator dit = a_list.dataIterator(); dit.ok(); ++dit)
    {
      m_rhsSetList[dit()] = a_list[dit()];
    }
}

void EBAMRPoissonOp::
getOpMatrix(const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_rhs)
{
  LevelData<EBCellFAB> column, phi01, rhs0;
  LevelData<BaseFab<int> > gids;

  IntVect idghosts = a_rhs.ghostVect();
  const DisjointBoxLayout &dbl = a_rhs.disjointBoxLayout();
  const EBLevelGrid eblg = this->getEBLG();
  const EBCellFactory ebcellfact(eblg.getEBISL() );
  column.define(dbl, 1, idghosts, ebcellfact);
  rhs0.define(dbl, 1, idghosts, ebcellfact);
  phi01.define(dbl, 1, idghosts, ebcellfact);
  gids.define(dbl,1,idghosts);

  EBLevelDataOps::setVal(rhs0,  0);
  EBLevelDataOps::setVal(phi01, 0);
  EBLevelDataOps::setVal(column,0);

  int data=0;
  for (DataIterator dit = a_rhs.dataIterator() ; dit.ok() ; ++dit )
    {
      const Box &box = dbl.get(dit());
      if (CH_SPACEDIM==3) data += box.size(0)*box.size(1)*box.size(2);
      else data += box.size(0)*box.size(1);
    }
  int gid0;
  int ierr;
#ifdef CH_MPI
  int result;
  MPI_Comm wcomm = Chombo_MPI::comm;
  ierr = MPI_Scan( &data, &result, 1, MPI_INT, MPI_SUM, wcomm );
  gid0 = result - data;
#ifdef CH_USE_PETSC
  PETSC_COMM_WORLD=Chombo_MPI::comm;
#endif
#else
  gid0 = 0;
  ierr = 0;
#endif
  int gid = gid0;

  for (DataIterator  dit = a_rhs.dataIterator() ; dit.ok() ; ++dit )
    {
      BaseFab<int> &gidsFab = gids[dit()];
      const Box& box = dbl.get(dit());
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next(), gid++ )
        {
          IntVect iv = bit();
          gidsFab(iv,0) = gid;
        }
    }

#ifdef CH_USE_PETSC
  Mat petsc_mat;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD,data,data,PETSC_DECIDE,PETSC_DECIDE,
                         8,PETSC_NULL,2,PETSC_NULL,&petsc_mat);
  //CHKERRQ(ierr);
  ierr = MatSetFromOptions( petsc_mat ); //CHKERRQ(ierr);
#endif

  IntVect lastIV = IntVect::Zero;
  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
      BaseFab<Real>& phi01FAB =  phi01[dit()].getSingleValuedFAB();
#ifdef CH_USE_PETSC
      BaseFab<Real>& columnFAB = column[dit()].getSingleValuedFAB();
      BaseFab<int> & gidsFab = gids[dit()];
#endif

      const Box& box = dbl.get(dit());
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next())
        {
          IntVect iv = bit();
          phi01FAB(lastIV,0) = 0;
          phi01FAB(iv,0) = 1;
          lastIV = iv;
          EBLevelDataOps::setVal(column,0);
          residual(column,phi01,rhs0); // results could be on all procs

#ifdef CH_USE_PETSC
          Real v;
          PetscInt  i,j;
          i = gidsFab(iv,0);
          for (DataIterator dit2 = a_rhs.dataIterator(); dit2.ok(); ++dit2)
            {
              const Box& box2 = dbl.get(dit2());

              for (BoxIterator bit2(box2); bit2.ok(); bit2.next())
                {
                  IntVect iv2 = bit2();
                  v = columnFAB(iv2,0);
                  j = gidsFab(iv2,0);
                  if (abs(v) > 1e-15)
                    {
                      ierr = MatSetValues(petsc_mat,1,&j,1,&i,&v,ADD_VALUES);
                      //CHKERRQ(ierr);
                    }
                }
            }

#endif
        }
    }
#ifdef CH_USE_PETSC
  ierr = MatAssemblyBegin(petsc_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(petsc_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  PetscViewer viewer;
  PetscViewerASCIIOpen( PETSC_COMM_WORLD, "A.m", &viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  MatView( petsc_mat, viewer );
  PetscViewerDestroy(&viewer);
  MPI_Barrier(PETSC_COMM_WORLD);
#endif
}

void EBAMRPoissonOp::
preCond(LevelData<EBCellFAB>&       a_lhs,
        const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIMERS("EBAMRPoissonOp::preCond");
  CH_TIMER("getInvDiagRHS", t1);
  CH_TIMER("relax", t2);

  // Recall that the operator is: alpha*phi + beta*lap(phi)
  // For isotropic-dx Poisson: alpha=0,beta=1 and the
  // diagonal term of this operator is: 4/h/h in 2D, 6/h/h in 3D.
  // Inverse of the diagonal term is our initial multiplier.

  CH_START(t1);
  getInvDiagRHS(a_lhs,a_rhs);
  CH_STOP(t1);
  CH_START(t2);
  relax(a_lhs, a_rhs, m_numPreCondIters);
  CH_STOP(t2);
}

void EBAMRPoissonOp::
getInvDiagRHS(LevelData<EBCellFAB>&       a_lhs,
              const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIMERS("EBAMRPoissonOp::getInvDiagRHS");
  CH_TIMER("regular", t1);
  CH_TIMER("irregular", t2);

  //this function computes: a_lhs = (1/diagonal)*a_rhs
  //use this to initialize the preconditioner

  Real scale = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      scale += -2.0 * m_beta * m_invDx2[idir];
    }
  scale = 1.0 / scale;

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      CH_START(t1);
      Box dblBox( m_eblg.getDBL().get(dit()) );
      EBCellFAB&                   lhs = a_lhs[dit()];
      BaseFab<Real>& lhsFAB = lhs.getSingleValuedFAB();
      const EBCellFAB&             rhs = a_rhs[dit()];
      const BaseFab<Real>& rhsFAB = rhs.getSingleValuedFAB();
      int ncomps = lhs.nComp();
      CH_assert(ncomps == rhs.nComp());

      FORT_GETINVDIAGRHS(CHF_FRA(lhsFAB),
                         CHF_CONST_FRA(rhsFAB),
                         CHF_CONST_REAL(scale),
                         CHF_CONST_INT(ncomps),
                         CHF_BOX(dblBox));
      CH_STOP(t1);

      CH_START(t2);
      const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[dit()];
      const BaseIVFAB<Real>& curBetaWeight =  m_betaDiagWeight[dit()];

      m_invDiagEBStencil[dit]->apply(lhs, rhs, 1., m_alpha, curAlphaWeight, m_beta, curBetaWeight, 1., false);
      CH_STOP(t2);
    }
}

void EBAMRPoissonOp::
applyOp(LevelData<EBCellFAB>&             a_opPhi,
        const LevelData<EBCellFAB>&       a_phi,
        bool                              a_homogeneousPhysBC)
{
  //homogeneous CFBCs because that is all we can do.
  applyOp(a_opPhi, a_phi, NULL, a_homogeneousPhysBC, true);
}
/**/
void EBAMRPoissonOp::
applyOp(LevelData<EBCellFAB>&                    a_opPhi,
        const LevelData<EBCellFAB>&              a_phi,
        const LevelData<EBCellFAB>* const        a_phiCoar,
        const bool&                              a_homogeneousPhysBC,
        const bool&                              a_homogeneousCFBC)
{
  //apply operator without input EB inhomogeneous BC leveldata
  applyOp(a_opPhi, a_phi, a_phiCoar, a_homogeneousPhysBC, a_homogeneousCFBC, NULL);
}
/**/
//get EB flux from leveldata if we have one.  otherwise use m_ebBC
//ignore input data in the case of homogeneous Phys BC
void EBAMRPoissonOp::
applyOp(LevelData<EBCellFAB>&                    a_opPhi,
        const LevelData<EBCellFAB>&              a_phi,
        const LevelData<EBCellFAB>* const        a_phiCoar,
        DataIterator&                            a_dit,
        const bool&                              a_homogeneousPhysBC,
        const bool&                              a_homogeneousCFBC,
        const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD //only non null in multifluid
        )
{
  CH_TIMERS("EBAMRPoissonOp::applyOp");
  CH_TIMER("coarse-fine_bcs", t6);
  CH_assert(a_opPhi.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  CH_assert( (! (a_phiCoar!=NULL && a_phiCoar->isDefined()) )
             ||  (a_phiCoar->ghostVect() == m_ghostCellsPhi) );

  CH_assert(a_phi.nComp() == a_opPhi.nComp());

  LevelData<EBCellFAB>& phi = const_cast<LevelData<EBCellFAB>&>(a_phi);
  if (m_hasCoar && (!s_turnOffBCs))
    {
      CH_START(t6);
      applyCFBCs(phi, a_phiCoar, a_homogeneousCFBC);
      CH_STOP(t6);
    }
  phi.exchange(phi.interval());

  applyOpNoCFBCs( a_opPhi,
                  a_phi,
                  a_phiCoar,
                  a_dit,
                  a_homogeneousPhysBC,
                  a_homogeneousCFBC,
                  a_ebFluxBCLD );
}

void EBAMRPoissonOp::
applyOp(LevelData<EBCellFAB>&                    a_opPhi,
        const LevelData<EBCellFAB>&              a_phi,
        const LevelData<EBCellFAB>* const        a_phiCoar,
        const bool&                              a_homogeneousPhysBC,
        const bool&                              a_homogeneousCFBC,
        const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD //only non null in multifluid
        )
{
  DataIterator dit = a_phi.dataIterator();
  applyOp(a_opPhi,
          a_phi,
          a_phiCoar,
          dit,
          a_homogeneousPhysBC,
          a_homogeneousCFBC,
          a_ebFluxBCLD);
}

void EBAMRPoissonOp::
applyOpNoCFBCs(LevelData<EBCellFAB>&                    a_opPhi,
               const LevelData<EBCellFAB>&              a_phi,
               const LevelData<EBCellFAB>* const        a_phiCoar,
               DataIterator&                            a_dit,
               const bool&                              a_homogeneousPhysBC,
               const bool&                              a_homogeneousCFBC,
               const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD //only non null in multifluid
               )
{
  CH_TIMERS("EBAMRPoissonOp::applyOpNoCFBCs");
  CH_TIMER("eb_bcs_apply", t3);
  CH_TIMER("regular_apply", t1);
  CH_TIMER("irregular_apply", t2);
  CH_TIMER("dom_bcs_apply", t4);
  CH_TIMER("alpha_apply", t5);
  int nComps = a_phi.nComp();
  LevelData<EBCellFAB>& phi = const_cast<LevelData<EBCellFAB>&>(a_phi);
  bool hasNoEBLevelData = (a_ebFluxBCLD==NULL);

  int ibox = 0;
  for (a_dit.reset(); a_dit.ok(); ++a_dit)
    {
      Box dblBox( m_eblg.getDBL().get(a_dit()) );
      const EBCellFAB& curPhiEBCellFAB = phi[a_dit()];
      Box curPhiBox = curPhiEBCellFAB.box();
      const BaseFab<Real>& curPhiFAB = curPhiEBCellFAB.getSingleValuedFAB();

      EBCellFAB& curOpPhiEBCellFAB = a_opPhi[a_dit()];
      BaseFab<Real>& curOpPhiFAB = curOpPhiEBCellFAB.getSingleValuedFAB();

      CH_START(t5);
      if (m_alpha == 0)
        {
          curOpPhiEBCellFAB.setVal(0.0);
        }
      else
        {
          curOpPhiEBCellFAB.copy(curPhiEBCellFAB);
          curOpPhiEBCellFAB.mult(m_alpha);
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
                             a_dit(),
                             m_beta);
      CH_STOP(t1);

      CH_START(t2);
      const BaseIVFAB<Real>& alphaWeight = m_alphaDiagWeight[a_dit()];
      m_opEBStencil[a_dit()]->apply(curOpPhiEBCellFAB, curPhiEBCellFAB, alphaWeight, m_alpha, m_beta, false);
      CH_STOP(t2);

      CH_START(t4);
      const BaseIVFAB<Real>& one = m_one[a_dit];
      Real alpha;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (a_homogeneousPhysBC)
            {
              alpha = 0.;
            }
          else
            {
              alpha = -m_beta*m_invDx[idir];
            }
          m_opEBStencilInhomDomLo[idir][a_dit]->apply(curOpPhiEBCellFAB,
                                                      (*m_cacheInhomDomBCLo[s_whichComp][idir])[a_dit()],
                                                      one,
                                                      alpha,
                                                      1.,
                                                      true);

          if (a_homogeneousPhysBC)
            {
              alpha = 0.;
            }
          else
            {
              alpha = m_beta*m_invDx[idir];
            }
          m_opEBStencilInhomDomHi[idir][a_dit]->apply(curOpPhiEBCellFAB,
                                                      (*m_cacheInhomDomBCHi[s_whichComp][idir])[a_dit()],
                                                      one,
                                                      alpha,
                                                      1.,
                                                      true);
        }
      CH_STOP(t4);

      const Real factor = m_beta/m_dx[0];
      CH_START(t3);
      if (hasNoEBLevelData)
        {
          //standard EB boundary conditions for inhomogeneous, single fluid
          if (!a_homogeneousPhysBC)
            {
              m_ebBC->applyEBFlux(curOpPhiEBCellFAB, curPhiEBCellFAB, m_vofItIrreg[a_dit()], (*m_eblg.getCFIVS()),
                                  a_dit(), m_origin, m_dx, factor,
                                  a_homogeneousPhysBC, s_time);
            }
        }
      else
        {
          //this stuff is for multifluid
          const EBISBox& ebisBox = m_eblg.getEBISL()[a_dit()];
          // Use vofit defined over EBISBox's boundary IVS because this mirrors
          // a_ebFluxBCLD's VoFs
          VoFIterator vofit(ebisBox.boundaryIVS(dblBox),ebisBox.getEBGraph());
          for (vofit.reset(); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const BaseIVFAB<Real>& ebInputFluxFAB = (*a_ebFluxBCLD)[a_dit()];
              Real areaFrac = ebisBox.bndryArea(vof);
              Real ebFlux = -ebInputFluxFAB(vof, 0);
              ebFlux *= areaFrac;
              // do area fraction scaling. (For single fluid, it was
              // incorporated into the Dirichlet EBBC stencil.)
              if (s_areaFracWeighted)
                {
                  ebFlux *= ebisBox.areaFracScaling(vof);
                }
              curOpPhiEBCellFAB(vof,0) += ebFlux * factor;
            }

        }
      CH_STOP(t3);

      ibox++;
    }
}

/***/
void
EBAMRPoissonOp::
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
  CH_TIME("EBAMRPoissonOp::applyOpRegularAllDirs");
  CH_assert(m_domainBC != NULL);

  //need to monkey with the ghost cells to account for boundary conditions
  if (!s_turnOffBCs)
    {
      BaseFab<Real>& phiFAB = (BaseFab<Real>&) a_curPhiFAB;
      applyDomainFlux(a_loBox, a_hiBox, a_hasLo, a_hasHi,
                      a_curDblBox, a_nComps, phiFAB,
                      a_homogeneousPhysBC, a_dit,m_beta);
    }

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
EBAMRPoissonOp::
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
  CH_TIME("EBAMRPoissonOp::applyDomainFlux");
  CH_assert(m_domainBC != NULL);

  for (int idir=0; idir<SpaceDim; idir++)
    {

      EBArith::loHi(a_loBox[idir], a_hasLo[idir],
                    a_hiBox[idir], a_hasHi[idir],
                    m_eblg.getDomain(),a_dblBox, idir);

      for (int comp = 0; comp<a_nComps; comp++)
        {

          if (a_hasLo[idir] == 1)
            {
              Box lbox=a_loBox[idir];
              lbox.shift(idir,-1);
              BaseFab<Real> loFaceFlux(a_loBox[idir],a_nComps);
              int side = -1;

              m_domainBC->getFaceFlux(loFaceFlux,a_phiFAB,m_origin,m_dx,idir,Side::Lo,a_dit,s_time,a_homogeneousPhysBC);

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

              m_domainBC->getFaceFlux(hiFaceFlux,a_phiFAB,m_origin,m_dx,idir,Side::Hi,a_dit,s_time,a_homogeneousPhysBC);

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

void EBAMRPoissonOp::applyOpNoBoundary(LevelData<EBCellFAB>&        a_opPhi,
                                       const LevelData<EBCellFAB>&  a_phi)
{
  s_turnOffBCs = true;
  applyOp(a_opPhi, a_phi, true);
  s_turnOffBCs = false;
}

void EBAMRPoissonOp::
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
}

void EBAMRPoissonOp::
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

void EBAMRPoissonOp::
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  EBLevelDataOps::assign(a_lhs,a_rhs);
}

Real EBAMRPoissonOp::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
 ProblemDomain domain;
 Real volume;

 return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);

  ///warning this will not include kappa and it will include values in covered cells
  //Real sum = 0;
  //for (DataIterator dit = a_1.dataIterator(); dit.ok(); ++dit)
  //  {
  //    FORT_EBAMRPDOTPROD(CHF_REAL(sum),
  //                       CHF_CONST_FRA1(a_1[dit()].getSingleValuedFAB(), 0),
  //                       CHF_CONST_FRA1(a_2[dit()].getSingleValuedFAB(), 0),
  //                       CHF_BOX(m_eblg.getDBL().get(dit())));
  //  }
  //Real volume = 0;
  //EBLevelDataOps::gatherBroadCast(sum, volume, 0);
  //return sum;
}

void EBAMRPoissonOp::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}

void EBAMRPoissonOp::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}

void EBAMRPoissonOp::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  EBLevelDataOps::scale(a_lhs,a_scale);
}

Real EBAMRPoissonOp::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
 CH_TIMERS("EBAMRPoissonOp::norm");
 CH_TIMER("mpi_allreduce",t1);

 Real maxNorm = 0.0;

 maxNorm = localMaxNorm(a_rhs);

 CH_START(t1);
#ifdef CH_MPI
       Real tmp = 1.;
       int result = MPI_Allreduce(&maxNorm, &tmp, 1, MPI_CH_REAL,
                         MPI_MAX, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
         { //bark!!!
           MayDay::Error("sorry, but I had a communcation error on norm");
         }
       maxNorm = tmp;
#endif
//  Real volume=1.;
//  EBLevelDataOps::gatherBroadCast(maxNorm, volume, 0);
 CH_STOP(t1);

 return maxNorm;
}

Real EBAMRPoissonOp::
localMaxNorm(const LevelData<EBCellFAB>& a_rhs)
{
 CH_TIME("EBAMRPoissonOp::localMaxNorm");
 return  EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);
//  ProblemDomain domain;
//  Real volume;
//  return EBLevelDataOps::kappaNorm(volume,a_rhs,EBLEVELDATAOPS_ALLVOFS,domain,0);
}
Real EBAMRPoissonOp::
staticMaxNorm(const LevelData<EBCellFAB>& a_rhs, const EBLevelGrid& a_eblg)
{

 Real maxNorm = 0.0;

 for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
   {
     int iRegIrregCovered;
     const BaseFab<int>& maskFAB = a_eblg.getEBISL()[dit()].getEBGraph().getMask(iRegIrregCovered);
     if (iRegIrregCovered != -1)//not all covered
       {
         if (iRegIrregCovered == 0)//has irreg
           {
             const Box& box = a_eblg.getDBL().get(dit());
             //             const EBISBox& ebisBox = a_eblg.getEBISL()[dit()];
             const BaseFab<Real>& rhsFAB = (a_rhs[dit()]).getSingleValuedFAB();
             FORT_MAXNORMMASK(CHF_REAL(maxNorm),
                              CHF_CONST_FRA1(rhsFAB,0),
                              CHF_BOX(box),
                              CHF_CONST_FRA1(maskFAB,0));

             //CP: this portion between the stars is the new faster code:
             //****************************
             int srccomp = 0;
             int ncomp   = 1;
             const BaseIVFAB<Real>& irrBFAB1 = a_rhs[dit()].getMultiValuedFAB();
             const Real* r = irrBFAB1.dataPtr(srccomp);
             int nvof    = irrBFAB1.numVoFs();
             for (int i=0; i<nvof*ncomp; i++)
               {
                 maxNorm = Max(maxNorm, Abs(r[i]));
               }
             //*****************************

             //CP: below and above else are the old codes that were slow
             // IntVectSet ivs = ebisBox.getMultiCells(box);
             // for (VoFIterator vofit(ivs, a_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
             //   {
             //     if (Abs(a_rhs[dit()](vofit(), 0)) > maxNorm) {maxNorm = Abs(a_rhs[dit()](vofit(), 0));}
             //   }
           }
         else//all reg
           {
             const Box& box = a_eblg.getDBL().get(dit());
             const BaseFab<Real>& rhsFAB = (a_rhs[dit()]).getSingleValuedFAB();
             FORT_MAXNORM(CHF_REAL(maxNorm),
                          CHF_CONST_FRA1(rhsFAB,0),
                          CHF_BOX(box));
           }
       }
   }

 return maxNorm;
}

void EBAMRPoissonOp::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  EBLevelDataOps::setToZero(a_lhs);
}

void EBAMRPoissonOp::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void  EBAMRPoissonOp::
setEBBC(const RefCountedPtr<BaseEBBC>&      a_ebBC)
{
  m_ebBC     = a_ebBC;
}

// MGLevelOp functions

void EBAMRPoissonOp::
createCoarser(LevelData<EBCellFAB>&       a_coar,
              const LevelData<EBCellFAB>& a_fine,
              bool                        a_ghosted)
{
  CH_assert(a_fine.nComp() == 1);
  const DisjointBoxLayout& dbl = m_eblgCoarMG.getDBL();
  EBISLayout coarEBISL = m_eblgCoarMG.getEBISL();
  /*
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), 2);

  int nghost = a_fine.ghostVect()[0];
  EBISLayout coarEBISL;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(coarEBISL,
                          dbl, coarDom, nghost);
  */
  EBCellFactory ebcellfact(coarEBISL);
  a_coar.define(dbl, 1,a_fine.ghostVect(),ebcellfact);
}

void EBAMRPoissonOp::
relax(LevelData<EBCellFAB>&       a_e,
      const LevelData<EBCellFAB>& a_residual,
      int                         a_iterations)
{
  CH_TIME("EBAMRPoissonOp::relax");

  CH_assert(a_e.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_e.nComp() == 1);
  CH_assert(a_residual.nComp() == 1);

  if (m_relaxType == 1)
    {
      for (int i = 0; i < a_iterations; i++)
        {
          levelMultiColorGS(a_e,a_residual);
        }
    }
  else  if (m_relaxType == 2)
    {
      for (int i = 0; i < a_iterations; i++)
        {
          levelGSRB(a_e,a_residual);
        }
    }
  else  if (m_relaxType == 3)
    {
      for (int i = 0; i < a_iterations; i++)
        {
          levelSlowRelax(a_e,a_residual);
        }
    }
  else if (m_relaxType == 4)
    {
      for (int i = 0; i < a_iterations; i++)
        {
          levelMultiColorGSClone(a_e,a_residual);
        }
    }
  else if (m_relaxType == 999)
    {
      //CP added, no relax. using this option to go directly into bottom solve
      //a_residual.copyTo(a_e); // This is wrong
      //do nothing
    }
  else
    {
      MayDay::Error("EBAMRPoissonOp::relax - invalid relaxation type");
    }
}
/****/
void EBAMRPoissonOp::
restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                 LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_rhsThisLevel)
{
  CH_TIME("EBAMRPoissonOp::restrictResidual");

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
void EBAMRPoissonOp::
prolongIncrement(LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_correctCoar)
{
  CH_TIME("EBAMRPoissonOp::prolongIncrement");
  Interval vars(0, 0);
  m_ebInterpMG.pwcInterp(a_phiThisLevel, a_correctCoar, vars);
}
//////////
void
EBAMRPoissonOp::
applyCFBCs(LevelData<EBCellFAB>&             a_phi,
           const LevelData<EBCellFAB>* const a_phiCoar,
           bool a_homogeneousCFBC,
           bool a_doOnlyRegularInterp) //=false by default
{
  CH_TIMERS("EBAMRPoissonOp::applyCFBCs");
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
          m_quadCFIWithCoar->interpolate(a_phi, *a_phiCoar, interv, a_doOnlyRegularInterp);
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
////
void
EBAMRPoissonOp::
applyHomogeneousCFBCs(LevelData<EBCellFAB>&   a_phi)
{
  CH_TIME("EBAMRPoissonOp::applyHomogeneousCFBCs");
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
/////////
void
EBAMRPoissonOp::
applyHomogeneousCFBCs(EBCellFAB&            a_phi,
                      const DataIndex&      a_datInd,
                      int                   a_idir,
                      Side::LoHiSide        a_hiorlo)
{
  CH_TIMERS("EBAMRPoissonOp::applyHomogeneousCFBCs2");
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
                      CHF_CONST_REAL(m_dx[a_idir]),
                      CHF_CONST_REAL(m_dxCoar[a_idir]),
                      CHF_CONST_INT(a_idir),
                      CHF_CONST_INT(ihiorlo));

      CH_STOP(t1);
    }
  else
    {
      if (!interpIVS.isEmpty())
        {
          CH_START(t2);
          Real halfdxcoar = m_dxCoar[a_idir]/2.0;
          Real halfdxfine = m_dx[a_idir]/2.0;
          Real xg = halfdxcoar -   halfdxfine;
          Real xc = halfdxcoar +   halfdxfine;
          Real xf = halfdxcoar + 3*halfdxfine;
          Real hf = m_dx[a_idir];
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
////
int EBAMRPoissonOp::
refToCoarser()
{
  return m_refToCoar;
}

////
int EBAMRPoissonOp::
refToFiner()
{
  return m_refToFine;
}

////
void EBAMRPoissonOp::
AMRResidual(LevelData<EBCellFAB>&       a_residual,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            const LevelData<EBCellFAB>& a_rhs,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("EBAMRPoissonOp::AMRResidual");
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


////
void EBAMRPoissonOp::
AMROperator(LevelData<EBCellFAB>&       a_LofPhi,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("EBAMRPoissonOp::AMROperator");
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



void EBAMRPoissonOp::
reflux(LevelData<EBCellFAB>& a_residual,
       const LevelData<EBCellFAB>& a_phiFine,
       const LevelData<EBCellFAB>& a_phi,
       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("EBAMRPoissonOp::reflux");
  this->fast_reflux(a_residual, a_phiFine, a_phi, a_finerOp);
}


void EBAMRPoissonOp::
fast_reflux(LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("EBAMRPoissonOp::fastReflux");
  CH_TIMER("define",t1);
  CH_TIMER("setToZero",t2);
  CH_TIMER("incrementCoar",t3);
  CH_TIMER("incrementFine",t4);
  CH_TIMER("reflux_from_reg",t5);
  Interval interv(0,0);
  CH_START(t1);
  CH_assert(a_phiFine.nComp() == 1);

  CH_STOP(t1);
  CH_START(t2);
  m_fastFR.setToZero();
  CH_STOP(t2);
  CH_START(t3);
  fast_incrementFRCoar(a_phiFine, a_phi);
  CH_STOP(t3);

  CH_START(t4);
  fast_incrementFRFine(a_phiFine, a_phi, a_finerOp);
  CH_STOP(t4);
  CH_START(t5);

  Real scale = 1.0/m_dx[0];
  m_fastFR.reflux(a_residual, interv, scale);

  CH_STOP(t5);
}

void EBAMRPoissonOp::
fast_incrementFRCoar(const LevelData<EBCellFAB>& a_phiFine,
                     const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBAMRPoissonOp::incrementFRCoar");
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
          Box ghostedBox = box;
          ghostedBox.grow(1);
          ghostedBox.grow(idir,-1);
          ghostedBox &= m_eblg.getDomain();

          EBFaceFAB coarflux(ebisBox, ghostedBox, idir, ncomp);
          //old way
          //getFlux(coarflux, coarfab, ghostedBox, box, m_eblg.getDomain(), ebisBox, m_dx, idir);

          //new way
          getFluxRegO(coarflux, coarfab, ghostedBox, m_dx);
          Real scale = 1.0;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex>*  faceit;
              Vector<VoFStencil>* stencil;
              int index = EBFastFR::index(idir, sit());
              if (m_hasEBCF)
                {
                  faceit  = &( m_faceitCoar[index][dit()]);
                  stencil = &(m_stencilCoar[index][dit()]);
                }
              getFluxEBCF(coarflux, coarfab, ghostedBox, *faceit, *stencil, m_dx);
            }
          for (SideIterator sit; sit.ok(); ++sit)
            {
              //both here means regular and irregular.   the name is a
              //relic of a former age
              m_fastFR.incrementCoarseBoth(coarflux, scale, dit(), interv, idir, sit());
            }
        }
    }
}

/***/
void EBAMRPoissonOp::
fast_incrementFRFine(const LevelData<EBCellFAB>& a_phiFine,
                     const LevelData<EBCellFAB>& a_phi,
                     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("EBAMRPoissonOp::incrementFRFine");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(m_hasFine);
  int ncomp = 1;
  Interval interv(0,0);
  EBAMRPoissonOp& finerEBAMROp = (EBAMRPoissonOp& )(*a_finerOp);

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
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Box fabBox = adjCellBox(boxFine, idir, sit(), 1);
              fabBox.shift(idir, -sign(sit()));

              Box ghostedBox = fabBox;
              ghostedBox.grow(1);
              ghostedBox.grow(idir,-1);
              ghostedBox &= m_eblgFine.getDomain();

              EBFaceFAB fluxFine(ebisBoxFine, ghostedBox, idir, ncomp);

              getFlux(fluxFine, phiFine, ghostedBox, fabBox, m_eblgFine.getDomain(), ebisBoxFine, m_dxFine, idir);

              Real scale = 1.0;
              //both here means regular and irregular.   the name is a
              //relic of a former age
              m_fastFR.incrementFineBoth(fluxFine, scale, ditf(), interv,  idir,    sit());
            }
        }
    }
}

/****/
void EBAMRPoissonOp::
getFlux(EBFaceFAB&                    a_fluxCentroid,
        const EBCellFAB&              a_phi,
        const Box&                    a_ghostedBox,
        const Box&                    a_fabBox,
        const ProblemDomain&          a_domain,
        const EBISBox&                a_ebisBox,
        const RealVect&               a_dx,
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

  {
    CH_TIME("regular stuff");
    FORT_REGGETFLUX(CHF_FRA(regFlux),
                    CHF_CONST_FRA(regPhi),
                    CHF_BOX(faceBox),
                    CHF_CONST_REAL(m_beta),
                    CHF_CONST_REALVECT(a_dx),
                    CHF_CONST_INT(a_idir),
                    CHF_CONST_INT(ncomp));
  }

  {
    CH_TIME("irregular stuff");
    EBFaceFAB fluxCenter(a_ebisBox, a_ghostedBox, a_idir,1);
    fluxCenter.copy(a_fluxCentroid);

    IntVectSet ivsCell = a_ebisBox.getIrregIVS(cellBox);
    if (!ivsCell.isEmpty())
      {
        FaceStop::WhichFaces stopCrit;
        if (m_eblg.getDomain().isPeriodic(a_idir))
          {
            stopCrit = FaceStop::SurroundingWithBoundary;
          }
        else
          {
            stopCrit = FaceStop::SurroundingNoBoundary;
          }
        for (FaceIterator faceit(ivsCell, a_ebisBox.getEBGraph(), a_idir,stopCrit);
            faceit.ok(); ++faceit)
          {
            const FaceIndex& face = faceit();
            Real phiHi = a_phi(face.getVoF(Side::Hi), 0);
            Real phiLo = a_phi(face.getVoF(Side::Lo), 0);
            Real fluxFace = m_beta*(phiHi - phiLo)/a_dx[a_idir];

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
  }
}

/****/
void EBAMRPoissonOp::
getFluxEBCF(EBFaceFAB&                    a_flux,
            const EBCellFAB&              a_phi,
            const Box&                    a_ghostedBox,
            Vector<FaceIndex>&            a_faceitEBCF,
            Vector<VoFStencil>&           a_stenEBCF,
            const RealVect&               a_dx)
{
  CH_TIME("EBAMRPoissonOp::getFluxEBCF");

  //only do the evil stuff if you have a coarse-fine /  EB crossing situation

  if (m_hasEBCF)
    {
      CH_TIME("EBCF stuff");
      for (int iface = 0; iface < a_faceitEBCF.size(); iface++)
        {
          const FaceIndex& face =     a_faceitEBCF[iface];
          const VoFStencil& stencil   = a_stenEBCF[iface];
          Real fluxval = 0;
          for (int isten = 0; isten < stencil.size(); isten++)
            {
              fluxval += stencil.weight(isten)*(a_phi(stencil.vof(isten), 0));
            }
          //note the last minute beta
          a_flux(face, 0) = m_beta*fluxval;
        }
    }
}

/****/
void EBAMRPoissonOp::
getFluxRegO(EBFaceFAB&                    a_flux,
            const EBCellFAB&              a_phi,
            const Box&                    a_ghostedBox,
            const RealVect&               a_dx)
{
  CH_TIME("EBAMRPoissonOp::getFluxEBCF");
  //has some extra cells so...
  a_flux.setVal(0.);
  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_flux.nComp());
  Box cellBox = a_ghostedBox;
  int idir = a_flux.direction();
  //want only interior faces
  cellBox.grow(idir, 1);
  cellBox &= m_eblg.getDomain();
  cellBox.grow(idir,-1);

  Box faceBox = surroundingNodes(cellBox, idir);

  //make a EBFaceFAB (including ghost cells) that will hold centered gradients
  BaseFab<Real>& regFlux  = a_flux.getSingleValuedFAB();
  const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();

  {
    CH_TIME("regular stuff");
    FORT_REGGETFLUX(CHF_FRA(regFlux),
                    CHF_CONST_FRA(regPhi),
                    CHF_BOX(faceBox),
                    CHF_CONST_REAL(m_beta),
                    CHF_CONST_REALVECT(a_dx),
                    CHF_CONST_INT(idir),
                    CHF_CONST_INT(ncomp));
  }
}

void EBAMRPoissonOp::dumpStencilMatrix()
{
  int phase = 0;
  map<StencilIndex, Real, StencilIndexComparator> mapper;

  // Get the Level and EB data
  const IntVect& domainSize = m_eblg.getDomain().size();

  // vector that will hold the volume fraction
  Vector<Real> kappa(0);

  // loop through VoFs adding non-covered ones to map, along
  //  with their volume fractions
  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box box = dbl[dit()];
      IntVectSet ivsall(box);
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const EBGraph& ebGraph = ebisBox.getEBGraph();
      for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
        {
          VolIndex srcVof(vofit());
          IntVect srciv = srcVof.gridIndex();
          if (!ebisBox.isCovered(srciv))
            {
              StencilIndex si(srcVof, phase);
              kappa.push_back(ebisBox.volFrac(srcVof));
              mapper[si] = ebisBox.areaFracScaling(srcVof);
            }
        }
    }

  // all the VoFs have been added to the map, and it is properly sorted.
  //  so, we just iterate through it setting the map_values sequentially.
  char charstr[100];
  map<StencilIndex, Real>::iterator mit;
  int imap = 0;
  for (mit = mapper.begin(); mit != mapper.end(); ++mit)
    {
      Real k = kappa[imap];
      Real alphaMax = mit->second;
      IntVect iv = mit->first.vof().gridIndex();
      mit->second = ++imap;
#if CH_SPACEDIM==2
      sprintf(charstr, "M_%03d(%03d,:) = [ %4d %4d %22.16e %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1],
              alphaMax, k);
#elif CH_SPACEDIM==3
      sprintf(charstr, "M_%03d(%03d,:) = [ %4d %4d %4d %22.16e %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1], iv[2],
              alphaMax, k);
#endif
      string str(charstr);
      pout() << str << endl;
    }

  int comp = 0;
  LayoutData< BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(comp);

  // loop through vofs outputting stencil weights
  for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box box = dbl[dit()];
      IntVectSet ivsall(box);
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const EBGraph& ebGraph = ebisBox.getEBGraph();
      BaseIVFAB<VoFStencil>& stenfab = (*fluxStencil)[dit()];
      const BaseIVFAB<Real>& alphaWeight = m_alphaDiagWeight[dit()];
      for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
        {
          VolIndex srcVof(vofit());
          IntVect srciv = srcVof.gridIndex();
          if (!ebisBox.isCovered(srciv))
            {
              StencilIndex si(srcVof, phase);
              map<StencilIndex, Real>::iterator srcit = mapper.find(si);
              if (srcit == mapper.end())
                {
                  pout() << "could not find " << srcVof << endl;
                  MayDay::Error("Source VoF not found in map");
                }
              else
                {
                  // mapping for srcVof
                  int is = (int) srcit->second;
                  // stencil for faces
                  VoFStencil sten = VoFStencil();
                  getDivFStencil(sten, srcVof, dit(), true);
                  {
                    VoFStencil domainFluxStencil = VoFStencil();
                    getDomainFluxStencil(domainFluxStencil,
                                         srcVof,
                                         comp,
                                         dit());
                    sten += domainFluxStencil;
                  }
                  if (ebisBox.isIrregular(srciv))
                    {
                      // stencil for EB
                      VoFStencil boundarySten = stenfab(srcVof, comp);
                      // add EB stencil to face stencil
                      sten += boundarySten;
                    }
                  // scale stencil by beta
                  sten *= m_beta;
                  // add alpha
                  if (m_alpha != 0)
                    {
                      if (ebisBox.isIrregular(srciv))
                        {
                          sten.add(srcVof, m_alpha*alphaWeight(srcVof, comp), comp);
                        }
                      else
                        {
                          sten.add(srcVof, m_alpha, comp);
                        }
                    }
                  for (int i=0; i<sten.size(); ++i)
                    {
                      VolIndex dstVof(sten.vof(i));
                      StencilIndex di(dstVof, phase);
                      map<StencilIndex, Real>::iterator dstit = mapper.find(di);
                      if (dstit == mapper.end())
                        {
                          pout() << "could not find " << dstVof << endl;
                          MayDay::Error("Destination VoF not found in map");
                        }
                      else if (sten.weight(i) != 0.)
                        {
                          // mapping for dstVof
                          int id = (int) dstit->second;
                          sprintf(charstr,"L_%03d(%03d, %03d) = %22.16e;",
                                  domainSize[0], is, id, sten.weight(i));
                          string str(charstr);
                          pout() << str << endl;
                        }
                    }
                }
            }
        }
    }
}

void EBAMRPoissonOp::getDomainFluxStencil(      VoFStencil& a_stencil,
                                          const VolIndex&   a_vof,
                                          const int         a_comp,
                                          const DataIndex&  a_dix)
{
  a_stencil.clear();
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_dix];
  for (int idir=0; idir<SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());
          for (int iface=0; iface<faces.size(); iface++)
            {
              if (faces[iface].isBoundary())
                {
                  VoFStencil thisStencil = VoFStencil();
                  m_domainBC->getFluxStencil(thisStencil,
                                             faces[iface],
                                             a_comp,
                                             m_dx,
                                             idir,
                                             sit(),
                                             ebisBox);
                  Real factor = m_invDx[idir];
                  if (s_areaFracWeighted)
                    {
                      factor *= ebisBox.areaFracScaling(a_vof);
                    }
                  thisStencil *= factor;
                  a_stencil += thisStencil;
                }
            }
        }
    }
}

void EBAMRPoissonOp::dumpReferenceStencilMatrix()
{
  int phase = 0;
  map<StencilIndex, Real, StencilIndexComparator> mapper;

  // Get the Level and EB data
  const IntVect& domainSize = m_eblg.getDomain().size();

  // vector that will hold the volume fraction
  Vector<Real> kappa(0);

  // loop through VoFs adding non-covered ones to map
  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box box = dbl[dit()];
      IntVectSet ivsall(box);
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const EBGraph& ebGraph = ebisBox.getEBGraph();
      for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
        {
          VolIndex srcVof(vofit());
          IntVect srciv = srcVof.gridIndex();
          if (!ebisBox.isCovered(srciv))
            {
              StencilIndex si(srcVof, phase);
              kappa.push_back(ebisBox.volFrac(srcVof));
              mapper[si] = ebisBox.areaFracScaling(srcVof);
            }
        }
    }

  // all the VoFs have been added to the map, and it is properly sorted.
  //  so, we just iterate through it setting the map_values sequentially.
  char charstr[100];
  map<StencilIndex, Real>::iterator mit;
  int imap = 0;
  for (mit = mapper.begin(); mit != mapper.end(); ++mit)
    {
      Real k = kappa[imap];
      Real alphaMax = mit->second;
      IntVect iv = mit->first.vof().gridIndex();
      mit->second = ++imap;
#if CH_SPACEDIM==2
      sprintf(charstr, "Mref_%03d(%03d,:) = [ %4d %4d %22.16e %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1],
              alphaMax, k);
#elif CH_SPACEDIM==3
      sprintf(charstr, "Mref_%03d(%03d,:) = [ %4d %4d %4d %22.16e %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1], iv[2],
              alphaMax, k);
#endif
      string str(charstr);
      pout() << str << endl;
    }

  // loop through vofs outputting stencil weights
  int nvar = 1;
  EBCellFactory factory(m_eblg.getEBISL());
  LevelData<EBCellFAB> phi(dbl, nvar, m_ghostCellsPhi, factory);
  LevelData<EBCellFAB> lhs(dbl, nvar, m_ghostCellsRHS, factory);

  for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box box = dbl[dit()];
      IntVectSet ivsall(box);
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const EBGraph& ebGraph = ebisBox.getEBGraph();
      EBCellFAB& data = phi[dit()];
      for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
        {
          VolIndex srcVof(vofit());
          IntVect srciv = srcVof.gridIndex();
          if (!ebisBox.isCovered(srciv))
            {
              StencilIndex si(srcVof, phase);
              map<StencilIndex, Real>::iterator srcit = mapper.find(si);
              if (srcit == mapper.end())
                {
                  pout() << "could not find " << srcVof << endl;
                  MayDay::Error("Source VoF not found in map");
                }
              else
                {
                  // mapping for srcVof
                  int is = (int) srcit->second;
                  // calculate weights for srcVof
                  setToZero(lhs);
                  setToZero(phi);
                  data(srcVof, 0) = 1;
                  applyOp(lhs, phi, true);
                  // loop through vofs looking for non-zero entries
                  for (DataIterator ditL=dbl.dataIterator(); ditL.ok(); ++ditL)
                    {
                      Box boxL = dbl[ditL()];
                      IntVectSet ivsL(boxL);
                      const EBISBox& ebisBoxL = m_eblg.getEBISL()[ditL()];
                      const EBGraph& ebGraphL = ebisBoxL.getEBGraph();
                      EBCellFAB& dataL = lhs[ditL()];
                      for (VoFIterator vofitL(ivsL, ebGraphL); vofitL.ok(); ++vofitL)
                        {
                          VolIndex lVof(vofitL());
                          IntVect liv = lVof.gridIndex();
                          if ( (!ebisBox.isCovered(liv)) &&
                               (dataL(lVof, 0) != 0.) )
                            {
                              StencilIndex li(lVof, phase);
                              map<StencilIndex, Real>::iterator lit = mapper.find(li);
                              if (lit == mapper.end())
                                {
                                  pout() << "could not find " << lVof << endl;
                                  MayDay::Error("Destination VoF not found in map");
                                }
                              else
                                {
                                  // mapping for lVof
                                  int il = (int) lit->second;
                                  sprintf(charstr,"Lref_%03d(%03d, %03d) = %22.16e;",
                                          domainSize[0], il, is, dataL(lVof, 0));
                                  string str(charstr);
                                  pout() << str << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/****/
void EBAMRPoissonOp::
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
void EBAMRPoissonOp::
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
void EBAMRPoissonOp::
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
//
void EBAMRPoissonOp::
AMROperatorNF(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoar,
              bool a_homogeneousPhysBC)
{
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);

  applyOp(a_LofPhi,a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);

}



/****/
void EBAMRPoissonOp::
AMRRestrict(LevelData<EBCellFAB>&       a_resCoar,
            const LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_correction,
            const LevelData<EBCellFAB>& a_coarCorrection, 
            bool a_skip_res )
{
  CH_TIME("EBAMRPoissonOp::AMRRestrict");
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_correction.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostCellsPhi);
  CH_assert(!a_skip_res);

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
///compute norm over all cells on coarse not covered by finer
Real
EBAMRPoissonOp::AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
                        const LevelData<EBCellFAB>& a_fineResid,
                        const int& a_refRat,
                        const int& a_ord)

{
  CH_TIME("EBAMRPoissonOp::AMRNorm");
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
///////////
void EBAMRPoissonOp::
AMRProlong(LevelData<EBCellFAB>&       a_correction,
           const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("EBAMRPoissonOp::AMRProlong");
  //use cached interpolation object
  Interval variables(0, 0);
  CH_assert(m_hasInterpAve);
  m_ebInterp.pwcInterp(a_correction, a_coarCorrection, variables);
}
///////////
void EBAMRPoissonOp::
AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                  const LevelData<EBCellFAB>& a_correction,
                  const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("EBAMRPoissonOp::AMRUpdateResidual");
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

///////////
void EBAMRPoissonOp::
AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                  const LevelData<EBCellFAB>& a_correction,
                  const LevelData<EBCellFAB>& a_coarCorrection,
                  const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD)
{
  CH_TIME("EBAMRPoissonOp::AMRUpdateResidual");
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_correction.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostCellsPhi);

  LevelData<EBCellFAB> lcorr;
  bool homogeneousPhys = true;
  bool homogeneousCF   = false;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  lcorr.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

  applyOp(lcorr, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF, a_ebFluxBCLD);

  incr(a_residual, lcorr, -1);
}

void EBAMRPoissonOp::
levelMultiColorGS(LevelData<EBCellFAB>&       a_phi,
                  const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBAMRPoissonOp::levelMultiColorGS");

  //this is a multigrid operator so only homogeneous CF BC and null coar level
  CH_assert(a_rhs.ghostVect()    == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect()    == m_ghostCellsPhi);

  int nComps = a_phi.nComp();
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox(m_eblg.getDBL().get(dit()));
      BaseFab<Real>& phiFAB       = (a_phi[dit()] ).getSingleValuedFAB();

      Box loBox[SpaceDim],hiBox[SpaceDim];
      int hasLo[SpaceDim],hasHi[SpaceDim];

      applyDomainFlux(loBox, hiBox, hasLo, hasHi,
                      dblBox, nComps, phiFAB,
                      true, dit(), m_beta);

    }

  // Loop over all possibilities (in all dimensions)
  for (int icolor=0; icolor < m_colors.size(); ++icolor)
    {
      if (m_hasCoar)
        {
          applyCFBCs(a_phi, NULL, true);
        }

      //when doLazyRelax==true, only relax if a_icolor==0
      //when doLazyRelax==false, relax every color
      if ((!s_doLazyRelax) || (icolor == 0))
        {
          a_phi.exchange(m_exchangeCopier);
        }

      colorGS(a_phi, a_rhs, icolor);
    }
}
void EBAMRPoissonOp::
colorGS(LevelData<EBCellFAB>&       a_phi,
        const LevelData<EBCellFAB>& a_rhs,
        const int&                  a_icolor)
{
  CH_TIME("EBAMRPoissonOp::colorGS");

  bool homogeneous = true;

  Real weight = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0 * m_beta * m_invDx2[idir];
    }
  weight = 1.0 / weight;

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& phifab = a_phi[dit()];
      const EBCellFAB& rhsfab = a_rhs[dit()];
      BaseFab<Real>& phiBaseFAB = (a_phi[dit()]).getSingleValuedFAB();
      const BaseFab<Real>& rhsBaseFAB = (a_rhs[dit()] ).getSingleValuedFAB();

      m_colorEBStencil[a_icolor][dit()]->cachePhi(phifab);

      GSColorAllRegular(phiBaseFAB, rhsBaseFAB, a_icolor, weight, homogeneous, dit());

      m_colorEBStencil[a_icolor][dit()]->uncachePhi(phifab);

      GSColorAllIrregular(phifab, rhsfab, a_icolor, homogeneous, dit());
    }
}

void EBAMRPoissonOp::
GSColorAllRegular(BaseFab<Real>&               a_phi,
                  const BaseFab<Real>&         a_rhs,
                  const int&                   a_icolor,
                  const Real&                  a_weight,
                  const bool&                  a_homogeneousPhysBC,
                  const DataIndex&             a_dit)
{
  CH_TIME("EBAMRPoissonOp::GSColorAllRegular");

  Box dblBox(m_eblg.getDBL().get(a_dit));

  IntVect loIV = dblBox.smallEnd();
  IntVect hiIV = dblBox.bigEnd();

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (loIV[idir] % 2 != m_colors[a_icolor][idir])
        {
          loIV[idir]++;
        }
    }

  if (loIV <= hiIV)
    {
      Box coloredBox(loIV, hiIV);

      for (int comp=0; comp<a_phi.nComp(); comp++)
        {
          FORT_DOALLREGULARMULTICOLOR(CHF_FRA1(a_phi,comp),
                                      CHF_CONST_FRA1(a_rhs,comp),
                                      CHF_CONST_REAL(a_weight),
                                      CHF_CONST_REAL(m_alpha),
                                      CHF_CONST_REAL(m_beta),
                                      CHF_CONST_REALVECT(m_dx),
                                      CHF_BOX(coloredBox));
        }
    }
}
void EBAMRPoissonOp::
GSColorAllIrregular(EBCellFAB&                   a_phi,
                    const EBCellFAB&             a_rhs,
                    const int&                   a_icolor,
                    const bool&                  a_homogeneousPhysBC,
                    const DataIndex&             a_dit)
{
  CH_TIMERS("EBAMRPoissonOp::GSColorAllIrregular");
  CH_TIMER("assignAlphaBetaWeights", t1);
  CH_TIMER("applyColorStencil", t2);

  if (m_vofItIrregColor[a_icolor][a_dit].size() != 0)
    {
      CH_START(t1);
      const BaseIVFAB<Real>& curAlphaWeight  = m_alphaDiagWeight[a_dit];
      const BaseIVFAB<Real>& curBetaWeight   = m_betaDiagWeight[a_dit];
      CH_STOP(t1);

      CH_START(t2);
      //phi = (I-lambda*L)phiOld
      Real safety = 1.0;
      m_colorEBStencil[a_icolor][a_dit]->relax(a_phi, a_rhs, curAlphaWeight, curBetaWeight, m_alpha, m_beta, safety);
      CH_STOP(t2);

    } //vofitIrregColor.size() != 0
}

void EBAMRPoissonOp::
levelGSRB(LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBAMRPoissonOp::levelGSRB");

  bool homogeneous = true;

  a_phi.exchange(m_exchangeCopier);
  //this is a multigrid operator so only homogeneous CF BC and null coar level
  CH_assert(a_rhs.ghostVect()    == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect()    == m_ghostCellsPhi);

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  int nComps = a_phi.nComp();
  int ibox = 0;
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox(m_eblg.getDBL().get(dit()));
      BaseFab<Real>& phiFAB       = (a_phi[dit()] ).getSingleValuedFAB();

      Box loBox[SpaceDim],hiBox[SpaceDim];
      int hasLo[SpaceDim],hasHi[SpaceDim];

      {
        CH_TIME("EBAMRPoissonOp::levelGSRB::applyDomainFlux");

        applyDomainFlux(loBox, hiBox, hasLo, hasHi,
                        dblBox, nComps, phiFAB,
                        homogeneous, dit(),m_beta);
      }
    }

  // do first red, then black passes
  for (int redBlack =0; redBlack <= 1; redBlack++)
    {
      CH_TIME("EBAMRPoissonOp::levelGSRB::Compute");

      //when doLazyRelax==true, only relax if a_icolor==0
      //when doLazyRelax==false, relax every color
      if ((!s_doLazyRelax) || (redBlack == 0))
        {
          CH_TIME("EBAMRPoissonOp::levelGSRB::ExchangeGSRB");
          a_phi.exchange(m_exchangeCopier);
        }

      if (m_hasCoar)
        {
        CH_TIME("EBAMRPoissonOp::levelGSRB::homogeneousCFInterp");
          applyCFBCs(a_phi, NULL, true);
        }

      Real weight = m_alpha;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          weight += -2.0 * m_beta * m_invDx2[idir];
        }
      weight = 1.0 / weight;

      for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB& phifab = a_phi[dit()];
          const EBCellFAB& rhsfab = a_rhs[dit()];

          //cache phi
          for (int c = 0; c < m_colors.size()/2; ++c)
            {
              m_colorEBStencil[m_colors.size()/2*redBlack+c][dit()]->cachePhi(phifab);
            }

          //reg cells
          const Box& region = dbl.get(dit());
          Box dblBox(m_eblg.getDBL().get(dit()));
          BaseFab<Real>& phiBaseFAB       = (a_phi[dit()] ).getSingleValuedFAB();
          const BaseFab<Real>& rhsBaseFAB = (a_rhs[dit()] ).getSingleValuedFAB();

          for (int comp = 0; comp < a_phi.nComp(); comp++)
            {
              FORT_DOALLREGULARGSRB(CHF_FRA1(phiBaseFAB,comp),
                                    CHF_CONST_FRA1(rhsBaseFAB,comp),
                                    CHF_CONST_REAL(weight),
                                    CHF_CONST_REAL(m_alpha),
                                    CHF_CONST_REAL(m_beta),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_BOX(region),
                                    CHF_CONST_INT(redBlack));
            }

          //uncache phi
          for (int c = 0; c < m_colors.size()/2; ++c)
            {
              m_colorEBStencil[m_colors.size()/2*redBlack+c][dit()]->uncachePhi(phifab);
            }

          for (int c = 0; c < m_colors.size()/2; ++c)
            {
              GSColorAllIrregular(phifab, rhsfab, m_colors.size()/2*redBlack+c, homogeneous, dit());
            }
          ibox++;
        }
    }
}

void EBAMRPoissonOp::
levelMultiColorGSClone(LevelData<EBCellFAB>&       a_phi,
                       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBAMRPoissonOp::levelMultiColorGS");

  //this is a multigrid operator so only homogeneous CF BC and null coar level
  CH_assert(a_rhs.ghostVect()    == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect()    == m_ghostCellsPhi);

  LevelData<EBCellFAB> phiOld;
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  phiOld.define( m_eblg.getDBL(), a_phi.nComp(),a_phi.ghostVect(), ebcellfact);

  // Loop over all possibilities (in all dimensions)
  for (int icolor=0; icolor < m_colors.size(); ++icolor)
    {
      if (m_hasCoar)
        {
          applyCFBCs(a_phi, NULL, true);
        }

      //when doLazyRelax==true, only relax if a_icolor==0
      //when doLazyRelax==false, relax every color
      if ((!s_doLazyRelax) || (icolor == 0))
        {
          a_phi.exchange(m_exchangeCopier);
        }

      EBLevelDataOps::clone(phiOld,a_phi);
      colorGSClone(a_phi, phiOld, a_rhs, icolor);
    }
}
void EBAMRPoissonOp::
colorGSClone(LevelData<EBCellFAB>&       a_phi,
             const LevelData<EBCellFAB>& a_phiOld,
             const LevelData<EBCellFAB>& a_rhs,
             const int&                  a_icolor)
{
  CH_TIME("EBAMRPoissonOp::colorGS");

  bool homogeneous = true;

  Real weight = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0 * m_beta * m_invDx2[idir];
    }
  weight = 1.0 / weight;

  GSColorAllRegularClone(a_phi, a_phiOld, a_rhs, a_icolor, weight, homogeneous);

  GSColorAllIrregularClone(a_phi, a_phiOld, a_rhs, a_icolor, homogeneous);

}
void EBAMRPoissonOp::
GSColorAllRegularClone(LevelData<EBCellFAB>&       a_phi,
                       const LevelData<EBCellFAB>& a_phiOld,
                       const LevelData<EBCellFAB>& a_rhs,
                       const int&                  a_icolor,
                       const Real&                 a_weight,
                       const bool&                 a_homogeneousPhysBC)
{
  CH_TIME("EBAMRPoissonOp::GSColorAllRegular");

  int nComps = a_phi.nComp();
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox(m_eblg.getDBL().get(dit()));
      BaseFab<Real>& phiFAB       = (a_phi[dit()] ).getSingleValuedFAB();
      const BaseFab<Real>& rhsFAB = (a_rhs[dit()] ).getSingleValuedFAB();

      EBCellFAB& phi = a_phi[dit()];
      m_colorEBStencil[a_icolor][dit()]->cachePhi(phi);

      Box loBox[SpaceDim],hiBox[SpaceDim];
      int hasLo[SpaceDim],hasHi[SpaceDim];

      applyDomainFlux(loBox, hiBox, hasLo, hasHi,
                      dblBox, nComps, phiFAB,
                      a_homogeneousPhysBC, dit(),m_beta);

      IntVect loIV = dblBox.smallEnd();
      IntVect hiIV = dblBox.bigEnd();

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (loIV[idir] % 2 != m_colors[a_icolor][idir])
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

void EBAMRPoissonOp::
GSColorAllIrregularClone(LevelData<EBCellFAB>&        a_phi,
                         const LevelData<EBCellFAB>&  a_phiOld,
                         const LevelData<EBCellFAB>&  a_rhs,
                         const int&                   a_icolor,
                         const bool&                  a_homogeneousPhysBC)
{
  CH_TIMERS("EBAMRPoissonOp::GSColorAllIrregular");
  CH_TIMER("applyColorStencil", t1);
  CH_TIMER("vofItDomLoBCS", t2);
  CH_TIMER("vofItDomHiBCS", t3);

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      if (m_vofItIrregColor[a_icolor][dit()].size() != 0)
        {
          Box dblBox( m_eblg.getDBL().get(dit()) );

          EBCellFAB&          phi = a_phi[   dit()];
          const EBCellFAB& phiOld = a_phiOld[dit()];
          const EBCellFAB& rhs = a_rhs[dit()];

          const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[dit()];
          const BaseIVFAB<Real>& curBetaWeight  = m_betaDiagWeight[dit()];
          int comp = 0;

          CH_START(t1);
          //phi = (I-lambda*L)phiOld
          Real safety = 1.0;
          m_colorEBStencil[a_icolor][dit()]->relaxClone(phi, phiOld, rhs, curAlphaWeight, curBetaWeight, m_alpha, m_beta, safety);
          CH_STOP(t1);

          //apply domain bcs to (I-lambda*L)phi (already done in colorStencil += fluxStencil, and hom only here))
          //apply domain bcs to (I-lambda*L)phi
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              CH_START(t2);
              for (m_vofItIrregColorDomLo[a_icolor][idir][dit()].reset(); m_vofItIrregColorDomLo[a_icolor][idir][dit()].ok();  ++m_vofItIrregColorDomLo[a_icolor][idir][dit()])
                {
                  Real flux;
                  const VolIndex& vof = m_vofItIrregColorDomLo[a_icolor][idir][dit()]();
                  Real weightIrreg = m_alpha*curAlphaWeight(vof,0) + m_beta*curBetaWeight(vof,0);
                  m_domainBC->getFaceFlux(flux,vof,comp,phiOld,
                                          m_origin,m_dx,idir,Side::Lo, dit(), s_time,
                                          a_homogeneousPhysBC);

                  phi(vof,comp) += (1./weightIrreg) * flux * m_beta*m_invDx[idir];
                }
              CH_STOP(t2);
              CH_START(t3);
              for (m_vofItIrregColorDomHi[a_icolor][idir][dit()].reset(); m_vofItIrregColorDomHi[a_icolor][idir][dit()].ok();  ++m_vofItIrregColorDomHi[a_icolor][idir][dit()])
                {
                  Real flux;
                  const VolIndex& vof = m_vofItIrregColorDomHi[a_icolor][idir][dit()]();
                  Real weightIrreg = m_alpha*curAlphaWeight(vof,0) + m_beta*curBetaWeight(vof,0);
                  m_domainBC->getFaceFlux(flux,vof,comp,phiOld,
                                          m_origin,m_dx,idir,Side::Hi,dit(),s_time,
                                          a_homogeneousPhysBC);

                  phi(vof,comp) -= (1./weightIrreg) * flux * m_beta*m_invDx[idir];
                }
              CH_STOP(t3);
            }
        }
    }
}

// Used in MFElliptic code...
void EBAMRPoissonOp::
levelGSRB(LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>& a_resid,
          const int                   a_color)
{
  CH_TIME("EBAMRPoissonOp::levelGSRB");

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  int ncomp = a_phi.nComp();

  Real weight = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0 * m_beta * m_invDx2[idir];
    }
  weight = 1.0 / weight;

  int ibox = 0;
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&       phifab = a_phi[dit()];
      const EBCellFAB& residfab = a_resid[dit()];

      //reg cells
      const Box&                 region =      dbl.get(dit());
      BaseFab<Real>&         phiBaseFAB =   phifab.getSingleValuedFAB();
      const BaseFab<Real>& residBaseFAB = residfab.getSingleValuedFAB();

      for (int icomp = 0; icomp < ncomp; icomp++)
        {
          FORT_REGGSRB(CHF_FRA1(phiBaseFAB,icomp),
                       CHF_CONST_FRA1(residBaseFAB,icomp),
                       CHF_CONST_REAL(weight),
                       CHF_BOX(region),
                       CHF_CONST_INT(a_color));
        }

      // Do the irregular cells
      const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[dit()];
      const BaseIVFAB<Real>& curBetaWeight  = m_betaDiagWeight[dit()];

      VoFIterator& vofit = m_vofItIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          const IntVect&  iv = VoF.gridIndex();

          if (iv.sum()%2 == a_color)
            {
              Real weightIrreg = m_alpha*curAlphaWeight(VoF,0)
                + m_beta*curBetaWeight(VoF, 0);
              for (int icomp=0; icomp<ncomp; icomp++)
                {
                  //have to subtract off what fortran did then add on irregular stuff
                  //only do subtraction on single-valued data holder
                  //this will save us from double counting on multivalued cells
                  //since the basefab data under multivalued cells is meaningless
                  Real residVal =  residfab(VoF,icomp);
                  phiBaseFAB(iv, icomp) -= weight*residVal;
                  //now do irregular increment
                  Real increment = 0;
                  if (Abs(weightIrreg) > 1.0e-12)
                    {
                      increment = (1.0/weightIrreg)*residVal;
                    }
                  phifab(VoF,icomp) +=  increment;
                }
            }
        }
      ibox++;
    }
}

// Still used in MFElliptic code...
void EBAMRPoissonOp::
levelMultiColorGS(LevelData<EBCellFAB>&       a_phi,
                  const LevelData<EBCellFAB>& a_resid,
                  const IntVect& color)
{
  int ncomp = a_phi.nComp();

  Real weight = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0 * m_beta * m_invDx2[idir];
    }
  weight = 1.0 / weight;

  int ibox = 0;
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& phi = a_phi[dit()];
      const Box& box   = m_eblg.getDBL().get(dit());
      const EBCellFAB& resid = a_resid[dit()];

      // Do the regular cells
      BaseFab<Real>& phiFAB         =   phi.getSingleValuedFAB();
      const BaseFab<Real>& residFAB = resid.getSingleValuedFAB();

      IntVect loIV = box.smallEnd();
      IntVect hiIV = box.bigEnd();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (loIV[idir] % 2 != color[idir])
            {
              loIV[idir]++;
            }
        }

      if (loIV <= hiIV)
        {
          Box coloredBox(loIV,hiIV);
          for (int comp=0; comp<ncomp; comp++)
            {
              FORT_REGMULTICOLORGS(CHF_FRA1(phiFAB,comp),
                                   CHF_CONST_REAL(weight),
                                   CHF_CONST_FRA1(residFAB,comp),
                                   CHF_BOX(coloredBox));
            }
        }

      // Do the irregular cells
      const BaseIVFAB<Real>& curAlphaWeight = m_alphaDiagWeight[dit()];
      const BaseIVFAB<Real>& curBetaWeight  = m_betaDiagWeight[dit()];

      VoFIterator& vofit = m_vofItIrreg[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          const IntVect&  iv = VoF.gridIndex();

          bool ok = true;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (iv[idir] % 2 != color[idir])
                {
                  ok = false;
                  break;
                }
            }

          if (ok)
            {
              Real weightIrreg = m_alpha*curAlphaWeight(VoF,0) + m_beta*curBetaWeight(VoF, 0);
              for (int comp=0; comp<ncomp; comp++)
                {
                  //have to subtract off what fortran did then add on irregular stuff
                  //only do subtraction on single-valued data holder
                  //this will save us from double counting on multivalued cells
                  //since the basefab data under multivalued cells is meaningless
                  Real residVal =  resid(VoF,comp);
                  phiFAB(iv, comp) -= weight*residVal;
                  //now do irregular increment
                  Real increment = 0;
                  if (Abs(weightIrreg ) > 1.0e-12)
                    {
                      increment = (1.0/weightIrreg)*residVal;
                    }
                  phi(VoF,comp) +=  increment;
                }
            }
        }
      ibox++;
    }
}

void EBAMRPoissonOp::
levelSlowRelax(LevelData<EBCellFAB>&       a_phi,
               const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBConductivityOp::levelSlowRelax");

  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);

  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  LevelData<EBCellFAB> lphi;
  create(lphi, a_rhs);
  for (int icolor = 0; icolor < m_colors.size(); icolor++)
    {
      if (m_hasCoar)
        {
          applyHomogeneousCFBCs(a_phi);
        }

      //after this lphi = L(phi)
      //this call contains bcs and exchange
      applyOp(  lphi,  a_phi, true);
      slowGSRBColor(a_phi, lphi, a_rhs, m_colors[icolor]);
    }
}

void
EBAMRPoissonOp::
slowGSRBColor(LevelData<EBCellFAB>&       a_phi,
              const LevelData<EBCellFAB>& a_lph,
              const LevelData<EBCellFAB>& a_rhs,
              const IntVect&              a_color)
{

  LevelData<EBCellFAB> relCoef;
  Real safety = 0.5;
  create(relCoef, a_rhs);
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  Real weight = m_alpha;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0 * m_beta * m_invDx2[idir];
    }
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      //set relaxation coeff to the regular value and fix
      // it up at irregular cells
      relCoef[dit()].setVal(safety/weight);
      const EBGraph& ebgraph = m_eblg.getEBISL()[dit()].getEBGraph();
      IntVectSet irregIVS = ebgraph.getIrregCells(dbl.get(dit()));

      for (VoFIterator vofit(irregIVS, ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          Real alphaWeight  = m_alphaDiagWeight[dit()](VoF,0);
          Real  betaWeight  =  m_betaDiagWeight[dit()](VoF,0);
          Real weightIrreg = m_alpha*alphaWeight + m_beta*betaWeight;
          if (weightIrreg > 1.0e-12)
            {
              relCoef[dit()](VoF,0) = safety/weightIrreg;
            }
          else
            {
              relCoef[dit()](VoF,0) = 0.0;
            }
        }
    }

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox  = dbl.get(dit());
      BaseFab<Real>&       regPhi =     a_phi[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regLph =     a_lph[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regRhs =     a_rhs[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regRel =   relCoef[dit()].getSingleValuedFAB();
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
          FORT_SLOWGSRBEBAMRPO(CHF_FRA1(regPhi,0),
                               CHF_CONST_FRA1(regLph,0),
                               CHF_CONST_FRA1(regRhs,0),
                               CHF_CONST_FRA1(regRel,0),
                               CHF_BOX(coloredBox));
        }

      const EBGraph& ebgraph = m_eblg.getEBISL()[dit()].getEBGraph();
      IntVectSet multiIVS = ebgraph.getMultiCells(dbl.get(dit()));

      for (VoFIterator vofit(multiIVS, ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
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
              Real lambda = relCoef[dit()](vof, 0);
              a_phi[dit()](vof, 0) += lambda*resid;
            }
        }
    }
}

#include "NamespaceFooter.H"
