#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBTensorCFInterp.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "EBAlias.H"
#include "EBArith.H"
#include "EBCellFactory.H"
#include "EBTensorCFInterp.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"

/***********************/
// default constructor
/***********************/
EBTensorCFInterp::
EBTensorCFInterp(const DisjointBoxLayout&       a_gridsFine,
                 const DisjointBoxLayout&       a_gridsCoar,
                 const EBISLayout&              a_ebislFine,
                 const EBISLayout&              a_ebislCoar,
                 const ProblemDomain&           a_domainCoar,
                 const int&                     a_nref,
                 const int&                     a_nvar,
                 const Real&                    a_dxFine,
                 const LayoutData<IntVectSet>&  a_cfivs,
                 const EBIndexSpace* const      a_ebisPtr,
                 bool                           a_doEBCFCrossing)
  :TensorCFInterp(a_gridsFine, &a_gridsCoar, a_dxFine, a_nref, a_nvar, refine(a_domainCoar, a_nref))
{
  m_nComp  = a_nvar;
  m_refRat = a_nref;
  m_domainCoar = a_domainCoar;
  m_domainFine = refine(m_domainCoar, a_nref);

  m_isDefined = true;

  m_ebquadcfi = RefCountedPtr<EBQuadCFInterp>
    (new EBQuadCFInterp(a_gridsFine, a_gridsCoar,
                        a_ebislFine, a_ebislCoar, a_domainCoar,
                        a_nref,a_nvar, a_cfivs,a_ebisPtr,a_doEBCFCrossing));

  m_ebcfdata = m_ebquadcfi->getEBCFData();
  m_doEBCFCrossing = m_ebcfdata->m_doEBCFCrossing;
  //define coarse ebis layouts buffers
  EBCellFactory coarFact(m_ebcfdata->m_ebislCoarsenedFine);
  m_ebBufferCoarsenedFine.define(m_ebcfdata->m_gridsCoarsenedFine, m_nComp, 4*IntVect::Unit, coarFact);
  if (m_doEBCFCrossing)
    {
      //build EBCF stencils
      buildEBCFCrossingStencils(a_cfivs);
    }
  buildEBCFCornerStencils(a_cfivs);
}
/***********************/
void
EBTensorCFInterp::
buildEBCFCrossingStencils(const LayoutData<IntVectSet>& a_cfivs)
{
  Real dxCoar = m_dxFine*m_refRatio;
  int numCompGrad = m_nComp*SpaceDim;
  for (int facedir = 0; facedir < SpaceDim; facedir++)
    {
      m_coarStencilLo[facedir].define(m_ebcfdata->m_gridsFine);
      m_coarStencilHi[facedir].define(m_ebcfdata->m_gridsFine);
      for (DataIterator dit = m_ebcfdata->m_gridsFine.dataIterator(); dit.ok(); ++dit)
        {
          //define all the holders
          const IntVectSet& loEBCFIVS = m_ebcfdata->m_ebcfivsLo[facedir][dit()];
          const IntVectSet& hiEBCFIVS = m_ebcfdata->m_ebcfivsHi[facedir][dit()];
          const EBISBox& ebisBoxCF = m_ebcfdata->m_ebislCoarsenedFine[dit()];

          m_coarStencilLo[facedir][dit()].define(loEBCFIVS, m_ebcfdata->m_ebislFine[dit()].getEBGraph(), numCompGrad);
          m_coarStencilHi[facedir][dit()].define(hiEBCFIVS, m_ebcfdata->m_ebislFine[dit()].getEBGraph(), numCompGrad);

          for (int ivar = 0; ivar < m_nComp; ivar++)
            {
              for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
                {
                  int icomp = TensorCFInterp::gradIndex(ivar, diffDir);

                  for (m_ebcfdata->m_vofItEBCFLo[facedir][dit()].reset(); m_ebcfdata->m_vofItEBCFLo[facedir][dit()].ok(); ++m_ebcfdata->m_vofItEBCFLo[facedir][dit()])
                    {
                      const VolIndex& vofFine = (m_ebcfdata->m_vofItEBCFLo[facedir][dit()])();
                      VolIndex vofCoar = m_ebcfdata->m_ebislFine.coarsen(vofFine, m_refRat, dit());
                      EBArith::getFirstDerivStencil(m_coarStencilLo[facedir][dit()](vofFine, icomp), vofCoar, ebisBoxCF,
                                                    diffDir, dxCoar, NULL, ivar);
                    }

                  for (m_ebcfdata->m_vofItEBCFHi[facedir][dit()].reset(); m_ebcfdata->m_vofItEBCFHi[facedir][dit()].ok(); ++m_ebcfdata->m_vofItEBCFHi[facedir][dit()])
                    {
                      const VolIndex& vofFine = (m_ebcfdata->m_vofItEBCFHi[facedir][dit()])();
                      VolIndex vofCoar = m_ebcfdata->m_ebislFine.coarsen(vofFine, m_refRat, dit());
                      EBArith::getFirstDerivStencil(m_coarStencilHi[facedir][dit()](vofFine, icomp), vofCoar, ebisBoxCF,
                                                    diffDir, dxCoar, NULL, ivar);
                    }
                }
            }
        }
    }
}
/***********************/

void
EBTensorCFInterp::
buildEBCFCornerStencils(const LayoutData<IntVectSet>& a_cfivs)
{
  int numCompGrad = m_nComp*SpaceDim;
  Real dxCoar = m_dxFine*m_refRatio;
  // let us get the corner IVS
  m_stencilCorners.define(m_ebcfdata->m_gridsFine);
  m_stencilEdges.define(m_ebcfdata->m_gridsFine);

  for (DataIterator dit = m_ebcfdata->m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      const IntVectSet& cornerIVS = m_ebcfdata->m_cornerIVS[dit()];
      const IntVectSet&   edgeIVS = m_ebcfdata->m_edgeIVS[dit()];
      const EBISBox& ebisBoxCF = m_ebcfdata->m_ebislCoarsenedFine[dit()];

      m_stencilCorners[dit()].define(cornerIVS, m_ebcfdata->m_ebislFine[dit()].getEBGraph(), numCompGrad);
      m_stencilEdges[dit()].define(    edgeIVS, m_ebcfdata->m_ebislFine[dit()].getEBGraph(), numCompGrad);

      for (int ivar = 0; ivar < m_nComp; ivar++)
        {
          for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
            {
              int icomp = TensorCFInterp::gradIndex(ivar, diffDir);
              for (m_ebcfdata->m_vofItCorners[dit()].reset(); m_ebcfdata->m_vofItCorners[dit()].ok(); ++m_ebcfdata->m_vofItCorners[dit()])
                {
                  const VolIndex& vofFine = m_ebcfdata->m_vofItCorners[dit()]();
                  VolIndex vofCoar = m_ebcfdata->m_ebislFine.coarsen(vofFine, m_refRat, dit());

                  EBArith::getFirstDerivStencil(m_stencilCorners[dit()](vofFine, icomp), vofCoar, ebisBoxCF,
                                                diffDir, dxCoar, NULL, ivar);
                }
              for (m_ebcfdata->m_vofItEdges[dit()].reset(); m_ebcfdata->m_vofItEdges[dit()].ok(); ++m_ebcfdata->m_vofItEdges[dit()])
                {
                  const VolIndex& vofFine = m_ebcfdata->m_vofItEdges[dit()]();
                  VolIndex vofCoar = m_ebcfdata->m_ebislFine.coarsen(vofFine, m_refRat, dit());
                  EBArith::getFirstDerivStencil(m_stencilEdges[dit()](vofFine, icomp), vofCoar, ebisBoxCF,
                                                diffDir, dxCoar, NULL, ivar);
                }
            }

        }
    }
}
/***********************/
void
EBTensorCFInterp::
coarseFineInterp(LevelData<EBCellFAB>&       a_fineData,
                 LevelData<EBCellFAB>&       a_tanGradF,
                 const LevelData<EBCellFAB>& a_coarData)
{

  CH_assert(a_fineData.nComp() == m_nComp);
  CH_assert(a_coarData.nComp() == m_nComp);
  CH_assert(a_tanGradF.nComp() == SpaceDim*m_nComp);
  LevelData<FArrayBox> fineDataLDFAB, coarDataLDFAB, tanGradFLDFAB;

  aliasEB(fineDataLDFAB, a_fineData);
  aliasEB(tanGradFLDFAB, a_tanGradF);
  aliasEB(coarDataLDFAB, (LevelData<EBCellFAB>&)a_coarData);
  TensorCFInterp::coarseFineInterp(fineDataLDFAB, tanGradFLDFAB, coarDataLDFAB);


  if (m_doEBCFCrossing)
    {
      interpEBCFCrossing(a_fineData, a_tanGradF, a_coarData);
    }

  interpEBCFCorners(     a_fineData, a_tanGradF, a_coarData);


}
void
EBTensorCFInterp::
coarseFineInterpH(LevelData<EBCellFAB>& a_fineData,
                  LevelData<EBCellFAB>& a_tanGradF)
{
  EBCellFactory        fact(m_ebcfdata->m_ebislCoar);
  LevelData<EBCellFAB> zero(m_ebcfdata->m_gridsCoar, a_fineData.nComp(), a_fineData.ghostVect(), fact);
  EBLevelDataOps::setVal(zero, 0.0);
  coarseFineInterp(a_fineData, a_tanGradF, zero);
}
/***********************/
EBTensorCFInterp::~EBTensorCFInterp()
{
}

void
EBTensorCFInterp::
interpEBCFCrossing(LevelData<EBCellFAB>&       a_fineData,
                   LevelData<EBCellFAB>&       a_tanGradF,
                   const LevelData<EBCellFAB>& a_coarData)
{
  //fill the data
  Interval interv(0, m_nComp-1);
  m_ebquadcfi->interpEBCFCrossing(a_fineData, a_coarData, interv);
  //fill the gradients
  a_coarData.copyTo(interv, m_ebBufferCoarsenedFine,interv);
  for (int ivar = 0; ivar < m_nComp; ivar++)
    {
      for (int idiff = 0; idiff < SpaceDim; idiff++)
        {
          int icomp = TensorCFInterp::gradIndex(ivar, idiff);

          for (int facedir = 0; facedir < SpaceDim; facedir++)
            {
              for (DataIterator dit = m_ebcfdata->m_gridsFine.dataIterator(); dit.ok(); ++dit)
                {
                  for (m_ebcfdata->m_vofItEBCFLo[facedir][dit()].reset(); m_ebcfdata->m_vofItEBCFLo[facedir][dit()].ok(); ++m_ebcfdata->m_vofItEBCFLo[facedir][dit()])
                    {
                      const   VolIndex& vofGhost =  (m_ebcfdata->m_vofItEBCFLo[facedir][dit()])();
                      const VoFStencil& coarSten = m_coarStencilLo[facedir][dit()](vofGhost, icomp);

                      Real coarContrib = applyVoFStencil(coarSten, m_ebBufferCoarsenedFine[dit()], icomp);
                      a_tanGradF[dit()](vofGhost, icomp) = coarContrib;

                    }
                  for (m_ebcfdata->m_vofItEBCFHi[facedir][dit()].reset(); m_ebcfdata->m_vofItEBCFHi[facedir][dit()].ok(); ++m_ebcfdata->m_vofItEBCFHi[facedir][dit()])
                    {
                      const   VolIndex& vofGhost =  (m_ebcfdata->m_vofItEBCFHi[facedir][dit()])();
                      const VoFStencil& coarSten = m_coarStencilHi[facedir][dit()](vofGhost, icomp);

                      Real coarContrib = applyVoFStencil(coarSten, m_ebBufferCoarsenedFine[dit()], icomp);
                      a_tanGradF[dit()](vofGhost, icomp) = coarContrib;
                    }
                }
            }
        }
    }
}
void
EBTensorCFInterp::
interpEBCFCorners(LevelData<EBCellFAB>&       a_fineData,
                  LevelData<EBCellFAB>&       a_tanGradF,
                  const LevelData<EBCellFAB>& a_coarData)
{
  Interval interv(0, m_nComp-1);
  //fill the data
  m_ebquadcfi->interpEBCFCorners(a_fineData, a_coarData, interv);
  //fill the gradients
  a_coarData.copyTo(interv, m_ebBufferCoarsenedFine, interv);

  for (int ivar = 0; ivar < m_nComp; ivar++)
    {
      for (int idiff = 0; idiff < SpaceDim; idiff++)
        {
          int icomp = TensorCFInterp::gradIndex(ivar, idiff);

          for (DataIterator dit = m_ebcfdata->m_gridsFine.dataIterator(); dit.ok(); ++dit)
            {
              for (m_ebcfdata->m_vofItCorners[dit()].reset(); m_ebcfdata->m_vofItCorners[dit()].ok(); ++m_ebcfdata->m_vofItCorners[dit()])
                {
                  const   VolIndex& vofGhost =  (m_ebcfdata->m_vofItCorners[dit()])();
                  const VoFStencil& coarSten = m_stencilCorners[dit()](vofGhost, icomp);

                  Real coarContrib = applyVoFStencil(coarSten, m_ebBufferCoarsenedFine[dit()], icomp);
                  a_tanGradF[dit()](vofGhost, icomp) = coarContrib;

                }
              for (m_ebcfdata->m_vofItEdges[dit()].reset(); m_ebcfdata->m_vofItEdges[dit()].ok(); ++m_ebcfdata->m_vofItEdges[dit()])
                {
                  const   VolIndex& vofGhost =  (m_ebcfdata->m_vofItEdges[dit()])();
                  const VoFStencil& coarSten = m_stencilEdges[dit()](vofGhost, icomp);

                  Real coarContrib = applyVoFStencil(coarSten, m_ebBufferCoarsenedFine[dit()], icomp);
                  a_tanGradF[dit()](vofGhost, icomp) = coarContrib;
                }
            }
        }
    }
}

#include "NamespaceFooter.H"
