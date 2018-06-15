#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBStencil.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "NamespaceHeader.H"
/**************/
/**************/
/**************/
EBStencil::EBStencil(const Vector<VolIndex>& a_srcVofs,
                     const BaseIVFAB<VoFStencil>& a_vofStencil,
                     const Box& a_box,
                     const EBISBox& a_ebisBox,
                     const IntVect& a_ghostVectPhi,
                     const IntVect& a_ghostVectLph,
                     int a_varDest,
                     bool a_doRelaxOpt,
                     int nComp, 
                     IntVectSet a_setIrreg,
                     bool       a_useInputSet)
  : m_box( a_box ),
    m_ebisBox( a_ebisBox ),
    m_ghostVectPhi( a_ghostVectPhi ),
    m_ghostVectLph( a_ghostVectLph ),
    m_destVar(a_varDest),
    m_doRelaxOpt(a_doRelaxOpt),
    m_nComp(nComp),
    m_setIrreg(a_setIrreg),
    m_useInputSets(a_useInputSet)
{
  CH_TIMERS("EBStencil::EBStencil");
  CH_TIMER("computeOffsets", t1);
  CH_START(t1);
  computeOffsets(a_srcVofs, a_vofStencil);
  CH_STOP(t1);
}
/**************/
/**************/
EBStencil::EBStencil(const Vector<VolIndex>& a_srcVofs,
                     const Vector<VoFStencil>& a_vofStencil,
                     const Box& a_boxLph,
                     const Box& a_boxPhi,
                     const EBISBox& a_ebisBoxLph,
                     const EBISBox& a_ebisBoxPhi,
                     const IntVect& a_ghostVectLph,
                     const IntVect& a_ghostVectPhi,
                     int a_varDest,
                     int nComp,
                     IntVectSet a_setIrreg,
                     bool       a_useInputSet)
  : m_ghostVectPhi( a_ghostVectPhi ),
    m_ghostVectLph( a_ghostVectLph ),
    m_destVar(a_varDest),
    m_doRelaxOpt(false),
    m_nComp(nComp),
    m_setIrreg(a_setIrreg),
    m_useInputSets(a_useInputSet)
{
//  CH_TIMERS("EBStencil::EBStencil");
//  CH_TIMER("computeOffsets", t1);
//  CH_START(t1);

  Box boxPhi = grow(a_boxPhi, a_ghostVectPhi);
  Box boxLph = grow(a_boxLph, a_ghostVectLph);
  m_lphBox = boxLph;
  m_phiBox = boxPhi;
  //debugging hook
  //m_srcVoFs = a_srcVofs;
  const IntVectSet& ivsPhi = a_ebisBoxPhi.getMultiCells(boxPhi);
  const IntVectSet& ivsLph = a_ebisBoxLph.getMultiCells(boxLph);
 
  const EBGraph& ebgraphPhi = a_ebisBoxPhi.getEBGraph();
  const EBGraph& ebgraphLph = a_ebisBoxLph.getEBGraph();
  BaseIVFAB<Real> baseivfabPhi(ivsPhi, ebgraphPhi, m_nComp);
  BaseIVFAB<Real> baseivfabLph(ivsLph, ebgraphLph, m_nComp);

  const IntVect& smallendPhi = boxPhi.smallEnd();
  const IntVect& smallendLph = boxLph.smallEnd();
  IntVect ncellsPhi = boxPhi.size();
  IntVect ncellsLph = boxLph.size();
  m_ebstencil.resize(a_srcVofs.size());
  m_destTerms.resize(a_srcVofs.size());
  m_cacheLph.resize(a_srcVofs.size());
  m_cachePhi.resize(a_srcVofs.size());

  for (int isrc = 0; isrc < a_srcVofs.size(); isrc++)
    {
      const VolIndex& srcVof = a_srcVofs[isrc];
      if (a_ebisBoxLph.numVoFs(srcVof.gridIndex()) > 1)
        {//multi-valued (the dataPtr(0) is correct--that is where we start from)
          m_destTerms[isrc].offset = baseivfabLph.getIndex(srcVof, m_destVar) - baseivfabLph.dataPtr(0);
          m_destTerms[isrc].multiValued = true;
        }
      else
        {//single-valued
          IntVect ivLph = srcVof.gridIndex()  - smallendLph;
          IntVect ivPhi = srcVof.gridIndex()  - smallendPhi;
          m_destTerms[isrc].offset = ivLph[0] + ivLph[1]*ncellsLph[0] ;
#if CH_SPACEDIM==3
          m_destTerms[isrc].offset +=  ivLph[2]*ncellsLph[0]*ncellsLph[1];
#endif

          //add in term due to variable number
#if CH_SPACEDIM==2
          m_destTerms[isrc].offset += m_destVar*ncellsLph[0]*ncellsLph[1];
#elif CH_SPACEDIM==3
          m_destTerms[isrc].offset += m_destVar*ncellsLph[0]*ncellsLph[1]*ncellsLph[2];
#else
          bogus_spacedim();
#endif
          m_destTerms[isrc].multiValued = false;
        }
      const VoFStencil& sten = a_vofStencil[isrc];
      for (int isten = 0; isten < sten.size(); isten++)
        {
          const VolIndex stencilVof = sten.vof(isten);
          int srcVar = sten.variable(isten);
          stencilTerm stenEntry;
          stenEntry.weight = sten.weight(isten);
          if (a_ebisBoxPhi.numVoFs(stencilVof.gridIndex()) > 1)
            {//multi-valued (the dataPtr(0) is correct--that is where we start from)
              stenEntry.offset = baseivfabPhi.getIndex(stencilVof, srcVar) - baseivfabPhi.dataPtr(0);
              m_ebstencil[isrc].multi.push_back(stenEntry);
            }
          else
            {//single-valued
              IntVect ivPhi = stencilVof.gridIndex()  - smallendPhi;
              stenEntry.offset = ivPhi[0] + ivPhi[1]*ncellsPhi[0] ;
#if CH_SPACEDIM==3
              stenEntry.offset +=  ivPhi[2]*ncellsPhi[0]*ncellsPhi[1];
#endif
              //add in term due to variable number
#if CH_SPACEDIM==2
              stenEntry.offset += srcVar*ncellsPhi[0]*ncellsPhi[1];
#elif CH_SPACEDIM==3
              stenEntry.offset += srcVar*ncellsPhi[0]*ncellsPhi[1]*ncellsPhi[2];
#else
              bogus_spacedim();
#endif
              m_ebstencil[isrc].single.push_back(stenEntry);
            }
        }
    }
  //  CH_STOP(t1);
}
/**************/
/**************/

void EBStencil::apply(EBCellFAB& a_lofphi, const EBCellFAB& a_phi, bool a_incrementOnly, int  a_ivar) const
{

  CH_TIMERS("EBStencil::apply");
  CH_TIMER("apply_loop_overvofs", t3);
  CH_TIMER("apply_header", t4);
  CH_START(t4);

  CH_assert(a_lofphi.getSingleValuedFAB().box() == m_lphBox);
  CH_assert(a_phi.getSingleValuedFAB().box()    == m_phiBox);

  const Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(a_ivar);
  Real*       singleValuedPtrLph = a_lofphi.getSingleValuedFAB().dataPtr(a_ivar);

  const Real*  multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(a_ivar);
  Real*        multiValuedPtrLph = a_lofphi.getMultiValuedFAB().dataPtr(a_ivar);

  CH_STOP(t4);
  CH_START(t3);
  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      const ebstencil_t& ebstencil = m_ebstencil[isrc];
      Real* lphiPtr = NULL;
      if (m_destTerms[isrc].multiValued)
        {
          lphiPtr  = multiValuedPtrLph + m_destTerms[isrc].offset;
        }
      else
        {
          lphiPtr  = singleValuedPtrLph + m_destTerms[isrc].offset;
        }

      Real& lphi = *lphiPtr;
      if (!a_incrementOnly)
        {
          lphi =  0.;
        }
      //single-valued
      for (int isingle = 0; isingle < ebstencil.single.size(); isingle++)
        {
          const int & offset = ebstencil.single[isingle].offset;
          const Real& phiVal = *(singleValuedPtrPhi + offset);
          const Real& weight = ebstencil.single[isingle].weight;
          lphi += phiVal*weight;
        }
      //multi-valued
      for (int imulti = 0; imulti < ebstencil.multi.size(); imulti++)
        {
          const int & offset = ebstencil.multi[imulti].offset;
          const Real& phiVal = *(multiValuedPtrPhi + offset);
          const Real& weight = ebstencil.multi[imulti].weight;
          lphi += phiVal*weight;
        }
    }
  CH_STOP(t3);
}

void EBStencil::apply(EBCellFAB&             a_lofphi,
                      const EBCellFAB&       a_phi,
                      const BaseIVFAB<Real>& a_alphaWeight,
                      Real                   a_alpha,
                      Real                   a_beta,
                      bool                   a_incrementOnly) const

{
  if (!m_doRelaxOpt)
    {
      MayDay::Error("this ebstencil was not configured to deal with baseivfab alpha");
    }
  CH_TIMERS("EBStencil::apply_alpha_beta");
  CH_TIMER("apply_alpha_beta_loop_overvofs", t3);
  CH_TIMER("apply_alpha_beta_header", t4);
  CH_START(t4);

  CH_assert(a_lofphi.getSingleValuedFAB().box() == m_lphBox);
  CH_assert(a_phi.getSingleValuedFAB().box()    == m_phiBox);

  const Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(0);
  Real*       singleValuedPtrLph = a_lofphi.getSingleValuedFAB().dataPtr(0);

  const Real* multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(0);
  Real*       multiValuedPtrLph = a_lofphi.getMultiValuedFAB().dataPtr(0);

  const Real* alphaWeightPtr = a_alphaWeight.dataPtr(0);

  CH_STOP(t4);
  CH_START(t3);
  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      const ebstencil_t& ebstencil = m_ebstencil[isrc];
      Real* lphiPtr = NULL;
      const Real* sourPtr = NULL;
      int alphaOffset = m_alphaBeta[isrc];
      const Real& alphaWeight = *(alphaWeightPtr + alphaOffset);
      if (m_destTerms[isrc].multiValued)
        {
          lphiPtr  = multiValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = multiValuedPtrPhi + m_sourTerms[isrc].offset;
        }
      else
        {
          lphiPtr  = singleValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = singleValuedPtrPhi + m_sourTerms[isrc].offset;
        }

      Real& lphi = *lphiPtr;
      const Real& sour = *sourPtr;
      if (!a_incrementOnly)
        {
          lphi =  0.;
        }
      //single-valued
      for (int isingle = 0; isingle < ebstencil.single.size(); isingle++)
        {
          const int & offset = ebstencil.single[isingle].offset;
          const Real& phiVal = *(singleValuedPtrPhi + offset);
          const Real& weight = ebstencil.single[isingle].weight;
          lphi += phiVal*weight;
        }
      //multi-valued
      for (int imulti = 0; imulti < ebstencil.multi.size(); imulti++)
        {
          const int & offset = ebstencil.multi[imulti].offset;
          const Real& phiVal = *(multiValuedPtrPhi + offset);
          const Real& weight = ebstencil.multi[imulti].weight;
          lphi += phiVal*weight;
        }
      lphi = a_alpha*alphaWeight*sour + a_beta*lphi;
    }
  CH_STOP(t3);
}

//For EB x domain where m_alpha, m_beta have changed since defineStencils
//we need both alphaWeight and betaWeight to calculate the relaxation parameter:
//lambdaDiagWeight = 1/(alpha*alphaWeight+beta*betaWeight)
void EBStencil::apply(EBCellFAB&             a_lofphi,
                      const EBCellFAB&       a_phi,
                      const Real             a_lambdaFactor,
                      const Real             a_alpha,
                      const BaseIVFAB<Real>& a_alphaWeight,
                      const Real             a_beta,
                      const BaseIVFAB<Real>& a_betaWeight,
                      Real                   a_one,
                      bool                   a_incrementOnly) const

{
  if (!m_doRelaxOpt)
    {
      MayDay::Error("this ebstencil was not configured to deal with baseivfab alpha");
    }
  CH_TIMERS("EBStencil::apply_alpha_beta");
  CH_TIMER("apply_alpha_beta_loop_overvofs", t3);
  CH_TIMER("apply_alpha_beta_header", t4);
  CH_START(t4);

  CH_assert(a_lofphi.getSingleValuedFAB().box() == m_lphBox);
  CH_assert(a_phi.getSingleValuedFAB().box()    == m_phiBox);

  const Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(0);
  Real*       singleValuedPtrLph = a_lofphi.getSingleValuedFAB().dataPtr(0);

  const Real* multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(0);
  Real*       multiValuedPtrLph = a_lofphi.getMultiValuedFAB().dataPtr(0);

  const Real* alphaWeightPtr = a_alphaWeight.dataPtr(0);
  const Real* betaWeightPtr = a_betaWeight.dataPtr(0);

  CH_STOP(t4);
  CH_START(t3);
  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];

      //const ebstencil_t& ebstencil = m_ebstencil[isrc];
      Real* lphiPtr = NULL;
      const Real* sourPtr = NULL;
      int alphaOffset = m_alphaBeta[isrc];
      int betaOffset = m_alphaBeta[isrc];
      const Real& alphaWeight = *(alphaWeightPtr + alphaOffset);
      const Real& betaWeight = *(betaWeightPtr + betaOffset);
      Real product = a_alpha*alphaWeight+a_beta*betaWeight;
      Real lambdaWeight = 0.;
      if (Abs(product) > 1.e-15) lambdaWeight = 1./product;
      if (m_destTerms[isrc].multiValued)
        {
          lphiPtr  = multiValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = multiValuedPtrPhi + m_sourTerms[isrc].offset;
        }
      else
        {
          lphiPtr  = singleValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = singleValuedPtrPhi + m_sourTerms[isrc].offset;
        }

      Real& lphi = *lphiPtr;
      const Real& sour = *sourPtr;
      if (!a_incrementOnly)
        {
          lphi =  0.;
        }
      // //single-valued
      // for (int isingle = 0; isingle < ebstencil.single.size(); isingle++)
      //   {
      //     const int & offset = ebstencil.single[isingle].offset;
      //     const Real& phiVal = *(singleValuedPtrPhi + offset);
      //     const Real& weight = ebstencil.single[isingle].weight;
      //     lphi += phiVal*weight;
      //   }
      // //multi-valued
      // for (int imulti = 0; imulti < ebstencil.multi.size(); imulti++)
      //   {
      //     const int & offset = ebstencil.multi[imulti].offset;
      //     const Real& phiVal = *(multiValuedPtrPhi + offset);
      //     const Real& weight = ebstencil.multi[imulti].weight;
      //     lphi += phiVal*weight;
      //   }

      // if (sour != 0.) pout() << srcVoF << endl;
      // if (sour != 0.) pout() << "            old phi " << lphi << endl;

      // lphi = a_lambdaFactor*lambdaWeight*sour + a_one*lphi;
      lphi += a_lambdaFactor*lambdaWeight*sour;
  }
  CH_STOP(t3);
}

void EBStencil::applyInhomDomBC(EBCellFAB&             a_lofphi,
                                const EBCellFAB&       a_phi,
                                const Real             a_factor) const

{
  CH_assert(a_lofphi.getSingleValuedFAB().box() == m_lphBox);
  CH_assert(a_phi.getSingleValuedFAB().box()    == m_phiBox);

  const Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(0);
  Real*       singleValuedPtrLph = a_lofphi.getSingleValuedFAB().dataPtr(0);

  const Real* multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(0);
  Real*       multiValuedPtrLph = a_lofphi.getMultiValuedFAB().dataPtr(0);

  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];

      //const ebstencil_t& ebstencil = m_ebstencil[isrc];
      Real* lphiPtr = NULL;
      const Real* sourPtr = NULL;
      if (m_destTerms[isrc].multiValued)
        {
          lphiPtr  = multiValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = multiValuedPtrPhi + m_sourTerms[isrc].offset;
        }
      else
        {
          lphiPtr  = singleValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = singleValuedPtrPhi + m_sourTerms[isrc].offset;
        }

      Real& lphi = *lphiPtr;
      const Real& sour = *sourPtr;
      lphi += a_factor*sour;
    }
}

void EBStencil::relax(EBCellFAB&             a_phi,
                      const EBCellFAB&       a_rhs,
                      const BaseIVFAB<Real>& a_alphaWeight,
                      const BaseIVFAB<Real>& a_betaWeight,
                      Real                   a_alpha,
                      Real                   a_beta,
                      Real                   a_safety) const
{
  if (!m_doRelaxOpt)
    {
      MayDay::Error("ebstencil::relax this ebstencil was not configured to deal with baseivfab alpha and beta");
    }

  CH_TIMERS("EBStencil::relax");
  CH_TIMER("relax_loop_overvofs", t3);
  CH_TIMER("relax_header", t4);
  CH_START(t4);

  CH_assert(a_rhs.getSingleValuedFAB().box() == m_lphBox);
  CH_assert(a_phi.getSingleValuedFAB().box() == m_phiBox);

  Real*       singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(0);
  const Real* singleValuedPtrLph =    a_rhs.getSingleValuedFAB().dataPtr(0);

  Real*       multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(0);
  const Real* multiValuedPtrLph =    a_rhs.getMultiValuedFAB().dataPtr(0);

  const Real* alphaWeightPtr = a_alphaWeight.dataPtr(0);
  const Real*  betaWeightPtr =  a_betaWeight.dataPtr(0);

  CH_STOP(t4);
  CH_START(t3);
  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      const ebstencil_t& ebstencil = m_ebstencil[isrc];
      const Real* rhsPtr = NULL;
      Real* sourPtr = NULL;
      //alpha and beta get the same offset
      int alphaOffset = m_alphaBeta[isrc];
      const Real& alphaWeight = *(alphaWeightPtr + alphaOffset);
      const Real&  betaWeight = *( betaWeightPtr + alphaOffset);
      if (m_destTerms[isrc].multiValued)
        {
          rhsPtr   = multiValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = multiValuedPtrPhi + m_sourTerms[isrc].offset;
        }
      else
        {
          rhsPtr   = singleValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = singleValuedPtrPhi + m_sourTerms[isrc].offset;
        }

      Real& phi = *sourPtr;
      const Real& rhs = *rhsPtr;

      Real denom = a_alpha*alphaWeight + a_beta*betaWeight;
      Real lambda = 0;
      if (Abs(denom) > 1.0e-12)
        lambda = a_safety/(denom);

      Real lphi = a_alpha*alphaWeight*phi;
      //single-valued
      for (int isingle = 0; isingle < ebstencil.single.size(); isingle++)
        {
          const int & offset = ebstencil.single[isingle].offset;
          const Real& phiVal = *(singleValuedPtrPhi + offset);
          const Real& weight = ebstencil.single[isingle].weight;
          lphi += phiVal*weight*a_beta;
        }
      //multi-valued
      for (int imulti = 0; imulti < ebstencil.multi.size(); imulti++)
        {
          const int & offset = ebstencil.multi[imulti].offset;
          const Real& phiVal = *(multiValuedPtrPhi + offset);
          const Real& weight = ebstencil.multi[imulti].weight;
          lphi += phiVal*weight*a_beta;
        }

      phi = phi + lambda * (rhs - lphi);
    }
  CH_STOP(t3);
}
void EBStencil::relaxClone(EBCellFAB&             a_phi,
                           const EBCellFAB&       a_phiOld,
                           const EBCellFAB&       a_rhs,
                           const BaseIVFAB<Real>& a_alphaWeight,
                           const BaseIVFAB<Real>& a_betaWeight,
                           Real a_alpha, Real a_beta, Real a_safety) const
{
  if (!m_doRelaxOpt)
    {
      MayDay::Error("ebstencil::relax this ebstencil was not configured to deal with baseivfab alpha and beta");
    }

  CH_TIMERS("EBStencil::relaxClone");
  CH_TIMER("relaxClone_loop_overvofs", t3);
  CH_TIMER("relaxClone_header", t4);
  CH_START(t4);

  CH_assert(a_rhs.getSingleValuedFAB().box() == m_lphBox);
  CH_assert(a_phi.getSingleValuedFAB().box() == m_phiBox);
  CH_assert(a_phiOld.getSingleValuedFAB().box() == m_phiBox);

  Real*       singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(0);
  const Real* singleValuedPtrPhiOld =    a_phiOld.getSingleValuedFAB().dataPtr(0);
  const Real* singleValuedPtrLph =    a_rhs.getSingleValuedFAB().dataPtr(0);

  Real*       multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(0);
  const Real*       multiValuedPtrPhiOld =    a_phiOld.getMultiValuedFAB().dataPtr(0);
  const Real* multiValuedPtrLph =    a_rhs.getMultiValuedFAB().dataPtr(0);

  const Real* alphaWeightPtr = a_alphaWeight.dataPtr(0);
  const Real*  betaWeightPtr =  a_betaWeight.dataPtr(0);

  CH_STOP(t4);
  CH_START(t3);
  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      const ebstencil_t& ebstencil = m_ebstencil[isrc];
      const Real* rhsPtr = NULL;
      Real* sourPtr = NULL;
      const Real* sourOldPtr = NULL;
      //alpha and beta get the same offset
      int alphaOffset = m_alphaBeta[isrc];
      const Real& alphaWeight = *(alphaWeightPtr + alphaOffset);
      const Real&  betaWeight = *( betaWeightPtr + alphaOffset);
      if (m_destTerms[isrc].multiValued)
        {
          rhsPtr   = multiValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = multiValuedPtrPhi + m_sourTerms[isrc].offset;
          sourOldPtr  = multiValuedPtrPhiOld + m_sourTerms[isrc].offset;
        }
      else
        {
          rhsPtr   = singleValuedPtrLph + m_destTerms[isrc].offset;
          sourPtr  = singleValuedPtrPhi + m_sourTerms[isrc].offset;
          sourOldPtr  = singleValuedPtrPhiOld + m_sourTerms[isrc].offset;
        }

      Real& phi = *sourPtr;
      const Real& phiOld = *sourOldPtr;
      const Real& rhs = *rhsPtr;

      Real denom = a_alpha*alphaWeight + a_beta*betaWeight;
      Real lambda = 0;
      if (Abs(denom) > 1.0e-12)
        lambda = a_safety/(denom);

      Real lphi = lambda*a_alpha*alphaWeight*phiOld;
      //single-valued
      for (int isingle = 0; isingle < ebstencil.single.size(); isingle++)
        {
          const int & offset = ebstencil.single[isingle].offset;
          const Real& phiVal = *(singleValuedPtrPhiOld + offset);
          const Real& weight = ebstencil.single[isingle].weight;
          lphi += lambda*phiVal*weight*a_beta;
        }
      //multi-valued
      for (int imulti = 0; imulti < ebstencil.multi.size(); imulti++)
        {
          const int & offset = ebstencil.multi[imulti].offset;
          const Real& phiVal = *(multiValuedPtrPhiOld + offset);
          const Real& weight = ebstencil.multi[imulti].weight;
          lphi += lambda*phiVal*weight*a_beta;
        }

      //lphi already has lambda multiplied in
      phi = phiOld + lambda*rhs - lphi;
    }
  CH_STOP(t3);
}
/**************/
/**************/
void
EBStencil::computeOffsets(const Vector<VolIndex>&      a_srcVofs,
                          const BaseIVFAB<VoFStencil>& a_vofStencil)
{
  CH_assert((SpaceDim == 2) || (SpaceDim == 3));
  Box boxPhi = grow(m_box, m_ghostVectPhi);
  Box boxLph = grow(m_box, m_ghostVectLph);
  m_lphBox = boxLph;
  m_phiBox = boxPhi;
  const IntVectSet& ivsPhi = m_ebisBox.getMultiCells(boxPhi);
  const IntVectSet& ivsLph = m_ebisBox.getMultiCells(boxLph);

  const EBGraph& ebgraph = m_ebisBox.getEBGraph();
  BaseIVFAB<Real> baseivfabPhi(ivsPhi, ebgraph, m_nComp);
  BaseIVFAB<Real> baseivfabLph(ivsLph, ebgraph, m_nComp);


  if (m_doRelaxOpt)
    {
      m_alphaBeta.resize(a_srcVofs.size());
      IntVectSet ivsIrr = m_ebisBox.getIrregIVS(m_box);
      if (m_useInputSets)
        {
          ivsIrr = m_setIrreg;
        }
      BaseIVFAB<Real> baseivfabIrr(ivsIrr, ebgraph, 1);
      for (int isrc = 0; isrc < a_srcVofs.size(); isrc++)
        {
          const VolIndex& srcVof = a_srcVofs[isrc];
          m_alphaBeta[isrc] = baseivfabIrr.getIndex(srcVof, 0) - baseivfabIrr.dataPtr(0);
        }
    }
  else
    {
      m_alphaBeta.resize(0);
    }

  const IntVect& smallendPhi = boxPhi.smallEnd();
  const IntVect& smallendLph = boxLph.smallEnd();

  IntVect ncellsPhi = boxPhi.size();
  IntVect ncellsLph = boxLph.size();
  m_ebstencil.resize(a_srcVofs.size());
  m_destTerms.resize(a_srcVofs.size());
  m_sourTerms.resize(a_srcVofs.size());

  m_cacheLph.resize(a_srcVofs.size());
  m_cachePhi.resize(a_srcVofs.size());
  //debugging hook
  //m_srcVoFs = a_srcVofs;
  for (int isrc = 0; isrc < a_srcVofs.size(); isrc++)
    {
      const VolIndex& srcVof = a_srcVofs[isrc];

      if (m_ebisBox.numVoFs(srcVof.gridIndex()) > 1)
        {//multi-valued (the dataPtr(0) is correct--that is where we start from)
          m_destTerms[isrc].offset = baseivfabLph.getIndex(srcVof, m_destVar) - baseivfabLph.dataPtr(0);
          m_destTerms[isrc].multiValued = true;
          m_sourTerms[isrc].offset = baseivfabPhi.getIndex(srcVof, m_destVar) - baseivfabPhi.dataPtr(0);
          m_sourTerms[isrc].multiValued = true;
        }
      else
        {//single-valued
          IntVect ivLph = srcVof.gridIndex()  - smallendLph;
          IntVect ivPhi = srcVof.gridIndex()  - smallendPhi;

          m_destTerms[isrc].offset = ivLph[0] + ivLph[1]*ncellsLph[0] ;
          m_sourTerms[isrc].offset = ivPhi[0] + ivPhi[1]*ncellsPhi[0] ;
#if CH_SPACEDIM==3
          m_destTerms[isrc].offset +=  ivLph[2]*ncellsLph[0]*ncellsLph[1];
          m_sourTerms[isrc].offset +=  ivPhi[2]*ncellsPhi[0]*ncellsPhi[1];
#endif
          //add in term due to variable number
#if CH_SPACEDIM==2
          m_destTerms[isrc].offset += m_destVar*ncellsLph[0]*ncellsLph[1];
          m_sourTerms[isrc].offset += m_destVar*ncellsPhi[0]*ncellsPhi[1];
#elif CH_SPACEDIM==3
          m_destTerms[isrc].offset += m_destVar*ncellsLph[0]*ncellsLph[1]*ncellsLph[2];
          m_sourTerms[isrc].offset += m_destVar*ncellsPhi[0]*ncellsPhi[1]*ncellsPhi[2];
#else
          bogus_spacedim();
#endif
          m_destTerms[isrc].multiValued = false;
        }
      const VoFStencil& sten = a_vofStencil(srcVof, 0);
      for (int isten = 0; isten < sten.size(); isten++)
        {
          const VolIndex stencilVof = sten.vof(isten);
          int srcVar = sten.variable(isten);
          stencilTerm stenEntry;
          stenEntry.weight = sten.weight(isten);
          //debugging hook
          //stenEntry.vof = stencilVof;
          if (m_ebisBox.numVoFs(stencilVof.gridIndex()) > 1)
            {//multi-valued(the dataPtr(0) is correct--that is where we start from)
              stenEntry.offset = baseivfabPhi.getIndex(stencilVof, srcVar) - baseivfabPhi.dataPtr(0);
              m_ebstencil[isrc].multi.push_back(stenEntry);
            }
          else
            {//single-valued
              IntVect ivPhi = stencilVof.gridIndex()  - smallendPhi;
              stenEntry.offset = ivPhi[0] + ivPhi[1]*ncellsPhi[0] ;
#if CH_SPACEDIM==3
              stenEntry.offset +=  ivPhi[2]*ncellsPhi[0]*ncellsPhi[1];
#endif
              //add in term due to variable number
#if CH_SPACEDIM==2
              stenEntry.offset += srcVar*ncellsPhi[0]*ncellsPhi[1];
#elif CH_SPACEDIM==3
              stenEntry.offset += srcVar*ncellsPhi[0]*ncellsPhi[1]*ncellsPhi[2];
#else
              bogus_spacedim();
#endif
              m_ebstencil[isrc].single.push_back(stenEntry);
            }
        }
    }
}

/**************/
/**************/
void
EBStencil::cachePhi(const EBCellFAB& a_phi, int a_ivar) const
{
//   CH_assert(m_cacheLph.size() == m_ebstencil.size());

  CH_assert(a_phi.getSingleValuedFAB().box()    == m_phiBox);
  const Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(a_ivar);
  const Real*  multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(a_ivar);

  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      if (m_sourTerms[isrc].multiValued)
        {
          m_cachePhi[isrc] = *(multiValuedPtrPhi + m_sourTerms[isrc].offset);
        }
      else
        {
          m_cachePhi[isrc] = *(singleValuedPtrPhi + m_sourTerms[isrc].offset);
        }
    }
}

void
EBStencil::cache(const EBCellFAB& a_lph, int a_ivar) const
{
//   CH_assert(m_cacheLph.size() == m_ebstencil.size());

  CH_assert(a_lph.getSingleValuedFAB().box()    == m_lphBox);
  const Real* singleValuedPtrLph =    a_lph.getSingleValuedFAB().dataPtr(a_ivar);
  const Real*  multiValuedPtrLph =    a_lph.getMultiValuedFAB().dataPtr(a_ivar);

  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      if (m_destTerms[isrc].multiValued)
        {
          m_cacheLph[isrc] = *(multiValuedPtrLph + m_destTerms[isrc].offset);
        }
      else
        {
          m_cacheLph[isrc] = *(singleValuedPtrLph + m_destTerms[isrc].offset);
        }
    }
}
/**************/
/**************/
void
EBStencil::uncachePhi(EBCellFAB& a_phi, int a_ivar) const
{
  CH_assert(a_phi.getSingleValuedFAB().box()    == m_phiBox);
  Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(a_ivar);
  Real*  multiValuedPtrPhi =    a_phi.getMultiValuedFAB().dataPtr(a_ivar);

  Real* phiPtr = NULL;

  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      if (m_sourTerms[isrc].multiValued)
        {
          phiPtr  = multiValuedPtrPhi + m_sourTerms[isrc].offset;
        }
      else
        {
          phiPtr  = singleValuedPtrPhi + m_sourTerms[isrc].offset;
        }

      *phiPtr = m_cachePhi[isrc];
    }
}
/**************/
void
EBStencil::uncache(EBCellFAB& a_lph, int a_ivar) const
{
  CH_assert(a_lph.getSingleValuedFAB().box()    == m_lphBox);
  Real* singleValuedPtrLph =    a_lph.getSingleValuedFAB().dataPtr(a_ivar);
  Real*  multiValuedPtrLph =    a_lph.getMultiValuedFAB().dataPtr(a_ivar);

  Real* lphiPtr = NULL;

  for (int isrc = 0; isrc < m_ebstencil.size(); isrc++)
    {
      if (m_destTerms[isrc].multiValued)
        {
          lphiPtr  = multiValuedPtrLph + m_destTerms[isrc].offset;
        }
      else
        {
          lphiPtr  = singleValuedPtrLph + m_destTerms[isrc].offset;
        }

      *lphiPtr = m_cacheLph[isrc];
    }
}

/**************/
/**************/
EBStencil::~EBStencil()
{
}
/**************/

#include "NamespaceFooter.H"
