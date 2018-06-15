#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBStenVarCoef.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "NamespaceHeader.H"
/**************/
EBStenVarCoef::
EBStenVarCoef(const Vector<VolIndex>& a_srcVofs,
              const BaseIVFAB<VoFStencil>& a_vofStencil,
              const Box& a_box,
              const EBISBox& a_ebisBox,
              const IntVect& a_ghostVect,
              int a_varDest)
  : m_box(       a_box       ),
    m_ebisBox(   a_ebisBox   ),
    m_ghostVect( a_ghostVect ),
    m_destVar(   a_varDest   )
{
  CH_TIME("EBStenVarCoef::EBStenVarCoef");
  computeOffsets(a_srcVofs, a_vofStencil);
}
/***/
void
EBStenVarCoef::
apply(EBCellFAB&             a_lofphi,
      const EBCellFAB&       a_phi,
      const EBCellFAB&       a_alphaWeight,
      const Real&            a_alpha,
      const EBCellFAB&       a_betaWeight,
      const Real&            a_beta)
{
  CH_TIME("EBStenVarCoef::apply");

  CH_assert(a_lofphi.getSingleValuedFAB().box()      == m_grownBox);
  CH_assert(a_phi.getSingleValuedFAB().box()         == m_grownBox);
  CH_assert(a_alphaWeight.getSingleValuedFAB().box() == m_grownBox);
  CH_assert(a_betaWeight.getSingleValuedFAB().box()  == m_grownBox);

  const Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(0);
  const Real*  multiValuedPtrPhi =    a_phi. getMultiValuedFAB().dataPtr(0);

  const Real* singleValuedPtrAlp =    a_alphaWeight.getSingleValuedFAB().dataPtr(0);
  const Real*  multiValuedPtrAlp =    a_alphaWeight. getMultiValuedFAB().dataPtr(0);

  const Real* singleValuedPtrBet =     a_betaWeight.getSingleValuedFAB().dataPtr(0);
  const Real*  multiValuedPtrBet =     a_betaWeight. getMultiValuedFAB().dataPtr(0);

  Real*        multiValuedPtrLph = a_lofphi. getMultiValuedFAB().dataPtr(m_destVar);
  Real*       singleValuedPtrLph = a_lofphi.getSingleValuedFAB().dataPtr(m_destVar);

  //plo is different from phi because of destvar
  //phi is for stencil evaluation.  plo is for the local value at this variable
  const Real*   multiValuedPtrPlo =    a_phi. getMultiValuedFAB().dataPtr(m_destVar);
  const Real*  singleValuedPtrPlo =    a_phi.getSingleValuedFAB().dataPtr(m_destVar);

  for (int isrc = 0; isrc < m_stencil.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      Real*       lphPtr  = NULL;
      const Real* phiPtr  = NULL;
      const Real* alpPtr  = NULL;
      const Real* betPtr  = NULL;
      if (m_sourTerms[isrc].multiValued)
        {
          lphPtr  = multiValuedPtrLph + m_sourTerms[isrc].offset;
          phiPtr  = multiValuedPtrPlo + m_sourTerms[isrc].offset;
          alpPtr  = multiValuedPtrAlp + m_sourTerms[isrc].offset;
          betPtr  = multiValuedPtrBet + m_sourTerms[isrc].offset;
        }
      else
        {
          lphPtr  = singleValuedPtrLph + m_sourTerms[isrc].offset;
          phiPtr  = singleValuedPtrPlo + m_sourTerms[isrc].offset;
          alpPtr  = singleValuedPtrAlp + m_sourTerms[isrc].offset;
          betPtr  = singleValuedPtrBet + m_sourTerms[isrc].offset;
        }

      Real& lph = *lphPtr;
      const Real&         phi = *phiPtr;
      const Real& alphaWeight = *alpPtr;
      const Real&  betaWeight = *betPtr;

      lph =  0.;
      const varcsten_t& stenpt = m_stencil[isrc];
      //single-valued
      for (int isingle = 0; isingle < stenpt.single.size(); isingle++)
        {
          const int & offset = stenpt.single[isingle].offset;
          const Real& phiVal = *(singleValuedPtrPhi + offset);
          const Real& weight = stenpt.single[isingle].weight;
          lph += phiVal*weight;
        }
      //multi-valued
      for (int imulti = 0; imulti < stenpt.multi.size(); imulti++)
        {
          const int & offset = stenpt.multi[imulti].offset;
          const Real& phiVal = *(multiValuedPtrPhi + offset);
          const Real& weight = stenpt.multi[imulti].weight;
          lph += phiVal*weight;
        }
      //at this point lph holds divF.  add in identity terms
      //and multiply by factors
      lph = a_alpha*alphaWeight*phi + a_beta*betaWeight*lph;
    }
}
/********/
void
EBStenVarCoef::
relax(EBCellFAB&  a_phi,
      const EBCellFAB&  a_rhs,
      const EBCellFAB&  a_alphaWeight,
      const EBCellFAB&  a_betaWeight,
      const EBCellFAB&  a_lambda,
      Real a_alpha, Real a_beta) const
{
  CH_TIME("EBStenVarCoef::relax");
  CH_assert(a_rhs.getSingleValuedFAB().box()         == m_grownBox);
  CH_assert(a_phi.getSingleValuedFAB().box()         == m_grownBox);
  CH_assert(a_alphaWeight.getSingleValuedFAB().box() == m_grownBox);
  CH_assert(a_betaWeight.getSingleValuedFAB().box()  == m_grownBox);
  CH_assert(a_lambda.getSingleValuedFAB().box()      == m_grownBox);


  Real*       singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(0);
  const Real* singleValuedPtrRhs =    a_rhs.getSingleValuedFAB().dataPtr(m_destVar);
  Real*       multiValuedPtrPhi  =     a_phi.getMultiValuedFAB().dataPtr(0);
  const Real* multiValuedPtrRhs  =     a_rhs.getMultiValuedFAB().dataPtr(m_destVar);

  const Real* singleValuedPtrAlp =    a_alphaWeight.getSingleValuedFAB().dataPtr(0);
  const Real*  multiValuedPtrAlp =    a_alphaWeight. getMultiValuedFAB().dataPtr(0);

  const Real* singleValuedPtrBet =      a_betaWeight.getSingleValuedFAB().dataPtr(0);
  const Real*  multiValuedPtrBet =      a_betaWeight. getMultiValuedFAB().dataPtr(0);

  const Real* singleValuedPtrLam =         a_lambda.getSingleValuedFAB().dataPtr(0);
  const Real*  multiValuedPtrLam =         a_lambda. getMultiValuedFAB().dataPtr(0);

  //plo is different from phi because of destvar
  //phi is for stencil evaluation.  plo is for the local value at this variable
  Real*        multiValuedPtrPlo =    a_phi. getMultiValuedFAB().dataPtr(m_destVar);
  Real*       singleValuedPtrPlo =    a_phi.getSingleValuedFAB().dataPtr(m_destVar);

  for (int isrc = 0; isrc < m_stencil.size(); isrc++)
    {

      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      Real*       phiPtr  = NULL;
      const Real* rhsPtr  = NULL;
      const Real* alpPtr  = NULL;
      const Real* betPtr  = NULL;
      const Real* lamPtr  = NULL;
      if (m_sourTerms[isrc].multiValued)
        {
          rhsPtr  = multiValuedPtrRhs + m_sourTerms[isrc].offset;
          phiPtr  = multiValuedPtrPlo + m_sourTerms[isrc].offset;
          alpPtr  = multiValuedPtrAlp + m_sourTerms[isrc].offset;
          betPtr  = multiValuedPtrBet + m_sourTerms[isrc].offset;
          lamPtr  = multiValuedPtrLam + m_sourTerms[isrc].offset;
        }
      else
        {
          rhsPtr  = singleValuedPtrRhs + m_sourTerms[isrc].offset;
          phiPtr  = singleValuedPtrPlo + m_sourTerms[isrc].offset;
          alpPtr  = singleValuedPtrAlp + m_sourTerms[isrc].offset;
          betPtr  = singleValuedPtrBet + m_sourTerms[isrc].offset;
          lamPtr  = singleValuedPtrLam + m_sourTerms[isrc].offset;
        }

      Real&       phi         = *phiPtr;
      const Real& lambda      = *lamPtr;
      const Real& rhs         = *rhsPtr;
      const Real& alphaWeight = *alpPtr;
      const Real&  betaWeight = *betPtr;


      Real lph =  0.;
      const varcsten_t& stenpt = m_stencil[isrc];
      //single-valued
      for (int isingle = 0; isingle < stenpt.single.size(); isingle++)
        {
          const int & offset = stenpt.single[isingle].offset;
          const Real& phiVal = *(singleValuedPtrPhi + offset);
          const Real& weight = stenpt.single[isingle].weight;
          lph += phiVal*weight;
        }
      //multi-valued
      for (int imulti = 0; imulti < stenpt.multi.size(); imulti++)
        {
          const int & offset = stenpt.multi[imulti].offset;
          const Real& phiVal = *(multiValuedPtrPhi + offset);
          const Real& weight = stenpt.multi[imulti].weight;
          lph += phiVal*weight;
        }
      //at this point lph holds divF.  add in identity terms
      //and multiply by factors
      lph = a_alpha*alphaWeight*phi + a_beta*betaWeight*lph;

      phi = phi + lambda*(rhs - lph);
    }
}
/**************/
void
EBStenVarCoef::
computeOffsets(const Vector<VolIndex>&      a_srcVofs,
               const BaseIVFAB<VoFStencil>& a_vofStencil)
{
  CH_assert((SpaceDim == 2) || (SpaceDim == 3));
  m_grownBox = grow(m_box, m_ghostVect);

  const IntVectSet& ivsPhi = m_ebisBox.getMultiCells(m_grownBox);
  const EBGraph& ebgraph = m_ebisBox.getEBGraph();
  BaseIVFAB<Real> baseivfabPhi(ivsPhi, ebgraph, m_destVar+1);
  const IntVect& smallendPhi = m_grownBox.smallEnd();

  IntVect ncellsPhi = m_grownBox.size();
  m_stencil.resize(  a_srcVofs.size());
  m_sourTerms.resize(a_srcVofs.size());

  m_cache.resize(a_srcVofs.size());
  //debugging hook
  m_srcVoFs = a_srcVofs;
  for (int isrc = 0; isrc < a_srcVofs.size(); isrc++)
    {
      const VolIndex& srcVof = a_srcVofs[isrc];

      if (m_ebisBox.numVoFs(srcVof.gridIndex()) > 1)
        {//multi-valued (the dataPtr(0) is correct--that is where we start from)
          m_sourTerms[isrc].offset = baseivfabPhi.getIndex(srcVof, 0) - baseivfabPhi.dataPtr(0);
          m_sourTerms[isrc].multiValued = true;
        }
      else
        {//single-valued
          IntVect ivPhi = srcVof.gridIndex()  - smallendPhi;

          m_sourTerms[isrc].offset = ivPhi[0] + ivPhi[1]*ncellsPhi[0] ;
#if CH_SPACEDIM==3
          m_sourTerms[isrc].offset +=  ivPhi[2]*ncellsPhi[0]*ncellsPhi[1];
#endif
          m_sourTerms[isrc].multiValued = false;
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
              m_stencil[isrc].multi.push_back(stenEntry);
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
              m_stencil[isrc].single.push_back(stenEntry);
            }
        }
    }
}
/**************/
void
EBStenVarCoef::
cache(const EBCellFAB& a_phi, int a_ivar)
{
  CH_assert(a_phi.getSingleValuedFAB().box()   == m_grownBox);
  const Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(a_ivar);
  const Real*  multiValuedPtrPhi =     a_phi.getMultiValuedFAB().dataPtr(a_ivar);

  for (int isrc = 0; isrc < m_sourTerms.size(); isrc++)
    {
      //debugging hook
      //const VolIndex& srcVoF = m_srcVoFs[isrc];
      if (m_sourTerms[isrc].multiValued)
        {
          m_cache[isrc] = *(multiValuedPtrPhi + m_sourTerms[isrc].offset);
        }
      else
        {
          m_cache[isrc] = *(singleValuedPtrPhi + m_sourTerms[isrc].offset);
        }
    }
}
/**************/
void
EBStenVarCoef::
uncache(EBCellFAB& a_phi, int a_ivar) const
{
  CH_assert(a_phi.getSingleValuedFAB().box()    == m_grownBox);

  Real* singleValuedPtrPhi =    a_phi.getSingleValuedFAB().dataPtr(a_ivar);
  Real*  multiValuedPtrPhi =     a_phi.getMultiValuedFAB().dataPtr(a_ivar);
  Real* phiPtr = NULL;

  for (int isrc = 0; isrc < m_sourTerms.size(); isrc++)
    {
      if (m_sourTerms[isrc].multiValued)
        {
          phiPtr  = multiValuedPtrPhi + m_sourTerms[isrc].offset;
        }
      else
        {
          phiPtr  = singleValuedPtrPhi + m_sourTerms[isrc].offset;
        }

      *phiPtr = m_cache[isrc];
    }
}
/**************/

#include "NamespaceFooter.H"
