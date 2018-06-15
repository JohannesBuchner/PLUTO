#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "EBLevelTGA.H"
#include "LevelTGAF_F.H"
#include "EBFluxFactory.H"
#include "NamespaceHeader.H"

void
EBLevelTGA::
setSourceGhostCells(LevelData<EBCellFAB>&    a_src,
                    const DisjointBoxLayout& a_grids,
                    int a_lev)
{
  int ncomp = a_src.nComp();
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid   = a_grids.get(dit());
      const Box& srcBox = a_src[dit()].box();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int iside = sign(sit());
              Box bc_box = adjCellBox(grid, idir, sit(), 1);

              for (int jdir = 0; jdir < SpaceDim; jdir++)
                {
                  //want corners too
                  if (jdir != idir)
                    {
                      bc_box.grow(jdir, 1);
                    }
                }

              //if fails might not have a ghost cell.
              bc_box &= m_eblg[a_lev].getDomain().domainBox();

              CH_assert(srcBox.contains(bc_box));
              if (grid.size(idir) >= 4)
                {
                  FORT_HORESGHOSTBC(CHF_FRA(a_src[dit()].getSingleValuedFAB()),
                                    CHF_BOX(bc_box),
                                    CHF_CONST_INT(idir),
                                    CHF_CONST_INT(iside),
                                    CHF_CONST_INT(ncomp));
                }
              else
                {
                  // valid region not wide enough to apply HOExtrap -- drop
                  // to linear extrap
                  FORT_RESGHOSTBC(CHF_FRA(a_src[dit()].getSingleValuedFAB()),
                                  CHF_BOX(bc_box),
                                  CHF_CONST_INT(idir),
                                  CHF_CONST_INT(iside),
                                  CHF_CONST_INT(ncomp));
                }
              IntVectSet ivs = m_eblg[a_lev].getEBISL()[dit()].getIrregIVS(bc_box);
              for (VoFIterator vofit(ivs, m_eblg[a_lev].getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
                {
                  for (int icomp = 0; icomp < ncomp; icomp++)
                    {
                      Real valNeigh = 0;
                      Vector<FaceIndex> faces = m_eblg[a_lev].getEBISL()[dit()].getFaces(vofit(), idir, flip(sit()));
                      for (int iface = 0; iface < faces.size(); iface++)
                        {
                          VolIndex vofNeigh = faces[iface].getVoF(flip(sit()));
                          valNeigh += a_src[dit()](vofNeigh, icomp);
                        }
                      if (faces.size() > 1) valNeigh /= faces.size();
                      a_src[dit()](vofit(), icomp) = valNeigh;
                    }
                }
            }
        }
    }
}


// in this function, we call the updateSoln function, compute
// the update, then subtract off the extra pieces to return the
// diffusive part of the update
void
EBLevelTGA::computeDiffusion(LevelData<EBCellFAB>&       a_diffusiveTerm,
                             LevelData<EBCellFAB>&       a_phiOld,
                             LevelData<EBCellFAB>&       a_src,
                             EBFluxRegister*             a_fineFluxRegPtr,
                             EBFluxRegister*             a_crseFluxRegPtr,
                             const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                             const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                             Real a_oldTime,
                             Real a_crseOldTime,
                             Real a_crseNewTime,
                             Real a_dt,
                             int  a_level,
                             bool a_zeroPhi,
                             bool a_rhsAlreadyKappaWeighted
                             )
{
  CH_assert(m_isEBLGSet);
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  //the IntVect::Zero is important here.
  EBFluxFactory fact(m_eblg[a_level].getEBISL());
  LevelData<EBFluxFAB>   flux(m_grids[a_level], ncomp,  IntVect::Unit, fact);

  BaseLevelTGA<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>::computeDiffusion
    (a_diffusiveTerm, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted);
}



void
EBLevelTGA::updateSoln(LevelData<EBCellFAB>&       a_phiNew,
                       LevelData<EBCellFAB>&       a_phiOld,
                       LevelData<EBCellFAB>&       a_src,
                       EBFluxRegister*             a_fineFluxRegPtr,
                       EBFluxRegister*             a_crseFluxRegPtr,
                       const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                       const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                       Real a_oldTime,
                       Real a_crseOldTime,
                       Real a_crseNewTime,
                       Real a_dt,
                       int  a_level,
                       bool a_zeroPhi,
                       bool a_rhsAlreadyKappaWeighted,
                       int  a_fluxStartComponent)
{
  CH_assert(m_isEBLGSet);
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  EBFluxFactory fact(m_eblg[a_level].getEBISL());
  //the IntVect::Zero is important here.
  LevelData<EBFluxFAB>   flux(m_grids[a_level], ncomp, IntVect::Zero, fact);

  BaseLevelTGA<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>::updateSoln
    (a_phiNew, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted, a_fluxStartComponent);
}

// in this function, we call the updateSoln function, compute
// the update, then subtract off the extra pieces to return the
// diffusive part of the update
void
EBLevelBackwardEuler::computeDiffusion(LevelData<EBCellFAB>&       a_diffusiveTerm,
                                       LevelData<EBCellFAB>&       a_phiOld,
                                       LevelData<EBCellFAB>&       a_src,
                                       EBFluxRegister*             a_fineFluxRegPtr,
                                       EBFluxRegister*             a_crseFluxRegPtr,
                                       const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                                       const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                                       Real a_oldTime,
                                       Real a_crseOldTime,
                                       Real a_crseNewTime,
                                       Real a_dt,
                                       int  a_level,
                                       bool a_zeroPhi,
                                       bool a_rhsAlreadyKappaWeighted)
{
  CH_assert(m_isEBLGSet);
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  //the IntVect::Zero is important here.
  EBFluxFactory fact(m_eblg[a_level].getEBISL());
  LevelData<EBFluxFAB>   flux(m_grids[a_level], ncomp,  IntVect::Unit, fact);

  BaseLevelBackwardEuler<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>::computeDiffusion
    (a_diffusiveTerm, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted);
}



void
EBLevelBackwardEuler::updateSoln(LevelData<EBCellFAB>&       a_phiNew,
                                 LevelData<EBCellFAB>&       a_phiOld,
                                 LevelData<EBCellFAB>&       a_src,
                                 EBFluxRegister*             a_fineFluxRegPtr,
                                 EBFluxRegister*             a_crseFluxRegPtr,
                                 const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                                 const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                                 Real a_oldTime,
                                 Real a_crseOldTime,
                                 Real a_crseNewTime,
                                 Real a_dt,
                                 int  a_level,
                                 bool a_zeroPhi,
                                 bool a_rhsAlreadyKappaWeighted,
                                 int  a_fluxStartComponent)
{
  CH_assert(m_isEBLGSet);
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  EBFluxFactory fact(m_eblg[a_level].getEBISL());
  //the IntVect::Zero is important here.
  LevelData<EBFluxFAB>   flux(m_grids[a_level], ncomp, IntVect::Zero, fact);

  BaseLevelBackwardEuler<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>::updateSoln
    (a_phiNew, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted, a_fluxStartComponent);
}


// in this function, we call the updateSoln function, compute
// the update, then subtract off the extra pieces to return the
// diffusive part of the update
void
EBLevelCrankNicolson::computeDiffusion(LevelData<EBCellFAB>&       a_diffusiveTerm,
                                       LevelData<EBCellFAB>&       a_phiOld,
                                       LevelData<EBCellFAB>&       a_src,
                                       EBFluxRegister*             a_fineFluxRegPtr,
                                       EBFluxRegister*             a_crseFluxRegPtr,
                                       const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                                       const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                                       Real a_oldTime,
                                       Real a_crseOldTime,
                                       Real a_crseNewTime,
                                       Real a_dt,
                                       int  a_level,
                                       bool a_zeroPhi,
                                       bool a_rhsAlreadyKappaWeighted)
{
  CH_assert(m_isEBLGSet);
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  //the IntVect::Zero is important here.
  EBFluxFactory fact(m_eblg[a_level].getEBISL());
  LevelData<EBFluxFAB>   flux(m_grids[a_level], ncomp,  IntVect::Unit, fact);

  BaseLevelCrankNicolson<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>::computeDiffusion
    (a_diffusiveTerm, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted);
}



void
EBLevelCrankNicolson::updateSoln(LevelData<EBCellFAB>&       a_phiNew,
                                 LevelData<EBCellFAB>&       a_phiOld,
                                 LevelData<EBCellFAB>&       a_src,
                                 EBFluxRegister*             a_fineFluxRegPtr,
                                 EBFluxRegister*             a_crseFluxRegPtr,
                                 const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                                 const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                                 Real a_oldTime,
                                 Real a_crseOldTime,
                                 Real a_crseNewTime,
                                 Real a_dt,
                                 int  a_level,
                                 bool a_zeroPhi,
                                 bool a_rhsAlreadyKappaWeighted,
                                 int  a_fluxStartComponent)
{
  CH_assert(m_isEBLGSet);
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  EBFluxFactory fact(m_eblg[a_level].getEBISL());
  //the IntVect::Zero is important here.
  LevelData<EBFluxFAB>   flux(m_grids[a_level], ncomp, IntVect::Zero, fact);

  BaseLevelCrankNicolson<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>::updateSoln
    (a_phiNew, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted, a_fluxStartComponent);
}

#include "NamespaceFooter.H"
