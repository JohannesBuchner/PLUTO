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

#include "LevelTGA.H"
#include "LevelTGAF_F.H"
#include "NamespaceHeader.H"

void
LevelTGA::
setSourceGhostCells(LevelData<FArrayBox>&    a_src,
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
              CH_assert(srcBox.contains(bc_box));
              if (grid.size(idir) >= 4)
                {
                  FORT_HORESGHOSTBC(CHF_FRA(a_src[dit()]),
                                    CHF_BOX(bc_box),
                                    CHF_CONST_INT(idir),
                                    CHF_CONST_INT(iside),
                                    CHF_CONST_INT(ncomp));
                }
              else
                {
                  // valid region not wide enough to apply HOExtrap -- drop
                  // to linear extrap
                  FORT_RESGHOSTBC(CHF_FRA(a_src[dit()]),
                                  CHF_BOX(bc_box),
                                  CHF_CONST_INT(idir),
                                  CHF_CONST_INT(iside),
                                  CHF_CONST_INT(ncomp));
                }

            }
        }
    }
}

// in this function, we call the updateSoln function, compute
// the update, then subtract off the extra pieces to return the
// diffusive part of the update
void
LevelTGA::computeDiffusion(LevelData<FArrayBox>& a_diffusiveTerm,
                           LevelData<FArrayBox>& a_phiOld,
                           LevelData<FArrayBox>& a_src,
                           LevelFluxRegister* a_fineFluxRegPtr,
                           LevelFluxRegister* a_crseFluxRegPtr,
                           const LevelData<FArrayBox>* a_crsePhiOldPtr,
                           const LevelData<FArrayBox>* a_crsePhiNewPtr,
                           Real a_oldTime,
                           Real a_crseOldTime,
                           Real a_crseNewTime,
                           Real a_dt,
                           int a_level)
{
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  //the IntVect::Zero is important here.
  LevelData<FluxBox>   flux(m_grids[a_level], ncomp,  IntVect::Zero);

  BaseLevelTGA<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>::computeDiffusion
    (a_diffusiveTerm, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level);
}

void
LevelTGA::updateSoln(LevelData<FArrayBox>& a_phiNew,
                     LevelData<FArrayBox>& a_phiOld,
                     LevelData<FArrayBox>& a_src,
                     LevelFluxRegister* a_fineFluxRegPtr,
                     LevelFluxRegister* a_crseFluxRegPtr,
                     const LevelData<FArrayBox>* a_crsePhiOldPtr,
                     const LevelData<FArrayBox>* a_crsePhiNewPtr,
                     Real a_oldTime,
                     Real a_crseOldTime,
                     Real a_crseNewTime,
                     Real a_dt,
                     int a_level,
                     bool a_zeroPhi,
                     int a_fluxStartComponent)
{
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  //the IntVect::Zero is important here.
  LevelData<FluxBox>   flux(m_grids[a_level], ncomp, IntVect::Zero);

  BaseLevelTGA<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>::updateSoln
    (a_phiNew, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_fluxStartComponent);
}

void
LevelBackwardEuler::
setSourceGhostCells(LevelData<FArrayBox>&    a_src,
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
              CH_assert(srcBox.contains(bc_box));
              if (grid.size(idir) >= 4)
                {
                  FORT_HORESGHOSTBC(CHF_FRA(a_src[dit()]),
                                    CHF_BOX(bc_box),
                                    CHF_CONST_INT(idir),
                                    CHF_CONST_INT(iside),
                                    CHF_CONST_INT(ncomp));
                }
              else
                {
                  // valid region not wide enough to apply HOExtrap -- drop
                  // to linear extrap
                  FORT_RESGHOSTBC(CHF_FRA(a_src[dit()]),
                                  CHF_BOX(bc_box),
                                  CHF_CONST_INT(idir),
                                  CHF_CONST_INT(iside),
                                  CHF_CONST_INT(ncomp));
                }

            }
        }
    }
}

// in this function, we call the updateSoln function, compute
// the update, then subtract off the extra pieces to return the
// diffusive part of the update
void
LevelBackwardEuler::computeDiffusion(LevelData<FArrayBox>& a_diffusiveTerm,
                           LevelData<FArrayBox>& a_phiOld,
                           LevelData<FArrayBox>& a_src,
                           LevelFluxRegister* a_fineFluxRegPtr,
                           LevelFluxRegister* a_crseFluxRegPtr,
                           const LevelData<FArrayBox>* a_crsePhiOldPtr,
                           const LevelData<FArrayBox>* a_crsePhiNewPtr,
                           Real a_oldTime,
                           Real a_crseOldTime,
                           Real a_crseNewTime,
                           Real a_dt,
                           int a_level)
{
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  //the IntVect::Zero is important here.
  LevelData<FluxBox>   flux(m_grids[a_level], ncomp,  IntVect::Zero);

  BaseLevelBackwardEuler<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>::computeDiffusion
    (a_diffusiveTerm, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level);
}

void
LevelBackwardEuler::updateSoln(LevelData<FArrayBox>& a_phiNew,
                     LevelData<FArrayBox>& a_phiOld,
                     LevelData<FArrayBox>& a_src,
                     LevelFluxRegister* a_fineFluxRegPtr,
                     LevelFluxRegister* a_crseFluxRegPtr,
                     const LevelData<FArrayBox>* a_crsePhiOldPtr,
                     const LevelData<FArrayBox>* a_crsePhiNewPtr,
                     Real a_oldTime,
                     Real a_crseOldTime,
                     Real a_crseNewTime,
                     Real a_dt,
                     int a_level,
                     bool a_zeroPhi,
                     int a_fluxStartComponent)
{
  // first compute updated solution
  int ncomp = a_phiOld.nComp();
  //the IntVect::Zero is important here.
  LevelData<FluxBox>   flux(m_grids[a_level], ncomp, IntVect::Zero);

  BaseLevelBackwardEuler<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>::updateSoln
    (a_phiNew, a_phiOld, a_src, flux,
     a_fineFluxRegPtr, a_crseFluxRegPtr,
     a_crsePhiOldPtr,  a_crsePhiNewPtr,
     a_oldTime, a_crseOldTime, a_crseNewTime, a_dt,
     a_level, a_zeroPhi, a_fluxStartComponent);
}

#include "NamespaceFooter.H"
