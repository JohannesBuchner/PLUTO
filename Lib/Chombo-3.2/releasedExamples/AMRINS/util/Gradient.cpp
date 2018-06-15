#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Gradient.H"
#include "GradientF_F.H"
#include "EdgeToCell.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "Mask.H"

#define NUMGRADGROW 1
// ----------------------------------------------------------
void
Gradient::levelGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                           LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>* a_phiCrsePtr,
                           const Real a_dx, const int a_nRefCrse,
                           const Box& a_dProblem)
{
  ProblemDomain physdomain(a_dProblem);
  levelGradientMAC(a_edgeGrad, a_phi, a_phiCrsePtr,
                   a_dx, a_nRefCrse, physdomain);
}

// ----------------------------------------------------------
void
Gradient::levelGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                           LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>* a_phiCrsePtr,
                           const Real a_dx, const int a_nRefCrse,
                           const ProblemDomain& a_dProblem)
{
  // this function is just a wrapper around the other
  // gradient function with predefined QuadCFInterp
  // and LayoutData<IntVectSet>

  QuadCFInterp cfInterpCrse;
  const DisjointBoxLayout& grids = a_phi.getBoxes();

  if (a_phiCrsePtr != NULL)
    {
      // define coarse-fine interpolation operator
      const DisjointBoxLayout& crseGrids = a_phiCrsePtr->getBoxes();

      int ncomp = a_phi.nComp();
      cfInterpCrse.define(grids, &crseGrids, a_dx, a_nRefCrse, ncomp,
                          a_dProblem);
    }
  else
    {
      // just send in undefined cfInterpCrse
    }

  // centered-difference stencil implies that nGrow = 1
  int nGrow = NUMGRADGROW;
  LayoutData<IntVectSet> gridsIVS;
  createGridIVS(gridsIVS, grids, nGrow);

  levelGradientMAC(a_edgeGrad,a_phi,a_phiCrsePtr,a_dx,
                   gridsIVS,
                   cfInterpCrse);
}

// ----------------------------------------------------------
void
Gradient::levelGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                           LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>* a_phiCrsePtr,
                           const Real a_dx,
                           const LayoutData<IntVectSet>& a_gridIVS,
                           QuadCFInterp& a_cfInterpCrse)
{
  CH_assert(a_edgeGrad.nComp() >= a_phi.nComp());

  // do coarse-fine BC's
  if (a_phiCrsePtr != NULL)
    {
      const LevelData<FArrayBox>& crseData = *a_phiCrsePtr;
      CH_assert(a_phi.nComp() == crseData.nComp());
      CH_assert(a_cfInterpCrse.isDefined());

      a_cfInterpCrse.coarseFineInterp(a_phi, crseData);
    }

  DataIterator dit = a_phi.dataIterator();

  if (a_edgeGrad.nComp() == a_phi.nComp())
    {
      // loop over boxes and compute gradient
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& thisPhi = a_phi[dit];
          FluxBox& thisEdgeGrad = a_edgeGrad[dit];

          // loop over directions
          for (int dir = 0; dir<SpaceDim; dir++)
            {
              FArrayBox& edgeGradDirFab = thisEdgeGrad[dir];

              // only do this in interior of grid
              Box edgeBox(a_edgeGrad.getBoxes()[dit]);
              edgeBox.surroundingNodes(dir);

              // for one-component gradient, only do normal gradients
              int edgeDir = dir;

              singleBoxMacGrad(edgeGradDirFab,
                               thisPhi, 0, 0, thisPhi.nComp(),
                               edgeBox, a_dx, dir, edgeDir,
                               a_gridIVS[dit]);
            } // end loop over dir
        } // end loop over grids
    } // end if only doing normal directions
  else
    {
      // multicomponent gradPhi means that we also need to do
      // transverse directions.
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& thisPhi = a_phi[dit];
          FluxBox& thisEdgeGrad = a_edgeGrad[dit];

          // loop over edges
          for (int edgeDir=0; edgeDir<SpaceDim; edgeDir++)
            {
              FArrayBox& thisEdgeGradDirFab = thisEdgeGrad[edgeDir];
              const Box& edgeBox = thisEdgeGradDirFab.box();

              // loop over component directions in edgeGrad
              // (this will be direction of gradient)
              for (int dir=0; dir<thisEdgeGrad.nComp(); dir++)
                {
                  int phiComp = 0;
                  singleBoxMacGrad(thisEdgeGradDirFab, thisPhi,
                                   dir, phiComp, 1,
                                   edgeBox, a_dx, dir, edgeDir,
                                   a_gridIVS[dit]);
                }
            } // end loop over edgeDir
        } // end loop over grids
    } // end if we're also computing transverse direcctions
}

// ----------------------------------------------------------
void
Gradient::levelGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                           const LevelData<FArrayBox>& a_phi,
                           const Real a_dx)

{
  CH_assert(a_edgeGrad.nComp() >= a_phi.nComp());

  // build IVS to handle grid detals
  LayoutData<IntVectSet> gridIVS;
  int nGrow = NUMGRADGROW;
  createGridIVS(gridIVS, a_edgeGrad.getBoxes(), nGrow);

  DataIterator dit = a_phi.dataIterator();

  if (a_edgeGrad.nComp() == a_phi.nComp())
    {
      // loop over boxes and compute gradient
      for (dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& thisPhi = a_phi[dit];
          FluxBox& thisEdgeGrad = a_edgeGrad[dit];

          // loop over directions
          for (int dir = 0; dir<SpaceDim; dir++)
            {
              FArrayBox& edgeGradDirFab = thisEdgeGrad[dir];
              const Box& edgeBox = edgeGradDirFab.box();
              // for one-component gradient, only do normal gradients
              int edgeDir = dir;
              singleBoxMacGrad(edgeGradDirFab,
                               thisPhi,
                               0, 0, thisPhi.nComp(),
                               edgeBox, a_dx, dir, edgeDir,
                               gridIVS[dit]);

            } // end loop over dir
        } // end loop over grids
    } // end if only doing normal dir
  else if (a_edgeGrad.nComp() == SpaceDim*(a_phi.nComp()))
    {
    // multicomponent gradPhi means that we also need to do
    // transverse directions.
      DataIterator dit = a_phi.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& thisPhi = a_phi[dit];
          FluxBox& thisEdgeGrad = a_edgeGrad[dit];

          // loop over edges
          for (int edgeDir=0; edgeDir<SpaceDim; edgeDir++)
            {
              FArrayBox& thisEdgeGradDirFab = thisEdgeGrad[edgeDir];
              const Box& edgeBox = thisEdgeGradDirFab.box();

              // loop over component directions in edgeGrad
              // (this will be direction of gradient)
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  // loop over components in phi
                  for (int phiComp=0; phiComp<a_phi.nComp(); ++phiComp)
                    {
                      int gradComp = phiComp*SpaceDim + dir;

                      singleBoxMacGrad(thisEdgeGradDirFab, thisPhi,
                                       gradComp, phiComp, 1,
                                       edgeBox, a_dx, dir, edgeDir,
                                       gridIVS[dit]);

                    } // end loop over components
                } // end loop over gradient directions
            } // end loop over edge orientation
        } // end loop over grids
    }
  else
    {
      // bad number of components in either phi or gradPhi
      MayDay::Error("levelGradientMAC: bad number of components!");
    }
}

// ----------------------------------------------------------
void
Gradient::compGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                          LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>* a_phiCrsePtr,
                          const LevelData<FArrayBox>* a_phiFinePtr,
                          const Real a_dx, const int a_nRefCrse,
                          const int a_nRefFine, const Box& a_dProblem)
{
  ProblemDomain physdomain(a_dProblem);
  levelGradientMAC(a_edgeGrad, a_phi, a_phiCrsePtr,
                   a_dx, a_nRefCrse, physdomain);
}

// ----------------------------------------------------------
void
Gradient::compGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                          LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>* a_phiCrsePtr,
                          const LevelData<FArrayBox>* a_phiFinePtr,
                          const Real a_dx, const int a_nRefCrse,
                          const int a_nRefFine,
                          const ProblemDomain& a_dProblem)
{
  levelGradientMAC(a_edgeGrad, a_phi, a_phiCrsePtr,
                   a_dx, a_nRefCrse, a_dProblem);
}

// ----------------------------------------------------------
void
Gradient::compGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                          LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>* a_phiCrsePtr,
                          const LevelData<FArrayBox>* a_phiFinePtr,
                          const Real a_dx, const int a_nRefFine,
                          const LayoutData<IntVectSet>& a_gridIVS,
                          QuadCFInterp& a_cfInterpCrse)
{
  levelGradientMAC(a_edgeGrad, a_phi, a_phiCrsePtr, a_dx,
                   a_gridIVS,
                   a_cfInterpCrse);
}

// ----------------------------------------------------------
void
Gradient::levelGradientCC(LevelData<FArrayBox>& a_grad,
                          LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>* a_phiCrsePtr,
                          const Real a_dx,
                          const int a_nRefCrse,
                          const Box& a_dProblem)
{
  ProblemDomain physdomain(a_dProblem);
  levelGradientCC(a_grad, a_phi, a_phiCrsePtr, a_dx, a_nRefCrse,
                  physdomain);
}

// ----------------------------------------------------------
void
Gradient::levelGradientCC(LevelData<FArrayBox>& a_grad,
                          LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>* a_phiCrsePtr,
                          const Real a_dx,
                          const int a_nRefCrse,
                          const ProblemDomain& a_dProblem)
{
  QuadCFInterp cfInterpCrse;

  if (a_phiCrsePtr != NULL)
    {
      const DisjointBoxLayout& grids = a_phi.getBoxes();
      const DisjointBoxLayout& crseGrids = a_phiCrsePtr->getBoxes();
      int ncomp = a_phi.nComp();
      cfInterpCrse.define(grids, &crseGrids, a_dx, a_nRefCrse, ncomp,
                          a_dProblem);

    }
  else
    {
      // just send in undefined cfInterpCrse
    }

  LayoutData<IntVectSet> gridIVS;
  int nGrow = NUMGRADGROW;
  createGridIVS(gridIVS, a_phi.getBoxes(), nGrow);

  levelGradientCC(a_grad, a_phi, a_phiCrsePtr, a_dx,
                  gridIVS,
                  cfInterpCrse);
}

// ----------------------------------------------------------
void
Gradient::levelGradientCC(LevelData<FArrayBox>& a_grad,
                          LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>* a_phiCrsePtr,
                          const Real a_dx,
                          const LayoutData<IntVectSet>& a_gridIVS,
                          QuadCFInterp& a_cfInterpCrse)
{
  // first compute edge-centered gradient
  const DisjointBoxLayout& grids = a_grad.getBoxes();
  int ncompgrad = a_grad.nComp();
  ncompgrad = ncompgrad/SpaceDim;
  LevelData<FluxBox> edgeGrad(grids, ncompgrad);

  levelGradientMAC(edgeGrad,a_phi, a_phiCrsePtr, a_dx,
                   a_gridIVS,
                   a_cfInterpCrse);

  // now average edges->cells
  EdgeToCell(edgeGrad, a_grad);

  // that was pretty easy, wasn't it?
}

// ----------------------------------------------------------
void
Gradient::levelGradientCC(LevelData<FArrayBox>& a_grad,
                          const LevelData<FArrayBox>& a_phi,
                          const Real a_dx)
{
  // in this case, assume that all BC's have already been set
  // can just loop over grids and call gradCC subroutine
  // directly
  DataIterator dit = a_grad.dataIterator();
  const DisjointBoxLayout& grids = a_grad.getBoxes();

  for (dit.reset(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_GRADCC(CHF_FRA1(a_grad[dit],dir),
                      CHF_CONST_FRA1(a_phi[dit],0),
                      CHF_BOX(grids[dit]),
                      CHF_CONST_REAL(a_dx),
                      CHF_INT(dir));
        }
    }
}

// ----------------------------------------------------------
void
Gradient::compGradientCC(LevelData<FArrayBox>& a_grad,
                         LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>* a_phiCrsePtr,
                         const LevelData<FArrayBox>* a_phiFinePtr,
                         const Real a_dx,
                         const int a_nRefCrse, const int a_nRefFine,
                         const Box& a_dProblem)
{
  ProblemDomain physdomain(a_dProblem);

  compGradientCC(a_grad, a_phi, a_phiCrsePtr, a_phiFinePtr, a_dx,
                 a_nRefCrse, a_nRefFine, physdomain);
}

// ----------------------------------------------------------
void
Gradient::compGradientCC(LevelData<FArrayBox>& a_grad,
                         LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>* a_phiCrsePtr,
                         const LevelData<FArrayBox>* a_phiFinePtr,
                         const Real a_dx,
                         const int a_nRefCrse, const int a_nRefFine,
                         const ProblemDomain& a_dProblem)
{
  QuadCFInterp cfInterpCrse;

  if (a_phiCrsePtr != NULL)
    {
      const DisjointBoxLayout grids = a_phi.getBoxes();
      const DisjointBoxLayout crseGrids = a_phiCrsePtr->getBoxes();
      int ncomp = a_phi.nComp();
      cfInterpCrse.define(grids, &crseGrids, a_dx, a_nRefCrse, ncomp,
                          a_dProblem);
    }
  else
    {
      // send in undefined cfInterpCrse
    }

  LayoutData<IntVectSet> gridIVS;
  int nComp = NUMGRADGROW;
  createGridIVS(gridIVS, a_phi.getBoxes(), nComp);

  compGradientCC(a_grad,a_phi,a_phiCrsePtr,a_phiFinePtr,a_dx,
                 a_nRefFine,a_dProblem,
                 gridIVS,
                 cfInterpCrse);
}

// ----------------------------------------------------------
void
Gradient::compGradientCC(LevelData<FArrayBox>& a_grad,
                         LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>* a_phiCrsePtr,
                         const LevelData<FArrayBox>* a_phiFinePtr,
                         const Real a_dx, const int a_nRefFine,
                         const Box& a_dProblem,
                         const LayoutData<IntVectSet>& a_gridIVS,
                         QuadCFInterp& cfInterpCrse)
{
  ProblemDomain physdomain(a_dProblem);

  compGradientCC(a_grad, a_phi, a_phiCrsePtr, a_phiFinePtr, a_dx,
                 a_nRefFine, physdomain,
                 a_gridIVS,
                 cfInterpCrse);
}

// ----------------------------------------------------------
void
Gradient::compGradientCC(LevelData<FArrayBox>& a_grad,
                         LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>* a_phiCrsePtr,
                         const LevelData<FArrayBox>* a_phiFinePtr,
                         const Real a_dx, const int a_nRefFine,
                         const ProblemDomain& a_dProblem,
                         const LayoutData<IntVectSet>& a_gridIVS,
                         QuadCFInterp& cfInterpCrse)
{
  // first compute edge-centered gradient
  const DisjointBoxLayout& grids = a_grad.getBoxes();
  int ncompGrad = a_grad.nComp();
  ncompGrad = ncompGrad/SpaceDim;

  // nGhost is necessary to ensure that all edges for
  // extrapolation of face-centered gradient at
  // coarse-fine interpolation are present on the same
  // grid
  int nGhost = 1;
  LevelData<FluxBox> edgeGrad(grids, ncompGrad, nGhost*IntVect::Unit);

  compGradientMAC(edgeGrad,a_phi,a_phiCrsePtr,a_phiFinePtr,
                  a_dx, a_nRefFine,
                  a_gridIVS,
                  cfInterpCrse);

  if (a_phiFinePtr != NULL)
    {
      // do one-sided differencing at interfaces w/ finer levels
      // loop over this level grids
      const DisjointBoxLayout& fineGrids = a_phiFinePtr->getBoxes();

      // do exchange to ensure that all edges are filled with
      // appropriate data
      edgeGrad.exchange(edgeGrad.interval());

      // will need a mask for this one
      IntVect maskGrow(IntVect::Unit);
      maskGrow *= 2;
      LevelData<BaseFab<int> > masks(grids, 1, maskGrow);
      Mask thisMask;
      thisMask.buildMasks(masks, a_dProblem, grids, &fineGrids, a_nRefFine);

      // will need coarsened fine grids
      BoxLayout coarsenedFineGrids;

      coarsen(coarsenedFineGrids, fineGrids, a_nRefFine);

      // now loop over this level's boxes:
      DataIterator dit = edgeGrad.dataIterator();

      for (dit.reset(); dit.ok(); ++dit)
        {
          FluxBox& thisEdgeGrad = edgeGrad[dit];
          //const Box& thisGradBox = thisEdgeGrad.box();
          const Box& thisGradBox = edgeGrad.getBoxes()[dit];
          BaseFab<int>& thisMask = masks[dit];

          // now loop over (coarsened) fine boxes
          LayoutIterator litFine = coarsenedFineGrids.layoutIterator();

          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlapBox(thisGradBox);
              // grow fine-grid box by one to make sure we catch case
              // where coarse-fine and coarse-coarse interfaces meet.
              const Box& crseFineBox = coarsenedFineGrids.get(litFine());
              Box testBox(crseFineBox);
              testBox.grow(1);
              overlapBox &= testBox;

              if (!overlapBox.isEmpty())
                {
                  // fine grid overlays coarse grid, so we need to modify grad.

                  // loop over directions
                  for (int dir = 0; dir<SpaceDim; dir++)
                    {
                      FArrayBox& thisGradDir = thisEdgeGrad[dir];

                      // figure out which edge faces we need to correct
                      Box loEdgeBox = adjCellLo(crseFineBox,dir,1);
                      Box hiEdgeBox = adjCellHi(crseFineBox,dir,1);

                      loEdgeBox &= thisGradBox;
                      hiEdgeBox &= thisGradBox;

                      loEdgeBox.shiftHalf(dir, 1);
                      hiEdgeBox.shiftHalf(dir,-1);

                      // check to see if we want to do correction
                      // in this direction, both for low and hi
                      int do_lo = 1;
                      int do_hi = 1;
                      // don't do correction if fine-grid and crse-grid
                      // edges coincide
                      if (overlapBox.smallEnd(dir)<=thisGradBox.smallEnd(dir))
                        do_lo = 0;
                      if (overlapBox.bigEnd(dir) >= thisGradBox.bigEnd(dir))
                        do_hi = 0;

                      FORT_CRSEONESIDEGRAD(CHF_FRA1(thisGradDir,0),
                                           CHF_FIA1(thisMask,0),
                                           CHF_BOX(loEdgeBox),
                                           CHF_BOX(hiEdgeBox),
                                           CHF_INT(dir),
                                           CHF_INT(do_lo),
                                           CHF_INT(do_hi));
                    } // end loop over directions
                } // end if fine grid overlays this grid
            } // end loop over fine grids
        } // end loop over this level's grids
    }  // end if a finer level exists

  // now average to cells
  EdgeToCell(edgeGrad, a_grad);
}

// ----------------------------------------------------------
void
Gradient::compGradientCC(LevelData<FArrayBox>& a_grad,
                         const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>* a_phiFinePtr,
                         const Real a_dx, const int a_nRefFine,
                         const Box& a_dProblem)
{
  ProblemDomain physdomain(a_dProblem);

  compGradientCC(a_grad,  a_phi, a_phiFinePtr, a_dx,
                 a_nRefFine, physdomain);
}

// ----------------------------------------------------------
void
Gradient::compGradientCC(LevelData<FArrayBox>& a_grad,
                         const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>* a_phiFinePtr,
                         const Real a_dx, const int a_nRefFine,
                         const ProblemDomain& a_dProblem)
{
  // first compute edge-centered gradient
  const DisjointBoxLayout& grids = a_grad.getBoxes();
  int ncompGrad = a_grad.nComp();
  ncompGrad = ncompGrad/SpaceDim;

  LevelData<FluxBox> edgeGrad(grids, ncompGrad);

  levelGradientMAC(edgeGrad,a_phi, a_dx);

  if (a_phiFinePtr != NULL)
    {
      // do one-sided differencing at interfaces w/ finer levels
      // loop over this level grids
      const DisjointBoxLayout& fineGrids = a_phiFinePtr->getBoxes();

      // will need a mask for this one
      IntVect maskGrow(IntVect::Unit);
      maskGrow *= 2;
      LevelData<BaseFab<int> > masks(grids, 1, maskGrow);
      Mask thisMask;
      thisMask.buildMasks(masks, a_dProblem, grids, &fineGrids, a_nRefFine);

      // will need coarsened fine grids
      BoxLayout coarsenedFineGrids;

      coarsen(coarsenedFineGrids, fineGrids, a_nRefFine);

      // now loop over this level's boxes:
      DataIterator dit = edgeGrad.dataIterator();

      for (dit.reset(); dit.ok(); ++dit)
        {
          FluxBox& thisEdgeGrad = edgeGrad[dit];
          const Box& thisGradBox = thisEdgeGrad.box();
          BaseFab<int>& thisMask = masks[dit];

          // now loop over (coarsened) fine boxes
          LayoutIterator litFine = coarsenedFineGrids.layoutIterator();

          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlapBox(thisGradBox);
              overlapBox &= coarsenedFineGrids.get(litFine());

              if (!overlapBox.isEmpty())
                {
                  // fine grid overlays coarse grid, so we need to modify grad.

                  // loop over directions
                  for (int dir = 0; dir<SpaceDim; dir++)
                    {
                      FArrayBox& thisGradDir = thisEdgeGrad[dir];

                      Box loEdgeBox = bdryLo(overlapBox,dir,1);
                      Box hiEdgeBox = bdryHi(overlapBox,dir,1);

                      // check to see if we want to do correction
                      // in this direction, both for low and hi
                      int do_lo = 1;
                      int do_hi = 1;
                      // don't do correction if fine-grid and crse-grid edges
                      // coincide
                      if (overlapBox.smallEnd(dir)<=thisGradBox.smallEnd(dir))
                        do_lo = 0;
                      if (overlapBox.bigEnd(dir) >= thisGradBox.bigEnd(dir))
                        do_hi = 0;

                      FORT_CRSEONESIDEGRAD(CHF_FRA1(thisGradDir,0),
                                           CHF_FIA1(thisMask,0),
                                           CHF_BOX(loEdgeBox),
                                           CHF_BOX(hiEdgeBox),
                                           CHF_INT(dir),
                                           CHF_INT(do_lo),
                                           CHF_INT(do_hi));
                    } // end loop over directions
                } // end if fine grid overlays this grid
            } // end loop over fine grids
        } // end loop over this level's grids
    }  // end if a finer level exists

  // now average to cells
  EdgeToCell(edgeGrad, a_grad);
}

// utility function to do the actual computin' (to reduce code duplication)
void
Gradient::singleBoxMacGrad(FArrayBox& a_gradFab,
                           const FArrayBox& a_phiFab,
                           int a_gradComp,
                           int a_phiComp,
                           int a_numComp,
                           const Box& a_edgeBox,
                           Real a_dx,
                           int a_dir,
                           int a_edgeDir,
                           const IntVectSet& a_gridIVS)
{
  // corner boxes identify which cells need special
  // stencil for transverse derivatives
  // "bottom", "top", "left", and "right" are relative
  // to edgeDir=0, dir=1
  Box BLcorner(a_edgeBox);
  Box BRcorner(a_edgeBox);
  Box ULcorner(a_edgeBox);
  Box URcorner(a_edgeBox);

  if (a_edgeDir != a_dir)
    {
      // now loop over directions and set box edge extents
      for (int tempDir=0; tempDir<SpaceDim; ++tempDir)
        {
          if ((tempDir==a_edgeDir))
            {
              BLcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              ULcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              BRcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
              URcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
            }
          else if (tempDir == a_dir)
            {
              BLcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              BRcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              ULcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
              URcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
            }
        }
    }

  // loop over components
  int phiComp = a_phiComp;
  for (int comp=a_gradComp; comp< a_gradComp+a_numComp; comp++)
    {
      // now compute gradient in direction dir on this grid
      FORT_NEWMACGRAD(CHF_FRA1(a_gradFab,comp),
                      CHF_CONST_FRA1(a_phiFab,phiComp),
                      CHF_BOX(a_edgeBox),
                      CHF_CONST_REAL(a_dx),
                      CHF_INT(a_dir),
                      CHF_INT(a_edgeDir));

      // now do corner cells specially (replace gradient values computed
      // in fortran with one-sided gradients where necessary)
      if (a_dir != a_edgeDir)
        {
          Real factor = 0.25/a_dx;
          IntVect dirBasis = BASISV(a_dir);
          IntVect edgeBasis = BASISV(a_edgeDir);

          // "lower left"
          {
            BoxIterator bit(BLcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV -dirBasis - edgeBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis+edgeBasis, phiComp)
                                                         -a_phiFab(thisIV-dirBasis+edgeBasis, phiComp)));
                  }
              } // end loop over cells in BLcorner
          }  // end BLcorner context

          // "lower right"
          {
            BoxIterator bit(BRcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV -dirBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis-edgeBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis-edgeBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis-2*edgeBasis, phiComp)
                                                         -a_phiFab(thisIV-dirBasis-2*edgeBasis, phiComp)));
                  }
              } // end loop over cells in BRcorner
          }  // end BRcorner context

          // "upper left"
          {
            BoxIterator bit(ULcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV +dirBasis - edgeBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis+edgeBasis, phiComp)
                                                       -a_phiFab(thisIV-dirBasis+edgeBasis, phiComp)));
                  }
              } // end loop over cells in ULcorner
          }  // end ULcorner context

          // "upper right"
          {
            BoxIterator bit(URcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV +dirBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis-edgeBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis-edgeBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis-2*edgeBasis, phiComp)
                                                         -a_phiFab(thisIV-dirBasis-2*edgeBasis, phiComp)));
                  }
              } // end loop over cells in URcorner
          }  // end URcorner context
        }

      phiComp++;
    }
}

void
Gradient::createGridIVS(LayoutData<IntVectSet>& a_gridsIVS,
                        const DisjointBoxLayout& a_grids,
                        const int a_nGrow)
{
  a_gridsIVS.define(a_grids);
  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      createGridIVS(a_gridsIVS[dit],
                    a_grids,
                    a_grids[dit],
                    a_nGrow);
    }
}

// simple little function to turn DBL->IVS in region around
// a given box.
void
Gradient::createGridIVS(IntVectSet& a_gridIVS,
                        const DisjointBoxLayout& a_grids,
                        const Box& a_localBox,
                        const int a_nGrow)
{
  LayoutIterator lit = a_grids.layoutIterator();
  // if this isn't already empty, then make it empty
  if (!a_gridIVS.isEmpty())
    {
      a_gridIVS.makeEmpty();
    }

  // this is the neigborhood around thisBox in which we're interested.
  Box checkBox(a_localBox);
  checkBox.grow(a_nGrow);

  bool done = false;
    for (lit.begin(); lit.ok() && !done; ++lit)
    {
      const Box& thisBox = a_grids.get(lit());
      Box intersectBox(thisBox);
      intersectBox &= checkBox;
      if (!intersectBox.isEmpty())
        {
          a_gridIVS |= intersectBox;
        }

      // to make sure we do boundaries correctly, also include
      // adjacent cells
      // also check for whether we can stop due to ordering
      done = true;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          Box lobox = adjCellLo(thisBox, dir, 1);
          Box hibox = adjCellHi(thisBox, dir, 1);

          lobox &= checkBox;
          hibox &= checkBox;
          if (!lobox.isEmpty())
            {
              a_gridIVS |= lobox;
            }
          if (!hibox.isEmpty())
            {
              a_gridIVS |= hibox;
            }
          if (lobox.smallEnd() <= checkBox.bigEnd())
            {
              done = false;
            }
        }
    }
}

