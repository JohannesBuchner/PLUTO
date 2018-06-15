#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoHiSide.H"

#include "PatchGodunov.H"
#include "NamespaceHeader.H"

// Flag everything as not defined or set
PatchGodunov::PatchGodunov()
{
  m_gdnvPhysics      = NULL;
  m_isDefined        = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet  = false;
}

PatchGodunov::~PatchGodunov()
{
  // Delete the initial/boundary condition object - if it exists
  if (m_gdnvPhysics != NULL)
    {
      delete m_gdnvPhysics;
      m_gdnvPhysics  = NULL;
    }
}

// Define this object and the boundary condition object
void PatchGodunov::define(const ProblemDomain&           a_domain,
                          const Real&                    a_dx,
                          const GodunovPhysics* const    a_physPtr,
                          const int&                     a_normalPredOrder,
                          const bool&                    a_useFourthOrderSlopes,
                          const bool&                    a_usePrimLimiting,
                          const bool&                    a_useCharLimiting,
                          const bool&                    a_useFlattening,
                          const bool&                    a_useArtificialViscosity,
                          const Real&                    a_artificialViscosity)
{
  // Store the domain and grid spacing
  m_domain = a_domain;
  m_dx     = a_dx;

  if (m_gdnvPhysics != NULL)
    {
      delete m_gdnvPhysics;
      m_gdnvPhysics  = NULL;
    }

  m_gdnvPhysics = a_physPtr->new_godunovPhysics();
  m_gdnvPhysics->define(m_domain, m_dx);

  m_util.define(m_domain, m_dx);

  m_normalPredOrder = a_normalPredOrder;

  // Slope flattening is only valid when 4th order slopes are being computed
  // and limited in some way
  CH_assert(!a_useFlattening || (a_useFourthOrderSlopes &&
                              (a_usePrimLimiting || a_useCharLimiting)));

  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_usePrimLimiting      = a_usePrimLimiting;
  m_useCharLimiting      = a_useCharLimiting;
  m_useFlattening        = a_useFlattening;

  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity    = a_artificialViscosity;

  m_isDefined = true;
}

// Set the current box.
void PatchGodunov::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime      = a_currentTime;
  m_isCurrentTimeSet = true;
}

// Set the current box.
void PatchGodunov::setCurrentBox(const Box& a_currentBox)
{
  m_currentBox      = a_currentBox;
  m_isCurrentBoxSet = true;

  m_gdnvPhysics->setCurrentBox(a_currentBox);
}

void PatchGodunov::updateState(FArrayBox&       a_U,
                               FluxBox&         a_F,
                               Real&            a_maxWaveSpeed,
                               const FArrayBox& a_S,
                               const Real&      a_dt,
                               const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_box == m_currentBox);

  int numPrim = m_gdnvPhysics->numPrimitives();
  int numFlux = m_gdnvPhysics->numFluxes();

  FluxBox whalf(a_box,numPrim);
  whalf.setVal(0.0);

  a_F.resize(a_box,numFlux);
  a_F.setVal(0.0);

  computeWHalf(whalf, a_U, a_S, a_dt, a_box);

  FArrayBox dU(a_U.box(),a_U.nComp());
  computeUpdate(dU, a_F, a_U, whalf, a_dt, a_box);

  a_U += dU;

  // Get and return the maximum wave speed on this patch/grid
  a_maxWaveSpeed = m_gdnvPhysics->getMaxWaveSpeed(a_U, m_currentBox);
}

void PatchGodunov::updateState(FArrayBox&       a_U,
                               FluxBox&         a_F,
                               FluxBox&         a_wHalf,
                               Real&            a_maxWaveSpeed,
                               const FArrayBox& a_S,
                               const Real&      a_dt,
                               const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_box == m_currentBox);

  //  int numPrim = m_gdnvPhysics->numPrimitives();
  int numFlux = m_gdnvPhysics->numFluxes();

  a_wHalf.setVal(0.0);

  a_F.resize(a_box,numFlux);
  a_F.setVal(0.0);

  computeWHalf(a_wHalf, a_U, a_S, a_dt, a_box);

  FArrayBox dU(a_U.box(),a_U.nComp());
  computeUpdate(dU, a_F, a_U, a_wHalf, a_dt, a_box);

  a_U += dU;

  // Get and return the maximum wave speed on this patch/grid
  a_maxWaveSpeed = m_gdnvPhysics->getMaxWaveSpeed(a_U, m_currentBox);
}

// Compute the time-centered values of the primitive variable on the
// interior cell faces of the cell-centered box a_box.
void PatchGodunov::computeWHalf(FluxBox&         a_WHalf,
                                const FArrayBox& a_U,
                                const FArrayBox& a_S,
                                const Real&      a_dt,
                                const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_box == m_currentBox);

  // Get the number of various variables
  int numPrim = m_gdnvPhysics->numPrimitives();

  // The current box of valid cells
  Box curBox = m_currentBox;

  // Boxes for face-centered state - used for the riemann() and
  // artificialViscosity() calls
  Box faceBox[SpaceDim];

  // Boxes for face-centered fluxes
  Box fluxBox[SpaceDim];

  // Boxes for cell-centered state - used for the updatePrim() calls
  Box ccBox[SpaceDim];

  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      // faceBox[dir1] is face-centered in direction "dir1",
      // is valid "curBox" grown by 1 in all directions except "dir1",
      // and stays one cell away, in "dir1", from the domain boundary.
      faceBox[dir1] = curBox; // valid cell-centered box
      faceBox[dir1].grow(1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1, -1);
      faceBox[dir1].surroundingNodes(dir1);

      // ccBox[dir1] is cell-centered,
      // is valid "curBox" grown by 1 in all directions except "dir1",
      // but only within the domain, "m_domain".
      ccBox[dir1] = curBox;
      ccBox[dir1].grow(1);
      ccBox[dir1].grow(dir1, -1);
      ccBox[dir1] &= m_domain;

      // fluxBox[dir1] is face-centered in direction "dir1",
      // consisting of all those faces of cells of "ccBox[dir1]".
      fluxBox[dir1] = ccBox[dir1];
      fluxBox[dir1].surroundingNodes(dir1);

      // The difference between faceBox[dir1] and fluxBox[dir1] is:
      // if curBox abuts the boundary of the domain in direction "dir1",
      // then fluxBox[dir1] contains faces along that boundary,
      // but faceBox[dir1] does not.
    }

  // cell-centered box:  but restrict it to within domain
  Box WBox = a_U.box();
  WBox &= m_domain;

  // Primitive variables
  FArrayBox W(WBox,numPrim);

  // Calculate the primitive variables W from the conserved variables a_U,
  // on cell-centered WBox.
  m_gdnvPhysics->consToPrim(W, a_U, WBox);

  // slopeBox is the cell-centered box where slopes will be needed:
  // it is one larger than the final update box "curBox" of valid cells,
  // but within domain.
  // On slopeBox we define the FABs flattening, WMinus, and WPlus.
  Box slopeBox = curBox;
  slopeBox.grow(1);
  slopeBox &= m_domain;

  // Compute flattening once for all slopes if needed
  FArrayBox flattening(slopeBox, 1); // cell-centered
  if (m_useFlattening)
    {
      Interval velInt    = m_gdnvPhysics->velocityInterval();
      int pressureIndex  = m_gdnvPhysics->pressureIndex();
      Real smallPressure = m_gdnvPhysics->smallPressure();
      int bulkIndex      = m_gdnvPhysics->bulkModulusIndex();

      m_util.computeFlattening(flattening,
                               W,
                               velInt,
                               pressureIndex,
                               smallPressure,
                               bulkIndex,
                               slopeBox);
    }

  // Intermediate, extrapolated primitive variables
  FArrayBox WMinus[SpaceDim];
  FArrayBox WPlus [SpaceDim];

  // Initial fluxes
  FArrayBox WHalf1[SpaceDim];

  // Source term computed from current state
  FArrayBox localSource;

  // If the source term is valid, make a local copy, increment it, and scale
  // it by 1/2 the timestep
  if (!a_S.box().isEmpty())
    {
      localSource.define(a_S.box(), a_S.nComp());
      localSource.copy(a_S);

      m_gdnvPhysics->incrementSource(localSource,W,slopeBox);

      localSource *= 0.5 * a_dt;
    }

  // Compute initial fluxes
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      // Size the intermediate, extrapolated primitive variables
      WMinus[dir1].resize(slopeBox,numPrim); // cell-centered
      WPlus [dir1].resize(slopeBox,numPrim); // cell-centered

      // Compute predictor step to obtain extrapolated primitive variables
      if (m_normalPredOrder == 0)
        {
          CTUNormalPred(WMinus[dir1],WPlus[dir1],a_dt,m_dx,W,
                        flattening,dir1,slopeBox);
        }
      else if (m_normalPredOrder == 1)
        {
          PLMNormalPred(WMinus[dir1],WPlus[dir1],a_dt,m_dx,W,
                        flattening,dir1,slopeBox);
        }
      else if (m_normalPredOrder == 2)
        {
          PPMNormalPred(WMinus[dir1],WPlus[dir1],a_dt,m_dx,W,
                        flattening,dir1,slopeBox);
        }
      else
        {
          MayDay::Error("PatchGodunov::computeWHalf: Normal predictor order must be 1 (PLM) or 2 (PPM)");
        }

      // If the source term is valid add it to the primitive quantities
      if (!localSource.box().isEmpty())
        {
          WMinus[dir1] += localSource;
          WPlus [dir1] += localSource;
        }

      // Solve the Riemann problem
      WHalf1[dir1].resize(fluxBox[dir1],numPrim); // face-centered
      m_gdnvPhysics->riemann(WHalf1[dir1],WPlus[dir1],WMinus[dir1],W,
                             m_currentTime,dir1,faceBox[dir1]);
    }

#if (CH_SPACEDIM == 3)
  // In 3D, compute some additional intermediate fluxes
  //
  // NOTE:  The diagonal entries of this array of fluxes are not
  // used and will not be defined.
  FArrayBox WHalf2[SpaceDim][SpaceDim];

  // Compute the intermediate, corrected fluxes in each direction
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      // Correct fluxes using fluxes from a different direction
      for (int dir2 = 0; dir2 < SpaceDim; dir2++)
        {
          // A different direction has been found
          if (dir2 != dir1)
            {
              // Temporary primitive variables
              FArrayBox WTempMinus(WMinus[dir1].box(),numPrim);
              FArrayBox WTempPlus (WPlus [dir1].box(),numPrim);
              FArrayBox AdWdx(WPlus[dir1].box(),numPrim);

              // preventing uninitialized memory reads which cause
              // FPE's on some machines
              // AdWdx.setVal(666.666);

              // Copy data for in place modification
              WTempMinus.copy(WMinus[dir1]);
              WTempPlus .copy(WPlus [dir1]);

              // Update the current, extrapolated primitive variable using a flux
              // in a different direction

              m_gdnvPhysics->quasilinearUpdate(AdWdx,
                                               WHalf1[dir2],
                                               WTempMinus,
                                               -(1.0/3.0) * a_dt / m_dx,
                                               dir2,
                                               ccBox[dir2]);
              WTempMinus += AdWdx;

              m_gdnvPhysics->quasilinearUpdate(AdWdx,
                                               WHalf1[dir2],
                                               WTempPlus,
                                               -(1.0/3.0) * a_dt / m_dx,
                                               dir2,
                                               ccBox[dir2]);
              WTempPlus  += AdWdx;

        // Solve the Riemann problem.

        WHalf2[dir1][dir2].resize(fluxBox[dir1],numPrim);
        m_gdnvPhysics->riemann(WHalf2[dir1][dir2],
                               WTempPlus,
                               WTempMinus,
                               W,
                               m_currentTime,
                               dir1,
                               faceBox[dir1]);
            }
        }
    }
#endif

  // faceBox and fluxBox are now a bit smaller for the final corrections
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      faceBox[dir1] = curBox;
      faceBox[dir1].grow(dir1,1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);

      fluxBox[dir1] = curBox;
      fluxBox[dir1].surroundingNodes(dir1);
    }

  // Do the final corrections to the fluxes
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      // Correct the flux using fluxes in the remaining direction(s)
      for (int dir2 = 0; dir2 < SpaceDim; dir2++)
        {
          // A different direction has been found
          if (dir2 != dir1)
            {
#if (CH_SPACEDIM == 2)
              // In 2D, the current primitive state is updated by a flux in
              // the other direction
              FArrayBox AdWdx(WPlus[dir1].box(),numPrim);

              m_gdnvPhysics->quasilinearUpdate(AdWdx,
                                               WHalf1[dir2],
                                               WMinus[dir1],
                                               -(1.0/2.0) * a_dt / m_dx,
                                               dir2,
                                               ccBox[dir2]);
              WMinus[dir1] += AdWdx;

              m_gdnvPhysics->quasilinearUpdate(AdWdx,
                                               WHalf1[dir2],
                                               WPlus [dir1],
                                               -(1.0/2.0) * a_dt / m_dx,
                                               dir2,
                                               ccBox[dir2]);
              WPlus[dir1] += AdWdx;

#elif (CH_SPACEDIM == 3)
              // In 3D, find a direction different from the two above
              int dir3 = 3 - dir1 - dir2;

              // Update the conservative state using both corrected fluxes in
              // the other two directions
              FArrayBox AdWdx(WPlus[dir1].box(),numPrim);

              m_gdnvPhysics->quasilinearUpdate(AdWdx,
                                               WHalf2[dir2][dir3],
                                               WMinus[dir1],
                                               -(1.0/2.0) * a_dt / m_dx,
                                               dir2,
                                               ccBox[dir2]);
              WMinus[dir1].plus(AdWdx,0,0,numPrim);

              m_gdnvPhysics->quasilinearUpdate(AdWdx,
                                               WHalf2[dir2][dir3],
                                               WPlus [dir1],
                                               -(1.0/2.0) * a_dt / m_dx,
                                               dir2,
                                               ccBox[dir2]);
              WPlus [dir1].plus(AdWdx,0,0,numPrim);
#else
              // Only 2D and 3D should be possible
              MayDay::Error("PatchGodunov::computeWHalf: CH_SPACEDIM not 2 or 3!");
#endif
            }
        }

      // Solve the Riemann problem to obtain time-centered face values.
      // be returned
      FArrayBox& WHalf = a_WHalf[dir1];
      m_gdnvPhysics->riemann(WHalf,
                             WPlus[dir1],
                             WMinus[dir1],
                             W,
                             m_currentTime,
                             dir1,
                             faceBox[dir1]);
    }
}

void PatchGodunov::computeUpdate(FArrayBox&       a_dU,
                                 FluxBox&         a_F,
                                 const FArrayBox& a_U,
                                 const FluxBox&   a_WHalf,
                                 const Real&      a_dt,
                                 const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_box == m_currentBox);

  m_gdnvPhysics->computeUpdate(a_dU,
                               a_F,
                               a_U,
                               a_WHalf,
                               m_useArtificialViscosity,
                               m_artificialViscosity,
                               m_currentTime,
                               m_dx,
                               a_dt,
                               a_box);
}

void PatchGodunov::CTUNormalPred(FArrayBox&       a_WMinus,
                                 FArrayBox&       a_WPlus,
                                 const Real&      a_dt,
                                 const Real&      a_dx,
                                 const FArrayBox& a_W,
                                 const FArrayBox& a_flat,
                                 const int&       a_dir,
                                 const Box&       a_box)
{
  // for CTU, increments are 0 -- straight copy from cells to faces
  a_WMinus.copy(a_W);
  a_WPlus.copy(a_W);

  // Apply a physics-dependent post-normal predictor step:
  // For example:
  //   - adjust/bound delta's so constraints on extrapolated primitive
  //     quantities are enforced (density and pressure > 0).
  //   - compute source terms that depend on the spatially varying
  //     coefficients.
  m_gdnvPhysics->postNormalPred(a_WMinus,a_WPlus,a_W,a_dt,a_dx,a_dir,a_box);
}


void PatchGodunov::PLMNormalPred(FArrayBox&       a_WMinus,
                                 FArrayBox&       a_WPlus,
                                 const Real&      a_dt,
                                 const Real&      a_dx,
                                 const FArrayBox& a_W,
                                 const FArrayBox& a_flat,
                                 const int&       a_dir,
                                 const Box&       a_box)
{
  int numprim = m_gdnvPhysics->numPrimitives();

  // This will hold 2nd or 4th order slopes
  FArrayBox dW(a_box,numprim);

  if (m_useFourthOrderSlopes)
    {
      // 2nd order slopes need to be computed over a larger box to accommodate
      // the 4th order slope computation
      Box boxVL = a_box;
      boxVL.grow(a_dir,1);
      boxVL &= m_domain;

      // Compute 2nd order (van Leer) slopes
      FArrayBox dWvL(boxVL, numprim);
      m_util.vanLeerSlopes(dWvL,a_W,numprim,
                           m_useCharLimiting || m_usePrimLimiting,
                           a_dir,boxVL);
      m_gdnvPhysics->getPhysIBC()->setBdrySlopes(dWvL,a_W,a_dir,m_currentTime);

      // Compute 4th order slopes, without limiting.
      m_util.fourthOrderSlopes(dW,a_W,dWvL, numprim, a_dir,a_box);
    }
  else
    {
      // Compute 2nd order (van Leer) slopes
      m_util.vanLeerSlopes(dW,a_W,numprim,
                           m_useCharLimiting || m_usePrimLimiting,
                           a_dir,a_box);
      m_gdnvPhysics->getPhysIBC()->setBdrySlopes(dW,a_W,a_dir,m_currentTime);
    }

  // To save on storage, we use the input values as temporaries for the
  // delta's
  a_WMinus.setVal(0.0);
  a_WPlus .setVal(0.0);

  if (m_useCharLimiting || m_useFourthOrderSlopes)
    {
      // Compute one-sided differences as inputs for limiting.
      m_util.oneSidedDifferences(a_WMinus,a_WPlus,a_W,a_dir,a_box);
    }

  FArrayBox lambda(a_box, numprim);
  m_gdnvPhysics->charValues(lambda, a_W, a_dir,a_box);

  if (m_useCharLimiting && m_usePrimLimiting)
    {
      MayDay::Error("PatchGodunov::PLMNormalPred: Attempt to limit slopes in primitive AND characteristic coordinates - not implemented");
    }

  // Apply limiter on characteristic or primitive variables. Either
  // way, must end up with characteristic variables to pass to to the
  // normal predictor utility.  Currently, cannot do both.

  // If doing characteristic limiting then transform before limiting
  if (m_useCharLimiting)
    {
      // Transform from primitive to characteristic variables
      m_gdnvPhysics->charAnalysis(a_WMinus,a_W,a_dir,a_box);
      m_gdnvPhysics->charAnalysis(a_WPlus ,a_W,a_dir,a_box);
      m_gdnvPhysics->charAnalysis(dW      ,a_W,a_dir,a_box);
    }

  if (m_useCharLimiting || m_usePrimLimiting)
    {
      // Limiting is already done for 2nd order slopes in primitive variables
      // so don't do it again
      if (m_useCharLimiting || m_useFourthOrderSlopes)
        {
          // Do slope limiting
          m_util.slopeLimiter(dW,a_WMinus,a_WPlus,numprim,a_box);
        }

      // Do slope flattening
      if (m_useFlattening)
        {
          m_util.applyFlattening(dW,a_flat,a_box);
        }
    }

  // If not doing characteristic limiting then transform after any limiting
  if (!m_useCharLimiting)
    {
      // Transform from primitive to characteristic variables
      m_gdnvPhysics->charAnalysis(dW,a_W,a_dir,a_box);
    }

  // To the normal prediction in characteristic variables
  m_util.PLMNormalPred(a_WMinus,a_WPlus,dW,lambda,a_dt / a_dx,a_box);

  // Construct the increments to the primitive variables
  m_gdnvPhysics->charSynthesis(a_WMinus,a_W,a_dir,a_box);
  m_gdnvPhysics->charSynthesis(a_WPlus ,a_W,a_dir,a_box);

  // Apply a physics-dependent post-normal predictor step:
  // For example:
  //   - adjust/bound delta's so constraints on extrapolated primitive
  //     quantities are enforced (density and pressure > 0).
  //   - compute source terms that depend on the spatially varying
  //     coefficients.
  m_gdnvPhysics->postNormalPred(a_WMinus,a_WPlus,a_W,a_dt,a_dx,a_dir,a_box);

  // Compute the state from the increments
  a_WMinus += a_W;
  a_WPlus  += a_W;
}

void PatchGodunov::PPMNormalPred(FArrayBox&       a_WMinus,
                                 FArrayBox&       a_WPlus,
                                 const Real&      a_dt,
                                 const Real&      a_dx,
                                 const FArrayBox& a_W,
                                 const FArrayBox& a_flat,
                                 const int&       a_dir,
                                 const Box&       a_box)
{
  int numprim = m_gdnvPhysics->numPrimitives();

  Box faceBox = a_box;
  // added by petermc, 22 Sep 2008:
  // for 4th order, need extra faces in all the directions
  if (m_highOrderLimiter) faceBox.grow(1);
  faceBox.surroundingNodes(a_dir);
  FArrayBox WFace(faceBox,numprim);

  // Return WFace on face-centered faceBox.
  m_util.PPMFaceValues(WFace,a_W,numprim,
                       m_useCharLimiting || m_usePrimLimiting,
                       a_dir,faceBox,m_currentTime,m_gdnvPhysics);

  // To save on storage, we use the input values as temporaries for the
  // delta's
  a_WMinus.setVal(0.0);
  a_WPlus .setVal(0.0);

  a_WMinus -= a_W;
  a_WPlus  -= a_W;

  WFace.shiftHalf(a_dir,1);
  a_WMinus += WFace;

  WFace.shift(a_dir,-1);
  a_WPlus  += WFace;

  FArrayBox lambda(a_box, numprim);
  m_gdnvPhysics->charValues(lambda, a_W, a_dir,a_box);

  if (m_useCharLimiting && m_usePrimLimiting)
    {
      MayDay::Error("PatchGodunov::PPMNormalPred: Attempt to limit slopes in primitive AND characteristic coordinates - not implemented");
    }

  // Apply limiter on characteristic or primitive variables. Either
  // way, must end up with characteristic variables to pass to to the
  // normal predictor utility.  Currently, cannot do both.

  // If doing characteristic limiting then transform before limiting
  if (m_useCharLimiting)
    {
      // Transform from primitive to characteristic variables
      m_gdnvPhysics->charAnalysis(a_WMinus,a_W,a_dir,a_box);
      m_gdnvPhysics->charAnalysis(a_WPlus ,a_W,a_dir,a_box);
    }

  if (m_useCharLimiting || m_usePrimLimiting)
    {
      // Do slope limiting
      // m_util.PPMLimiter(a_WMinus,a_WPlus,numprim,a_box);

      // petermc, 4 Sep 2008:  included a_W and a_dir in argument list
      m_util.PPMLimiter(a_WMinus, a_WPlus, a_W, numprim, a_dir, a_box);

      // Do slope flattening
      if (m_useFlattening)
        {
          m_util.applyFlattening(a_WMinus,a_flat,a_box);
          m_util.applyFlattening(a_WPlus ,a_flat,a_box);
        }
    }

  // If not doing characteristic limiting then transform after any limiting
  if (!m_useCharLimiting)
    {
      // Transform from primitive to characteristic variables
      m_gdnvPhysics->charAnalysis(a_WMinus,a_W,a_dir,a_box);
      m_gdnvPhysics->charAnalysis(a_WPlus ,a_W,a_dir,a_box);
    }

  // To the normal prediction in characteristic variables
  m_util.PPMNormalPred(a_WMinus,a_WPlus,lambda,a_dt / a_dx,numprim,a_box);

  // Construct the increments to the primitive variables
  m_gdnvPhysics->charSynthesis(a_WMinus,a_W,a_dir,a_box);
  m_gdnvPhysics->charSynthesis(a_WPlus ,a_W,a_dir,a_box);

  // Apply a physics-dependent post-normal predictor step:
  // For example:
  //   - adjust/bound delta's so constraints on extrapolated primitive
  //     quantities are enforced (density and pressure > 0).
  //   - compute source terms that depend on the spatially varying
  //     coefficients.
  m_gdnvPhysics->postNormalPred(a_WMinus,a_WPlus,a_W,a_dt,a_dx,a_dir,a_box);

  // Compute the state from the increments
  a_WMinus += a_W;
  a_WPlus  += a_W;
}

void PatchGodunov::highOrderLimiter(bool a_highOrderLimiter)
{
  CH_assert(m_isDefined);
  m_highOrderLimiter = a_highOrderLimiter;
  m_util.highOrderLimiter(a_highOrderLimiter);
}

GodunovPhysics* PatchGodunov::getGodunovPhysicsPtr()
{
  return m_gdnvPhysics;
}

// Return true if everything is defined and setup
bool PatchGodunov::isDefined() const
{
  return m_isDefined        &&
         m_isCurrentTimeSet &&
         m_isCurrentBoxSet  ;
}

Real PatchGodunov::dx() const
{
  CH_assert(m_isDefined);
  return m_dx;
}

//-----------------------------------------------------------------------
const ProblemDomain&
PatchGodunov::problemDomain() const
{
  CH_assert(m_isDefined);
  return m_domain;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
