#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <string>
using std::string;

#include "OldPatchGodunov.H"
#include "OldLoHiCenter.H"
#include "LoHiSide.H"

#include "OldPatchGodunovF_F.H"
#include "UsingNamespace.H"

// Flag everything as not defined or set
OldPatchGodunov::OldPatchGodunov()
{
  m_isDefined = false;
  m_isBCSet = false;
  m_isSlopeSet = false;
  m_isArtViscSet = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet = false;

  m_bc = NULL;
}

OldPatchGodunov::~OldPatchGodunov()
{
  // Delete the initial/boundary condition object - if it exists
  if (m_bc != NULL)
  {
    delete m_bc;
  }
}

// Define this object and the boundary condition object
void OldPatchGodunov::define(const ProblemDomain& a_domain,
                             const Real&    a_dx)
{
  CH_assert(m_isBCSet);

  // Store the domain and grid spacing
  m_domain = a_domain;
  m_dx = a_dx;

  // Set the domain and grid spacing in the boundary condition object
  m_bc->define(m_domain,m_dx);

  m_isDefined = true;
}

// Set the boundary condition object
void OldPatchGodunov::setPhysIBC(OldPhysIBC* a_bc)
{
  // Delete old boundary condition object - if any
  if (m_bc != NULL)
  {
    delete m_bc;
  }

  // Store new boundary condition object
  m_bc = a_bc->new_physIBC();

  // just in case we're re-defining the BC
  if (m_isDefined)
  {
    m_bc->define(m_domain, m_dx);
  }

  m_isBCSet = true;
}

// Return the current boundary condition object
OldPhysIBC* OldPatchGodunov::getPhysIBC() const
{
  CH_assert(m_isBCSet);

  return m_bc;
}

// Set the parameters for slope computation
void OldPatchGodunov::setSlopeParameters(bool a_fourthOrderSlopes,
                                         bool a_flattening,
                                         bool a_limitSlopes)
{
  // Slope flattening is only allowed with 4th order slopes
  CH_assert(a_fourthOrderSlopes || !a_flattening);

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_fourthOrderSlopes;
  m_useFlattening = a_flattening;
  m_limitSlopes = a_limitSlopes;

  m_isSlopeSet = true;
}

// Return true if 4th order slopes are being computed, false if 2nd order
bool OldPatchGodunov::useFourthOrderSlopes()
{
  CH_assert(m_isSlopeSet);

  return m_useFourthOrderSlopes;
}

// Return true if slope flattening is being applied
bool OldPatchGodunov::useFlattening()
{
  CH_assert(m_isSlopeSet);

  return m_useFlattening;
}

// Return true if slopes should be limited
bool OldPatchGodunov::limitSlopes()
{
  // CH_assert(m_isSlopeSet);

  return m_limitSlopes;
}

// Flag the use of artificial viscosity and set the coefficient
void OldPatchGodunov::setArtificialViscosity(bool a_useArtificialViscosity,
                                             Real a_artificialViscosity)
{
  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity > 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity = a_artificialViscosity;

  m_isArtViscSet = true;
}

// Return true if artificial viscosity is being applied
bool OldPatchGodunov::useArtificialViscosity()
{
  CH_assert(m_isArtViscSet);

  return m_useArtificialViscosity;
}

// Return the artificial viscosity coefficient
Real OldPatchGodunov::artificialViscosityCoefficient()
{
  CH_assert(m_isArtViscSet);
  CH_assert(m_useArtificialViscosity);

  return m_artificialViscosity;
}

// Set the current physical time - used for time dependent boundary conditions
void OldPatchGodunov::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime = a_currentTime;
  m_isCurrentTimeSet = true;
}

// Set the current box for the updateState() call
void OldPatchGodunov::setCurrentBox(const Box& a_currentBox)
{
  m_currentBox = a_currentBox;
  m_isCurrentBoxSet = true;
}

// Update the conserved variables, return the final fluxes used for this,
// and return the maximum wave speed on this patch/grid
void OldPatchGodunov::updateState(FArrayBox&       a_U,
                                  FArrayBox        a_F[CH_SPACEDIM],
                                  Real&            a_maxWaveSpeed,
                                  const FArrayBox& a_S,
                                  const Real&      a_dt,
                                  const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_box == m_currentBox);

  // Get the number of various variables
  int numFlux = numFluxes();

  // Create a temp storage for a_U - RS 06/10/02
  FArrayBox Uold(a_U.box(),a_U.nComp());
  Uold.copy(a_U);

  // The current box
  const Box& curBox = m_currentBox;

  // allocate flux storage
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      // This box is face centered in direction "dir1", is one bigger
      // than the input box in all directions except "dir1"
      Box fluxBox;

      fluxBox = curBox;
      fluxBox.grow(1);
      fluxBox.grow(dir1,-1);
      fluxBox &= m_domain;
      fluxBox.surroundingNodes(dir1);

      FArrayBox& F = a_F[dir1];
      F.resize(fluxBox,numFlux);
    }

  computeFluxes(a_U,a_F, a_S, a_dt, a_box);

  // Update conserved variables to the next time step using the final
  // flux in each direction
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      finalUpdate(a_U,a_F[dir1],a_dt/m_dx,dir1,curBox);
    }

  // Post processing on the conserved variables
  postUpdateCons(a_U,Uold,a_dt,m_dx,curBox);

  // Get and return the maximum wave speed on this patch/grid
  a_maxWaveSpeed = getMaxWaveSpeed(a_U,m_currentBox);
}

void OldPatchGodunov::computeFluxes(FArrayBox&       a_U,
                                    FArrayBox        a_F[CH_SPACEDIM],
                                    const FArrayBox& a_S,
                                    const Real&      a_dt,
                                    const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_box == m_currentBox);

  // Get the number of various variables
  int numFlux = numFluxes();
  int numPrim = numPrimitives();

  // The current box
  Box curBox = m_currentBox;

  // Boxes for face centered state - used for the riemann() and
  // artificialViscosity() calls
  Box faceBox[SpaceDim];

  // Boxes for face centered fluxes
  Box fluxBox[SpaceDim];

  // Boxes for cell centered state - used for the updatePrim() calls
  Box ccBox[SpaceDim];

  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
  {
    // This box is face centered in direction "dir1", is one bigger than the
    // input box in all directions except "dir1" and stays one cell away from
    // the domain boundary in "dir1"
    faceBox[dir1] = curBox;
    faceBox[dir1].grow(1);
    faceBox[dir1] &= m_domain;
    faceBox[dir1].grow(dir1,-1);
    faceBox[dir1].surroundingNodes(dir1);

    // This box is face centered in direction "dir1", is one bigger than the
    // input box in all directions except "dir1"
    fluxBox[dir1] = curBox;
    fluxBox[dir1].grow(1);
    fluxBox[dir1].grow(dir1,-1);
    fluxBox[dir1] &= m_domain;
    fluxBox[dir1].surroundingNodes(dir1);

    CH_assert(a_F[dir1].box().contains(fluxBox[dir1]));

    // This is the cell centered analog of "fluxBox"
    ccBox[dir1] = curBox;
    ccBox[dir1].grow(1);
    ccBox[dir1].grow(dir1,-1);
    ccBox[dir1] &= m_domain;
  }

  // Primitive variables
  Box WBox = a_U.box();
  WBox &= m_domain;

  FArrayBox W(WBox,numPrim);

  // Calculate the primitive variables from the conserved variables
  consToPrim(W,a_U,WBox);

  // Define the box where slopes will be needed (one larger than the final
  // update box)
  Box slopeBox = curBox;
  slopeBox.grow(1);
  slopeBox &= m_domain;

  // Compute flattening once for all slopes if needed
  FArrayBox flattening;
  if (useFlattening())
  {
    flattening.define(slopeBox,1);
    computeFlattening(flattening,W,slopeBox);
  }

  // Intermediate, extrapolated primitive variables
  FArrayBox WMinus[SpaceDim];
  FArrayBox WPlus [SpaceDim];

  // Initial fluxes
  FArrayBox F1[SpaceDim];

  // Compute initial fluxes
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
  {
    // Size the intermediate, extrapolated primitive variables
    WMinus[dir1].resize(slopeBox,numPrim);
    WPlus [dir1].resize(slopeBox,numPrim);

    // Compute slopes
    int numSlope = numSlopes();
    FArrayBox dW(slopeBox,numSlope);

    slope(dW,W,flattening,dir1,slopeBox);

    // Compute predictor step to obtain extrapolated primitive variables
    normalPred(WMinus[dir1],WPlus[dir1],W,dW,a_dt/m_dx,dir1,slopeBox);

    // If the source term is valid add it to the primitive quantities
    if (!a_S.box().isEmpty())
    {
      incrementWithSource(WMinus[dir1],a_S,0.5*a_dt,slopeBox);
      incrementWithSource(WPlus [dir1],a_S,0.5*a_dt,slopeBox);
    }

    // Solve the Riemann problem and get fluxes
    F1[dir1].resize(fluxBox[dir1],numFlux);
    F1[dir1].setVal(0.0);
    riemann(F1[dir1],WPlus[dir1],WMinus[dir1],dir1,faceBox[dir1]);

    // Use the user supplied PhysBC object to obtain boundary fluxes
    m_bc->fluxBC(F1[dir1],W,WMinus[dir1],dir1,Side::Lo,m_currentTime);
    m_bc->fluxBC(F1[dir1],W,WPlus [dir1],dir1,Side::Hi,m_currentTime);
  }

#if (CH_SPACEDIM == 3)
  // In 3D, compute some additional intermediate fluxes
  //
  // NOTE:  The diagonal entries of this array of fluxes are not
  // used and will not be defined.
  FArrayBox F2[SpaceDim][SpaceDim];

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
        FArrayBox WMinusTemp(WMinus[dir1].box(),numPrim);
        FArrayBox WPlusTemp (WPlus [dir1].box(),numPrim);

        // Copy data for in place modification
        WMinusTemp.copy(WMinus[dir1]);
        WPlusTemp .copy(WPlus [dir1]);

        // Update the current, extrapolated primitive variable using a flux
        // in a different direction
        updatePrim(WMinusTemp,WPlusTemp,F1[dir2],(1.0/3.0)*a_dt/m_dx,dir2,ccBox[dir2]);

        // Solve the Riemann problem and get fluxes
        F2[dir1][dir2].resize(fluxBox[dir1],numFlux);
        F2[dir1][dir2].setVal(0.0);

        riemann(F2[dir1][dir2],WPlusTemp,WMinusTemp,dir1,faceBox[dir1]);

        // Use the user supplied PhysBC object to obtain boundary fluxes
        m_bc->fluxBC(F2[dir1][dir2],W,WMinusTemp,dir1,Side::Lo,m_currentTime);
        m_bc->fluxBC(F2[dir1][dir2],W,WPlusTemp ,dir1,Side::Hi,m_currentTime);
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
        updatePrim(WMinus[dir1],WPlus[dir1],F1[dir2],(1.0/2.0)*a_dt/m_dx,dir2,ccBox[dir2]);
#elif (CH_SPACEDIM == 3)
        // In 3D, find a direction different from the two above
        int dir3 = 3 - dir1 - dir2;

        // Update the conservative state using both corrected fluxes in
        // the other two directions
        updatePrim(WMinus[dir1],WPlus[dir1],F2[dir2][dir3],(1.0/2.0)*a_dt/m_dx,dir2,ccBox[dir2]);
#else
        // Only 2D and 3D should be possible
        MayDay::Error("OldPatchGodunov::computeFluxes - CH_SPACEDIM not 2 or 3!");
#endif
      }
    }

    // Solve the Riemann problem and get fluxes, these final fluxes will
    // be returned
    FArrayBox& F = a_F[dir1];
    // do setval here to prevent uninitialized memory read on some machines
    F.setVal(666.666);
    riemann(F,WPlus[dir1],WMinus[dir1],dir1,faceBox[dir1]);

    // Use the user supplied PhysBC object to obtain boundary fluxes
    m_bc->fluxBC(F,W,WMinus[dir1],dir1,Side::Lo,m_currentTime);
    m_bc->fluxBC(F,W,WPlus [dir1],dir1,Side::Hi,m_currentTime);

    // Apply artificial viscosity if reqested
    if (useArtificialViscosity())
    {
      // Compute the divergence of the velocity for artificial viscosity
      FArrayBox divVelDir;

      divVelDir.define(fluxBox[dir1],1);
      divVelDir.setVal(0.0);

      divVel(divVelDir,W,dir1,slopeBox);

      // Apply artificial viscosity to the fluxes
      artificialViscosity(F,a_U,divVelDir,dir1,faceBox[dir1]);

      // Use the user supplied PhysBC object to apply artificial viscosity
      // to the boundary fluxes
      m_bc->artViscBC(F,a_U,divVelDir,dir1,m_currentTime);
    }
  }
}

// Generate default names for the conserved variables, "variable#"
Vector<string> OldPatchGodunov::stateNames()
{
  Vector<string> retval;

  int cnum = numConserved();

  for (int ivar = 0; ivar < cnum; ivar++)
  {
    char varNameChar[80];
    sprintf(varNameChar,"variable%d",ivar);
    retval.push_back(string(varNameChar));
  }

  return retval;
}

// Return true if everything is defined and setup
bool OldPatchGodunov::isDefined() const
{
  return m_isDefined        &&
         m_isBCSet          &&
         m_isSlopeSet       &&
         m_isArtViscSet     &&
         m_isCurrentTimeSet &&
         m_isCurrentBoxSet;
}

// Compute the flattening coefficients from the primitive variables
void OldPatchGodunov::computeFlattening(FArrayBox&       a_flattening,
                                        const FArrayBox& a_W,
                                        const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_W.box().contains(a_box));

  // The current direction
  int idir;

  // The directional flattening coefficients
  FArrayBox zetaDir(a_box,SpaceDim);

  // The divergence of the velocity
  FArrayBox dVel(a_box,SpaceDim);

  // The interval of the primitive variables corresponding to the velocity
  Interval velInterval= velocityInterval();
  int v0index = velInterval.begin();

  // Get the directional flattening coefficients in each direction
  for (idir = 0; idir < SpaceDim; idir++)
  {
    // A box one larger (in direction "idir") than the final result box
    Box box1 = a_box;
    box1.grow(idir,1);

    // A box two larger (in direction "idir") than the final result box
    Box box2 = a_box;
    box2.grow(idir,2);

    // A box three larger (in direction "idir") than the final result box
    Box box3 = a_box;
    box3.grow(idir,3);

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on "box3"
    Box loBox,hiBox,centerBox,entireBox;
    int hasLo,hasHi;

    oldLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                  box3,m_domain,idir);

    // Compute the first differences in "pressure"
    FArrayBox delta1p(entireBox,1);

    int pressIndex = pressureIndex();

    FORT_OLDGETGRADF(CHF_FRA1(delta1p,0),
                     CHF_CONST_FRA1(a_W,pressIndex),
                     CHF_CONST_INT(idir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on "box2"
    oldLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                  box2,m_domain,idir);

    // Compute the second differences in "pressure"
    FArrayBox delta2p(entireBox,1);

    FORT_OLDGETDPTWOF(CHF_FRA1(delta2p,0),
                      CHF_CONST_FRA1(delta1p,0),
                      CHF_CONST_INT(idir),
                      CHF_BOX(loBox),
                      CHF_CONST_INT(hasLo),
                      CHF_BOX(hiBox),
                      CHF_CONST_INT(hasHi),
                      CHF_BOX(centerBox));

    // Compute a 3-way minimum of the "bulk modulus"
    FArrayBox bulkMin(entireBox,1);

    int bulkIndex = bulkModulusIndex();

    FORT_OLDMIN3PTSF(CHF_FRA1(bulkMin,0),
                     CHF_CONST_FRA1(a_W,bulkIndex),
                     CHF_CONST_INT(idir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));

    // Use the first and second differences normalized by the 3-way minimum
    // computed above to generate directional flattening coefficients
    FArrayBox zetaTwiddleDir(entireBox,1);

    FORT_OLDGETFLATF(CHF_FRA1(zetaTwiddleDir,0),
                     CHF_CONST_FRA1(delta1p,0),
                     CHF_CONST_FRA1(delta2p,0),
                     CHF_CONST_FRA1(bulkMin,0),
                     CHF_BOX(entireBox));

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on "box1"
    oldLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                  box1,m_domain,idir);

    // Take a 3-way minimum of the directional flattening coefficients
    FORT_OLDMIN3PTSF(CHF_FRA1(zetaDir,idir),
                     CHF_CONST_FRA1(zetaTwiddleDir,0),
                     CHF_CONST_INT(idir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));

    // Compute each component of the divergence of the velocity
    FORT_OLDGETGRADF(CHF_FRA1(dVel,idir),
                     CHF_CONST_FRA1(a_W,v0index+idir),
                     CHF_CONST_INT(idir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));
  }

  // At each point, set the flattening coefficient to the minimum of all
  // the directional flattening coefficients if the divergence of the velocity
  // is negative, otherwise set it to 1 (no flattening).
  FORT_OLDMINFLATF(CHF_FRA1(a_flattening,0),
                   CHF_CONST_FRA(zetaDir),
                   CHF_CONST_FRA(dVel),
                   CHF_BOX(a_box));
}

// Compute the directional slopes (2nd or 4th order) of the primitive
// variables applying slope limiting and, possibly, slope flattening.
// Slopes are only computed for primitive variables with indices 0 through
// numSlopes()-1.
void OldPatchGodunov::slope(FArrayBox&       a_dW,
                            const FArrayBox& a_W,
                            const FArrayBox& a_flattening,
                            const int&       a_dir,
                            const Box&       a_box)
{
  // Number of slopes to compute
  int numSlope = numSlopes();

  CH_assert(a_dW.nComp() == numSlope);
  CH_assert(a_W.nComp() >= numSlope);

  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // A box two larger (in direction "a_dir") than the final result box
  Box box2 = a_box;
  box2.grow(a_dir,2);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.  The data used for the differences is
  // defined on "box2" for 4th order slope and "box1" for 2nd order slopes
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  if (useFourthOrderSlopes())
  {
    oldLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                  box2,m_domain,a_dir);
  }
  else
  {
    oldLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                  box1,m_domain,a_dir);
  }

  // Compute 2nd order slopes - including one sided differences
  FArrayBox delta2W(entireBox,numSlope);
  FArrayBox deltaWL(entireBox,numSlope);
  FArrayBox deltaWR(entireBox,numSlope);

  FORT_OLDSECONDSLOPEDIFFSF(CHF_FRA(delta2W),
                            CHF_FRA(deltaWL),
                            CHF_FRA(deltaWR),
                            CHF_CONST_FRA(a_W),
                            CHF_CONST_INT(numSlope),
                            CHF_CONST_INT(a_dir),
                            CHF_BOX(loBox),
                            CHF_CONST_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_CONST_INT(hasHi),
                            CHF_BOX(centerBox));

  // Apply the slope limiter, if desired
  if (limitSlopes())
    {
      applyLimiter(delta2W,deltaWL,deltaWR,a_dir,centerBox);
    }

  // Use the user supplied PhysBC object to obtain slopes at the boundary
  m_bc->setBdrySlopes(delta2W,a_W,a_dir,m_currentTime);

  // 4th order slopes requested
  if (useFourthOrderSlopes())
  {
    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on "box1"
    oldLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                  box1,m_domain,a_dir);

    // Compute 4th order slopes
    FArrayBox delta4W(entireBox,numSlope);
    CH_assert(delta4W.box().contains(entireBox));

    FORT_OLDFOURTHSLOPEDIFFSF(CHF_FRA(delta4W),
                              CHF_CONST_FRA(a_W),
                              CHF_CONST_FRA(delta2W),
                              CHF_CONST_INT(numSlope),
                              CHF_CONST_INT(a_dir),
                              CHF_BOX(loBox),
                              CHF_CONST_INT(hasLo),
                              CHF_BOX(hiBox),
                              CHF_CONST_INT(hasHi),
                              CHF_BOX(centerBox));

    // Apply the slope limiter
    if (limitSlopes())
      {
        applyLimiter(delta4W,deltaWL,deltaWR,a_dir,centerBox);
      }

    // Apply slope flattening if requested
    if (useFlattening())
    {
      FORT_OLDAPPLYFLATF(CHF_FRA(delta4W),
                         CHF_CONST_FRA1(a_flattening,0),
                         CHF_CONST_INT(numSlope),
                         CHF_BOX(entireBox));
    }

    // Copy the 4th order slopes for return
    a_dW.copy(delta4W);
  }
  else
  {
    // Copy the 2nd order slopes for return
    a_dW.copy(delta2W);
  }
}

// Increment the primitive variables with sources corresponding to the
// conserved variables - the default implementation does nothing (it is
// implemented to satisfy compilation for users without source terms).
void OldPatchGodunov::incrementWithSource(FArrayBox&       a_W,
                                          const FArrayBox& a_S,
                                          const Real&      a_scale,
                                          const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_S.box().contains(a_box));
  MayDay::Error("You have provided a source term without providing an instantiation of the incrementWithSource function in your OldPatchGodunov-derived class. Please provide one");
}

// Apply artificial viscosity to the fluxes
void OldPatchGodunov::artificialViscosity(FArrayBox&       a_F,
                                          const FArrayBox& a_U,
                                          const FArrayBox& a_divVel,
                                          const int&       a_dir,
                                          const Box&       a_box)
{
  // Get the artificial viscosity coefficient
  Real coeff = artificialViscosityCoefficient();

  // Apply the artificial viscosity
  FORT_OLDARTVISCF(CHF_FRA(a_F),
                   CHF_CONST_FRA(a_U),
                   CHF_CONST_FRA1(a_divVel,0),
                   CHF_CONST_REAL(coeff),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(a_box));
}

// Apply a van Leer limiter directly to the slopes.  No characteristic
// analysis/projection is done as this is specific to each application.
// The user should implement their own applyLimiter() in their derived class
// if characteristic analysis/projection is needed.
void OldPatchGodunov::applyLimiter(FArrayBox&       a_dW,
                                   const FArrayBox& a_dWLeft,
                                   const FArrayBox& a_dWRight,
                                   const int&       a_dir,
                                   const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_OLDVANLEERLIMITERF(CHF_FRA(a_dW),
                          CHF_CONST_FRA(a_dWLeft),
                          CHF_CONST_FRA(a_dWRight),
                          CHF_CONST_INT(a_dir),
                          CHF_BOX(a_box));
}

// Compute a face centered divergence of the velocity
void OldPatchGodunov::divVel(FArrayBox&       a_divVel,
                             const FArrayBox& a_W,
                             const int        a_dir,
                             const Box&       a_box)
{
  Box dveltanBox = a_divVel.box();
  dveltanBox.enclosedCells(a_dir);
  dveltanBox.grow(1);
  dveltanBox &= m_domain;

  // First, we need to calculate the directional derivatives of
  // the tangential components of velocity at the cell-centers.
  FArrayBox dveltan(dveltanBox,SpaceDim-1);

  // Get the interval of the primitive variables corresponding to the velocity
  Interval velInterval= velocityInterval();
  int v0index = velInterval.begin();

  // Go through the tangential directions
  for (int i = 0, dir = 0; dir < SpaceDim; ++dir)
  {
    if (dir != a_dir)
    {
      // This velocity component is tangential.  Build the box in which
      // d(v[dir])/d(x[dir]) is to be computed.
      Box primBox = a_divVel.box();
      primBox.enclosedCells(a_dir);
      primBox.grow(dir,1).grow(a_dir,1);

      // Compute where centered differences can be used and where one sided
      // differences need to be used.  The data used for the differences is
      // defined on "dveltanBox"
      Box gradBox,hiBox,loBox,centerBox;
      int hasLo,hasHi;

      oldLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,gradBox,
                    primBox,m_domain,dir);

      // Compute d(v[dir])/d(x[dir]).
      FORT_OLDGETGRADF(CHF_FRA1(dveltan,i),
                       CHF_CONST_FRA1(a_W,v0index+dir),
                       CHF_CONST_INT(dir),
                       CHF_BOX(loBox),
                       CHF_CONST_INT(hasLo),
                       CHF_BOX(hiBox),
                       CHF_CONST_INT(hasHi),
                       CHF_BOX(centerBox));
      ++i;
    }
  }

  // Now, we can calculate the divergence of the normal velocity
  // at the center normal-direction edges. To do this, we determine
  // which edges at which we have sufficient data to compute centered
  // estimates of h*(div(u)). At the remaining edges. i.e. those
  // corresponding to the physical boundaries, we use zeroth-order
  // extrapolation.

  Box divBox = a_divVel.box();
  divBox.enclosedCells(a_dir);
  divBox.grow(a_dir,1);

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  oldLoHiCenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                    divBox,m_domain,a_dir);

  // All of the boxes computed above are shifted so as to be cell-centered,
  // with the index of the cell center being identified with the low edge.
  // We then shift a_divVel to be compatible with that convention on input,
  // and undo the shift on output.  Basically, everything is made compatible
  // a_W (which is cell-centered).

  loBox.shiftHalf(a_dir,1);
  centerBox.shiftHalf(a_dir,1);
  hiBox.shiftHalf(a_dir,1);

  a_divVel.shiftHalf(a_dir,1);

  FORT_OLDDIVUEDGEF(CHF_FRA1(a_divVel,0),
                    CHF_CONST_FRA1(a_W,v0index+a_dir),
                    CHF_CONST_FRA(dveltan),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(centerBox));

  a_divVel.shiftHalf(a_dir,-1);
}

// Empty function - up to the application to provide this if
// necessary
void OldPatchGodunov::postUpdateCons(FArrayBox&       a_U,
                                     const FArrayBox& a_Uold,
                                     const Real&      a_dt,
                                     const Real&      a_dx,
                                     const Box&       a_box)
{

}
