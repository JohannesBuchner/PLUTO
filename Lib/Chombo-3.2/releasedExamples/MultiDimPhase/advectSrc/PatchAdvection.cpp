#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PatchAdvection.H"
#include "PatchAdvectionF_F.H"

#include "OldLoHiCenter.H"
#include "LoHiSide.H"
#include "UsingNamespace.H"

PatchAdvection::PatchAdvection():OldPatchGodunov()
{
}

// Factory method - this object is its own factory.  It returns a pointer
// to new OldPatchGodunov object with its initial and boundary condtions, slope
// parameters, and artificial viscosity information defined.
OldPatchGodunov* PatchAdvection::new_patchGodunov() const
{
  // Make the new object
  PatchAdvection* newPtr = new PatchAdvection();

  // Pass the current initial and boundary condition (IBC) object to this new
  // patch integrator so it can create its own IBC object and define it as it
  // needs to (i.e. new domain and grid spacing if necessary)
  newPtr->setPhysIBC(getPhysIBC());

  // Set the slope and artificial viscosity parameters
  newPtr->setSlopeParameters(m_useFourthOrderSlopes,
                             m_useFlattening,
                             m_limitSlopes);
  newPtr->setArtificialViscosity(m_useArtificialViscosity,
                                 m_artificialViscosity);

  // copy state info
  newPtr->m_numCons = m_numCons;
  newPtr->m_numPrim = m_numPrim;
  newPtr->m_numFlux = m_numFlux;
  newPtr->m_numSlope = m_numSlope;

  OldPatchGodunov* retval =
    static_cast<OldPatchGodunov*>(newPtr);

  // Return the new object
  return retval;
}



// Compute the maximum wave speed
Real PatchAdvection::getMaxWaveSpeed(const FArrayBox& a_U,
                                     const Box&       a_box)
{
  CH_assert(isDefined());

  // note that for this advection problem, the max wavespeed
  // is based entirely on the cell-centered velocity
  Real maxSpeed, speed = 0.0;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& vel = *m_cellVelPtr;
      speed = vel.norm(0,dir,1);
      if (maxSpeed < speed) maxSpeed = speed;
    }

  return maxSpeed;
}

// Number of conserved variables
int PatchAdvection::numConserved()
{
  return m_numCons;
}

// Names of the conserved variables
Vector<string> PatchAdvection::stateNames()
{
  // just hack around this for now
  Vector<string> retval(m_numPrim, "scalar");

  return retval;
}

// Number of flux variables
int PatchAdvection::numFluxes()
{
  // In some computations there may be more fluxes than conserved variables
  return m_numFlux;
}

// Number of primitive variables
int PatchAdvection::numPrimitives()
{
  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return m_numPrim;
}

//  Number of primitive variables for which slopes are computed
int PatchAdvection::numSlopes()
{
  // This may be less than the number of primitive variables for the
  // reason given in numPrimitives() comments
  return m_numSlope;
}

// Compute the primitive variables from the conserved variables
void PatchAdvection::consToPrim(FArrayBox&       a_W,
                                const FArrayBox& a_U,
                                const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  // since the primitive and conserved variables are the same in
  // this set of equations, even I can do this one...
  a_W.copy(a_U, a_box);


}

// Compute the conserved variables from the primitive variables
void PatchAdvection::primToCons(FArrayBox&       a_U,
                                const FArrayBox& a_W,
                                const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  // since primitive variables are conserved in this problem,
  // this is easy
  a_U.copy(a_W, a_box);


}

// Extrapolate the primitive variables to the cell faces
void PatchAdvection::normalPred(FArrayBox&       a_WMinus,
                                FArrayBox&       a_WPlus,
                                const FArrayBox& a_W,
                                const FArrayBox& a_dW,
                                const Real&      a_scale,
                                const int&       a_dir,
                                const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert (m_cellVelPtr != NULL);
  FArrayBox& cellVel = *m_cellVelPtr;
  CH_assert(cellVel.box().contains(a_box));

  FORT_NORMALPREDF(CHF_FRA(a_WMinus),
                   CHF_FRA(a_WPlus),
                   CHF_CONST_FRA(a_W),
                   CHF_CONST_FRA(a_dW),
                   CHF_CONST_FRA(cellVel),
                   CHF_CONST_REAL(a_scale),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(a_box));
}

// increment with source term
void PatchAdvection::incrementWithSource(FArrayBox& a_W,
                                         const FArrayBox& a_S,
                                         const Real& a_scale,
                                         const Box& a_box)
{
  CH_assert (isDefined());

  CH_assert (a_W.box().contains(a_box));
  CH_assert (a_S.box().contains(a_box));
  // may eventually want to have version which only acts on a subset
  // of components, but for now, this is required (can always set non-active
  // components to 0):
  CH_assert (a_S.nComp() == a_W.nComp());

  // now we just add in the source term, right? (do i need to do
  // cons->prim, increment w/ source, prim->cons? -- in this case, since
  // cons == prim, it's not important, but as an example, it would be
  // useful to actually have it correct in the code.
  //

  // for the moment, skip the cons->prim->cons conversion
  // would be nice if FArrayBox had a scaleMult function which
  // worked on a sub-box (scale-mult functions only work on
  // entire domain) -- this is a problem because of the const-ness of
  // a_S
  // this is downright silly
  FArrayBox tempSrc(a_S.box(), a_S.nComp());
  tempSrc.copy(a_S);
  tempSrc *= a_scale;
  int startComp = 0;
  int numComp = a_S.nComp();
  a_W.plus(tempSrc, a_box, startComp, startComp, numComp);


}


// Compute a Riemann problem and generate fluxes at the faces
void PatchAdvection::riemann(FArrayBox&       a_F,
                             const FArrayBox& a_WLeft,
                             const FArrayBox& a_WRight,
                             const int&       a_dir,
                             const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_F.box().contains(a_box));

  CH_assert(m_advVelPtr != NULL);
  FluxBox& advVel = *m_advVelPtr;
  //CH_assert(advVel.box().contains(a_box));


  // Get the numbers of relevant variables
  int numPrim = numPrimitives();
  int numFlux = numFluxes();

  CH_assert(a_F.nComp() == numFlux);
  CH_assert(a_WLeft.nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Solution to the Riemann problem
  FArrayBox Wgdnv(a_box,numPrim);

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft.box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));

  FArrayBox& advVelDir = advVel[a_dir];
  CH_assert(advVelDir.box().contains(a_box));

  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.
  FORT_RIEMANNF(CHF_FRA(Wgdnv),
                CHF_CONST_FRA(shiftWLeft),
                CHF_CONST_FRA(shiftWRight),
                CHF_CONST_FRA1(advVelDir,0),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));

  // Shift the left and right primitive variable boxes back to their original
  // position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);

  a_F.copy(Wgdnv);
}

// Update the primitive variables using face-centered variables
void PatchAdvection::updatePrim(FArrayBox&       a_WMinus,
                                FArrayBox&       a_WPlus,
                                const FArrayBox& a_F,
                                const Real&      a_scale,
                                const int&       a_dir,
                                const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_WPlus.box().contains(a_box));
  CH_assert(a_WMinus.box().contains(a_box));

  // Number of conserved variables
  int numCons = numConserved();

  // Temporary storage for converting primitive variables to conserved
  // variables and back again (after update with the fluxes)
  FArrayBox U(a_box,numCons);

  // for W+-, convert primitive quantities to conservative quantities,
  // update with the fluxes, and convert back
  primToCons(U,a_WMinus,a_box);
  updateCons(U,a_F,a_scale,a_dir,a_box);
  consToPrim(a_WMinus,U,a_box);

  primToCons(U,a_WPlus,a_box);
  updateCons(U,a_F,a_scale,a_dir,a_box);
  consToPrim(a_WPlus,U,a_box);
}

// Update the conserved variable using fluxes and a scaling factor
void PatchAdvection::updateCons(FArrayBox&       a_U,
                                const FArrayBox& a_F,
                                const Real&      a_scale,
                                const int&       a_dir,
                                const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert (cellVelPtr() != NULL);

  const FArrayBox& cellVel = *cellVelPtr();
  CH_assert (cellVel.box().contains(a_box));

  FORT_UPDATECONSADF(CHF_FRA(a_U),
                     CHF_CONST_FRA(a_F),
                     CHF_CONST_FRA(cellVel),
                     CHF_CONST_REAL(a_scale),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(a_box));
}

// Update the conserved variable using fluxes and a scaling factor
void PatchAdvection::finalUpdate(FArrayBox&       a_U,
                                 const FArrayBox& a_F,
                                 const Real&      a_scale,
                                 const int&       a_dir,
                                 const Box&       a_box)
{
  // in this case, we just call updateCons..
  updateCons(a_U, a_F, a_scale, a_dir, a_box);
}



// Interval within the primitive variables corresponding to the velocities
// since this doesn't apply to this set of equations, return a bogus interval
Interval PatchAdvection::velocityInterval()
{
  Interval retval(-1,-1);

  return retval;
}

// Component index within the primitive variables of the pressure
// since this doesn't apply to this set of equations, return a bogus value
int PatchAdvection::pressureIndex()
{
  return -1;
}

// Component index within the primitive variables of the bulk modulus
// since this doesn't apply to this set of equations, return a bogus value
int PatchAdvection::bulkModulusIndex()
{
  return -1;
}

// set advection velocity pointer
void
PatchAdvection::setAdvVelPtr(FluxBox* a_advVelPtr)
{
  m_advVelPtr = a_advVelPtr;
}


// access advection velocity
const FluxBox*
PatchAdvection::advectionVelPtr() const
{
  return m_advVelPtr;
}

// set cell-centered velocity pointer
void
PatchAdvection::setCellVelPtr(FArrayBox* a_cellVelPtr)
{
  m_cellVelPtr = a_cellVelPtr;
}


// access advection velocity
const FArrayBox*
PatchAdvection::cellVelPtr() const
{
  return m_cellVelPtr;
}


// set number of variables
void
PatchAdvection::setNumVar(int a_numVar)
{
  // note that numCons == numPrim == numFlux == numVar
  // in this case.
  m_numCons = a_numVar;
  m_numPrim = a_numVar;
  m_numCons = a_numVar;
  m_numFlux = a_numVar;
  m_numSlope = a_numVar;
}

