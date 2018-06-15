#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectionPhysics.H"
#include "AdvectionPhysicsF_F.H"

#include "LoHiSide.H"

AdvectionPhysics::AdvectionPhysics():GodunovPhysics()
{
  // Bogus values
  m_numCons = -1;
  m_numPrim = -1;
  m_numFlux = -1;

  // Undefined until these are defined
  m_nCompDefined = false;
}

AdvectionPhysics::~AdvectionPhysics()
{
}

void AdvectionPhysics::setNComp(int nComp)
{
  CH_assert(nComp > 0);

  // All these are the same variable for advection
  m_numCons = nComp;
  m_numPrim = nComp;
  m_numFlux = nComp;

  m_nCompDefined = true;
}

// Compute the maximum wave speed
Real AdvectionPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
                                       const Box&       a_box)
{
  CH_assert(isDefined());

  // note that for this advection problem, the max wavespeed
  // is based entirely on the cell-centered velocity
  Real maxSpeed, speed = 0.0;
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox& vel = *m_cellVelPtr;
      speed = vel.norm(0,dir,1);
      if (maxSpeed < speed) maxSpeed = speed;
    }

  return maxSpeed;
}

// Factory method - this object is its own factory.  It returns a pointer
// to new AdvectionPhysics object with the same definition as this object.
GodunovPhysics* AdvectionPhysics::new_godunovPhysics() const
{
  // Make the new object
  AdvectionPhysics* newPtr = new AdvectionPhysics();

  // Define the same domain and dx
  newPtr->define(m_domain,m_dx);

  // Pass the current initial and boundary condition (IBC) object to this new
  // patch integrator so it can create its own IBC object and define it as it
  // needs to (i.e. new domain and grid spacing if necessary)
  newPtr->setPhysIBC(getPhysIBC());

  newPtr->m_numCons = m_numCons;
  newPtr->m_numPrim = m_numPrim;
  newPtr->m_numFlux = m_numFlux;

  newPtr->m_nCompDefined = m_nCompDefined;

  GodunovPhysics* retval = static_cast<GodunovPhysics*>(newPtr);

  // Return the new object
  return retval;
}

// Number of conserved variables
int AdvectionPhysics::numConserved()
{
  return m_numCons;
}

// Names of the conserved variables
Vector<string> AdvectionPhysics::stateNames()
{
  Vector<string> retval(m_numPrim, "scalar");

  return retval;
}

// Number of flux variables
int AdvectionPhysics::numFluxes()
{
  // In some computations there may be more fluxes than conserved variables
  return m_numFlux;
}

void AdvectionPhysics::getFlux(FArrayBox&       a_flux,
                               const FArrayBox& a_WHalf,
                               const int&       a_dir,
                               const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_flux.nComp() == a_WHalf.nComp());

  FArrayBox& advVelDir = (*m_advVelPtr)[a_dir];
  // flux is WHalf*advVel
  a_flux.copy(a_WHalf, a_box);
  // componentwise multiplication, since all components of flux
  // get multiplied by 0th componenent of advVel
  for (int comp=0; comp<a_flux.nComp(); comp++)
    {
      a_flux.mult(advVelDir, a_box, 0, comp, 1);
    }
}


// Is everthing defined
bool AdvectionPhysics::isDefined() const
{
  return (m_isDefined && m_nCompDefined);
}

// Number of primitive variables
int AdvectionPhysics::numPrimitives()
{
  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return m_numPrim;
}

void AdvectionPhysics::charAnalysis(FArrayBox&       a_dW,
                                    const FArrayBox& a_W,
                                    const int&       a_dir,
                                    const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_dW.nComp() == a_W.nComp());

  // This is the identity mapping for advection - do nothing
}

void AdvectionPhysics::charSynthesis(FArrayBox&       a_dW,
                                     const FArrayBox& a_W,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  CH_assert(a_dW.nComp() == a_W.nComp());

  // This is the identity mapping for advection - do nothing
}

void AdvectionPhysics::charValues(FArrayBox&       a_lambda,
                                  const FArrayBox& a_W,
                                  const int&       a_dir,
                                  const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_lambda.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  CH_assert(m_cellVelPtr != NULL);
  FArrayBox& cellVel = *m_cellVelPtr;

  int nComp = a_lambda.nComp();
  int copyComp = 1;

  for (int comp = 0; comp < nComp; comp++)
  {
    a_lambda.copy(cellVel,a_box,a_dir,a_box,comp,copyComp);
  }
}

// Increment the source term using the primitive state (not needed)
void AdvectionPhysics::incrementSource(FArrayBox&       a_S,
                                       const FArrayBox& a_W,
                                       const Box&       a_box)
{
  CH_assert (isDefined());

  CH_assert (a_S.box().contains(a_box));
  CH_assert (a_W.box().contains(a_box));

  CH_assert (a_S.nComp() == a_W.nComp());

  // The source term does not explicitly depend on the current primitive state
}

// Compute a Riemann problem and generate fluxes at the faces
void AdvectionPhysics::riemann(FArrayBox&       a_WGdnv,
                               const FArrayBox& a_WLeft,
                               const FArrayBox& a_WRight,
                               const FArrayBox& a_W,
                               const Real&      a_time,
                               const int&       a_dir,
                               const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_WGdnv.box().contains(a_box));

  CH_assert(m_advVelPtr != NULL);
  FluxBox& advVel = *m_advVelPtr;
  // CH_assert(advVel.box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();

  CH_assert(a_WGdnv .nComp() == numPrim);
  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Solution to the Riemann problem

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft .box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));

  FArrayBox& advVelDir = advVel[a_dir];
  CH_assert(advVelDir.box().contains(a_box));

  // Riemann solver computes WGdnv all edges that are not on the physical
  // boundary.
  FORT_RIEMANNF(CHF_FRA(a_WGdnv),
                CHF_CONST_FRA(shiftWLeft),
                CHF_CONST_FRA(shiftWRight),
                CHF_CONST_FRA1(advVelDir,0),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

void AdvectionPhysics::postNormalPred(FArrayBox&       a_WMinus,
                                      FArrayBox&       a_WPlus,
                                      const FArrayBox& a_W,
                                      const Real&      a_dt,
                                      const Real&      a_dx,
                                      const int&       a_dir,
                                      const Box&       a_box)
{
  CH_assert(isDefined());

  // Nothing needs to be done here
}

void AdvectionPhysics::quasilinearUpdate(FArrayBox&       a_AdWdx,
                                         const FArrayBox& a_WHalf,
                                         const FArrayBox& a_W,
                                         const Real&      a_scale,
                                         const int&       a_dir,
                                         const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_AdWdx.box().contains(a_box));
  CH_assert(a_W    .box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();

  CH_assert(a_AdWdx.nComp() == numPrim);
  CH_assert(a_WHalf.nComp() == numPrim);
  CH_assert(a_W    .nComp() == numPrim);

  const FArrayBox& cellVel = *m_cellVelPtr;

  FORT_QUASILINEARUPDATE(CHF_FRA(a_AdWdx),
                         CHF_CONST_FRA(a_WHalf),
                         CHF_CONST_FRA1(cellVel,a_dir),
                         CHF_CONST_REAL(a_scale),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_box));
}

// Compute the primitive variables from the conserved variables
void AdvectionPhysics::consToPrim(FArrayBox&       a_W,
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

// Interval within the primitive variables corresponding to the velocities
// since this doesn't apply to this set of equations, return a bogus interval
Interval AdvectionPhysics::velocityInterval()
{
  MayDay::Error("AdvectionPhysics::velocityInterval - not defined");

  Interval retval(-1,-1);
  return retval;
}

// Component index within the primitive variables of the pressure
// since this doesn't apply to this set of equations, return a bogus value
int AdvectionPhysics::pressureIndex()
{
  MayDay::Error("AdvectionPhysics::pressureIndex - not defined");

  return -1;
}

// Used to limit the absolute value of a "pressure" difference
// since this doesn't apply to this set of equations, return a bogus value
Real AdvectionPhysics::smallPressure()
{
  MayDay::Error("AdvectionPhysics::smallPressure - not defined");

  return -1.0;
}

// Component index within the primitive variables of the bulk modulus
// since this doesn't apply to this set of equations, return a bogus value
int AdvectionPhysics::bulkModulusIndex()
{
  MayDay::Error("AdvectionPhysics::bulkModulusIndex - not defined");

  return -1;
}

// Set face-centered advection velocity pointer
void
AdvectionPhysics::setAdvVelPtr(FluxBox* a_advVelPtr)
{
  m_advVelPtr = a_advVelPtr;
}

// Get face-centered advection velocity pointer
const FluxBox*
AdvectionPhysics::advectionVelPtr() const
{
  return m_advVelPtr;
}

// Set cell-centered advection velocity pointer
void
AdvectionPhysics::setCellVelPtr(FArrayBox* a_cellVelPtr)
{
  m_cellVelPtr = a_cellVelPtr;
}

// Get cell-centered advection velocity pointer
const FArrayBox*
AdvectionPhysics::cellVelPtr() const
{
  return m_cellVelPtr;
}
