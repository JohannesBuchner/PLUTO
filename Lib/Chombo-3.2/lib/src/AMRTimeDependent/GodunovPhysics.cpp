#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GodunovPhysics.H"
#include "GodunovPhysicsF_F.H"
#include "NamespaceHeader.H"

GodunovPhysics::GodunovPhysics()
{
  m_isDefined = false;
  m_isBCSet = false;
  m_bc = NULL;
  m_useFourthOrderArtificialViscosity = false;
}

PhysIBC* GodunovPhysics::getPhysIBC() const
{
  CH_assert(m_isBCSet);
  return m_bc;
}

void GodunovPhysics::setPhysIBC(PhysIBC* a_bc)
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

GodunovPhysics::~GodunovPhysics()
{
  if (m_bc != NULL)
  {
    delete m_bc;
  }
}

void GodunovPhysics::define(const ProblemDomain& a_domain,
                            const Real&          a_dx)
{
  m_domain    = a_domain;
  m_dx        = a_dx;
  m_isDefined = true;

  if (m_bc != NULL)
  {
    m_bc->define(m_domain, m_dx);
  }

  m_util.define(m_domain, m_dx);
}

void GodunovPhysics::setCurrentBox(const Box& a_currentBox)
{
  // Do nothing but assert the object has been defined
  CH_assert(isDefined());
}
void GodunovPhysics::setFourthOrderArtificialViscosityParameter(const Real& a_M0sq)
{
  m_M0sq = a_M0sq;
  m_useFourthOrderArtificialViscosity = true;
}
bool GodunovPhysics::fourthOrderArtificialViscosityIsDefined() const
{
  return m_useFourthOrderArtificialViscosity;
}
Real GodunovPhysics::getFourthOrderArtificialViscosityParameter() const
{
  return m_M0sq;
}
int GodunovPhysics::densityIndex()
{
  return -1;
}
void GodunovPhysics::computeUpdate(FArrayBox&       a_dU,
                                   FluxBox&         a_F,
                                   const FArrayBox& a_U,
                                   const FluxBox&   a_WHalf,
                                   const bool&      a_useArtificialViscosity,
                                   const Real&      a_artificialViscosity,
                                   const Real&      a_currentTime,
                                   const Real&      a_dx,
                                   const Real&      a_dt,
                                   const Box&       a_box)
{
  CH_assert(isDefined());

  a_dU.setVal(0.0);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    // Get flux from WHalf
    getFlux(a_F[idir],a_WHalf[idir],idir,a_F[idir].box());

    if (a_useArtificialViscosity)
    {
      artVisc(a_F[idir],a_U,
              a_artificialViscosity,a_currentTime,
              idir,a_box);
    }

    // Compute flux difference fHi-fLo
    FArrayBox diff(a_dU.box(), a_dU.nComp());
    diff.setVal(0.0);

    FORT_FLUXDIFFF(CHF_FRA(diff),
                   CHF_CONST_FRA(a_F[idir]),
                   CHF_CONST_INT(idir),
                   CHF_BOX(a_box));

    // Add flux difference to dU
    a_dU += diff;
  }

  // Multiply dU by dt/dx because that is what the output expects
  a_dU *= -a_dt / a_dx;
}

void GodunovPhysics::getFlux(FArrayBox&       a_flux,
                             const FArrayBox& a_WHalf,
                             const int&       a_dir,
                             const Box&       a_box)
{
  MayDay::Error("GodunovPhysics::getFlux:  Default implementation called - this should never happen");
}

bool GodunovPhysics::isDefined() const
{
  return m_isDefined;
}

void GodunovPhysics::artVisc(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const Real&      a_artificialViscosity,
                             const Real&      a_currentTime,
                             const int&       a_dir,
                             const Box&       a_box)
{
  CH_assert(a_U.box().contains(a_box)); // a_box:  valid cells

  // Take the cell-centered box, a_box, and derive a face-centered box in
  // direction a_dir.
  // faceBox consists of all a_dir-faces of valid cells.
  Box faceBox = a_box;
  faceBox &= m_domain;
  faceBox.surroundingNodes(a_dir);

  CH_assert(a_F.box().contains(faceBox));

  // Derive a cell-centered box where the primitive variables are needed.
  // wBox consists of valid cells grown by 1, intersected with m_domain.
  Box wBox = faceBox;
  wBox.enclosedCells(a_dir);
  wBox.grow(1);
  wBox &= m_domain;

  // Get the primitive variables from the conserved variables (as needed).
  int numPrim = numPrimitives();
  FArrayBox W(wBox, numPrim);
  consToPrim(W, a_U, wBox); // W on valid cells + 1 layer of ghosts

  // Compute the divergence of the velocity
  FArrayBox divu(faceBox, 1);
  Interval velInt = velocityInterval();
  m_util.divVel(divu, W, velInt, a_dir, faceBox);
  // If using fourth-order artificial viscosity, apply the nonlinear operator to
  // divu.
  if (m_useFourthOrderArtificialViscosity)
    m_util.divVelHO(divu, W, a_dir, faceBox, this);

  // Change fluxes due to artificial viscosity on the interior faces
  m_util.artificialViscosity(a_F, a_U,
                             divu, a_artificialViscosity, a_dir, faceBox);

  // Change fluxes due to artificial viscosity on the boundary faces
  m_bc->artViscBC(a_F, a_U, divu, a_dir, a_currentTime);
}

void GodunovPhysics::soundSpeed(FArrayBox& a_speed,
                                const FArrayBox& a_U,
                                const Box&       a_box)
{
  MayDay::Abort("GodunovPhysics::soundSpeed function not implemented");
}

#include "NamespaceFooter.H"
