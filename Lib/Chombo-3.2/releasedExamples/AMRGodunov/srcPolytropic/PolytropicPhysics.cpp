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

#include "PolytropicPhysics.H"

#include "PolytropicPhysicsF_F.H"
#include "LGintegrator.H"

#include "NamespaceHeader.H"

PolytropicPhysics::PolytropicPhysics(const Real& a_smallPressure)
  :
  GodunovPhysics()
{
  m_smallPressure = a_smallPressure;
}

PolytropicPhysics::~PolytropicPhysics()
{
}

// Compute the maximum wave speed
Real PolytropicPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
                                        const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));

  Real speed = 0.0;

  FORT_MAXWAVESPEEDF(CHF_REAL(speed),
                     CHF_CONST_FRA(a_U),
                     CHF_BOX(a_box));

  return speed;
}

// Compute the speed of sound
void PolytropicPhysics::soundSpeed(FArrayBox& a_speed,
                                   const FArrayBox& a_U,
                                   const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));
  CH_assert(a_speed.contains(a_box));

  FORT_SOUNDSPEEDF(CHF_CONST_FRA1(a_speed, 0),
                   CHF_CONST_FRA(a_U),
                   CHF_BOX(a_box));
}

GodunovPhysics* PolytropicPhysics::new_godunovPhysics() const
{
  CH_assert(m_isBCSet);

  GodunovPhysics* retval = static_cast<GodunovPhysics*>
                           (new PolytropicPhysics(m_smallPressure));
  retval->setPhysIBC(m_bc);
  // added by petermc, 8 Jul 2009
  if (fourthOrderArtificialViscosityIsDefined())
    retval->setFourthOrderArtificialViscosityParameter(getFourthOrderArtificialViscosityParameter());

  return retval;
}

// Number of conserved variables
int PolytropicPhysics::numConserved()
{
  CH_assert(isDefined());

  return UNUM;
}

// Names of the conserved variables
Vector<string> PolytropicPhysics::stateNames()
{
  CH_assert(isDefined());

  Vector<string> retval;

  retval.push_back("density");

  retval.push_back("x-momentum");
  if (SpaceDim >= 2)
  {
    retval.push_back("y-momentum");

    if (SpaceDim >= 3)
    {
      retval.push_back("z-momentum");
    }
  }

  retval.push_back("energy-density");

  return retval;
}

// Number of flux variables
int PolytropicPhysics::numFluxes()
{
  CH_assert(isDefined());

  // In some computations there may be more fluxes than conserved variables
  return UNUM;
}

int PolytropicPhysics::densityIndex()
{
  return WRHO;
}

void PolytropicPhysics::getFlux(FArrayBox&       a_flux,
                                const FArrayBox& a_whalf,
                                const int&       a_dir,
                                const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_GETFLUXF(CHF_FRA(a_flux),
                CHF_CONST_FRA(a_whalf),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));
}

// Number of primitive variables
int PolytropicPhysics::numPrimitives()
{
  CH_assert(isDefined());

  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return WNUM;
}

// Compute expansion amplitudes of dW in right eigenvectors.
void PolytropicPhysics::charAnalysis(FArrayBox&       a_dW,
                                     const FArrayBox& a_W,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARANALYSISF(CHF_FRA(a_dW),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(a_box));
}

void PolytropicPhysics::charSynthesis(FArrayBox&       a_dW,
                                      const FArrayBox& a_W,
                                      const int&       a_dir,
                                      const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARSYNTHESISF(CHF_FRA(a_dW),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));
}

void PolytropicPhysics::charValues(FArrayBox&       a_lambda,
                                   const FArrayBox& a_W,
                                   const int&       a_dir,
                                   const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARVALUESF(CHF_FRA(a_lambda),
                   CHF_CONST_FRA(a_W),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(a_box));
}

// Add to (increment) the source terms given the current state
// Empty implementation - does nothing.
void PolytropicPhysics::incrementSource(FArrayBox&       a_S,
                                        const FArrayBox& a_W,
                                        const Box&       a_box)
{
}

// Compute a Riemann problem and generate fluxes at the faces.
void PolytropicPhysics::riemann(FArrayBox&       a_WGdnv,
                                const FArrayBox& a_WLeft,
                                const FArrayBox& a_WRight,
                                const FArrayBox& a_W,
                                const Real&      a_time,
                                const int&       a_dir,
                                const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_WGdnv.box().contains(a_box));

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

  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.
  FORT_RIEMANNF(CHF_FRA(a_WGdnv),
                CHF_CONST_FRA(shiftWLeft),
                CHF_CONST_FRA(shiftWRight),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

void PolytropicPhysics::postNormalPred(FArrayBox&       a_dWMinus,
                                       FArrayBox&       a_dWPlus,
                                       const FArrayBox& a_W,
                                       const Real&      a_dt,
                                       const Real&      a_dx,
                                       const int&       a_dir,
                                       const Box&       a_box)
{
  // Bound "a_dWMinus" and "a_dWPlus" so that density and pressure will never
  // be less than "smallr" and "smallp", respectively
  FORT_POSTNORMALPREDF(CHF_FRA(a_dWMinus),
                       CHF_FRA(a_dWPlus),
                       CHF_CONST_FRA(a_W),
                       CHF_BOX(a_box));
}

// Compute increment of primitive variables.
void PolytropicPhysics::quasilinearUpdate(FArrayBox&       a_dWdx,
                                          const FArrayBox& a_WHalf,
                                          const FArrayBox& a_W,
                                          const Real&      a_scale,
                                          const int&       a_dir,
                                          const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_dWdx.box().contains(a_box));

  FORT_GETADWDXF(CHF_FRA(a_dWdx),
                 CHF_CONST_FRA(a_WHalf),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_REAL(a_scale),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(a_box));
}

// Compute the primitive variables from the conserved variables
void PolytropicPhysics::consToPrim(FArrayBox&       a_W,
                                   const FArrayBox& a_U,
                                   const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  FORT_CONSTOPRIMF(CHF_FRA(a_W),
                   CHF_CONST_FRA(a_U),
                   CHF_BOX(a_box));
}

// Interval within the primitive variables corresponding to the velocities
Interval PolytropicPhysics::velocityInterval()
{
  CH_assert(isDefined());

#if CH_SPACEDIM==1
  Interval retval(WVELX,WVELX);
#elif CH_SPACEDIM==2
  Interval retval(WVELX,WVELY);
#elif CH_SPACEDIM==3
  Interval retval(WVELX,WVELZ);
#else
  // bogus_spacedim();
  Interval retval(WRHO+1, WPRES-1);
#endif

  return retval;
}

/// Interval within the flux variables corresponding to vector flux
Interval PolytropicPhysics::vectorFluxInterval()
{
  CH_assert(isDefined());

#if CH_SPACEDIM==1
  Interval retval(UMOMX,UMOMX);
#elif CH_SPACEDIM==2
  Interval retval(UMOMX,UMOMY);
#elif CH_SPACEDIM==3
  Interval retval(UMOMX,UMOMZ);
#else
  // bogus_spacedim();
  Interval retval(URHO+1, UENG-1);
#endif

  return retval;
}

// Component index within the primitive variables of the pressure
int PolytropicPhysics::pressureIndex()
{
  CH_assert(isDefined());

  return WPRES;
}

// Used to limit the absolute value of a "pressure" difference
Real PolytropicPhysics::smallPressure()
{
  CH_assert(isDefined());

  return m_smallPressure;
}

// Component index within the primitive variables of the bulk modulus
int PolytropicPhysics::bulkModulusIndex()
{
  CH_assert(isDefined());

  return WPRES;
}

#ifdef CH_USE_HDF5
void PolytropicPhysics::expressions(HDF5HeaderData& a_expressions) const
{
  a_expressions.m_string["scalar gamma"] = "1.4";

  a_expressions.m_string["vector velocity"] = "momentum/density";
  a_expressions.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  a_expressions.m_string["scalar pressure"] = "(gamma-1)*(<energy-density>-kinetic_energy*density)";
  a_expressions.m_string["scalar soundspeed"] = "sqrt(gamma*(pressure/density))";
  a_expressions.m_string["scalar log10entropy"] = "log10(pressure) - gamma*log10(<density>)";

  a_expressions.m_string["scalar x-velocity"] = "<x-momentum>/density";
  if (SpaceDim >= 2)
  {
    a_expressions.m_string["scalar y-velocity"] = "<y-momentum>/density";
    if (SpaceDim >= 3)
    {
      a_expressions.m_string["scalar z-velocity"] = "<z-momentum>/density";
    }
  }
  if (SpaceDim == 2)
    {
      a_expressions.m_string["scalar vorticity"] = "curl(velocity)";
    }
}

#include "NamespaceFooter.H"

#endif
