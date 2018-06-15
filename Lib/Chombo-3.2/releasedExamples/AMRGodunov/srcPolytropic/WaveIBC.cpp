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
#include "LoHiCenter.H"

#include "WaveIBC.H"
#include "WaveIBCF_F.H"
#include "NamespaceHeader.H"

// Null constructor
WaveIBC::WaveIBC()
{
  m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines
WaveIBC::WaveIBC(Real&           a_smallPressure,
                 const Real&     a_gamma,
                 const Real&     a_ambientDensity,
                 const Real&     a_deltaDensity,
                 const int&      a_pressure,
                 const IntVect&  a_waveNumber,
                 const RealVect& a_center,
                 const RealVect& a_velocity,
                 const Real&     a_artvisc)
{
  setFortranCommon(a_smallPressure,
                   a_gamma,
                   a_ambientDensity,
                   a_deltaDensity,
                   a_pressure,
                   a_waveNumber,
                   a_center,
                   a_velocity,
                   a_artvisc);
}

WaveIBC::~WaveIBC()
{
}

// Sets parameters in a common block used by Fortran routines:
//   a_smallPressure  - Lower limit for pressure (returned)
//   a_gamma          - Gamma for polytropic, gamma-law gas
//   a_ambientDensity - Ambient density add to the plane-wave
//   a_deltaDensity   - Mean of the plane-wave
//   a_pressure       - If 0, use isentropic pressure
//                      if 1, use the constant pressure of 1.0
//   a_waveNumber     - The wave number of the plane-wave
//   a_velocity       - Initial velocity of the gas
//   a_center         - Position of a maximum of the plane-wave
//   a_artvisc        - Artificial viscosity coefficient
void WaveIBC::setFortranCommon(Real&           a_smallPressure,
                               const Real&     a_gamma,
                               const Real&     a_ambientDensity,
                               const Real&     a_deltaDensity,
                               const int&      a_pressure,
                               const IntVect&  a_waveNumber,
                               const RealVect& a_center,
                               const RealVect& a_velocity,
                               const Real&     a_artvisc)
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_WAVESETF(CHF_REAL(a_smallPressure),
                CHF_CONST_REAL(a_gamma),
                CHF_CONST_REAL(a_ambientDensity),
                CHF_CONST_REAL(a_deltaDensity),
                CHF_CONST_INT(a_pressure),
                CHF_CONST_INTVECT(a_waveNumber),
                CHF_CONST_REALVECT(a_center),
                CHF_CONST_REALVECT(a_velocity),
                CHF_CONST_REAL(a_artvisc));

  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* WaveIBC::new_physIBC()
{
  WaveIBC* retval = new WaveIBC();
  retval->m_isFortranCommonSet = m_isFortranCommonSet;

  return static_cast<PhysIBC*>(retval);
}

// Set up initial conditions
void WaveIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;

    // Set up initial condition in this grid
    FORT_WAVEINITF(CHF_FRA(U),
                   CHF_CONST_REAL(m_dx),
                   CHF_BOX(uBox));
  }
}

// Set boundary fluxes
void WaveIBC::primBC(FArrayBox&            a_WGdnv,
                     const FArrayBox&      a_Wextrap,
                     const FArrayBox&      a_W,
                     const int&            a_dir,
                     const Side::LoHiSide& a_side,
                     const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    int lohisign;
    Box tmp = a_WGdnv.box();

    // Determine which side and thus shifting directions
    lohisign = sign(a_side);
    tmp.shiftHalf(a_dir,lohisign);

    // Is there a domain boundary next to this grid
    if (!m_domain.contains(tmp))
    {
      tmp &= m_domain;

      Box boundaryBox;

      // Find the strip of cells next to the domain boundary
      if (a_side == Side::Lo)
      {
        boundaryBox = bdryLo(tmp,a_dir);
      }
      else
      {
        boundaryBox = bdryHi(tmp,a_dir);
      }

      // Set the boundary fluxes
      FORT_SOLIDBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_Wextrap),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
    }
  }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void WaveIBC::setBdrySlopes(FArrayBox&       a_dW,
                            const FArrayBox& a_W,
                            const int&       a_dir,
                            const Real&      a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    Box loBox,hiBox,centerBox,domain;
    int hasLo,hasHi;
    Box slopeBox = a_dW.box();
    slopeBox.grow(a_dir,1);

    // Generate the domain boundary boxes, loBox and hiBox, if there are
    // domain boundarys there
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
               slopeBox,m_domain,a_dir);

    // Set the boundary slopes if necessary
    if ((hasLo != 0) || (hasHi != 0))
    {
      FORT_SLOPEBCSF(CHF_FRA(a_dW),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi));
    }
  }
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void WaveIBC::artViscBC(FArrayBox&       a_F,
                        const FArrayBox& a_U,
                        const FArrayBox& a_divVel,
                        const int&       a_dir,
                        const Real&      a_time)
{
}
#include "NamespaceFooter.H"
