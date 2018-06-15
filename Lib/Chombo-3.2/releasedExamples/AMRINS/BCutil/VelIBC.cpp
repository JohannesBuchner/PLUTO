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

#include "SolidF_F.H"
#include "SolidVelBCF_F.H"
#include "VelIBC.H"

// Null constructor
VelIBC::VelIBC()
{
  setDefaultValues();
}

void
VelIBC::setDefaultValues()
{
  // default values
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_bcVal[dir][0] = 0;
      m_bcVal[dir][1] = 0;
      m_slopeVal[dir][0] = 0;
      m_slopeVal[dir][1] = 0;
    }
  m_isBCvalSet = false;
  m_isSlopeValSet = false;
}

// Factory method - this object is its own factory:
// Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
// its define() must be called before it is used) and m_isFortranCommonSet
// set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* VelIBC::new_physIBC()
{
  VelIBC* retval = new VelIBC();
  retval->m_velocity = m_velocity;
  retval->m_probtype = m_probtype;
  retval->m_isBCvalSet = m_isBCvalSet;
  retval->m_isSlopeValSet = m_isSlopeValSet;

  for (int dir=0; dir<SpaceDim; dir++)
    {
      retval->m_bcVal[dir][0] = m_bcVal[dir][0];
      retval->m_bcVal[dir][1] = m_bcVal[dir][1];
      retval->m_slopeVal[dir][0] = m_slopeVal[dir][0];
      retval->m_slopeVal[dir][1] = m_slopeVal[dir][1];
  }

  return static_cast<PhysIBC*>(retval);
}

// Set boundary fluxes
void VelIBC::primBC(FArrayBox&            a_WGdnv,
                    const FArrayBox&      a_Wextrap,
                    const FArrayBox&      a_W,
                    const int&            a_dir,
                    const Side::LoHiSide& a_side,
                    const Real&           a_time)
{
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
    {
      // This needs to be fixed
      // CH_assert(m_isBCvalSet);

      int lohisign;
      Box tmp = a_WGdnv.box() & m_domain;
      Real bcVal;

      // Determine which side and thus shifting directions
      if (a_side == Side::Lo)
        {
          lohisign = -1;
          bcVal = m_bcVal[a_dir][0];
        }
      else
        {
          lohisign = 1;
          bcVal = m_bcVal[a_dir][1];
        }

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
          FORT_SOLIDVELBCF(CHF_FRA(a_WGdnv),
                           CHF_CONST_FRA(a_Wextrap),
                           CHF_CONST_REAL(bcVal),
                           CHF_CONST_INT(lohisign),
                           CHF_CONST_REAL(m_dx),
                           CHF_CONST_INT(a_dir),
                           CHF_BOX(boundaryBox));
        }
    }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void VelIBC::setBdrySlopes(FArrayBox&       a_dW,
                           const FArrayBox& a_W,
                           const int&       a_dir,
                           const Real&      a_time)
{
  //  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
    {
      // This needs to be fixed
      // CH_assert(m_isSlopeValSet);

      Box loBox,hiBox,centerBox,domain;
      int hasLo,hasHi;
      Box slopeBox = a_dW.box()&m_domain;

      Real loVal = m_slopeVal[a_dir][0];
      Real hiVal = m_slopeVal[a_dir][1];

      // Generate the domain boundary boxes, loBox and hiBox, if there are
      // domain boundaries there
      loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
                 slopeBox,m_domain,a_dir);

      // Set the boundary slopes if necessary
      if ((hasLo != 0) || (hasHi != 0))
        {
          FORT_SLOPEBCSF(CHF_FRA(a_dW),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_REAL(m_dx),
                         CHF_CONST_INT(a_dir),
                         CHF_CONST_REAL(loVal),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_CONST_REAL(hiVal),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi));
        }
    }
}

// Set up initial conditions
void VelIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be in this function -- we're required to
  // define it to conform to the PhysIBC interface, but data initialization
  // is handled differently for the AMRINS code
  MayDay::Error("Shouldn't be in VelIBC::initialize");
}

/// set velocity
void
VelIBC::advectionVel(const RealVect& a_advVel)
{
  m_velocity = a_advVel;
}

///
const RealVect&
VelIBC::advectionVel() const
{
  return m_velocity;
}

// set probtype
void
VelIBC::probType(const int a_probtype)
{
  m_probtype = a_probtype;
}

//
int
VelIBC::probType() const
{
  return m_probtype;
}

void
VelIBC::setNormalWallVel(Real a_bcVal, int a_dir,
                         Side::LoHiSide a_hiLo)
{
  if (a_hiLo == Side::Lo)
    {
      m_bcVal[a_dir][0] = a_bcVal;
    }
  else
    {
      m_bcVal[a_dir][1] = a_bcVal;
    }
  // this is good enough to say that bcvals have been set
  m_isBCvalSet = true;
}

Real
VelIBC::getNormalWallVel(int a_dir, Side::LoHiSide a_hiLo) const
{
  Real bcval;
  if (a_hiLo == Side::Lo)
    {
      bcval = m_bcVal[a_dir][0];
    }
  else
    {
      bcval = m_bcVal[a_dir][1];
    }

  return bcval;
}

void
VelIBC::setSlopeValue(Real a_slopeVal, int a_dir,
                      Side::LoHiSide a_hiLo)
{
  if (a_hiLo == Side::Lo)
    {
      m_slopeVal[a_dir][0] = a_slopeVal;
    }
  else
    {
      m_slopeVal[a_dir][1] = a_slopeVal;
    }

  m_isSlopeValSet = true;
}

Real
VelIBC::getSlopeValue(int a_dir, Side::LoHiSide a_hiLo) const
{
  Real bcval;
  if (a_hiLo == Side::Lo)
    {
      bcval = m_slopeVal[a_dir][0];
    }
  else
    {
      bcval = m_slopeVal[a_dir][1];
    }

  return bcval;
}

void
VelIBC::artViscBC(FArrayBox&       a_F,
                  const FArrayBox& a_U,
                  const FArrayBox& a_divVel,
                  const int&       a_dir,
                  const Real&      a_time)
{
  MayDay::Error("VelIBC::artViscBC - not implemented");
}
