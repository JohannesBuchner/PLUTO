#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTG 7-31-2001

#include "EBFluxFAB.H"
#include "NamespaceHeader.H"

void
EBFluxFAB::
clone(const EBFluxFAB& a_input)
{
  define(a_input.m_ebisBox, a_input.m_region, a_input.m_nComp);
  this->setVal(0.);
  (*this) += a_input;
}

void
EBFluxFAB::
alias(Vector<EBFaceFAB*> a_inputFAB)
{
  clear();
  m_nComp = a_inputFAB[0]->nComp();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_fluxes[idir] = a_inputFAB[idir];
      CH_assert(m_nComp == a_inputFAB[idir]->nComp());
    }
  m_aliased  = true;
}
int
EBFluxFAB::
size(const Box& R, const Interval& comps) const
{
  int retval = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      retval += m_fluxes[idir]->size(R, comps);
    }

  return retval;
}

void
EBFluxFAB::
linearOut(void* buf, const Box& R, const Interval& comps) const
{
  unsigned char* charbuf = (unsigned char*) buf;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_fluxes[idir]->linearOut(charbuf, R, comps);
      int dirsize =  m_fluxes[idir]->size(R, comps);
      charbuf += dirsize;
    }
}

void
EBFluxFAB::
linearIn(void* buf, const Box& R, const Interval& comps)
{
  unsigned char* charbuf = (unsigned char*) buf;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_fluxes[idir]->linearIn(charbuf, R, comps);
      int dirsize =  m_fluxes[idir]->size(R, comps);
      charbuf += dirsize;
    }
}
// first do simple access functions
// ---------------------------------------------------------
int
EBFluxFAB::
nComp() const
{
  return m_nComp;
}

// ---------------------------------------------------------
const Box&
EBFluxFAB::getRegion() const
{
  return m_region;
}

// ---------------------------------------------------------
EBFaceFAB&
EBFluxFAB::operator[] (const int dir)
{
  CH_assert(m_nComp >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_fluxes[dir] != NULL);

  return *m_fluxes[dir];
}

// ---------------------------------------------------------
const EBFaceFAB&
EBFluxFAB::operator[] (const int dir)  const
{
  CH_assert(m_nComp >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_fluxes[dir] != NULL);

  return *m_fluxes[dir];
}

// ---------------------------------------------------------
// constructors and destructors
// ---------------------------------------------------------
EBFluxFAB::EBFluxFAB()
{
  setDefaultValues();
}

// ---------------------------------------------------------
EBFluxFAB::EBFluxFAB(const EBISBox& a_ebisBox,
                     const Box& a_region, int a_nComp)
{
  setDefaultValues();
  define(a_ebisBox, a_region, a_nComp);
}
void
EBFluxFAB::define(const EBISBox& a_ebisBox,
                  const Box& a_region, int a_nComp)
{
  clear();
  m_isDefined = true;
  m_aliased = false;
  m_ebisBox = a_ebisBox;
  m_region = a_region & a_ebisBox.getRegion();
  m_nComp = a_nComp;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_fluxes[idir] != NULL)
        {
          delete m_fluxes[idir];
          m_fluxes[idir] = NULL;
        }
      m_fluxes[idir] = new EBFaceFAB(a_ebisBox, a_region, idir, a_nComp);
    }
}
// ---------------------------------------------------------
bool
EBFluxFAB::isDefined() const
{
  return m_isDefined;
}
// ---------------------------------------------------------
EBFluxFAB::~EBFluxFAB()
{

  clear();
}

// ---------------------------------------------------------
void
EBFluxFAB::clear()
{
  // first delete storage
  if (!m_aliased)
    {
      for (int dir=0; dir < SpaceDim; dir++)
        {
          if (m_fluxes[dir] != NULL)
            {
              delete m_fluxes[dir];
            }
        }
    }

  // now reset all other variables
  setDefaultValues();
}

void
EBFluxFAB::setDefaultValues()
{
  for (int dir=0; dir < SpaceDim; dir++)
    {
      m_fluxes[dir] = NULL;
    }

  // set the box to the empty box...
  m_region = Box();
  m_aliased = false;
}

// ---------------------------------------------------------
void
EBFluxFAB::setVal(const Real& val)
{
  CH_assert(m_nComp > 0);

  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->setVal(val);
    }
}

// ---------------------------------------------------------
void
EBFluxFAB::copy(const Box& Rfrom, const Interval& Cdest,
                const Box& Rto, const EBFluxFAB& src,
                const Interval& Csrc)
{
  //CH_assert(Rfrom == Rto);
  //Box R = Rto;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      const EBFaceFAB& srcFab = src[dir];
      m_fluxes[dir]->copy(Rfrom,Cdest,Rto, srcFab,Csrc);
    }

}
#include "NamespaceFooter.H"
