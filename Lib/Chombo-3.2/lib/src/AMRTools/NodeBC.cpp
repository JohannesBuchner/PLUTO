#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeBC.cpp
// adapted from GhostBC by DTGraves, Mon, July 19, 1999
// petermc, 13 Feb 2001

#include "NodeBC.H"
#include "NodeBCF_F.H"
#include "NamespaceHeader.H"

using std::cerr;
using std::endl;

// ---------------------------------------------------------
FaceNodeBC::FaceNodeBC(int a_dir, Side::LoHiSide a_sd)
  : m_components(-1,-1)
{
  Interval comps(0,0);
  define(a_dir, a_sd, comps);
}

// ---------------------------------------------------------
FaceNodeBC::FaceNodeBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps)
  : m_components(-1,-1)
{
  define(a_dir, a_sd, a_comps);
}

// ---------------------------------------------------------
void
FaceNodeBC::define(int a_dir, Side::LoHiSide a_sd)
{
  Interval comps(0,0);
  define(a_dir, a_sd, comps);
}


// ---------------------------------------------------------
void
FaceNodeBC::define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps)
{
  m_direction = a_dir;
  m_side = a_sd;
  m_components = a_comps;
  // m_inhomogeneous = false;
}


// ---------------------------------------------------------
void
FaceNodeBC::applyEitherBCs(FArrayBox& a_state,
                           const Box& a_domain,
                           Real a_dx,
                           bool a_homogeneous) const
{
  ProblemDomain probdomain(a_domain);
  applyEitherBCs(a_state, probdomain, a_dx, a_homogeneous);
}


// ---------------------------------------------------------
void
FaceNodeBC::applyEitherBCs(FArrayBox& a_state,
                           const ProblemDomain& a_domain,
                           Real a_dx,
                           bool a_homogeneous) const
{
  // apply bc to regular data at box boundary

  // periodic BC's trump all other BCs -- do nothing in that case
  if (a_domain.isPeriodic(m_direction)) return;

  // bx is NODE-centered
  const Box& bx = a_state.box();

  // Set domainNodes to be all the NODEs of a_domain.
  // petermc, 28 Mar 2001:  don't exclude corners.  Previously had set
  // domainNodes to be all the NODEs of a_domain except for those on
  // faces in other dimensions.
  Box domainNodes = surroundingNodes(a_domain.domainBox());
  // domainNodes.grow(-1); // shrink by 1 in all directions
  if (domainNodes.isEmpty()) return;
  // domainNodes.grow(m_direction, 1); // grow by 1 in m_direction only

  // First check (set isbc) whether a_state's box is on
  // side m_side, m_direction of a_domain.
  bool isbc;
  // endCoordinate is coordinate in m_direction on m_side
  int endCoordinate;
  switch (m_side)
    {
    case Side::Lo:
      {
        endCoordinate = domainNodes.smallEnd(m_direction);
        isbc = (bx.smallEnd(m_direction) <= endCoordinate);
        break;
      }
    case Side::Hi:
      {
        endCoordinate = domainNodes.bigEnd(m_direction);
        isbc = (endCoordinate <= bx.bigEnd(m_direction));
        break;
      }
    default:
      {
        cerr << "FaceNodeBC::applyEitherBCs():  bogus side" << endl;
        abort();
      }
    }
  // If a_state has no nodes on this face, then return.
  if (! isbc) return;

  // Set bc_box to be the nodes on this boundary face.
  Box bc_box = bx & domainNodes;
  bc_box.setRange(m_direction, endCoordinate, 1);

  // Get the arrays that describe the boundary condition.
  int ncomps = a_state.nComp();
  FArrayBox neumfac(bc_box, ncomps);
  FArrayBox dircfac(bc_box, ncomps);
  FArrayBox inhmval(bc_box, ncomps);
  fillBCValues(neumfac, dircfac, inhmval, a_dx, a_domain);
  if (a_homogeneous) inhmval.setVal(0.0);

  // Apply boundary condition.
  applyBCs(bc_box, a_state, neumfac, dircfac, inhmval, a_dx);
}


// ---------------------------------------------------------
void
FaceNodeBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                  const Box& a_domain,
                                  Real a_dx) const
{
  ProblemDomain probdomain(a_domain);
  applyEitherBCs(a_state, probdomain, a_dx, false);
}


// ---------------------------------------------------------
void
FaceNodeBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                  const ProblemDomain& a_domain,
                                  Real a_dx) const
{
  applyEitherBCs(a_state, a_domain, a_dx, false);
}


// ---------------------------------------------------------
void
FaceNodeBC::applyHomogeneousBCs(FArrayBox& a_state,
                                const Box& a_domain,
                                Real a_dx) const
{
  ProblemDomain probdomain(a_domain);
  applyEitherBCs(a_state, probdomain, a_dx, true);
}


// ---------------------------------------------------------
void
FaceNodeBC::applyHomogeneousBCs(FArrayBox& a_state,
                                const ProblemDomain& a_domain,
                                Real a_dx) const
{
  applyEitherBCs(a_state, a_domain, a_dx, true);
}


// ---------------------------------------------------------
void
FaceNodeBC::applyBCs(const Box& a_bcBox,
                     FArrayBox& a_state,
                     const FArrayBox& a_neumfac,
                     const FArrayBox& a_dircfac,
                     const FArrayBox& a_inhmval,
                     Real a_dx) const
{
#ifndef NDEBUG
  // a_state must contain nodes one level in from boundary, because
  // you need them for the Neumann condition.
  // Set biggerbox to be a_bcBox together with an adjacent layer of nodes
  // just inside the boundary, and check that a_state contains biggerbox.
  Box biggerbox = a_bcBox;
  switch (m_side)
    {
    case Side::Lo:
      {
        biggerbox.growHi(m_direction, 1);
        break;
      }
    case Side::Hi:
      {
        biggerbox.growLo(m_direction, 1);
        break;
      }
    default:
      {
        cerr << "FaceNodeBC::applyBCs(): bogus side" << endl;
        abort();
      }
    }

  CH_assert(a_state.box().contains(biggerbox));
  CH_assert(a_dircfac.box().contains(a_bcBox));
  CH_assert(a_neumfac.box().contains(a_bcBox));
  CH_assert(a_inhmval.box().contains(a_bcBox));

  CH_assert(m_components.begin() >= 0);
  CH_assert(m_components.end() < a_state.nComp());
#endif
  int iSide = sign(m_side);
  int compbegin = m_components.begin();
  int compend   = m_components.end();
  FORT_FACENODEBC(CHF_FRA(a_state),
                  CHF_CONST_FRA(a_neumfac),
                  CHF_CONST_FRA(a_dircfac),
                  CHF_CONST_FRA(a_inhmval),
                  CHF_BOX(a_bcBox),
                  CHF_CONST_INT(m_direction),
                  CHF_CONST_INT(iSide),
                  CHF_CONST_REAL(a_dx),
                  CHF_CONST_INT(compbegin),
                  CHF_CONST_INT(compend));
}


// ---------------------------------------------------------
DomainNodeBC::DomainNodeBC()
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loBC[idir] = NULL;
      m_hiBC[idir] = NULL;
    }
}


// ---------------------------------------------------------
DomainNodeBC::~DomainNodeBC()
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_loBC[idir] != NULL) delete m_loBC[idir];
      if (m_hiBC[idir] != NULL) delete m_hiBC[idir];
      m_loBC[idir] = NULL;
      m_hiBC[idir] = NULL;
    }
}


// ---------------------------------------------------------
void
DomainNodeBC::setFaceNodeBC(const FaceNodeBC& a_BC)
{
  int direction = a_BC.m_direction;
  Interval components = a_BC.m_components;
  Side::LoHiSide side = a_BC.m_side;
  CH_assert ((side == Side::Lo) || (side == Side::Hi));

  switch (side)
    {
    case Side::Lo:
      {
        if (m_loBC[direction] != NULL) delete m_loBC[direction];
        m_loBC[direction] = a_BC.new_boxBC();
        m_loBC[direction]->define(direction, side, components);
        break;
      }
    case Side::Hi:
      {
        if (m_hiBC[direction] != NULL) delete m_hiBC[direction];
        m_hiBC[direction] = a_BC.new_boxBC();
        m_hiBC[direction]->define(direction, side, components);
        break;
      }
    default:
      {
        cerr << "DomainNodeBC::setFaceNodeBC(): bogus side" << endl;
        abort();
      }
    }
}


// ---------------------------------------------------------
const FaceNodeBC&
DomainNodeBC::operator() (int a_direction, Side::LoHiSide a_side) const
{
  CH_assert ((a_side == Side::Lo) || (a_side == Side::Hi));
  FaceNodeBC* retPtr = NULL;
  switch (a_side)
    {
    case Side::Lo:
      {
        CH_assert(m_loBC[a_direction] != NULL);
        retPtr = m_loBC[a_direction];
        break;
      }
    case Side::Hi:
      {
        CH_assert(m_hiBC[a_direction] != NULL);
        retPtr = m_hiBC[a_direction];
        break;
      }
    default:
      {
        cerr << "DomainNodeBC::operator(): bogus side" << endl;
        abort();
      }
    }
  return  *retPtr;
}



// ---------------------------------------------------------
void
DomainNodeBC::applyHomogeneousBCs(NodeFArrayBox& a_state,
                                  const Box& a_domain,
                                  Real a_dx) const
{
  ProblemDomain probdomain(a_domain);
  applyHomogeneousBCs(a_state, probdomain, a_dx);
}


// ---------------------------------------------------------
void
DomainNodeBC::applyHomogeneousBCs(NodeFArrayBox& a_state,
                                  const ProblemDomain& a_domain,
                                  Real a_dx) const
{
  FArrayBox& stateFab = a_state.getFab();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // m_{lo|hi}BC[idir] is a pointer to a FaceNodeBC.
      m_loBC[idir]->applyHomogeneousBCs(stateFab, a_domain, a_dx);
      m_hiBC[idir]->applyHomogeneousBCs(stateFab, a_domain, a_dx);
    }
}


// ---------------------------------------------------------
void
DomainNodeBC::applyInhomogeneousBCs(NodeFArrayBox& a_state,
                                    const Box& a_domain,
                                    Real a_dx) const
{
  ProblemDomain probdomain(a_domain);
  applyInhomogeneousBCs(a_state, probdomain, a_dx);
}

// ---------------------------------------------------------
void
DomainNodeBC::applyInhomogeneousBCs(NodeFArrayBox& a_state,
                                    const ProblemDomain& a_domain,
                                    Real a_dx) const
{
  FArrayBox& stateFab = a_state.getFab();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // m_{lo|hi}BC[idir] is a pointer to a FaceNodeBC.
      m_loBC[idir]->applyInhomogeneousBCs(stateFab, a_domain, a_dx);
      m_hiBC[idir]->applyInhomogeneousBCs(stateFab, a_domain, a_dx);
    }
}


// ---------------------------------------------------------
DomainNodeBC::DomainNodeBC(const DomainNodeBC& a_dgbcin)
{

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loBC[idir] = NULL;
      m_hiBC[idir] = NULL;
    }
  *this = a_dgbcin;
}


// ---------------------------------------------------------
DomainNodeBC&
DomainNodeBC::operator=(const DomainNodeBC& a_dgbcin)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
        {
          if (a_dgbcin.isBCDefined(idir,sit()))
            {
              setFaceNodeBC(a_dgbcin(idir,sit()));
            }
          else
            {
              // set pointer to NULL
              resetFaceNodeBC(idir,sit());
            }
        }
    }
  return *this;
}


// ---------------------------------------------------------
bool
DomainNodeBC::isBCDefined(int a_dir, const Side::LoHiSide a_side) const
{
  bool isdefined;

  CH_assert ((a_side == Side::Lo) || (a_side == Side::Hi));

  switch (a_side)
    {
    case Side::Lo:
      {
        isdefined = (m_loBC[a_dir] != NULL);
        break;
      }
    case Side::Hi:
      {
        isdefined = (m_hiBC[a_dir] != NULL);
        break;
      }
    default:
      {
        cerr << "DomainNodeBC::isDefined(): bogus side" << endl;
        abort();
      }
    }
  return isdefined;
}


// ---------------------------------------------------------
void
DomainNodeBC::resetFaceNodeBC(const int a_dir, const Side::LoHiSide a_side)
{

  CH_assert ((a_side == Side::Lo) || (a_side == Side::Hi));

  switch (a_side)
    {
    case Side::Lo:
      {
        if (m_loBC[a_dir] != NULL) delete m_loBC[a_dir];
        m_loBC[a_dir] = NULL;
        break;
      }
    case Side::Hi:
      {
        if (m_hiBC[a_dir] != NULL) delete m_hiBC[a_dir];
        m_hiBC[a_dir] = NULL;
        break;
      }
    default:
      {
        cerr << "DomainNodeBC::resetFaceNodeBC(): bogus side" << endl;
        abort();
      }
    }
}

#include "NamespaceFooter.H"
