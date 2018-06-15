#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Mon, July 19, 1999

#include "GhostBC.H"
#include "GhostBCF_F.H"
#include "NamespaceHeader.H"

BoxGhostBC::BoxGhostBC(int             a_dir,
                       Side::LoHiSide  a_sd,
                       const Interval& a_comps)
  : m_components(-1,-1)
{
  define(a_dir, a_sd, a_comps);
}

BoxGhostBC::BoxGhostBC(int            a_dir,
                       Side::LoHiSide a_sd)
  : m_components(-1,-1)
{
  Interval comps(0,0);
  define(a_dir,a_sd, comps);
}

void BoxGhostBC::define(int            a_dir,
                        Side::LoHiSide a_sd)
{
  Interval comps(0,0);
  define(a_dir, a_sd, comps);
}

void BoxGhostBC::define(int             a_dir,
                        Side::LoHiSide  a_sd,
                        const Interval& a_comps)
{
  m_direction = a_dir;
  m_side = a_sd;
  m_components = a_comps;
}

void BoxGhostBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                       const Box& a_domain,
                                       Real       a_dx) const
{
  ProblemDomain probdomain(a_domain);
  applyInhomogeneousBCs(a_state, probdomain, a_dx);
}

void BoxGhostBC::applyInhomogeneousBCs(FArrayBox&           a_state,
                                       const ProblemDomain& a_domain,
                                       Real                 a_dx) const
{
  // apply ghost bc to regular data at box boundary
  const Box& bx= a_state.box();

  // Set bcBox to be the ghost cells of a_state's box that lie
  // outside a_domain, adjacent to the face associated with this object.
  Box bcBox;
  bool isbc;
  int idir = m_direction;

  // periodic BC's trump all other BCs -- do nothing in that case
  if (!a_domain.isPeriodic(m_direction))
  {
    const Box& domBox = a_domain.domainBox();
    if (m_side == Side::Lo)
      {
        isbc = (bx.smallEnd(idir) < domBox.smallEnd(idir));
        if (isbc)
          {
            int ichop = domBox.smallEnd(m_direction);
            bcBox = a_state.box();
            bcBox.chop(m_direction, ichop);
          }
      }
    else
    if (m_side == Side::Hi)
      {
        isbc = (bx.bigEnd(idir) > domBox.bigEnd(idir));
        if (isbc)
          {
            int ichop = domBox.bigEnd(m_direction)+1;
            Box chopBox = a_state.box();
            bcBox = chopBox.chop(m_direction,ichop);
          }
      }
    else
      {
        cerr << "BoxGhostBC::applyghostbc: bogus side" << endl;
        abort();
      }

    if (isbc)
      {
        int ncomps = a_state.nComp();
        FArrayBox neumfac(bcBox, ncomps);
        FArrayBox dircfac(bcBox, ncomps);
        FArrayBox inhmval(bcBox, ncomps);

        fillBCValues(neumfac, dircfac, inhmval, a_dx, domBox);

        applyBCs(bcBox, a_state, neumfac, dircfac, inhmval, a_dx);
      }
  } // end if not periodic
}

void BoxGhostBC::applyHomogeneousBCs(FArrayBox& a_state,
                                     const Box& a_domain,
                                     Real       a_dx) const
{
  ProblemDomain probdomain(a_domain);
  applyHomogeneousBCs(a_state, probdomain, a_dx);
}

void BoxGhostBC::applyHomogeneousBCs(FArrayBox&           a_state,
                                     const ProblemDomain& a_domain,
                                     Real                 a_dx) const
{
  // apply ghost bc to regular data at box boundary
  const Box& bx= a_state.box();

  // Set bcBox to be the ghost cells of a_state's box that lie
  // outside a_domain, adjacent to the face associated with this object.
  Box bcBox;
  bool isbc;
  int idir = m_direction;

  if (!a_domain.isPeriodic(m_direction))
    {
      const Box& domBox = a_domain.domainBox();
      if (m_side == Side::Lo)
        {
          isbc = (bx.smallEnd(idir) <domBox.smallEnd(idir));
          if (isbc)
            {
              int ichop = domBox.smallEnd(m_direction);
              bcBox = a_state.box();
              bcBox.chop(m_direction, ichop);
            }
        }
      else
      if (m_side == Side::Hi)
        {
          isbc = (bx.bigEnd(idir) >domBox.bigEnd(idir));
          if (isbc)
            {
              int ichop = domBox.bigEnd(m_direction)+1;
              Box chopBox = a_state.box();
              bcBox = chopBox.chop(m_direction,ichop);
            }
        }
      else
        {
          cerr << "BoxGhostBC::applyghostbc: bogus side" << endl;
          abort();
        }

      if (isbc)
        {
          int ncomps = a_state.nComp();
          FArrayBox neumfac(bcBox, ncomps);
          FArrayBox dircfac(bcBox, ncomps);
          FArrayBox inhmval(bcBox, ncomps);

          fillBCValues(neumfac, dircfac, inhmval, a_dx, domBox);
          inhmval.setVal(0.0);

          applyBCs(bcBox, a_state, neumfac, dircfac, inhmval, a_dx);
        }
    }
}

void BoxGhostBC::applyBCs(const Box&       a_bcBox,
                          FArrayBox&       a_state,
                          const FArrayBox& a_neumfac,
                          const FArrayBox& a_dircfac,
                          const FArrayBox& a_inhmval,
                          Real             a_dx) const
{
  int iSide = sign(m_side);
  int iDir = m_direction;

  CH_assert((m_side == Side::Hi) || (m_side == Side::Lo));

  Box biggerbox = a_bcBox;

  // state has to contain one cell on side of ghost cell
  if (m_side == Side::Lo)
    {
      biggerbox.growHi(iDir, 1);
    }
  else
    {
      biggerbox.growLo(iDir, 1);
    }

  // in order to deal with nghost > 1, need some measure of
  // where interior domain is. for now, just look at cells
  // adjacent to bcBox
  Box interiorBox;
  if (m_side == Side::Lo)
    {
      interiorBox = adjCellHi(a_bcBox,iDir,2);
    }
  else
    {
      interiorBox = adjCellLo(a_bcBox,iDir,2);
    }

  CH_assert(a_state.box().contains(biggerbox));
  CH_assert(a_dircfac.box().contains(a_bcBox));
  CH_assert(a_neumfac.box().contains(a_bcBox));
  CH_assert(a_inhmval.box().contains(a_bcBox));
  CH_assert(m_components.begin() >= 0);
  CH_assert(m_components.end() < a_state.nComp());

  int compbegin = m_components.begin();
  int compend   = m_components.end();

  FORT_BOXGHOSTBC(CHF_FRA(a_state),
                  CHF_CONST_FRA(a_neumfac),
                  CHF_CONST_FRA(a_dircfac),
                  CHF_CONST_FRA(a_inhmval),
                  CHF_BOX(a_bcBox),
                  CHF_BOX(interiorBox),
                  CHF_CONST_INT(iDir),
                  CHF_CONST_INT(iSide),
                  CHF_CONST_REAL(a_dx),
                  CHF_CONST_INT(compbegin),
                  CHF_CONST_INT(compend));
}

DomainGhostBC::DomainGhostBC()
{

}

DomainGhostBC::~DomainGhostBC()
{
  clear();
}


void
DomainGhostBC::setBoxGhostBC(const BoxGhostBC& a_ghostBC)
{
  int direction = a_ghostBC.m_direction;
  Interval components = a_ghostBC.m_components;
  Side::LoHiSide side = a_ghostBC.m_side;

  CH_assert ((side == Side::Lo) || (side == Side::Hi));

  BoxGhostBC* newbc = a_ghostBC.new_boxghostbc();
  newbc->define(direction,side, components);

  switch (side)
    {
      case Side::Lo:
        {
          m_loGhostBC[direction].append(newbc);

          break;
        }

      case Side::Hi:
        {
          m_hiGhostBC[direction].append(newbc);

          break;
        }

      default:
        {
          MayDay::Error("DomainGhostBC::setboxghostbc: bogus side");
        }
    }
}

const List<BoxGhostBC*>& DomainGhostBC::operator()(int            a_direction,
                                                   Side::LoHiSide a_side) const
{
  CH_assert ((a_side == Side::Lo) || (a_side == Side::Hi));

  if (a_side == Side::Lo)
    {
      return  m_loGhostBC[a_direction];
    }
  else
  if (a_side == Side::Hi)
    {
      return m_hiGhostBC[a_direction];
    }
  else
    {
      MayDay::Error("DomainGhostBC::operator(): bogus side");
    }

  // we should never reach this point -- just putting this in
  // to avoid compiler warning
  MayDay::Error("DomainGhostBC::operator(): reached inaccessible point");
  return m_loGhostBC[0];

}

void DomainGhostBC::applyHomogeneousBCs(FArrayBox& a_state,
                                        const Box& a_domain,
                                        Real       a_dx) const
{
  ProblemDomain probdomain(a_domain);

  applyHomogeneousBCs(a_state, probdomain, a_dx);
}

void DomainGhostBC::applyHomogeneousBCs(FArrayBox&           a_state,
                                        const ProblemDomain& a_domain,
                                        Real                 a_dx) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ListIterator<BoxGhostBC*> litLo(m_loGhostBC[idir]);
      for (litLo.rewind(); litLo.ok(); ++litLo)
        {
          m_loGhostBC[idir][litLo]->applyHomogeneousBCs(a_state,
                                                        a_domain, a_dx);
        }
      ListIterator<BoxGhostBC*> litHi(m_hiGhostBC[idir]);
      for (litHi.rewind(); litHi.ok(); ++litHi)
        {
          m_hiGhostBC[idir][litHi]->applyHomogeneousBCs(a_state,
                                                        a_domain, a_dx);
        }
    }
}

void DomainGhostBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                          const Box& a_domain,
                                          Real       a_dx) const
{
  ProblemDomain probdomain(a_domain);

  applyInhomogeneousBCs(a_state, probdomain,a_dx);
}

void DomainGhostBC::applyInhomogeneousBCs(FArrayBox&           a_state,
                                          const ProblemDomain& a_domain,
                                          Real                 a_dx) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ListIterator<BoxGhostBC*> litLo(m_loGhostBC[idir]);
      for (litLo.rewind(); litLo.ok(); ++litLo)
        {
          m_loGhostBC[idir][litLo]->applyInhomogeneousBCs(a_state,
                                                          a_domain, a_dx);
        }
      ListIterator<BoxGhostBC*> litHi(m_hiGhostBC[idir]);
      for (litHi.rewind(); litHi.ok(); ++litHi)
        {
          m_hiGhostBC[idir][litHi]->applyInhomogeneousBCs(a_state,
                                                          a_domain, a_dx);
        }

    }
}

DomainGhostBC::DomainGhostBC(const DomainGhostBC& a_dgbcin)
{
  *this = a_dgbcin;
}

DomainGhostBC& DomainGhostBC::operator=(const DomainGhostBC& a_dgbcin)
{
  // first clear ourselves, just in case
  clear();

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
        {
          if (a_dgbcin.isBCDefined(idir,sit()))
            {
              const List<BoxGhostBC*>& thisList = a_dgbcin(idir, sit());
              ListIterator<BoxGhostBC*> lit(thisList);
              for (lit.rewind(); lit.ok(); ++lit)
                {
                  CH_assert (thisList[lit] != NULL);
                  setBoxGhostBC(*thisList[lit]);
                }
            }
          else
            {
              // set pointer to NULL
              resetBoxGhostBC(idir,sit());
            }
        }
    }

  return *this;
}

bool DomainGhostBC::isBCDefined(const int            a_dir,
                                const Side::LoHiSide a_side) const
{
  bool isdefined;

  CH_assert ((a_side == Side::Lo) || (a_side == Side::Hi));

  if (a_side == Side::Lo)
    {
      isdefined = (m_loGhostBC[a_dir].isNotEmpty());
    }
  else
  if (a_side == Side::Hi)
    {
      isdefined = (m_hiGhostBC[a_dir].isNotEmpty());
    }
  else
    {
      MayDay::Error("DomainGhostBC::isDefined(): bogus side");
    }

  return isdefined;
}

void DomainGhostBC::resetBoxGhostBC(const int            a_dir,
                                    const Side::LoHiSide a_side)
{
  CH_assert ((a_side == Side::Lo) || (a_side == Side::Hi));

  if (a_side == Side::Lo)
    {
      if (m_loGhostBC[a_dir].isNotEmpty())
        {
          ListIterator<BoxGhostBC*> lit(m_loGhostBC[a_dir]);
          for (lit.rewind(); lit.ok(); ++lit)
            {
              if (m_loGhostBC[a_dir][lit] != NULL)
                {
                  delete m_loGhostBC[a_dir][lit];
                  m_loGhostBC[a_dir][lit] = NULL;
                }
            } // end loop over list

          m_loGhostBC[a_dir].clear();
        } // end if list isn't empty
    } // end lo side
  else
    if (a_side == Side::Hi)
    {
      if (m_hiGhostBC[a_dir].isNotEmpty())
        {
          ListIterator<BoxGhostBC*> lit(m_hiGhostBC[a_dir]);
          for (lit.rewind(); lit.ok(); ++lit)
            {
              if (m_hiGhostBC[a_dir][lit] != NULL)
                {
                  delete m_hiGhostBC[a_dir][lit];
                  m_hiGhostBC[a_dir][lit] = NULL;
                }
            } // end hiop over list

          m_hiGhostBC[a_dir].clear();
        } // end if list isn't empty
    }
  else
    {
      MayDay::Error("DomainGhostBC::resetBoxGhostBC(): bogus side");
    }
}


void
DomainGhostBC::clear()
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ListIterator<BoxGhostBC*> litLo(m_loGhostBC[idir]);
      for (litLo.begin(); litLo.ok(); ++litLo)
        {
          if (m_loGhostBC[idir][litLo] != NULL)
            {
              delete m_loGhostBC[idir][litLo];
              m_loGhostBC[idir][litLo] = NULL;
            }
        } // end pass through items in list
      // now that we've cleaned up to avoid memory leaks, can clear list.
      m_loGhostBC[idir].clear();

      ListIterator<BoxGhostBC*> litHi(m_hiGhostBC[idir]);
      for (litHi.begin(); litHi.ok(); ++litHi)
        {
          if (m_hiGhostBC[idir][litHi] != NULL)
            {
              delete m_hiGhostBC[idir][litHi];
              m_hiGhostBC[idir][litHi] = NULL;
            }
        } // end pass through items in list
      // now that we've cleaned up to avoid memory leaks, can clear list.
      m_hiGhostBC[idir].clear();
    }
}

#include "NamespaceFooter.H"
