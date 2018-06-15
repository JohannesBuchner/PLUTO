#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Dan Martin, Fri, Jan 14, 2000

#include "MayDay.H"
#include "HOExtrapBC.H"
#include "BCF_F.H"

// -------------------------------------------------------------
HOExtrapBC::HOExtrapBC() : BoxGhostBC()
{
}

// -------------------------------------------------------------
HOExtrapBC::~HOExtrapBC()
{
}

// -------------------------------------------------------------
HOExtrapBC::HOExtrapBC(int dir, Side::LoHiSide sd,
                       const Interval& a_comps) :BoxGhostBC(dir,sd,a_comps)
{
}

// -------------------------------------------------------------
HOExtrapBC::HOExtrapBC(int dir, Side::LoHiSide sd)
  :BoxGhostBC(dir,sd)
{
}

// -------------------------------------------------------------
void
HOExtrapBC::fillBCValues(FArrayBox& a_neumfac,
                         FArrayBox& a_dircfac,
                         FArrayBox& a_inhmval,
                         Real a_dx,
                         const ProblemDomain& a_domain) const
{

  // this doesn't do anything because we don't need these arrays...
  // this function is just needed to complete the instantiation.
}

// -------------------------------------------------------------
void
HOExtrapBC::fillBCValues(FArrayBox& a_neumfac,
                         FArrayBox& a_dircfac,
                         FArrayBox& a_inhmval,
                         Real a_dx,
                         const Box& a_domain) const
{

  // this doesn't do anything because we don't need these arrays...
  // this function is just needed to complete the instantiation.
}

// -------------------------------------------------------------
void
HOExtrapBC::applyHomogeneousBCs(FArrayBox& a_state,
                                const ProblemDomain& a_domain,
                                Real a_dx) const
{
  applyHOExtrapBCs(a_state, a_domain, a_dx);
}

// -------------------------------------------------------------
void
HOExtrapBC::applyHomogeneousBCs(FArrayBox& a_state,
                              const Box& a_domain,
                              Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyHOExtrapBCs(a_state, physdomain, a_dx);
}

// -------------------------------------------------------------
void
HOExtrapBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                  const ProblemDomain& a_domain,
                                  Real a_dx) const
{
  applyHOExtrapBCs(a_state, a_domain, a_dx);
}

// -------------------------------------------------------------
void
HOExtrapBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                  const Box& a_domain,
                                  Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyHOExtrapBCs(a_state, physdomain, a_dx);
}

// -------------------------------------------------------------
void
HOExtrapBC::applyHOExtrapBCs(FArrayBox& a_state, const Box& a_domain,
                             Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyHOExtrapBCs(a_state, physdomain, a_dx);
}

// -------------------------------------------------------------
void
HOExtrapBC::applyHOExtrapBCs(FArrayBox& a_state,
                             const ProblemDomain& a_domain,
                             Real a_dx) const
{

  // only do things in non-periodic case
  if (!a_domain.isPeriodic(m_direction))
    {

      // first construct BC box

      Box bx= a_state.box();
      bx.grow(-1);
      Box bc_box;
      bool isbc;
      int idir = m_direction;
      if (m_side == Side::Lo)
        {
          isbc = (bx.smallEnd(idir) <= a_domain.domainBox().smallEnd(idir));
          if (isbc)
            {
              int ichop = a_domain.domainBox().smallEnd(m_direction);
              bc_box = a_state.box();
              bc_box.chop(m_direction, ichop);
            }
        }
      else if (m_side == Side::Hi)
        {
          isbc = (bx.bigEnd(idir) >= a_domain.domainBox().bigEnd(idir));
          if (isbc)
            {
              int ichop = a_domain.domainBox().bigEnd(m_direction)+1;
              Box chop_box = a_state.box();
              bc_box = chop_box.chop(m_direction,ichop);
            }
        }
      else
        {
          cerr << "DomainGhostBC::applyghostbc: bogus side" << endl;
          abort();
        }

      if (isbc)
      {
        int iSide = sign(m_side);
        int iDir = m_direction;

#ifndef NDEBUG
        CH_assert((m_side == Side::Hi)||(m_side == Side::Lo));
        Box biggerbox = bc_box;
        //state has to contain one cell on side of ghost cell
        if (m_side == Side::Lo)
          biggerbox.growHi(iDir, 1);
        else
          biggerbox.growLo(iDir, 1);

        CH_assert(a_state.box().contains(biggerbox));

        CH_assert(m_components.begin() >= 0);
        CH_assert(m_components.end() < a_state.nComp());
#endif

        // may need to drop order if box is not big enough
        // (assume one row of ghost cells)
        if (a_state.box().size(iDir) > 4)
        {
          int startcomp = m_components.begin();
          int endcomp = m_components.end();
          FORT_HOEXTRAPGHOSTBC(CHF_FRA(a_state),
                               CHF_BOX(bc_box),
                               CHF_CONST_INT(iDir),
                               CHF_CONST_INT(iSide),
                               CHF_CONST_REAL(a_dx),
                               CHF_CONST_INT(startcomp),
                               CHF_CONST_INT(endcomp));
        }
        else
        {
          // valid region not wide enough to apply HOExtrap -- drop
          // to linear extrap
          int startcomp = m_components.begin();
          int endcomp = m_components.end();
          FORT_EXTRAPGHOSTBC(CHF_FRA(a_state),
                             CHF_BOX(bc_box),
                             CHF_CONST_INT(iDir),
                             CHF_CONST_INT(iSide),
                             CHF_CONST_REAL(a_dx),
                             CHF_CONST_INT(startcomp),
                             CHF_CONST_INT(endcomp));
        }

      } // end if this is a bc

    }
}

// -------------------------------------------------------------
BoxGhostBC*
HOExtrapBC::new_boxghostbc() const
{
  HOExtrapBC* newop = new HOExtrapBC();
  if (newop == NULL)
  {
    MayDay::Error("Out of memory in HOExtrapBC::new_boxghostbc");
  }

  return static_cast<BoxGhostBC*>(newop);

}
