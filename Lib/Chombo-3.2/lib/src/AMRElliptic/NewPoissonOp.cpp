#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewPoissonOp.H"
#include "FORT_PROTO.H"
#include "AMRPoissonOpF_F.H"
#include "NamespaceHeader.H"

void NewPoissonOp::define( const RealVect& a_dx,
                           const ProblemDomain& a_domain,
                           BCFunc a_bc)
{
  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
}

void NewPoissonOp::residual(  FArrayBox& a_lhs, const FArrayBox& a_phi,
                              const FArrayBox& a_rhs, bool a_homogeneous)
{
  applyOp(a_lhs, a_phi, a_homogeneous);
  a_lhs*=-1.0;
  incr(a_lhs, a_rhs, 1);
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void NewPoissonOp::preCond(   FArrayBox& a_phi, const FArrayBox& a_rhs)
{
    // diagonal term of this operator is 4/h/h in 2D, 6/h/h in 3D,
  // so inverse of this is our initial multiplier
  Real mult = -m_dx[0]*m_dx[0]/(2.0*SpaceDim);
  Interval comps = a_phi.interval();
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  a_phi.copy(a_rhs.box(), comps,
             a_rhs.box(), a_rhs,
             comps);

  a_phi *= mult;

  relax(a_phi, a_rhs, 2);
}

void NewPoissonOp::applyOp(   FArrayBox& a_lhs, const FArrayBox& a_phi,
                              bool a_homogeneous )
{

  FArrayBox& phi = (FArrayBox&)a_phi;
  Real dx = m_dx[0];

  m_bc(phi, m_domain.domainBox(), m_domain, dx, a_homogeneous);

  Real alpha = 0.0;
  Real beta = 1.0;
  const Box& region = a_lhs.box();
  FORT_OPERATORLAP(CHF_FRA(a_lhs),
                   CHF_CONST_FRA(phi),
                   CHF_BOX(region),
                   CHF_CONST_REAL(dx),
                   CHF_CONST_REAL(alpha),
                   CHF_CONST_REAL(beta));

}

void NewPoissonOp::create(    FArrayBox& a_lhs, const FArrayBox& a_rhs)
{
  a_lhs.define(a_rhs.box(), a_rhs.nComp());
}

void NewPoissonOp::assign(    FArrayBox& a_lhs, const FArrayBox& a_rhs)
{
  a_lhs.copy(a_rhs);
}

Real NewPoissonOp::dotProduct(const FArrayBox& a_1, const FArrayBox& a_2)
{
  return a_1.dotProduct(a_2);
}
void NewPoissonOp::incr( FArrayBox& a_lhs, const FArrayBox& a_x, Real a_scale)
{
  a_lhs.plus(a_x, a_scale);
}
void NewPoissonOp::axby( FArrayBox& a_lhs, const FArrayBox& a_x,
                         const FArrayBox& a_y, Real a_a, Real a_b)
{
  a_lhs.copy(a_x);
  a_lhs*=a_a;
  a_lhs.plus(a_y, a_b);
}
void NewPoissonOp::scale(FArrayBox& a_lhs, const Real& a_scale)
{
  a_lhs*=a_scale;
}
void NewPoissonOp::setToZero(FArrayBox& a_lhs)
{
  a_lhs.setVal(0.0);
}

Real NewPoissonOp::norm(const FArrayBox& a_x, int a_ord)
{
  return a_x.norm(a_ord);
}

void NewPoissonOp::relax(FArrayBox& a_e,
                         const FArrayBox& a_residual,
                         int a_iterations)
{
  for (int i=0; i<a_iterations; i++)
  {
    levelGSRB(a_e, a_residual);
  }
}

void NewPoissonOp::createCoarser(FArrayBox& a_coarse,
                                 const FArrayBox& a_fine,
                                 bool ghosted)
{
  Box c = a_fine.box();
  if (ghosted)
    {
      c.grow(-1);
    }
  c.coarsen(2);
  c.refine(2);
  if (ghosted) c.grow(1);
  CH_assert(c == a_fine.box());
  if (ghosted) c.grow(-1);
  c.coarsen(2);
  if (ghosted)
    {
      c.grow(1);
    }
  a_coarse.define(c, a_fine.nComp());

}

void NewPoissonOp::restrictResidual(FArrayBox& a_resCoarse,
                                    FArrayBox& a_phiFine,
                                    const FArrayBox& a_rhsFine)
{

  Real dx = m_dx[0];

  FArrayBox& phi = a_phiFine;
  m_bc(phi,  m_domain.domainBox(), m_domain,  dx, true);

  const FArrayBox& rhs = a_rhsFine;
  FArrayBox& res = a_resCoarse;

  Box region = rhs.box();

  res.setVal(0.0);
  Real alpha = 0;
  Real beta = 1.0;
  FORT_RESTRICTRES(CHF_FRA(res),
                   CHF_CONST_FRA(phi),
                   CHF_CONST_FRA(rhs),
                   CHF_CONST_REAL(alpha),
                   CHF_CONST_REAL(beta),
                   CHF_BOX(region),
                   CHF_CONST_REAL(dx));

}

void NewPoissonOp::prolongIncrement(FArrayBox& a_phiThisLevel,
                                    const FArrayBox& a_correctCoarse)
{

  FArrayBox& phi =  a_phiThisLevel;
  const FArrayBox& coarse = a_correctCoarse;
  //FArrayBox& c = (FArrayBox&)coarse;
  Box region = a_phiThisLevel.box();
  region.grow(-1);
  Box cBox = a_correctCoarse.box();
  cBox.grow(-1);
  CH_assert(cBox == coarsen(region, 2));

  int r = 2;

  FORT_PROLONG(CHF_FRA(phi),
               CHF_CONST_FRA(coarse),
               CHF_BOX(region),
               CHF_CONST_INT(r));
}

/***/
void NewPoissonOp::
colorGS(FArrayBox&       a_phi,
        const FArrayBox& a_rhs,
        const IntVect&   color)
{
  bool homogeneous = true;

  Real weight = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0  / (m_dx[idir]*m_dx[idir]);
    }
  weight = 1.0 / weight;

  levelRelaxColor(a_phi, a_rhs, color, weight, homogeneous);
}

/***/
bool NewPoissonOp::
nextColor(IntVect& color,
          const IntVect& limit)
{
  color[0]++;

  for (int i=0; i<CH_SPACEDIM-1; ++i)
    {
      if (color[i] > limit[i])
        {
          color[i] = 0;
          color[i+1]++;
        }
    }
  if (color[CH_SPACEDIM-1] > limit[CH_SPACEDIM-1])
    return false;

  return true;
}
/***/
void NewPoissonOp::
levelGSRB(FArrayBox&       a_phi,
          const FArrayBox& a_rhs)
{
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  Real dx = m_dx[0];

  // do first red, then black passes
  for (int whichPass =0; whichPass <= 1; whichPass++)
    {

      m_bc(a_phi,  m_domain.domainBox(), m_domain,  dx, true);

      // now step through grids...
      //fill in intersection of ghostcells and a_phi's boxes
      // dfm -- for a 5 point stencil, this should not be necessary
      //a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      // invoke physical BC's where necessary
      //m_bc(a_phi,  m_domain,  dx, true); //new approach to help with checker
#ifndef NDEBUG

#endif

      FORT_GSRBLAPLACIAN(CHF_FRA(a_phi),
                         CHF_CONST_FRA(a_rhs),
                         CHF_BOX(a_rhs.box()),
                         CHF_CONST_REAL(dx),
                         CHF_CONST_INT(whichPass));

    } // end loop through red-black
}

void NewPoissonOp::
levelRelaxColor(FArrayBox&       a_phi,
                const FArrayBox& a_rhs,
                const IntVect&   a_color,
                const Real&      a_weight,
                const bool&      a_homogeneousPhysBC)
{
  CH_assert(a_phi.nComp() == 1);

  Box dblBox = a_rhs.box();

  //apply domain boundary condtions
  m_bc(a_phi, m_domain.domainBox(), m_domain, m_dx[0], a_homogeneousPhysBC);

  //find box over which we are iterating
  IntVect loIV = dblBox.smallEnd();
  IntVect hiIV = dblBox.bigEnd();

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (loIV[idir] % 2 != a_color[idir])
        {
          loIV[idir]++;
        }
    }

  if (loIV <= hiIV)
    {
      Box coloredBox(loIV, hiIV);
      int comp = 0;
      Real alpha = 0;
      Real beta = 1.;
      FORT_AMRPMULTICOLOR(CHF_FRA1(a_phi,comp),
                          CHF_CONST_FRA1(a_rhs,comp),
                          CHF_CONST_REAL(a_weight),
                          CHF_CONST_REAL(alpha),
                          CHF_CONST_REAL(beta),
                          CHF_CONST_REALVECT(m_dx),
                          CHF_BOX(coloredBox));
    }
}

NewPoissonOpFactory::NewPoissonOpFactory()
{
  m_dx  = RealVect::Unit;
  m_bc  = NULL;
}

NewPoissonOpFactory::NewPoissonOpFactory(RealVect& a_dx, BCFunc a_bc)
  :m_dx(a_dx), m_bc(a_bc)
{
}

void NewPoissonOpFactory::define(const RealVect& a_dx,  BCFunc a_bc)
{
  m_dx=a_dx;
  m_bc=a_bc;
}

void NewPoissonOp::
createCoarsened(FArrayBox&                  a_lhs,
                const FArrayBox&            a_rhs,
                const int &                 a_refRat)
{
  int ncomp = a_rhs.nComp();
  //fill ebislayout
  Box coarsenedBox = coarsen(a_rhs.box(), a_refRat);
  a_lhs.define(coarsenedBox, ncomp);
}
NewPoissonOp* NewPoissonOpFactory::MGnewOp(const ProblemDomain& a_FineindexSpace,
                                           int   a_depth,
                                           bool  a_homoOnly)
{
  CH_assert(a_depth >= 0 );
  CH_assert(m_bc != NULL);
  NewPoissonOp* newOp = new NewPoissonOp();
  RealVect dx = m_dx;
  ProblemDomain domain = a_FineindexSpace;
  for (int i=0; i<a_depth; i++)
    {
      Box d = domain.domainBox();
      d.coarsen(8);
      d.refine(8);
      if (domain.domainBox() == d)
        {
          dx*=2;
          domain.coarsen(2);
        }
      else
        {
          return NULL;
        }
    }
  newOp->define(dx, domain, m_bc);
  return newOp;
}

void NewPoissonOpFactory::MGreclaim(NewPoissonOp* a_reclaim)
{
  delete a_reclaim;
}
#include "NamespaceFooter.H"
