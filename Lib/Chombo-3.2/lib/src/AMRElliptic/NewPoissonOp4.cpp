#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewPoissonOp4.H"
#include "NewPoissonOp4F_F.H"
#include "FORT_PROTO.H"
#include "AMRPoissonOpF_F.H"
#include "NamespaceHeader.H"

void NewPoissonOp4::define( const RealVect& a_dx,
                           const ProblemDomain& a_domain,
                           BCFunc a_bc)
{
  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
}

void NewPoissonOp4::residual(  FArrayBox& a_lhs, const FArrayBox& a_phi,
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
void NewPoissonOp4::preCond(   FArrayBox& a_phi, const FArrayBox& a_rhs)
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

void NewPoissonOp4::applyOp(   FArrayBox& a_lhs, const FArrayBox& a_phi,
                              bool a_homogeneous )
{

  FArrayBox& phi = (FArrayBox&)a_phi;
  Real dx = m_dx[0];

  m_bc(phi, m_domain.domainBox(), m_domain, dx, a_homogeneous);

  Real alpha = 0.0;
  Real beta = 1.0;
  const Box& region = a_lhs.box();
  FORT_OPERATORLAP4(CHF_FRA(a_lhs),
                    CHF_CONST_FRA(phi),
                    CHF_BOX(region),
                    CHF_CONST_REAL(dx),
                    CHF_CONST_REAL(alpha),
                    CHF_CONST_REAL(beta));

}

void NewPoissonOp4::create(    FArrayBox& a_lhs, const FArrayBox& a_rhs)
{
  a_lhs.define(a_rhs.box(), a_rhs.nComp());
}

void NewPoissonOp4::assign(    FArrayBox& a_lhs, const FArrayBox& a_rhs)
{
  a_lhs.copy(a_rhs);
}

Real NewPoissonOp4::dotProduct(const FArrayBox& a_1, const FArrayBox& a_2)
{
  return a_1.dotProduct(a_2);
}

void NewPoissonOp4::incr( FArrayBox& a_lhs, const FArrayBox& a_x, Real a_scale)
{
  a_lhs.plus(a_x, a_scale);
}

void NewPoissonOp4::axby( FArrayBox& a_lhs, const FArrayBox& a_x,
                         const FArrayBox& a_y, Real a_a, Real a_b)
{
  a_lhs.copy(a_x);
  a_lhs*=a_a;
  a_lhs.plus(a_y, a_b);
}

void NewPoissonOp4::scale(FArrayBox& a_lhs, const Real& a_scale)
{
  a_lhs*=a_scale;
}

void NewPoissonOp4::setToZero(FArrayBox& a_lhs)
{
  a_lhs.setVal(0.0);
}

Real NewPoissonOp4::norm(const FArrayBox& a_x, int a_ord)
{
  return a_x.norm(a_ord);
}

void NewPoissonOp4::relax(FArrayBox& a_e,
                         const FArrayBox& a_residual,
                         int a_iterations)
{
  for (int i=0; i<a_iterations; i++)
  {
    levelGSRB(a_e, a_residual);
  }
}

void NewPoissonOp4::createCoarser(FArrayBox& a_coarse,
                                 const FArrayBox& a_fine,
                                 bool ghosted)
{
  Box c = a_fine.box();
  // check the box
  if (ghosted) c.grow(-m_nGhost);
  c.coarsen(2);
  c.refine(2);
  if (ghosted) c.grow( m_nGhost);
  CH_assert(c == a_fine.box());
  // do the coarsening
  if (ghosted) c.grow(-m_nGhost);
  c.coarsen(2);
  if (ghosted) c.grow( m_nGhost);
  a_coarse.define(c, a_fine.nComp());
}

void NewPoissonOp4::restrictResidual(FArrayBox& a_resCoarse,
                                    FArrayBox& a_phiFine,
                                    const FArrayBox& a_rhsFine)
{
  m_bc(a_phiFine, m_domain.domainBox(), m_domain, m_dx[0], true);

  a_resCoarse.setVal(0.0);
  Real alpha = 0.0, beta = 1.0;
  FORT_RESTRICTRES4(CHF_FRA(a_resCoarse),
                    CHF_CONST_FRA(a_phiFine),
                    CHF_CONST_FRA(a_rhsFine),
                    CHF_CONST_REAL(alpha),
                    CHF_CONST_REAL(beta),
                    CHF_BOX(a_rhsFine.box()),
                    CHF_CONST_REAL(m_dx[0]));
}

void NewPoissonOp4::prolongIncrement(FArrayBox& a_phiThisLevel,
                                    const FArrayBox& a_correctCoarse)
{
  const int r = 2;

  FArrayBox& phi =  a_phiThisLevel;
  const FArrayBox& coarse = a_correctCoarse;
  //FArrayBox& c = (FArrayBox&)coarse;
  Box region = a_phiThisLevel.box();
  region.grow(-m_nGhost);
  Box cBox = a_correctCoarse.box();
  cBox.grow(-m_nGhost);
  CH_assert(cBox == coarsen(region, r));

  int prolongOrder = 0;
  if (prolongOrder==0)
  {
    FORT_PROLONG(CHF_FRA(phi),
                 CHF_CONST_FRA(coarse),
                 CHF_BOX(region),
                 CHF_CONST_INT(r));
  }
  else if (prolongOrder==1)
  {
    // set boundary condition for coarse-level correction
    FORT_PROLONGLINEAR(CHF_FRA(phi),
                       CHF_CONST_FRA(coarse),
                       CHF_BOX(region),
                       CHF_BOX(cBox),
                       CHF_CONST_INT(r));
  }
}

/***/
void NewPoissonOp4::
levelGSRB(FArrayBox&       a_phi,
          const FArrayBox& a_rhs)
{
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  Real dx = m_dx[0];
  FArrayBox tmp;
  create(tmp, a_phi);
  assign(tmp, a_phi);
  for (int whichpass=0; whichpass<2; whichpass++)
    {
      if (whichpass==0)
        m_bc(a_phi,  m_domain.domainBox(), m_domain,  dx, true);
      else
        m_bc(tmp,  m_domain.domainBox(), m_domain,  dx, true);
      FORT_GSRBLAPLACIAN4(CHF_FRA(a_phi),
                          CHF_CONST_FRA(a_rhs),
                          CHF_BOX(a_rhs.box()),
                          CHF_CONST_REAL(dx),
                          CHF_FRA(tmp),
                          CHF_CONST_INT(whichpass));
    }
}

NewPoissonOp4Factory::NewPoissonOp4Factory()
{
  m_dx  = RealVect::Unit;
  m_bc  = NULL;
}

NewPoissonOp4Factory::NewPoissonOp4Factory(RealVect& a_dx, BCFunc a_bc)
  :m_dx(a_dx), m_bc(a_bc)
{
}

void NewPoissonOp4Factory::define(const RealVect& a_dx,  BCFunc a_bc)
{
  m_dx=a_dx;
  m_bc=a_bc;
}

void NewPoissonOp4::
createCoarsened(FArrayBox&                  a_lhs,
                const FArrayBox&            a_rhs,
                const int &                 a_refRat)
{
  int ncomp = a_rhs.nComp();
  Box coarsenedBox = coarsen(a_rhs.box(), a_refRat);
  a_lhs.define(coarsenedBox, ncomp);
}

NewPoissonOp4* NewPoissonOp4Factory::MGnewOp(const ProblemDomain& a_FineindexSpace,
                                           int   a_depth,
                                           bool  a_homoOnly)
{
  CH_assert(a_depth >= 0 );
  CH_assert(m_bc != NULL);
  NewPoissonOp4* newOp = new NewPoissonOp4();
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

void NewPoissonOp4Factory::MGreclaim(NewPoissonOp4* a_reclaim)
{
  delete a_reclaim;
}

#include "NamespaceFooter.H"
