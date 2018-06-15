#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*****************/
/*****************/

#include "VelocityAMRPoissonOp.H"
#include "AMRPoissonOpF_F.H"

/***********************/
// default constructor
/***********************/
// ---------------------------------------------------------
VelocityAMRPoissonOp::VelocityAMRPoissonOp() : AMRPoissonOp()
{
}

// destructor
// ---------------------------------------------------------
VelocityAMRPoissonOp::~VelocityAMRPoissonOp()
{
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
/***********************/
// ---------------------------------------------------------
void
VelocityAMRPoissonOp::applyOpI(LevelData<FArrayBox>&        a_lhs,
                               const LevelData<FArrayBox>&  a_phi,
                               bool                         a_homogeneous,
                               VelBCHolder&                 a_velBC)
{
  CH_TIME("VelocityAMRPoissonOp::applyOpI");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  {
    CH_TIME("apply Op BC's");
    a_velBC.applyBCs(phi, dbl,
                     m_domain, m_dx,
                     a_homogeneous);
  }

  if (s_exchangeMode == 0)
    phi.exchange(phi.interval(), m_exchangeCopier);
  else if (s_exchangeMode == 1)
    phi.exchangeNoOverlap(m_exchangeCopier);
  else
    MayDay::Error("exchangeMode");

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit];
      FORT_OPERATORLAP(CHF_FRA(a_lhs[dit]),
                       CHF_CONST_FRA(phi[dit]),
                       CHF_BOX(region),
                       CHF_CONST_REAL(m_dx),
                       CHF_CONST_REAL(m_alpha),
                       CHF_CONST_REAL(m_beta));
    }
}

// ---------------------------------------------------------
void VelocityAMRPoissonOp::AMROperatorNF(LevelData<FArrayBox>&        a_LofPhi,
                                         const LevelData<FArrayBox>&  a_phi,
                                         const LevelData<FArrayBox>&  a_phiCoarse,
                                         bool                         a_homogeneousPhysBC,
                                         VelBCHolder&                 a_velBC)
{
  CH_TIME("VelocityAMRPoissonOp::AMROperatorNF");

  CH_assert(a_phi      .isDefined());
  CH_assert(a_phiCoarse.isDefined());

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  // m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  for (int comp = 0; comp < SpaceDim; comp++)
    {
      Interval intvl(comp, comp);

      LevelData<FArrayBox> phiComp;
      aliasLevelData(phiComp, &phi, intvl);

      LevelData<FArrayBox> phiCoarseComp;
      aliasLevelData(phiCoarseComp, (LevelData<FArrayBox>*) &a_phiCoarse, intvl);

      m_interpWithCoarser.coarseFineInterp(phiComp, phiCoarseComp);
    }
  //apply physical boundary conditions in applyOp
  applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC, a_velBC);
}


// ---------------------------------------------------------
void VelocityAMRPoissonOp::applyOp(LevelData<FArrayBox>&        a_lhs,
                                   const LevelData<FArrayBox>&  a_phi,
                                   const LevelData<FArrayBox>*  a_phiCoarsePtr,
                                   bool                         a_homogeneous,
                                   VelBCHolder&                 a_velBC)
{
  if (a_phiCoarsePtr == NULL)
    {
      applyOpI(a_lhs, a_phi, a_homogeneous, a_velBC);
    }
  else
    {
      AMROperatorNF(a_lhs, a_phi, *a_phiCoarsePtr, a_homogeneous, a_velBC);
    }
}
