#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "InflowOutflowIBC.H"
#include "InflowOutflowAdvectBC.H"
#include "InflowOutflowPoissonDomainBC.H"
#include "VoFIterator.H"
#include "NeumannPoissonEBBC.H"
#include "InflowOutflowIBC.H"
#include "EBLevelDataOps.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

/***/
RefCountedPtr<BaseEBBCFactory>
InflowOutflowIBC::getPressureEBBC()  const
{
  NeumannPoissonEBBCFactory* neumptr = new NeumannPoissonEBBCFactory();

  neumptr->setValue(0.);

  RefCountedPtr<BaseEBBCFactory>  retval = RefCountedPtr<BaseEBBCFactory>(neumptr);
  return retval;
}
/***/
RefCountedPtr<BaseDomainBCFactory>
InflowOutflowIBC::getVelBC(int a_icomp) const
{
  //this is only used in the helmholtz solve so always no slip
  RefCountedPtr<BaseDomainBCFactory> helmBC = RefCountedPtr<BaseDomainBCFactory>(new InflowOutflowHelmholtzDomainBCFactory(m_flowDir,
                                                                                                                           m_inflowVel,
                                                                                                                           m_doPoiseInflow,
                                                                                                                           m_doSlipWallsHi,
                                                                                                                           m_doSlipWallsLo,
                                                                                                                           m_poiseInflowFunc,
                                                                                                                           m_doWomersleyInflow));
  return helmBC;
}

///
RefCountedPtr<BaseEBBCFactory>
InflowOutflowIBC::getVelocityEBBC(int a_velComp) const
{
  RefCountedPtr<BaseEBBCFactory>  retval;
  //slip walls only apply to the domain boundary
  DirichletPoissonEBBCFactory* dirptr = new DirichletPoissonEBBCFactory();
  Real val = 0.0;
  dirptr->setValue(val);
  dirptr->setOrder(m_orderEBBC);
  retval = RefCountedPtr<BaseEBBCFactory>(dirptr);

  return retval;
}
///
RefCountedPtr<BaseDomainBCFactory>
InflowOutflowIBC::
getPressBC() const
{
  RefCountedPtr<BaseDomainBCFactory> poisBC
    = RefCountedPtr<BaseDomainBCFactory>
        (new InflowOutflowPoissonDomainBCFactory(m_flowDir,
                                                 m_inflowVel,
                                                 m_doPoiseInflow,
                                                 m_doSlipWallsHi,
                                                 m_doSlipWallsLo,
                                                 m_poiseInflowFunc,
                                                 m_doWomersleyInflow));
  return poisBC;
}
///
RefCountedPtr<EBPhysIBCFactory>
InflowOutflowIBC::
getVelAdvectBC(int a_velComp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval
    = RefCountedPtr<EBPhysIBCFactory>
        (new InflowOutflowAdvectBCFactory(m_flowDir,
                                          m_inflowVel,
                                          a_velComp,
                                          m_doPoiseInflow,
                                          m_doSlipWallsHi,
                                          m_doSlipWallsLo,
                                          m_poiseInflowFunc,
                                          m_doWomersleyInflow));
  return retval;
}
///
RefCountedPtr<EBPhysIBCFactory>
InflowOutflowIBC::
getScalarAdvectBC(const int&  a_comp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval = RefCountedPtr<EBPhysIBCFactory>(new ExtrapAdvectBCFactory());
  return retval;
}

///
RefCountedPtr<BaseDomainBCFactory>
InflowOutflowIBC::
getMACVelBC() const
{
  RefCountedPtr<BaseDomainBCFactory> poisBC
    = RefCountedPtr<BaseDomainBCFactory>
        (new InflowOutflowPoissonDomainBCFactory(m_flowDir,
                                                 m_inflowVel,
                                                 m_doPoiseInflow,
                                                 m_doSlipWallsHi,
                                                 m_doSlipWallsLo,
                                                 m_poiseInflowFunc,
                                                 m_doWomersleyInflow));
  return poisBC;
}

///
void
InflowOutflowIBC::
initializeVelocity(LevelData<EBCellFAB>&    a_velocity,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
  // EBLevelDataOps::setVal(a_velocity, m_inflowVel, m_flowDir);
  if (m_initPoiseData)
    {
      // if (!m_isPoiseDefined)
      //   {
      //     poiseuilleDefine(a_domain, a_dx);
      //   }
      for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB& velFAB = a_velocity[dit()];
          const Box& grid   = a_grids.get(dit());
          const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
          IntVectSet ivsTot(grid);
          for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              //set coordinates based on slip walls
              RealVect prob_lo = RealVect::Zero;
              //assumes inflow boundary is always on lo side of domain
              const RealVect loc  = EBArith::getVofLocation(vof, a_dx, prob_lo);
              Real radius = m_poiseInflowFunc->getRadius(loc);
              velFAB(vof, m_flowDir) = m_poiseInflowFunc->getVel(radius)[m_flowDir];
            }
        }
    }
}
///
void
InflowOutflowIBC::
initializePressureGradient(LevelData<EBCellFAB>&    a_gradient,
                           const DisjointBoxLayout& a_grids,
                           const EBISLayout&        a_ebisl,
                           const ProblemDomain&     a_domain,
                           const RealVect&          a_origin,
                           const Real&              a_time,
                           const RealVect&          a_dx) const
{
  EBLevelDataOps::setToZero(a_gradient);
  if (m_initPoiseData)
    {
      // if (!m_isPoiseDefined)
      //   {
      //     poiseuilleDefine(a_domain, a_dx);
      //   }
      for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB& gradFAB = a_gradient[dit()];
          const Box& grid   = a_grids.get(dit());
          const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
          IntVectSet ivsTot(grid);
          for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              gradFAB(vof, m_flowDir) = m_poiseInflowFunc->getGradP()[m_flowDir];
            }
        }
    }
}
///
void
InflowOutflowIBC::
initializePressure(LevelData<EBCellFAB>&    a_pressure,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
  EBLevelDataOps::setToZero(a_pressure);
}
///
void
InflowOutflowIBC::
initializeScalar ( LevelData<EBCellFAB>&    a_scalar,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
  EBLevelDataOps::setToZero(a_scalar);
}

#include "NamespaceFooter.H"
