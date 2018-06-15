#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PoisselleTube.H"
#include "PoisselleTubeVelBC.H"
#include "PoisselleTubePoissonDomainBC.H"
#include "VoFIterator.H"
#include "NeumannPoissonEBBC.H"
#include "EBLevelDataOps.H"
#include "ParmParse.H"

#include "UsingNamespace.H"

/***/
RefCountedPtr<BaseEBBCFactory>
PoisselleTube::getPressureEBBC()  const
{
  NeumannPoissonEBBCFactory* neumptr = new NeumannPoissonEBBCFactory();

  neumptr->setValue(0.);

  RefCountedPtr<BaseEBBCFactory>  retval = RefCountedPtr<BaseEBBCFactory>(neumptr);
  return retval;
}
/***/
RefCountedPtr<BaseDomainBCFactory>
PoisselleTube::getVelBC(int a_icomp) const
{
  //this is only used in the helmholtz solve so always no slip
  RefCountedPtr<BaseDomainBCFactory> helmBC = RefCountedPtr<BaseDomainBCFactory>(new PoisselleTubeHelmholtzDomainBCFactory(m_bcval[a_icomp]));
  return helmBC;
}

///
RefCountedPtr<BaseEBBCFactory>
PoisselleTube::getVelocityEBBC(int a_velComp) const
{
  DirichletPoissonEBBCFactory* dirptr = new DirichletPoissonEBBCFactory();
  Real val = 0.0;
  dirptr->setValue(val);
  int orderEBBC = 1;
  ParmParse pp;
  pp.get("order_ebbc", orderEBBC);
  dirptr->setOrder(orderEBBC);
  RefCountedPtr<BaseEBBCFactory>  retval = RefCountedPtr<BaseEBBCFactory>(dirptr);
  return retval;
}
///
RefCountedPtr<BaseDomainBCFactory>
PoisselleTube::
getPressBC() const
{
  RefCountedPtr<BaseDomainBCFactory> poisBC = RefCountedPtr<BaseDomainBCFactory>(new PoisselleTubePoissonDomainBCFactory(m_bcval));
  return poisBC;
}
///
RefCountedPtr<EBPhysIBCFactory>
PoisselleTube::
getVelAdvectBC(int a_velComp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval = RefCountedPtr<EBPhysIBCFactory>(new PoisselleTubeVelBCFactory(m_bcval[a_velComp], a_velComp));
  return retval;
}
///
RefCountedPtr<EBPhysIBCFactory>
PoisselleTube::
getScalarAdvectBC(const int&  a_comp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval = RefCountedPtr<EBPhysIBCFactory>(new ExtrapAdvectBCFactory());
  return retval;
}

///
RefCountedPtr<BaseDomainBCFactory>
PoisselleTube::
getMACVelBC() const
{
  RefCountedPtr<BaseDomainBCFactory> poisBC = RefCountedPtr<BaseDomainBCFactory>(new PoisselleTubePoissonDomainBCFactory(m_bcval));
  return poisBC;
}

///
void
PoisselleTube::
initializeVelocity(LevelData<EBCellFAB>&    a_velocity,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_grids.get(dit()));
      for (VoFIterator vofit(ivsBox, a_ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          RealVect xval = EBArith::getVofLocation(vofit() , a_dx, RealVect::Zero);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              //bogus values for normal and time because they do not matter
              Real velComp = m_bcval[idir].value(xval, RealVect::Zero, a_time, idir);
              a_velocity[dit()](vofit(), idir) = velComp;
            }
        }
    }
}
///
void
PoisselleTube::
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
PoisselleTube::
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

