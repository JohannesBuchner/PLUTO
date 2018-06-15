#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NoFlowVortex.H"
#include "NoFlowAdvectBC.H"
#include "ExtrapAdvectBC.H"
#include "VoFIterator.H"
#include "NeumannPoissonEBBC.H"

#include "NamespaceHeader.H"

RefCountedPtr<BaseEBBCFactory>
NoFlowVortex::
getPressureEBBC() const
{
  NeumannPoissonEBBCFactory* neumptr = new NeumannPoissonEBBCFactory();

  neumptr->setValue(0.);
  return RefCountedPtr<BaseEBBCFactory>(neumptr );
}

///
RefCountedPtr<BaseDomainBCFactory>
NoFlowVortex::getVelBC(int a_icomp) const
{
  //this is only used in the helmholtz solve so always no slip
  DirichletPoissonDomainBCFactory* dirptr = new DirichletPoissonDomainBCFactory();
  Real val = 0.0;
  dirptr->setValue(val);
  RefCountedPtr<BaseDomainBCFactory>  retval( dirptr );
  return retval;
}

///
RefCountedPtr<BaseEBBCFactory>
NoFlowVortex::getVelocityEBBC(int a_velComp) const
{
  DirichletPoissonEBBCFactory* dirptr = new DirichletPoissonEBBCFactory();
  Real val = 0.0;
  dirptr->setValue(val);
  int orderEBBC = 2;
  dirptr->setOrder(orderEBBC);
  RefCountedPtr<BaseEBBCFactory>  retval( dirptr );
  return retval;
}
///
NoFlowVortex::NoFlowVortex(const RealVect& a_center,
                           const Real&     a_coreRadius,
                           const Real&     a_strength)
{
  m_center = a_center;
  m_coreRadius = a_coreRadius;
  m_strength = a_strength;
  Real domVal = 0.0;

  m_domBCPress  = RefCountedPtr<BaseDomainBCFactory>(
    new NeumannPoissonDomainBCFactory());
  m_domBCMACVel = RefCountedPtr<BaseDomainBCFactory>(
    new DirichletPoissonDomainBCFactory());

  NeumannPoissonDomainBCFactory*   bcPress = (NeumannPoissonDomainBCFactory*)   (&(*m_domBCPress));
  DirichletPoissonDomainBCFactory* bcMACVel =(DirichletPoissonDomainBCFactory*) (&(*m_domBCMACVel));
  bcPress->setValue(domVal);
  bcMACVel->setValue(domVal);
}
///
/**
   Return pressure boundary conditions for domain.
*/
RefCountedPtr<BaseDomainBCFactory>
NoFlowVortex::
getPressBC() const
{
  return m_domBCPress;
}

RefCountedPtr<EBPhysIBCFactory>
NoFlowVortex::
getVelAdvectBC(int a_velComp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval( new NoFlowAdvectBCFactory(a_velComp));
  return retval;
}
RefCountedPtr<EBPhysIBCFactory>
NoFlowVortex::
getScalarAdvectBC(const int&  a_comp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval( new ExtrapAdvectBCFactory() );
  return retval;
}

///
RefCountedPtr<BaseDomainBCFactory>
NoFlowVortex::
getMACVelBC() const
{
  return m_domBCMACVel;
}

///
void
NoFlowVortex::getRadius(Real&           a_radius,
                        RealVect&       a_xdiff,
                        const RealVect& a_xval) const
{
  a_radius = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_xdiff[idir] = a_xval[idir] - m_center[idir];
      a_radius += a_xdiff[idir]*a_xdiff[idir];
    }

  a_radius = sqrt(a_radius);
}
///
void
NoFlowVortex::
getVelPt(RealVect& a_vel, const RealVect& a_xval, Real a_pi) const
{
  Real radius;
  RealVect xdiff;
  getRadius(radius, xdiff, a_xval);

  radius /= m_coreRadius;
  xdiff  /= m_coreRadius;

  Real vmag = cos(a_pi*radius/2.0)*sin(a_pi*radius/2.0);

#if CH_SPACEDIM == 2

  a_vel[0] = -vmag*xdiff[1]/radius;
  a_vel[1] =  vmag*xdiff[0]/radius;

#elif CH_SPACEDIM == 3

  Real xplusy = sqrt(xdiff[0]*xdiff[0] + xdiff[1]*xdiff[1]);
  Real xplusz = sqrt(xdiff[0]*xdiff[0] + xdiff[2]*xdiff[2]);
  Real yplusz = sqrt(xdiff[1]*xdiff[1] + xdiff[2]*xdiff[2]);

  a_vel[0] = -vmag*yplusz/radius;
  a_vel[1] =  vmag*xplusz/radius;
  a_vel[2] =  vmag*xplusy/radius;

#else
  bogus_spacedim();
#endif

}

void
NoFlowVortex::
getScalarPt(Real& a_scal, const RealVect& a_xval) const
{
  RealVect xdiff;
  Real radius;
  getRadius(radius, xdiff, a_xval);
  if (radius < m_coreRadius)
    {
      a_scal = radius;
    }
  else
    {
      a_scal = 0.0;
    }
}
void
NoFlowVortex::
getXVal(RealVect& a_xval,
        const RealVect& a_origin,
        const VolIndex& a_vof,
        const RealVect& a_dx) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_xval[idir] = a_origin[idir] + a_dx[idir]*(Real(a_vof.gridIndex()[idir]) + 0.5);
    }
}
///
void
NoFlowVortex::
initializeVelocity(LevelData<EBCellFAB>&    a_velocity,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
  Real pi = 4.0*atan(1.0);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_grids.get(dit()));
      for (VoFIterator vofit(ivsBox, a_ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          RealVect xval, velpt;
          getXVal(xval, a_origin,vofit(), a_dx);
          getVelPt(velpt, xval, pi);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_velocity[dit()](vofit(), idir) = velpt[idir];
            }
        }
    }
}
///
void
NoFlowVortex::
initializePressure(LevelData<EBCellFAB>&    a_pressure,
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
          RealVect xval, pt;
          getXVal(xval, a_origin,vofit(), a_dx);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_pressure[dit()](vofit(), idir) = 1.e99;
            }
        }
    }
}
///
void
NoFlowVortex::
initializeScalar ( LevelData<EBCellFAB>&    a_scalar,
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
          RealVect xval;
          Real scal;
          getXVal(xval, a_origin,vofit(), a_dx);
          getScalarPt(scal, xval);
          a_scalar[dit()](vofit(), 0) = scal;
        }
    }
}

#include "NamespaceFooter.H"
