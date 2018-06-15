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

#include "AMRNavierStokesFactory.H"

// ---------------------------------------------------------------
AMRNavierStokesFactory::AMRNavierStokesFactory()
{
  setDefaultValues();
}

// ---------------------------------------------------------------
AMRNavierStokesFactory::~AMRNavierStokesFactory()
{
  if (m_physBCPtr != NULL)
    {
      delete m_physBCPtr;
      m_physBCPtr = NULL;
    }
}

// ---------------------------------------------------------------
void
AMRNavierStokesFactory::CFL(const Real a_cfl)
{
  m_cfl = a_cfl;
}

// ---------------------------------------------------------------
void
AMRNavierStokesFactory::refinementThreshold(const Real a_refine_threshold)
{
  m_refine_thresh = a_refine_threshold;
}

// ---------------------------------------------------------------
void
AMRNavierStokesFactory::setPhysBC(const PhysBCUtil& a_physBC)
{
  if (m_physBCPtr != NULL)
    {
      delete m_physBCPtr;
      m_physBCPtr = NULL;
    }
  m_physBCPtr = a_physBC.newPhysBCUtil();
}

// ---------------------------------------------------------------
void
AMRNavierStokesFactory::setLimitSolverCoarsening(bool a_limitSolverCoarsening)
{
  m_limitSolverCoarsening = a_limitSolverCoarsening;
}

// ---------------------------------------------------------------
// "virtual constructor"
AMRLevel*
AMRNavierStokesFactory::new_amrlevel() const
{
  AMRNavierStokes* amrns_ptr = new AMRNavierStokes();
  amrns_ptr->CFL(m_cfl);
  amrns_ptr->refinementThreshold(m_refine_thresh);
  CH_assert (m_physBCPtr != NULL);
  amrns_ptr->setPhysBC(*m_physBCPtr);
  amrns_ptr->limitSolverCoarsening(m_limitSolverCoarsening);

  return (static_cast <AMRLevel*> (amrns_ptr));
}

// ---------------------------------------------------------------
void
AMRNavierStokesFactory::setDefaultValues()
{
  m_cfl = 0.5;
  // don't really know what default value to put here
  m_refine_thresh = -1.0;
  // set default Ptr?
  m_physBCPtr = new PhysBCUtil;
  // default is to not limit coarsening
  m_limitSolverCoarsening = false;

}
