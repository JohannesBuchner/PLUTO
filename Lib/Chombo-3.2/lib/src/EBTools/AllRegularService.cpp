#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#include "AllRegularService.H"
#include "PolyGeom.H"
#include "NamespaceHeader.H"
/*******************/
/*******************/
AllRegularService::AllRegularService()
{
  PolyGeom::setVectDx(RealVect::Unit);
}

/*******************/
/*******************/
AllRegularService::~AllRegularService()
{
}

/*******************/
/*******************/
bool
AllRegularService::isRegular(const Box& a_region,
                             const ProblemDomain& a_domain,
                             const RealVect& a_origin,
                             const Real& a_dx) const
{
  return true;
}

/*******************/
/*******************/
bool
AllRegularService::isCovered(const Box& a_region,
                             const ProblemDomain& a_domain,
                             const RealVect& a_origin,
                             const Real& a_dx) const
{
  return false;
}

/*******************/
/*******************/
void
AllRegularService::fillGraph(BaseFab<int>&        a_regIrregCovered,
                             Vector<IrregNode>&   a_nodes,
                             const Box&           a_validRegion,
                             const Box&           a_ghostRegion,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx) const
{
  PolyGeom::setVectDx(RealVect::Unit);
  //set all cells to regular
  a_regIrregCovered.setVal(1);
}
/*******************/
/*******************/
#include "NamespaceFooter.H"
