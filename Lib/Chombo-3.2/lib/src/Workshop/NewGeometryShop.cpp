#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "GeometryShop.H"
#include "NewGeometryShop.H"
#include "ReferenceHeightIF.H"
#include "NoRefinement.H"
#include "FixedRefinement.H"
#include "ComputeCutCellMoments.H"

#include "PolyGeom.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

bool NewGeometryShop::s_verbose = false;

/**********************************************/
/*********************************************/
NewGeometryShop::NewGeometryShop(const BaseIF& a_baseIF,
                               const RealVect         & a_origin,
                               const RealVect         & a_vectDx,
                               const ProblemDomain    & a_domain,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                               const int              & a_order,
                               const int              & a_degreeP,
#else
                               const int              & a_orderPmax,
                               const int              & a_degreePmax,
#endif
                               const bool             & a_useConstraints)
  :m_phase(-1),
   m_baseIF(a_baseIF.newImplicitFunction())

{
  CH_TIME("NewGeometryShop::NewGeometryShop");

  m_domain = a_domain;
  m_origin = a_origin;
  m_vectDx = a_vectDx;

#if RECURSIVE_GEOMETRY_GENERATION == 0
  m_order          = a_order;
  m_degreeP        = a_degreeP;
#else
  m_orderPmax      = a_orderPmax;
  m_degreePmax     = a_degreePmax;
#endif
  m_useConstraints = a_useConstraints;

  //this is used for scaling face areas
  m_volScaleFactor = 1.0;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      m_volScaleFactor *= m_vectDx[idir];
    }
  //this is used for scaling boundary area
  Real maxDxComponent = 0;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if (m_vectDx[idir] > maxDxComponent)
        {
          maxDxComponent = m_vectDx[idir];
          m_bndryAreaScaleFactor = pow(m_vectDx[idir],SpaceDim -1);
        }
    }

  //records vector dx as an IndexTM
  //m_dxVect = convertRealVect2IndexTM(m_vectDx);
  for (int idir = 0 ; idir < SpaceDim ; idir++)
    {
      m_dxVect[idir] = m_vectDx[idir];
    }
  if (GLOBALDIM != SpaceDim)
    {
      const ReferenceHeightIF *refIF = dynamic_cast<const ReferenceHeightIF *>(m_baseIF);
      if (refIF == NULL)
      {
        MayDay::Error("NewGeomteryShop constructor: Couldn't cast BaseIF pointer to ReferenceHeightIF");
      }

      m_dxVect[SpaceDim] = refIF->getReferenceHeight();
    }

// Only create these in serial since they aren't distributed in parallel and
// this causes parallel runs to fail when the domain gets large enough.  This
// needs a better long-term fix.
#ifndef CH_MPI
  //initialize the FAB containing the residuals
  int nComp = 15;
  m_residuals.define(m_domain.domainBox(),nComp);
  m_residuals.setVal(0.0);

  int nGradComp = 15;
  m_gradNormal.define(m_domain.domainBox(),nGradComp);
  m_gradNormal.setVal(0.0);
#endif

  m_threshold = 1.0e-15;
  m_numCellsClipped = 0;
}

/**********************************************/
/*********************************************/
NewGeometryShop::~NewGeometryShop()
{
  CH_TIME("NewGeometryShop::~NewGeometryShop");

  delete(m_baseIF);
}

/**********************************************/
/*********************************************/



bool NewGeometryShop::isRegular(const Box           & a_region,
                               const ProblemDomain & a_domain,
                               const RealVect      & a_origin,
                               const Real          & a_dx) const
{
  CH_TIME("NewGeometryShop::isRegular");
  //set a vectDx
  RealVect vectDx = m_vectDx;

  // first check any of the Box corners is outside, and return false
  // right away. (bvs)
  RealVect physCorner;
  IntVect lo = a_region.smallEnd();
  IntVect len = a_region.size();
  Box unitBox(IntVect::Zero, IntVect::Unit);
  for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
          physCorner[idir] = vectDx[idir]*(current[idir]) + a_origin[idir];
        }
      Real functionValue = m_baseIF->value(physCorner);
      if (functionValue > 0.0 )
        {
          return false;
        }
    }

  return isRegularEveryPoint(a_region, a_domain, a_origin, a_dx);
}

bool NewGeometryShop::isRegularEveryPoint(const Box&           a_region,
                                         const ProblemDomain& a_domain,
                                         const RealVect&      a_origin,
                                         const Real&          a_dx) const
{
  CH_TIME("NewGeometryShop::isRegularEveryPoint");

  //set a vectDx
  RealVect vectDx = m_vectDx;

  // All corner indices for the current box
  Box allCorners(a_region);
  allCorners.surroundingNodes();

  RealVect physCorner;
  BoxIterator bit(allCorners);
  // If every corner is inside, the box is regular
  for (int i=0; i<2; i++)
    {
      for (; bit.ok(); ++bit, ++bit)
        {
          // Current corner
          const IntVect& corner = bit();

          // Compute physical coordinate of corner


          for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
              physCorner[idir] = vectDx[idir]*(corner[idir]) + a_origin[idir];
            }

          // If the implicit function value is positive then the current corner is
          // outside the domain
          Real functionValue = m_baseIF->value(physCorner);

          if (functionValue > 0.0 )
            {
              return false;
            }
        }
      bit.reset();
      ++bit;
    }

  return true;
}

/**********************************************/
/*********************************************/
bool NewGeometryShop::isCovered(const Box           & a_region,
                               const ProblemDomain & a_domain,
                               const RealVect      & a_origin,
                               const Real          & a_dx) const
{
  CH_TIME("NewGeometryShop::isCovered");

  //set a vectDx
  RealVect vectDx = m_vectDx;

  // first check any of the Box corners is outside, and return false
  // right away. (bvs)
  RealVect physCorner;
  IntVect lo = a_region.smallEnd();
  IntVect len = a_region.size();
  Box unitBox(IntVect::Zero, IntVect::Unit);
  for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
          physCorner[idir] = vectDx[idir]*(current[idir]) + a_origin[idir];
        }
      Real functionValue = m_baseIF->value(physCorner);
      if (functionValue < 0.0 )
        {
          return false;
        }
    }

  return isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx);
}

bool NewGeometryShop::isCoveredEveryPoint(const Box&           a_region,
                                         const ProblemDomain& a_domain,
                                         const RealVect&      a_origin,
                                         const Real&          a_dx) const
{
  CH_TIME("NewGeometryShop::isCoveredEveryPoint");

  //set a vectDx
  RealVect vectDx = m_vectDx;

  // All corner indices for the current box
  Box allCorners(a_region);
  allCorners.surroundingNodes();

  RealVect physCorner;
  BoxIterator bit(allCorners);
  // If every corner is inside, the box is regular
  for (int i=0; i<2; i++)
    {
      for (; bit.ok(); ++bit, ++bit)
        {
          // Current corner
          IntVect corner = bit();

          // Compute physical coordinate of corner

          for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
              physCorner[idir] = vectDx[idir]*(corner[idir]) + a_origin[idir];
            }

          // If the implicit function value is positive then the current corner is
          // outside the domain
          Real functionValue = m_baseIF->value(physCorner);

          if (functionValue < 0.0 )
            {
              return false;
            }
        }
      bit.reset();
      ++bit;
    }

  return true;
}

bool NewGeometryShop::isIrregular(const Box           & a_region,
                                 const ProblemDomain & a_domain,
                                 const RealVect      & a_origin,
                                 const Real          & a_dx) const
{
  CH_TIME("NewGeometryShop::isIrregular");

  //set a vectDx
  RealVect vectDx = m_vectDx;

  // first check any of the Box corners is outside, and return false
  // right away. (bvs)
  RealVect physCorner;
  IntVect lo = a_region.smallEnd();
  IntVect len = a_region.size();
  for (int idir = 0; idir < CH_SPACEDIM; ++idir)
    {
      physCorner[idir] = vectDx[idir]*(lo[idir]) + a_origin[idir];
    }
  Real originVal = m_baseIF->value(physCorner);

  Box unitBox(IntVect::Zero, IntVect::Unit);
  for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
          physCorner[idir] = vectDx[idir]*(current[idir]) + a_origin[idir];
        }
      Real functionValue = m_baseIF->value(physCorner);
      if (functionValue * originVal < 0.0 )
        {
          return true;
        }
    }

  // return isIrregularEveryPoint(a_region, a_domain, a_origin, a_dx, originVal);
  return !(isRegularEveryPoint(a_region, a_domain, a_origin, a_dx) ||
           isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx));
}

bool NewGeometryShop::isIrregularEveryPoint(const Box&           a_region,
                                         const ProblemDomain& a_domain,
                                         const RealVect&      a_origin,
                                         const Real&          a_dx,
                                         const Real&          a_originVal) const
{
  CH_TIME("NewGeometryShop::isIrregularEveryPoint");

  //set a vectDx
  RealVect vectDx = m_vectDx;

  // All corner indices for the current box
  Box allCorners(a_region);
  allCorners.surroundingNodes();

  RealVect physCorner;
  BoxIterator bit(allCorners);
  // If every corner is inside, the box is regular
  for (int i=0; i<2; i++)
    {
      for (; bit.ok(); ++bit, ++bit)
        {
          // Current corner
          IntVect corner = bit();

          // Compute physical coordinate of corner

          for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
              physCorner[idir] = vectDx[idir]*(corner[idir]) + a_origin[idir];
            }

          // If the implicit function value is positive then the current corner is
          // outside the domain
          Real functionValue = m_baseIF->value(physCorner);

          if (functionValue * a_originVal < 0.0 )
            {
              return true;
            }
        }
      bit.reset();
      ++bit;
    }

  return false;
}

/**********************************************/
/*********************************************/
void NewGeometryShop::fillGraph(BaseFab<int>&       a_regIrregCovered,
                               Vector<IrregNode>&  a_nodes,
                               const Box&          a_validRegion,
                               const Box&          a_ghostRegion,
                               const ProblemDomain & a_domain,
                               const RealVect      & a_origin,
                               const Real          & a_dx) const
{
  CH_TIME("NewGeometryShop::fillGraph");

  CH_assert(a_domain.contains(a_ghostRegion));

  IntVectSet ivsirreg = IntVectSet(DenseIntVectSet(a_ghostRegion, false));

  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
    {
      const IntVect iv =bit();
      Box miniBox(iv, iv);

      RvgDim cellCenter;
      for (int idir = 0;idir < SpaceDim; ++idir)
        {
          cellCenter[idir] = m_dxVect[idir]*(iv[idir] +0.5) + a_origin[idir];
        }
      if (GLOBALDIM != SpaceDim)
        {
          const ReferenceHeightIF *refIF = dynamic_cast<const ReferenceHeightIF *>(m_baseIF);
          if (refIF == NULL)
          {
            MayDay::Error("NewGeomteryShop constructor: Couldn't cast BaseIF pointer to ReferenceHeightIF");
          }

          cellCenter[SpaceDim] = refIF->getOrigin()[SpaceDim] + (refIF->getReferenceHeight() / 2.0);
        }

#if RECURSIVE_GEOMETRY_GENERATION != 0
      // P+R+D = m_degreePmax + m_orderPmax + GLOBALDIM and
      // R = Rmax when P = 0, so Rmax = m_degreePmax + m_orderPmax
      int maxOrder = m_degreePmax + m_orderPmax;
#endif

      //member data: sign(chosen from -1,0,1) of each vertex,
      //location of each edge intersection, cellCenter,normal and gradNormal
#if RECURSIVE_GEOMETRY_GENERATION == 0
      IFData<GLOBALDIM> edgeData(*m_baseIF,m_dxVect,cellCenter,m_order);
#else
      IFData<GLOBALDIM> edgeData(*m_baseIF,m_dxVect,cellCenter,maxOrder);
#endif

      //create a CutCellMoment object, in order to detect whether any face coincides with the interface
      CutCellMoments <GLOBALDIM> globalDimCutCell(edgeData);

#if (USING_TOP_FACE_MOMENTS)
      Iv2 bdId;
      bdId[BDID_DIR]   = 2;
      bdId[BDID_HILO]  = 1;
      CutCellMoments<GLOBALDIM-1> cutCell;
      cutCell = globalDimCutCell.m_bdCutCellMoments[bdId];
#else
      CutCellMoments<GLOBALDIM> cutCell;
      cutCell = globalDimCutCell;
#endif
      if (cutCell.isCovered())
        {
          //set covered cells to -1
          a_regIrregCovered(iv, 0) = -1;
        }
      else if (cutCell.isRegular())
        {
          //set regular cells to 1
          a_regIrregCovered(iv, 0) =  1;
        }
      else
        {
          //set irregular cells to 0
          //irregular if any face coincides with interface and edgeData.m_allVerticesIn = true
          a_regIrregCovered(iv, 0) =  0;
          if (a_validRegion.contains(iv))
            {
              ivsirreg |= iv;
            }
        }
    }

  //now loop through irregular cells and make nodes for each  one.
  for (IVSIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
    {
      VolIndex vof(ivsit(), 0);
      Real     volFrac, bndryArea;
      RealVect normal, volCentroid, bndryCentroid;
      Vector<int> loArc[SpaceDim];
      Vector<int> hiArc[SpaceDim];
      Vector<Real> loAreaFrac[SpaceDim];
      Vector<Real> hiAreaFrac[SpaceDim];
      Vector<RealVect> loFaceCentroid[SpaceDim];
      Vector<RealVect> hiFaceCentroid[SpaceDim];

      computeVoFInternals(volFrac,
                          loArc,
                          hiArc,
                          loAreaFrac,
                          hiAreaFrac,
                          bndryArea,
                          normal,
                          volCentroid,
                          bndryCentroid,
                          loFaceCentroid,
                          hiFaceCentroid,
                          ivsirreg,
                          vof,
                          a_domain,
                          a_origin,
                          a_dx,
                          m_vectDx,
                          ivsit());
      {
        CH_TIME("fillGraph::endOfirregularCellLoop");
        IrregNode newNode;
        newNode.m_cell          = ivsit();
        newNode.m_volFrac       = volFrac;
        newNode.m_cellIndex     = 0;
        newNode.m_volCentroid   = volCentroid;
        newNode.m_bndryCentroid = bndryCentroid;
        for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          {
            int loNodeInd = newNode.index(faceDir, Side::Lo);
            int hiNodeInd = newNode.index(faceDir, Side::Hi);
            newNode.m_arc[loNodeInd]          = loArc[faceDir];
            newNode.m_arc[hiNodeInd]          = hiArc[faceDir];
            newNode.m_areaFrac[loNodeInd]     = loAreaFrac[faceDir];
            newNode.m_areaFrac[hiNodeInd]     = hiAreaFrac[faceDir];
            newNode.m_faceCentroid[loNodeInd] = loFaceCentroid[faceDir];
            newNode.m_faceCentroid[hiNodeInd] = hiFaceCentroid[faceDir];
          }
        a_nodes.push_back(newNode);
      }
    } //end loop over cells in the box
}

/**********************************************/
/*********************************************/
void NewGeometryShop::computeVoFInternals(Real&               a_volFrac,
                                         Vector<int>         a_loArc[SpaceDim],
                                         Vector<int>         a_hiArc[SpaceDim],
                                         Vector<Real>        a_loAreaFrac[SpaceDim],
                                         Vector<Real>        a_hiAreaFrac[SpaceDim],
                                         Real&               a_bndryArea,
                                         RealVect&           a_normal,
                                         RealVect&           a_volCentroid,
                                         RealVect&           a_bndryCentroid,
                                         Vector<RealVect>    a_loFaceCentroid[SpaceDim],
                                         Vector<RealVect>    a_hiFaceCentroid[SpaceDim],
                                         const IntVectSet&   a_ivsIrreg,
                                         const VolIndex&     a_vof,
                                         const ProblemDomain&a_domain,
                                         const RealVect&     a_origin,
                                         const Real&         a_dx,
                                         const RealVect&     a_vectDx,
                                         const IntVect&      a_iv)const
{
  CH_TIME("GeometryShop::ComputeVofInternals");

  //assigns a_iv to m_currIv
  settCurrIv(a_iv);

  //for each CutCellMoments<dim>, we record the cell Center
  //(in physical coordinates at the global dimension)
  RvgDim cellCenter;
  for (int idir = 0;idir < SpaceDim; ++idir)
    {
      cellCenter[idir] = m_dxVect[idir]*(m_currIv[idir] +0.5) + m_origin[idir];
    }
  if (GLOBALDIM != SpaceDim)
    {
      const ReferenceHeightIF *refIF = dynamic_cast<const ReferenceHeightIF *>(m_baseIF);
      if (refIF == NULL)
      {
        MayDay::Error("NewGeomteryShop constructor: Couldn't cast BaseIF pointer to ReferenceHeightIF");
      }

      cellCenter[SpaceDim] = refIF->getOrigin()[SpaceDim] + (refIF->getReferenceHeight() / 2.0);

    }

#if RECURSIVE_GEOMETRY_GENERATION != 0
  // P+R+D = m_degreePmax + m_orderPmax + GLOBALDIM and
  // R = Rmax when P = 0, so Rmax = m_degreePmax + m_orderPmax
  int maxOrder = m_degreePmax + m_orderPmax;
#endif

  // member data: sign (chosen from -1,0,1) of each vertex,
  // location of each edge intersection, cellCenter, normal and gradNormal
#if RECURSIVE_GEOMETRY_GENERATION == 0
  IFData<GLOBALDIM> edgeData(*m_baseIF,m_dxVect,cellCenter,m_order);
#else
  IFData<GLOBALDIM> edgeData(*m_baseIF,m_dxVect,cellCenter,maxOrder);
#endif

  //construct data holders for all moments
  ComputeCutCellMoments<GLOBALDIM> computeThisVof(edgeData);

  //refines if constraints violated or normal = (0,0)
  int maxNumberRefinements = 1;
  RefinementCriterion<GLOBALDIM> refinementCriterion(maxNumberRefinements);

  //compute the moments and save answers in thisVof
#if RECURSIVE_GEOMETRY_GENERATION == 0
  computeThisVof.computeMoments(m_order,m_degreeP,m_useConstraints,refinementCriterion);
#else
  computeThisVof.computeMoments(m_orderPmax,m_degreePmax,m_useConstraints,refinementCriterion);
#endif
  CutCellMoments<GLOBALDIM> thisVof = computeThisVof.m_cutCellMoments;

#if (USING_TOP_FACE_MOMENTS)
    {
      Iv2 bdId;
      bdId[BDID_DIR]   = 2;
      bdId[BDID_HILO]  = 1;
      ((NewGeometryShop *)this)->m_cutCellMoments = thisVof.m_bdCutCellMoments[bdId];
    }
#else
    {
      ((NewGeometryShop *)this)->m_cutCellMoments = thisVof;
    }
#endif

  //fill a gradNormal check FAB
  //fillGradNormFAB();

  //fillresiduals
  //fillResiduals(degreeP);

  //fillVolFrac:
  a_volFrac = fillVolFrac();

  //fillloArc
  fillLoArc(a_loArc,
            a_ivsIrreg);

  //fillHiArc
  fillHiArc(a_hiArc,
            a_ivsIrreg);

  //face area:lo
  fillLoAreaFrac(a_loAreaFrac);

  //face area:hi
  fillHiAreaFrac(a_hiAreaFrac);

  //fill boundary area
  a_bndryArea = fillBndryArea();

  //fillNormal
  a_normal = fillNormal();

  //fillvolCentroid:
  a_volCentroid = fillvolCentroid();

  //fillbndryCentroid
  a_bndryCentroid = fillBndryCentroid();

  //fillLoFaceCentroid
  fillLoFaceCentroid(a_loFaceCentroid);

  //fillHiFaceCentroid
  fillHiFaceCentroid(a_hiFaceCentroid);

  //sanity check for computed quantities
  clipComputedVal(a_volFrac,
                  a_loAreaFrac,
                  a_hiAreaFrac,
                  a_bndryArea,
                  a_volCentroid,
                  a_bndryCentroid,
                  a_loFaceCentroid,
                  a_hiFaceCentroid,
                  m_currIv);
}

void NewGeometryShop::settCurrIv(const IntVect& a_iv)const
{
  m_currIv = a_iv;
}

#if RECURSIVE_GEOMETRY_GENERATION == 0
void NewGeometryShop::fillResiduals(int & a_degreeP)const
#else
void NewGeometryShop::fillResiduals(int & a_degreePmax)const
#endif
{
// Only create these in serial since they aren't distributed in parallel and
// this causes parallel runs to fail when the domain gets large enough.  This
// needs a better long-term fix.
#ifndef CH_MPI
#if RECURSIVE_GEOMETRY_GENERATION == 0
  for (int iDegree = 0 ; iDegree < a_degreeP + 1 ; iDegree++)
#else
  for (int iDegree = 0 ; iDegree <= a_degreePmax; iDegree++)
#endif
    {
      for (int normJ = 0 ; normJ < 3 ; normJ++)
        {
          m_residuals(m_currIv,iDegree * 3 + normJ) = m_cutCellMoments.getResidual(iDegree,normJ);
        }
    }
#endif
}

//fillVolFrac
Real NewGeometryShop::fillVolFrac()const
{
  EBorVol momentMap = VolMoment;

 //when called with VolMoment , get Vol returns the volume fraction
  Real volume = m_cutCellMoments.getVol(momentMap);

  return volume/m_volScaleFactor;
}

//records connectivity between vofs
void NewGeometryShop::fillArc(Vector<int>        a_arc[SpaceDim],
                             const int        & a_hilo,
                             const IntVectSet & a_ivsIrreg)const
{
  Iv2 bdId;
  //a_hilo is 0 or 1
  bdId[BDID_HILO] = a_hilo;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      bdId[BDID_DIR] = idir;
      bool covered =  m_cutCellMoments.getBdCutCellMoments(bdId).isCovered();

      if (covered)
        {
          a_arc[idir].resize(0);
        }
      else
        {
          a_arc[idir].resize(1);

          //otherIV is the iv in the idir direction on the a_hilo side
          IntVect otherIV = m_currIv;
          otherIV[idir] += (a_hilo*2) - 1;

          if (m_domain.contains(otherIV))
            {
              int otherCellIndex;
              if (a_ivsIrreg.contains(otherIV))
                {
                  otherCellIndex = 0;
                }
              else
                {
                  //arc to regular cell
                  otherCellIndex = -2;
                }
              a_arc[idir][0]=otherCellIndex;
            }
          else if (!m_domain.contains(otherIV))
            {
              //boundary arcs always -1
              a_arc[idir][0] = -1;
            }
        }
    }
}
  //records connectivity between vofs
  void NewGeometryShop::fillLoArc(Vector<int>        a_loArc[SpaceDim],
                                 const IntVectSet & a_ivsIrreg)const
{
   int hilo = 0;
   fillArc(a_loArc,
           hilo,
           a_ivsIrreg);
}

//records connectivity between vofs
void NewGeometryShop::fillHiArc(Vector<int>        a_hiArc[SpaceDim],
                               const IntVectSet & a_ivsIrreg)const
{
  int hilo = 1;
  fillArc(a_hiArc,
          hilo,
          a_ivsIrreg);
}

//face area
void NewGeometryShop::fillAreaFrac(Vector<Real>  a_areaFrac[SpaceDim],
                                  const int   & a_hilo)const
{
  //volume moments of coordinate aligned faces
  EBorVol momentMap = VolMoment;

  Iv2 bdId;
  bdId[BDID_HILO] = a_hilo;
  for (int idir = 0; idir < SpaceDim ;++ idir)
    {
      bdId[BDID_DIR] = idir;
      bool covered = m_cutCellMoments.getBdCutCellMoments(bdId).isCovered();
      //when called with VolMoment,get Vol returns the volume fraction
      if (!covered)
        {
          a_areaFrac[idir].resize(1);
          a_areaFrac[idir][0] = m_cutCellMoments.getBdCutCellMoments(bdId).getVol(momentMap);
          //scale area fraction
          Real scaleFactor = m_vectDx[idir]/m_volScaleFactor;
          a_areaFrac[idir][0] *= scaleFactor;
        }
      else
        {
          a_areaFrac[idir].resize(0);
        }
    }

}
//face area:lo
void NewGeometryShop::fillLoAreaFrac(Vector<Real>a_loAreaFrac[SpaceDim])const
{
  int hilo = 0;
  fillAreaFrac(a_loAreaFrac,hilo);
}
//face area:hi
void NewGeometryShop::fillHiAreaFrac(Vector<Real>a_loAreaFrac[SpaceDim])const
{
  int hilo = 1;
  fillAreaFrac(a_loAreaFrac,hilo);
}

//fill boundary area
Real NewGeometryShop::fillBndryArea()const
{
  EBorVol momentMap = EBMoment;
  Real bndryArea =  m_cutCellMoments.getVol(momentMap);
  bndryArea /= m_bndryAreaScaleFactor;

  return bndryArea;
}

//fillNormal
RealVect NewGeometryShop::fillNormal() const
{
  RealVect normal;

  IndexTM<int,SpaceDim> zeroDerivative = IndexTM<int,SpaceDim>::Zero;

  map<IndexTM<int,SpaceDim>,
      IndexTM<Real,SpaceDim>,
      LexLT<IndexTM<int,SpaceDim> > >::const_iterator iter =
        m_cutCellMoments.m_IFData.m_normalDerivatives.find(zeroDerivative);

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      normal[idir] = iter->second[idir];
    }

 return normal;
}

//fillvolCentroid
RealVect NewGeometryShop::fillvolCentroid()const
{
  EBorVol momentMap = VolMoment;

  //returns IndexTM physical coordinates
  RealVect cutCellCentroidPhysCoord = m_cutCellMoments.getCentroid(momentMap);

  //returns RealVect coordinates relative to cell center
  RealVect centroid = convert2RelativeCoord(cutCellCentroidPhysCoord);

  return centroid;
}

//fillbndryCentroid
RealVect NewGeometryShop::fillBndryCentroid()const
{
  EBorVol momentMap = EBMoment;

  //returns IndexTM physical coordinates
  RealVect cutCellCentroidPhysCoord = m_cutCellMoments.getCentroid(momentMap);

  //returns realvect relative to m_currIv
  RealVect centroid = convert2RelativeCoord(cutCellCentroidPhysCoord);
  return(centroid);
}

//fillLoFaceCentroid
void NewGeometryShop::fillFaceCentroid(Vector<RealVect>  a_faceCentroid[SpaceDim],
                                      const int       & a_hilo)const
{
//volume moments of coordinate aligned faces
  EBorVol momentMap = VolMoment;

  Iv2 bdId;
  bdId[BDID_HILO] = a_hilo;
  for (int idir = 0; idir < SpaceDim ;++ idir)
    {
      bdId[BDID_DIR] = idir;
      bool covered = m_cutCellMoments.getBdCutCellMoments(bdId).isCovered();
       if (!covered)
        {
          a_faceCentroid[idir].resize(1);

          IndexTM<Real,SpaceDim - 1>centroidRel2FacePhysCoord;
          //returns IndexTM in physical (SpaceDim-1) coordinates
          centroidRel2FacePhysCoord = m_cutCellMoments.getBdCutCellMoments(bdId).getCentroid(momentMap);
          //assigns SpaceDim-1 values to the appropriate components of a RealVect
          RealVect centroidPhysCoord;

          centroidPhysCoord[idir] = 0.0;
          for (int jdir = 0; jdir < SpaceDim; ++jdir)
            {
              if (jdir < idir)
                {
                  centroidPhysCoord[jdir] = centroidRel2FacePhysCoord[jdir];
                }
              else if (jdir > idir)
                {
                  centroidPhysCoord[jdir] = centroidRel2FacePhysCoord[jdir - 1];
                }
            }
          //returns values relative to the cell center
          RealVect centroid = convert2RelativeCoord(centroidPhysCoord);
          a_faceCentroid[idir][0] = centroid;
        }
      else
        {
          a_faceCentroid[idir].resize(0);
        }
    }
}
//fillLoFaceCentroid
void NewGeometryShop::fillLoFaceCentroid(Vector<RealVect>a_loFaceCentroid[SpaceDim])const
{
  int hilo = 0;
  fillFaceCentroid(a_loFaceCentroid,hilo);
}

  //fillHiFaceCentroid
void NewGeometryShop::fillHiFaceCentroid(Vector<RealVect>a_hiFaceCentroid[SpaceDim])const
{
 int hilo = 1;
 fillFaceCentroid(a_hiFaceCentroid,hilo);
}


//takes RealVect physical coordinates to coordinates relative to m_currIv
RealVect NewGeometryShop::convert2RelativeCoord(const RealVect& a_rVect)const
{
  //rvect[idir] = dx*(iv[idir] + 0.5 + relCoord[idir]) + origin[idir]
  RealVect retval;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      // retval[idir] = a_rVect[idir] - m_origin[idir];
      retval[idir] = a_rVect[idir] /  m_vectDx[idir];
      //retval[idir] -= (m_currIv[idir] + 0.5);
    }

  return retval;
}

//takes IndexTM physical coordinates to coordinates relative to m_currIv
RealVect NewGeometryShop::convert2RelativeCoord(const IndexTM<Real,SpaceDim>& a_rVect)const
{
  //rvect[idir] = dx*(iv[idir] + 0.5 + relCoord[idir]) + origin[idir]
  RealVect retval;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      // retval[idir] = a_rVect[idir] - m_origin[idir];
      retval[idir] = a_rVect[idir] / m_vectDx[idir];
      //      retval[idir] -= (m_currIv[idir] + 0.5);
    }
  return retval;
}

//next four functions take RealVect and IntVect to IndexTM and viseVersa
IndexTM<Real,SpaceDim>NewGeometryShop::convertRealVect2IndexTM(const RealVect& a_realVect)const
{
  IndexTM<Real,SpaceDim> retval;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      retval[idir] = a_realVect[idir];
    }
  return retval;
}
//2 of 4
RealVect NewGeometryShop::convertIndexTM2RealVect(const IndexTM<Real,SpaceDim>& a_indexTm)const
{
  RealVect retval;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      retval[idir] = a_indexTm[idir];
    }
  return retval;
}
//3 of 4
IndexTM<int,SpaceDim>NewGeometryShop::convertIntVect2IndexTM(const IntVect& a_intVect)const
{
  IndexTM<int,SpaceDim> retval;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      retval[idir] = a_intVect[idir];
    }
  return retval;
}

//4 of 4
IntVect NewGeometryShop::convertIndexTM2IntVect(const IndexTM<int,SpaceDim>& a_indexTm)const
{
  IntVect retval;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      retval[idir] = a_indexTm[idir];
    }
  return retval;
}
void NewGeometryShop::clipComputedVal(Real            &    a_volFrac,
                                     Vector<Real>         a_loAreaFrac[SpaceDim],
                                     Vector<Real>         a_hiAreaFrac[SpaceDim],
                                     Real            &    a_bndryArea,
                                     RealVect        &    a_volCentroid,
                                     RealVect        &    a_bndryCentroid,
                                     Vector<RealVect>     a_loFaceCentroid[SpaceDim],
                                     Vector<RealVect>     a_hiFaceCentroid[SpaceDim],
                                     const IntVect   &    a_iv)const

{

  //clipping
  //only report adjustments when discrepancy is above the threshold
  bool thisVofClipped = false;
  Real discrepancy = 0.0;
  Real volDiscrepancy = 0.0;
  //volFrac out of bounds
  if (a_volFrac < 0.0)
    {
      volDiscrepancy = Abs(a_volFrac);
      char message[1024];
      if (SpaceDim ==2)
        {
          sprintf(message,"vol fraction (%e) out of bounds. Clipping: (%d,%d)",
                  a_volFrac,a_iv[0],a_iv[1]);
        }
      else if (SpaceDim == 3)
        {
          sprintf(message,"vol frac (%e) out of bounds. Clipping: (%d,%d,%d)",
                  a_volFrac,a_iv[0],a_iv[1],a_iv[2]);
        }
      else
        {
          sprintf(message,"SpaceDim not 2 or 3");
        }

      if (volDiscrepancy > m_threshold)
        {
          pout() << message << endl;
        }
      //do the clipping
      thisVofClipped = true;
      a_volFrac = -a_volFrac;
    }

  if (a_volFrac > 1.0)
    {
      volDiscrepancy = Abs(1.0 - a_volFrac);
      char message[1024];
      if (SpaceDim ==2)
        {
          sprintf(message,"vol fraction (%e) out of bounds. Clipping: (%d,%d)",
                  a_volFrac,a_iv[0],a_iv[1]);
        }
      else if (SpaceDim == 3)
        {
          sprintf(message,"vol frac (%e) out of bounds. Clipping: (%d,%d,%d)",
                  a_volFrac,a_iv[0],a_iv[1],a_iv[2]);
        }
      else
        {
          sprintf(message,"SpaceDim not 2 or 3");
        }
      if (volDiscrepancy>m_threshold)
        {
          pout() << message << endl;
        }
      //do the clipping
      thisVofClipped = true;
      a_volFrac = 1.0;
    }

  //area frac out of bounds
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_loAreaFrac[idir].size();num ++)
        {
          //lo frac too high
          if (a_loAreaFrac[idir][num] > 1.0)
            {
              discrepancy = Abs(1 - a_loAreaFrac[idir][num]);
              char message[1024];
              if (SpaceDim ==2)
                {
                  sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d)",
                          a_loAreaFrac[idir][num],a_iv[0],a_iv[1]);
                }
              else if (SpaceDim == 3)
                {
                  sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                          a_loAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                }
              else
                {
                  sprintf(message,"SpaceDim not 2 or 3");
                }
              if (discrepancy>m_threshold)
                {
                  pout() << message << endl;
                }

              //do the clipping
              thisVofClipped = true;
              a_loAreaFrac[idir][num] = 1.0;
            }
          //lo frac too low
          if (a_loAreaFrac[idir][num] < 0.0)
            {
              discrepancy = Abs(a_loAreaFrac[idir][num]);
              char message[1024];
              if (SpaceDim ==2)
                {
                  sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d)",
                          a_loAreaFrac[idir][num],a_iv[0],a_iv[1]);
                }
              else if (SpaceDim == 3)
                {
                  sprintf(message,"lo area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                          a_loAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                }
              else
                {
                  sprintf(message,"SpaceDim not 2 or 3");
                }
              if (discrepancy>m_threshold)
                {
                  pout() << message << endl;
                }

              //do the clipping
              thisVofClipped = true;
              a_loAreaFrac[idir][num] = -a_loAreaFrac[idir][num];
            }

        }
    }

  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_hiAreaFrac[idir].size();num ++)
        {
          //hi frac too high
          if (a_hiAreaFrac[idir][num] > 1.0)
            {
              discrepancy = Abs(1 - a_hiAreaFrac[idir][num]);
              char message[1024];
              if (SpaceDim ==2)
                {
                  sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d)",
                          a_hiAreaFrac[idir][num],a_iv[0],a_iv[1]);
                }
              else if (SpaceDim == 3)
                {
                  sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                          a_hiAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                }
              else
                {
                  sprintf(message,"SpaceDim not 2 or 3");
                }
              if (discrepancy>m_threshold)
                {
                  pout() << message << endl;
                }

              //do the clipping
              thisVofClipped = true;
              a_hiAreaFrac[idir][num] = 1.0;
            }
          //hi frac too low
          if (a_hiAreaFrac[idir][num] < 0.0)
            {
              discrepancy = Abs(a_hiAreaFrac[idir][num]);
              char message[1024];
              if (SpaceDim ==2)
                {
                  sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d)",
                          a_hiAreaFrac[idir][num],a_iv[0],a_iv[1]);
                }
              else if (SpaceDim == 3)
                {
                  sprintf(message,"hi area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                          a_hiAreaFrac[idir][num],a_iv[0],a_iv[1],a_iv[2]);
                }
              else
                {
                  sprintf(message,"SpaceDim not 2 or 3");
                }
              if (discrepancy>m_threshold)
                {
                  pout() << message << endl;
                }

              //do the clipping
              thisVofClipped = true;
              a_hiAreaFrac[idir][num] = -a_hiAreaFrac[idir][num];
            }

        }

    }

  //bndry area out of bounds
  if (a_bndryArea < 0.0)
    {
      discrepancy = Abs(a_bndryArea);
      char message[1024];
      if (SpaceDim ==2)
        {
          sprintf(message,"boundary area fraction (%e) out of bounds. Clipping: (%d,%d)",
                  a_bndryArea,a_iv[0],a_iv[1]);
        }
      else if (SpaceDim == 3)
        {
          sprintf(message,"boundary area fraction (%e) out of bounds. Clipping: (%d,%d,%d)",
                  a_bndryArea,a_iv[0],a_iv[1],a_iv[2]);
        }
      else
        {
          sprintf(message,"SpaceDim not 2 or 3");
        }
      if (discrepancy>m_threshold)
        {
          pout() << message << endl;
        }

      //do the clipping
      thisVofClipped = true;
      a_bndryArea = 0.0;
    }

  //volCentroid out of bounds
  for (int idir =0;idir<SpaceDim;++idir)
    {
      if (a_volCentroid[idir] > 0.5)
        {
          discrepancy = Abs(0.5 - a_volCentroid[idir]);
          char message[1024];
          if (SpaceDim ==2)
            {
              sprintf(message,"volCentroid (%e) out of bounds. Clipping: (%d,%d)",
                      a_volCentroid[idir],a_iv[0],a_iv[1]);
            }
          else if (SpaceDim == 3)
            {
              sprintf(message,"volCentroid(%e) out of bounds. Clipping: (%d,%d,%d)",
                      a_volCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
            }
          else
            {
              sprintf(message,"SpaceDim not 2 or 3");
            }
          if (discrepancy>m_threshold)
            {
              pout() << message << endl;
            }
          //do the clipping
          thisVofClipped = true;
          a_volCentroid[idir] = 0.5;
        }
      if (a_volCentroid[idir] < -0.5)
        {
          discrepancy = Abs(-0.5 - a_volCentroid[idir]);
          char message[1024];
          if (SpaceDim ==2)
            {
              sprintf(message,"volCentroid (%e) out of bounds. Clipping: (%d,%d)",
                      a_volCentroid[idir],a_iv[0],a_iv[1]);
            }
          else if (SpaceDim == 3)
            {
              sprintf(message,"volCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                      a_volCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
            }
          else
            {
              sprintf(message,"SpaceDim not 2 or 3");
            }
          if (discrepancy > m_threshold)
            {
              pout() << message << endl;
            }
          //do the clipping
          thisVofClipped = true;
          a_volCentroid[idir] = -0.5;
        }

      //boundary centroid out of bounds
      if (a_bndryCentroid[idir] > 0.5)
        {
          discrepancy = Abs(0.5 - a_bndryCentroid[idir]);
          char message[1024];
          if (SpaceDim ==2)
            {
              sprintf(message,"bndryCentroid (%e) out of bounds. Clipping: (%d,%d)",
                      a_bndryCentroid[idir],a_iv[0],a_iv[1]);
            }
          else if (SpaceDim == 3)
            {
              sprintf(message,"boundary Centroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                      a_bndryCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
            }
          else
            {
              sprintf(message,"SpaceDim not 2 or 3");
            }
          if (discrepancy>m_threshold)
            {
              pout() << message << endl;
            }
          //do the clipping
          thisVofClipped = true;
          a_bndryCentroid[idir] = 0.5;
        }

      if (a_bndryCentroid[idir] < -0.5)
        {
          discrepancy = Abs(-0.5 - a_bndryCentroid[idir]);
          char message[1024];
          if (SpaceDim ==2)
            {
              sprintf(message,"bndryCentroid (%e) out of bounds. Clipping: (%d,%d)",
                      a_bndryCentroid[idir],a_iv[0],a_iv[1]);
            }
          else if (SpaceDim == 3)
            {
              sprintf(message,"bndryCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                      a_bndryCentroid[idir],a_iv[0],a_iv[1],a_iv[2]);
            }
          else
            {
              sprintf(message,"SpaceDim not 2 or 3");
            }
          if (discrepancy > m_threshold)
            {
              pout() << message << endl;
            }
          //do the clipping
          thisVofClipped = true;
          a_bndryCentroid[idir] = -0.5;
        }
    }
  //loFaceCentroid out of bounds
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_loFaceCentroid[idir].size();num ++)
        {
          for (int jdir = 0; jdir<SpaceDim; ++jdir)
            {
              if (jdir != idir)
                {
                  if (a_loFaceCentroid[idir][num][jdir] > 0.5)
                    {
                      discrepancy = Abs(0.5 - a_loFaceCentroid[idir][num][jdir]);
                      char message[1024];
                      if (SpaceDim ==2)
                        {
                          sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                  a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                        }
                      else if (SpaceDim == 3)
                        {
                          sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                  a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                        }
                      else
                        {
                          sprintf(message,"SpaceDim not 2 or 3");
                        }
                      if (discrepancy > m_threshold)
                        {
                          pout() << message << endl;
                        }
                      //do the clipping
                      thisVofClipped = true;
                      a_loFaceCentroid[idir][num][jdir] = 0.5;
                    }
                  if (a_loFaceCentroid[idir][num][jdir] < -0.5)
                    {
                      discrepancy = Abs(-0.5 - a_loFaceCentroid[idir][num][jdir]);
                      char message[1024];
                      if (SpaceDim ==2)
                        {
                          sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                  a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                        }
                      else if (SpaceDim == 3)
                        {
                          sprintf(message,"loFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                  a_loFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                        }
                      else
                        {
                          sprintf(message,"SpaceDim not 2 or 3");
                        }
                      if (discrepancy > m_threshold)
                        {
                          pout() << message << endl;
                        }
                      //do the clipping
                      thisVofClipped = true;
                      a_loFaceCentroid[idir][num][jdir] = -0.5;
                    }
                }
            }
        }
    }
  //hiFaceCentroid out of bounds
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_hiFaceCentroid[idir].size();num ++)
        {
          for (int jdir = 0; jdir<SpaceDim; ++jdir)
            {
              if (jdir != idir)
                {
                  if (a_hiFaceCentroid[idir][num][jdir] > 0.5)
                    {
                      discrepancy = Abs(0.5 - a_hiFaceCentroid[idir][num][jdir]);
                      char message[1024];
                      if (SpaceDim ==2)
                        {
                          sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                  a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                        }
                      else if (SpaceDim == 3)
                        {
                          sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                  a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                        }
                      else
                        {
                          sprintf(message,"SpaceDim not 2 or 3");
                        }
                      if (discrepancy > m_threshold)
                        {
                          pout() << message << endl;
                        }
                      //do the clipping
                      thisVofClipped = true;
                      a_hiFaceCentroid[idir][num][jdir] = 0.5;
                    }
                  if (a_hiFaceCentroid[idir][num][jdir] < -0.5)
                    {
                      discrepancy = Abs(0.5 - a_hiFaceCentroid[idir][num][jdir]);
                      char message[1024];
                      if (SpaceDim ==2)
                        {
                          sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d)",
                                  a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1]);
                        }
                      else if (SpaceDim == 3)
                        {
                          sprintf(message,"hiFaceCentroid (%e) out of bounds. Clipping: (%d,%d,%d)",
                                  a_hiFaceCentroid[idir][num][jdir],a_iv[0],a_iv[1],a_iv[2]);
                        }
                      else
                        {
                          sprintf(message,"SpaceDim not 2 or 3");
                        }
                      if (discrepancy > m_threshold)
                        {
                          pout() << message << endl;
                        }
                      //do the clipping
                      thisVofClipped = true;
                      a_hiFaceCentroid[idir][num][jdir] = -0.5;
                    }
                }

              if (thisVofClipped)
                {
                  NewGeometryShop* changedThis = (NewGeometryShop *) this;
                  changedThis->m_numCellsClipped += 1;
                }
            }
        }
    }
}
int NewGeometryShop::getNumCellsClipped()
{
  return m_numCellsClipped;
}

#if RECURSIVE_GEOMETRY_GENERATION == 0
void NewGeometryShop::outputResidual(int & type,int & a_degreeP)const
#else
void NewGeometryShop::outputResidual(int & type,int & a_degreePmax)const
#endif
{
// Only create these in serial since they aren't distributed in parallel and
// this causes parallel runs to fail when the domain gets large enough.  This
// needs a better long-term fix.
#ifndef CH_MPI
  if (type == 0)
    {
#ifdef CH_USE_HDF5
      const char name[] = "Residual.hdf5";
      int nComp = 15;
      Vector<string> labels(nComp);
#if RECURSIVE_GEOMETRY_GENERATION == 0
      for (int iDegree = 0 ; iDegree <= a_degreeP ; iDegree++)
#else
      for (int iDegree = 0 ; iDegree <= a_degreePmax; iDegree++)
#endif
        {
          for (int iRes = 0 ; iRes < 3 ; iRes ++)
            {
              char labelChSt[80];
              sprintf(labelChSt, "ResidualsDegree_%d Norm_%d", iDegree,iRes);
              string label(labelChSt);
              labels[iDegree*3 + iRes] = label;
            }
        }
      writeFABname(&m_residuals,name,labels);
#endif
    }
  else
    {
#if RECURSIVE_GEOMETRY_GENERATION == 0
      for (int iDegree = 0 ; iDegree <=a_degreeP; iDegree++)
#else
      for (int iDegree = 0 ; iDegree <= a_degreePmax; iDegree++)
#endif
        {
          for (int iRes = 0 ; iRes < 3; iRes++)
            {
              Real res = m_residuals.norm(iRes,iDegree*3+iRes);
              pout()<<"Residual["<<iDegree*3+iRes<<"]="<<res<<endl;
            }
        }
    }
#endif
}

void NewGeometryShop::outputGradNormal()const
{
// Only create these in serial since they aren't distributed in parallel and
// this causes parallel runs to fail when the domain gets large enough.  This
// needs a better long-term fix.
#ifndef CH_MPI
#ifdef CH_USE_HDF5
  const char name[] = "GradNormal.hdf5";
  writeFABname(&m_gradNormal,name);
#endif
#endif
}

#include "NamespaceFooter.H"
