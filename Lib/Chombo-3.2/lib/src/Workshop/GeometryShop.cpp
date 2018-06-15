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

#include "GeometryService.H"
#include "GeometryShop.H"

#include "PolyGeom.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

bool GeometryShop::s_verbose = false;

GeometryShop::GeometryShop(const BaseIF& a_localGeom,
                           int           a_verbosity,
                           RealVect      a_vectDx,
                           Real          a_thrshdVoF)
  :m_phase(-1)
{
  CH_TIME("GeometryShop::GeometryShop");

  m_vectDx = a_vectDx;

  RealVect vectDx;
  vectDx = RealVect::Unit;
  PolyGeom::setVectDx(vectDx);

  m_implicitFunction = a_localGeom.newImplicitFunction();

  // See if this is an STL description - m_stlIF will be NULL if it isn't
  m_stlIF = dynamic_cast<const STLIF *>(m_implicitFunction);

  m_verbosity = a_verbosity;

  Real arg1 = 10.0;
  Real arg2 = -m_verbosity;

  m_threshold = 1.0e-15*pow(arg1, arg2);

  m_numCellsClipped = 0;

  m_thrshdVoF = a_thrshdVoF;

  m_STLBoxSet = false;
}

GeometryShop::~GeometryShop()
{
  delete(m_implicitFunction);
}

bool GeometryShop::twoEdgeIntersections(edgeMo a_edges[4])const
{
  bool retval;
  int count = 0;
  for (int iedge = 0;iedge<4;++iedge)
    {
      if (!a_edges[iedge].isCovered() && a_edges[iedge].getEdgeLength()<1.0)
        {
          count += 1;
        }
    }
  if (count >= 2 && count !=0)
    {
      retval = true;
    }
  else
    {
      retval = false;
    }
  return retval;
}

bool GeometryShop::isRegular(const Box&           a_region,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx) const
{
  CH_TIME("GeometryShop::isRegular");
  // set a vectDx
  RealVect vectDx = m_vectDx;

  // If this isn't an STL description then the value function works
  if (m_stlIF == NULL)
  {
    // first check any of the Box corners is outside, and return false
    // right away. (bvs)
    IntVect lo = a_region.smallEnd();
    IntVect len = a_region.size();
    Box unitBox(IntVect::Zero, IntVect::Unit);
    for (BoxIterator bit(unitBox); bit.ok(); ++bit)
      {
        IntVect current = lo + len*bit();
        RealVect physCorner;
        for (int idir = 0; idir < CH_SPACEDIM; ++idir)
          {
            physCorner[idir] = vectDx[idir]*current[idir] + a_origin[idir];
          }
        Real functionValue = m_implicitFunction->value(physCorner);
        if (functionValue > 0.0 )
          {
            return false;
          }
      }
  }

  return isRegularEveryPoint(a_region, a_domain, a_origin, a_dx);
}

bool GeometryShop::isRegularEveryPoint(const Box&           a_region,
                                       const ProblemDomain& a_domain,
                                       const RealVect&      a_origin,
                                       const Real&          a_dx) const
{
  CH_TIME("GeometryShop::isRegularEveryPoint");
  // set a vectDx
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
              physCorner[idir] = vectDx[idir]*corner[idir] + a_origin[idir];
            }

          if (m_stlIF == NULL)
          {
            // If the implicit function value is positive then the current
            // corner is outside the domain
            Real functionValue = m_implicitFunction->value(physCorner);

            if (functionValue > 0.0 )
              {
                return false;
              }
          }
          else
          {
            if (!m_STLBoxSet)
            {
              m_stlIF->getExplorer()->Explore(a_domain.domainBox(),a_domain,a_origin,vectDx);
              m_STLBoxSet = true;
            }

            MayDay::Error("STL not implemented");
          }
        }
      bit.reset();
      ++bit;
    }

  return true;
}

bool GeometryShop::isIrregular(const Box&           a_region,
                               const ProblemDomain& a_domain,
                               const RealVect&      a_origin,
                               const Real&          a_dx) const
{

  CH_TIME("GeometryShop::isIrregular");
  // set a vectDx
  RealVect vectDx = m_vectDx;

  // first check any of the Box corners is outside, and return false
  // right away. (bvs)
  RealVect physCorner;
  IntVect lo = a_region.smallEnd();
  IntVect len = a_region.size();
  for (int idir = 0; idir < CH_SPACEDIM; ++idir)
    {
      physCorner[idir] = vectDx[idir]*lo[idir] + a_origin[idir];
    }
  Real originVal = m_implicitFunction->value(physCorner);

  Box unitBox(IntVect::Zero, IntVect::Unit);
  for (BoxIterator bit(unitBox); bit.ok(); ++bit)
    {
      IntVect current = lo + len*bit();
      for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
          physCorner[idir] = vectDx[idir]*current[idir] + a_origin[idir];
        }

      if (m_stlIF == NULL)
      {
        Real functionValue = m_implicitFunction->value(physCorner);
        if (functionValue * originVal < 0.0 )
          {
            return true;
          }
      }
      else
      {
        if (!m_STLBoxSet)
        {
          m_stlIF->getExplorer()->Explore(a_domain.domainBox(),a_domain,a_origin,vectDx);
          m_STLBoxSet = true;
        }

        MayDay::Error("STL not implemented");
      }
    }

  // return isIrregularEveryPoint(a_region, a_domain, a_origin, a_dx, originVal);
  return !(isRegularEveryPoint(a_region, a_domain, a_origin, a_dx) ||
           isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx));
}

bool GeometryShop::isIrregularEveryPoint(const Box&           a_region,
                                         const ProblemDomain& a_domain,
                                         const RealVect&      a_origin,
                                         const Real&          a_dx,
                                         const Real&          a_originVal) const
{
  CH_TIME("GeometryShop::isIrregularEveryPoint");
  // set a vectDx
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
              physCorner[idir] = vectDx[idir]*corner[idir] + a_origin[idir];
            }

          if (m_stlIF == NULL)
          {
            // If the implicit function value is positive then the current
            // corner is outside the domain
            Real functionValue = m_implicitFunction->value(physCorner);

            if (functionValue * a_originVal < 0.0 )
              {
                return true;
              }
          }
          else
          {
            if (!m_STLBoxSet)
            {
              m_stlIF->getExplorer()->Explore(a_domain.domainBox(),a_domain,a_origin,vectDx);
              m_STLBoxSet = true;
            }

            MayDay::Error("STL not implemented");
          }
        }
      bit.reset();
      ++bit;
    }

  return false;
}

bool GeometryShop::isCovered(const Box&           a_region,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx) const
{
  CH_TIME("GeometryShop::isCovered");

  // set a vectDx
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
          physCorner[idir] = vectDx[idir]*current[idir] + a_origin[idir];
        }

      if (m_stlIF == NULL)
      {
        Real functionValue = m_implicitFunction->value(physCorner);
        if (functionValue < 0.0 )
          {
            return false;
          }
      }
      else
      {
        if (!m_STLBoxSet)
        {
          m_stlIF->getExplorer()->Explore(a_domain.domainBox(),a_domain,a_origin,vectDx);
          m_STLBoxSet = true;
        }

        MayDay::Error("STL not implemented");
      }
    }

  return isCoveredEveryPoint(a_region, a_domain, a_origin, a_dx);
}

bool GeometryShop::isCoveredEveryPoint(const Box&           a_region,
                                       const ProblemDomain& a_domain,
                                       const RealVect&      a_origin,
                                       const Real&          a_dx) const
{
  CH_TIME("GeometryShop::isCoveredEveryPoint");
  // set a vectDx
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
              physCorner[idir] = vectDx[idir]*corner[idir] + a_origin[idir];
            }

          if (m_stlIF == NULL)
          {
            // If the implicit function value is positive then the current
            // corner is outside the domain
            Real functionValue = m_implicitFunction->value(physCorner);

            if (functionValue < 0.0 )
              {
                return false;
              }
          }
          else
          {
            if (!m_STLBoxSet)
            {
              m_stlIF->getExplorer()->Explore(a_domain.domainBox(),a_domain,a_origin,vectDx);
              m_STLBoxSet = true;
            }

            MayDay::Error("STL not implemented");
          }
        }
      bit.reset();
      ++bit;
    }

  return true;
}

GeometryService::InOut GeometryShop::InsideOutside(const Box&           a_region,
                                                   const ProblemDomain& a_domain,
                                                   const RealVect&      a_origin,
                                                   const Real&          a_dx) const
{
  CH_TIME("GeometryShop::InsideOutside");

  GeometryService::InOut rtn;

  if (m_implicitFunction->fastIntersection(a_region, a_domain, a_origin, a_dx))
    {
      rtn = m_implicitFunction->InsideOutside(a_region, a_domain, a_origin, a_dx);
    }
  else
    {
      RealVect vectDx;

      // All corner indices for the current box
      Box allCorners(a_region);
      allCorners.surroundingNodes();

      RealVect physCorner(allCorners.smallEnd());

      if (m_vectDx[0] != 0.0)
      {
        vectDx[0] = a_dx;
        for (int idir = 1; idir < SpaceDim; idir++)
        {
          vectDx[idir] = vectDx[0] * m_vectDx[idir] / m_vectDx[0];
        }
      }
      else
      {
        vectDx = a_dx * RealVect::Unit;
      }

      STLExplorer* stlExplorer = NULL;
      if (m_stlIF != NULL)
      {
        if (!m_STLBoxSet)
        {
          m_stlIF->getExplorer()->Explore(a_domain.domainBox(),a_domain,a_origin,vectDx);
          m_STLBoxSet = true;
        }

        stlExplorer = m_stlIF->getExplorer();
      }

      physCorner *= vectDx;
      physCorner += a_origin;

      Real firstValue;
      Real firstSign;

      if (m_stlIF == NULL)
      {
        firstValue = m_implicitFunction->value(physCorner);
      }
      else
      {
        bool in;
        stlExplorer->GetPointInOut(allCorners.smallEnd(),in);

        if (in)
        {
          firstValue = -1.0;
        }
        else
        {
          firstValue =  1.0;
        }
      }
      
      firstSign  = copysign(1.0, firstValue);

      if ( firstSign < 0 )
        {
          rtn = GeometryService::Regular;
        }
      else
        {
          rtn = GeometryService::Covered;
        }

      BoxIterator bit(allCorners);

      for (; bit.ok(); ++bit)
        {
          // Current corner
          IntVect corner = bit();

          Real functionValue;
          Real functionSign;

          if (m_stlIF == NULL)
          {
            // Compute physical coordinate of corner
            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
              {
                physCorner[idir] = vectDx[idir]*corner[idir] + a_origin[idir];
              }

            // If the implicit function value is positive then the current
            // corner is outside the domain
            functionValue = m_implicitFunction->value(physCorner);
          }
          else
          {
            bool in;
            stlExplorer->GetPointInOut(corner,in);

            if (in)
            {
              functionValue = -1.0;
            }
            else
            {
              functionValue =  1.0;
            }
          }

          functionSign = copysign(1.0, functionValue);

          if (functionValue == 0 || firstValue == 0)
            {
              if (functionSign * firstSign < 0)
                {
                  rtn = GeometryService::Irregular;
                  return rtn;
                }
            }
          if (functionValue * firstValue < 0.0 )
            {
              rtn = GeometryService::Irregular;
              return rtn;
            }
        }
    }

  return rtn;
}

/**********************************************/
/*********************************************/
void
GeometryShop::fillGraph(BaseFab<int>        & a_regIrregCovered,
                        Vector<IrregNode>   & a_nodes,
                        const Box           & a_validRegion,
                        const Box           & a_ghostRegion,
                        const ProblemDomain & a_domain,
                        const RealVect      & a_origin,
                        const Real          & a_dx) const
{
  CH_TIMERS("GeometryShop::fillGraph");
  CH_TIMER("part1",p1);
  CH_TIMER("part2",p2);
  CH_TIMER("part3",p3);
  CH_TIMER("part4",p4);

  CH_START(p1);
  CH_assert(a_domain.contains(a_ghostRegion));
  RealVect vectDx;
  Real thrshd = m_thrshdVoF;
  // if (thrshd > 0)
  //   pout() << "GeometryShop:: Using thrshd: " << thrshd << endl;
  if (m_vectDx == RealVect::Zero)
    {
      vectDx = a_dx*RealVect::Unit;
    }
  else
    {
      vectDx = m_vectDx;
    }

  PolyGeom::setVectDx(vectDx);
  IntVectSet ivsirreg = IntVectSet(DenseIntVectSet(a_ghostRegion, false));
  IntVectSet ivsdrop  = IntVectSet(DenseIntVectSet(a_ghostRegion, false));// CP
  long int numCovered=0, numReg=0, numIrreg=0;
  CH_STOP(p1);

  CH_START(p2);
  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
    {
      const IntVect iv =bit();
      Box miniBox(iv, iv);
      GeometryService::InOut inout = InsideOutside(miniBox, a_domain, a_origin, a_dx);

      if (inout == GeometryService::Covered)
        {
          // set covered cells to -1
          a_regIrregCovered(iv, 0) = -1;
          numCovered++;
        }
      else if (inout == GeometryService::Regular)
        {
          // set regular cells to 1
          a_regIrregCovered(iv, 0) =  1;
          numReg++;
        }
      else
        {
          // set irregular cells to 0
          a_regIrregCovered(iv, 0) =  0;
          if (a_validRegion.contains(iv))
            {
              ivsirreg |= iv;
              numIrreg++;
            }
        }
    }
  // pout()<< "GeometryShop:: Counting cells:  " << numCovered<< "  "<< numReg<< "  "<< numIrreg  <<endl;
  CH_STOP(p2);

  CH_START(p3);
  // now loop through irregular cells and make nodes for each  one.
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
                          vectDx,
                          ivsit());


          IrregNode newNode;
          newNode.m_cell          = ivsit();
          newNode.m_volFrac       = volFrac;
          newNode.m_cellIndex     = 0;
          newNode.m_volCentroid   = volCentroid;
          newNode.m_bndryCentroid = bndryCentroid;
          // if (thrshd == 0.)//begin treb
          //   {
          //     //this piece of code successfully removes volFrac=1 cells where EB cuts the vertex
          //     //this piece of code cannot be used with CP's small volFrac removal below
          //     //this piece of code does not work for removal of volFrac << 1 because of regular next to covered
          //     bool isIrregNode;
          //     if ((volFrac < -thrshd || volFrac > thrshd) && (volFrac < 1.-thrshd || volFrac > 1.+thrshd))
          //       // if (volFrac < 1.)
          //       {
          //         isIrregNode = true;
          //       }
          //     else
          //       {
          //         isIrregNode = false;
          //       }
      
          //     for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          //       {
          //         int loNodeInd = newNode.index(faceDir, Side::Lo);
          //         int hiNodeInd = newNode.index(faceDir, Side::Hi);
          //         newNode.m_arc[loNodeInd]          = loArc[faceDir];
          //         newNode.m_arc[hiNodeInd]          = hiArc[faceDir];
          //         newNode.m_areaFrac[loNodeInd]     = loAreaFrac[faceDir];
          //         newNode.m_areaFrac[hiNodeInd]     = hiAreaFrac[faceDir];
          //         newNode.m_faceCentroid[loNodeInd] = loFaceCentroid[faceDir];
          //         newNode.m_faceCentroid[hiNodeInd] = hiFaceCentroid[faceDir];
          //         if (!isIrregNode)//only go through this logic if a covered or regular cell has chance of being irregular????
          //           {
          //             //covered cell with no arcs in all directions
          //             if ((volFrac == 0.) && (loArc[faceDir].size() == 0) && (hiArc[faceDir].size() == 0))
          //               {
          //                 isIrregNode = false;
          //               }
          //             //regular cell with exactly one arc one lo and hi sides in all directions
          //             else if ((volFrac > 1.-thrshd && volFrac < 1.+thrshd) && (loArc[faceDir].size() == 1) && (hiArc[faceDir].size() == 1))
          //               // else if (volFrac == 1. && (loArc[faceDir].size() == 1) && (hiArc[faceDir].size() == 1))
          //               {
          //                 for (int numFace = 0; numFace < loAreaFrac[faceDir].size(); numFace++)
          //                   {
          //                     if (loAreaFrac[faceDir][numFace] > 1.-thrshd && loAreaFrac[faceDir][numFace] < 1.+thrshd)
          //                       {
          //                         isIrregNode = false;
          //                       }
          //                     else
          //                       {
          //                         isIrregNode = true;
          //                       }
          //                   }
          //                 for (int numFace = 0; numFace < hiAreaFrac[faceDir].size(); numFace++)
          //                   {
          //                     if (hiAreaFrac[faceDir][numFace] > 1.-thrshd && hiAreaFrac[faceDir][numFace] < 1.+thrshd)
          //                       {
          //                         isIrregNode = false;
          //                       }
          //                     else
          //                       {
          //                         isIrregNode = true;
          //                       }
          //                   }
          //               }
          //             //if none of those fit for all directions, then irregular, and don't come back through this logic
          //             else
          //               {
          //                 isIrregNode = true;
          //               }
          //           }
          //       }

          //     if (isIrregNode)
          //       {
          //         a_nodes.push_back(newNode);
          //       }
          //     else
          //       {
          //         if ((volFrac > 1.-thrshd && volFrac < 1.+thrshd))
          //           // if (volFrac == 1.)
          //           {
          //             a_regIrregCovered(ivsit(), 0) =  1;
          //             pout() << "Removing regular vof " << vof << " from irreg node list" << endl;
          //           }
          //       }
          //   }//end treb
          // else//begin CP
          //   {
              if (thrshd > 0. && volFrac < thrshd)
                {
                  ivsdrop |= ivsit();
                  a_regIrregCovered(ivsit(), 0) = -1;
                  if (m_verbosity > 2)
                    {
                      pout() << "Removing vof " << vof << " with volFrac " << volFrac << endl;
                    }
                }//CP record these nodes to be removed
              else
                {
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
            // }//end CP
    } // end loop over cells in the box
  CH_STOP(p3);

  CH_START(p4);
  // CP: fix sweep that removes cells with volFrac less than a certain threshold
  for (IVSIterator ivsit(ivsdrop); ivsit.ok(); ++ivsit)
    {
      VolIndex vof(ivsit(), 0);
      IntVect iv = vof.gridIndex();

      // multiple nodes in a gridcell location?
      // where is this guy in a_nodes?--search in m_cell?
      // how to access this guy's neighbor?

      // newNode.m_cell          = ivsit();
      // newNode.m_volFrac       = volFrac;
      // newNode.m_cellIndex     = 0;
      // newNode.m_volCentroid   = volCentroid;
      // newNode.m_bndryCentroid = bndryCentroid;
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int isign = sign(sit());
              IntVect otherIV = iv + isign*BASISV(faceDir);
              if (a_validRegion.contains(otherIV))
                {
                  if (a_regIrregCovered(otherIV,0) == 0)
                    {
                      // i am in the case where the other cell
                      // is also irregular.   I just made a previously
                      // irregular cell covered so I have to check to
                      // see if it had any faces pointed this way.
                      int inode = -1;
                      bool found = false;
                      for (int ivec = 0; ivec < a_nodes.size() && ! found; ivec++)
                        {
                          if (a_nodes[ivec].m_cell == otherIV)
                            {
                              inode = ivec;
                              found = true;
                            }
                        }
                      if (!found && a_validRegion.contains(otherIV))
                        {
                          MayDay::Error("something wrong in our logic");
                        }
                      if (found)
                        {
                          int arcindex = a_nodes[inode].index(faceDir, flip(sit()));
                          a_nodes[inode].m_arc[         arcindex].resize(0);
                          a_nodes[inode].m_areaFrac[    arcindex].resize(0);
                          a_nodes[inode].m_faceCentroid[arcindex].resize(0);
                        }
                    }
                  else if (a_regIrregCovered(otherIV,0) == 1)
                    // CP This is a very peculiar case which has so far
                    // not happened. A irregular cell with tiny cell volFrac
                    // is connected to a regular cell
                    // may not work well. Need to be debugged before it is used

                    {
                      // MayDay::Error("need to be debugged before this branch is used!!!");
                      VolIndex vof(otherIV, 0);
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
                                          ivsirreg, // CP this one????
                                          vof,
                                          a_domain,
                                          a_origin,
                                          a_dx,
                                          vectDx,
                                          otherIV); // CP
                      IrregNode newNode;
                      // case where neighbor is regular.  need to make
                      // a new node with IrregNode newNode;

                      newNode.m_cell          = otherIV;
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
                      a_regIrregCovered(otherIV,0) = 0;

                    }//else if
                }//valid region
            }//sit
        }//facedir
    }//ivsdrop
  CH_STOP(p4);
}

/**********************************************/
/*********************************************/
void
GeometryShop::computeVoFInternals(Real&               a_volFrac,
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
  CH_TIME("GeometryShop::computeVoFInternals");

  // need maxDx to properly scale a_bndryArea
  Real maxDx = 0.0;

  for (int idir = 0; idir <SpaceDim; ++idir)
    {
      CH_assert(a_vectDx[idir] > 0);
      if (a_vectDx[idir] >maxDx)
        {
          maxDx = a_vectDx[idir];
        }
    }

  if (SpaceDim == 2)
    {
      // In 2D a vof is a faceMo = edgeMo[4],boundary length and normal vector
      faceMo Face;
      edgeMo edges[4];

      bool faceCovered;
      bool faceRegular;
      bool faceDontKnow;
      int faceNormal = 2;

      // get edgeType and intersection points
      edgeData2D(edges,
                 faceCovered,
                 faceRegular,
                 faceDontKnow,
                 a_dx,
                 a_vectDx,
                 a_iv,
                 a_domain,
                 a_origin);

      CH_assert(faceRegular || faceCovered || faceDontKnow);
      CH_assert((!(faceRegular && faceCovered)) && (!(faceRegular && faceDontKnow)) && (!(faceDontKnow && faceCovered)));
      // define the faceMo
      Face.define(edges,faceNormal,faceCovered,faceRegular,faceDontKnow);

      Moments geom;

      // answer0,answer1 are vectors whose components contain geometric information
      Vector<Real> answer0;
      Vector<Real> answer1;

      int order = 0;
      answer0 = geom.momentCalc2D(order,Face);

      // extract the info from answer0 and answer1

      // volfrac
      a_volFrac = answer0[0];

      order = 1;
      answer1 = geom.momentCalc2D(order,Face);

      // centroid
      for (int idir=0; idir<SpaceDim;++idir)
        {
          a_volCentroid[idir] = answer1[SpaceDim-1-idir];
        }

      if (a_volFrac <= 0.0)
        {
          a_volCentroid = RealVect::Zero;
        }
      else
        {
          a_volCentroid /= a_volFrac;
        }

      // normal
      Real normalVec[SpaceDim];
      Face.getNormal(normalVec);
      for (int idir = 0;idir < SpaceDim;++idir)
        {
          a_normal[idir] = normalVec[idir];
        }

      for (int idir = 0;idir < SpaceDim;++idir)
        {
          // (nx,ny)->(nxdy,nydx)
          a_normal[idir] = normalVec[idir]*a_vectDx[1 - idir];
        }
      Real anisBd = 0.0;
      for (int idir = 0;idir < SpaceDim;++idir)
        {
          anisBd += a_normal[idir]*a_normal[idir];
        }
      anisBd = sqrt(anisBd);
      if (anisBd !=0.0)
        {
          a_normal /= anisBd;
        }

      // compute bndryArea and bndryCentroid
      a_bndryArea = Face.getBdLength();
      if (a_bndryArea <= 0.0)
        {
          a_bndryCentroid = RealVect::Zero;
          a_bndryArea = 0.0;
        }
      else
        {
          for (int idir = 0;idir < SpaceDim; ++idir)
            {
              a_bndryCentroid[idir] = answer0[SpaceDim-idir]/a_bndryArea;
            }
          a_bndryArea *= anisBd;
          a_bndryArea /= maxDx;
        }

      for (int edgeNormal = 0; edgeNormal < SpaceDim; ++edgeNormal)
        {
          // loside
          bool coveredLo = edges[edgeNormal*2].isCovered();
          if (coveredLo)
            {
              a_loArc[edgeNormal].resize(0);
              a_loFaceCentroid[edgeNormal].resize(0);
              a_loAreaFrac[edgeNormal].resize(0);
            }
          else if (!coveredLo)
            {
              a_loArc[edgeNormal].resize(1);
              a_loFaceCentroid[edgeNormal].resize(1);
              a_loAreaFrac[edgeNormal].resize(1);
              IntVect otherIV = a_iv;
              otherIV[edgeNormal] -= 1;
              if (a_domain.contains(otherIV))
                {
                  int otherCellIndex;
                  if (a_ivsIrreg.contains(otherIV))
                    {
                      otherCellIndex = 0;
                    }
                  else
                    {
                      // arc to regular cell
                      otherCellIndex = -2;
                    }
                  a_loArc[edgeNormal][0]=otherCellIndex;
                }
              else if (!a_domain.contains(otherIV))
                {
                  // boundary arcs always -1
                  a_loArc[edgeNormal][0] = -1;
                }
              a_loFaceCentroid[edgeNormal][0] = edges[edgeNormal*2].getEdgeCentroid();
              a_loAreaFrac[edgeNormal][0] = edges[edgeNormal*2].getEdgeLength();
            }

          // hiside
          bool coveredHi = edges[edgeNormal*2+1].isCovered();
          if (coveredHi)
            {
              a_hiArc[edgeNormal].resize(0);
              a_hiFaceCentroid[edgeNormal].resize(0);
              a_hiAreaFrac[edgeNormal].resize(0);
            }
          else if (!coveredHi)
            {
              a_hiArc[edgeNormal].resize(1);
              a_hiFaceCentroid[edgeNormal].resize(1);
              a_hiAreaFrac[edgeNormal].resize(1);
              IntVect otherIV = a_iv;
              otherIV[edgeNormal] += 1;
              if (a_domain.contains(otherIV))
                {
                  int otherCellIndex;
                  if (a_ivsIrreg.contains(otherIV))
                    {
                      otherCellIndex = 0;
                    }
                  else
                    {
                      // arc to regular cell
                      otherCellIndex = -2;
                    }
                  a_hiArc[edgeNormal][0] = otherCellIndex;
                }
              else if (!a_domain.contains(otherIV))
                {
                  // boundary arcs always -1
                  a_hiArc[edgeNormal][0] = -1;
                }
              a_hiFaceCentroid[edgeNormal][0] = edges[edgeNormal*2+1].getEdgeCentroid();
              a_hiAreaFrac[edgeNormal][0] = edges[edgeNormal*2+1].getEdgeLength();

            }

        }
    }

  if (SpaceDim==3)
    {
      // 1) using the intvect, build up the classes in Moments: edgeMO,
      //    faceMo and finally vofMo
      // 2) check for covered or regular faces
      // 3) call  momentCalc3D
      // 4) keep track of what the output means and fill the variables
      //    requested

      faceMo Faces[6];
      int index = -1;

      for (int faceNormal = 0;faceNormal < SpaceDim;++faceNormal)
        {
          for (int hiLoFace = 0;hiLoFace < 2;++hiLoFace)
            {
              index += 1;
              edgeMo edges[4];
              bool faceCovered;
              bool faceRegular;
              bool faceDontKnow;
              edgeData3D(edges,
                         faceCovered,
                         faceRegular,
                         faceDontKnow,
                         hiLoFace,
                         faceNormal,
                         a_dx,
                         a_vectDx,
                         a_iv,
                         a_domain,
                         a_origin);

              CH_assert(faceRegular || faceCovered || faceDontKnow);
              CH_assert((!(faceRegular && faceCovered)) && (!(faceRegular && faceDontKnow)) && (!(faceDontKnow && faceCovered)));
              Faces[index].define(edges,faceNormal,faceCovered,faceRegular,faceDontKnow);

              // if the face is covered we will deal with it later
              if (!Faces[index].isCovered())
                {
                  Moments geom;
                  // answer0 and answer1 have all the geometric facts
                  Vector<Real> answer0;
                  Vector<Real> answer1;

                  int order = 0;
                  answer0 = geom.momentCalc2D(order,Faces[index]);

                  // area of this face
                  Real area=answer0[0];

                  order = 1;
                  answer1 = geom.momentCalc2D(order,Faces[index]);

                  // first we get the centroid of the face
                  RealVect faceCentroid;
                  for (int idir = 0; idir<SpaceDim;++idir)
                    {
                      faceCentroid[idir] = answer1[SpaceDim-1-idir];
                    }

                  if (area > 0.0)
                    {
                      faceCentroid /= area;
                    }
                  else
                    {
                      faceCentroid = RealVect::Zero;
                    }

                  // record these facts in member data
                  Faces[index].setFaceCentroid(faceCentroid);
                  Faces[index].setFaceArea(area);
                }

              else if (Faces[index].isCovered())
                {
                  RealVect faceCentroid = RealVect::Zero;
                  Faces[index].setFaceCentroid(faceCentroid);
                  Real area = 0.0;
                  Faces[index].setFaceArea(area);
                }
            }
        }
#if 1
      // iterate over faces recalculating face area
      // face order is xLo,xHi,yLo,yHi,zLo,zHi
      for (int iFace = 0; iFace < 2*SpaceDim; ++iFace)
        {
          // recalculate face area
          faceMo& face = Faces[iFace];

          // collect exactly two irregular edges or mayday.
          Vector<RealVect> crossingPt;

          Vector<int> cPtHiLo;
          Vector<int> edgeHiLo;
          Vector<int> edgeDir;

          if (!(face.isCovered()) && !(face.isRegular()))
            {
              for ( int iEdge = 0; iEdge < 4; ++iEdge)
                {
                  const edgeMo& curEdge = face.retrieveEdge(iEdge);

                  if (curEdge.dontKnow())
                    {
                      // loPt of edge
                      RealVect loPt = curEdge.getLo();

                      // hiPt of edge
                      RealVect hiPt = curEdge.getHi();

                      // direction that varies over edge
                      int direction = curEdge.direction();
                      bool intersectLo = curEdge.getIntersectLo();

                      // for irregular or regular edges at most one pt away from corner
                      if (intersectLo)
                        {
                          crossingPt.push_back(loPt);
                          cPtHiLo.push_back(-1);

                          int hilo = iEdge % 2;
                          edgeHiLo.push_back(hilo);

                          edgeDir.push_back(direction);
                        }
                      else // forced by dontknow
                        {
                          crossingPt.push_back(hiPt);
                          cPtHiLo.push_back(1);

                          int hilo = iEdge % 2;
                          edgeHiLo.push_back(hilo);

                          edgeDir.push_back(direction);
                        }
                    }
                }
            }

          if (crossingPt.size() == 2)
            {
              // get midpoint of line connecting intersection points
              RealVect midPt = crossingPt[0];
              midPt += crossingPt[1];
              midPt *= 0.5;

              // faceNormal
              int faceNormal = Faces[iFace].getFaceNormal();

              // directions over which the face varies
              int dir1 = (faceNormal + 1) % SpaceDim;
              int dir2 = (faceNormal + 2) % SpaceDim;

              // find max of (deltaDir1,deltaDir2)
              Real deltaDir1 = Abs(crossingPt[0][dir1] - crossingPt[1][dir1]);
              Real deltaDir2 = Abs(crossingPt[0][dir2] - crossingPt[1][dir2]);

              // maxDir will be the direction of integration (independent variable)
              // minDir will the direction in which the integrand varies(dependent variable)
              int maxDir;
              int minDir;
              if (deltaDir1 > deltaDir2)
                {
                  maxDir = dir1;
                  minDir = dir2;
                }
              else
                {
                  maxDir = dir2;
                  minDir = dir1;
                }
              // flip area?
              bool complementArea;

              CH_assert (cPtHiLo.size() == 2);

              // loEdge-loEdge
              if (edgeHiLo[0] == 0 && edgeHiLo[1] == 0)
                {
                  // both crossingPts must be Hi or both must be Lo
                  CH_assert((cPtHiLo[0] == 1 && cPtHiLo[1] == 1) ||
                         (cPtHiLo[0] == -1 && cPtHiLo[1] == -1));

                  // prismArea gives triangle

                  // if the cPtHiLo[0]= hiPt then one wants triangle
                  if (cPtHiLo[0] == 1)
                    {
                      complementArea = false;
                    }
                  else
                    {
                      complementArea = true;
                    }
                }

              // hiEdge-hiEdge
              else if (edgeHiLo[0] == 1 && edgeHiLo[1] == 1)
                {
                  // both crossingPts must be Hi or both must be Lo
                  CH_assert((cPtHiLo[0] == 1 && cPtHiLo[1] == 1) ||
                         (cPtHiLo[0] == -1 && cPtHiLo[1] == -1));
                  // prismArea gives the trapezoid

                  // cPtHiLo[0] == Lo => one wants the triangle area
                  if (cPtHiLo[0] == -1)
                    {
                      complementArea = true;
                    }
                  else
                    {
                      complementArea = false;
                    }
                }

              // hiEdge-loEdge
              else if (edgeHiLo[0] == 1 && edgeHiLo[1] == 0)
                {
                  // cpPtHiLo must be the same or the opposite of edgeHiLo
                  CH_assert((cPtHiLo[0] == 1 && cPtHiLo[1] == -1) ||
                         (cPtHiLo[0] == -1 && cPtHiLo[1] == 1));
                  // maxDir > minDir =>prismArea gives trapezoid
                  // maxDir < minDir =>prismArea gives triangle

                  // (cPtHiLo[1] == 1) => one wants trapezoid
                  if (cPtHiLo[1] == 1)
                    {
                      if (maxDir < minDir)
                        {
                          complementArea = true;
                        }
                      else
                        {
                          complementArea = false;
                        }
                    }
                  else
                    // (cPtHiLo[1] == -1) => one wants triangle
                    {
                      if (maxDir < minDir)
                        {
                          complementArea = false;
                        }
                      else
                        {
                          complementArea = true;
                        }
                    }
                }

              // loEdge-hiEdge
              else if (edgeHiLo[0] == 0 && edgeHiLo[1] == 1)
                {
                  // triangle + triangle complement or two trapezoids comprise this case

                  // two trapezoids
                  if (cPtHiLo[1] == 1 && cPtHiLo[0] == 1)
                    {
                      CH_assert(edgeDir[0] == edgeDir[1]);
                      complementArea = false;
                    }
                  else if (cPtHiLo[1] == -1 && cPtHiLo[0] == -1)
                    {
                      CH_assert(edgeDir[0] == edgeDir[1]);
                      complementArea = true;
                    }
                  // triangle + triangle complement
                  // if the cPtHiLo[0]= loPt then one wants triangle
                  else if (cPtHiLo[0] == -1 && cPtHiLo[1] == 1 )
                    {
                      CH_assert(edgeDir[0] != edgeDir[1]);
                      // if maxDir < minDir prismArea gives trapezoid
                      // if maxDir > minDir prismArea gives triangle
                      if (maxDir < minDir)
                        {
                          complementArea = true;
                        }
                      else
                        {
                          complementArea = false;
                        }
                    }

                  // cPtHiLo[0]= hiPt => one wants trapezoid
                  else if (cPtHiLo[0] == 1 && cPtHiLo[1] == -1)
                    {
                      if (maxDir < minDir)
                        {
                          complementArea = false;
                        }
                      else
                        {
                          complementArea = true;
                        }
                    }
                }

              else
                {
                  MayDay::Abort("cPtDir or mindir or maxDir not set correctly");
                }

              // segLo is the lo end of segment within face[iEdge] for Brent Rootfinder
              RealVect segLo = midPt;
              segLo[minDir] = -0.5;

              // segHi is the hi end of segment within face[iEdge] for Brent Rootfinder
              RealVect segHi = midPt;
              segHi[minDir] = 0.5;

              // put segLo and segHi in physical coordinates
              RealVect physSegLo;
              RealVect physSegHi;
              RealVect physMidPt;
              for (int idir = 0; idir < SpaceDim; ++idir)
                {
                  physSegLo[idir] = a_vectDx[idir]*(segLo[idir] + a_iv[idir] + 0.5) + a_origin[idir];
                  physSegHi[idir] = a_vectDx[idir]*(segHi[idir] + a_iv[idir] + 0.5) + a_origin[idir];
                  physMidPt[idir] = a_vectDx[idir]*(midPt[idir] + a_iv[idir] + 0.5) + a_origin[idir];
                }

              // find upDir
              pair<int,Side::LoHiSide> upDir;

              // physIntercept is along the segment[physSegLo,physSegHi]
              // this segment passes through midPt with direction minDir
              Real physIntercept;
              bool dropOrder = false;


              if (m_stlIF == NULL)
              {
                Real fLo = m_implicitFunction->value(physSegLo);
                Real fHi = m_implicitFunction->value(physSegHi);

                // This guards against the "root must be bracketed" error
                // by dropping order
                if (fLo*fHi > 0.0)
                  {
                    dropOrder = true;
                  }
                else
                  {
                    physIntercept = BrentRootFinder(physSegLo, physSegHi, minDir);
                  }
              }
              else
              {
                dropOrder = true;
              }

              if (!dropOrder)
                {
                  // put physIntercept into relative coordinates
                  Real intercept = physIntercept - a_origin[minDir];
                  intercept  /= a_vectDx[minDir];
                  intercept -= (a_iv[minDir]+0.5);

                  // push_back third pt onto crossingPt
                  crossingPt.push_back(midPt);
                  crossingPt[2][minDir] = intercept;

                  // integrate w.r.t xVec using Prismoidal Rule
                  RealVect xVec;
                  RealVect yVec;

                  // the order of (xVec,yVec) will be sorted out in PrismoidalAreaCalc
                  xVec[0] = crossingPt[0][maxDir];
                  xVec[1] = crossingPt[2][maxDir];
                  xVec[2] = crossingPt[1][maxDir];

                  yVec[0] = crossingPt[0][minDir];
                  yVec[1] = crossingPt[2][minDir];
                  yVec[2] = crossingPt[1][minDir];

                  // Prismoidal's rule
                  Real area = PrismoidalAreaCalc(xVec,yVec);

                  // Only use area if it is valid
                  if (area >= 0.0 && area <= 1.0)
                  {
                    // assign area to this value or (1 - this value)
                    if (complementArea)
                      {
                        area = 1.0 - area;
                      }

                    Faces[iFace].setFaceArea(area);
                  }
                }
            }
        }
#endif
      // fill in some arguments of computeVofInternals for the faces
      for (int faceNormal = 0;faceNormal < SpaceDim;++faceNormal)
        {
          bool coveredLo = Faces[faceNormal*2].isCovered();
          if (coveredLo)
            {
              a_loArc[faceNormal].resize(0);
              a_loFaceCentroid[faceNormal].resize(0);
              a_loAreaFrac[faceNormal].resize(0);
            }
          else if (!coveredLo)
            {
              a_loArc[faceNormal].resize(1);
              a_loFaceCentroid[faceNormal].resize(1);
              a_loAreaFrac[faceNormal].resize(1);
              IntVect otherIV = a_iv;
              otherIV[faceNormal] -= 1;

              if (a_domain.contains(otherIV))
                {
                  int otherCellIndex;
                  if (a_ivsIrreg.contains(otherIV))
                    {
                      otherCellIndex = 0;
                    }
                  else
                    {
                      // arc to regular cell
                      otherCellIndex = -2;
                    }
                  a_loArc[faceNormal][0] = otherCellIndex;
                }
              else if (!a_domain.contains(otherIV))
                {
                  // boundary arcs always -1
                  a_loArc[faceNormal][0] = -1;
                }

              a_loFaceCentroid[faceNormal][0] = Faces[faceNormal*2].getFaceCentroid();

              a_loAreaFrac[faceNormal][0] = Faces[faceNormal*2].getFaceArea();

            }
          else
            {
              MayDay::Abort("is it coveredLo?");
            }
          bool coveredHi = Faces[faceNormal*2+1].isCovered();
          if (coveredHi)
            {
              a_hiArc[faceNormal].resize(0);
              a_hiFaceCentroid[faceNormal].resize(0);
              a_hiAreaFrac[faceNormal].resize(0);
            }
          else if (!coveredHi)
            {
              a_hiArc[faceNormal].resize(1);
              a_hiFaceCentroid[faceNormal].resize(1);
              a_hiAreaFrac[faceNormal].resize(1);
              IntVect otherIV = a_iv;
              otherIV[faceNormal] += 1;
              if (a_domain.contains(otherIV))
                {
                  int otherCellIndex;
                  if (a_ivsIrreg.contains(otherIV))
                    {
                      otherCellIndex = 0;
                    }
                  else
                    {
                      // arc to regular cell
                      otherCellIndex = -2;
                    }
                  a_hiArc[faceNormal][0] = otherCellIndex;
                }
              else if (!a_domain.contains(otherIV))
                {
                  // boundaryArcs always -1
                  a_hiArc[faceNormal][0] = -1;
                }
              a_hiFaceCentroid[faceNormal][0] = Faces[faceNormal*2+1].getFaceCentroid();
              a_hiAreaFrac[faceNormal][0] = Faces[faceNormal*2+1].getFaceArea();
            }
          else
            {
              MayDay::Abort("is it coveredHi?");
            }
        }

      // We have enough face data to construct the vof
      vofMo Vof;
      Vof.define(Faces);

      Moments geom;
      int order = 0;
      Vector<Real> answer0;
      answer0 = geom.momentCalc3D(order,Vof);

      a_volFrac = answer0[0];

      Vector<Real> answer1;
      order = 1;
      answer1 = geom.momentCalc3D(order,Vof);

      for (int idir=0; idir<SpaceDim;++idir)
        {
          a_volCentroid[idir] = answer1[SpaceDim-1-idir];
        }

      if (a_volFrac == 0.0)
        {
          a_volCentroid = RealVect::Zero;
        }
      else
        {
          a_volCentroid /= a_volFrac;

        }

      Real normalVec[SpaceDim];
      Vof.getNormal(normalVec);
      for (int idir = 0;idir < SpaceDim;++idir)
        {
          a_normal[idir] = normalVec[idir];
        }

      for (int idir = 0;idir < SpaceDim;++idir)
        {
          // (nx,ny,nz)->(nxdydz,nydxdz,nzdxdy)
          a_normal[idir] = normalVec[idir]*a_vectDx[((idir-1)*(idir-2))/2]*a_vectDx[2 - ((idir-1)*idir)/2];
        }
      Real anisBd = 0.0;
      for (int idir = 0;idir < SpaceDim;++idir)
        {
          anisBd += a_normal[idir]*a_normal[idir];
        }
      anisBd = sqrt(anisBd);
      if (anisBd !=0.0)
        {
          a_normal /= anisBd;
        }

      a_bndryArea = Vof.getBdArea();
      if (a_bndryArea > 0.0)
        {
          for (int idir = 0;idir < SpaceDim;++idir)
            {
              a_bndryCentroid[idir] = answer0[SpaceDim-idir]/a_bndryArea;
            }
          a_bndryArea *= anisBd;
          a_bndryArea /= (maxDx*maxDx);
        }
      else
        {
          a_bndryCentroid = RealVect::Zero;
          a_bndryArea = 0.0;
        }
    }

  // clipping
  // only report adjustments when discrepancy is above the threshold
  bool thisVofClipped = false;
  Real discrepancy = 0.0;
  Real volDiscrepancy = 0.0;
  // volFrac out of bounds
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
          MayDay::Warning(message);
        }
      // do the clipping
      thisVofClipped = true;
      a_volFrac = 0.0;
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
          MayDay::Warning(message);
        }
      // do the clipping
      thisVofClipped = true;
      a_volFrac = 1.0;
    }

  // area frac out of bounds
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_loAreaFrac[idir].size();num ++)
        {
          // lo frac too high
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
              if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                {
                  MayDay::Warning(message);
                }

              // do the clipping
              thisVofClipped = true;
              a_loAreaFrac[idir][num] = 1.0;
            }
          // lo frac too low
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
              if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                {
                  MayDay::Warning(message);
                }

              // do the clipping
              thisVofClipped = true;
              a_loAreaFrac[idir][num] = 0.0;
            }

        }
    }

  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_hiAreaFrac[idir].size();num ++)
        {
          // hi frac too high
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
              if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                {
                  MayDay::Warning(message);
                }

              // do the clipping
              thisVofClipped = true;
              a_hiAreaFrac[idir][num] = 1.0;
            }
          // hi frac too low
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
              if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
                {
                  MayDay::Warning(message);
                }

              // do the clipping
              thisVofClipped = true;
              a_hiAreaFrac[idir][num] = 0.0;
            }

        }

    }

  // bndry area out of bounds
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
      if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
        {
          MayDay::Warning(message);
        }

      // do the clipping
      thisVofClipped = true;
      a_bndryArea = 0.0;
    }

  if (a_bndryArea > sqrt(2.0))
    {
      discrepancy = Abs(sqrt(2.0) - a_bndryArea);
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
      if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
        {
          MayDay::Warning(message);
        }
      // do the clipping
      thisVofClipped = true;
      a_bndryArea = sqrt(2.0);
    }

  // volCentroid out of bounds
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
          if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
            {
              MayDay::Warning(message);
            }
          // do the clipping
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
          if (discrepancy > m_threshold&& volDiscrepancy>m_threshold)
            {
              MayDay::Warning(message);
            }
          // do the clipping
          thisVofClipped = true;
          a_volCentroid[idir] = -0.5;
        }

      // boundary centroid out of bounds
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
          if (discrepancy>m_threshold && volDiscrepancy>m_threshold)
            {
              MayDay::Warning(message);
            }
          // do the clipping
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
          if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
            {
              MayDay::Warning(message);
            }
          // do the clipping
          thisVofClipped = true;
          a_bndryCentroid[idir] = -0.5;
        }
    }
  // loFaceCentroid out of bounds
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_loFaceCentroid[idir].size();num ++)
        {
          for (int jdir = 0; jdir<SpaceDim; ++jdir)
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
                  if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
                    {
                      MayDay::Warning(message);
                    }
                  // do the clipping
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
                  if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
                    {
                      MayDay::Warning(message);
                    }
                  // do the clipping
                  thisVofClipped = true;
                  a_loFaceCentroid[idir][num][jdir] = -0.5;
                }
            }
        }
    }
  // hiFaceCentroid out of bounds
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      for (int num = 0; num < a_hiFaceCentroid[idir].size();num ++)
        {
          for (int jdir = 0; jdir<SpaceDim; ++jdir)
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
                  if (discrepancy > m_threshold && volDiscrepancy>m_threshold)
                    {
                      MayDay::Warning(message);
                    }
                  // do the clipping
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
                  if (discrepancy > m_threshold && volDiscrepancy > m_threshold)
                    {
                      MayDay::Warning(message);
                    }
                  // do the clipping
                  thisVofClipped = true;
                  a_hiFaceCentroid[idir][num][jdir] = -0.5;
                }
            }
        }
    }

  if (thisVofClipped)
    {
      GeometryShop* changedThis = (GeometryShop *) this;
      changedThis->m_numCellsClipped += 1;
    }
}

int GeometryShop::getNumCellsClipped()
{
  return m_numCellsClipped;
}

void GeometryShop::edgeData3D(edgeMo a_edges[4],
                              bool& a_faceCovered,
                              bool& a_faceRegular,
                              bool& a_faceDontKnow,
                              const int a_hiLoFace,
                              const int a_faceNormal,
                              const Real& a_dx,
                              const RealVect& a_vectDx,
                              const IntVect& a_iv,
                              const ProblemDomain& a_domain,
                              const RealVect& a_origin) const
{
  CH_TIME("GeometryShop::edgeData3D");
  a_faceRegular = true;
  a_faceCovered = true;
  a_faceDontKnow = false;

  int index = -1;

  // edge order is lexigraphic xLo,xHi,yLo,yHi
  for (int dom = 0; dom < 3; ++dom)
    {
      if (dom != a_faceNormal)
        {
          for (int lohi = 0; lohi < 2; ++lohi)
            {
              // which edge 0,1,2, or 3 in lexigraphic order is given by index
              index += 1;
              // range is the direction along which the edge varies
              int range = 3 - a_faceNormal - dom;

              RealVect LoPt;
              bool LoPtChanged = false;
              Real funcLo;

              RealVect HiPt;
              bool HiPtChanged = false;
              Real funcHi;

              // put LoPt in physical coordinates
              LoPt[a_faceNormal] = a_origin[a_faceNormal]+
                (a_iv[a_faceNormal]+a_hiLoFace)*a_vectDx[a_faceNormal];
              LoPt[dom] = a_origin[dom] + (a_iv[dom]+lohi)*a_vectDx[dom];
              LoPt[range] = a_origin[range] + (a_iv[range])*a_vectDx[range];

              // put HiPt in physical coordinates
              HiPt[a_faceNormal] = a_origin[a_faceNormal] +
                (a_iv[a_faceNormal]+a_hiLoFace)*a_vectDx[a_faceNormal];
              HiPt[dom] = a_origin[dom] + (a_iv[dom]+lohi)*a_vectDx[dom];
              HiPt[range] = a_origin[range] + (a_iv[range]+1)*a_vectDx[range];

              // find the midpoint
              RealVect MidPt = LoPt;
              MidPt += HiPt;
              MidPt /= 2.0;

              Real signHi;
              Real signLo;

              // placeholders for edgeType
              bool covered  = false;
              bool regular  = false;
              bool dontKnow = false;

              RealVect interceptPt = RealVect::Zero;

              if (m_stlIF == NULL)
              {
                funcHi = m_implicitFunction->value(HiPt);
                funcLo = m_implicitFunction->value(LoPt);
              }
              else
              {
                IntVect loIV(a_iv);
                loIV[a_faceNormal] += a_hiLoFace;
                loIV[dom]          += lohi;

                IntVect hiIV(a_iv);
                hiIV[a_faceNormal] += a_hiLoFace;
                hiIV[dom]          += lohi;
                hiIV[range]        += 1;

                CellEdge curEdge(loIV,hiIV);

                bool loIn,hiIn;
                m_stlIF->getExplorer()->GetCellEdgeIntersection(curEdge,interceptPt,loIn,hiIn);

                if (loIn)
                {
                  funcLo = -1.0;
                }
                else
                {
                  funcLo = 1.0;
                }

                if (hiIn)
                {
                  funcHi = -1.0;
                }
                else
                {
                  funcHi = 1.0;
                }
              }

              // For level set data negative -> in the fluid
              //                    positive -> out of the fluid
              signHi = -funcHi;
              signLo = -funcLo;

              edgeType(covered,regular,dontKnow,signHi,signLo);

              // now we know the boolean values so we can set the edge Hi  and Lo pts
              if (covered)
                {
                  a_faceRegular=false;

                  LoPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_vectDx[range];
                  LoPtChanged = true;

                  HiPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_vectDx[range];
                  HiPtChanged = true;
                }
              else if (regular)
                {
                  a_faceCovered = false;
                }
              else if (dontKnow)
                {
                  a_faceRegular = false;
                  a_faceCovered = false;
                  a_faceDontKnow = true;

                  Real intercept;
                  if (m_stlIF == NULL)
                  {
                    // find where the surface intersects the edge
                    intercept = BrentRootFinder(LoPt, HiPt, range);
                  }
                  else
                  {
                    intercept = interceptPt[range];
                  }

                  if (funcHi >= 0 && funcLo*funcHi <= 0)
                    {
                      HiPt[range] = intercept;
                      HiPtChanged = true;
                    }
                  else if (funcLo >= 0 && funcLo*funcHi <= 0)
                    {
                      LoPt[range] = intercept;
                      LoPtChanged = true;
                    }
                  else
                    {
                      MayDay::Abort("Bogus intersection calculated");
                    }
                }

              // put LoPt and HiPt in local coordinates
              if (a_hiLoFace == 0)
                {
                  LoPt[a_faceNormal] = -0.5;
                  HiPt[a_faceNormal] = -0.5;
                }
              else
                {
                  LoPt[a_faceNormal] =  0.5;
                  HiPt[a_faceNormal] =  0.5;
                }

              if (lohi == 0)
                {
                  LoPt[dom] = -0.5;
                  HiPt[dom] = -0.5;
                }
              else
                {
                  LoPt[dom] =  0.5;
                  HiPt[dom] =  0.5;
                }

              if (LoPtChanged)
                {
                  LoPt[range] -= a_origin[range];
                  LoPt[range] /= a_vectDx[range];
                  LoPt[range] -= (a_iv[range] + 0.5);
                }
              else
                {
                  LoPt[range] = -0.5;
                }

              if (HiPtChanged)
                {
                  HiPt[range] -= a_origin[range];
                  HiPt[range] /= a_vectDx[range];
                  HiPt[range] -= (a_iv[range] + 0.5);
                }
              else
                {
                  HiPt[range] =  0.5;
                }

              CH_assert((!(regular && covered)) && (!(regular && dontKnow)) && (!(dontKnow && covered)));
              CH_assert(regular || covered || (!(LoPtChanged && HiPtChanged)));
              CH_assert(regular || covered || dontKnow);
              CH_assert(regular || covered || LoPtChanged || HiPtChanged);
              bool intersectLo = LoPtChanged;
              edgeMo Edge;
              // range means the coordinate direction that varies over the length of the edge
              Edge.define(LoPt,HiPt,intersectLo, range,covered,regular,dontKnow);
              a_edges[index] = Edge;
            }
        }
    }
  return;
}
void GeometryShop::edgeData2D(edgeMo a_edges[4],
                              bool& a_faceCovered,
                              bool& a_faceRegular,
                              bool& a_faceDontKnow,
                              const Real& a_dx,
                              const RealVect& a_vectDx,
                              const IntVect& a_iv,
                              const ProblemDomain& a_domain,
                              const RealVect& a_origin) const
{
  CH_TIME("GeometryShop::edgeData2D");
  // index counts which edge:xLo=0,xHi=1,yLo=2,yHi=3
  int index = -1;

  a_faceRegular = true;
  a_faceCovered = true;
  a_faceDontKnow = false;

  // domain means the direction normal to the edge
  for (int domain = 0; domain < 2;++domain)
    {
      for (int lohi = 0; lohi < 2;++lohi)
        {
          index += 1;

          // range is the direction along the edge
          int range = 1-domain;

          // Express HiPt. LoPt and MidPt in physical coordinates
          RealVect LoPt;
          bool LoPtChanged = false;
          LoPt[domain] = a_origin[domain] + (a_iv[domain]+lohi)*a_vectDx[domain];
          LoPt[range]  = a_origin[range]  + (a_iv[range])*a_vectDx[range];

          RealVect HiPt;
          bool HiPtChanged = false;
          HiPt[domain] = a_origin[domain] + (a_iv[domain]+lohi)*a_vectDx[domain];
          HiPt[range] =  a_origin[range] + (a_iv[range]+1)*a_vectDx[range];

          RealVect MidPt = HiPt;
          MidPt += LoPt;
          MidPt /= 2.0;

          // decide which type of edge
          bool covered;
          bool regular;
          bool dontKnow;

          // function value
          Real funcHi = m_implicitFunction->value(HiPt);
          Real funcLo = m_implicitFunction->value(LoPt);

          // the sign of signHi and signLo determine edgetype
          Real signHi;
          Real signLo;

          // For level set data negative -> in the fluid
          //                    positive -> out of the fluid
          signHi = -funcHi;
          signLo = -funcLo;

          edgeType(covered,regular,dontKnow,signHi,signLo);

          // Given edgeType, set the edge Hi  and Lo pts
          if (covered)
            {
              a_faceRegular=false;

              LoPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_vectDx[range];
              LoPtChanged = true;

              HiPt[range] = a_origin[range] + (a_iv[range]+0.5)*a_vectDx[range];
              HiPtChanged = true;
            }
          else if (regular)
            {
              a_faceCovered = false;
            }
          else if (dontKnow)
            {
              a_faceRegular  = false;
              a_faceCovered  = false;
              a_faceDontKnow = true;

              // find where the surface intersects the edge
              Real intercept;

              intercept = BrentRootFinder(LoPt, HiPt, range);

              // choose the midpoint for an ill-conditioned problem
              if (intercept<LoPt[range] || intercept>HiPt[range])
                {
                  pout()<<"GeometryShop::edgeData: Ill-conditioned edge data"<<endl;
                  intercept = (LoPt[range]+HiPt[range])/2.0;
                }

              if (funcHi >= 0 && funcLo*funcHi <= 0)
                {
                  HiPt[range] = intercept;
                  HiPtChanged = true;
                }
              else if (funcLo >= 0 && funcLo*funcHi <= 0)
                {
                  LoPt[range] = intercept;
                  LoPtChanged = true;
                }
              else
                {
                  MayDay::Abort("Bogus intersection calculated");
                }
            }

          // express the answer relative to dx and cell-center
          if (lohi == 0)
            {
              LoPt[domain] = -0.5;
              HiPt[domain] = -0.5;
            }
          else
            {
              LoPt[domain] =  0.5;
              HiPt[domain] =  0.5;
            }

          if (LoPtChanged)
            {
              LoPt[range] -= a_origin[range];
              LoPt[range] /= a_vectDx[range];
              LoPt[range] -= (a_iv[range] + 0.5);
            }
          else
            {
              LoPt[range] = -0.5;
            }

          if (HiPtChanged)
            {
              HiPt[range] -= a_origin[range];
              HiPt[range] /= a_vectDx[range];
              HiPt[range] -= (a_iv[range] + 0.5);
            }
          else
            {
              HiPt[range] =  0.5;
            }

          CH_assert(regular || covered || (!(LoPtChanged && HiPtChanged)));
          CH_assert(regular || covered || dontKnow);
          CH_assert((!(regular && covered)) && (!(regular && dontKnow)) && (!(dontKnow && covered)));
          CH_assert(regular || covered || LoPtChanged || HiPtChanged);
          // default is something invalid
          bool intersectLo = LoPtChanged;

          // define this edge
          // Note we have some irregular edges of 0 length
          a_edges[index].define(LoPt,HiPt,intersectLo, range,covered,regular,dontKnow);
        }
    }
}

void GeometryShop::edgeType(bool& a_covered,
                            bool& a_regular,
                            bool& a_dontKnow,
                            Real& a_signHi,
                            Real& a_signLo) const
{
  // if signHi and signLo are both positive
  if (a_signHi > 0.0 && a_signLo > 0.0)
    {
      a_covered  = false;
      a_regular  = true;
      a_dontKnow = false;
    }

  // if signHi and signLo are both negative
  else if (a_signHi <= 0.0 && a_signLo <= 0.0)
    {
      a_covered  = true;
      a_regular  = false;
      a_dontKnow = false;
    }

  // if signHi or signLo both are zero
  else if (a_signHi == 0.0 && a_signLo == 0.0)
    {
      a_covered  = true;
      a_regular  = false;
      a_dontKnow = false;
    }

  // otherwise signLo*signHi <= 0
  // in this case we will look for an intersection point
  else
    {
      a_covered  = false;
      a_regular  = false;
      a_dontKnow = true;
    }

  return;
}

Real GeometryShop::Min(const Real x, const Real y)const
{
  Real retval;
  if (x < y)
    {
      retval = x;
    }
  else
    {
      retval = y;
    }
  return retval;
}

//  The following is an implementation of "Brent's Method"
//    for one-dimensional root finding. Pseudo-code for this
//    algorithm can be found on p. 253 of "Numerical Recipes"
//    ISBN 0-521-30811-9
Real GeometryShop::BrentRootFinder(const RealVect& a_x1,
                                   const RealVect& a_x2,
                                   const int&      a_range) const
{
  const Real tol = PolyGeom::getTolerance();

  //  Max allowed iterations and floating point precision
  const unsigned int  MAXITER = 100;
#if defined(CH_USE_DOUBLE)
  const Real      EPS   = 3.0e-15;
#elif defined(CH_USE_FLOAT)
  const Real      EPS   = 3.0e-7;
#else
#error Unknown Chombo precision
#endif
  unsigned int i;
  RealVect aPt;
  RealVect bPt;
  Real c, fa, fb, fc;
  Real d, e;
  Real tol1, xm;
  Real p, q, r, s;

  aPt = a_x1;
  bPt = a_x2;

  fa = -m_implicitFunction->value(aPt);
  fb = -m_implicitFunction->value(bPt);

  //  Init these to be safe
  c = d = e = 0.0;

  if (fb*fa > 0)
    {
      pout() << "fa " << fa << " fb " << fb <<endl;
      MayDay::Abort("GeometryShop::BrentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
    }

  fc = fb;

  for (i = 0; i < MAXITER; i++)
    {
      if (fb*fc > 0)
        {
          //  Rename a, b, c and adjust bounding interval d
          c = aPt[a_range];
          fc  = fa;
          d = bPt[a_range] - aPt[a_range];
          e = d;
        }

      if (Abs(fc) < Abs(fb))
        {
          aPt[a_range] = bPt[a_range];
          bPt[a_range] = c;
          c = aPt[a_range];
          fa  = fb;
          fb  = fc;
          fc  = fa;
        }

      //  Convergence check
      tol1  = 2.0 * EPS * Abs(bPt[a_range]) + 0.5 * tol;
      xm    = 0.5 * (c - bPt[a_range]);

      if (Abs(xm) <= tol1 || fb == 0.0)
        {
          break;
        }

      if (Abs(e) >= tol1 && Abs(fa) > Abs(fb))
        {
          //  Attempt inverse quadratic interpolation
          s = fb / fa;
          if (aPt[a_range] == c)
            {
              p = 2.0 * xm * s;
              q = 1.0 - s;
            }
          else
            {
              q = fa / fc;
              r = fb / fc;
              p = s * (2.0 * xm * q * (q-r) - (bPt[a_range]-aPt[a_range]) * (r-1.0));
              q = (q-1.0) * (r-1.0) * (s-1.0);
            }

          //  Check whether in bounds
          if (p > 0) q = -q;

          p = Abs(p);

          if (2.0 * p < Min(3.0*xm*q-Abs(tol1*q), Abs(e*q)))
            {
              //  Accept interpolation
              e = d;
              d = p / q;
            }
          else
            {
              //  Interpolation failed, use bisection
              d = xm;
              e = d;
            }
        }
      else
        {
          //  Bounds decreasing too slowly, use bisection
          d = xm;
          e = d;
        }

      //  Move last best guess to a
      aPt[a_range] = bPt[a_range];
      fa  = fb;

      //  Evaluate new trial root
      if (Abs(d) > tol1)
        {
          bPt[a_range] = bPt[a_range] + d;
        }
      else
        {
          if (xm < 0) bPt[a_range] = bPt[a_range] - tol1;
          else        bPt[a_range] = bPt[a_range] + tol1;
        }

      fb = -m_implicitFunction->value(bPt);
    }

  if (i >= MAXITER)
    {
      cerr  << "BrentRootFinder: exceeding maximum iterations: "
            << MAXITER << endl;
    }
  //  //  Keep statistics
  //     statCount++;
  //     statSum += i;
  //     statSum2  += i*i;

  return bPt[a_range];
}

// Coordinates are in the box: [-0.5,0.5] x [-0.5,0.5]
// It is assumed that the spacing in x is uniform
Real GeometryShop::PrismoidalAreaCalc(RealVect& a_xVec,
                                      RealVect& a_yVec) const
{
  // The area of the parabola determined by the three input points
  Real retval;

  // See if parabolic approximation, a*x^2 + b*x + c, stays in bounds
  Real aScale;  // a times 2h^2
  Real bScale;  // b times 2h
  Real cScale;  // c

  aScale = a_yVec[2] - 2*a_yVec[1] + a_yVec[0];
  bScale = a_yVec[2] - a_yVec[0];
  cScale = a_yVec[1];

  // If the parabolic extreme isn't between a_xVec[0] and a_xVec[2]
  // or the parabolic extreme is inside the current cell
  //
  //   Note:  The tests below were designed to avoid any division operations
  //          to avoid numeric problems.  What follows is an explanation of
  //          the tests.
  //
  //   If we put a_xVec[0] at x = -h, a_xVec[1] at x = 0, and a_xVec[2] at
  //   x = h then the parabolic fit, f(x), to the data is:
  //
  //       f(x) = a*x^2 + b*x + c
  //
  //   Where:
  //
  //       a = aScale / (2*h^2)
  //       b = bScale / (2*h)
  //       c = cScale
  //
  //   The location of the parabolic extreme is:
  //
  //       xExtreme = -b / (2*a) = h * -bScale / (2*aScale)
  //
  //   This isn't between a_xVec[0] and a_xVec[2] if |x| > h and this implies:
  //
  //       | h * -bScale / (2*aScale) | > h  =>
  //       | bScale | > 2 * | aScale |
  //
  //   Which is the first test below.
  //
  //   The value at the parabolic extreme is:
  //
  //       f(xExtreme) = (4*a*c - b^2) / (4*a)
  //                   = (8*aScale*cScale - bScale^2) / (8*aScale)
  //
  //   And this is inbounds if it is between -1/2 and 1/2 implying:
  //
  //       | f(xExtreme) | <= 1/2                                =>
  //       | (8*aScale*cScale - bScale^2) / (8*aScale) | <= 1/2  =>
  //       | (8*aScale*cScale - bScale^2) | <= | 4*aScale |
  //
  //   Which is the second test below.
  if ((Abs(bScale) > 2*Abs(aScale)) ||
      (Abs(8*aScale*cScale - bScale*bScale) <= Abs(4*aScale)))
    {
      // Compute the area of the parabola

      // Integrate w.r.t. x, y is a function of x
      // (x_vec[1],y_vec[1]) is the middle point

      // Learn whether x_vec[0] or x_vec[2] is smaller (with origin = cell center)
      Real largeX = a_xVec[2];
      Real smallX = a_xVec[0];

      // Y's that go with the largeX and smallX
      Real largeX_Y = a_yVec[2];
      Real smallX_Y = a_yVec[0];

      // Swap large for small, if necessary
      if (a_xVec[0] > a_xVec[2])
        {
          largeX = a_xVec[0];
          smallX = a_xVec[2];

          largeX_Y = a_yVec[0];
          smallX_Y = a_yVec[2];
        }

      // Find h
      Real h = 0.5 * (largeX - smallX);
      CH_assert (h >= 0.0);

      // Prismoidal rule, which is exact on cubics.
      retval = (h/3.0) * (a_yVec[0] + 4.0 * a_yVec[1] + a_yVec[2]);

      // yVec relative to cell center. Hence add 0.5, but
      // adding 0.5 to each y in the line above adds (h/3) * 3 = h)
      retval += h;

      // Unaccounted for rectangle at the high end of the face.
      retval += (0.5 - largeX) * (largeX_Y + 0.5);

      // Similarly, the unaccounted for rectangle may be at the small end.
      retval += (smallX + 0.5) * (smallX_Y + 0.5);
    }
  else
    {
      // If the extreme of the parabola occurs between a_xVec[0] and a_xVec[2]
      // and the lies outside the current cell then return a negative area so
      // the linear approximation is used
      retval = -1.0;
    }

  return retval;
}

#include "NamespaceFooter.H"
