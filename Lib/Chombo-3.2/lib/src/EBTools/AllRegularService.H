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

#ifndef _ALLREGULARSERVICE_H_
#define _ALLREGULARSERVICE_H_

#include "GeometryService.H"
#include "NamespaceHeader.H"

/// A geometry service that is always regular.
/**
 */
class AllRegularService: public GeometryService
{
public:
  AllRegularService();

  virtual ~AllRegularService();

  ///
  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  virtual bool isRegular(const Box& a_region,
                         const ProblemDomain& a_domain,
                         const RealVect& a_origin,
                         const Real& a_dx) const;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
   virtual bool isCovered(const Box& a_region,
                          const ProblemDomain& a_domain,
                          const RealVect& a_origin,
                          const Real& a_dx) const;

  ///
  /**
   */
  virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                         Vector<IrregNode>&   a_nodes,
                         const Box&           a_validRegion,
                         const Box&           a_ghostRegion,
                         const ProblemDomain& a_domain,
                         const RealVect&      a_origin,
                         const Real&          a_dx) const;
};
#include "NamespaceFooter.H"
#endif
