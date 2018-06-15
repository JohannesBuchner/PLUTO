#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _IFINSPECTOR_H_
#define _IFINSPECTOR_H_

#include "BaseIF.H"
#include "parstream.H"
#include "RefCountedPtr.H"
#include "MayDay.H"
#include "Notation.H"
#include "EBCellFAB.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

///
/**
  fill all the values of cell-centered FAB with given origin
  (corner of box) and dx
*/


class IFInspector
{

  void fillValues(EBCellFAB     & a_ebFab,
                  const BaseIF  & a_imIF,
                  const RealVect& a_origin,
                  const RealVect& a_dx);

  void fillValues(LevelData<EBCellFAB> & a_levelFab,
                  const BaseIF         & a_imIF,
                  const RealVect       & a_origin,
                  const RealVect       & a_dx);

  void fillNormal(EBCellFAB   & a_ebFab,
                  const BaseIF& a_imIF,
                  const RealVect& a_origin,
                  const RealVect& a_dx);

  void fillNormal(LevelData<EBCellFAB> & a_levelFab,
                  const BaseIF         & a_imIF,
                  const RealVect       & a_origin,
                  const RealVect       & a_dx);

  void fillGradNormal(EBCellFAB     & a_ebFab,
                      const BaseIF  & a_imIF,
                      const RealVect& a_origin,
                      const RealVect& a_dx);

  void fillGradNormal(LevelData<EBCellFAB>& a_levelFab,
                      const BaseIF        & a_imIF,
                      const RealVect      & a_origin,
                      const RealVect      & a_dx);

};

#include "NamespaceFooter.H"

#endif
