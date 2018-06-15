#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRGODUNOVFACTORY_H_
#define _EBAMRGODUNOVFACTORY_H_

#include "EBAMRGodunov.H"
#include "AMRLevelFactory.H"
#include "EBPatchGodunovFactory.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBAMRGodunovFactory : public AMRLevelFactory
{
public:
  ///
  AMRLevel* new_amrlevel() const;

  ///
  EBAMRGodunovFactory(const Real&                        a_initialDtMultiplier,
                      const Real&                        a_cfl,
                      const int &                        a_redistRad,
                      const RealVect&                    a_domainLength,
                      const Real&                        a_refineThresh,
                      const int &                        a_tagBufferSize,
                      const int &                        a_verbosity,
                      const bool&                        a_useMassRedist,
                      const bool&                        a_doSmushing,
                      const bool&                        a_doRZCoords,
                      const bool&                        a_hasSourceTerm,
                      const EBPatchGodunovFactory* const a5_patchGodunov,
                      bool                               a_tagAll = false);

  ///
  virtual ~EBAMRGodunovFactory();

protected:
  bool                   m_tagAll;
  Real                   m_initialDtMultiplier;
  Real                   m_cfl;
  int                    m_redistRad;
  RealVect               m_domainLength;
  Real                   m_refineThresh;
  int                    m_tagBufferSize;
  int                    m_verbosity;
  bool                   m_useMassRedist;
  bool                   m_doSmushing;
  bool                   m_doRZCoords;
  bool                   m_hasSourceTerm;
  const EBPatchGodunovFactory* m_patchGodunov;


private:
  EBAMRGodunovFactory()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBAMRGodunovFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRGodunovFactory(const EBAMRGodunovFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
