#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include "EBArith.H"
#include "EBFluxRegister.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBIndexSpace.H"
#include "parstream.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"
#include "CH_Timer.H"
/*******************/
IntVect ivdebugfr(D_DECL(18, 6, 0));
EBFluxRegister::~EBFluxRegister()
{
  //  pout() << "in ebfluxreg destructor" << endl;
}
/*******************/
EBFluxRegister::
EBFluxRegister(const DisjointBoxLayout& a_dblFine,
               const DisjointBoxLayout& a_dblCoar,
               const EBISLayout&        a_ebislFine,
               const EBISLayout&        a_ebislCoar,
               const Box&               a_domainCoar,
               const int&               a_refRat,
               const int&               a_nvar,
               const EBIndexSpace*      ebisPtr,
               bool a_forceNoEBCF)
{
  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_refRat, a_nvar, ebisPtr, a_forceNoEBCF);
}
/*******************/
void
EBFluxRegister::
define(const DisjointBoxLayout& a_gridsFine,
       const DisjointBoxLayout& a_gridsCoar,
       const EBISLayout&        a_ebislFine,
       const EBISLayout&        a_ebislCoar,
       const ProblemDomain&     a_domainCoar,
       const int&               a_refRat,
       const int&               a_nvar,
       const EBIndexSpace*      a_ebisPtr,
       bool a_forceNoEBCF)
{
  ProblemDomain domainFine = refine(a_domainCoar, a_refRat);
  EBLevelGrid eblgCoar(a_gridsCoar, a_ebislCoar, a_domainCoar);
  EBLevelGrid eblgFine(a_gridsFine, a_ebislFine,   domainFine);
  EBFastFR::define(eblgFine, eblgCoar, a_refRat, a_nvar, a_forceNoEBCF);
}
/*******************/
bool
EBFluxRegister::
copyBIFFToEBFF(EBFaceFAB&             a_dst,
               const BaseIFFAB<Real>& a_src,
               const Box            & a_box,
               const EBISBox&         a_ebisBox)
{
  IntVectSet ivs = a_src.getIVS();
  ivs &= a_box;

  bool hasCells = !(ivs.isEmpty());
  if (hasCells)
    {
      a_dst.define(a_ebisBox, a_box, a_src.direction(), a_src.nComp());
      a_dst.setVal(0.);
      for (FaceIterator faceit(ivs, a_ebisBox.getEBGraph(),
                              a_src.direction(),
                              FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          for (int icomp = 0; icomp < a_src.nComp(); icomp++)
            {
              a_dst(faceit(), icomp) = a_src(faceit(), icomp);
            }
        }
    }
  return (hasCells);
}
/*******************/
void
EBFluxRegister::
incrementCoarseIrregular(const BaseIFFAB<Real>& a_coarFlux,
                         const Real&            a_scale,
                         const DataIndex&       a_coarDatInd,
                         const Interval&        a_variables,
                         const int&             a_dir)
{
  if (m_hasEBCF)
    {
      EBFaceFAB facefab;
      bool hasCells = copyBIFFToEBFF(facefab, a_coarFlux, m_eblgCoar.getDBL()[a_coarDatInd], m_eblgCoar.getEBISL()[a_coarDatInd]);
      if (hasCells)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              incrementCoarIrreg(facefab, a_scale, a_coarDatInd, a_variables, a_dir, sit());
            }
        }
    }
}
/*******************/
void
EBFluxRegister::
incrementFineIrregular(const BaseIFFAB<Real>& a_fineFlux,
                       const Real&            a_scale,
                       const DataIndex&       a_fineDatInd,
                       const Interval&        a_variables,
                       const int&             a_dir,
                       const Side::LoHiSide&  a_sd)
{
  if (m_hasEBCF)
    {
      EBFaceFAB facefab;
      bool hasCells = copyBIFFToEBFF(facefab, a_fineFlux, m_eblgFine.getDBL()[a_fineDatInd], m_eblgFine.getEBISL()[a_fineDatInd]);
      if (hasCells)
        {
          incrementFineIrreg(facefab, a_scale, a_fineDatInd, a_variables, a_dir, a_sd);
        }
    }
}
/*******************/
void
EBFluxRegister::
incrementRedistRegister(EBCoarToCoarRedist& a_register,
                        const Interval&     a_variables,
                        const Real&         a_scale)
{
  if (m_hasEBCF)
    {
      LevelData<BaseIVFAB<Real> >& registerMass = a_register.m_regsCoar;
      LayoutData<IntVectSet>&      registerSets = a_register.m_setsCoar;
      //reflux into an empty LevelData<EBCellFAB>
      //Multiply this by (kappa)(1-kappa)
      //add result into register mass
      EBCellFactory ebcfCoar(m_eblgCoar.getEBISL());
      LevelData<EBCellFAB> increment(m_eblgCoar.getDBL(), m_nComp, m_saveCoar.ghostVect(), ebcfCoar);
      EBLevelDataOps::clone (increment, m_saveCoar);
      EBLevelDataOps::setVal(increment, 0.0);

      reflux(increment, a_variables, a_scale, true);

      for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int iindex = index(idir, sit());
                  Vector<IntVectSet> setsCoar = (m_setsCoar[iindex])[dit()];
                  for (int iset = 0; iset < setsCoar.size(); iset++)
                    {
                      const IntVectSet& setCoa  = registerSets[dit()];
                      const IntVectSet& setReg  = setsCoar[iset];
                      IntVectSet set = setCoa;
                      set &= setReg;
                      const EBISBox& ebisBox =m_eblgCoar.getEBISL()[dit()];
                      for (VoFIterator vofit(set, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                        {
                          const VolIndex& vof = vofit();

                          int ibleck = 0;
                          if ((vof.gridIndex() == ivdebugfr) && EBFastFR::s_verbose)
                            {
                              ibleck = 1;
                              pout()    << setprecision(10)
                                        << setiosflags(ios::showpoint)
                                        << setiosflags(ios::scientific);
                              pout() << "incrcotoco:" << endl;

                            }

                          for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                            {
                              Real extraMass = increment[dit()](vofit(), icomp);
                              Real oldMass   =  registerMass[dit()](vofit(), icomp);
                              Real newMass   = oldMass + extraMass;
                              Real diff = extraMass;
                              if (ibleck == 1)
                                {
                                  pout() << "( " << oldMass << ", " << newMass << ", "  << diff << ")";
                                }

                              registerMass[dit()](vofit(), icomp) += extraMass;

                              //set increment to zero in case it gets
                              //hit twice (more than one direction or
                              //whatever
                              increment[dit()](vof, icomp) = 0;

                            } //loop over comps
                          if (ibleck == 1)
                            {
                              ibleck = 0;
                              pout() << endl;
                            }
                        }//loop over vofs in the set
                    }//loop over sets in this coarse box
                } //loop over sides
            } //loop over directions
        } //dataiterator loop
    } //You are using Bonetti's defense against me, uh?
} //I thought it fitting, considering the rocky terrain.

/*******************/
void
EBFluxRegister::
incrementRedistRegister(EBCoarToFineRedist& a_register,
                        const Interval&     a_variables,
                        const Real&         a_scale)
{
  if (m_hasEBCF)
    {
      LevelData<BaseIVFAB<Real> >& registerMass = a_register.m_regsCoar;
      LayoutData<IntVectSet>&      registerSets = a_register.m_setsCoar;
      //reflux into an empty LevelData<EBCellFAB>
      //Multiply this by (kappa)(1-kappa)
      //add result into register mass
      EBCellFactory ebcfCoar(m_eblgCoar.getEBISL());
      LevelData<EBCellFAB> increment(m_eblgCoar.getDBL(), m_nComp, m_saveCoar.ghostVect(), ebcfCoar);
      EBLevelDataOps::clone(increment, m_saveCoar);
      EBLevelDataOps::setVal(increment, 0.0);
      reflux(increment, a_variables, a_scale, true);

      for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int iindex = index(idir, sit());
                  Vector<IntVectSet> setsCoar = (m_setsCoar[iindex])[dit()];
                  for (int iset = 0; iset < setsCoar.size(); iset++)
                    {
                      const IntVectSet& setCoa  = registerSets[dit()];
                      const IntVectSet& setReg  = setsCoar[iset];
                      IntVectSet set = setCoa;
                      set &= setReg;
                      const EBISBox& ebisBox =m_eblgCoar.getEBISL()[dit()];
                      for (VoFIterator vofit(set, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                        {
                          const VolIndex vof = vofit();
                          int ibleck = 0;
                          if ((vof.gridIndex() == ivdebugfr) && EBFastFR::s_verbose)
                            {
                              ibleck = 1;
                              pout()    << setprecision(10)
                                        << setiosflags(ios::showpoint)
                                        << setiosflags(ios::scientific);
                              pout() << "incrcotofi:" << endl;
                            }

                          for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                            {
                              Real extraMass = increment[dit()](vofit(), icomp);

                              if ((ibleck == 1) && icomp == 0)
                                {
                                  //incrredist
                                  Real soluOld = registerMass[dit()](vof, icomp);
                                  Real soluNew = soluOld + extraMass;

                                  Real fluxDif = -extraMass;
                                  pout() << "(" << extraMass << ", " << soluOld << ", " << soluNew << ",  " << fluxDif << ")   " ;
                                }

                              registerMass[dit()](vof, icomp) += extraMass;

                              //set increment to zero in case it gets
                              //hit twice (more than one direction or
                              //whatever
                              increment[dit()](vof, icomp) = 0;

                            }
                          if (ibleck == 1)
                            {
                              pout() << endl;
                              ibleck = 0;
                            }

                        }
                    }
                }
            }
        }
    }
}
/*******************/

#include "NamespaceFooter.H"
