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
#include "EBAMRDataOps.H"
#include "CH_Timer.H"
#include "EBMGAverage.H"
#include "EBArith.H"
#include <iostream>
#include <iomanip>

#include "NamespaceHeader.H"

void
EBAMRDataOps::
getErrorFromCoarseAndFine(Vector< LevelData<EBCellFAB>* >&           a_errorCoar,
                          const Vector< LevelData<EBCellFAB>* >&     a_solnCoar,
                          const Vector< DisjointBoxLayout >&         a_gridsCoar,
                          const Vector< EBISLayout >&                a_ebislCoar,
                          const ProblemDomain&                       a_level0DomainCoar,
                          const Vector< LevelData<EBCellFAB>* >&     a_solnFine,
                          const Vector< DisjointBoxLayout >&         a_gridsFine,
                          const Vector< EBISLayout >&                a_ebislFine,
                          const ProblemDomain&                       a_level0DomainFine,
                          const Vector<int>&                         a_refRat,
                          bool                                       a_kappaAlreadyMultipliedIn)
{

  int nlevels = a_solnFine.size();
  a_errorCoar.resize(nlevels);
  int nref = 2;  //nothing to do with param refinement ratio. this is the refinement between the two solutions
  int nvar = a_solnFine[0]->nComp();
  Interval interv(0, nvar-1);

  ProblemDomain domLevCoar = a_level0DomainCoar;
  const EBIndexSpace* ebis = a_ebislFine[0].getEBIS();
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebislCoar[ilev]);
      int nvar = a_solnCoar[ilev]->nComp();

      a_errorCoar[ilev] = new LevelData<EBCellFAB>(a_gridsCoar[ilev], nvar,   a_solnFine[ilev]->ghostVect(), ebcellfact);

      if (a_kappaAlreadyMultipliedIn)
        {
          //use simple averaging
          IntVect gv = a_solnFine[ilev]->ghostVect();
          EBMGAverage averageOp(a_gridsFine[ilev],
                                a_gridsCoar[ilev],
                                a_ebislFine[ilev],
                                a_ebislCoar[ilev],
                                domLevCoar, nref, nvar, ebis, gv);
          //here make error = Ave(fine)
          averageOp.average(*a_errorCoar[ilev], *a_solnFine[ilev], interv);
        }
      else
        {
          //use kappa weighted averaging
          EBCoarseAverage averageOp(a_gridsFine[ilev],
                                    a_gridsCoar[ilev],
                                    a_ebislFine[ilev],
                                    a_ebislCoar[ilev],
                                    domLevCoar, nref, nvar, ebis);
          //here make error = Ave(fine)
          averageOp.average(*a_errorCoar[ilev], *a_solnFine[ilev], interv);
        }
      //now subtract off coarse so error= Ave(Fine) - coar
      for (DataIterator dit = a_gridsCoar[ilev].dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB&       errorFAB = (*a_errorCoar[ilev])[dit()];
          const EBCellFAB& solnFAB  = (*a_solnCoar[ilev])[dit()];

          errorFAB -= solnFAB;
        }
      domLevCoar.refine(a_refRat[ilev]);
    }
}

void
EBAMRDataOps::
getErrorFromCoarseAndFine(LevelData<EBCellFAB>&           a_errorCoar,
                          const LevelData<EBCellFAB>&     a_solnCoar,
                          const DisjointBoxLayout &       a_gridsCoar,
                          const EBISLayout &              a_ebislCoar,
                          const ProblemDomain&            a_level0DomainCoar,
                          const LevelData<EBCellFAB>&     a_solnFine,
                          const DisjointBoxLayout &       a_gridsFine,
                          const EBISLayout &              a_ebislFine,
                          const ProblemDomain&            a_level0DomainFine,
                          bool a_kappaAlreadyMultipliedIn)
{

  int nref = 2;  //nothing to do with param refinement ratio. this is the refinement between the two solutions
  int nvar = a_solnFine.nComp();
  Interval interv(0, nvar-1);

  ProblemDomain domLevCoar = a_level0DomainCoar;
  const EBIndexSpace* ebis = a_ebislFine.getEBIS();
  EBCellFactory ebcellfact(a_ebislCoar);

  IntVect gv = a_solnFine.ghostVect();
  a_errorCoar.define(a_gridsCoar, nvar,   gv, ebcellfact);

  if (a_kappaAlreadyMultipliedIn)
    {
      //use simple averaging
      EBMGAverage averageOp(a_gridsFine,
                            a_gridsCoar,
                            a_ebislFine,
                            a_ebislCoar,
                            domLevCoar, nref, nvar, ebis, gv);
      //here make error = Ave(fine)
      averageOp.average(a_errorCoar, a_solnFine, interv);
    }
  else
    {
      //use kappa weighted averaging
      EBCoarseAverage averageOp(a_gridsFine,
                                a_gridsCoar,
                                a_ebislFine,
                                a_ebislCoar,
                                domLevCoar, nref, nvar, ebis);
      //here make error = Ave(fine)
      averageOp.average(a_errorCoar, a_solnFine, interv);
    }
  //now subtract off coarse so error= Ave(Fine) - coar
  for (DataIterator dit = a_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&       errorFAB = a_errorCoar[dit()];
      const EBCellFAB& solnFAB  =  a_solnCoar[dit()];
      errorFAB -= solnFAB;
    }
}
int EBAMRDataOps::countVoF(const Vector<DisjointBoxLayout>& a_dbl,
                           const Vector<EBISLayout>&        a_ebisl,
                           const Vector<ProblemDomain>&     a_domain)
{
  //do some counting of boxes and iv's
  pout() << "   Here is an accounting of boxes, ivs, and vofs:" << endl;

  int numLevels = a_domain.size();
  int numBoxAMR     = 0;
  int numRegIVAMR   = 0;
  int numRegInIVAMR = 0;
  int numIrregAMR   = 0;
  int numMultiAMR   = 0;
  int numVofsAMR    = 0;

  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      int numBox     = 0;
      int numRegIV   = 0;
      int numRegInIV = 0;
      int numIrreg   = 0;
      int numMulti   = 0;
      int numVofs    = 0;
      int minBoxes   = a_domain[ilev].domainBox().numPts();//an arbitrary big number
      int maxBoxes   = 0;
      int numDomIV   = a_domain[ilev].domainBox().numPts();

      const DisjointBoxLayout &dbl = a_dbl[ilev];

      //loop over boxes
      for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
        {
          //count this box
          numBox++;

          //count these reg iv's (including covered iv's)
          const Box& dblBox = dbl[dit()];
          numRegIV += dblBox.numPts();

          //count these irreg ivs
          const EBISBox&     ebisBox = a_ebisl[ilev][dit()];
          const IntVectSet& ivsIrreg = ebisBox.getIrregIVS(dblBox);
          const IntVectSet& ivsMulti = ebisBox.getMultiCells(dblBox);
          numIrreg += ivsIrreg.numPts();

          //count all vofs in this dblBox
          IntVectSet ivsBox(dblBox);
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              numVofs++;
              const VolIndex& vof = vofit();
              if (ebisBox.isRegular(vof.gridIndex()))
                {
                  numRegInIV++;
                }
            }

          //count all multi vofs in this dblBox
          for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              numMulti++;
            }
        }

      //figure out the minBoxes and maxBoxes for all procs (measures load balance)
      minBoxes = EBLevelDataOps::parallelMin(numBox);
      maxBoxes = EBLevelDataOps::parallelMax(numBox);

      //make sure we have the right sum over all processors
      Real domainNumBox     = EBLevelDataOps::parallelSum(numBox);
      Real domainNumRegIV   = EBLevelDataOps::parallelSum(numRegIV);
      Real domainNumRegInIV = EBLevelDataOps::parallelSum(numRegInIV);
      Real domainNumIrreg   = EBLevelDataOps::parallelSum(numIrreg);
      Real domainNumMulti   = EBLevelDataOps::parallelSum(numMulti);
      Real domainNumVofs    = EBLevelDataOps::parallelSum(numVofs);

      pout() << "     Stats for Level " << ilev << ":" << endl;
      pout() << "       Domain = " << a_domain[ilev] << endl;
      pout() << "       The number of               boxes = " << setw(12) << left << domainNumBox     << "(processor: "<< numBox     <<")"<<endl;
      pout() << "       The number of        regular iv's = " << setw(12) << left << domainNumRegIV   << "(processor: "<< numRegIV   <<")"<<endl;
      pout() << "       The number of  uncov regular iv's = " << setw(12) << left << domainNumRegInIV << "(processor: "<< numRegInIV <<")"<<endl;
      pout() << "       The number of      irregular iv's = " << setw(12) << left << domainNumIrreg   << "(processor: "<< numIrreg   <<")"<<endl;
      pout() << "       The number of          multi vofs = " << setw(12) << left << domainNumMulti   << "(processor: "<< numMulti   <<")"<<endl;
      pout() << "       The number of                vofs = " << setw(12) << left << domainNumVofs    << "(processor: "<< numVofs    <<")"<<endl;
      pout() << "       The number of        domain cells = " << setw(12) << left << numDomIV         <<endl;
#ifdef CH_MPI
      pout() << "       The min number   boxes for a proc = " << minBoxes  << endl;
      pout() << "       The max number   boxes for a proc = " << maxBoxes  << endl;
#endif
      Real frac = (Real)numRegIV/(Real)(a_domain[ilev].domainBox().numPts());
      pout() << "       (#regular cells)/(#domain cells)  = " << frac << endl;

      //increment the AMR counters
      numBoxAMR     += domainNumBox;
      numRegIVAMR   += domainNumRegIV;
      numRegInIVAMR += domainNumRegInIV;
      numIrregAMR   += domainNumIrreg;
      numMultiAMR   += domainNumMulti;
      numVofsAMR    += domainNumVofs;
    }

  if (numLevels > 1)
    {
      pout() << "     Summary Stats for the " << numLevels << " level AMR Hierarchy:" << endl;
      pout() << "       The number of               boxes = " << numBoxAMR     << endl;
      pout() << "       The number of        regular iv's = " << numRegIVAMR   << endl;
      pout() << "       The number of  uncov regular iv's = " << numRegInIVAMR << endl;
      pout() << "       The number of      irregular iv's = " << numIrregAMR   << endl;
      pout() << "       The number of          multi vofs = " << numMultiAMR   << endl;
      pout() << "       The number of                vofs = " << numVofsAMR    << endl;
      Real amrFrac = (Real)numRegIVAMR/(Real)(a_domain[numLevels-1].domainBox().numPts());
      pout() << "       (#regular cells)/(#finest domain cells) = " << amrFrac << endl;
    }
  return numVofsAMR;
}

void EBAMRDataOps::averageCellToFaces(Vector<LevelData<EBFluxFAB>* >&        a_dataFlux,
                                      const Vector<LevelData<EBCellFAB>* >&  a_dataCell,
                                      const Vector<EBLevelGrid>&             a_eblg,
                                      const Vector<int>&                     a_refRatio,
                                      const int&                             a_comp)
{
  int numLevels = a_dataCell.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::averageCellToFaces(*a_dataFlux[ilev],
                                         *a_dataCell[ilev],
                                         a_eblg[ilev].getDBL(),
                                         a_eblg[ilev].getEBISL(),
                                         a_eblg[ilev].getDomain(),
                                         a_comp);
    }

  averageDown(a_dataFlux,
              a_eblg,
              a_refRatio);

}

void EBAMRDataOps::averageCellToFacesMAC(Vector<LevelData<EBFluxFAB>* >&        a_dataFlux,
                                         const Vector<LevelData<EBCellFAB>* >&  a_dataCell,
                                         const Vector<EBLevelGrid>&             a_eblg,
                                         const Vector<int>&                     a_refRatio)
{
  int numLevels = a_dataCell.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::averageCellToFacesMAC(*a_dataFlux[ilev],
                                            *a_dataCell[ilev],
                                            a_eblg[ilev].getDBL(),
                                            a_eblg[ilev].getEBISL(),
                                            a_eblg[ilev].getDomain());
    }

  averageDown(a_dataFlux,
              a_eblg,
              a_refRatio);

}

void EBAMRDataOps::averageCellToFaces(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_dataFlux,
                                      const Vector<LevelData<EBCellFAB>* >&          a_dataCell,
                                      const Vector<EBLevelGrid>&                     a_eblg,
                                      const Vector<int>&                             a_refRatio,
                                      const int&                                     a_comp)
{
  int numLevels = a_dataCell.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::averageCellToFaces(*a_dataFlux[ilev],
                                         *a_dataCell[ilev],
                                         a_eblg[ilev].getDBL(),
                                         a_eblg[ilev].getEBISL(),
                                         a_eblg[ilev].getDomain(),
                                         a_comp);
    }

  averageDown(a_dataFlux,
              a_eblg,
              a_refRatio);

}

void EBAMRDataOps::quadCFInterpAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<EBLevelGrid>&       a_eblg,
                                   const Vector<int>&               a_refRatio)
{
  int numLevels = a_data.size();
  for (int ilev = 1; ilev < numLevels; ilev++)
    {
      quadCFInterpOne(a_data,
                      a_eblg,
                      a_refRatio,
                      ilev);
    }
}
void EBAMRDataOps::quadCFInterpAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<DisjointBoxLayout>& a_dbl,
                                   const Vector<EBISLayout>&        a_ebisl,
                                   const Vector<ProblemDomain>&     a_domain,
                                   const Vector<int>&               a_refRatio)
{
  int numLevels = a_data.size();
  for (int ilev = 1; ilev < numLevels; ilev++)
    {
      quadCFInterpOne(a_data,
                      a_dbl,
                      a_ebisl,
                      a_domain,
                      a_refRatio,
                      ilev);
    }
}

void EBAMRDataOps::quadCFInterpOne(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<EBLevelGrid>&       a_eblg,
                                   const Vector<int>&               a_refRatio,
                                   const int&                       a_fineLevel)
{
  int coarseLevel = a_fineLevel-1;
  CH_assert(coarseLevel > -1);
  int ncomp = a_data[0]->nComp();
  EBQuadCFInterp quadCFIWithFine(a_eblg[a_fineLevel].getDBL(),
                                 a_eblg[coarseLevel].getDBL(),
                                 a_eblg[a_fineLevel].getEBISL(),
                                 a_eblg[coarseLevel].getEBISL(),
                                 a_eblg[coarseLevel].getDomain(),
                                 a_refRatio[coarseLevel],
                                 ncomp, *a_eblg[a_fineLevel].getCFIVS());

  Interval  intval(0,ncomp-1);
  quadCFIWithFine.interpolate(*a_data[a_fineLevel],
                              *a_data[coarseLevel],
                              intval);
  a_data[a_fineLevel]->exchange(intval);
}
void EBAMRDataOps::quadCFInterpOne(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<DisjointBoxLayout>& a_dbl,
                                   const Vector<EBISLayout>&        a_ebisl,
                                   const Vector<ProblemDomain>&     a_domain,
                                   const Vector<int>&               a_refRatio,
                                   const int&                       a_fineLevel)
{
  int coarseLevel = a_fineLevel-1;
  CH_assert(coarseLevel > -1);
  int ncomp = a_data[0]->nComp();
  LayoutData<IntVectSet> cfivs;
  EBArith::defineCFIVS(cfivs, a_dbl[a_fineLevel], a_domain[a_fineLevel]);
  EBQuadCFInterp quadCFIWithFine(a_dbl[a_fineLevel],
                                 a_dbl[coarseLevel],
                                 a_ebisl[a_fineLevel],
                                 a_ebisl[coarseLevel],
                                 a_domain[coarseLevel],
                                 a_refRatio[coarseLevel],
                                 ncomp, cfivs);

  Interval  intval(0,ncomp-1);
  quadCFIWithFine.interpolate(*a_data[a_fineLevel],
                              *a_data[coarseLevel],
                              intval);
  a_data[a_fineLevel]->exchange(intval);
}

void EBAMRDataOps::quadCFInterpOne(LevelData<EBCellFAB>&       a_dataFine,
                                   const LevelData<EBCellFAB>& a_dataCoar,
                                   const DisjointBoxLayout&    a_dblFine,
                                   const DisjointBoxLayout&    a_dblCoar,
                                   const EBISLayout&           a_ebislFine,
                                   const EBISLayout&           a_ebislCoar,
                                   const ProblemDomain&        a_domainCoar,
                                   const int&                  a_refRatio)
{
  int ncomp = a_dataCoar.nComp();
  LayoutData<IntVectSet> cfivs;
  ProblemDomain domFine = refine(a_domainCoar, a_refRatio);
  EBArith::defineCFIVS(cfivs, a_dblFine, domFine);
  EBQuadCFInterp quadCFIWithFine(a_dblFine,
                                 a_dblCoar,
                                 a_ebislFine,
                                 a_ebislCoar,
                                 a_domainCoar,
                                 a_refRatio,
                                 ncomp, cfivs);

  Interval  intval(0,ncomp-1);
  quadCFIWithFine.interpolate(a_dataFine,
                              a_dataCoar,
                              intval);
  a_dataFine.exchange(intval);
}

void EBAMRDataOps::pwlFillPatchAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<EBLevelGrid>&       a_eblg,
                                   const Vector<int>&               a_refRatio)
{
  CH_TIME("EBAMRDataOps::pwlFillPatchAll(lg)");
  int numLevels = a_data.size();
  for (int ilev = 1; ilev < numLevels; ilev++)
    {
      pwlFillPatchOne(*a_data[ilev],
                      *a_data[ilev-1],
                      a_eblg[ilev  ].getDBL(),
                      a_eblg[ilev-1].getDBL(),
                      a_eblg[ilev  ].getEBISL(),
                      a_eblg[ilev-1].getEBISL(),
                      a_eblg[ilev-1].getDomain(),
                      a_refRatio[ilev-1]);
    }
}
void EBAMRDataOps::pwlFillPatchAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<DisjointBoxLayout>& a_dbl,
                                   const Vector<EBISLayout>&        a_ebisl,
                                   const Vector<ProblemDomain>&     a_domain,
                                   const Vector<int>&               a_refRatio)
{
  CH_TIME("EBAMRDataOps::pwlFillPatchAll");
  int numLevels = a_data.size();
  for (int ilev = 1; ilev < numLevels; ilev++)
    {
      pwlFillPatchOne(*a_data[ilev],
                      *a_data[ilev-1],
                      a_dbl[ilev],
                      a_dbl[ilev-1],
                      a_ebisl[ilev],
                      a_ebisl[ilev-1],
                      a_domain[ilev-1],
                      a_refRatio[ilev-1]);
    }
}

void EBAMRDataOps::pwlFillPatchOne(LevelData<EBCellFAB>&       a_dataFine,
                                   LevelData<EBCellFAB>&       a_dataCoar,
                                   const DisjointBoxLayout&    a_dblFine,
                                   const DisjointBoxLayout&    a_dblCoar,
                                   const EBISLayout&           a_ebislFine,
                                   const EBISLayout&           a_ebislCoar,
                                   const ProblemDomain&        a_domainCoar,
                                   const int&                  a_refRatioCoar)
{
  CH_TIME("EBAMRDataOps::pwlFillPatchOne");
  int nComp = a_dataCoar.nComp();
  IntVect ghostVect = a_dataCoar.ghostVect();

  Interval  intval(0,nComp-1);
  EBPWLFillPatch patcher(a_dblFine,
                         a_dblCoar,
                         a_ebislFine,
                         a_ebislCoar,
                         a_domainCoar,
                         a_refRatioCoar,
                         nComp,
                         ghostVect[0]);

  a_dataCoar.exchange(intval);

  Real coarTimeOld = 0.0;
  Real coarTimeNew = 1.0;
  Real fineTime    = 0.0;
  patcher.interpolate(a_dataFine,
                      a_dataCoar,
                      a_dataCoar,
                      coarTimeOld,
                      coarTimeNew,
                      fineTime,
                      intval);

  a_dataFine.exchange(intval);
}

void EBAMRDataOps::exchangeAll(Vector<LevelData<EBCellFAB>* >& a_data)
{
  CH_TIME("EBAMRDataOps::exchangeAll(cell)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeAll(*a_data[ilev]);
    }
}
void EBAMRDataOps::exchangeCorners(Vector<LevelData<EBCellFAB>* >& a_data,
                                   const ProblemDomain&            a_domain)
{
  CH_TIME("EBAMRDataOps::exchangeCorners(cell)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeCorners(*a_data[ilev],a_domain);
    }
}

void EBAMRDataOps::exchangeComp(Vector<LevelData<EBCellFAB>* >& a_data,
                                const int&                      a_comp)
{
  CH_TIME("EBAMRDataOps::exchangeComp(cell,comp)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeComp(*a_data[ilev],
                                   a_comp);
    }
}

void EBAMRDataOps::exchangeAll(Vector<LevelData<EBFluxFAB>* >& a_data)
{
  CH_TIME("EBAMRDataOps::exchangeAll(flux)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeAll(*a_data[ilev]);
    }
}

void EBAMRDataOps::exchangeComp(Vector<LevelData<EBFluxFAB>* >& a_data,
                                const int&                      a_comp)
{
  CH_TIME("EBAMRDataOps::exchangeComp(flux,comp)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeComp(*a_data[ilev],
                                   a_comp);
    }
}
void EBAMRDataOps::exchangeAll(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_data)
{
  CH_TIME("EBAMRDataOps::exchangeAll(cellRef)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeAll(*a_data[ilev]);
    }
}

void EBAMRDataOps::exchangeComp(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_data,
                                const int&                      a_comp)
{
  CH_TIME("EBAMRDataOps::exchangeComp(cellRef,comp)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeComp(*a_data[ilev],
                                   a_comp);
    }
}

void EBAMRDataOps::exchangeAll(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_data)
{
  CH_TIME("EBAMRDataOps::exchangeAll(fluxRef)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeAll(*a_data[ilev]);
    }
}

void EBAMRDataOps::exchangeComp(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_data,
                                const int&                      a_comp)
{
  CH_TIME("EBAMRDataOps::exchangeComp(fluxRef,comp)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::exchangeComp(*a_data[ilev],
                                   a_comp);
    }
}

void EBAMRDataOps::coarsenDown(Vector<LevelData<EBCellFAB>* >&  a_data,
                               const Vector<EBLevelGrid>&       a_eblg,
                               const Vector<int>&               a_refRatio)

{
  // do coarsen down here
  int numLevels = a_data.size();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      //define as one-variable object and call for both velocity and scalars
      EBCoarsen ebcoarsen(a_eblg[ilev  ],
                          a_eblg[ilev-1],
                          a_refRatio[ilev-1],
                          ncomp);

      Interval interval(0, ncomp-1);
      ebcoarsen.coarsenFine(*a_data[ilev-1],
                            *a_data[ilev  ],
                            interval);
    } // end loop over levels
  exchangeAll(a_data);
}
void EBAMRDataOps::averageDown(Vector<LevelData<EBCellFAB>* >&  a_data,
                               const Vector<EBISLayout>&        a_ebisl,
                               const Vector<DisjointBoxLayout>& a_dbl,
                               const Vector<ProblemDomain>&     a_domain,
                               const Vector<int>&               a_refRatio)

{
  CH_TIME("EBAMRDataOps::averageDown(data,ebisl,dbl,domain,refrat)");
  // do average down here
  int numLevels = a_data.size();
  const EBIndexSpace* ebis = a_ebisl[0].getEBIS();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      EBCoarseAverage ebaverage(a_dbl[ilev  ],
                                a_dbl[ilev-1],
                                a_ebisl[ilev  ],
                                a_ebisl[ilev-1],
                                a_domain[ilev-1],
                                a_refRatio[ilev-1],
                                ncomp,
                                ebis);

      Interval interval(0, ncomp-1);
      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  exchangeAll(a_data);
}
void EBAMRDataOps::averageDown(Vector<LevelData<EBCellFAB>* >&  a_data,
                               const Vector<EBLevelGrid>&       a_eblg,
                               const Vector<int>&               a_refRatio)

{
  CH_TIME("EBAMRDataOps::averageDown(data,eblg,refrat)");
  // do average down here
  int numLevels = a_data.size();
  const EBIndexSpace* ebis = a_eblg[0].getEBIS();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      EBCoarseAverage ebaverage(a_eblg[ilev  ].getDBL(),
                                a_eblg[ilev-1].getDBL(),
                                a_eblg[ilev  ].getEBISL(),
                                a_eblg[ilev-1].getEBISL(),
                                a_eblg[ilev-1].getDomain(),
                                a_refRatio[ilev-1],
                                ncomp,
                                ebis);

      Interval interval(0, ncomp-1);
      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  exchangeAll(a_data);
}
void EBAMRDataOps::averageDown(LevelData<EBCellFAB>&       a_dataCoar,
                               const LevelData<EBCellFAB>& a_dataFine,
                               const EBISLayout&           a_ebislCoar,
                               const EBISLayout&           a_ebislFine,
                               const DisjointBoxLayout&    a_dblCoar,
                               const DisjointBoxLayout&    a_dblFine,
                               const ProblemDomain&        a_domainCoar,
                               const int&                  a_refRatio)

{
  CH_TIME("EBAMRDataOps::averageDown(dataC,dataF,ebislC,ebislF,dblC,dblF,domC,refrat)");
  // do average down here
  int ncomp = a_dataFine.nComp();
  const EBIndexSpace* ebis = a_ebislFine.getEBIS();
  EBCoarseAverage ebaverage(a_dblFine,
                            a_dblCoar,
                            a_ebislFine,
                            a_ebislCoar,
                            a_domainCoar,
                            a_refRatio,
                            ncomp,
                            ebis);

  Interval interval(0, ncomp-1);
  ebaverage.average(a_dataCoar,
                    a_dataFine,
                    interval);
  a_dataCoar.exchange();
}

void EBAMRDataOps::averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data,
                               const Vector<EBISLayout>&        a_ebisl,
                               const Vector<DisjointBoxLayout>& a_dbl,
                               const Vector<ProblemDomain>&     a_domain,
                               const Vector<int>&               a_refRatio)
{
  CH_TIME("EBAMRDataOps::averageDown(fluxData,ebisl,dbl,dom,refrat");
  // do average down here
  int numLevels = a_data.size();
  const EBIndexSpace* ebis = a_ebisl[0].getEBIS();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      EBCoarseAverage ebaverage(a_dbl[ilev  ],
                                a_dbl[ilev-1],
                                a_ebisl[ilev  ],
                                a_ebisl[ilev-1],
                                a_domain[ilev-1],
                                a_refRatio[ilev-1],
                                ncomp,
                                ebis);

      Interval interval(0, ncomp-1);
      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  exchangeAll(a_data);
}
void EBAMRDataOps::averageDown(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_data,
                               const Vector<EBISLayout>&                      a_ebisl,
                               const Vector<DisjointBoxLayout>&               a_dbl,
                               const Vector<ProblemDomain>&                   a_domain,
                               const Vector<int>&                             a_refRatio)
{
  CH_TIME("EBAMRDataOps::averageDown(dataRef,ebisl,dbl,dom,refrat");
  // do average down here
  int numLevels = a_data.size();
  const EBIndexSpace* ebis = a_ebisl[0].getEBIS();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      EBCoarseAverage ebaverage(a_dbl[ilev  ],
                                a_dbl[ilev-1],
                                a_ebisl[ilev  ],
                                a_ebisl[ilev-1],
                                a_domain[ilev-1],
                                a_refRatio[ilev-1],
                                ncomp,
                                ebis);

      Interval interval(0, ncomp-1);
      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  exchangeAll(a_data);
}

void EBAMRDataOps::averageDown(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_data,
                               const Vector<EBISLayout>&                      a_ebisl,
                               const Vector<DisjointBoxLayout>&               a_dbl,
                               const Vector<ProblemDomain>&                   a_domain,
                               const Vector<int>&                             a_refRatio)
{
  CH_TIME("EBAMRDataOps::averageDown(fluxDataRef,ebisl,dbl,dom,refrat");
  // do average down here
  int numLevels = a_data.size();
  const EBIndexSpace* ebis = a_ebisl[0].getEBIS();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      EBCoarseAverage ebaverage(a_dbl[ilev  ],
                                a_dbl[ilev-1],
                                a_ebisl[ilev  ],
                                a_ebisl[ilev-1],
                                a_domain[ilev-1],
                                a_refRatio[ilev-1],
                                ncomp,
                                ebis);

      Interval interval(0, ncomp-1);
      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  exchangeAll(a_data);
}

void EBAMRDataOps::averageDown(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_data,
                               const Vector<EBLevelGrid>&                     a_eblg,
                               const Vector<int>&                             a_refRatio)
{
  CH_TIME("EBAMRDataOps::averageDown(fluxDataRef,eblg,refrat");
  // do average down here
  int numLevels = a_data.size();
  const EBIndexSpace* ebis = a_eblg[0].getEBIS();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      EBCoarseAverage ebaverage(a_eblg[ilev  ].getDBL(),
                                a_eblg[ilev-1].getDBL(),
                                a_eblg[ilev  ].getEBISL(),
                                a_eblg[ilev-1].getEBISL(),
                                a_eblg[ilev-1].getDomain(),
                                a_refRatio[ilev-1],
                                ncomp,
                                ebis);

      Interval interval(0, ncomp-1);
      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  exchangeAll(a_data);
}

void EBAMRDataOps::averageDown(Vector<LevelData<EBFluxFAB>* > & a_data,
                               const Vector<EBLevelGrid>&       a_eblg,
                               const Vector<int>&               a_refRatio)
{
  CH_TIME("EBAMRDataOps::averageDown(fluxData,eblg,refrat)");
  // do average down here
  int numLevels = a_data.size();
  const EBIndexSpace* ebis = a_eblg[0].getEBIS();
  for (int ilev = (numLevels - 1); ilev>0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      EBCoarseAverage ebaverage(a_eblg[ilev  ].getDBL(),
                                a_eblg[ilev-1].getDBL(),
                                a_eblg[ilev  ].getEBISL(),
                                a_eblg[ilev-1].getEBISL(),
                                a_eblg[ilev-1].getDomain(),
                                a_refRatio[ilev-1],
                                ncomp,
                                ebis);

      Interval interval(0, ncomp-1);
      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  exchangeAll(a_data);
}


void EBAMRDataOps::setCoveredAMRVal(Vector<LevelData<EBCellFAB>* >& a_data,
                                    const Vector<EBISLayout>&       a_ebisl,
                                    const Vector<int>&              a_refRat,
                                    const Real&                     a_value)
{
  CH_TIME("EBAMRDataOps::setCoveredAMRVal");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels-1; ilev++)
    {
      int coarLev = ilev;
      int fineLev = ilev+1;

      int ncomp = (*a_data[ilev]).nComp();
      Interval interv(0, ncomp-1);
      const DisjointBoxLayout& coarGrids = (*a_data[coarLev]).disjointBoxLayout();
      const DisjointBoxLayout& fineGrids = (*a_data[fineLev]).disjointBoxLayout();

      //zero out under finer grids
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB&     coarFAB = (*a_data[coarLev])[dit()];
          const EBISBox& ebisBox = a_ebisl[coarLev][dit()];

          LayoutIterator litFine = fineGrids.layoutIterator();
          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlayBox = coarFAB.box();
              Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat[coarLev]);

              overlayBox &= coarsenedGrid;
              if (!overlayBox.isEmpty())
                {
                  IntVectSet ivsZero(overlayBox);

                  for (VoFIterator vofit(ivsZero, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                    {
                      for (int ivar =0; ivar < ncomp; ivar++)
                        {
                          coarFAB(vofit(), ivar) = a_value;
                        }
                    }
                }
            }
        }
    }
}

void EBAMRDataOps::setCoveredAMRVal(Vector<LevelData<EBCellFAB>* >& a_data,
                                    const Vector<EBLevelGrid>&      a_eblg,
                                    const Vector<int>&              a_refRat,
                                    const Real&                     a_value)
{
  CH_TIME("EBAMRDataOps::setCoveredAMRVal(eblg)");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels-1; ilev++)
    {
      int coarLev = ilev;
      int fineLev = ilev+1;

      int ncomp = (*a_data[ilev]).nComp();
      Interval interv(0, ncomp-1);
      const DisjointBoxLayout& coarGrids = (*a_data[coarLev]).disjointBoxLayout();
      const DisjointBoxLayout& fineGrids = (*a_data[fineLev]).disjointBoxLayout();

      //zero out under finer grids
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB&     coarFAB = (*a_data[coarLev])[dit()];
          const EBISBox& ebisBox = a_eblg[coarLev].getEBISL()[dit()];

          LayoutIterator litFine = fineGrids.layoutIterator();
          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlayBox = coarFAB.box();
              Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat[coarLev]);

              overlayBox &= coarsenedGrid;
              if (!overlayBox.isEmpty())
                {
                  IntVectSet ivsZero(overlayBox);

                  for (VoFIterator vofit(ivsZero, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                    {
                      for (int ivar =0; ivar < ncomp; ivar++)
                        {
                          coarFAB(vofit(), ivar) = a_value;
                        }
                    }
                }
            }
        }
    }
}



void EBAMRDataOps::setCoveredVal(Vector<LevelData<EBCellFAB>* >&a_data,
                                 const Real& a_value)
{
  CH_TIME("EBAMRDataOps::setCoveredVal");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::setCoveredVal(*a_data[ilev],
                                    a_value);
    }

}

void EBAMRDataOps::defineAMRData(Vector<LevelData<EBCellFAB>* >&  a_amrData,
                                 const Vector<EBLevelGrid>&       a_eblg,
                                 const IntVect&                   a_ghosts,
                                 const int&                       a_nComp,
                                 const int&                       a_numLevels)
{
  a_amrData.resize(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      EBCellFactory ebcellfact(a_eblg[ilev].getEBISL());
      a_amrData[ilev] = new LevelData<EBCellFAB>(a_eblg[ilev].getDBL(),
                                                 a_nComp,
                                                 a_ghosts,
                                                 ebcellfact);
    }
}

void EBAMRDataOps::deleteAMRData(Vector<LevelData<EBCellFAB>* >&  a_amrData)
{
  for (int ilev = 0; ilev < a_amrData.size(); ilev++)
    {
      delete a_amrData[ilev];
    }
}
void EBAMRDataOps::defineAMRData(Vector<LevelData<EBFluxFAB>* >&  a_amrData,
                                 const Vector<EBLevelGrid>&       a_eblg,
                                 const IntVect&                   a_ghosts,
                                 const int&                       a_nComp,
                                 const int&                       a_numLevels)
{
  a_amrData.resize(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      EBFluxFactory ebfluxfact(a_eblg[ilev].getEBISL());
      a_amrData[ilev] = new LevelData<EBFluxFAB>(a_eblg[ilev].getDBL(),
                                                 a_nComp,
                                                 a_ghosts,
                                                 ebfluxfact);
    }
}
void EBAMRDataOps::defineAMRData(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_amrData,
                                 const Vector<EBLevelGrid>&                     a_eblg,
                                 const IntVect&                                 a_ghosts,
                                 const int&                                     a_nComp,
                                 const int&                                     a_numLevels)
{
  a_amrData.resize(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      EBCellFactory ebcellfact(a_eblg[ilev].getEBISL());
      a_amrData[ilev] = RefCountedPtr<LevelData<EBCellFAB> >(
        new LevelData<EBCellFAB>(a_eblg[ilev].getDBL(),
                                 a_nComp,
                                 a_ghosts,
                                 ebcellfact) );
    }
}
void EBAMRDataOps::defineAMRData(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_amrData,
                                 const Vector<EBLevelGrid>&                     a_eblg,
                                 const IntVect&                                 a_ghosts,
                                 const int&                                     a_nComp,
                                 const int&                                     a_numLevels)
{
  a_amrData.resize(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      EBFluxFactory ebfluxfact(a_eblg[ilev].getEBISL());
      a_amrData[ilev] = RefCountedPtr<LevelData<EBFluxFAB> >(
            new LevelData<EBFluxFAB>(a_eblg[ilev].getDBL(),
                                     a_nComp,
                                     a_ghosts,
                                     ebfluxfact) );
    }
}

void EBAMRDataOps::setToZero(Vector<LevelData<EBCellFAB>* >& a_result)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::setToZero(*a_result[ilev]);
    }
}
void EBAMRDataOps::setVal(Vector<LevelData<EBCellFAB>* >& a_result,
                          const Real&                     a_value)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::setVal(*a_result[ilev],
                             a_value);
    }
}
void EBAMRDataOps::setVal(Vector<LevelData<EBCellFAB>* >& a_result,
                          const Real&                     a_value,
                          const int&                      a_comp)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::setVal(*a_result[ilev],
                             a_value,
                             a_comp);
    }
}
void EBAMRDataOps::setToZero(Vector<LevelData<EBFluxFAB>* >& a_result)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::setToZero(*a_result[ilev]);
    }
}
void EBAMRDataOps::setVal(Vector<LevelData<EBFluxFAB>* >& a_result,
                          const Real&                     a_value)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::setVal(*a_result[ilev],
                             a_value);
    }
}

void EBAMRDataOps::axby(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                        const Vector<LevelData<EBCellFAB>* >& a_x,
                        const Vector<LevelData<EBCellFAB>* >& a_y,
                        const Real& a_a,
                        const Real& a_b)
{
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::axby(*a_lhs[ilev],
                           *a_x[ilev],
                           *a_y[ilev],
                           a_a,
                           a_b);
    }
}
void EBAMRDataOps::axby(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                        const Vector<LevelData<EBCellFAB>* >& a_x,
                        const Vector<LevelData<EBCellFAB>* >& a_y,
                        const Real& a_a,
                        const Real& a_b,
                        const int&  a_lhsComp,
                        const int&  a_xComp,
                        const int&  a_yComp)
{
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::axby(*a_lhs[ilev],
                           *a_x[ilev],
                           *a_y[ilev],
                           a_a,
                           a_b,
                           a_lhsComp,
                           a_xComp,
                           a_yComp);
    }
}
void EBAMRDataOps::assign(Vector<LevelData<EBCellFAB>* >&       a_to,
                          const Vector<LevelData<EBCellFAB>* >& a_from,
                          const Interval&                       a_toInterval,
                          const Interval&                       a_fromInterval)
{
  int numLevels = a_to.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::assign(*a_to[ilev],
                             *a_from[ilev],
                             a_toInterval,
                             a_fromInterval);
    }
}
void EBAMRDataOps::assign(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&       a_to,
                          const Vector<LevelData<EBCellFAB>* >&                a_from,
                          const Interval&                                      a_toInterval,
                          const Interval&                                      a_fromInterval)
{
  CH_TIME("EBAMRDataOps::assign(cell_to,cell_from,int_to,int_from)");
  int numLevels = a_to.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::assign(*a_to[ilev],
                             *a_from[ilev],
                             a_toInterval,
                             a_fromInterval);
    }
}
void EBAMRDataOps::assign(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                          const Vector<LevelData<EBCellFAB>* >& a_rhs)
{
  CH_TIME("EBAMRDataOps::assign(cell_lhs,cell_rhs)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::assign(*a_lhs[ilev],
                             *a_rhs[ilev]);
    }
}
void EBAMRDataOps::assign(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&       a_lhs,
                          const Vector<LevelData<EBCellFAB>* >&                a_rhs)
{
  CH_TIME("EBAMRDataOps::assign(cellRef_lhs,cellRef_rhs)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::assign(*a_lhs[ilev],
                             *a_rhs[ilev]);
    }
}
void EBAMRDataOps::assign(Vector<LevelData<EBFluxFAB>* >&       a_lhs,
                          const Vector<LevelData<EBFluxFAB>* >& a_rhs)
{
  CH_TIME("EBAMRDataOps::assign(flux_lhs,flux_rhs)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::assign(*a_lhs[ilev],
                             *a_rhs[ilev]);
    }
}
void EBAMRDataOps::assign(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_lhs,
                          const Vector<LevelData<EBFluxFAB>* >&          a_rhs)
{
  CH_TIME("EBAMRDataOps::assign(fluxRef_lhs,fluxRef_rhs)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::assign(*a_lhs[ilev],
                             *a_rhs[ilev]);
    }
}
void EBAMRDataOps::incr(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                        const Vector<LevelData<EBCellFAB>* >& a_rhs,
                        const Real& a_scale)
{
  CH_TIME("EBAMRDataOps::incr(lhs,rhs,scale)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::incr(*a_lhs[ilev],
                           *a_rhs[ilev],
                           a_scale);
    }
}

void EBAMRDataOps::incr(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                        const Real& a_scale)
{
  CH_TIME("EBAMRDataOps::incr(lhs,scale)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::incr(*a_lhs[ilev],
                           a_scale);
    }
}

void EBAMRDataOps::scale(Vector<LevelData<EBFluxFAB>* >& a_result,
                         const Real&                     a_value)
{
  CH_TIME("EBAMRDataOps::scale(flux_result,scale)");
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::scale(*a_result[ilev],
                            a_value);
    }
}

void EBAMRDataOps::scale(Vector<LevelData<EBCellFAB>* >& a_result,
                         const Real&                     a_value)
{
  CH_TIME("EBAMRDataOps::scale(cell_result,scale)");
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::scale(*a_result[ilev],
                            a_value);
    }
}

void EBAMRDataOps::scale(Vector<LevelData<EBCellFAB>* >& a_result,
                         const Real&                     a_value,
                         const int&                      a_comp)
{
  CH_TIME("EBAMRDataOps::scale(flux_result,scale,comp)");
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::scale(*a_result[ilev],
                            a_value,
                            a_comp);
    }
}

void EBAMRDataOps::sum(Vector<LevelData<EBCellFAB>* >&       a_result,
                       const Vector<LevelData<EBCellFAB>* >& a_in1,
                       const Vector<LevelData<EBCellFAB>* >& a_in2)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::sum(*a_result[ilev],
                          *a_in1[ilev],
                          *a_in2[ilev]);
    }
}

void EBAMRDataOps::addConstant(Vector<LevelData<EBCellFAB>* >& a_data,
                               const Real&           a_constant)
{
  for (int ilev = 0; ilev < a_data.size(); ilev++)
    {
      EBLevelDataOps::addConstant(*a_data[ilev], a_constant);
    }
}

void EBAMRDataOps::product(Vector<LevelData<EBCellFAB>* >&       a_result,
                           const Vector<LevelData<EBCellFAB>* >& a_in1,
                           const Vector<LevelData<EBCellFAB>* >& a_in2)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::product(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev]);
    }
}
void EBAMRDataOps::product(Vector<LevelData<EBCellFAB>* >&       a_result,
                           const Vector<LevelData<EBCellFAB>* >& a_in1,
                           const Vector<LevelData<EBCellFAB>* >& a_in2,
                           const int&                            a_rComp,
                           const int&                            a_1Comp,
                           const int&                            a_2Comp)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::product(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev],
                              a_rComp,
                              a_1Comp,
                              a_2Comp);
    }
}

void EBAMRDataOps::divideVectorByScalar(Vector<LevelData<EBCellFAB>* >&       a_vectorOut,
                                        const Vector<LevelData<EBCellFAB>* >& a_vectorIn,
                                        const Vector<LevelData<EBCellFAB>* >& a_scalar)
{
  //a_vector /= a_scalar
  int numLevels = a_vectorOut.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::divideVectorByScalar(*a_vectorOut[ilev],
                                           *a_vectorIn[ ilev],
                                           *a_scalar[   ilev]);
    }
}

void EBAMRDataOps::divide(Vector<LevelData<EBCellFAB>* >&       a_result,
                          const Vector<LevelData<EBCellFAB>* >& a_in1,
                          const Vector<LevelData<EBCellFAB>* >& a_in2)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::divide(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev]);
    }
}
void EBAMRDataOps::divide(Vector<LevelData<EBCellFAB>* >&       a_result,
                          const Vector<LevelData<EBCellFAB>* >& a_in1,
                          const Vector<LevelData<EBCellFAB>* >& a_in2,
                          const int&                            a_rComp,
                          const int&                            a_1Comp,
                          const int&                            a_2Comp)

{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::divide(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev],
                             a_rComp,
                             a_1Comp,
                             a_2Comp);
    }
}
void EBAMRDataOps::product(Vector<LevelData<EBFluxFAB>* >&       a_result,
                           const Vector<LevelData<EBFluxFAB>* >& a_in1,
                           const Vector<LevelData<EBFluxFAB>* >& a_in2)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::product(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev]);
    }
}

void EBAMRDataOps::kappaWeight(Vector<LevelData<EBCellFAB>* >& a_data)
{
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::kappaWeight(*a_data[ilev]);
    }
}

void EBAMRDataOps::kappaScale(Vector<LevelData<EBCellFAB>* >& a_data,
                              const Real&                     a_value)
{
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBLevelDataOps::kappaScale(*a_data[ilev],
                                 a_value);
    }
}


Real EBAMRDataOps::kappaNorm(Real&                                 a_volume,
                             const Vector<LevelData<EBCellFAB>* >& a_data,
                             int                                   a_which,
                             const Vector<ProblemDomain>&          a_domain,
                             int                                   a_p)
{
  MayDay::Error("NOT CODED YET");
  Real val = 1.e99;
  return val;
}

Real EBAMRDataOps::noKappaNorm(Real&                                 a_volume,
                               const Vector<LevelData<EBCellFAB>* >& a_data,
                               int                                   a_which,
                               const Vector<ProblemDomain>&          a_domain,
                               int                                   a_p)
{
  MayDay::Error("NOT CODED YET");
  Real val = 1.e99;
  return val;
}

Real EBAMRDataOps::kappaDotProduct(Real&                                 a_volume,
                                   const Vector<LevelData<EBCellFAB>* >& a_data1,
                                   const Vector<LevelData<EBCellFAB>* >& a_data2,
                                   int                                   a_which,
                                   const Vector<ProblemDomain>&          a_domain)
{
  MayDay::Error("NOT CODED YET");
  Real val = 1.e99;
  return val;
}

Real EBAMRDataOps::noKappaDotProduct(Real&                       a_volume,
                                     const Vector<LevelData<EBCellFAB>* >& a_data1,
                                     const Vector<LevelData<EBCellFAB>* >& a_data2,
                                     int                         a_which,
                                     const Vector<ProblemDomain>&        a_domain)
{
  MayDay::Error("NOT CODED YET");
  Real val = 1.e99;
  return val;
}

Real
EBAMRDataOps::sum(const Vector<LevelData<EBCellFAB>* >& a_data,
                  const Vector<DisjointBoxLayout>&      a_grids,
                  const Vector<EBISLayout>&             a_ebisl,
                  const Vector<int> &                   a_refRat,
                  int   a_comp,
                  bool  a_multiplyByKappa)
{
  Real retval = 0.0;
  int maxlev = a_data.size()-1;

  for (int ilev  = 0; ilev <= maxlev; ilev++)
    {
      //don't count stuff covered by finer levels
      IntVectSet ivsExclude;
      if (ilev < maxlev)
        {
          //put next finer grids into an IVS
          for (LayoutIterator lit = a_grids[ilev+1].layoutIterator(); lit.ok(); ++lit)
            {
              ivsExclude |= a_grids[ilev+1].get(lit());
            }
          //coarsen back down to this level.
          //ivs will now hold the image of the next finer level
          ivsExclude.coarsen(a_refRat[ilev]);
        }

      Real sumLev;
      sumLev = EBLevelDataOps::sum(*(a_data[ilev]),  a_grids[ilev],
                                   a_ebisl[ilev],  ivsExclude,
                                   a_comp, a_multiplyByKappa);
      retval += sumLev;
    }

  return retval;
}

void EBAMRDataOps::setMaxMin(Vector<LevelData<EBCellFAB>* >& a_data,
                             const Real&                     a_maxVal,
                             const Real&                     a_minVal,
                             const int&                      a_comp)
{
  //this function sets the max and min values
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      for (DataIterator dit=a_data[ilev]->dataIterator();dit.ok();++dit)
        {
          EBCellFAB& dataEBFAB = (*a_data[ilev])[dit()];
          const Box& region = dataEBFAB.getRegion();
          const IntVectSet ivsBox(region);
          const EBISBox& ebisBox = dataEBFAB.getEBISBox();
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Real& val = dataEBFAB(vof,a_comp);
              val = Max(a_minVal,val);
              val = Min(a_maxVal,val);
            }
        }
    }
}

void EBAMRDataOps::getMaxMin(Real&                                 a_maxVal,
                             Real&                                 a_minVal,
                             const Vector<LevelData<EBCellFAB>* >& a_data,
                             const int&                            a_comp)
{
  //this function gets the max and min values
  a_maxVal = -1.e99;
  a_minVal =  1.e99;
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      Real maxThis = -1.e99;
      Real minThis =  1.e99;
      EBLevelDataOps::getMaxMin(maxThis,minThis,*a_data[ilev],a_comp);
      a_maxVal = Max(maxThis,a_maxVal);
      a_minVal = Min(minThis,a_minVal);
    }
}
Real EBAMRDataOps::subtractOffMean(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<DisjointBoxLayout>& a_grids,
                                   const Vector<EBISLayout>&        a_ebisl,
                                   const Vector<int> &              a_refRat)
{
  int pmode = -2;  //norm = (1/v)(sum(a_data dv)) ---no absolute values and multiply kappa as you go
  const Real mean = EBArith::norm(a_data, a_grids, a_ebisl, a_refRat, 0, pmode);

  //subtract off mean of data
  incr(a_data, -mean);
  return mean;
}

Real EBAMRDataOps::subtractOffMean(Vector<LevelData<EBCellFAB>* >&  a_data,
                                   const Vector<EBLevelGrid>&       a_eblg,
                                   const Vector<int> &              a_refRat)
{
  Vector<EBISLayout>        ebisl(a_eblg.size());
  Vector<DisjointBoxLayout> grids(a_eblg.size());
  Vector<ProblemDomain>     domain(a_eblg.size());
  for (int ilev = 0; ilev < a_eblg.size(); ilev++)
    {
      ebisl[ilev]  = a_eblg[ilev].getEBISL();
      grids[ilev]  = a_eblg[ilev].getDBL();
      domain[ilev] = a_eblg[ilev].getDomain();
    }

  int pmode = -2;  //norm = (1/v)(sum(a_data dv)) ---no absolute values and multiply kappa as you go
  const Real mean = EBArith::norm(a_data, grids, ebisl, a_refRat, 0, pmode);

  //subtract off mean of data
  incr(a_data, -mean);
  return mean;
}

void EBAMRDataOps::checkThisData(const Vector<LevelData<EBCellFAB>* >& a_data,
                                 const string&                         a_name,
                                 const IntVect&                        a_iv1,
                                 const IntVect&                        a_iv2,
                                 const Real&                           a_shift)
{
  bool isNANINF = false;

  isNANINF = checkNANINF(a_data,a_iv1,a_iv2,a_shift);
  if (isNANINF)
    {
      pout() << "failed naninf check, data name = " << a_name << std::endl;
      CH_assert(false);//this is easier to catch in the debugger
      MayDay::Error("failed naninf check  ");
    }
}

void EBAMRDataOps::checkThisData(const Vector<LevelData<EBFluxFAB>* >& a_data,
                                 const string&                         a_name)
{
  bool isNANINF = false;

  isNANINF = checkNANINF(a_data);
  if (isNANINF)
    {
      pout() << "failed naninf check, data name = " << a_name << std::endl;
      CH_assert(false);//this is easier to catch in the debugger
      MayDay::Error("failed naninf check  ");
    }
}

bool EBAMRDataOps::checkNANINF(const Vector<LevelData<EBCellFAB>* >& a_data,
                               const IntVect&                        a_iv1,
                               const IntVect&                        a_iv2,
                               const Real&                           a_shift)
{
  const bool checkSymmetry = (a_iv1 != a_iv2);

  //this function checks for nans and infs
  bool dataIsNANINF = false;
  int numLevels = a_data.size();
  int ncomp = a_data[0]->nComp();
  for (int icomp = 0; icomp < ncomp;++icomp)
    {
      for (int ilev = 0; ilev < numLevels; ilev++)
        {
          Real val1 = 1.e27;
          Real val2 = 1.e37;

          for (DataIterator dit=a_data[ilev]->dataIterator();dit.ok();++dit)
            {
              const EBCellFAB& dataEBFAB = (*a_data[ilev])[dit()];
              const Box& region = dataEBFAB.getRegion();
              IntVectSet ivsBox(region);
              const EBISBox& ebisBox = dataEBFAB.getEBISBox();
              for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
                {
                  const VolIndex& vof = vofit();
                  Real val = dataEBFAB(vof,icomp);
                  if (isnan(val) || isinf(val) || Abs(val)>1.e40)
                    {
                      pout() << "     ilev = " << ilev << " icomp = " << icomp << " vof = " << vof << " val = " << val << std::endl;
                      dataIsNANINF = true;
                    }
                  if (vof.gridIndex()==a_iv1)
                    {
                      val1 = Abs(val) - a_shift;
                    }
                  if (vof.gridIndex()==a_iv2)
                    {
                      val2 = Abs(val) - a_shift;
                    }
                }
            }
          //check symmetry
          if (checkSymmetry)
            {
              if ((val1-val2) > 1.e-8)
                {
                  pout() << "symmetry break: icomp = " << icomp << "; iv1 = " << a_iv1 << "; iv2 = " << a_iv2 << std::endl;
                  pout() << "val1 = " << val1 << std::endl;
                  pout() << "val2 = " << val2 << std::endl;
                  dataIsNANINF = true;
                }
            }
        }
      //pout the range
      const bool verbosity = false;
      if (verbosity)
        {
          Real maxVal,minVal;
          getMaxMin(maxVal,minVal,a_data,icomp);
          pout() << "max value = " << maxVal << " for comp = " << icomp << std::endl;
          pout() << "min value = " << minVal << " for comp = " << icomp << std::endl;
        }
    }

  if (dataIsNANINF)
    {
#ifdef CH_USE_HDF5
      char* fname = tempnam(NULL,NULL);
      writeEBAMRname(&a_data, fname);
      pout() << "   bad data written to an hdf5 file: " << fname << std::endl;
#else
      pout() << "   found a NaN or Infinity.  Recompile with HDF5 to get a data dump." << std::endl ;
#endif
    }

  return dataIsNANINF;
}

bool EBAMRDataOps::checkNANINF(const Vector<LevelData<EBFluxFAB>* >& a_data)
{
  //this function checks for nans and infs
  bool dataIsNANINF = false;
  int numLevels = a_data.size();
  int ncomp = a_data[0]->nComp();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      const DisjointBoxLayout& dbl = a_data[ilev]->getBoxes();
      for (DataIterator dit=a_data[ilev]->dataIterator();dit.ok();++dit)
        {
          const EBFluxFAB& dataFluxFAB = (*a_data[ilev])[dit()];
          const Box& dblBox = dbl[dit()];
          for (int idir=0;idir<SpaceDim;idir++)
            {
              const EBFaceFAB& dataFaceFAB = dataFluxFAB[idir];
              IntVectSet ivsBox(dblBox);
              const EBISBox& ebisBox = dataFaceFAB.getEBISBox();
              const EBGraph& ebGraph = ebisBox.getEBGraph();
              FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
              for (FaceIterator faceit(ivsBox, ebGraph, idir, stopCritGrid); faceit.ok(); ++faceit)
                {
                  const FaceIndex& face = faceit();
                  for (int icomp = 0; icomp < ncomp;++icomp)
                    {
                      Real val = dataFaceFAB(face,icomp);
                      if (isnan(val) || isinf(val) || Abs(val)>1.e40)
                        {
                          pout() << "     ilev = " << ilev << " icomp = " << icomp << " face = " << face << " val = " << val << std::endl;
                          dataIsNANINF = true;
                        }
                    }
                }
            }
        }
    }

  if (dataIsNANINF)
    {
#ifdef CH_USE_HDF5
      //char* fname = tempnam(NULL,NULL);
      //writeEBAMRname(&a_data, fname);
      //pout() << "   bad data written to an hdf5 file: " << fname << std::endl;
      pout() << "   bad face data...not sure how to output this data type...anyone know how?" << std::endl;
#else
      pout() << "   found a NaN or Infinity.  Recompile with HDF5 to get a data dump." << std::endl ;
#endif
    }

  return dataIsNANINF;
}

#include "NamespaceFooter.H"
