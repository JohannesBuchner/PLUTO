#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include "QuadCFInterp.H"
#include "BoxIterator.H"
#include "QuadCFInterpF_F.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "CH_Timer.H"
#include "CFIVS.H"
#include "TensorCFInterp.H"
#include "NamespaceHeader.H"
using std::endl;
/***********************/
// default constructor
/***********************/

bool QuadCFInterp::newCFInterMode = true;

/**/
void
QuadCFInterp::
interpPhiOnIVS(LevelData<FArrayBox>& a_phif,
               const FArrayBox& a_phistar,
               const DataIndex& a_datInd,
               const int a_idir,
               const Side::LoHiSide a_hiorlo,
               const IntVectSet& a_interpIVS,
               Real a_dxLevel,
               Real a_dxCrse,
               int a_ncomp)
{
  IVSIterator fine_ivsit(a_interpIVS);
  FArrayBox& a_phi = a_phif[a_datInd];

  Real x1 = a_dxLevel;
  Real x2 = 0.5*(3.*a_dxLevel+a_dxCrse);
  Real denom = 1.0-((x1+x2)/x1);
  Real idenom = 1.0/(denom); // divide is more expensive usually
  Real x = 2.*a_dxLevel;
  Real xsquared = x*x;

  Real m1 = 1/(x1*x1);
  Real m2 = 1/(x1*(x1-x2));

  Real q1 = 1/(x1-x2);
  Real q2 = x1+x2;

  int ihilo = sign(a_hiorlo);
  IntVect ai = -2*ihilo*BASISV(a_idir);
  IntVect bi = -  ihilo*BASISV(a_idir);
  IntVect ci =    ihilo*BASISV(a_idir);


  for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
    {
      IntVect ivf = fine_ivsit();
      // quadratic interpolation
      for (int ivar = 0; ivar < a_ncomp; ivar++)
        {
          Real pa =      a_phi(ivf + ai, ivar);
          Real pb =      a_phi(ivf + bi, ivar);
          Real pc =  a_phistar(ivf + ci, ivar);
          //phi = ax**2 + bx + c, x = 0 at pa
          Real a = (pb-pa)*m1 - (pb-pc)*m2;
          a *= idenom;
          Real b = (pb-pc)*q1 - a*q2;
          Real c = pa;
          a_phi(ivf,ivar) = a*xsquared + b*x + c;
        } //end loop over components
    } //end loop over fine intvects
}
/**/
void
QuadCFInterp::
homogeneousCFInterpPhi(LevelData<FArrayBox>& a_phif,
                       const DataIndex& a_datInd,
                       int a_idir,
                       Side::LoHiSide a_hiorlo,
                       LayoutData<CFIVS> a_loCFIVS[SpaceDim],
                       LayoutData<CFIVS> a_hiCFIVS[SpaceDim],
                       Real a_dxLevel,
                       Real a_dxCrse,
                       int a_ncomp)
{
  const CFIVS* cfivs_ptr = NULL;
  if (a_hiorlo == Side::Lo)
    cfivs_ptr = &a_loCFIVS[a_idir][a_datInd];
  else
    cfivs_ptr = &a_hiCFIVS[a_idir][a_datInd];

  const IntVectSet& interp_ivs = cfivs_ptr->getFineIVS();
  if (!interp_ivs.isEmpty())
    {
      int ihilo = sign(a_hiorlo);
      Box phistarbox = interp_ivs.minBox();
      phistarbox.shift(a_idir, ihilo);
      FArrayBox phistar(phistarbox, a_ncomp);
      //hence the homogeneous...
      phistar.setVal(0.);

      //given phistar, interpolate on fine ivs to fill ghost cells for phi
      interpPhiOnIVS(a_phif, phistar, a_datInd, a_idir, a_hiorlo,
                     interp_ivs, a_dxLevel, a_dxCrse, a_ncomp);
    }
}
/**/
void
QuadCFInterp::
homogeneousCFInterpTanGrad(LevelData<FArrayBox>& a_tanGrad,
                           const LevelData<FArrayBox>& a_phi,
                           const DataIndex& a_DatInd,
                           int a_idir,
                           Side::LoHiSide a_hiorlo,
                           Real a_dxLevel,
                           Real a_dxCrse,
                           int a_ncomp,
                           LayoutData<TensorFineStencilSet> a_loTanStencilSets[SpaceDim],
                           LayoutData<TensorFineStencilSet> a_hiTanStencilSets[SpaceDim])
{

  const TensorFineStencilSet* cfstencil_ptr = NULL;
  if (a_hiorlo == Side::Lo)
    cfstencil_ptr = &a_loTanStencilSets[a_idir][a_DatInd];
  else
    cfstencil_ptr = &a_hiTanStencilSets[a_idir][a_DatInd];

  Real x1 = a_dxLevel;
  Real x2 = 0.5*(3.*a_dxLevel+a_dxCrse);
  Real denom = 1.0-((x1+x2)/x1);
  Real idenom = 1.0/(denom); // divide is more expensive usually
  Real x = 2.*a_dxLevel;
  Real xsquared = x*x;

  Real m1 = 1/(x1*x1);
  Real m2 = 1/(x1*(x1-x2));

  Real q1 = 1/(x1-x2);
  Real q2 = x1+x2;

  const FArrayBox& phi = a_phi[a_DatInd];
  FArrayBox& tanGrad = a_tanGrad[a_DatInd];

  // loop over gradient directions
  for (int gradDir = 0; gradDir<SpaceDim; gradDir++)
    {
      if (gradDir != a_idir)
        {

          // first do centered stencil
          const IntVectSet& centeredIVS =
            cfstencil_ptr->getCenteredStencilSet(gradDir);

          int ihilo = sign(a_hiorlo);

          if (!centeredIVS.isEmpty())
            {
              // do centered computation
              IVSIterator cntrd_ivs(centeredIVS);
              // want to average fine-grid gradient with coarse
              // grid gradient, which is 0 (which is where the
              // extra factor of one-half comes from)

              Real gradMult = (0.5/a_dxLevel);
              for (cntrd_ivs.begin(); cntrd_ivs.ok(); ++cntrd_ivs)
                {
                  IntVect ivf = cntrd_ivs();
                  IntVect finePhiLoc = ivf - ihilo*BASISV(a_idir);
                  IntVect finePhiLoc2 = finePhiLoc - ihilo*BASISV(a_idir);
                  // loop over variables
                  for (int ivar = 0; ivar<a_phi.nComp(); ivar++)
                    {
                      Real fineHi = phi(finePhiLoc2+BASISV(gradDir),ivar);
                      Real fineLo = phi(finePhiLoc2-BASISV(gradDir),ivar);
                      Real fineGrada = gradMult*(fineHi-fineLo);

                      fineHi = phi(finePhiLoc+BASISV(gradDir),ivar);
                      fineLo = phi(finePhiLoc-BASISV(gradDir),ivar);
                      Real fineGradb = gradMult*(fineHi-fineLo);

                      // homogeneous interp implies that gradc is 0
                      Real gradc = 0;

                      int gradComp = TensorCFInterp::gradIndex(ivar,gradDir);

                      Real a = (fineGradb-fineGrada)*m1 - (fineGradb-gradc)*m2;
                      a *= idenom;
                      Real b = (fineGradb-gradc)*q1 - a*q2;
                      Real c = fineGrada;

                      tanGrad(ivf,gradComp) = a*xsquared + b*x + c;

                    }
                } // end loop over centered difference cells
            } // end if there are centered cells

          // now do forward-difference cells
          const IntVectSet& forwardIVS =
            cfstencil_ptr->getForwardStencilSet(gradDir);

          if (!forwardIVS.isEmpty())
            {
              // do forward-difference computations
              IVSIterator fwd_ivs(forwardIVS);
              // set up multipliers for gradient; since we want to average
              // fine-grid gradient with coarse-grid gradient (which is 0),
              // include an extra factor of one-half here.
              Real mult0 = -1.5/a_dxLevel;
              Real mult1 = 2.0/a_dxLevel;
              Real mult2 = -0.5/a_dxLevel;

              for (fwd_ivs.begin(); fwd_ivs.ok(); ++fwd_ivs)
                {
                  IntVect ivf = fwd_ivs();
                  IntVect finePhiLoc = ivf - ihilo*BASISV(a_idir);
                  IntVect finePhiLoc2 = finePhiLoc - ihilo*BASISV(a_idir);
                  //now loop overvariables
                  for (int var= 0; var<a_phi.nComp(); var++)
                    {
                      Real fine0 = phi(finePhiLoc2,var);
                      Real fine1 = phi(finePhiLoc2+BASISV(gradDir),var);
                      Real fine2 = phi(finePhiLoc2+2*BASISV(gradDir),var);
                      Real fineGrada = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      fine0 = phi(finePhiLoc,var);
                      fine1 = phi(finePhiLoc+BASISV(gradDir),var);
                      fine2 = phi(finePhiLoc+2*BASISV(gradDir),var);
                      Real fineGradb = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      Real gradc = 0.0;
                      int gradComp = TensorCFInterp::gradIndex(var,gradDir);
                      // now compute gradient

                      Real a = (fineGradb-fineGrada)*m1 - (fineGradb-gradc)*m2;
                      a *= idenom;
                      Real b = (fineGradb-gradc)*q1 - a*q2;
                      Real c = fineGrada;

                      tanGrad(ivf,gradComp) = a*xsquared + b*x + c;

                    } // end loop over variables
                } // end loop over forward-difference locations
            } // end if there are forward-difference cells

          // now do backward-difference cells
          const IntVectSet& backwardIVS =
            cfstencil_ptr->getBackwardStencilSet(gradDir);

          if (!backwardIVS.isEmpty())
            {
              IVSIterator back_ivs(backwardIVS);
              // set up multipliers for gradient -- since we want to average
              // fine-grid gradient with coarse-grid gradient (which is 0),
              // include an extra factor of one-half here.
              Real mult0 = -1.5/a_dxLevel;
              Real mult1 = 2.0/a_dxLevel;
              Real mult2 = -0.5/a_dxLevel;

              for (back_ivs.begin(); back_ivs.ok(); ++back_ivs)
                {
                  IntVect ivf = back_ivs();
                  IntVect finePhiLoc = ivf - ihilo*BASISV(a_idir);
                  IntVect finePhiLoc2 = finePhiLoc - ihilo*BASISV(a_idir);
                  // now loop over variables
                  for (int var=0; var<a_phi.nComp(); var++)
                    {
                      Real fine0 = phi(finePhiLoc2,var);
                      Real fine1 = phi(finePhiLoc2-BASISV(gradDir),var);
                      Real fine2 = phi(finePhiLoc2-2*BASISV(gradDir),var);
                      Real fineGrada = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      fine0 = phi(finePhiLoc,var);
                      fine1 = phi(finePhiLoc-BASISV(gradDir),var);
                      fine2 = phi(finePhiLoc-2*BASISV(gradDir),var);
                      Real fineGradb = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      Real gradc = 0.0;

                      int gradComp = TensorCFInterp::gradIndex(var,gradDir);


                      Real a = (fineGradb-fineGrada)*m1 - (fineGradb-gradc)*m2;
                      a *= idenom;
                      Real b = (fineGradb-gradc)*q1 - a*q2;
                      Real c = fineGrada;

                      tanGrad(ivf,gradComp) = a*xsquared + b*x + c;

                    } // end loop over variables
                } // end loop over backward-difference cells
            } // end if there are backward-difference cells


        } // end if gradDir is a tangential direction
    } // end loop over gradient directions

}
/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void
QuadCFInterp::
homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                    LevelData<FArrayBox>& a_tanGrad,
                    LayoutData<CFIVS> a_loCFIVS[SpaceDim],
                    LayoutData<CFIVS> a_hiCFIVS[SpaceDim],
                    Real a_dxLevel,
                    Real a_dxCrse,
                    int a_ncomp,
                    LayoutData<TensorFineStencilSet> a_loTanStencilSets[SpaceDim],
                    LayoutData<TensorFineStencilSet> a_hiTanStencilSets[SpaceDim])
{

  // need to do this to be sure that tangential derivatives are computed
  // correctly
  a_phif.exchange(a_phif.interval());
  DataIterator dit = a_phif.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();

      // first fill in cells for phi
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterpPhi(a_phif,datInd,idir,sit(),
                                     a_loCFIVS, a_hiCFIVS,
                                     a_dxLevel, a_dxCrse, a_ncomp);
            }
        }

      // now fill in tangential gradient cells
      for (int idir = 0; idir<SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterpTanGrad(a_tanGrad, a_phif,
                                         datInd,idir,sit(),
                                         a_dxLevel, a_dxCrse, a_ncomp,
                                         a_loTanStencilSets, a_hiTanStencilSets);
            }
        }
    }
}
void
QuadCFInterp::clear()
{
  m_isDefined = false;
  m_level = -1;
  m_dxFine = -1;

}

QuadCFInterp::QuadCFInterp()
{
  clear();
}

/***********************/
/***********************/
QuadCFInterp::QuadCFInterp(
                           const DisjointBoxLayout& a_fineBoxes,
                           const DisjointBoxLayout* a_coarBoxes,
                           Real  a_dxFine,
                           int a_refRatio, int a_nComp,
                           const Box& a_domf)
{
  ProblemDomain fineProbDomain(a_domf);
  define(a_fineBoxes,a_coarBoxes, a_dxFine,a_refRatio,a_nComp, fineProbDomain);
}

/***********************/
/***********************/
QuadCFInterp::QuadCFInterp(
                           const DisjointBoxLayout& a_fineBoxes,
                           const DisjointBoxLayout* a_coarBoxes,
                           Real  a_dxFine,
                           int a_refRatio, int a_nComp,
                           const ProblemDomain& a_domf)
{
  define(a_fineBoxes,a_coarBoxes, a_dxFine,a_refRatio,a_nComp, a_domf);
}

/***********************/
/***********************/

void
QuadCFInterp::define(
                     const DisjointBoxLayout& a_fineBoxes,
                     const DisjointBoxLayout* a_coarBoxesPtr,
                     Real  a_dxLevel,
                     int a_refRatio, int a_nComp,
                     const ProblemDomain& a_domf)
{
  CH_TIME("QuadCFInterp::define");
  clear();
  m_isDefined = true;
  CH_assert(a_nComp > 0);

  CH_assert (!a_domf.isEmpty());
  // consistency check
  CH_assert (a_fineBoxes.checkPeriodic(a_domf));

  m_domainFine = a_domf;
  m_dxFine = a_dxLevel;
  m_refRatio = a_refRatio;
  m_nComp = a_nComp;
  m_inputFineLayout = a_fineBoxes;
  bool fineCoversCoarse = false;
  if (a_coarBoxesPtr != NULL)
    {
      int factor = D_TERM6(a_refRatio, *a_refRatio, *a_refRatio,
                           *a_refRatio, *a_refRatio, *a_refRatio);
      long long numPts = a_fineBoxes.numCells()/factor;
      numPts -= a_coarBoxesPtr->numCells();
      if (numPts == 0) fineCoversCoarse = true;
    }
  m_fineCoversCoarse = fineCoversCoarse;
  if (a_coarBoxesPtr == NULL || fineCoversCoarse)
    m_level = 0;
  else
    m_level = 1;

  if (m_level > 0)
    {
      // (DFM) only check for valid refRatio if a coarser level exists
      CH_assert(a_refRatio >= 1);
      const DisjointBoxLayout& coarBoxes = *a_coarBoxesPtr;
      m_inputCoarLayout = coarBoxes;
      CH_assert (coarBoxes.checkPeriodic(coarsen(a_domf,a_refRatio)));
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_loQCFS[idir].define(a_fineBoxes);
          m_hiQCFS[idir].define(a_fineBoxes);

        }



      //locoarboxes and hicoarboxes are now open
      //and have same processor mapping as a_fineboxes
      //make boxes for coarse buffers
      m_coarBoxes.deepCopy(a_fineBoxes);
      m_coarBoxes.coarsen(m_refRatio);
      m_coarBoxes.grow(2);

      m_coarBoxes.close();
      m_coarBuffer.define(m_coarBoxes, m_nComp);
      m_copier.define(coarBoxes, m_coarBoxes);


      if (!newCFInterMode) //old n^2 algorithm (bvs)
        {
          //make cfstencils and boxes for coarse buffers
          DataIterator dit = a_fineBoxes.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box& fineBox = a_fineBoxes[dit()];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  //low side cfstencil
                  m_loQCFS[idir][dit()].define(a_domf,
                                               fineBox,
                                               a_fineBoxes,
                                               coarBoxes,
                                               a_refRatio,
                                               idir,
                                               Side::Lo);
                  //high side cfstencil
                  m_hiQCFS[idir][dit()].define(a_domf,
                                               fineBox,
                                               a_fineBoxes,
                                               coarBoxes,
                                               a_refRatio,
                                               idir,
                                               Side::Hi);

                }
            }
        } else
        {
          //new "moving window" version of CF stencil building
          Vector<Box> periodicFine;
          CFStencil::buildPeriodicVector(periodicFine, a_domf, a_fineBoxes);
          Vector<Box> coarsenedFine(periodicFine);
          for (int i=0; i<coarsenedFine.size(); ++i)
            {
              coarsenedFine[i].coarsen(a_refRatio);
            }
          DataIterator dit = a_fineBoxes.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box& fineBox = a_fineBoxes[dit()];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  //low side cfstencil
                  m_loQCFS[idir][dit()].define(a_domf,
                                               fineBox,
                                               periodicFine,
                                               coarsenedFine,
                                               coarBoxes,
                                               a_refRatio,
                                               idir,
                                               Side::Lo);
                  //high side cfstencil
                  m_hiQCFS[idir][dit()].define(a_domf,
                                               fineBox,
                                               periodicFine,
                                               coarsenedFine,
                                               coarBoxes,
                                               a_refRatio,
                                               idir,
                                               Side::Hi);

                }
            }
        }



    }
}

/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void
QuadCFInterp::coarseFineInterp(BaseFab<Real> & a_phif,
                               const BaseFab<Real> & a_phic,
                               const QuadCFStencil& a_qcfs,
                               const Side::LoHiSide a_hiorlo,
                               const int a_idir,
                               const Interval& a_variables) const
{
  CH_TIME("QuadCFInterp::coarseFineInterp(BaseFab<Real> & a_phif,...)");
  CH_assert(isDefined());
  //nothing happens if m_level == 0
  if (m_level > 0)
    {
      if (!a_qcfs.isEmpty())
        {
          BaseFab<Real>  phistar;

          //first find extended value phistar
          //includes finding slopes of coarse solution bar
          getPhiStar(phistar, a_phic, a_qcfs, a_hiorlo, a_idir, a_variables);

          //given phistar, interpolate on fine ivs
          interpOnIVS(a_phif, phistar, a_qcfs, a_hiorlo, a_idir, a_variables);
        }
    }
}

/***********************/
//get extended phi (lives next to interpivs)
/***********************/
void
QuadCFInterp::getPhiStar(BaseFab<Real> & a_phistar,
                         const BaseFab<Real> & a_phic,
                         const QuadCFStencil& a_qcfs,
                         const Side::LoHiSide a_hiorlo,
                         const int a_idir,
                         const Interval& a_variables) const
{
  CH_TIMERS("QuadCFInterp::getPhiStar");
  //CH_TIMER("QuadCFInterp::computeFirstDerivative",  t1st);
  //CH_TIMER("QuadCFInterp::computesecondDerivative", t2nd);
  //CH_TIMER("QuadCFInterp::computemixedDerivative",  tmixed);
  CH_TIMER("QuadCFInterp::slopes", tslopes);
  CH_TIMER("QuadCFInterp::notPacked", tnp);
  CH_TIMER("QuadCFInterp::preamble", tpreamble);
  CH_assert(isDefined());
  CH_assert(a_qcfs.isDefined());
#if (CH_SPACEDIM > 1)
  Real dxf = m_dxFine;
  Real dxc = m_refRatio*dxf;
#endif

  // if we think of a_idir as the "me" direction, then
  // the other directions can be "you1" and "you2"
#if (CH_SPACEDIM == 3)
  int you1, you2;
  if (a_idir == 0)
    {
      you1 = 1;
      you2 = 2;
    }
  else if (a_idir == 1)
    {
      you1 = 0;
      you2 = 2;
    }
  else
    {
      you1 = 0;
      you2 = 1;
    }
#else // (CH_SPACEDIM == 2)
  int you1;
  if (a_idir == 0)
    {
      you1 = 1;
    }
  else
    {
      you1 = 0;
    }
#endif

  //if cfsten is empty, nothing to interpolate.
  if (!a_qcfs.isEmpty())
    {
      CH_START(tpreamble);
      CH_assert(m_level > 0);
      const IntVectSet& interp_ivs = a_qcfs.getFineIVS();
      const IntVectSet& coarsl_ivs = a_qcfs.getCoarIVS();
      if (!coarsl_ivs.isDense())
        {
          MayDay::Error("What the hell?? TreeIntVectSet ???");
        }
      if (!interp_ivs.isDense())
        {
          MayDay::Error("What the hell?? TreeIntVectSet ???");
        }

      Box coarinterpbox = coarsl_ivs.minBox();
      int ncomp = a_phic.nComp();
      CH_assert(ncomp == m_nComp);
      CH_assert(a_phic.box().contains((coarinterpbox)));

      // allocate phistar here
      int ihilo = sign(a_hiorlo);
      Box phistarbox = interp_ivs.minBox();
      phistarbox.shift(a_idir, ihilo);
      a_phistar.resize(phistarbox, ncomp);
      CH_STOP(tpreamble);
      for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
        {

          CH_START(tslopes);
          //phi = phino + slope*x + half*x*x*curvature
          BaseFab<Real>  coarslope(coarinterpbox, SpaceDim);
          BaseFab<Real>  coarcurva(coarinterpbox, SpaceDim);
#if (CH_SPACEDIM == 3)
          BaseFab<Real>  coarmixed(coarinterpbox, 1);
#endif

          // coarslope.setVal(0.);

          //first find extended value phistar. get slopes of coarse solution
          IVSIterator coar_ivsit(coarsl_ivs);

          for (coar_ivsit.begin(); coar_ivsit.ok(); ++coar_ivsit)
            {
              // this isn't relevant for 1D
#if (CH_SPACEDIM > 1)
              const IntVect& coariv = coar_ivsit();

              // coarslope(coariv, a_idir) = 0.0;
              // coarcurva(coariv, a_idir) = 0.0;

              coarslope(coariv, you1) =
                a_qcfs.computeFirstDerivative (a_phic, you1, ivar, coariv, dxc);
              coarcurva(coariv, you1) =
                a_qcfs.computeSecondDerivative(a_phic, you1, ivar, coariv, dxc);

#endif

#if (CH_SPACEDIM == 3)
              coarslope(coariv, you2) =
                a_qcfs.computeFirstDerivative (a_phic, you2, ivar, coariv, dxc);
              coarcurva(coariv, you2) =
                a_qcfs.computeSecondDerivative(a_phic, you2, ivar, coariv, dxc);
              coarmixed(coariv) =
                a_qcfs.computeMixedDerivative(a_phic, ivar, coariv, dxc);

#endif
            } //end loop over coarse intvects
          CH_STOP(tslopes);
          if (a_qcfs.finePacked() && CH_SPACEDIM==3)
            {
              const IntVect& iv = phistarbox.smallEnd();
              IntVect civ(iv);
              civ.coarsen(m_refRatio);
              Box region = a_qcfs.packedBox();
#if (CH_SPACEDIM == 3)
              FORT_PHISTAR(CHF_FRA_SHIFT(a_phistar, iv),
                           CHF_BOX_SHIFT(region, iv),
                           CHF_CONST_FRA_SHIFT(a_phic, civ),
                           CHF_FRA_SHIFT(coarslope, civ),
                           CHF_FRA_SHIFT(coarcurva, civ),
                           CHF_FRA_SHIFT(coarmixed, civ),
                           CHF_CONST_REAL(dxf),
                           CHF_CONST_INT(ivar),
                           CHF_CONST_INT(a_idir),
                           CHF_CONST_INT(ihilo),
                           CHF_CONST_INT(m_refRatio));
#endif
            }
          else
            {
              CH_START(tnp);
              IntVect ivf, ivc, ivstar;
              // ifdef is here to prevent unused variable wasrnings in 1D
#if (CH_SPACEDIM > 1)
              int jf, jc;
              Real xf, xc, x1;
#endif
              Real pc, update1=0, update2=0, update3=0;
              IVSIterator fine_ivsit(interp_ivs);
              for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
                {
                  ivf = fine_ivsit();
                  ivc = coarsen(ivf, m_refRatio);
                  ivstar = ivf;
                  ivstar.shift(a_idir, ihilo);
                  pc = a_phic(ivc,ivar);

                  // for 1D, none of this is necessary -- just copy
                  // coarse value into phiStar
#if (CH_SPACEDIM > 1)
                  jf = ivf[you1];
                  jc = ivc[you1];
                  xf = (jf+0.5)*dxf;
                  xc = (jc+0.5)*dxc;
                  x1 = xf-xc;
                  update1= x1*coarslope(ivc, you1) + 0.5*x1*x1*coarcurva(ivc, you1);

#endif

#if (CH_SPACEDIM==3)
                  Real x2;
                  jf = ivf[you2];
                  jc = ivc[you2];
                  xf = (jf+0.5)*dxf;
                  xc = (jc+0.5)*dxc;
                  x2 = xf-xc;
                  update2 =  x2*coarslope(ivc, you2) + 0.5*x2*x2*coarcurva(ivc, you2);

                  //add in mixed derivative component
                  update3 =  x1*x2*coarmixed(ivc);
#endif
                  a_phistar(ivstar, ivar) = pc+update1+update2+update3;
                } //end loop over fine intvects
              CH_STOP(tnp);
            } // end if for not packed optimization

        }//end loop over variables


    } //end if (level>0 && !hocfs.isempty())
} //end function getphistar

/***********************/
/***********************/
bool
QuadCFInterp::isDefined() const
{
  return m_isDefined;
}

void
QuadCFInterp::interpOnIVS(BaseFab<Real> & a_phif,
                          const BaseFab<Real> & a_phistar,
                          const QuadCFStencil& a_qcfs,
                          const Side::LoHiSide a_hiorlo,
                          const int a_idir,
                          const Interval& a_variables) const
{
  CH_TIME("QuadCFInterp::interpOnIVS");
  CH_assert(isDefined());
  CH_assert(a_qcfs.isDefined());

  //if cfsten is empty, nothing to interpolate.
  if (!a_qcfs.isEmpty())
    {
      //if there IS something to interpolate, the level ident
      //had better be greater than zero.  Otherwise a null
      //was sent in as coarse grids on construction
      CH_assert(m_level > 0);
      const IntVectSet& interp_ivs = a_qcfs.getFineIVS();

      int ihilo = sign(a_hiorlo);
      int nref = m_refRatio;
      if (!a_qcfs.finePacked())
        {
          IVSIterator fine_ivsit(interp_ivs);
          CH_assert(a_phistar.nComp() == a_phif.nComp());
          CH_assert(a_phistar.nComp() == m_nComp);



          for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
            {
              IntVect ivf = fine_ivsit();
              // quadratic interpolation
              for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                {
                  Real pa =    a_phif (ivf -2*ihilo*BASISV(a_idir), ivar);
                  Real pb =    a_phif (ivf -  ihilo*BASISV(a_idir), ivar);
                  Real ps = a_phistar(ivf +  ihilo*BASISV(a_idir), ivar);
                  //phi = ax**2 + bx + c, x = 0 at pa
                  Real h = m_dxFine;
                  Real a = (2./h/h)*(2.*ps + pa*(nref+1.0) -pb*(nref+3.0))/
                    (nref*nref + 4*nref + 3.0);
                  Real b = (pb-pa)/h - a*h;
                  Real c = pa;
                  Real x = 2.*h;
                  a_phif (ivf,ivar) = a*x*x + b*x + c;
                } //end loop over components
            } //end loop over fine intvects
        } else
        { // data is packed, just call Fortran
          int b=a_variables.begin();
          int e=a_variables.end();
          FORT_QUADINTERP(CHF_FRA(a_phif),
                          CHF_CONST_FRA(a_phistar),
                          CHF_BOX(a_qcfs.packedBox()),
                          CHF_CONST_INT(ihilo),
                          CHF_CONST_REAL(m_dxFine),
                          CHF_CONST_INT(a_idir),
                          CHF_CONST_INT(b),
                          CHF_CONST_INT(e),
                          CHF_CONST_INT(nref));
        }
    } //end if (level>0 && !oscfs.isempty())
} //end function interponivs

/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void
QuadCFInterp::coarseFineInterp(LevelData<FArrayBox>& a_phif,
                               const LevelData<FArrayBox>& a_phic)
{
  CH_TIME("QuadCFInterp::coarseFineInterp");
  CH_assert(isDefined());


  Interval variables = a_phic.interval();

  if (m_level > 0)
    {
      CH_assert(a_phif.nComp() == m_nComp);
      CH_assert(a_phic.nComp() == m_nComp);
      CH_assert(a_phif.ghostVect() >= IntVect::Unit);
      CH_assert(a_phic.boxLayout() == m_inputCoarLayout);
      CH_assert(a_phif.boxLayout() == m_inputFineLayout);
      a_phic.copyTo(a_phic.interval(), m_coarBuffer, m_coarBuffer.interval(),
                    m_copier);
      for (int idir = 0; idir < SpaceDim; idir++)
        {


          DataIterator ditFine = a_phif.dataIterator();
          for (ditFine.begin(); ditFine.ok(); ++ditFine)
            {
              DataIndex datIndGlo =ditFine();
              BaseFab<Real> & phif = a_phif[datIndGlo];
              const BaseFab<Real> & phiC = m_coarBuffer[datIndGlo];
              //lo side cfinterp
              //recall that buffers have fine processor mapping
              {
                const QuadCFStencil& loQCFS = m_loQCFS[idir][datIndGlo];
                coarseFineInterp(phif, phiC ,loQCFS, Side::Lo, idir, variables);

              }
              //hi side cfinterp
              {
                const QuadCFStencil& hiQCFS = m_hiQCFS[idir][datIndGlo];
                coarseFineInterp(phif, phiC, hiQCFS, Side::Hi, idir, variables);
              }
            }//end iteration over boxes in fine grid
        } //end iteration over directions
    }
}
/***********************/
/***********************/
QuadCFInterp::~QuadCFInterp()
{
  clear();
}
#include "NamespaceFooter.H"
