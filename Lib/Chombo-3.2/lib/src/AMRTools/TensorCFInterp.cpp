#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TensorCFInterp.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "NamespaceHeader.H"

/***********************/
// default constructor
/***********************/
void
TensorCFInterp::clear()
{
  m_isDefined = false;
  m_level = -1;
  m_dxFine = -1;
}

TensorCFInterp::TensorCFInterp()
{
  clear();
}

/***********************/
/***********************/
TensorCFInterp::TensorCFInterp(
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
TensorCFInterp::TensorCFInterp(
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
TensorCFInterp::define(
                       const DisjointBoxLayout& a_fineBoxes,
                       const DisjointBoxLayout* a_coarBoxesPtr,
                       Real  a_dxLevel,
                       int a_refRatio, int a_nComp,
                       const Box& a_domf)
{
  ProblemDomain fineProbDomain(a_domf);
  define(a_fineBoxes, a_coarBoxesPtr, a_dxLevel, a_refRatio, a_nComp,
         a_domf);

}


/***********************/
/***********************/

void
TensorCFInterp::define(
                       const DisjointBoxLayout& a_fineBoxes,
                       const DisjointBoxLayout* a_coarBoxesPtr,
                       Real  a_dxLevel,
                       int a_refRatio, int a_nComp,
                       const ProblemDomain& a_domf)
{
  clear();
  m_isDefined = true;
  CH_assert(a_nComp > 0);

  CH_assert (!a_domf.isEmpty());
  // consistency check
  CH_assert (a_fineBoxes.checkPeriodic(a_domf));

  m_dxFine = a_dxLevel;
  m_refRatio = a_refRatio;
  m_nComp = a_nComp;
  m_probDomain = a_domf;

  m_inputFineLayout = a_fineBoxes;
  bool fineCoversCoarse = false;

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

      coarsen(m_coarsenedFineBoxes, m_inputFineLayout, a_refRatio);

      m_coarsenedFineBuffer.define(m_coarsenedFineBoxes, a_nComp, IntVect::Unit);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_loQCFS[idir].define(a_fineBoxes);
          m_hiQCFS[idir].define(a_fineBoxes);
          m_loFineStencilSets[idir].define(a_fineBoxes);
          m_hiFineStencilSets[idir].define(a_fineBoxes);
        }

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


              // now define stencils for fine-level gradient computation
              const IntVectSet& loFineIVS = m_loQCFS[idir][dit()].getFineIVS();
              m_loFineStencilSets[idir][dit()].define(loFineIVS,
                                                      a_domf,
                                                      idir);

              const IntVectSet& hiFineIVS = m_hiQCFS[idir][dit()].getFineIVS();
              m_hiFineStencilSets[idir][dit()].define(hiFineIVS,
                                                      a_domf,
                                                      idir);

            }
        }
    }
}


/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void
TensorCFInterp::coarseFineInterp(BaseFab<Real> & a_phif,
                                 BaseFab<Real> & a_gradf,
                                 BaseFab<Real> & a_tanGradStar,
                                 const BaseFab<Real> & a_phic,
                                 const QuadCFStencil& a_qcfs,
                                 const Side::LoHiSide a_hiorlo,
                                 const int a_idir,
                                 const Interval& a_variables) const
{
  CH_assert(isDefined());
  //nothing happens if m_level == 0
  if (m_level > 0)
    {
      if (!a_qcfs.isEmpty())
        {
          BaseFab<Real>  phistar;

          //first find extended value phistar
          //includes finding slopes of coarse solution bar
          // tanGradStar will also be computed here
          getPhiStar(phistar, a_tanGradStar, a_phic, a_qcfs,
                     a_hiorlo, a_idir, a_variables);

          //given phistar, interpolate on fine ivs
          interpOnIVS(a_phif, a_gradf, phistar, a_qcfs, a_hiorlo,
                      a_idir, a_variables);
        }
    }
}

/***********************/
//get extended phi (lives next to interpivs)
/***********************/
void
TensorCFInterp::getPhiStar(BaseFab<Real> & a_phistar,
                           BaseFab<Real> & a_tanGradStar,
                           const BaseFab<Real> & a_phic,
                           const QuadCFStencil& a_qcfs,
                           const Side::LoHiSide a_hiorlo,
                           const int a_idir,
                           const Interval& a_variables) const
{

  CH_assert(isDefined());
  CH_assert(a_qcfs.isDefined());
  Real dxf = m_dxFine;
  Real dxc = m_refRatio*dxf;

  //if cfsten is empty, nothing to interpolate.
  if (!a_qcfs.isEmpty())
    {
      //if there IS something to interpolate, the level ident
      //had better be greater than zero.  Otherwise a null
      //was sent in as coarse grids on construction
      CH_assert(m_level > 0);
      const IntVectSet& interp_ivs = a_qcfs.getFineIVS();
      const IntVectSet& coarsl_ivs = a_qcfs.getCoarIVS();
      Box coarinterpbox = coarsl_ivs.minBox();
      int ncomp = a_phic.nComp();
      CH_assert(ncomp == m_nComp);
      CH_assert(a_phic.box().contains((coarinterpbox)));

      // allocate phistar here
      int ihilo = sign(a_hiorlo);
      Box phistarbox = interp_ivs.minBox();
      phistarbox.shift(a_idir, ihilo);
      a_phistar.resize(phistarbox, ncomp);
      a_phistar.setVal(666.666);
      // tanGradStar has same location as phiStar
      // notice that we're allocating SpaceDim components for
      // each component in phi (for each coordinate direction).
      // possible memory optimization would only allocate SpaceDim-1
      // components (won't need normal dir)
      a_tanGradStar.resize(phistarbox, ncomp*SpaceDim);
      a_tanGradStar.setVal(666.666);

      for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
        {
          //phi = phino + slope*x + half*x*x*curvature
          BaseFab<Real>  coarslope(coarinterpbox, SpaceDim);
          BaseFab<Real>  coarcurva(coarinterpbox, SpaceDim);

#if (CH_SPACEDIM == 3)
          //directions != idir,  needed for mixed derivative in 3d
          //add to phi x*y*mixed_derivative
          Vector<int> xydir;
          for (int isldir = 0; isldir < SpaceDim; isldir++)
            if (isldir != a_idir)
              xydir.push_back(isldir);
          BaseFab<Real>  coarmixed(coarinterpbox, 1);
#endif

          coarslope.setVal(0.);

          //first find extended value phistar
          //get slopes of coarse solution
          IVSIterator coar_ivsit(coarsl_ivs);
          for (coar_ivsit.begin(); coar_ivsit.ok(); ++coar_ivsit)
            {
              IntVect coariv = coar_ivsit ();
              for (int isldir = 0; isldir < SpaceDim; isldir++)
                {
                  if (isldir == a_idir)
                    {
                      coarslope(coariv, isldir) = 0.0;
                      coarcurva(coariv, isldir) = 0.0;
                    }
                  else
                    {
                      coarslope(coariv, isldir) =
                        a_qcfs.computeFirstDerivative(a_phic, isldir,
                                                      ivar, coariv, dxc);
                      coarcurva(coariv, isldir) =
                        a_qcfs.computeSecondDerivative(a_phic, isldir,
                                                       ivar, coariv, dxc);
                    } //end if (isldir == idir)
                } // end loop over slope dir
#if (CH_SPACEDIM == 3)
              coarmixed(coariv) =
                a_qcfs.computeMixedDerivative(a_phic, ivar, coariv, dxc);
#endif
            } //end loop over coarse intvects

          IVSIterator fine_ivsit(interp_ivs);
          for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
            {
              IntVect ivf = fine_ivsit();
              IntVect ivc = coarsen(ivf,m_refRatio);
              IntVect ivstar = ivf;
              ivstar.shift(a_idir, ihilo);
              a_phistar(ivstar, ivar)  = a_phic(ivc,ivar);
              // initialize tanGrad to derivatives at coarse point
              for (int gradDir=0; gradDir<SpaceDim; gradDir++)
                {
                  if (gradDir != a_idir)
                    {
                      // grad will be stored so that dimensions
                      // vary fastest
                      int gradComp = gradIndex(ivar, gradDir);
                      a_tanGradStar(ivstar, gradComp) = coarslope(ivc,gradDir);
                    }
                }


              for (int isldir = 0; isldir < SpaceDim; isldir++)
                {
                  if (isldir != a_idir)
                    {
                      int jf = ivf[isldir];
                      int jc = ivc[isldir];
                      Real xf = (jf+0.5)*dxf;
                      Real xc = (jc+0.5)*dxc;
                      Real x = xf-xc;
                      a_phistar(ivstar, ivar) +=
                        (      x*coarslope(ivc, isldir)) +
                        (0.5*x*x*coarcurva(ivc, isldir));
                      // increment tangential gradient
                      int gradComp = gradIndex(ivar, isldir);
                      a_tanGradStar(ivstar, gradComp) +=
                        x*coarcurva(ivc,isldir);

                    }//end if (isldir != idir)
                } //end loop over slope directions
#if (CH_SPACEDIM==3)
              //add in mixed derivative component
              int dif = ivf[xydir[0]];
              int dic = ivc[xydir[0]];
              int djf = ivf[xydir[1]];
              int djc = ivc[xydir[1]];

              Real xf = (dif+0.5)*dxf;
              Real xc = (dic+0.5)*dxc;
              Real yf = (djf+0.5)*dxf;
              Real yc = (djc+0.5)*dxc;
              Real x = xf-xc;
              Real y = yf-yc;
              a_phistar(ivstar, ivar) += x*y*coarmixed(ivc);
              // mixed derivative gets added to each component of
              // tanGrad.  this is simple enough to do in 3d --
              // we can just hardwire it.  if higher dimensional
              // problem, this needs to be rethought.
              CH_assert (SpaceDim <= 3);

              // first do direction 0
              int gradComp = gradIndex(ivar,xydir[0]);
              a_tanGradStar(ivstar, gradComp) += y*coarmixed(ivc);
              // now do direction 1
              gradComp = gradIndex(ivar,xydir[1]);
              a_tanGradStar(ivstar, gradComp) += x*coarmixed(ivc);


#endif
            } //end loop over fine intvects
        } //end loop over variables
    } //end if (level>0 && !hocfs.isempty())
} //end function getphistar

/***********************/
/***********************/
bool
TensorCFInterp::isDefined() const
{
  return m_isDefined;
}

void
TensorCFInterp::interpOnIVS(BaseFab<Real> & a_phif,
                            BaseFab<Real> & a_gradf,
                            const BaseFab<Real> & a_phistar,
                            const QuadCFStencil& a_qcfs,
                            const Side::LoHiSide a_hiorlo,
                            const int a_idir,
                            const Interval& a_variables) const
{
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
      IVSIterator fine_ivsit(interp_ivs);
      CH_assert(a_phistar.nComp() == a_phif.nComp());
      CH_assert(a_phistar.nComp() == m_nComp);

      int nref = m_refRatio;
      int ihilo = sign(a_hiorlo);
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

              int gradVar = gradIndex(ivar,a_idir);
              a_gradf(ivf, gradVar) = ihilo*(2.*a*x + b);
            } //end loop over components
        } //end loop over fine intvects
    } //end if (level>0 && !oscfs.isempty())
} //end function interponivs


/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void
TensorCFInterp::coarseFineInterp(LevelData<FArrayBox>& a_phif,
                                 LevelData<FArrayBox>& a_gradf,
                                 const LevelData<FArrayBox>& a_phic)
{
  CH_assert(isDefined());
  CH_assert(a_phif.nComp() == m_nComp);
  CH_assert(a_phic.nComp() == m_nComp);
  CH_assert(a_phif.ghostVect() >= IntVect::Unit);
  CH_assert(a_phic.boxLayout() == m_inputCoarLayout);
  CH_assert(a_phif.boxLayout() == m_inputFineLayout);

  Interval variables = a_phic.interval();

  if (m_level > 0)
    {
      LevelData<FArrayBox>& phif = (LevelData<FArrayBox>&)(a_phif);
      LevelData<FArrayBox>& phic = (LevelData<FArrayBox>&)(a_phic);
      phif.exchange();
      phic.exchange();
      Interval interv(0, m_nComp-1);
      a_phic.copyTo(interv, m_coarsenedFineBuffer, interv);

      DataIterator ditFine = a_phif.dataIterator();
      for (ditFine.begin(); ditFine.ok(); ++ditFine)
        {
          DataIndex datIndGlo =ditFine();
          BaseFab<Real> & phif = a_phif[datIndGlo];
          BaseFab<Real> & gradf = a_gradf[datIndGlo];

          // this will contain the coarse-grid gradients,
          // centered in the same locations as phiStar.
          // this is necessary because cannot compute interpolated
          // grad(phi) until ghost cell values are computed.
          // these will be allocated at the same time as phiStar
          BaseFab<Real> gradStar[SpaceDim][2];
          // loop over directions
          for (int idir=0; idir<SpaceDim; idir++)
            {
              //lo side cfinterp
              //recall that buffers have fine processor mapping
              {
                const BaseFab<Real> & loPhiC = m_coarsenedFineBuffer[datIndGlo];
                const QuadCFStencil& loQCFS = m_loQCFS[idir][datIndGlo];
                BaseFab<Real>& thisGradStar = gradStar[idir][0];
                coarseFineInterp(phif, gradf, thisGradStar, loPhiC,loQCFS,
                                 Side::Lo, idir, variables);
              }
              //hi side cfinterp
              {
                const BaseFab<Real> & hiPhiC = m_coarsenedFineBuffer[datIndGlo];
                const QuadCFStencil& hiQCFS = m_hiQCFS[idir][datIndGlo];
                BaseFab<Real>& thisGradStar = gradStar[idir][1];
                coarseFineInterp(phif, gradf, thisGradStar, hiPhiC, hiQCFS,
                                 Side::Hi, idir, variables);
              }
            } //end iteration over directions
          // now can compute tangential gradients
          BaseFab<Real> & tanGrad = a_gradf[ditFine()];

          for (int idir=0; idir<SpaceDim; idir++)
            {
              // lo side
              BaseFab<Real>& thisGradStarLo = gradStar[idir][0];
              const QuadCFStencil& loQCFS = m_loQCFS[idir][datIndGlo];
              const TensorFineStencilSet& loTFSS
                = m_loFineStencilSets[idir][datIndGlo];
              computeTanGrad(tanGrad, phif, thisGradStarLo,
                             loTFSS, loQCFS, Side::Lo, idir, variables);

              // hi side
              BaseFab<Real>& thisGradStarHi = gradStar[idir][1];
              const QuadCFStencil& hiQCFS = m_hiQCFS[idir][datIndGlo];
              const TensorFineStencilSet& hiTFSS
                = m_hiFineStencilSets[idir][datIndGlo];
              computeTanGrad(tanGrad, phif, thisGradStarHi,
                             hiTFSS, hiQCFS,Side::Hi, idir, variables);
            }

        }//end iteration over boxes in fine grid
    }
}


/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void
TensorCFInterp::coarseFineInterpH(LevelData<FArrayBox>& a_phif,
                                  LevelData<FArrayBox>& a_gradf)
{
  CH_assert(isDefined());
  CH_assert(a_phif.nComp() == m_nComp);
  CH_assert(a_phif.ghostVect() >= IntVect::Unit);
  CH_assert(a_phif.boxLayout() == m_inputFineLayout);

  Interval variables = a_phif.interval();

  if (m_level > 0)
    {
      DataIterator ditFine = a_phif.dataIterator();
      for (ditFine.begin(); ditFine.ok(); ++ditFine)
        {
          m_coarsenedFineBuffer[ditFine()].setVal(0.);
        }

      for (ditFine.begin(); ditFine.ok(); ++ditFine)
        {
          DataIndex datIndGlo =ditFine();
          BaseFab<Real> & phif = a_phif[datIndGlo];
          BaseFab<Real> & gradf = a_gradf[datIndGlo];

          // this will contain the coarse-grid gradients,
          // centered in the same locations as phiStar.
          // this is necessary because cannot compute interpolated
          // grad(phi) until ghost cell values are computed.
          // these will be allocated at the same time as phiStar
          BaseFab<Real> gradStar[SpaceDim][2];
          // loop over directions
          for (int idir=0; idir<SpaceDim; idir++)
            {
              //lo side cfinterp
              //recall that buffers have fine processor mapping
              {
                const BaseFab<Real> & loPhiC = m_coarsenedFineBuffer[datIndGlo];
                const QuadCFStencil& loQCFS = m_loQCFS[idir][datIndGlo];
                BaseFab<Real>& thisGradStar = gradStar[idir][0];
                coarseFineInterp(phif, gradf, thisGradStar, loPhiC,loQCFS,
                                 Side::Lo, idir, variables);
              }
              //hi side cfinterp
              {
                const BaseFab<Real> & hiPhiC = m_coarsenedFineBuffer[datIndGlo];
                const QuadCFStencil& hiQCFS = m_hiQCFS[idir][datIndGlo];
                BaseFab<Real>& thisGradStar = gradStar[idir][1];
                coarseFineInterp(phif,gradf, thisGradStar, hiPhiC, hiQCFS,
                                 Side::Hi, idir, variables);
              }
            } //end iteration over directions
          // now can compute tangential gradients
          BaseFab<Real> & tanGrad = a_gradf[ditFine()];

          for (int idir=0; idir<SpaceDim; idir++)
            {
              // lo side
              BaseFab<Real>& thisGradStarLo = gradStar[idir][0];
              const QuadCFStencil& loQCFS = m_loQCFS[idir][datIndGlo];
              const TensorFineStencilSet& loTFSS
                = m_loFineStencilSets[idir][datIndGlo];
              computeTanGrad(tanGrad, phif, thisGradStarLo,
                             loTFSS, loQCFS, Side::Lo, idir, variables);

              // hi side
              BaseFab<Real>& thisGradStarHi = gradStar[idir][1];
              const QuadCFStencil& hiQCFS = m_hiQCFS[idir][datIndGlo];
              const TensorFineStencilSet& hiTFSS
                = m_hiFineStencilSets[idir][datIndGlo];
              computeTanGrad(tanGrad, phif, thisGradStarHi,
                             hiTFSS, hiQCFS,Side::Hi, idir, variables);
            }

        }//end iteration over boxes in fine grid
    }
}

void
TensorCFInterp::computeTanGrad(BaseFab<Real>& a_tanGradf,
                               const BaseFab<Real>& a_phif,
                               const BaseFab<Real>& a_tanGradStar,
                               const TensorFineStencilSet& a_TFSS,
                               const QuadCFStencil& a_qcfs,
                               const Side::LoHiSide& a_hiorlo,
                               const int a_idir,
                               const Interval& a_vars) const
{

  CH_assert(isDefined());
  CH_assert(a_TFSS.isDefined());

  // if cfsten is empty, nothing to interpolate
  if (!a_TFSS.isEmpty())
    {
      // if there IS something to interpolate, the level ident
      // had better be greater than zero.  Otherwise, a null
      // was sent in as coarse grids on construction
      CH_assert(m_level >0);

      // first do centered cells -- since extra work from computing
      // all cells as centered, then going back and fixing the ones
      // which need to be one-sided is probably minimal, it's probably
      // faster to loop over all cells in the qcfStencil once here,
      // then go back and fix up the cells which abut physical boundaries,
      // since that will result in traversing the intVect set only once,
      // rather than SpaceDim-1 times (yeah, i know that 2-1 is 1, so
      // no real benefit in 2d, but should help in 3d
      const IntVectSet& interp_ivs = a_qcfs.getFineIVS();
      IVSIterator fine_ivsit(interp_ivs);
      CH_assert(a_tanGradf.nComp() == a_tanGradStar.nComp());
      CH_assert(a_tanGradf.nComp() == SpaceDim*m_nComp);


      int ihilo = sign(a_hiorlo);
      Real gradFactor = (0.5/m_dxFine);
      // Real crseMult = 2.0/(1+m_refRatio);
      // Real fineMult = 1.0 - crseMult;
      int nref = m_refRatio;

      for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
        {
          IntVect ivf = fine_ivsit();
          IntVect finePhiLoc = ivf - ihilo*BASISV(a_idir);
          IntVect finePhiLoc2 = finePhiLoc - ihilo*BASISV(a_idir);
          // loop over variables
          for (int ivar=a_vars.begin(); ivar<=a_vars.end(); ivar++)
            {
              // now need to loop over tangential directions to compute
              // derivatives.
              for (int gradDir=0; gradDir<SpaceDim; gradDir++)
                {
                  if (gradDir != a_idir)
                    {
                      Real fineHi = a_phif (finePhiLoc+BASISV(gradDir),ivar);
                      Real fineLo = a_phif (finePhiLoc-BASISV(gradDir),ivar);
                      Real fineGradb = gradFactor*(fineHi-fineLo);

                      Real fineHi2 = a_phif (finePhiLoc2+BASISV(gradDir),ivar);
                      Real fineLo2 = a_phif (finePhiLoc2-BASISV(gradDir),ivar);
                      Real fineGrada = gradFactor*(fineHi2-fineLo2);


                      int gradComp = gradIndex(ivar,gradDir);
                      Real grads = a_tanGradStar(ivf+ihilo*BASISV(a_idir),
                                                 gradComp);

                      Real h = m_dxFine;
                      Real a = (2./h/h)*(2.*grads
                                         +fineGrada*(nref+1.0)
                                         -fineGradb*(nref+3.0))/
                        (nref*nref + 4*nref + 3.0);
                      Real b = (fineGradb-fineGrada)/h - a*h;
                      Real c = fineGrada;
                      Real x = 2.*h;
                      a_tanGradf(ivf,gradComp) = a*x*x + b*x + c;

#if 0
                      // now average fine-grid gradient with coarse-grid
                      // interpolated gradient
                      a_tanGradf(ivf, gradComp) = fineMult*fineGradb
                        +crseMult*a_tanGradStar(ivf+ihilo*BASISV(a_idir),
                                                gradComp);
#endif
                    }
                } // end loop over gradient directions
            } // end loop over variables
        } // end loop over fine intVectSet

      // now go back and do shifted stencils, if appropriate
      for (int gradDir=0; gradDir<SpaceDim; gradDir++)
        {
          if (gradDir != a_idir)
            {
              const IntVectSet& forIVS = a_TFSS.getForwardStencilSet(gradDir);
              if (!forIVS.isEmpty())
                {
                  IVSIterator fwd_ivsit(forIVS);

                  Real mult0 = -1.5/m_dxFine;
                  Real mult1 = 2.0/m_dxFine;
                  Real mult2 = -0.5/m_dxFine;
                  for (fwd_ivsit.begin(); fwd_ivsit.ok(); ++fwd_ivsit)
                    {
                      IntVect ivf = fwd_ivsit();
                      IntVect finePhiLoc = ivf - ihilo*BASISV(a_idir);
                      IntVect finePhiLoc2 = finePhiLoc - ihilo*BASISV(a_idir);
                      // loop over variables
                      for (int var = a_vars.begin(); var<=a_vars.end(); ++var)
                      {
                        Real fine0 = a_phif (finePhiLoc2,var);
                        Real fine1 = a_phif (finePhiLoc2+BASISV(gradDir),var);
                        Real fine2 = a_phif (finePhiLoc2+2*BASISV(gradDir),var);
                        Real fineGrada = mult0*fine0 +mult1*fine1 +mult2*fine2;

                        fine0 = a_phif (finePhiLoc,var);
                        fine1 = a_phif (finePhiLoc+BASISV(gradDir),var);
                        fine2 = a_phif (finePhiLoc+2*BASISV(gradDir),var);
                        Real fineGradb = mult0*fine0 +mult1*fine1 +mult2*fine2;

                        int gradComp = gradIndex(var,gradDir);
                        Real grads = a_tanGradStar(ivf+ihilo*BASISV(a_idir),
                                                   gradComp);



                        Real h = m_dxFine;
                        Real a = (2./h/h)*(2.*grads
                                           +fineGrada*(nref+1.0)
                                           -fineGradb*(nref+3.0))/
                          (nref*nref + 4*nref + 3.0);
                        Real b = (fineGradb-fineGrada)/h - a*h;
                        Real c = fineGrada;
                        Real x = 2.*h;
                        a_tanGradf(ivf,gradComp) = a*x*x + b*x + c;


                        // now average fine-grid 2nd-order gradient with
                        // coarse-grid interpolated gradient
#if 0
                        a_tanGradf(ivf, gradComp) =
                          fineMult*(mult0*fine0 + mult1*fine1 + mult2*fine2)
                          +crseMult*a_tanGradStar(ivf+ihilo*BASISV(a_idir),
                                                  gradComp);
#endif
                      } // end loop over variables
                    } // end loop over forward-difference cells
                } // end if there are forward-difference gradients to compute

              // now do backward-difference cells
              const IntVectSet& backIVS =a_TFSS.getBackwardStencilSet(gradDir);
              if (!backIVS.isEmpty())
                {
                  IVSIterator back_ivsit(backIVS);

                  Real mult0 = 1.5/m_dxFine;
                  Real mult1 = -2.0/m_dxFine;
                  Real mult2 = 0.5/m_dxFine;
                  for (back_ivsit.begin(); back_ivsit.ok(); ++ back_ivsit)
                    {
                      IntVect ivf = back_ivsit();
                      IntVect finePhiLoc = ivf - ihilo*BASISV(a_idir);
                      IntVect finePhiLoc2 = finePhiLoc - ihilo*BASISV(a_idir);
                      // loop over variables
                      for (int var=a_vars.begin(); var<=a_vars.end(); ++var)
                        {
                          Real fine0 = a_phif (finePhiLoc2,var);
                          Real fine1 = a_phif (finePhiLoc2-BASISV(gradDir),var);
                          Real fine2 = a_phif (finePhiLoc2-2*BASISV(gradDir),var);
                          Real fineGrada = mult0*fine0
                            + mult1*fine1 + mult2*fine2;


                          fine0 = a_phif (finePhiLoc,var);
                          fine1 = a_phif (finePhiLoc-BASISV(gradDir),var);
                          fine2 = a_phif (finePhiLoc-2*BASISV(gradDir),var);
                          Real fineGradb = mult0*fine0 +mult1*fine1 +mult2*fine2;


                          int gradComp = gradIndex(var,gradDir);
                          Real grads = a_tanGradStar(ivf+ihilo*BASISV(a_idir),
                                                  gradComp);

                          Real h = m_dxFine;
                          Real a = (2./h/h)*(2.*grads
                                             +fineGrada*(nref+1.0)
                                             -fineGradb*(nref+3.0))/
                            (nref*nref + 4*nref + 3.0);
                          Real b = (fineGradb-fineGrada)/h - a*h;
                          Real c = fineGrada;
                          Real x = 2.*h;
                          a_tanGradf(ivf,gradComp) = a*x*x + b*x + c;

                          // now average fine-grid 2nd-order gradient with
                          // coarse-grid interpolated gradient
#if 0
                          a_tanGradf(ivf, gradComp) =
                            fineMult*(mult0*fine0 + mult1*fine1 + mult2*fine2)
                            +crseMult*a_tanGradStar(ivf+ihilo*BASISV(a_idir),
                                                    gradComp);
#endif

                      } // end loop over variables
                    } // end loop over backward-difference cells
                } // end if there are backward-difference gradients to compute

            } // end if gradDir is a tangential direction
        } // end loop over gradDir


    } // end if there is anything to do here
}

/***********************/
/***********************/
TensorCFInterp::~TensorCFInterp()
{
  clear();
}











#include "NamespaceFooter.H"
