#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "QuadCFStencil.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "Tuple.H"
#include "DataIterator.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

bool
QuadCFStencil::isDefined() const
{
  return m_isDefined;
}

/**
   compute second derivative in devdir direction
   at coarse point ivin
   Uses centered finite diff approximation if
   m_ivsStandard.contains(ivin ).
   Otherwise, it uses the stencil from m_second
   Asserts that derivdir != direction of face
*/
Real
QuadCFStencil::computeSecondDerivative(
                                       const BaseFab<Real> & a_phic,
                                       int a_derivdir,
                                       int a_ivar,
                                       const IntVect& a_ivin,
                                       Real a_dx) const
{
  CH_assert(isDefined());
  CH_assert(a_ivar >= 0);
  CH_assert(a_ivar < a_phic.nComp());
  CH_assert(!a_phic.box().isEmpty());
  CH_assert(a_derivdir >= 0);
  CH_assert(a_derivdir != m_direction);
  CH_assert(a_derivdir < SpaceDim);
  CH_assert(a_dx > 0);
  //use standard centered diff approx if
  //a_ivin is contained in m_ivsStandard
  Real secondd;
  if (m_ivsStandard.contains(a_ivin))
    {
      IntVect basisv = BASISV(a_derivdir);
      IntVect ivcc = a_ivin;
      IntVect ivlo = a_ivin - basisv;
      IntVect ivhi = a_ivin + basisv;
      CH_assert(a_phic.box().contains(ivcc));
      CH_assert(a_phic.box().contains(ivlo));
      CH_assert(a_phic.box().contains(ivhi));
      Real phicc = a_phic(ivcc, a_ivar);
      Real philo = a_phic(ivlo, a_ivar);
      Real phihi = a_phic(ivhi, a_ivar);
      secondd = (phihi + philo - 2.0*phicc)/(a_dx*a_dx);
    }
  else
    {
      CH_assert(m_ivsQuadd.contains(a_ivin));
      CH_assert(m_dropOrd.box().contains(a_ivin));
      if (m_dropOrd(a_ivin))
        {
          secondd = 0.0;
        }
      else
        {
          const BaseFab<DerivStencil>& stenfab
            = m_second[a_derivdir];
          CH_assert(stenfab.box().contains(a_ivin));
          const DerivStencil& dersten = stenfab(a_ivin);
          secondd = 0.0;
          bool keepzer = false;
          for (int isten = 0; isten < dersten.size(); isten++)
            {
              Real wgt = dersten.getWeight(isten);
              IntVect ivphi = dersten.getIndex(isten);
              //can't really assert this because blocking factor
              //could screw it up (grid one cell wide and all that)
              // CH_assert(a_phic.box().contains(ivphi));
              if (a_phic.box().contains(ivphi))
                secondd += wgt*a_phic(ivphi, a_ivar);
              else
                keepzer = true;
            }
          if (keepzer) secondd = 0.0;
          secondd /= (a_dx*a_dx);
        }
    }
  return secondd;
}

/**
   compute first derivative in devdir direction
   at coarse point a_ivin
   Uses centered finite diff approximation if
   m_ivsStandard.contains(a_ivin ).
   Otherwise, it uses the stencil from m_firstD
   Asserts that a_derivdir != direction of face
*/

Real
QuadCFStencil::computeFirstDerivative(
                                      const BaseFab<Real> & a_phic,
                                      int a_derivdir,
                                      int a_ivar,
                                      const IntVect& a_ivin,
                                      Real a_dx) const
{
  CH_assert(isDefined());
  CH_assert(a_ivar >= 0);
  CH_assert(a_ivar < a_phic.nComp());
  CH_assert(!a_phic.box().isEmpty());
  CH_assert(a_derivdir >= 0);
  CH_assert(a_derivdir < SpaceDim);
  CH_assert(a_derivdir != m_direction);
  CH_assert(a_dx > 0);
  //use standard centered diff approx if
  //a_ivin is contained in m_ivsStandard
  Real firstd;
  if (m_ivsStandard.contains(a_ivin))
    {
      IntVect basisv = BASISV(a_derivdir);
      IntVect ivlo = a_ivin - basisv;
      IntVect ivhi = a_ivin + basisv;
      CH_assert(a_phic.box().contains(ivlo));
      CH_assert(a_phic.box().contains(ivhi));
      Real philo = a_phic(ivlo, a_ivar);
      Real phihi = a_phic(ivhi, a_ivar);
      firstd = (phihi - philo)/(2.0*a_dx);
    }
  else
    {
      const BaseFab<DerivStencil>& stenfab
        = m_firstD[a_derivdir];
      CH_assert(stenfab.box().contains(a_ivin));
      const DerivStencil& dersten = stenfab(a_ivin);
      firstd = 0.0;
      bool keepzer = false;
      for (int isten = 0; isten < dersten.size(); isten++)
        {
          Real wgt = dersten.getWeight(isten);
          const IntVect& ivphi= dersten.getIndex(isten);
          //can't really assert this because blocking factor
          //could screw it up (grid one cell wide and all that)
          // CH_assert(a_phic.box().contains(ivphi));
          if (a_phic.box().contains(ivphi))
            firstd += wgt*a_phic(ivphi, a_ivar);
          else
            keepzer = true;
        }
      if (keepzer) firstd = 0.0;
      firstd /= (a_dx);
    }
  return firstd;
}

/**
   compute mixed derivative (direction is unambiguous.
   x and y are the two directions tangential to face)
   at coarse point a_ivin. \\
   In two dimensions, always returns 0. \\
   Uses centered finite diff approximation if
   m_ivsStandard.contains(a_ivin ).
   It uses the stencil from m_mixedSten.\\
   Returns 0 if SpaceDim != 3
*/
Real
QuadCFStencil::computeMixedDerivative(
                                      const BaseFab<Real> & a_phic,
                                      int a_ivar,
                                      const IntVect& a_ivin,
                                      Real a_dx) const
{
#if (CH_SPACEDIM == 3)
  CH_assert(isDefined());
  CH_assert(a_ivar >= 0);
  CH_assert(a_ivar < a_phic.nComp());
  CH_assert(!a_phic.box().isEmpty());
  CH_assert(a_dx > 0);
  //compute derivative directions
  //(defined as normal to m_direction)
  int itran1,itran2;
  if (m_direction==0)
    {
      itran2 =1;
      itran1 =2;
    }
  else if (m_direction==1)
    {
      itran2 =0;
      itran1 =2;
    }
  else if (m_direction==2)
    {
      itran2 =0;
      itran1 =1;
    }
  else
    {
      std::cerr << "Quadcfstencil::define-bogus m_direction" << endl;
      abort();
    }

  //use standard centered diff approx if
  //a_ivin is contained in m_ivsStandard
  Real mixedd;
  if (m_ivsStandard.contains(a_ivin))
    {
      IntVect basex = BASISV(itran2);
      IntVect basey = BASISV(itran1);
      IntVect ivur = a_ivin + basex + basey;
      IntVect ivul = a_ivin - basex + basey;
      IntVect ivlr = a_ivin + basex - basey;
      IntVect ivll = a_ivin - basex - basey;
      CH_assert(a_phic.box().contains(ivur));
      CH_assert(a_phic.box().contains(ivul));
      CH_assert(a_phic.box().contains(ivlr));
      CH_assert(a_phic.box().contains(ivll));
      Real phiur = a_phic(ivur, a_ivar);
      Real phiul = a_phic(ivul, a_ivar);
      Real philr = a_phic(ivlr, a_ivar);
      Real phill = a_phic(ivll, a_ivar);
      mixedd = (phiur+phill-philr-phiul)/(4.*a_dx*a_dx);
    }
  else
    {
      CH_assert(m_ivsQuadd.contains(a_ivin));
      CH_assert(m_dropOrd.box().contains(a_ivin));
      if (m_dropOrd(a_ivin))
        {
          mixedd = 0.0;
        }
      else
        {
          CH_assert(m_mixedSten.box().contains(a_ivin));
          const DerivStencil& dersten = m_mixedSten(a_ivin);
          mixedd = 0.0;
          bool keepzer = false;
          for (int isten = 0; isten < dersten.size(); isten++)
            {
              Real wgt = dersten.getWeight(isten);
              const IntVect& ivphi =  dersten.getIndex(isten);
              //can't really assert this because blocking factor
              //could screw it up (grid one cell wide and all that)
              // CH_assert(a_phic.box().contains(ivphi));
              if (a_phic.box().contains(ivphi))
                mixedd += wgt*a_phic(ivphi, a_ivar);
              else
                keepzer = true;
            }
          if (keepzer) mixedd = 0.0;
          mixedd /= (a_dx*a_dx);
        }
    }
  return mixedd;
#else
  return 0.0;
#endif
}

void
QuadCFStencil::setDefaultValues()
{
  m_isDefined = false;
  m_direction = -1;
}

QuadCFStencil::QuadCFStencil()
{
  setDefaultValues();
}

QuadCFStencil::~QuadCFStencil()
{

}

QuadCFStencil::QuadCFStencil(
                             const Box& a_fine_domain,
                             const Box& a_grid,
                             const DisjointBoxLayout& a_fineBoxes,
                             const DisjointBoxLayout& a_coarBoxes,
                             int a_refRatio, int a_direction,
                             Side::LoHiSide a_hiorlo)
{
  setDefaultValues();
  ProblemDomain fineProbDomain(a_fine_domain);
  define(fineProbDomain,  a_grid, a_fineBoxes, a_coarBoxes,
         a_refRatio,  a_direction, a_hiorlo);
}

QuadCFStencil::QuadCFStencil(
                             const ProblemDomain& a_fine_domain,
                             const Box& a_grid,
                             const DisjointBoxLayout& a_fineBoxes,
                             const DisjointBoxLayout& a_coarBoxes,
                             int a_refRatio, int a_direction,
                             Side::LoHiSide a_hiorlo)
{
  setDefaultValues();
  define(a_fine_domain, a_grid, a_fineBoxes, a_coarBoxes,
         a_refRatio,  a_direction, a_hiorlo);
}

void
QuadCFStencil::define(
                      const Box& a_fine_domain,
                      const Box& a_grid,
                      const DisjointBoxLayout& a_fineBoxes,
                      const DisjointBoxLayout& a_coarBoxes,
                      int a_refRatio,
                      int a_direction,
                      Side::LoHiSide a_hiorlo)
{
  ProblemDomain fineProbDomain(a_fine_domain);
  define(fineProbDomain, a_grid, a_fineBoxes, a_coarBoxes,
         a_refRatio, a_direction, a_hiorlo);

}

void
QuadCFStencil::define(
                      const ProblemDomain& a_fine_domain,
                      const Box& a_grid,
                      const DisjointBoxLayout& a_fineBoxes,
                      const DisjointBoxLayout& a_coarBoxes,
                      int a_refRatio,
                      int a_direction,
                      Side::LoHiSide a_hiorlo)
{

  CH_assert(!a_fine_domain.isEmpty());

  CH_assert(a_refRatio >= 1);
  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert(!a_grid.isEmpty());
  CH_assert(a_fine_domain.contains(a_grid));
  CH_assert((a_hiorlo == Side::Lo) ||
         (a_hiorlo == Side::Hi));
  m_isDefined = true;

  m_direction = a_direction;

  //define base cfstencil
  m_baseCFS.define(a_fine_domain,
                   a_grid,
                   a_fineBoxes, a_coarBoxes,
                   a_refRatio, a_direction,
                   a_hiorlo);

  const IntVectSet& baseIVS = m_baseCFS.getCoarIVS();

  if (baseIVS.isEmpty()) return;

  //these are the directions tangential to the face
#if CH_SPACEDIM == 1
  // this shouldn't ever be used, but it's simpler
  // to keep it around rather than block out every
  // place it exists...
  Tuple<int,1> vinttran;
#else
  Tuple<int,SpaceDim-1> vinttran;
#endif

  int iint = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != m_direction)
        {
          vinttran[iint] = idir;
          iint++;
        }
    }

  //compute all cells necessary for both stencils
  //and interpolation points
  Box coarseGrid = coarsen(a_grid, a_refRatio);
  ProblemDomain coarseDomain = coarsen(a_fine_domain, a_refRatio);

  // create a box which will indicate whether a cell is adjacent
  // to a periodic boundary or not (all cells contained in the box
  // are _not_ adjacent to a periodic boundary)
  Box periodicTestBox(coarseDomain.domainBox());
  if (coarseDomain.isPeriodic())
    for (int idir=0; idir<SpaceDim; idir++)
      if (coarseDomain.isPeriodic(idir))
        periodicTestBox.grow(idir,-1);

  Box edgeBoxG1, edgeBoxG2;
  if (a_hiorlo == Side::Lo)
    {
      edgeBoxG2 = adjCellLo(coarseGrid,m_direction,1);
      edgeBoxG1 = adjCellLo(coarseGrid,m_direction,1);
    }
  else
    {
      edgeBoxG2 = adjCellHi(coarseGrid,m_direction,1);
      edgeBoxG1 = adjCellHi(coarseGrid,m_direction,1);
    }

  for (int iint = 0; iint < SpaceDim-1; iint++)
    {
      int itran = vinttran[iint];
      edgeBoxG2.grow(itran, 2);
      edgeBoxG1.grow(itran, 1);
    }
  edgeBoxG2 &= coarseDomain;
  edgeBoxG1 &= coarseDomain;
  //this is the coarsened version of the fine union
  //of rectangles
  DisjointBoxLayout coarsenedFineBoxes;
  coarsen(coarsenedFineBoxes, a_fineBoxes, a_refRatio);
  IntVectSet ivsAllGood(edgeBoxG2);
  LayoutIterator litFine = coarsenedFineBoxes.layoutIterator();
  for (litFine.reset(); litFine.ok(); ++litFine)
    {
      ivsAllGood -= coarsenedFineBoxes.get(litFine());
      // need to subtract periodic images as well
      if (coarseDomain.isPeriodic()
          && !periodicTestBox.contains(coarsenedFineBoxes.get(litFine()))
          && !periodicTestBox.contains(coarseGrid))
        {
          ShiftIterator shiftIt = coarseDomain.shiftIterator();
          IntVect shiftMult = coarseDomain.domainBox().size();
          Box shiftedBox = coarsenedFineBoxes.get(litFine());
          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
            {
              IntVect shiftVect = shiftMult*shiftIt();
              shiftedBox.shift(shiftVect);
              ivsAllGood -= shiftedBox;
              shiftedBox.shift(-shiftVect);
            } // end loop over periodic shift directions
        } // end if is periodic
    }
  //create temp IVS which has all points in stencils
  //which fell off coarse grid and subtract those
  IntVectSet ivsEdgeOfEarth = ivsAllGood;
  LayoutIterator litCoar =a_coarBoxes.layoutIterator();
  for (litCoar.reset(); litCoar.ok(); ++litCoar)
    {
      ivsEdgeOfEarth -= a_coarBoxes.get(litCoar());
      // if is periodic, need to subtract off periodic images here as well
      if (coarseDomain.isPeriodic()
          && !periodicTestBox.contains(a_coarBoxes.get(litCoar()))
          && !periodicTestBox.contains(coarseGrid))
        {
          ShiftIterator shiftIt = coarseDomain.shiftIterator();
          IntVect shiftMult = coarseDomain.domainBox().size();
          Box shiftedBox = a_coarBoxes.get(litCoar());
          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
            {
              IntVect shiftVect = shiftMult*shiftIt();
              shiftedBox.shift(shiftVect);
              ivsEdgeOfEarth -= shiftedBox;
              shiftedBox.shift(-shiftVect);
            } // end loop over periodic shift directions
        } // end if is periodic

    }
  ivsAllGood -= ivsEdgeOfEarth;

  //compute cells where pure centered diffs can be taken
  //start with the set of all coarse intvects
  //on the coarse-fine interface and compute
  //the set of all cells for which all the
  //centered-difference stencil is within
  //the coarse grid == m_ivsStandard
  m_ivsStandard = ivsAllGood;
  //ivsStandard = grow((ivsAllGood & edgeBoxG1), itran,-1)
  m_ivsStandard &= edgeBoxG1;

  for (int iint = 0; iint < SpaceDim-1; iint++)
    {
      int itran = vinttran[iint];
      m_ivsStandard.grow(itran, -1);
    }
  //the points where something special =
  //all points- standard points
  //by conservation
  m_ivsQuadd = baseIVS;
  m_ivsQuadd -= m_ivsStandard;

  if (!m_ivsQuadd.isEmpty())
    {
      buildStencils(ivsAllGood);
    }
}

void
QuadCFStencil::define(
                      const ProblemDomain& a_fine_domain,
                      const Box& a_grid,
                      const Vector<Box>& a_periodicBoxes,
                      const Vector<Box>& a_coarsenedPeriodicBoxes,
                      const DisjointBoxLayout& a_coarBoxes,
                      int a_refRatio,
                      int a_direction,
                      Side::LoHiSide a_hiorlo)
{

  CH_assert(!a_fine_domain.isEmpty());

  CH_assert(a_refRatio >= 1);
  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert(!a_grid.isEmpty());
  CH_assert(a_fine_domain.contains(a_grid));
  CH_assert((a_hiorlo == Side::Lo) ||
         (a_hiorlo == Side::Hi));
  m_isDefined = true;

  m_direction = a_direction;

  //define base cfstencil
  m_baseCFS.define(a_fine_domain,
                   a_grid,
                   a_periodicBoxes,
                   a_refRatio, a_direction,
                   a_hiorlo);

  const IntVectSet& baseIVS = m_baseCFS.getCoarIVS();

  if (baseIVS.isEmpty()) return;

  //these are the directions tangential to the face
#if CH_SPACEDIM == 1
  // this shouldn't ever be used, but it's simpler
  // to keep it around rather than block out every
  // place it exists...
  Tuple<int,1> vinttran;
#else
  Tuple<int,SpaceDim-1> vinttran;
#endif

  int iint = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != m_direction)
        {
          vinttran[iint] = idir;
          iint++;
        }
    }

  //compute all cells necessary for both stencils
  //and interpolation points
  Box coarseGrid = coarsen(a_grid, a_refRatio);
  ProblemDomain coarseDomain = coarsen(a_fine_domain, a_refRatio);

  // create a box which will indicate whether a cell is adjacent
  // to a periodic boundary or not (all cells contained in the box
  // are _not_ adjacent to a periodic boundary)
  Box periodicTestBox(coarseDomain.domainBox());
  if (coarseDomain.isPeriodic())
    for (int idir=0; idir<SpaceDim; idir++)
      if (coarseDomain.isPeriodic(idir))
        periodicTestBox.grow(idir,-1);

  Box edgeBoxG1, edgeBoxG2;
  if (a_hiorlo == Side::Lo)
    {
      edgeBoxG2 = adjCellLo(coarseGrid,m_direction,1);
      edgeBoxG1 = adjCellLo(coarseGrid,m_direction,1);
    }
  else
    {
      edgeBoxG2 = adjCellHi(coarseGrid,m_direction,1);
      edgeBoxG1 = adjCellHi(coarseGrid,m_direction,1);
    }

  for (int iint = 0; iint < SpaceDim-1; iint++)
    {
      int itran = vinttran[iint];
      edgeBoxG2.grow(itran, 2);
      edgeBoxG1.grow(itran, 1);
    }
  edgeBoxG2 &= coarseDomain;
  edgeBoxG1 &= coarseDomain;
  //this is the coarsened version of the fine union
  //of rectangles
  if (edgeBoxG2.isEmpty()) return;

  IntVectSet ivsAllGood(edgeBoxG2);

  int w1 = edgeBoxG2.smallEnd()[0];
  int w2 = edgeBoxG2.bigEnd()[0];

  for (int i=0; i<a_coarsenedPeriodicBoxes.size(); ++i)
    {
      const Box& box = a_coarsenedPeriodicBoxes[i];
      if (box.bigEnd()[0] >= w1)
        {
          ivsAllGood -= box;
          if (box.smallEnd()[0] > w2) i=a_coarsenedPeriodicBoxes.size();
        }
    }

  //create temp IVS which has all points in stencils
  //which fell off coarse grid and subtract those
  IntVectSet ivsEdgeOfEarth = ivsAllGood;
  LayoutIterator litCoar =a_coarBoxes.layoutIterator();
  for (litCoar.reset(); litCoar.ok(); ++litCoar)
    {
      ivsEdgeOfEarth -= a_coarBoxes.get(litCoar());
      // if is periodic, need to subtract off periodic images here as well
      if (coarseDomain.isPeriodic()
          && !periodicTestBox.contains(a_coarBoxes.get(litCoar()))
          && !periodicTestBox.contains(coarseGrid))
        {
          ShiftIterator shiftIt = coarseDomain.shiftIterator();
          IntVect shiftMult = coarseDomain.domainBox().size();
          Box shiftedBox = a_coarBoxes.get(litCoar());
          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
            {
              IntVect shiftVect = shiftMult*shiftIt();
              shiftedBox.shift(shiftVect);
              ivsEdgeOfEarth -= shiftedBox;
              shiftedBox.shift(-shiftVect);
            } // end loop over periodic shift directions
        } // end if is periodic

    }
  ivsAllGood -= ivsEdgeOfEarth;

  //compute cells where pure centered diffs can be taken
  //start with the set of all coarse intvects
  //on the coarse-fine interface and compute
  //the set of all cells for which all the
  //centered-difference stencil is within
  //the coarse grid == m_ivsStandard
  m_ivsStandard = ivsAllGood;
  //ivsStandard = grow((ivsAllGood & edgeBoxG1), itran,-1)
  m_ivsStandard &= edgeBoxG1;

  for (int iint = 0; iint < SpaceDim-1; iint++)
    {
      int itran = vinttran[iint];
      m_ivsStandard.grow(itran, -1);
    }
  //the points where something special =
  //all points- standard points
  //by conservation
  m_ivsQuadd = baseIVS;
  m_ivsQuadd -= m_ivsStandard;

  if (!m_ivsQuadd.isEmpty())
    {
      buildStencils(ivsAllGood);
    }
}

void QuadCFStencil::buildStencils(const IntVectSet& ivsAllGood)
{

 #if CH_SPACEDIM == 1
  // this shouldn't ever be used, but it's simpler
  // to keep it around rather than block out every
  // place it exists...
  Tuple<int,1> vinttran;
#else
  Tuple<int,SpaceDim-1> vinttran;
#endif

  int iint = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != m_direction)
        {
          vinttran[iint] = idir;
          iint++;
        }
    }
  const Box& osbox = m_ivsQuadd.minBox();
  m_dropOrd.resize(osbox);
  m_dropOrd.setVal(false);

#if (CH_SPACEDIM == 3)
  //first calculate mixed derivative stencils
  //only in the case of 3 dimensions

  m_mixedSten.resize(osbox);
  //create four boxes of stencil
  int itran1 = vinttran[0];
  int itran2 = vinttran[1];
  Box stenbox(IntVect::Zero, BASISV(itran1)+BASISV(itran2));
  //fabs for four corners of stencil
  //0 shifts shown as comments for pedagogical reasons

  BaseFab<Real>  fabur(stenbox, 1);
  fabur.setVal(-1.0);
  fabur(BASISV(itran1)) = 1.0;
  fabur(BASISV(itran2)) = 1.0;
  //fabur.shift(itran1, 0);
  //fabur.shift(itran1, 0);

  BaseFab<Real>  fabll(stenbox, 1);
  fabll.copy(fabur);
  fabll.shift(itran1, -1);
  fabll.shift(itran2, -1);

  BaseFab<Real>  fablr(stenbox, 1);
  fablr.copy(fabur);
  //fablr.shift(itran1, 0);
  fablr.shift(itran2, -1);

  BaseFab<Real>  fabul(stenbox, 1);
  fabul.copy(fabur);
  fabur.shift(itran1, -1);
  //fabll.shift(itran2, 0);
  IVSIterator irregit(m_ivsQuadd);
  for (irregit.begin(); irregit.ok(); ++irregit)
    {
      IntVect ivirr = irregit();
      fabur.shift(ivirr);
      fablr.shift(ivirr);
      fabul.shift(ivirr);
      fabll.shift(ivirr);
      Box bur(fabur.box());
      Box bul(fabul.box());
      Box blr(fablr.box());
      Box bll(fabll.box());
      DerivStencil& dersten= m_mixedSten(ivirr);
      dersten.define();
      bool hasamixed = false;
      int imixed = 0;
      //add each stencil whose box is
      //contained in ivsAllGood and
      //average the result
      if (ivsAllGood.contains(bur))
        {
          hasamixed = true;
          addFabToSten(fabur, dersten);
          imixed++;
        }
      if (ivsAllGood.contains(bul))
        {
          hasamixed = true;
          addFabToSten(fabul, dersten);
          imixed++;
        }
      if (ivsAllGood.contains(blr))
        {
          hasamixed = true;
          addFabToSten(fablr, dersten);
          imixed++;
        }
      if (ivsAllGood.contains(bll))
        {
          hasamixed = true;
          addFabToSten(fabll, dersten);
          imixed++;
        }
      if (imixed > 0)
        {
          Real rnmixed = imixed;
          dersten /= rnmixed;
        }
      m_dropOrd(ivirr) = (!hasamixed);
      //shift fabs back so that stencils work
      fabur.shift(-ivirr);
      fablr.shift(-ivirr);
      fabul.shift(-ivirr);
      fabll.shift(-ivirr);
    }
#endif
  //now calculate all non-mixed derivative stencils
  for (int ivec = 0; ivec < SpaceDim-1; ivec++)
    {
      int idir = vinttran[ivec];
      IntVect basisv = BASISV(idir);
      BaseFab<DerivStencil>& second_ds = m_second [idir];
      BaseFab<DerivStencil>& firstd_ds = m_firstD[idir];
      second_ds.resize(osbox);
      firstd_ds.resize(osbox);
      IVSIterator irregit(m_ivsQuadd);
      for (irregit.begin(); irregit.ok(); ++irregit)
        {
          IntVect iv = irregit();
          if (!m_dropOrd(iv))
            {
              Box bsten(iv, iv);
              bsten.grow(idir, 1);
              if (ivsAllGood.contains(bsten))
                {
                  //here centered diffs in this direction
                  //apply.  this should never get hit in
                  //two dimensions.  In 2d, if iv is part of
                  //irregit, then it should need one
                  //sided stuff
                  CH_assert(SpaceDim == 3);
                  BaseFab<Real>  fabsten(bsten,1);
                  fabsten(iv-basisv) =  1.0;
                  fabsten(iv+basisv) =  1.0;
                  fabsten(iv)        = -2.0;
                  DerivStencil& secdersten= second_ds(iv);
                  secdersten.define();
                  addFabToSten(fabsten, secdersten);
                  fabsten(iv-basisv) = -0.5;
                  fabsten(iv+basisv) =  0.5;
                  fabsten(iv)        =  0.0;
                  DerivStencil& firdersten= firstd_ds(iv);
                  firdersten.define();
                  addFabToSten(fabsten, firdersten);
                }
              else
                {
                  Box bshift(bsten);
                  bshift.shift(basisv);
                  if (ivsAllGood.contains(bshift))
                    {
                      //point on the negative side is covered but
                      //the other one is not so introduce shifted
                      //stencil
                      BaseFab<Real>  fabsten(bshift,1);
                      fabsten(iv)          =  1.0;
                      fabsten(iv+  basisv) = -2.0;
                      fabsten(iv+2*basisv) =  1.0;
                      DerivStencil& secdersten= second_ds(iv);
                      secdersten.define();
                      addFabToSten(fabsten, secdersten);
                      //need three point stencil for derivative
                      //to get second-order, shifted
                      fabsten(iv)          = -3.0/2.0;
                      fabsten(iv+  basisv) =  4.0/2.0;
                      fabsten(iv+2*basisv) = -1.0/2.0;
                      DerivStencil& firdersten= firstd_ds(iv);
                      firdersten.define();
                      addFabToSten(fabsten, firdersten);
                    }
                  else
                    {
                      bshift.shift(-2*basisv);
                      if (ivsAllGood.contains(bshift))
                        {
                          //point on the positive side is covered but
                          //the other one is not so introduce shifted
                          //stencil
                          BaseFab<Real>  fabsten(bshift,1);
                          fabsten(iv)          =  1.0;
                          fabsten(iv-  basisv) = -2.0;
                          fabsten(iv-2*basisv) =  1.0;
                          DerivStencil& secdersten= second_ds(iv);
                          secdersten.define();
                          addFabToSten(fabsten, secdersten);
                          //need three point stencil for derivative
                          //to get second-order, shifted
                          fabsten(iv)          =  3.0/2.0;
                          fabsten(iv-  basisv) = -4.0/2.0;
                          fabsten(iv-2*basisv) =  1.0/2.0;
                          DerivStencil& firdersten= firstd_ds(iv);
                          firdersten.define();
                          addFabToSten(fabsten, firdersten);
                        }
                      else
                        {
                          //both sides have a covered point.
                          //need to drop order and use
                          //one-sided first derivative
                          m_dropOrd(iv) = true;
                          IntVect ivhi = iv + basisv;
                          if (ivsAllGood.contains(ivhi))
                            {
                              //derivative taken on positive side
                              Box bfirst(iv, ivhi);
                              BaseFab<Real>  fabsten(bfirst, 1);
                              fabsten(iv  ) = -1.0;
                              fabsten(ivhi) =  1.0;
                              DerivStencil& firdersten= firstd_ds(iv);
                              firdersten.define();
                              addFabToSten(fabsten, firdersten);
                            }
                          else
                            {
                              //derivative taken on negative side
                              IntVect ivlo = iv - basisv;
                              if (ivsAllGood.contains(ivlo))
                                {
                                  Box bfirst(ivlo, iv);
                                  BaseFab<Real>  fabsten(bfirst, 1);
                                  fabsten(iv  ) =  1.0;
                                  fabsten(ivlo) = -1.0;
                                  DerivStencil& firdersten= firstd_ds(iv);
                                  firdersten.define();
                                  addFabToSten(fabsten, firdersten);
                                }
                              else
                                {
                                  //neither is there.  set fab sten to zero.
                                  Box bfirst(iv, iv);
                                  BaseFab<Real>  fabsten(bfirst, 1);
                                  fabsten(iv  ) =  0.0;
                                  DerivStencil& firdersten= firstd_ds(iv);
                                  firdersten.define();
                                  addFabToSten(fabsten, firdersten);
                                } //end error message
                            }//end positive point covered
                        }//end drop order else
                    }//end point on positive side covered else
                }//end if (there needs one side in THIS direction)
            } //end if (!m_dropOrd) (if mixed is calculable)
        } //end loop over one sided cells
    }//end loop over tangential directions
}

/**
   add fab to the stencil.
   the weight equal to the value of the
   of the fab data at the box
*/
void
QuadCFStencil::addFabToSten(
                            const BaseFab<Real> & a_fabin,
                            DerivStencil&  a_sten
                            )
{
  CH_assert(isDefined());
  CH_assert(!a_fabin.box().isEmpty());
  CH_assert(a_fabin.nComp() == 1);
  BoxIterator bit(a_fabin.box());
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      Real weight = a_fabin(iv);
      a_sten.accumulate(iv, weight);
    }
}
#include "NamespaceFooter.H"
