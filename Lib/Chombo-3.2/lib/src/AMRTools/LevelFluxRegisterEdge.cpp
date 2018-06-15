#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelFluxRegisterEdge.H"
#include "LayoutIterator.H"
#include "DebugOut.H"
#include "Copier.H"
#include "parstream.H"

#include "NamespaceHeader.H"

void
LevelFluxRegisterEdge::setDefaultValues()
{
  m_isDefined = false;
  m_nComp = -1;
  m_nRefine = -1;
}

int LevelFluxRegisterEdge::index(int dir, Side::LoHiSide side)
{
  return side*SpaceDim+dir;
}

int LevelFluxRegisterEdge::getRegComp(const int& faceDir,
                                      const int& edgeDir)
{
  int regcomp;
  // need to think of a better way to do this.
  // regcomp is the remaining dimension after discounting the
  // two arguments
  CH_assert (SpaceDim == 3);
  if (faceDir == 0)
  {
    regcomp = 1;
    if (edgeDir == 1)
    {
      regcomp = 2;
    }
  }
  else if (faceDir == 1)
  {
    regcomp = 0;
    if (edgeDir == 0)
    {
      regcomp = 2;
    }
  }
  else if (faceDir == 2)
  {
    regcomp = 0;
    if (edgeDir == 0)
    {
      regcomp = 1;
    }
  }

  return regcomp;
}



LevelFluxRegisterEdge::LevelFluxRegisterEdge(const DisjointBoxLayout& a_dblFine,
                                             const DisjointBoxLayout& a_dblCoar,
                                             const Box& a_dProblem,
                                             int a_nRefine,
                                             int a_nComp)
{
  setDefaultValues();
  ProblemDomain physdomain(a_dProblem);
  define(a_dblFine, a_dblCoar, physdomain, a_nRefine, a_nComp);
}

LevelFluxRegisterEdge::LevelFluxRegisterEdge(const DisjointBoxLayout& a_dblFine,
                                             const DisjointBoxLayout& a_dblCoar,
                                             const ProblemDomain& a_dProblem,
                                             int a_nRefine,
                                             int a_nComp)
{
  setDefaultValues();
  define(a_dblFine, a_dblCoar, a_dProblem, a_nRefine, a_nComp);
}





void
LevelFluxRegisterEdge::define(
                              const DisjointBoxLayout& a_dbl,
                              const DisjointBoxLayout& a_dblCoarse,
                              const Box& a_dProblem,
                              int a_nRefine,
                              int a_nComp)
{
  ProblemDomain physdomain(a_dProblem);
  define(a_dbl, a_dblCoarse, physdomain, a_nRefine, a_nComp);
}


// new define
void
LevelFluxRegisterEdge::define(
                              const DisjointBoxLayout& a_dbl,
                              const DisjointBoxLayout& a_dblCoarse,
                              const ProblemDomain& a_dProblem,
                              int a_nRefine,
                              int a_nComp)
{
  m_isDefined = true;
  CH_assert(a_nRefine > 0);
  CH_assert(a_nComp > 0);
  CH_assert(!a_dProblem.isEmpty());
  m_nComp = a_nComp;
  m_nRefine = a_nRefine;
  m_domainCoarse = coarsen(a_dProblem, a_nRefine);
  CH_assert (a_dblCoarse.checkPeriodic(m_domainCoarse));

  // allocate copiers
  m_crseCopiers.resize(SpaceDim*2);

  SideIterator side;

  // create a Vector<Box> of the fine boxes which also includes periodic images,
  // since we don't really care about the processor layouts, etc
  Vector<Box> periodicFineBoxes;
  CFStencil::buildPeriodicVector(periodicFineBoxes, a_dProblem, a_dbl);
  // now coarsen these boxes...
  for (int i=0; i<periodicFineBoxes.size(); i++)
    {
      periodicFineBoxes[i].coarsen(m_nRefine);
    }

  for (int idir=0 ; idir<SpaceDim; ++idir)
  {
    for (side.begin(); side.ok(); ++side)
    {
      // step one, build fineBoxes, flux register boxes
      // indexed by the fine level but in the coarse index
      // space
      DisjointBoxLayout fineBoxes,tmp;
      // first create coarsened dbl, then compute flux register boxes
      // adjacent to coarsened fine boxes
      coarsen(tmp, a_dbl, m_nRefine);
      if (side() == Side::Lo)
        {
          adjCellLo(fineBoxes, tmp, idir,1);
        }
      else
        {
          adjCellHi(fineBoxes, tmp, idir,1);
        }

      // now define the FluxBoxes of fabFine on this DisjointBoxLayout
      m_fabFine[index(idir, side())].define(fineBoxes, a_nComp);



      LayoutData<Vector<Vector<IntVectSet> > >& ivsetsVect
        = m_refluxLocations[index(idir, side())];
      ivsetsVect.define(a_dblCoarse);

      LayoutData<Vector<DataIndex> >& mapsV =
        m_coarToCoarMap[index(idir, side())];
      mapsV.define(a_dblCoarse);

      DisjointBoxLayout coarseBoxes = a_dblCoarse;
      DataIterator dit = a_dblCoarse.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          unsigned int thisproc = a_dblCoarse.procID(dit());
          if (thisproc == procID())
          {
            ivsetsVect[DataIndex(dit())].resize(SpaceDim);
          }
          const Box& coarseBox = a_dblCoarse[dit()];
          int count = 0;
          for (int i=0; i<periodicFineBoxes.size(); i++)
            {
              Box regBox;
              if (side() == Side::Lo)
                {
                  regBox = adjCellLo(periodicFineBoxes[i], idir, 1);
                }
              else
                {
                  regBox = adjCellHi(periodicFineBoxes[i], idir, 1);
                }

              // do this little dance in order to ensure that
              // we catch corner cells which might be in different
              // boxes.
              Box testBox(regBox);
              testBox.grow(1);
              testBox.grow(idir,-1);
              if (testBox.intersectsNotEmpty(coarseBox))
                {
                  testBox &= coarseBox;
                  ++count;
                  unsigned int proc = a_dblCoarse.procID(dit());
                  const DataIndex index = DataIndex(dit());
                  if (proc == procID())
                  {
                    mapsV[DataIndex(dit())].push_back(index);
                    // loop over face directions here
                    for (int faceDir=0; faceDir<SpaceDim; faceDir++)
                    {
                      // do nothing in normal direction
                      if (faceDir != idir)
                      {
                        // this should give us the face indices for the
                        // faceDir-centered faces adjacent to the coarse-fine
                        // interface which are contained in the current
                        // coarse box
                        Box intersectBox(regBox);
                        Box coarseEdgeBox(coarseBox);
                        coarseEdgeBox.surroundingNodes(faceDir);
                        intersectBox.surroundingNodes(faceDir);
                        intersectBox &= coarseEdgeBox;
                        intersectBox.shiftHalf(faceDir,1);
                        IntVectSet localIV(intersectBox);
                        ivsetsVect[DataIndex(dit())][faceDir].push_back(localIV);
                      }
                    }
                  }
                }
            } // end loop over boxes on coarse level
        }
      m_regCoarse.define(coarseBoxes, a_nComp, IntVect::Unit);

      // last thing to do is to define copiers
      m_crseCopiers[index(idir, side())].define(fineBoxes, coarseBoxes,
                                                IntVect::Unit);
    }
  }
}

LevelFluxRegisterEdge::LevelFluxRegisterEdge()
{
  setDefaultValues();
}

LevelFluxRegisterEdge::~LevelFluxRegisterEdge()
{
}

bool
LevelFluxRegisterEdge::isDefined() const
{
  return m_isDefined;
}

void
LevelFluxRegisterEdge::undefine()
{
  m_isDefined = false;
}

void
LevelFluxRegisterEdge::setToZero()
{

  for (DataIterator dit = m_regCoarse.dataIterator(); dit.ok(); ++dit)
    m_regCoarse[dit()].setVal(0.0);


  SideIterator side;
  for (int idir=0 ; idir<SpaceDim; ++idir)
  {
    for (side.begin(); side.ok(); ++side)
    {
      LevelData<FluxBox>& fineReg = m_fabFine[index(idir, side())];


      for (DataIterator dit = fineReg.dataIterator(); dit.ok(); ++dit)
        fineReg[dit()].setVal(0.0);

    }
  }
}

void
LevelFluxRegisterEdge::incrementCoarse(FArrayBox& a_coarseFlux,
                                       Real a_scale,
                                       const DataIndex& a_coarseDataIndex,
                                       const Interval& a_srcInterval,
                                       const Interval& a_dstInterval)
{
  CH_assert(isDefined());
  CH_assert(!a_coarseFlux.box().isEmpty());

  CH_assert(a_srcInterval.size() == a_dstInterval.size());
  CH_assert(a_srcInterval.begin() >= 0);
  CH_assert(a_srcInterval.end() < a_coarseFlux.nComp());
  CH_assert(a_dstInterval.begin() >= 0);
  CH_assert(a_dstInterval.end() < m_nComp);

  // get edge-centering of coarseFlux
  const Box& edgeBox = a_coarseFlux.box();
  int edgeDir = -1;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (edgeBox.type(dir) == IndexType::CELL)
        {
          if (edgeDir == -1)
            {
              edgeDir = dir;
            }
          else
            {
              // already found a cell-centered direction (should only be
              // one for edge-centering)
              MayDay::Error("LevelFluxRegisterEdge::incrementCoarse -- e-field not edge-centered");
            }
        }
    } // end loop over directions
  CH_assert(edgeDir != -1);

  FArrayBox& thisCrseReg = m_regCoarse[a_coarseDataIndex][edgeDir];

  thisCrseReg.plus(a_coarseFlux, -a_scale, a_srcInterval.begin(),
                   a_dstInterval.begin(), a_srcInterval.size());

}

void
LevelFluxRegisterEdge::incrementFine(
                                     FArrayBox& a_fineFlux,
                                     Real a_scale,
                                     const DataIndex& a_fineDataIndex,
                                     const Interval& a_srcInterval,
                                     const Interval& a_dstInterval)
{
  CH_assert(isDefined());
  CH_assert(!a_fineFlux.box().isEmpty());
  CH_assert(a_srcInterval.size() == a_dstInterval.size());
  CH_assert(a_srcInterval.begin() >= 0);
  CH_assert(a_srcInterval.end() < a_fineFlux.nComp());
  CH_assert(a_dstInterval.begin() >= 0);
  CH_assert(a_dstInterval.end() < m_nComp);

  int edgeDir = -1;
  for (int sideDir = 0; sideDir<SpaceDim; sideDir++)
    {
      if (a_fineFlux.box().type(sideDir) == IndexType::CELL)
        {
          edgeDir = sideDir;
        }
    }
  CH_assert(edgeDir >= 0);
  CH_assert(edgeDir < SpaceDim);

  for (int faceDir=0; faceDir<SpaceDim; faceDir++)
    {
      if (faceDir != edgeDir)
        {

          SideIterator sit;
          for (sit.begin(); sit.ok(); ++sit)
            {
              incrementFine(a_fineFlux,
                            a_scale,
                            a_fineDataIndex,
                            a_srcInterval,
                            a_dstInterval,
                            faceDir,
                            sit());
            }
        }
    }
}

void
LevelFluxRegisterEdge::incrementFine(
                                     FArrayBox& a_fineFlux,
                                     Real a_scale,
                                     const DataIndex& a_fineDataIndex,
                                     const Interval& a_srcInterval,
                                     const Interval& a_dstInterval,
                                     int a_dir,
                                     Side::LoHiSide a_sd)
{
  CH_assert(isDefined());
  CH_assert(!a_fineFlux.box().isEmpty());
  CH_assert(a_srcInterval.size() == a_dstInterval.size());
  CH_assert(a_srcInterval.begin() >= 0);
  CH_assert(a_srcInterval.end() < a_fineFlux.nComp());
  CH_assert(a_dstInterval.begin() >= 0);
  CH_assert(a_dstInterval.end() < m_nComp);

  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert((a_sd == Side::Lo)||(a_sd == Side::Hi));
  //
  //
  //denom is the number of fine faces per coarse face
  //this is intrinsically dimension-dependent
#if (CH_SPACEDIM == 2)
  Real denom = 1;
#elif (CH_SPACEDIM == 3)
  Real denom = m_nRefine;
#else
  // This code doesn't make any sense in 1D, and hasn't been implemented
  // for DIM > 3
  Real denom = -1.0;
  MayDay::Error("LevelFluxRegisterEdge -- bad SpaceDim");
#endif

  Real scale = a_scale/denom;

  // need which fluxbox face we're doing this for
  Box thisBox = a_fineFlux.box();
  int fluxComp = -1;
  for (int sideDir=0; sideDir<SpaceDim; sideDir++)
  {
    // we do nothing in the direction normal to face
    if (sideDir != a_dir)
    {
      if (thisBox.type(sideDir) == IndexType::CELL)
      {
        fluxComp = sideDir;
      }
    }
  }
  CH_assert (fluxComp >= 0);
  int regcomp = getRegComp(a_dir, fluxComp);

  FluxBox& thisReg = m_fabFine[index(a_dir, a_sd)][a_fineDataIndex];
  FArrayBox& reg = thisReg[regcomp];

  a_fineFlux.shiftHalf(a_dir, sign(a_sd));

  // this is a way of geting a face-centered domain
  // box which we can then use to intersect with things
  // to screen out cells outside the physical domain
  // (nothing is screened out in periodic case)
  Box shiftedValidDomain = m_domainCoarse.domainBox();
  shiftedValidDomain.grow(2);
  shiftedValidDomain &= m_domainCoarse;
  shiftedValidDomain.surroundingNodes(regcomp);

  BoxIterator regIt(reg.box() & shiftedValidDomain);
  for (regIt.begin(); regIt.ok(); ++regIt)
    {
      const IntVect& coarseIndex = regIt();
      // create a cell-centered box, then shift back to face-centered
      Box box(coarseIndex, coarseIndex);
      box.shiftHalf(regcomp,-1);
      // to avoid adding in edges which do not overlie coarse-grid
      // edges, will refine only in non-fluxComp directions to
      // determine box from which to grab fluxes.
      IntVect refineVect(m_nRefine*IntVect::Unit);
      //refineVect.setVal(fluxComp,1);
      box.refine(refineVect);
      if (a_sd == Side::Lo) box.growLo(a_dir,-(m_nRefine-1));
      else                 box.growHi(a_dir,-(m_nRefine-1));
      BoxIterator fluxIt(box);
      for (fluxIt.begin(); fluxIt.ok(); ++fluxIt)
        {
          int src = a_srcInterval.begin();
          int dest =  a_dstInterval.begin();
          for ( ; src <=a_srcInterval.end(); ++src,++dest)
            reg(coarseIndex, dest) += scale*a_fineFlux(fluxIt(), src);
        }
    }
  a_fineFlux.shiftHalf(a_dir, -sign(a_sd));
}

void
LevelFluxRegisterEdge::refluxCurl(LevelData<FluxBox>& a_uCoarse,
                                  Real a_scale)
{
  CH_assert(isDefined());
  CH_assert(a_uCoarse.nComp() == m_nComp);

  SideIterator side;
  // idir is the normal direction to the coarse-fine interface
  for (int idir=0 ; idir<SpaceDim; ++idir)
  {
    for (side.begin(); side.ok(); ++side)
    {
      LevelData<FluxBox>& fineReg = m_fabFine[index(idir, side())];

      // first, create temp LevelData<FluxBox> to hold "coarse flux"
      const DisjointBoxLayout coarseBoxes = m_regCoarse.getBoxes();

      // this fills the place of what used to be m_fabCoarse in the old
      // implementation
      LevelData<FluxBox> coarReg(coarseBoxes, m_nComp, IntVect::Unit);


      // now fill the coarReg with the curl of the stored coarse-level
      // edge-centered flux
      DataIterator crseDit = coarseBoxes.dataIterator();
      for (crseDit.begin(); crseDit.ok(); ++crseDit)
        {
          FluxBox& thisCoarReg = coarReg[crseDit];
          thisCoarReg.setVal(0.0);

          EdgeDataBox& thisEdgeData = m_regCoarse[crseDit];

          for (int edgeDir=0; edgeDir<SpaceDim; edgeDir++)
            {
              if (idir != edgeDir)
                {
                  FArrayBox& crseEdgeDataDir = thisEdgeData[edgeDir];
                  for (int faceDir = 0; faceDir<SpaceDim; faceDir++)
                    {
                      if (faceDir != edgeDir)
                        {
                          FArrayBox& faceData = thisCoarReg[faceDir];
                          int shiftDir = -1;
                          for (int i=0; i<SpaceDim; i++)
                            {
                              if ((i != faceDir) && (i != edgeDir) )
                                {
                                  shiftDir = i;
                                }
                            }
                          CH_assert(shiftDir >= 0);
                          crseEdgeDataDir.shiftHalf(shiftDir, sign(side()));
                          // scaling already taken care of in incrementCrse
                          Real scale = 1.0;
                          faceData.plus(crseEdgeDataDir, scale, 0, 0, faceData.nComp());
                          crseEdgeDataDir.shiftHalf(shiftDir, -sign(side()));
                        } // end if not normal direction
                    } // end loop over face directions
                } // end if edgeDir != idir
            } // end loop over edge directions
        } // end loop over crse boxes


      // first, we need to create a temp LevelData<FluxBox>
      // to make a local copy in the coarse layout space of
      // the fine register increments

      LevelData<FluxBox> fineRegLocal(coarReg.getBoxes(), m_nComp, IntVect::Unit);

      fineReg.copyTo(fineReg.interval(), fineRegLocal,
                     fineRegLocal.interval(),
                     m_crseCopiers[index(idir,side())]);

      for (DataIterator it = a_uCoarse.dataIterator(); it.ok(); ++it)
        {
          // loop over flux components here
          for (int fluxComp=0; fluxComp < SpaceDim; fluxComp++)
            {
              // we don't do anything in the normal direction
              if (fluxComp != idir)
                {
                  // fluxDir is the direction of the face-centered flux
                  FArrayBox& U = a_uCoarse[it()][fluxComp];
                  // set up IntVectSet to avoid double counting of updates
                  Box coarseGridBox = U.box();
                  // transfer to Cell-centered, then create IVS
                  coarseGridBox.shiftHalf(fluxComp,1);
                  IntVectSet nonUpdatedEdges(coarseGridBox);

                  // remember, we want to take the curl here
                  // also recall that fluxComp is the component
                  // of the face-centered curl (not the edge-centered
                  // vector field that we're refluxing, which is why
                  // the sign may seem like it's the opposite of what
                  // you might expect!
                  Real local_scale = -sign(side())*a_scale;
                  //int testDir = (fluxComp+1)%(SpaceDim);
                  if (((fluxComp+1)%(SpaceDim)) == idir)  local_scale *= -1;
                  Vector<IntVectSet>& ivsV =
                    m_refluxLocations[index(idir, side())][it()][fluxComp];
                  Vector<DataIndex>&  indexV =
                    m_coarToCoarMap[index(idir, side())][it()];
                  IVSIterator iv;

                  for (int i=0; i<ivsV.size(); ++i)
                    {
                      iv.define(ivsV[i]);
                      const FArrayBox& coar = coarReg[indexV[i]][fluxComp];
                      const FArrayBox& fine = fineRegLocal[indexV[i]][fluxComp];
                      for (iv.begin(); iv.ok(); ++iv)
                        {
                          IntVect thisIV = iv();
                          if (nonUpdatedEdges.contains(thisIV))
                          {
                            for (int comp=0; comp <m_nComp; ++comp)
                              {
                                //Real coarVal = coar(thisIV, comp);
                                //Real fineVal = fine(thisIV, comp);
                                U(thisIV, comp) -=
                                local_scale*(coar(thisIV, comp)
                                             +fine(thisIV, comp));
                              }
                            nonUpdatedEdges -= thisIV;
                          }
                        }

                    }
                } // end if not normal face
            } // end loop over fluxbox directions
        } // end loop over coarse boxes
    } // end loop over sides
  } // end loop over directions
}


// write the contents of all the registers
void
LevelFluxRegisterEdge::dump()
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      dumpLoCoar(idir);
      dumpLoFine(idir);
      dumpHiCoar(idir);
      dumpHiFine(idir);
    }
}
void
LevelFluxRegisterEdge::dumpLoCoar(int a_idir)
{
  // doesn't do anything right now
}
void
LevelFluxRegisterEdge::dumpHiCoar(int a_idir)
{
  // doesn't actually do anything at the moment
}

void
LevelFluxRegisterEdge::dumpLoFine(int a_idir)
{
  // doesn't actually do anything at the moment
}
void
LevelFluxRegisterEdge::dumpHiFine(int a_idir)
{
  // doesn't actually do anything at the moment
}

#include "NamespaceFooter.H"
