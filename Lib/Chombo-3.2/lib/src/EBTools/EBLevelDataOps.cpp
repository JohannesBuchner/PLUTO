#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "SPMD.H"

#include "EBLevelDataOps.H"
#include "EBLevelDataOpsF_F.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include "BoxIterator.H"
#include "CH_Timer.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "EBArith.H"
#include "EBIndexSpace.H"
#include "EBLoadBalance.H"
#include "CornerCopier.H"
#include "NamespaceHeader.H"

/*****/
void
EBLevelDataOps::pruneCoveredBoxes(Vector<Box>&              a_boxes,
                                  const ProblemDomain&      a_domain,
                                  const EBIndexSpace*       a_ebisPtr)
{
  Vector<int> procs;
  EBLoadBalance(procs,  a_boxes, a_domain);
  DisjointBoxLayout dbl(a_boxes, procs);
  EBISLayout ebisl;
  a_ebisPtr->fillEBISLayout(ebisl, dbl, a_domain, 0);

  /**
     Algorithm:
     loop through boxes and figure out which are covered
     gather them all together and popd them out of the a_boxes
     We COULD have done this by just gathering all the uncovered boxes
     and doing gather-broadcast on that list but that would not preserve
     the box order even in the case of there being no covered boxes.
   */
  Vector<Box> coveredBoxesLocal;
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = dbl.get(dit());
      if (ebisl[dit()].isCovered(grid))
        {
          coveredBoxesLocal.push_back(grid);
        }
    }

  Vector<Vector<Box> > allCoveredBoxes;
  int baseproc = 0;
  gather(allCoveredBoxes, coveredBoxesLocal,  baseproc);
  Vector<Box> newBoxes;
  if (procID() == baseproc)
    {
      for (int ibox = 0; ibox < a_boxes.size(); ibox++)
        {
          const Box& oldBox = a_boxes[ibox];
          bool isCovered = false;
          for (int iproc = 0; iproc < allCoveredBoxes.size();  iproc++)
            {
              for (int icov = 0; icov < allCoveredBoxes[iproc].size();  icov++)
                {
                  if (oldBox == allCoveredBoxes[iproc][icov])
                    {
                      isCovered = true;
                    }
                }
            }
          if (!isCovered)
            {
              newBoxes.push_back(oldBox);
            }
        }
    }
  broadcast(newBoxes, baseproc);
  a_boxes = newBoxes;
}
/*****/
Real EBLevelDataOps::parallelSum(const Real& a_value)
{
  // Find the sum of all a_value's
  Real sum = a_value;
#ifdef CH_MPI
  Real sendBuf = a_value;
  int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_CH_REAL,MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelSum");
    }
#endif
  return sum;
}

/*****/
int  EBLevelDataOps::parallelSum(const int& a_value)
{
  // Find the sum of all a_value's
  int sum = a_value;
#ifdef CH_MPI
  int sendBuf = a_value;
  int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_INT,MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelSum");
    }
#endif
  return sum;
}

/*****/
long long EBLevelDataOps::parallelSum(const long long& a_value)
{
  // Find the sum of all a_value's
  long sum = a_value;
#ifdef CH_MPI
  long long sendBuf = a_value;
  int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_LONG_LONG,MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelSum");
    }
#endif
  return sum;
}

/*****/
int EBLevelDataOps::parallelMin(const int& a_value)
{
  // Find the minimum of a_value's
  int val = a_value;
#ifdef CH_MPI
  int sendBuf = a_value;
  int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_INT,MPI_MIN, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMin");
    }
#endif
  return val;
}

/*****/
int  EBLevelDataOps::parallelMax(const int& a_value)
{
  // Find the maximum of a_value's
  int val = a_value;
#ifdef CH_MPI
  int sendBuf = a_value;
  int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_INT,MPI_MAX, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMax");
    }
#endif
  return val;
}
/*****/
Real EBLevelDataOps::parallelMin(const Real& a_value)
{
  // Find the minimum of a_value's
  Real val = a_value;
#ifdef CH_MPI
  Real sendBuf = a_value;
  int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_CH_REAL,MPI_MIN, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMin");
    }
#endif
  return val;
}

/*****/
Real EBLevelDataOps::parallelMax(const Real& a_value)
{
  // Find the maximum of a_value's
  Real val = a_value;
#ifdef CH_MPI
  Real sendBuf = a_value;
  int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_CH_REAL,MPI_MAX, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMax");
    }
#endif
  return val;
}

/*****/
void EBLevelDataOps::averageCellToFace(EBFaceFAB&             a_faceData,
                                       const EBCellFAB&       a_cellData,
                                       const EBGraph&         a_ebGraph,
                                       const Box&             a_dblBox,
                                       const int&             a_ghostFluxTan,
                                       const int&             a_idir,
                                       const ProblemDomain&   a_domain,
                                       const int&             a_cellComp,
                                       const int&             a_faceComp)
{
  //so covered data won't be something awful
  a_faceData.setVal(0.0);
  //this function fills all faces in and surrounding the dblBox,
  //  ghost faces tangential to the dblBox are filled

  BaseFab<Real>&       regFaceData = a_faceData.getSingleValuedFAB();
  const BaseFab<Real>& regCellData = a_cellData.getSingleValuedFAB();

  //do faces inside the domain
  {
    //make a cell centered box that has non-domain faces surrounding it ghost cells tangential to a_idir
    Box cellBox = a_dblBox;
    cellBox.grow(a_ghostFluxTan);
    cellBox &= a_domain;
    cellBox.grow(a_idir,-a_ghostFluxTan);

    //get the surrounding faces box
    Box faceBox = cellBox;
    faceBox.surroundingNodes(a_idir);

    //do regular cells in fortran for the faces in faceBox
    FORT_EBAVECELLTOFACE(CHF_FRA1(      regFaceData, a_faceComp),
                       CHF_CONST_FRA1(regCellData, a_cellComp),
                       CHF_CONST_INT(a_idir),
                       CHF_BOX(faceBox));

    //fix up irregular faces
    IntVectSet ivsIrreg = a_ebGraph.getIrregCells(cellBox);
    FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingNoBoundary;
    for (FaceIterator faceit(ivsIrreg, a_ebGraph, a_idir, stopCritGrid); faceit.ok(); ++faceit)
      {
        const FaceIndex& face = faceit();

        a_faceData(face, a_faceComp) = 0.5*(a_cellData(face.getVoF(Side::Hi), a_cellComp) +
                                            a_cellData(face.getVoF(Side::Lo), a_cellComp));
      }
  }

  //fix up domain faces
  for (SideIterator outSide; outSide.ok(); ++outSide)
    {
      Box edge = adjCellBox(a_dblBox,a_idir,outSide(),1);
      if (!a_domain.contains(edge))
        {
          Side::LoHiSide inSide = flip(outSide());
          edge.shift(a_idir,sign(inSide));
          IntVectSet ivsDomainEdge(edge);

          //put reasonable values in covered faces at domain boundaries
          for (BoxIterator boxit(edge); boxit.ok(); ++boxit)
            {
              IntVect ivCell = boxit();
              IntVect ivFace = boxit();
              if (outSide() == Side::Hi)
                ivFace.shift(BASISV(a_idir));
              regFaceData(ivFace, a_faceComp) = regCellData(ivCell, a_cellComp);
            }

          FaceStop::WhichFaces stopDomain =  FaceStop::AllBoundaryOnly;
          for (FaceIterator faceit(ivsDomainEdge, a_ebGraph, a_idir, stopDomain); faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();

              //get the vof just inside of the domain on this face
              const VolIndex&   vof = face.getVoF(inSide);
              const Real&     value = a_cellData(vof, a_cellComp);

              //see if there are faces further inside of vof
              const Vector<FaceIndex> faces = a_ebGraph.getFaces(vof,a_idir,inSide);
              if (faces.size() == 1)
                {//linear extrapolation
                  const VolIndex&   vofNeigh = faces[0].getVoF(inSide);
                  const Real&     valueNeigh = a_cellData(vofNeigh, a_cellComp);
                  a_faceData(face, a_faceComp) = 0.5*(3.0*value - valueNeigh);
                }
              else
                {//zero order extrapolation for != 1 inside faces
                  a_faceData(face, a_faceComp) = value;
                }
            }
        }
    }
}
/*****/
//first do cell centered average on grown box.   then interpolate to centroids
void EBLevelDataOps::averageCellToFace(EBFaceFAB           &      a_fluxData,
                                       const EBCellFAB     &      a_cellData,
                                       const Box           &      a_grid,
                                       const EBISBox       &      a_ebisBox,
                                       const ProblemDomain &      a_domain,
                                       int isrc, int idst, int inco,
                                       bool a_interpolateToCentroid)
{
  int idir = a_fluxData.direction();
  EBFaceFAB cellCenteredFlux;
  Box ccFluxBox = a_grid;
  if (a_interpolateToCentroid)
    {
      ccFluxBox.grow(1);
      ccFluxBox &= a_domain;
    }

  cellCenteredFlux.define(a_ebisBox, ccFluxBox, idir, a_fluxData.nComp());
  cellCenteredFlux.setVal(0.);

  faceCenteredAverageCellsToFaces(cellCenteredFlux, a_cellData, ccFluxBox, a_ebisBox, a_domain, isrc, idst, inco);
  //first copy (this does all the regular cells)
  Interval srcInt(isrc, isrc+inco-1);
  Interval dstInt(idst, idst+inco-1);
  a_fluxData.copy(a_grid, dstInt, a_grid,  cellCenteredFlux, srcInt);

  if (a_interpolateToCentroid)
    {
      //if required, do the fancy interpolation to centroids.
      IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(a_grid);
      IntVectSet cfivs;
      FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ivsIrreg, a_ebisBox.getEBGraph(), idir, stopCritGrid); faceit.ok(); ++faceit)
        {
          FaceStencil sten = EBArith::getInterpStencil(faceit(), cfivs, a_ebisBox, a_domain);
          for (int icomp = 0; icomp < inco; icomp++)
            {
              sten.setAllVariables(isrc+icomp);
              a_fluxData(faceit(), idst+icomp) = applyFaceStencil(sten, cellCenteredFlux, 0);
            }
        }
    }
}
/****/
void EBLevelDataOps::faceCenteredAverageCellsToFaces(EBFaceFAB           &      a_faceData,
                                                     const EBCellFAB     &      a_cellData,
                                                     const Box           &      ccFluxBox,
                                                     const EBISBox       &      a_ebisBox,
                                                     const ProblemDomain &      a_domain,
                                                     int isrc, int idst, int inco)
{
  //get the surrounding faces box
  Box faceBox = ccFluxBox;
  int idir = a_faceData.direction();
  faceBox.surroundingNodes(idir);
  //do regular cells in fortran for the faces in faceBox
  for (int icomp = 0; icomp < inco; icomp++)
    {
      int faceComp= idst + icomp;
      int cellComp= isrc + icomp;

      FORT_EBAVECELLTOFACE(CHF_FRA1(      a_faceData.getSingleValuedFAB(), faceComp),
                         CHF_CONST_FRA1(a_cellData.getSingleValuedFAB(), cellComp),
                         CHF_CONST_INT(idir),
                         CHF_BOX(faceBox));
    }

  //fix up irregular faces and boundary faces
  IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(ccFluxBox);
  Box interiorBox = a_domain.domainBox();
  interiorBox.grow(idir, 1);
  interiorBox &= a_domain.domainBox();
  interiorBox.grow(idir, -1);

  IntVectSet ivsEdge(a_domain.domainBox());
  ivsEdge -= interiorBox;
  ivsEdge &= ccFluxBox;
  ivsIrreg  |= ivsEdge;

  FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
  for (FaceIterator faceit(ivsIrreg, a_ebisBox.getEBGraph(), idir, stopCritGrid); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      if (!face.isBoundary())
        {
          for (int icomp = 0; icomp < inco; icomp++)
            {
              a_faceData(face, idst + icomp) = 0.5*(a_cellData(face.getVoF(Side::Hi), isrc + icomp) +
                                                    a_cellData(face.getVoF(Side::Lo), isrc + icomp));
            }
        }
      else
        {

          for (int icomp = 0; icomp < inco; icomp++)
            {
              VolIndex whichVoF;
              const Box& domainBox = a_domain.domainBox();
              if (domainBox.contains(face.gridIndex(Side::Lo)))
                {
                  whichVoF = face.getVoF(Side::Lo);
                }
              else if (domainBox.contains(face.gridIndex(Side::Hi)))
                {
                  whichVoF = face.getVoF(Side::Hi);
                }
              else
                {
                  MayDay::Error("face and domain inconsistent:  logic bust in average cells to faces");
                }

              a_faceData(face, idst + icomp) = a_cellData(whichVoF, isrc + icomp);
            }
        }
    }

}

/*****/
void EBLevelDataOps::averageCellToFace(LevelData<EBFluxFAB>&         a_fluxData,
                                       const LevelData<EBCellFAB>&   a_cellData,
                                       const DisjointBoxLayout&      a_grids,
                                       const EBISLayout&             a_ebisl,
                                       const ProblemDomain&          a_domain,
                                       int isrc, int idst, int inco,
                                       bool a_interpolateToCentroid)
{

  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          averageCellToFace(a_fluxData[dit()][idir],
                            a_cellData[dit()],
                            a_grids[dit()],
                            a_ebisl[dit()],
                            a_domain, isrc, idst, inco,
                            a_interpolateToCentroid);
        }
    }
}

/*****/
void EBLevelDataOps::averageCellToFaces(LevelData<EBFluxFAB>&         a_fluxData,
                                        const LevelData<EBCellFAB>&   a_cellData,
                                        const DisjointBoxLayout&      a_grids,
                                        const EBISLayout&             a_ebisl,
                                        const ProblemDomain&          a_domain,
                                        const int&                    a_comp)
{//fills ghost faces only if they are tangential to the face dir
  IntVect ghostCellIn = a_cellData.ghostVect();
  IntVect ghostFluxIn = a_fluxData.ghostVect();
  CH_assert(ghostCellIn[0]==ghostCellIn[1]);
  CH_assert(ghostCellIn[0]==ghostCellIn[SpaceDim-1]);
  CH_assert(ghostFluxIn[0]==ghostFluxIn[1]);
  CH_assert(ghostFluxIn[0]==ghostFluxIn[SpaceDim-1]);
  int ghostFluxTan   = ghostCellIn[0];
  ghostFluxTan = std::min(ghostFluxIn[0],ghostFluxTan);
  CH_assert(ghostFluxIn[0] >= ghostFluxTan);

  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& cellData = a_cellData[dit()];
      const EBGraph&    ebgraph = a_ebisl[dit()].getEBGraph();
      const Box&         dblBox = a_grids.get(dit());
      EBFluxFAB&       fluxData = a_fluxData[dit()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB&       faceData = fluxData[idir];
          averageCellToFace(faceData, cellData, ebgraph, dblBox, ghostFluxTan, idir, a_domain, a_comp, a_comp);
        }
    }
  exchangeComp(a_fluxData,a_comp);
}

void EBLevelDataOps::averageCellToFacesMAC(LevelData<EBFluxFAB>&         a_fluxData,
                                           const LevelData<EBCellFAB>&   a_cellData,
                                           const DisjointBoxLayout&      a_grids,
                                           const EBISLayout&             a_ebisl,
                                           const ProblemDomain&          a_domain)
{
  int cellNComp = a_cellData.nComp();
  int faceNComp = a_fluxData.nComp();
  CH_assert(cellNComp == SpaceDim);
  CH_assert(faceNComp == 1);
  IntVect ghostCellIn = a_cellData.ghostVect();
  IntVect ghostFluxIn = a_fluxData.ghostVect();
  IntVect ghostFlux   = ghostCellIn - IntVect::Unit;
  ghostFlux = min(ghostFluxIn,ghostFlux);
  ghostFlux = min(ghostFlux,IntVect::Unit);
  CH_assert(ghostFluxIn >= ghostFlux);

  int faceComp = 0;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& cellData = a_cellData[dit()];
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      const Box&      dblBox = a_grids.get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB&       faceData = a_fluxData[dit()][idir];
          averageCellToFace(faceData, cellData, ebgraph, dblBox, ghostFlux[0], idir, a_domain, idir, faceComp);
        }
    }
  exchangeComp(a_fluxData,faceComp);
}

void EBLevelDataOps::exchangeAll(LevelData<EBCellFAB>& a_data)

{
  const int ncomp = a_data.nComp();
  Interval interval(0,ncomp-1);
  a_data.exchange(interval);
  return;
}
void EBLevelDataOps::exchangeCorners(LevelData<EBCellFAB>& a_data,
                                     const ProblemDomain&  a_domain)

{
  const int ncomp = a_data.nComp();
  Interval interval(0,ncomp-1);
  CornerCopier corners(a_data.disjointBoxLayout(),a_data.disjointBoxLayout(),a_domain,IntVect::Unit);
  a_data.exchange(interval,corners);
  return;
}

void EBLevelDataOps::exchangeComp(LevelData<EBCellFAB>& a_data,
                                  const int&            a_comp)
{
  const int ncomp = a_data.nComp();
  CH_assert(a_comp<ncomp);
  Interval interval(a_comp,a_comp);
  a_data.exchange(interval);
  return;
}

void EBLevelDataOps::exchangeAll(LevelData<EBFluxFAB>& a_data)
{
  const int ncomp = a_data.nComp();
  Interval interval(0,ncomp-1);
  a_data.exchange(interval);
  return;
}

void EBLevelDataOps::exchangeComp(LevelData<EBFluxFAB>& a_data,
                                  const int&            a_comp)
{
  const int ncomp = a_data.nComp();
  CH_assert(a_comp<ncomp);
  Interval interval(a_comp,a_comp);
  a_data.exchange(interval);
  return;
}


void EBLevelDataOps::setCoveredVal(LevelData<EBCellFAB>&a_data,
                                   const Real&          a_value)
{
  CH_TIME("EBLevelDataOps::setCoveredVal(cell)");
  int ncomp = a_data.nComp();
  for (DataIterator dit = a_data.dataIterator();dit.ok();++dit)
    {
      for (int icomp = 0; icomp < ncomp;++icomp)
        {
          a_data[dit()].setCoveredCellVal(a_value,icomp);
        }
    }
}
void EBLevelDataOps::setCoveredVal(LevelData<EBFluxFAB>&a_data,
                                   const Real&          a_value)
{
  CH_TIME("EBLevelDataOps::setCoveredVal(face)");
  int ncomp = a_data.nComp();
  for (DataIterator dit = a_data.dataIterator();dit.ok();++dit)
    {
      for (int icomp = 0; icomp < ncomp;++icomp)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              EBFaceFAB& faceFAB = a_data[dit()][idir];
              faceFAB.setCoveredFaceVal(a_value,icomp);
            }
        }
    }
}

void EBLevelDataOps::setCoveredVal(LevelData<EBFluxFAB>&a_data,
                                   const int&           a_comp,
                                   const Real&          a_value)
{
  CH_TIME("EBLevelDataOps::setCoveredVal(face,comp)");
  for (DataIterator dit = a_data.dataIterator();dit.ok();++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_data[dit()][idir].setCoveredFaceVal(a_value,a_comp);
        }
    }
}
void EBLevelDataOps::setIrregVal(LevelData<EBCellFAB>&    a_data,
                                 const DisjointBoxLayout& a_dbl,
                                 const EBISLayout&        a_ebisl,
                                 const Real&              a_value)
{
  int ncomp = a_data.nComp();
  for (DataIterator dit = a_data.dataIterator();dit.ok();++dit)
    {
      for (int icomp = 0; icomp < ncomp;++icomp)
        {
          const Box& grid = a_dbl.get(dit());
          const EBISBox& ebisBox = a_ebisl[dit()];
          const IntVectSet& ivsIrreg = ebisBox.getIrregIVS(grid);
          for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex vof = vofit();
              a_data[dit()](vof,icomp) = a_value;
            }
        }
    }
}

void EBLevelDataOps::setCoveredVal(LevelData<EBCellFAB>&a_data,
                                   const int&           a_comp,
                                   const Real&          a_value)
{
  CH_TIME("EBLevelDataOps::setCoveredVal(cell,comp)");
  for (DataIterator dit = a_data.dataIterator();dit.ok();++dit)
    {
      a_data[dit()].setCoveredCellVal(a_value,a_comp);
    }
}
bool EBLevelDataOps::checkNANINF(const LevelData<EBCellFAB>&a_data,
                                 const IntVect&             a_iv1,
                                 const IntVect&             a_iv2,
                                 const Real&                a_shift)
{
  const bool checkSymmetry = (a_iv1 != a_iv2);

  //this function checks for nans and infs
  bool dataIsNANINF = false;
  int ncomp = a_data.nComp();
  for (int icomp = 0; icomp < ncomp;++icomp)
    {
      Real val1 = 1.e27;
      Real val2 = 1.e37;

      for (DataIterator dit=a_data.dataIterator();dit.ok();++dit)
        {
          const EBCellFAB& dataEBFAB = a_data[dit()];
          const Box& region = dataEBFAB.getRegion();
          IntVectSet ivsBox(region);
          const EBISBox& ebisBox = dataEBFAB.getEBISBox();
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Real val = dataEBFAB(vof,icomp);
              if (isnan(val) || isinf(val) || Abs(val)>1.e40)
                {
                  pout() << "      icomp = " << icomp << " vof = " << vof << " val = " << val << std::endl;
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
      MayDay::Warning("Found a NaN or Infinity.");
    }
  return dataIsNANINF;
}
void EBLevelDataOps::getMaxMin(Real&                       a_maxVal,
                               Real&                       a_minVal,
                               const LevelData<EBCellFAB>& a_data,
                               const int&                  a_comp,
                               const bool&                 a_doAbs)
{
  CH_TIME("EBLevelDataOps::getMaxMin");
  //this function gets the max and min (valid) values
  a_maxVal = -1.e99;
  a_minVal =  1.e99;

  const DisjointBoxLayout& dbl = a_data.disjointBoxLayout();
  for (DataIterator dit=a_data.dataIterator();dit.ok();++dit)
    {
      const EBCellFAB& dataEBFAB = a_data[dit()];
      const Box&        dblBox   = dbl.get(dit());
      const IntVectSet ivsBox(dblBox);
      const EBISBox& ebisBox = dataEBFAB.getEBISBox();
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const Real& val = dataEBFAB(vof,a_comp);
          if (a_doAbs)
            {
              a_maxVal = Max(a_maxVal,Abs(val));
              a_minVal = Min(a_minVal,Abs(val));
            }
          else
            {
              a_maxVal = Max(a_maxVal,val);
              a_minVal = Min(a_minVal,val);
            }
        }
    }
  a_minVal = EBLevelDataOps::parallelMin(a_minVal);
  a_maxVal = EBLevelDataOps::parallelMax(a_maxVal);
}

void EBLevelDataOps::averageMultiVofsToRegFAB(LevelData<EBCellFAB>&    a_data,
                                              const DisjointBoxLayout& a_dbl,
                                              const EBISLayout&        a_ebisl)
{
  const IntVect& ghostVect = a_data.ghostVect();
  const int&         nComp = a_data.nComp();
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&         data     = a_data[dit()];
      BaseFab<Real>&     regData  = data.getSingleValuedFAB();
      const EBISBox&     ebisBox  = a_ebisl[dit()];
      const Box&         dblBox   = a_dbl.get(dit());

      Box bigBox = dblBox;
      bigBox.grow(ghostVect);
      const IntVectSet&  multiIVS = ebisBox.getMultiCells(bigBox);
      for (IVSIterator ivsit(multiIVS); ivsit.ok(); ++ivsit)
        {
          const IntVect&            iv = ivsit();
          const Vector<VolIndex>& vofs = ebisBox.getVoFs(iv);

          for (int comp = 0; comp < nComp; comp++)
            {
              Real sumData = 0.0;
              Real sumFrac = 0.0;
              for (int i = 0; i < vofs.size(); i++)
                {
                  const VolIndex& thisVof  = vofs[i];
                  const Real&     kappa    = ebisBox.volFrac(thisVof);
                  sumFrac += kappa;
                  sumData += kappa*data(thisVof,comp);
                }
              sumData /= sumFrac;
              regData(iv, comp) = sumData;
            }
        }
    }
}

void EBLevelDataOps::copyToMultiVofsFromRegFAB(LevelData<EBCellFAB>&    a_data,
                                               const DisjointBoxLayout& a_dbl,
                                               const EBISLayout&        a_ebisl)
{
  const IntVect& ghostVect = a_data.ghostVect();
  const int&         nComp = a_data.nComp();
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&         data     = a_data[dit()];
      BaseFab<Real>&     regData  = data.getSingleValuedFAB();
      const EBISBox&     ebisBox  = a_ebisl[dit()];
      const Box&         dblBox   = a_dbl.get(dit());

      Box bigBox = dblBox;
      bigBox.grow(ghostVect);
      const IntVectSet&  multiIVS = ebisBox.getMultiCells(bigBox);
      for (IVSIterator ivsit(multiIVS); ivsit.ok(); ++ivsit)
        {
          const IntVect&            iv = ivsit();
          const Vector<VolIndex>& vofs = ebisBox.getVoFs(iv);

          for (int comp = 0; comp < nComp; comp++)
            {
              const Real& rData = regData(iv,comp);
              for (int i = 0; i < vofs.size(); i++)
                {
                  const VolIndex& thisVof  = vofs[i];
                  data(thisVof,comp) = rData;
                }
            }
        }
    }
}


void EBLevelDataOps::defineLevelData(LevelData<EBCellFAB>&    a_levelData,
                                     const EBISLayout&        a_ebisl,
                                     const DisjointBoxLayout& a_dbl,
                                     const IntVect&           a_ghosts,
                                     const int&               a_nComp)
{
  EBCellFactory ebcellfact(a_ebisl);
  a_levelData.define(a_dbl,a_nComp,a_ghosts,ebcellfact);
}

void EBLevelDataOps::defineLevelData(LevelData<EBFluxFAB>&    a_levelData,
                                     const EBISLayout&        a_ebisl,
                                     const DisjointBoxLayout& a_dbl,
                                     const IntVect&           a_ghosts,
                                     const int&               a_nComp)
{
  EBFluxFactory ebfluxfact(a_ebisl);
  a_levelData.define(a_dbl,a_nComp,a_ghosts,ebfluxfact);
}

void EBLevelDataOps::setToZero(LevelData<EBCellFAB>& a_result)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      a_result[dit()].setVal(0.0);
    }
}

void EBLevelDataOps::setToZero(LevelData<EBFluxFAB>& a_result)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      EBFluxFAB& result = a_result[dit()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          result[idir].setVal(0.0);
        }
    }
}


void EBLevelDataOps::setVal(LevelData<EBCellFAB>& a_result,
                            const Real&           a_value)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& result = a_result[d];

      result.setVal(a_value);
    }
}

void EBLevelDataOps::setVal(LevelData<BaseIVFAB<Real> >& a_result,
                            const Real&           a_value)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      BaseIVFAB<Real>& result = a_result[d];

      result.setVal(a_value);
    }
}
void EBLevelDataOps::setVal(LevelData<EBCellFAB>& a_result,
                            const Real&           a_value,
                            const int&            a_comp)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& result = a_result[d];

      result.setVal(a_comp,a_value);
    }
}

void EBLevelDataOps::setVal(LevelData<EBFluxFAB>& a_result,
                            const Real&           a_value)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBFluxFAB& result = a_result[d];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          result[idir].setVal(a_value);
        }
    }
}

void EBLevelDataOps::axby( LevelData<EBCellFAB>&       a_lhs,
                           const LevelData<EBCellFAB>& a_x,
                           const LevelData<EBCellFAB>& a_y,
                           const Real& a,
                           const Real& b)
{
  //  CH_assert(a_lhs.disjointBoxLayout() == a_x.disjointBoxLayout());

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& data = a_lhs[d];
      data.axby(a_x[d], a_y[d], a, b);
      //data.copy(a_x[d]);
      //data.mult(a);
      //data.plus(a_y[d], b);
    }
}
void EBLevelDataOps::axby( LevelData<EBCellFAB>&       a_lhs,
                           const LevelData<EBCellFAB>& a_x,
                           const LevelData<EBCellFAB>& a_y,
                           const Real& a,
                           const Real& b,
                           const int&  a_lhsComp,
                           const int&  a_xComp,
                           const int&  a_yComp)

{
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& data = a_lhs[d];
      data.axby(a_x[d], a_y[d], a, b,a_lhsComp,a_xComp,a_yComp);
    }
}

void EBLevelDataOps::assign(LevelData<EBCellFAB>&       a_to,
                            const LevelData<EBCellFAB>& a_from,
                            const Interval&             a_toInterval,
                            const Interval&             a_fromInterval)
{
  CH_TIME("EBLevelDataOps::assign(to,from,toInterval,fromInterval)");
  a_from.copyTo(a_fromInterval, a_to, a_toInterval);
}

void EBLevelDataOps::assign(LevelData<EBCellFAB>& a_lhs,
                            const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBLevelDataOps::assign(to,from)");
  Interval interv(0, a_rhs.nComp()-1);
  a_rhs.copyTo(interv, a_lhs, interv);
}
void EBLevelDataOps::clone(LevelData<EBCellFAB>& a_lhs,
                           const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBLevelDataOps::clone");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB&       lhs = a_lhs[d];
      const EBCellFAB& rhs = a_rhs[d];
      lhs.copy(rhs);
    }
}
void EBLevelDataOps::assign(LevelData<EBFluxFAB>& a_lhs,
                            const LevelData<EBFluxFAB>& a_rhs)
{
  CH_TIME("EBLevelDataOps::assign(to,from)");
  Interval interv(0, a_rhs.nComp()-1);
  a_rhs.copyTo(interv, a_lhs, interv);
}
void EBLevelDataOps::incr( LevelData<EBCellFAB>& a_lhs,
                           const LevelData<EBCellFAB>& a_rhs,
                           const Real& a_scale)
{
  //  CH_assert(a_lhs.disjointBoxLayout() == a_rhs.disjointBoxLayout());
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      a_lhs[d].plus(a_rhs[d], a_scale);
    }
}

void EBLevelDataOps::incr( LevelData<EBCellFAB>& a_lhs,
                           const Real& a_scale)
{
  //  CH_assert(a_lhs.disjointBoxLayout() == a_rhs.disjointBoxLayout());
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      a_lhs[d] += a_scale;
    }
}

void EBLevelDataOps::scale(LevelData<EBFluxFAB>& a_result,
                           const Real&           a_value)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBFluxFAB& result = a_result[d];
      for (int idir = 0;idir<SpaceDim;idir++)
        {
          result[idir] *= a_value;
        }
    }
}

void EBLevelDataOps::scale(LevelData<EBCellFAB>& a_result,
                           const LevelData<EBCellFAB>&           a_value)
{
  CH_assert(a_value.nComp() == 1);
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& result = a_result[d];
      for (int icomp = 0; icomp < a_result.nComp(); icomp++)
        {
          int isrc = 0; int idst = icomp; int nco = 1;
          result.mult(a_value[dit()], isrc, idst, nco);
        }
    }
}

void EBLevelDataOps::scale(LevelData<EBCellFAB>& a_result,
                           const Real&           a_value)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& result = a_result[d];

      result.mult(a_value);
    }
}

void EBLevelDataOps::scale(LevelData<EBCellFAB>& a_result,
                           const Real&           a_value,
                           const int&            a_comp)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& result = a_result[d];

      result.axby(result,result,0.0,a_value,a_comp,a_comp,a_comp);
    }
}

void EBLevelDataOps::sum(LevelData<EBCellFAB>&       a_result,
                         const LevelData<EBCellFAB>& a_in1,
                         const LevelData<EBCellFAB>& a_in2)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& result = a_result[d];

      const EBCellFAB& in1 = a_in1[d];
      const EBCellFAB& in2 = a_in2[d];

      result.copy(in1);
      result += in2;
    }
}

void EBLevelDataOps::addConstant(LevelData<EBCellFAB>& a_data,
                                 const Real&           a_constant)
{
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {

      DataIndex d = dit();
      EBCellFAB& data = a_data[d];

      data += a_constant;
    }
}
void EBLevelDataOps::power(LevelData<EBCellFAB>& a_data,
                           const Real&           a_exponent,
                           const int&            a_comp)
{
  CH_TIME("EBLevelDataOps::power");
  CH_assert(a_comp<a_data.nComp());

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& data = a_data[dit()];
      const Box& region = data.getRegion();
      const EBISBox& dataEBISBox = data.getEBISBox();
      const EBGraph& dataEBGraph = dataEBISBox.getEBGraph();

      const IntVectSet dataIVS = IntVectSet(region);
      for (VoFIterator vofit(dataIVS,dataEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const Real& dataIn = data(vof,a_comp);
          data(vof,a_comp) = pow(dataIn,a_exponent);
        }
    }
}

void EBLevelDataOps::product(LevelData<EBCellFAB>&       a_result,
                             const LevelData<EBCellFAB>& a_in1,
                             const LevelData<EBCellFAB>& a_in2)
{

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      EBCellFAB& result = a_result[d];

      const EBCellFAB& in1 = a_in1[d];
      const EBCellFAB& in2 = a_in2[d];

      result.copy(in1);
      result *= in2;
    }
}
void EBLevelDataOps::product(LevelData<EBCellFAB>&       a_result,
                             const LevelData<EBCellFAB>& a_in1,
                             const LevelData<EBCellFAB>& a_in2,
                             const int&                  a_rComp,
                             const int&                  a_1Comp,
                             const int&                  a_2Comp)
{
  const DisjointBoxLayout& dbl = a_result.disjointBoxLayout();
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      const Box& dblBox = dbl.get(dit());

      EBCellFAB& result = a_result[d];

      const EBCellFAB& in1 = a_in1[d];
      const EBCellFAB& in2 = a_in2[d];

      const Box& regionr = dblBox;
      const Box& region1 = dblBox;

      const Interval intR(a_rComp,a_rComp);
      const Interval int1(a_1Comp,a_1Comp);
//       const Interval int2(a_2Comp,a_2Comp);

      result.copy(region1,intR,regionr,in1,int1);
      result.mult(in2,a_2Comp,a_rComp,1);
    }
}

void EBLevelDataOps::divideVectorByScalar(LevelData<EBCellFAB>&       a_vectorOut,
                                          const LevelData<EBCellFAB>& a_vectorIn,
                                          const LevelData<EBCellFAB>& a_scalar)
{
  //a_vector /= a_scalar
  int nCompVectorOut = a_vectorOut.nComp();
  int nCompVectorIn  = a_vectorIn.nComp();
  int nCompScalar = a_scalar.nComp();
  CH_assert(nCompScalar==1);
  CH_assert(nCompVectorOut==nCompVectorIn);

  for (int icomp = 0;icomp<nCompVectorOut;++icomp)
    {
      divide(a_vectorOut,a_vectorIn,a_scalar,icomp,icomp,0);
    }
}

void EBLevelDataOps::divide(LevelData<EBCellFAB>&       a_result,
                            const LevelData<EBCellFAB>& a_in1,
                            const LevelData<EBCellFAB>& a_in2)
{

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      EBCellFAB& result = a_result[d];

      const EBCellFAB& in1 = a_in1[d];
      const EBCellFAB& in2 = a_in2[d];

      result.copy(in1);
      result /= in2;
    }
}
void EBLevelDataOps::divide(LevelData<EBCellFAB>&       a_result,
                            const LevelData<EBCellFAB>& a_in1,
                            const LevelData<EBCellFAB>& a_in2,
                            const int&                  a_rComp,
                            const int&                  a_1Comp,
                            const int&                  a_2Comp)
{
  const DisjointBoxLayout& dbl = a_result.disjointBoxLayout();
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      const Box& dblBox = dbl.get(dit());

      EBCellFAB& result = a_result[d];
      const EBCellFAB& in1 = a_in1[d];
      const EBCellFAB& in2 = a_in2[d];

      const Box& regionr = dblBox;
      const Box& region1 = dblBox;

      const Interval intR(a_rComp,a_rComp);
      const Interval int1(a_1Comp,a_1Comp);
//       const Interval int2(a_2Comp,a_2Comp);

      result.copy(region1,intR,regionr,in1,int1);
      result.divide(in2,a_2Comp,a_rComp,1);
    }
}
void EBLevelDataOps::invert(LevelData<EBCellFAB>&       a_result,
                            const LevelData<EBCellFAB>& a_in1)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      EBCellFAB& result = a_result[d];

      const EBCellFAB& in1 = a_in1[d];
      EBCellFAB one(in1.getEBISBox(),
                    in1.box(),
                    in1.nComp());
      one.setVal(1.0);

      one /= in1;
      result.copy(one);
    }
}
void EBLevelDataOps::product(LevelData<EBFluxFAB>&       a_result,
                             const LevelData<EBFluxFAB>& a_in1,
                             const LevelData<EBFluxFAB>& a_in2)
{

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB& result = a_result[d][idir];
          const EBFaceFAB& in1 = a_in1[d][idir];
          const EBFaceFAB& in2 = a_in2[d][idir];

          result.copy(in1);
          result *= in2;
        }
    }
}

//-----------------------------------------------------------------------
void EBLevelDataOps::kappaWeight(LevelData<EBCellFAB>& a_data)
{
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& data = a_data[d];

      kappaWeight(data);
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBLevelDataOps::kappaWeight(EBCellFAB& a_data)
{
  int nComp = a_data.nComp();

  const Box& dataBox = a_data.box();
  const EBISBox& dataEBISBox = a_data.getEBISBox();
  const EBGraph& dataEBGraph = dataEBISBox.getEBGraph();

  const IntVectSet& dataIrreg = dataEBISBox.getIrregIVS(dataBox);

  for (VoFIterator vofit(dataIrreg,dataEBGraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < nComp; comp++)
        {
          a_data(vof,comp) *= dataEBISBox.volFrac(vof);
        }
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBLevelDataOps::areaFracScalingWeight(LevelData<EBCellFAB>& a_data)
{
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      EBCellFAB& data = a_data[d];

      areaFracScalingWeight(data);
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBLevelDataOps::areaFracScalingWeight(EBCellFAB& a_data)
{
  int nComp = a_data.nComp();

  const Box& dataBox = a_data.box();
  const EBISBox& dataEBISBox = a_data.getEBISBox();
  const EBGraph& dataEBGraph = dataEBISBox.getEBGraph();

  const IntVectSet& dataIrreg = dataEBISBox.getIrregIVS(dataBox);

  for (VoFIterator vofit(dataIrreg,dataEBGraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < nComp; comp++)
        {
          a_data(vof,comp) *= dataEBISBox.areaFracScaling(vof);
        }
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBLevelDataOps::kappaScale(LevelData<EBCellFAB>& a_data,
                                const Real&           a_scale)
{
  kappaWeight(a_data);
  scale(a_data,a_scale);
}
//-----------------------------------------------------------------------

void EBLevelDataOps::gatherBroadCast(Real& a_accum, Real& a_volume, const int& a_p)
{
  //   Vector<Real> accum(1,a_accum);
//   gatherBroadCast(accum, a_volume, a_p);
//   a_accum = accum[0];
#ifdef CH_MPI
  Real tmp=a_volume;
  MPI_Allreduce(&tmp, &a_volume, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  tmp = a_accum;
  if (a_p==0)
    {
      MPI_Allreduce(&tmp, &a_accum, 1, MPI_CH_REAL,
                    MPI_MAX, Chombo_MPI::comm);
    }
  else
    {
      MPI_Allreduce(&tmp, &a_accum, 1, MPI_CH_REAL,
                    MPI_SUM, Chombo_MPI::comm);
    }
#endif

}

void EBLevelDataOps::gatherBroadCast(Vector<Real>& a_accum, Real& a_volume, const int& a_p)
{
#ifdef CH_MPI
  {
    CH_TIME("MPI_Barrier EBLevelDataOps::gatherBroadcast");
      MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  //gather what each processor thinks is the accum and volume
  int ncomp = a_accum.size();
  Real volumeLoc = a_volume;
  Vector<Real> accumLoc  = a_accum;

  Vector<Real> volumeVec;
  Vector<Vector<Real> > accumVec;
  int baseproc = 0;
  gather(volumeVec, volumeLoc, baseproc);
  gather( accumVec,  accumLoc, baseproc);

  a_volume = 0.;
  for (int i=0; i<ncomp; i++) a_accum[i]=0.0;
  if (procID() == baseproc)
    {
      for (int ivec = 0; ivec < numProc(); ivec++)
        {
          a_volume += volumeVec[ivec];
          Vector<Real> cur =   accumVec[ivec];
          for (int i=0; i<ncomp; i++)
          {
            if (a_p == 0)
              {
                if (cur[i] > a_accum[i])
                  {
                    a_accum[i] = cur[i];
                  }
              }
            else
              {
                a_accum[i]  += cur[i];
              }
          }
        }
    }

  //broadcast the sum to all processors.
  broadcast( a_accum, baseproc);
  broadcast(a_volume, baseproc);
}
Real EBLevelDataOps::kappaNorm(Real&                       a_volume,
                               const LevelData<EBCellFAB>& a_data,
                               int                         a_which,
                               const ProblemDomain&        a_domain,
                               int                         a_p)
{
  Real accum = 0.0;
  a_volume = 0.0;
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      const EBCellFAB& data = a_data[d];
      const Box& box = a_data.getBoxes().get(d);
      Real cur;
      Real curVolume;

      cur = sumKappaPow(curVolume,data,box,a_which,a_domain,a_p);

      a_volume += curVolume;

      if (a_p == 0)
        {
          if (cur > accum)
            {
              accum = cur;
            }
        }
      else
        {
          accum += cur;
        }
    }

  //do broadcast-gather thing on volume and accumulation
  gatherBroadCast(accum, a_volume, a_p);

  if (a_p != 0)
    {
      if (a_volume > 0.0)
        {
          accum = accum / a_volume;
        }

      accum = pow(accum,Real(1.0)/a_p);
    }


  return accum;
}

Vector<Real> EBLevelDataOps::vectorKappaNorm(Real&                       a_volume,
                                     const LevelData<EBCellFAB>& a_data,
                                     int                         a_which,
                                     const ProblemDomain&        a_domain,
                                     int                         a_p)
{
  CH_TIME("EBLevelDataOps::vectorKappaNorm");

  int ncomp = a_data.nComp();
  Vector<Real> accum(ncomp,0.0);
  Vector<Real> cur;
  a_volume = 0.0;
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      const EBCellFAB& data = a_data[d];
      const Box& box = a_data.getBoxes().get(d);
      Real curVolume;

      cur = vectorSumKappaPow(curVolume,data,box,a_which,a_domain,a_p);

      a_volume += curVolume;

      for (int i=0; i<ncomp; i++)
      {
        if (a_p == 0)
          {
            if (cur[i] > accum[i])
              {
                accum[i] = cur[i];
              }
          }
        else
          {
            accum[i] += cur[i];
          }
      }
    }

  //AttachDebugger(3);
  //do broadcast-gather thing on volume and accumulation
  gatherBroadCast(accum, a_volume, a_p);

  for (int i=0; i<ncomp; i++)
  {
    if (a_p != 0)
      {
        if (a_volume > 0.0)
          {
            accum[i] = accum[i] / a_volume;
          }

        accum[i] = pow(accum[i], Real(1.0)/a_p);
      }

  }
  return accum;
}

Real EBLevelDataOps::noKappaNorm(Real&                       a_volume,
                                 const LevelData<EBCellFAB>& a_data,
                                 int                         a_which,
                                 const ProblemDomain&        a_domain,
                                 int                         a_p)
{
  Real accum = 0.0;

  a_volume = 0.0;

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      const EBCellFAB& data = a_data[d];
      const Box& box = a_data.getBoxes().get(d);
      Real cur;
      Real curVolume;

      cur = sumNoKappaPow(curVolume,data,box,a_which,a_domain,a_p);

      a_volume += curVolume;

      if (a_p == 0)
        {
          if (cur > accum)
            {
              accum = cur;
            }
        }
      else
        {
          accum += cur;
        }
    }

  gatherBroadCast(accum, a_volume, a_p);
  if (a_p != 0)
    {
      if (a_volume > 0.0)
        {
          accum = accum / a_volume;
        }

      accum = pow(accum,Real(1.0)/a_p);
    }

  return accum;
}

Real EBLevelDataOps::kappaDotProduct(Real&                       a_volume,
                                     const LevelData<EBCellFAB>& a_data1,
                                     const LevelData<EBCellFAB>& a_data2,
                                     int                         a_which,
                                     const ProblemDomain&        a_domain)
{
  CH_TIME("EBLevelDataOps::kappaDotProduct");
  Real accum = 0.0;

  a_volume = 0.0;

  for (DataIterator dit = a_data1.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      const EBCellFAB& data1 = a_data1[d];
      const EBCellFAB& data2 = a_data2[d];
      const Box& box = a_data1.getBoxes().get(d);
      Real cur;
      Real curVolume;

      cur = sumKappaDotProductAllCells(curVolume,data1,data2,box,a_which,a_domain); // This is the new, fast version
      //cur = sumKappaDotProduct(curVolume,data1,data2,box,a_which,a_domain);
      a_volume += curVolume;
      accum += cur;
    }

  gatherBroadCast(accum, a_volume, 1);

  if (a_volume > 0.0)
    {
      accum = accum / a_volume;
    }

  return accum;
}

Real EBLevelDataOps::noKappaDotProduct(Real&                       a_volume,
                                       const LevelData<EBCellFAB>& a_data1,
                                       const LevelData<EBCellFAB>& a_data2,
                                       int                         a_which,
                                       const ProblemDomain&        a_domain)
{
  Real accum = 0.0;

  a_volume = 0.0;

  for (DataIterator dit = a_data1.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      const EBCellFAB& data1 = a_data1[d];
      const EBCellFAB& data2 = a_data2[d];
      const Box& box = a_data1.getBoxes().get(d);
      Real cur;
      Real curVolume;

      cur = sumNoKappaDotProduct(curVolume,data1,data2,box,a_which,a_domain);

      a_volume += curVolume;
      accum += cur;
    }

  gatherBroadCast(accum, a_volume, 1);

  if (a_volume > 0.0)
    {
      accum = accum / a_volume;
    }

  return accum;
}

Real EBLevelDataOps::kappaSumLevel(Real&                       a_volume,
                                   const LevelData<EBCellFAB>& a_data,
                                   int                         a_which,
                                   const ProblemDomain&        a_domain)
{
  Real accum = 0.0;

  a_volume = 0.0;

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      const EBCellFAB& data = a_data[d];
      const Box& box = a_data.getBoxes().get(d);
      Real cur;
      Real curVolume;

      cur = sumKappa(curVolume,data,box,a_which,a_domain);

      a_volume += curVolume;
      accum += cur;
    }

  gatherBroadCast(accum, a_volume, 1);

  if (a_volume > 0.0)
    {
      accum = accum / a_volume;
    }

  return accum;
}

Real EBLevelDataOps::noKappaSumLevel(Real&                       a_volume,
                                     const LevelData<EBCellFAB>& a_data,
                                     int                         a_which,
                                     const ProblemDomain&        a_domain)
{
  Real accum = 0.0;

  a_volume = 0.0;
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      const EBCellFAB& data = a_data[d];
      const Box& box = a_data.getBoxes().get(d);
      Real cur;
      Real curVolume;

      cur = sumNoKappa(curVolume,data,box,a_which,a_domain);

      a_volume += curVolume;
      accum += cur;
    }

  gatherBroadCast(accum, a_volume, 1);

  if (a_volume > 0.0)
    {
      accum = accum / a_volume;
    }

  return accum;
}

// Protected functions - for internal use
Real EBLevelDataOps::sumKappaPow(Real&                a_volume,
                                 const EBCellFAB&     a_data, const Box& curBox,
                                 int                  a_which,
                                 const ProblemDomain& a_domain,
                                 int                  a_p)
{
  Vector<Real> rtn;
  rtn = vectorSumKappaPow(a_volume, a_data, curBox, a_which, a_domain, a_p);
  return rtn[0];
}

Vector<Real> EBLevelDataOps::vectorSumKappaPow(Real&                a_volume,
                                               const EBCellFAB&     a_data, const Box& curBox,
                                               int                  a_which,
                                               const ProblemDomain& a_domain,
                                               int                  a_p)
{
  int ncomp = a_data.nComp();

  Vector<Real> val(ncomp,0);
  Vector<Real> sum(ncomp,0);
  //const Box& curBox = a_data.box();
  const EBISBox& curEBISBox = a_data.getEBISBox();


  a_volume = 0.0;


  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) == EBLEVELDATAOPS_INTERIORREGVOFS &&
     (a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) == EBLEVELDATAOPS_BOUNDARYREGVOFS &&
     (a_which & EBLEVELDATAOPS_INTERIORIRREGVOFS) == EBLEVELDATAOPS_INTERIORIRREGVOFS &&
     (a_which & EBLEVELDATAOPS_BOUNDARYIRREGVOFS) == EBLEVELDATAOPS_BOUNDARYIRREGVOFS)
    {
      // ok, do everything norm optimization
      for (BoxIterator bit(curBox); bit.ok();  ++bit)
        {
          const IntVect& iv = bit();
          if (curEBISBox.isCovered(iv))
          {
            // do nothing
          }
          else if (curEBISBox.numVoFs(iv) == 1)
            {
              VolIndex vof(iv,0);
              Real volFrac = curEBISBox.volFrac(vof);
              a_volume += volFrac;

              for (int comp=0; comp<ncomp; comp++)
              {
                val[comp]     = volFrac * Abs(a_data(vof,comp));
              }
              if (a_p == 0)
                {
                  for (int comp=0; comp<ncomp; comp++)
                  {
                    sum[comp] = Max(sum[comp],val[comp]);
                  }
                }
              else
                {
                  for (int comp=0; comp<ncomp; comp++)
                  {
                    Real p = a_p;
                    Real integrand = pow(val[comp],p);
                    sum[comp] += integrand;
                  }
                }

            }
          else
            {
              //multi-valued
              Vector<VolIndex> vofs = curEBISBox.getVoFs(iv);
              for (int v=0; v<vofs.size();v++)
              {
                const VolIndex& vof = vofs[v];
                Real volFrac = curEBISBox.volFrac(vof);
                a_volume += volFrac;
                for (int comp=0; comp<ncomp; comp++)
                {
                  val[comp]     = volFrac * Abs(a_data(vof,comp));
                }
                if (a_p == 0)
                  {
                    for (int comp=0; comp<ncomp; comp++)
                    {
                      sum[comp] = Max(sum[comp],val[comp]);
                    }
                  }
                else
                  {
                    for (int comp=0; comp<ncomp; comp++)
                    {
                      Real p = a_p;
                      Real integrand = pow(val[comp],p);
                      sum[comp] += integrand;
                    }
                  }
              }
            }
        }
      return sum;
    }



  IntVectSet regIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) == EBLEVELDATAOPS_INTERIORREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      regIvsRegion = IntVectSet(interior);
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) == EBLEVELDATAOPS_BOUNDARYREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      regIvsRegion |= boundary;
    }

  for (VoFIterator vofit(regIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isRegular(iv))
        {
          for (int comp=0; comp<ncomp; comp++)
          {
            val[comp] = Abs(a_data(vof,comp));
          }
          a_volume += 1.0;

          if (a_p == 0)
            {
              for (int comp=0; comp<ncomp; comp++)
              {
                sum[comp] = Max(sum[comp],val[comp]);
              }
            }
          else
            {
              for (int comp=0; comp<ncomp; comp++)
              {
                Real p = a_p;
                Real integrand = pow(val[comp],p);

                sum[comp] += integrand;
              }
            }
        }
    }

  IntVectSet irregIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORIRREGVOFS) == EBLEVELDATAOPS_INTERIORIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      irregIvsRegion = IntVectSet(interior);
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYIRREGVOFS) == EBLEVELDATAOPS_BOUNDARYIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      irregIvsRegion |= boundary;
    }

  for (VoFIterator vofit(irregIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isIrregular(iv))
        {
          Real volfrac = curEBISBox.volFrac(vof);
          for (int comp=0; comp<ncomp; comp++)
          {
            val[comp]     = volfrac * Abs(a_data(vof,comp));
          }

          a_volume += volfrac;

          if (a_p == 0)
            {
              for (int comp=0; comp<ncomp; comp++)
              {
                sum[comp] = Max(sum[comp],val[comp]);
              }
            }
          else
            {
              for (int comp=0; comp<ncomp; comp++)
              {
                Real p = a_p;
                Real integrand = pow(val[comp],p);

                sum[comp] += integrand;
              }
            }
        }
    }

  return sum;
}

Real EBLevelDataOps::sumNoKappaPow(Real&                a_volume,
                                   const EBCellFAB&     a_data, const Box& curBox,
                                   int                  a_which,
                                   const ProblemDomain& a_domain,
                                   int                  a_p)
{
  CH_assert(a_data.nComp() == 1);

  //const Box& curBox = a_data.box();
  const EBISBox& curEBISBox = a_data.getEBISBox();
  int comp = 0;

  Real sum = 0.0;

  a_volume = 0.0;

  IntVectSet regIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) == EBLEVELDATAOPS_INTERIORREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      regIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) == EBLEVELDATAOPS_BOUNDARYREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      regIvsRegion |= boundary;
    }

  for (VoFIterator vofit(regIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isRegular(iv))
        {
          Real val = Abs(a_data(vof,comp));

          a_volume += 1.0;

          if (a_p == 0)
            {
              sum = Max(sum,val);
            }
          else
            {
              Real p = a_p;
              Real integrand = pow(val,p);

              sum += integrand;
            }
        }
    }

  IntVectSet irregIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORIRREGVOFS) == EBLEVELDATAOPS_INTERIORIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      irregIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYIRREGVOFS) == EBLEVELDATAOPS_BOUNDARYIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      irregIvsRegion |= boundary;
    }

  for (VoFIterator vofit(irregIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isIrregular(iv))
        {
          Real volfrac = curEBISBox.volFrac(vof);
          Real val     = Abs(a_data(vof,comp));

          a_volume += volfrac;

          if (a_p == 0)
            {
              sum = Max(sum,val);
            }
          else
            {
              Real p = a_p;
              Real integrand = pow(val,p);

              sum += integrand;
            }
        }
    }

  return sum;
}

Real EBLevelDataOps::sumKappaDotProduct(Real&                a_volume,
                                        const EBCellFAB&     a_data1,
                                        const EBCellFAB&     a_data2,
                                        const Box& curBox,
                                        int                  a_which,
                                        const ProblemDomain& a_domain)
{


  // const Box& curBox1 = a_data1.box();
  const EBISBox& curEBISBox1 = a_data1.getEBISBox();

  int ncomp = a_data1.nComp();

  Real sum = 0.0;

  a_volume = 0.0;
  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) != EBLEVELDATAOPS_INTERIORREGVOFS)
    {
      MayDay::Error("This code has been optimized to not handle the boundary-only case");
    }


  Box region = curBox;
  if ((a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) != EBLEVELDATAOPS_BOUNDARYREGVOFS)
    {
      region.grow(1);
      region &= a_domain;
      region.grow(-1);
    }

  for (BoxIterator bit(region); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if (curEBISBox1.isCovered(iv))
      {
        // do nothing
      }
      else if (curEBISBox1.numVoFs(iv) == 1)
        {
          VolIndex vof(iv,0);
          Real volFrac = curEBISBox1.volFrac(vof);
          a_volume += volFrac;
          for (int comp=0; comp<ncomp; comp++)
          {
            Real val1 = volFrac * a_data1(vof,comp);
            Real val2 = volFrac * a_data2(vof,comp);
            sum += val1*val2;
          }
        }
      else
        {
          //multi-valued
          Vector<VolIndex> vofs = curEBISBox1.getVoFs(iv);
          for (int v=0; v<vofs.size();v++)
          {
            const VolIndex& vof = vofs[v];
            Real volFrac = curEBISBox1.volFrac(vof);
            a_volume += volFrac;
            for (int comp=0; comp<ncomp; comp++)
            {
              Real val1 = volFrac * a_data1(vof,comp);
              Real val2 = volFrac * a_data2(vof,comp);
              sum += val1*val2;
            }
          }
        }
    }



  return sum;
}

Real EBLevelDataOps::sumNoKappaDotProduct(Real&                a_volume,
                                          const EBCellFAB&     a_data1,
                                          const EBCellFAB&     a_data2,const Box& curBox1,
                                          int                  a_which,
                                          const ProblemDomain& a_domain)
{
  CH_assert(a_data1.nComp() == 1);
  CH_assert(a_data2.nComp() == 1);

  const EBISBox& curEBISBox1 = a_data1.getEBISBox();
  int comp1 = 0;

  int comp2 = 0;



  Real sum = 0.0;

  a_volume = 0.0;

  IntVectSet regIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) == EBLEVELDATAOPS_INTERIORREGVOFS)
    {
      Box interior = curBox1;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      regIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) == EBLEVELDATAOPS_BOUNDARYREGVOFS)
    {
      Box interior = curBox1;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox1);
      boundary -= interior;

      regIvsRegion |= boundary;
    }

  for (VoFIterator vofit(regIvsRegion,curEBISBox1.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox1.isRegular(iv))
        {
          Real val1 = a_data1(vof,comp1);
          Real val2 = a_data2(vof,comp2);

          a_volume += 1.0;
          sum += val1*val2;
        }
    }

  IntVectSet irregIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORIRREGVOFS) == EBLEVELDATAOPS_INTERIORIRREGVOFS)
    {
      Box interior = curBox1;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      irregIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYIRREGVOFS) == EBLEVELDATAOPS_BOUNDARYIRREGVOFS)
    {
      Box interior = curBox1;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox1);
      boundary -= interior;

      irregIvsRegion |= boundary;
    }

  for (VoFIterator vofit(irregIvsRegion,curEBISBox1.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox1.isIrregular(iv))
        {
          Real volfrac = curEBISBox1.volFrac(vof);
          Real val1    = a_data1(vof,comp1);
          Real val2    = a_data2(vof,comp2);

          a_volume += volfrac;
          sum += val1*val2;
        }
    }

  return sum;
}

Real EBLevelDataOps::sumKappa(Real&                a_volume,
                              const EBCellFAB&     a_data, const Box& curBox,
                              int                  a_which,
                              const ProblemDomain& a_domain)
{
  CH_assert(a_data.nComp() == 1);

  //const Box& curBox = a_data.box();
  const EBISBox& curEBISBox = a_data.getEBISBox();
  int comp = 0;

  Real sum = 0.0;

  a_volume = 0.0;

  IntVectSet regIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) == EBLEVELDATAOPS_INTERIORREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      regIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) == EBLEVELDATAOPS_BOUNDARYREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      regIvsRegion |= boundary;
    }

  for (VoFIterator vofit(regIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isRegular(iv))
        {
          Real val = a_data(vof,comp);

          a_volume += 1.0;
          sum += val;
        }
    }

  IntVectSet irregIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORIRREGVOFS) == EBLEVELDATAOPS_INTERIORIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      irregIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYIRREGVOFS) == EBLEVELDATAOPS_BOUNDARYIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      irregIvsRegion |= boundary;
    }

  for (VoFIterator vofit(irregIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isIrregular(iv))
        {
          Real volfrac = curEBISBox.volFrac(vof);
          Real val     = volfrac * a_data(vof,comp);

          a_volume += volfrac;
          sum += val;
        }
    }

  return sum;
}


Real EBLevelDataOps::sumKappaDotProductAllCells(Real&        a_volume,
                                        const EBCellFAB&     a_data1,
                                        const EBCellFAB&     a_data2,
                                        const Box& curBox,
                                        int                  a_which,
                                        const ProblemDomain& a_domain)
{
  //CP: This is a version of sumKappaDotProduct that assumes a_data1 and a_data2
  // are covering the same region, and the data in the two match up comformably
  CH_TIMERS("EBLevelDataOps::sumKappaDotProduct");
  CH_TIMER("EBLevelDataOps::sumKappaDotProduct::part1",part1);
  CH_TIMER("EBLevelDataOps::sumKappaDotProduct::part2a",part2a);
  CH_TIMER("EBLevelDataOps::sumKappaDotProduct::part2b",part2b);


  // const Box& curBox1 = a_data1.box();
  const EBISBox& curEBISBox1 = a_data1.getEBISBox();

  int ncomp = a_data1.nComp();

  Real sum = 0.0;
  // Real sum1= 0.0;

  a_volume = 0.0;
  // Real volume1 = 0.0;
  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) != EBLEVELDATAOPS_INTERIORREGVOFS)
    {
      MayDay::Error("This code has been optimized to not handle the boundary-only case");
    }


  Box region = curBox;
  if ((a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) != EBLEVELDATAOPS_BOUNDARYREGVOFS)
    {
      region.grow(1);
      region &= a_domain;
      region.grow(-1);
    }
  CH_START(part1);
  const BaseFab<Real>& regFAB1 = a_data1.getSingleValuedFAB();
  const BaseFab<Real>& regFAB2 = a_data2.getSingleValuedFAB();

  for (BoxIterator bit(region); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if (curEBISBox1.isCovered(iv))
      {
        // do nothing
      }
      else if (curEBISBox1.numVoFs(iv) == 1)
        {

          VolIndex vof(iv,0);
          Real volFrac = curEBISBox1.volFrac(vof);
          a_volume += volFrac;
          //volume1  += volFrac;
          for (int comp=0; comp<ncomp; comp++)
          {
            Real val1 = volFrac * regFAB1(iv,comp);
            Real val2 = volFrac * regFAB2(iv,comp);
            sum += val1*val2;
            // sum1+= val1*val2;
          }
        }
      // else
      //   {
      //     //multi-valued
      //     Vector<VolIndex> vofs = curEBISBox1.getVoFs(iv);
      //     for (int v=0; v<vofs.size();v++){
      //       const VolIndex& vof = vofs[v];
      //       Real volFrac = curEBISBox1.volFrac(vof);
      //       volume1 += volFrac;
      //       for (int comp=0; comp<ncomp; comp++){
      //         Real val1 = volFrac * a_data1(vof,comp);
      //         Real val2 = volFrac * a_data2(vof,comp);
      //         sum1 += val1*val2;
      //       }
      //     }
      //   }
    }
  CH_STOP(part1);
  // Having dealt with single valued cells, now we deal with multi-valued ones
  // we take advantage of the fact that (tested)
  // curEBISBox1.getVoFs(iv) is pretty fast, while a_data1(vof,comp) is agonizingly slow
  CH_START(part2a);
  const BaseIVFAB<Real>& irrBFAB1 = a_data1.getMultiValuedFAB();
  const BaseIVFAB<Real>& irrBFAB2 = a_data2.getMultiValuedFAB();
  const MiniIVFAB<Real>& irrFAB1 = static_cast<const MiniIVFAB<Real>& >(irrBFAB1);
  const MiniIVFAB<Real>& irrFAB2 = static_cast<const MiniIVFAB<Real>& >(irrBFAB2);

  int nvof1 = irrFAB1.numVoFs();
  int nvof  = irrFAB2.numVoFs();
  CH_assert(nvof1 == nvof);
  int srccomp = 0;
  int destcomp = 0;
  const Real* d1 = irrFAB1.dataPtr(destcomp);
  const Real* d2 = irrFAB2.dataPtr(srccomp);
  const Vector<VolIndex>& VoFs = irrFAB1.getVoFs();
  //const Vector<VolIndex>& VoFs2 = irrFAB2.getVoFs();

  CH_STOP(part2a);
  CH_START(part2b);
  for (int i=0; i<nvof; i++)
  //for (int i=nvof-1; i>=0; i--)
    {
      const VolIndex& vof = VoFs[i];
      Real volFrac = curEBISBox1.volFrac(vof);
      if (volFrac > 0)
        {
          const IntVect& iv = vof.gridIndex();
          if (region.contains(iv))
            {
              a_volume += volFrac;
              for (int j=0; j<ncomp; j++)
                {
                  int id = i + nvof*j;
                  sum += d1[id] * volFrac * d2[id] * volFrac;
                }
            }
        }
    }

  CH_STOP(part2b);
  return sum;
}



/***/
Real
EBLevelDataOps::sum(const LevelData<EBCellFAB> &   a_data,
                    const DisjointBoxLayout &      a_grids,
                    const EBISLayout &             a_ebisl,
                    const IntVectSet &             a_ivsExclude,
                    int   a_comp,
                    bool  a_multiplyByKappa)
{
  //sum up stuff in on this processor
  Real sumLocal = 0.0 ;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      IntVectSet ivsGrid(grid);
      ivsGrid -= a_ivsExclude;
      const EBISBox& ebisBox = a_ebisl[dit()];
      for (VoFIterator vofit(ivsGrid, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex vof = vofit();
          Real sumPt = a_data[dit()](vof, a_comp);
          if (a_multiplyByKappa)
            {
              Real kappa = ebisBox.volFrac(vof);
              sumPt *= kappa;
            }
          sumLocal += sumPt;
        }
    }

  //gather everything to proc 0
  int baseproc = 0;
  Vector<Real> sumVec;

  gather(sumVec, sumLocal, baseproc);

  //add up all the sums
  Real sumTot = 0.;
  if (procID() == baseproc)
    {
      for (int iproc = 0; iproc < sumVec.size(); iproc++)
        {
          sumTot += sumVec[iproc];
        }
    }
  //broadcast the sum to all procs
  broadcast(sumTot, baseproc);

  return sumTot;

}
/***/
Real EBLevelDataOps::sumNoKappa(Real&                a_volume,
                                const EBCellFAB&     a_data, const Box& curBox,
                                int                  a_which,
                                const ProblemDomain& a_domain)
{
  CH_assert(a_data.nComp() == 1);

  //const Box& curBox = a_data.box();
  const EBISBox& curEBISBox = a_data.getEBISBox();
  int comp = 0;

  Real sum = 0.0;

  a_volume = 0.0;

  IntVectSet regIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORREGVOFS) == EBLEVELDATAOPS_INTERIORREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      regIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYREGVOFS) == EBLEVELDATAOPS_BOUNDARYREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      regIvsRegion |= boundary;
    }

  for (VoFIterator vofit(regIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isRegular(iv))
        {
          Real val = a_data(vof,comp);

          a_volume += 1.0;
          sum += val;
        }
    }

  IntVectSet irregIvsRegion;

  if ((a_which & EBLEVELDATAOPS_INTERIORIRREGVOFS) == EBLEVELDATAOPS_INTERIORIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      irregIvsRegion |= interior;
    }

  if ((a_which & EBLEVELDATAOPS_BOUNDARYIRREGVOFS) == EBLEVELDATAOPS_BOUNDARYIRREGVOFS)
    {
      Box interior = curBox;

      interior.grow(1);
      interior &= a_domain;
      interior.grow(-1);

      IntVectSet boundary(curBox);
      boundary -= interior;

      irregIvsRegion |= boundary;
    }

  for (VoFIterator vofit(irregIvsRegion,curEBISBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      if (curEBISBox.isIrregular(iv))
        {
          Real volfrac = curEBISBox.volFrac(vof);
          Real val     = a_data(vof,comp);

          a_volume += volfrac;
          sum += val;
        }
    }

  return sum;
}
#include "NamespaceFooter.H"
