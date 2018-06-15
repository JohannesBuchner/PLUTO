#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Thurs, Aug 16, 2001

#include <fstream>
#include <string>
using std::fstream;
using std::string;
#include <cstdio>
#include <cmath>
#include "CH_HDF5.H"
#include "AMRIO.H"
#include "EBAMRIO.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "PolyGeom.H"
#include "EBAlias.H"
#include "EBIndexSpace.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "VisItChomboDriver.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5

static int g_whichCellIndex = 0;
static Real g_coveredCellValue = -98.7654321;

void
writeEBHDF5(const string& a_filename,
            const Vector<DisjointBoxLayout>& a_vectGrids,
            const Vector<LevelData<EBCellFAB>* > & a_vectData,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const Vector<int>& a_vectRatio,
            const int& a_numLevels,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect)
{
  writeEBHDF5(a_filename, a_vectGrids, &a_vectData, NULL, NULL,
              a_vectNames, a_domain, a_dx,  a_dt,  a_time,
              a_vectRatio, a_numLevels,a_replaceCovered,  a_coveredValues,
              a_ghostVect);
}


void
writeEBHDF5(const string& a_filename,
            const Vector<DisjointBoxLayout>& a_vectGrids,
            const Vector<LevelData<EBCellFAB>* > * a_phase1,
            const Vector<LevelData<EBCellFAB>* > * a_phase2,
            const Vector<LevelData<FArrayBox>* > * a_levelset,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const Vector<int>& a_vectRatio,
            const int& a_numLevels,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect)
{
  CH_TIME("EBAMRIO::writeEBHDF5");
  CH_assert(a_numLevels > 0);
  CH_assert(a_phase1->size() >= a_numLevels);
  CH_assert(a_vectRatio.size() >= a_numLevels-1);

  IntVect ghostIV = a_ghostVect;
  Vector<LevelData<FArrayBox>* > chomboData(a_numLevels, NULL);

  int ncomp = (*(*a_phase1)[0]).nComp();
  int ncomp2 = 0;

  int nnames = a_vectNames.size();

  Real dx = a_dx;

  int indexVolFrac = ncomp;

  if (a_phase2 != NULL)
  {
    ncomp2 = (*(*a_phase2)[0]).nComp();
  }

  indexVolFrac += ncomp2;
  int indexBoundaryArea = indexVolFrac+1;
  int indexAreaFrac = indexBoundaryArea+1;
  int indexNormal = indexAreaFrac+2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int indexLevelSet = indexDist + 1;

  int ncompTotal = indexDist + 1;

  if (a_phase2 != NULL)
  {
    ncompTotal = indexLevelSet+1;
  }


  CH_assert(nnames >= ncomp);

  Vector<string> names(ncompTotal);

  for (int i = 0; i < a_vectNames.size(); i++)
    {
      names[i] = a_vectNames[i];
    }

  string volFracName("fraction-0");
  string boundaryAreaName("boundaryArea-0");

  Vector<string> areaName(6);
  areaName[0] = "xAreafractionLo-0";
  areaName[1] = "xAreafractionHi-0";
  areaName[2] = "yAreafractionLo-0";
  areaName[3] = "yAreafractionHi-0";
  areaName[4] = "zAreafractionLo-0";
  areaName[5] = "zAreafractionHi-0";

  Vector<string> normName(3);
  normName[0] = "xnormal-0";
  normName[1] = "ynormal-0";
  normName[2] = "znormal-0";

  string distName("distance-0");

  names[indexVolFrac] = volFracName;
  names[indexBoundaryArea] = boundaryAreaName;

  for (int i = 0; i < 2*SpaceDim; i++)
    {
      names[indexAreaFrac+i] = areaName[i];
    }

  for (int i = 0; i < SpaceDim; i++)
    {
      names[indexNormal+i] = normName[i];
    }

  names[indexDist] = distName;

  if (a_phase2 != NULL)
    {
      names[indexLevelSet] = "LevelSet";

    }

  // set things up for each level
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      const DisjointBoxLayout& grids = a_vectGrids[ilev];
      const LevelData<EBCellFAB>& ebcfData1 = *(*a_phase1)[ilev];



      if (a_replaceCovered)
        {
          CH_assert(a_coveredValues.size() == ncomp);
        }

      //copy data into something writeAMRHierarchy can grok
      chomboData[ilev] =
        new LevelData<FArrayBox>(grids, ncompTotal, ghostIV);
      LevelData<FArrayBox>& fabData = *chomboData[ilev];

      int ibox = 0;
      // go through all the grids on this level
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          //const EBCellFAB& ebcellfab1 = ebcfData1[dit()];
          const EBCellFAB& ebcf = ebcfData1[dit()];
          EBISBox ebisBox1 = ebcf.getEBISBox();
          EBCellFAB ebcellfab1(ebisBox1, ebcf.getRegion(), ebcf.nComp());
          ebcellfab1.copy(ebcf);

          for (int icomp = 0; icomp < ebcellfab1.nComp(); icomp++)
            {
              ebcellfab1.setInvalidData(0.0, icomp);
            }
          const EBISBox& ebisbox = ebcellfab1.getEBISBox();
          FArrayBox& currentFab = fabData[dit()];
          currentFab.setVal(0.);
          // copy regular data
          currentFab.copy(ebcellfab1.getSingleValuedFAB(),0,0,ncomp);

          // copy the multicell data
          {
            IntVectSet ivsMulti = ebisBox1.getMultiCells(currentFab.box());

            for (VoFIterator vofit(ivsMulti, ebisBox1.getEBGraph()); vofit.ok(); ++vofit)
              {
                if (vofit().cellIndex() == g_whichCellIndex)
                  {
                    IntVect iv = vofit().gridIndex();
                    for (int ivar = 0; ivar < ncomp; ++ivar)
                      {
                        currentFab(iv,ivar) = ebcellfab1(vofit(),ivar);
                      }
                  }
              }
          }

          // set default volume fraction
          currentFab.setVal(1.0,indexVolFrac);

          // set default boundary area
          currentFab.setVal(0.0,indexBoundaryArea);

          // set default area fractions
          for (int i = 0; i < 2*SpaceDim; i++)
            {
              currentFab.setVal(1.0,indexAreaFrac+i);
            }

          // set default normal
          for (int i = 0; i < SpaceDim; i++)
            {
              currentFab.setVal(0.0,indexNormal+i);
            }

          // set default distance of EB from corner
          currentFab.setVal(0.0,indexDist);

          // set special values
          // iterate through the current grid
          // NOTE:  this is probably an inefficient way to do this!
          //can be diff than ebcellfab because of ghost
          //Box bregion = ebcellfab.getRegion();
          Box bregion = currentFab.box();
          bregion &= ebisbox.getDomain();
          for (BoxIterator bit(bregion); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();

              // set special values for covered cells
              if (ebisbox.isCovered(iv))
                {
                  // replace regular data if flagged
                  if (a_replaceCovered)
                    {
                      for (int icomp = 0; icomp < ncomp; icomp++)
                        {
                          Real cval = a_coveredValues[icomp];

                          currentFab(iv,icomp) = cval;
                        }
                    }

                  // volume fraction is zero
                  currentFab(iv,indexVolFrac) = 0.0;

                  // boundary area is zero
                  currentFab(iv,indexBoundaryArea) = 0.0;

                  // area fractions are zero
                  for (int i = 0; i < 2*SpaceDim; i++)
                    {
                      currentFab(iv,indexAreaFrac+i) = 0.0;
                    }
                }

              // set special values for irregular cells
              if (ebisbox.isIrregular(iv))
                {
                  Vector<VolIndex> vofs = ebisbox.getVoFs(iv);

                  //it is the last vof that ends up in the data
                  //we put the corresponding volfrac with it
                  //
                  const VolIndex& vofSpec = vofs[vofs.size()-1];
                  Real volFrac = ebisbox.volFrac(vofSpec);
                  RealVect normal = ebisbox.normal(vofSpec);
                  Real bndryArea = ebisbox.bndryArea(vofSpec);

                  if (bndryArea == 0.0)
                    {
                      if (volFrac > 0.5)
                        {
                          volFrac = 1.0;
                        }
                      else
                        {
                          volFrac = 0.0;
                        }

                      normal = RealVect::Zero;
                    }

                  // set volume fraction
                  currentFab(iv,indexVolFrac) = volFrac;

                  // set boundary area
                  currentFab(iv,indexBoundaryArea) = bndryArea;

                  // set area fractions
                  for (int i = 0; i < SpaceDim; i++)
                    {
                      Vector<FaceIndex> faces;

                      faces = ebisbox.getFaces(vofSpec,i,Side::Lo);
                      if (faces.size() == 0)
                        {
                          currentFab(iv,indexAreaFrac+2*i) = 0.0;
                        }
                      else
                        {
                          currentFab(iv,indexAreaFrac+2*i) = ebisbox.areaFrac(faces[0]);
                        }

                      faces = ebisbox.getFaces(vofSpec,i,Side::Hi);
                      if (faces.size() == 0)
                        {
                          currentFab(iv,indexAreaFrac+2*i+1) = 0.0;
                        }
                      else
                        {
                          currentFab(iv,indexAreaFrac+2*i+1) = ebisbox.areaFrac(faces[0]);
                        }
                    }

                  // set normal
                  for (int i = 0; i < SpaceDim; i++)
                    {
                      currentFab(iv,indexNormal+i) = normal[i];
                    }

                  // set distance unless the length of the normal is zero
                  Real length = PolyGeom::dot(normal,normal);

                  if (length > 0)
                    {
                      Real dist = PolyGeom::computeAlpha(volFrac,normal) * dx;
                      currentFab(iv,indexDist) = -dist;
                    }
                }
            }
          if (a_phase2 != NULL)
            {
              const LevelData<EBCellFAB>& ebcfData2 = *(*a_phase2)[ilev];
              const EBCellFAB& ebcellfab2 = ebcfData2[dit()];
              const EBISBox& ebisbox2 = ebcellfab2.getEBISBox();
              const LevelData<FArrayBox>& lls = *(*a_levelset)[ilev];
              const FArrayBox& ls = lls[dit()];

              // copy regular data
              currentFab.copy(ebcellfab2.getSingleValuedFAB(),0, ncomp, ncomp2);

              // copy the multicell data
              EBISBox ebisBox2 = ebcellfab2.getEBISBox();
              {
                IntVectSet ivsMulti = ebisBox2.getMultiCells(currentFab.box());

                for (VoFIterator vofit(ivsMulti, ebisBox2.getEBGraph()); vofit.ok(); ++vofit)
                  {
                    if (vofit().cellIndex() == g_whichCellIndex)
                      {
                        IntVect iv = vofit().gridIndex();
                        for (int ivar = 0; ivar < ncomp2; ++ivar)
                          {
                            currentFab(iv,ivar+ncomp) = ebcellfab2(vofit(),ivar);
                          }
                      }
                  }
              }

              currentFab.copy(ls, 0, indexLevelSet, 1);

              if (a_replaceCovered)
              {
                Box bregion = currentFab.box();
                bregion &= ebisbox2.getDomain();
                for (BoxIterator bit(bregion); bit.ok(); ++bit)
                  {
                    const IntVect& iv = bit();

                    // set special values for covered cells
                    if (ebisbox2.isCovered(iv))
                      {
                        // replace regular data if flagged

                        for (int icomp = 0; icomp < ncomp2; icomp++)
                          {
                            Real cval = a_coveredValues[icomp];

                            currentFab(iv,icomp+ncomp) = cval;
                          }

                      }
                  }
              }
            }
          ibox ++;
        }//end loop over boxes
      //i got sick of infs in ghost cells.
      //fabData.exchange(Interval(0, fabData.nComp()-1));

      // Only require a_numLevels-1 refine ratios - this is all that are
      // needed
      if (ilev < a_numLevels-1)
      {
        dx /= a_vectRatio[ilev];
      }
    } //end loop over levels

  // write the data
  WriteAMRHierarchyHDF5(a_filename,
                        a_vectGrids,
                        chomboData,
                        names,
                        a_domain.domainBox(),
                        a_dx,
                        a_dt,
                        a_time,
                        a_vectRatio,
                        a_numLevels);

  //clean up memory
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      delete chomboData[ilev];
    }
}

void
writeEBHDF5(const string& a_filename,
            const DisjointBoxLayout& a_grids,
            const LevelData<EBCellFAB>& a_data,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect)
{
  int numLevels = 1;
  Vector<LevelData<EBCellFAB>* > vectData(numLevels,
                                          (LevelData<EBCellFAB>*)&a_data);
  Vector<DisjointBoxLayout> vectGrids(numLevels, a_grids);
  Vector<int> refRat(numLevels,2);

  writeEBHDF5(a_filename,
              vectGrids,
              vectData,
              a_vectNames,
              a_domain,
              a_dx,
              a_dt,
              a_time,
              refRat,
              numLevels,
              a_replaceCovered,
              a_coveredValues,
              a_ghostVect);
}

class LocalEBCellFactory
  : public DataFactory<EBCellFAB>
{
public:
  /// factory function.
  /**
      Creates a new baseebcellfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual EBCellFAB* create(const Box& a_box, int a_ncomps,
                            const DataIndex& a_dit) const;

  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  LocalEBCellFactory(const EBCellFAB* a_unmanagedPtr)
  {
    CH_assert (a_unmanagedPtr != NULL);
    m_nonmanagedPtr = a_unmanagedPtr;
  }

  ///
  virtual ~LocalEBCellFactory()
  {
  }

private:

  const EBCellFAB* m_nonmanagedPtr;
};

EBCellFAB*
LocalEBCellFactory::create(const Box& a_box, int a_ncomps,
                           const DataIndex& a_dit) const
{

//  check is broken by EBCellFAB's BaseFab's no longer being
//  intersected with the domain
//  if (a_box != m_nonmanagedPtr->getRegion())
//    {
//      MayDay::Error("Incorrect use of LocalEBCellFactory, try again");
//    }
  return new EBCellFAB(m_nonmanagedPtr->getEBISBox(), a_box, a_ncomps);
}

void
writeEBHDF5(const string& a_filename,
            const Box& a_grid,
            const EBCellFAB& a_data,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect)
{
  Vector<int> procID(1, 0);
  Vector<Box> vecbox(1, a_grid);
  DisjointBoxLayout dbl(vecbox, procID);
  int ncomp = a_data.nComp();
  Interval interv(0, ncomp-1);
  LocalEBCellFactory factory(&a_data);
  LevelData<EBCellFAB> ldf(dbl, ncomp, IntVect::Zero, factory);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      ldf[dit()].copy(a_grid, interv, a_grid, a_data, interv);
    }

  writeEBHDF5(a_filename,
              dbl,
              ldf,
              a_vectNames,
              a_domain,
              a_dx,
              a_dt,
              a_time,
              a_replaceCovered,
              a_coveredValues,
              a_ghostVect);
}

// put simple debugging functions here, at the end of the hdf5 stuff
static void
ChomboVisVisualizeFile(const char *fname)
{
  char command[2000];
  sprintf(command,"$CHOMBOVIS_HOME/bin/chombovis debug_level=0 %s &",fname);
  int ret;
  ret = system(command);
}

static void
ChomboBrowserBrowseFile(const char *fname)
{
  char command[2000];
  sprintf(command,"$CHOMBOVIS_HOME/bin/chombobrowser debug_level=0 %s &",fname);
  int ret;
  ret = system(command);
}

static VisItChomboDriver visit;
static void
VisItVisualizeFile(const char *fname)
{
  visit.VisualizeFile(fname);
}

static void
VisItBrowseFile(const char *fname)
{
  visit.BrowseFile(fname);
}


static void
VisualizeFile(const char *fname)
{
  const char *use_visit = getenv("CHOMBO_USE_CHOMBOVIS");
  if (use_visit)
      ChomboVisVisualizeFile(fname);
  else
      VisItVisualizeFile(fname);
}

static void
BrowseFile(const char *fname)
{
  const char *use_visit = getenv("CHOMBO_USE_CHOMBOVIS");
  if (use_visit)
      ChomboBrowserBrowseFile(fname);
  else
      VisItBrowseFile(fname);
}

void
writeEBFAB(const EBCellFAB* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = "EBfab.hdf5";
  writeEBFABname(a_dataPtr, fname);
}

void
viewEBFAB(const EBCellFAB* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeEBFABname(a_dataPtr, fname);
  VisualizeFile(fname);

}


void
browseEBFAB(const EBCellFAB* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeEBFABname(a_dataPtr, fname);
  BrowseFile(fname);

}

void
writeEBFABname(const EBCellFAB* a_dataPtr,
               const char*      a_filename)
{
  if (a_dataPtr == NULL) return;

  const BaseFab<Real>& data = a_dataPtr->getSingleValuedFAB();

  int nComp = data.nComp();

  Vector<string> names(nComp);
  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      names[ivar] = label;
    }

  Box domain = data.box();

  // put bogus numbers here
  Real dx = 1.0;
  Real dt = 1.0;
  Real time = 1.0;

  bool replaceCovered = false;
  Vector<Real> coveredValues;

  string fileString(a_filename);

  writeEBHDF5(fileString,
              domain,
              *a_dataPtr,
              names,
              domain,
              dx,
              dt,
              time,
              replaceCovered,
              coveredValues);
}

void multiFaceValues(const EBFaceFAB* a_face,
                     const int        a_side,
                     const int        a_iv0,
                     const int        a_iv1,
                     const int        a_iv2)
{
  IntVect iv(D_DECL(a_iv0,a_iv1,a_iv2));
  EBISBox ebisBox = a_face->getEBISBox();
  Vector< FaceIndex> problemFaces;
  if (a_side == 0)
    {
      problemFaces = ebisBox.getAllFaces(iv,a_face->direction(),Side::Lo);
    }
  else if (a_side ==1)
    {
      problemFaces = ebisBox.getAllFaces(iv,a_face->direction(),Side::Hi);
    }
  else
    {
      cout<<" Third argument must be 0 for low side or 1 for high side";
      return;
    }
  int numFaces = problemFaces.size();
  cout << "In index ("<<iv[0]<<","<<iv[1]<<"), there are "<<numFaces;
  cout <<" faces on the "<<a_side<<" side, in the "<<a_face->direction()<<" direction"<<endl;

  Vector<Real> faceValue(a_face->nComp());
  for (int iFace = 0; iFace<numFaces; ++iFace)
    {
      cout<< "faceFab["<<iFace<<",*) = (";
      for (int ivar = 0;ivar <a_face->nComp() ; ++ivar)
        {
          if (ivar > 0)
            {
              cout <<", ";
            }
          cout<<(*a_face)(problemFaces[iFace],ivar);
        }
      cout <<")"<<endl;
    }
}
void multiCellValues(const EBCellFAB* a_ebCellFab,
                     const int        a_iv0,
                     const int        a_iv1,
                     const int        a_iv2)
{
  IntVect iv(D_DECL(a_iv0,a_iv1,a_iv2));
  EBISBox ebisBox = a_ebCellFab->getEBISBox();

  Vector< VolIndex> vofs = ebisBox.getVoFs(iv);

  int numVofs = vofs.size();
  cout << "In index ("<<iv[0]<<","<<iv[1]<<"), there are "<<numVofs;
  cout <<" vofs."<<endl;

  Vector<Real> vofValue(a_ebCellFab->nComp());
  for (int iVof = 0; iVof<numVofs; ++iVof)
    {
      cout<< "ebCellFab["<<vofs[iVof].cellIndex()<<",*) = (";
      for (int ivar = 0;ivar <a_ebCellFab->nComp() ; ++ivar)
        {
          if (ivar > 0)
            {
              cout <<", ";
            }
          cout<<(*a_ebCellFab)(vofs[iVof],ivar);
        }
      cout <<")"<<endl;
    }
}

void viewEBFluxLD(const LevelData<EBFluxFAB>* a_fluxLD, int a_dir)
{

  int nComp = a_fluxLD->nComp();
  const DisjointBoxLayout& dbl = a_fluxLD->disjointBoxLayout();
  LevelData<FArrayBox> fabLD(dbl, nComp);

  for (DataIterator dit = a_fluxLD->dataIterator(); dit.ok(); ++dit)
    {
      const EBFluxFAB& fluxFAB = (*a_fluxLD)[dit()];
      const EBFaceFAB& faceFAB = fluxFAB[a_dir];
      EBISBox ebisBox = faceFAB.getEBISBox();
      Box faceBox = faceFAB.getCellRegion();
      EBFaceFAB  tempFaceFAB(ebisBox, faceBox, a_dir, nComp);
      tempFaceFAB.copy(faceFAB);
      FArrayBox& tempFAB = tempFaceFAB.getFArrayBox();

      const IntVectSet& ivsMulti = ebisBox.getEBGraph().getMultiCells(faceBox);
      for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          IntVect ivHi = iv + BASISV(a_dir);
          const Vector<FaceIndex> loFaces = ebisBox.getFaces(vof, a_dir, Side::Lo);
          const Vector<FaceIndex> hiFaces = ebisBox.getFaces(vof, a_dir, Side::Hi);
          int iLoFace = 999;
          if (g_whichCellIndex + 1 <= loFaces.size())
            {
              iLoFace = g_whichCellIndex;
            }
          else
            {
              iLoFace = loFaces.size() - 1;
            }
          int iHiFace = 999;
          if (g_whichCellIndex + 1<= hiFaces.size())
            {
              iHiFace = g_whichCellIndex;
            }
          else
            {
              iHiFace = hiFaces.size() - 1;
            }
          if (loFaces.size() > 0)
            {
              for (int ivar = 0; ivar < nComp; ++ivar)
                {
                  Real loTemp = tempFaceFAB(loFaces[iLoFace], ivar);
                  tempFAB(iv, ivar) = loTemp;
                }
            }
          if (hiFaces.size() > 0)
            {
              for (int ivar = 0; ivar < nComp; ++ivar)
                {
                  Real hiTemp = tempFaceFAB(hiFaces[iHiFace], ivar);
                  tempFAB(ivHi, ivar) = hiTemp;
                }
            }
        }
      int comp = 0;
      tempFAB.shiftHalf(a_dir, -1);
      fabLD[dit()].copy(tempFAB, faceBox, comp, faceBox, comp, nComp);
    }

  viewLevel(&fabLD);
}


void writeEBFluxLDname(const LevelData<EBFluxFAB>* a_fluxLD, int a_dir, 
                   const char* a_filename)
{

  int nComp = a_fluxLD->nComp();
  const DisjointBoxLayout& dbl = a_fluxLD->disjointBoxLayout();
  LevelData<FArrayBox> fabLD(dbl, nComp);

  for (DataIterator dit = a_fluxLD->dataIterator(); dit.ok(); ++dit)
    {
      const EBFluxFAB& fluxFAB = (*a_fluxLD)[dit()];
      const EBFaceFAB& faceFAB = fluxFAB[a_dir];
      EBISBox ebisBox = faceFAB.getEBISBox();
      Box faceBox = faceFAB.getCellRegion();
      EBFaceFAB  tempFaceFAB(ebisBox, faceBox, a_dir, nComp);
      tempFaceFAB.copy(faceFAB);
      FArrayBox& tempFAB = tempFaceFAB.getFArrayBox();

      const IntVectSet& ivsMulti = ebisBox.getEBGraph().getMultiCells(faceBox);
      for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          IntVect ivHi = iv + BASISV(a_dir);
          const Vector<FaceIndex> loFaces = ebisBox.getFaces(vof, a_dir, Side::Lo);
          const Vector<FaceIndex> hiFaces = ebisBox.getFaces(vof, a_dir, Side::Hi);
          int iLoFace = 999;
          if (g_whichCellIndex + 1 <= loFaces.size())
            {
              iLoFace = g_whichCellIndex;
            }
          else
            {
              iLoFace = loFaces.size() - 1;
            }
          int iHiFace = 999;
          if (g_whichCellIndex + 1<= hiFaces.size())
            {
              iHiFace = g_whichCellIndex;
            }
          else
            {
              iHiFace = hiFaces.size() - 1;
            }
          if (loFaces.size() > 0)
            {
              for (int ivar = 0; ivar < nComp; ++ivar)
                {
                  Real loTemp = tempFaceFAB(loFaces[iLoFace], ivar);
                  tempFAB(iv, ivar) = loTemp;
                }
            }
          if (hiFaces.size() > 0)
            {
              for (int ivar = 0; ivar < nComp; ++ivar)
                {
                  Real hiTemp = tempFaceFAB(hiFaces[iHiFace], ivar);
                  tempFAB(ivHi, ivar) = hiTemp;
                }
            }
        }
      int comp = 0;
      tempFAB.shiftHalf(a_dir, -1);
      fabLD[dit()].copy(tempFAB, faceBox, comp, faceBox, comp, nComp);
    }
  writeLevelname(&fabLD, a_filename);
  // viewLevel(&fabLD);
}


void viewEBFace(const EBFaceFAB* a_face)
{
  EBISBox ebisBox = a_face->getEBISBox();
  Box box = a_face->getCellRegion();
  int faceDir = a_face->direction();
  int nComp = a_face->nComp();
  EBFaceFAB tempFaceFab(ebisBox,box, faceDir,nComp);
  tempFaceFab.copy((*a_face));
  for (int ivar = 0; ivar < tempFaceFab.nComp(); ++ivar)
    {
      tempFaceFab.setCoveredFaceVal(g_coveredCellValue,ivar);
    }
  FArrayBox & regWithMultiCells = tempFaceFab.getFArrayBox();

  const IntVectSet& ivsMulti = ebisBox.getEBGraph().getMultiCells(box);

  for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      IntVect ivHi = iv + BASISV(a_face->direction());
      const Vector<FaceIndex> loFaces = ebisBox.getFaces(vof,a_face->direction(),Side::Lo);
      const Vector<FaceIndex> hiFaces = ebisBox.getFaces(vof,a_face->direction(),Side::Hi);
      int iLoFace = 999;

      if (g_whichCellIndex + 1 <= loFaces.size())
        {
          iLoFace = g_whichCellIndex;
        }
      else
        {
          iLoFace = loFaces.size() - 1;
        }
      int iHiFace = 999;
      if (g_whichCellIndex + 1<= hiFaces.size())
        {
          iHiFace = g_whichCellIndex;
        }
      else
        {
          iHiFace = hiFaces.size() - 1;
        }
      if (loFaces.size() > 0)
        {
          for (int ivar = 0; ivar < a_face->nComp(); ++ivar)
            {
              Real loTemp = (*a_face)(loFaces[iLoFace],ivar);
              regWithMultiCells(iv,ivar) = loTemp;
            }
        }

      if (hiFaces.size() > 0)
        {
          for (int ivar = 0; ivar < a_face->nComp(); ++ivar)
            {
              Real hiTemp = (*a_face)(hiFaces[iHiFace],ivar);
              regWithMultiCells(ivHi,ivar) = hiTemp;
            }
        }
    }

  viewFAB(&regWithMultiCells);
}


void browseEBFace(const EBFaceFAB* a_face)
{
  const FArrayBox& reg = a_face->getFArrayBox();
  browseFAB(&reg);
}

void
writeEBLevel(const LevelData<EBCellFAB>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = "EBLDF.hdf5";
  writeEBLevelname(a_dataPtr, fname);
}

void
viewEBLevel(const LevelData<EBCellFAB>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeEBLevelname(a_dataPtr, fname);
  VisualizeFile(fname);

}

void
browseEBLevel(const LevelData<EBCellFAB>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeEBLevelname(a_dataPtr, fname);
  BrowseFile(fname);
}

void
writeEBLevelname(const LevelData<EBCellFAB>* a_dataPtr,
                 const char*                 a_filename)
{
  if (a_dataPtr == NULL) return;

  const LevelData<EBCellFAB>& data = *a_dataPtr;

  int nComp = data.nComp();
  IntVect ghostVect = data.ghostVect();

  Vector<string> names(nComp);
  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      names[ivar] = label;
    }

  // get the domain box
  const DisjointBoxLayout& levelBoxes = data.getBoxes();
  ProblemDomain domain = levelBoxes.physDomain();

  if (domain.domainBox().isEmpty())
  {
    IntVectSet domainIVS;
    for (DataIterator dit = data.dataIterator(); dit.ok(); ++dit)
      {
        const EBCellFAB& fab    = data[dit()];
        const EBISBox&  ebisBox = fab.getEBISBox();
        const Box& curBox = ebisBox.getDomain().domainBox();
        domainIVS |= curBox;
      }
    Box domainBox = domainIVS.minBox();
    domain.define(domainBox);
  }

  // put bogus numbers here
  Real dx = 1.0;
  Real dt = 1.0;
  Real time = 1.0;

  bool replaceCovered = false;
  Vector<Real> coveredValues;

  string fileString(a_filename);
  writeEBHDF5(fileString,
              levelBoxes,
              *a_dataPtr,
              names,
              domain,
              dx,
              dt,
              time,
              replaceCovered,
              coveredValues,
              ghostVect);
}

void
writeEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = "EBAMR.hdf5";
  writeEBAMRname(a_dataPtr, fname);
}

void
viewEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeEBAMRname(a_dataPtr, fname);
  VisualizeFile(fname);

}

void
visEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  string fname = tempnam(NULL,NULL);
  fname.append(".hdf5"); // do this so visit will recognize the file
  writeEBAMRname(a_dataPtr, fname.c_str());

  char command[2000];
  sprintf(command,"visit -o %s &",fname.c_str());
  int ret;
  ret = system(command);
}

void
browseEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeEBAMRname(a_dataPtr, fname);
  BrowseFile(fname);
}

void
writeEBAMRname(const Vector<LevelData<EBCellFAB>* >* a_dataPtr,
               const char*                           a_filename)
{
  if (a_dataPtr == NULL) return;

  const LevelData<EBCellFAB>& coarData = *((*a_dataPtr)[0]);

  int nComp = coarData.nComp();
  IntVect ghostVect = coarData.ghostVect();

  Vector<string> names(nComp);
  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      names[ivar] = label;
    }

  // create a vector of DisjointBoxLayout's
  Vector<DisjointBoxLayout> amrDBL;
  const int numLevels = a_dataPtr->size();
  for (int ilevel = 0; ilevel < numLevels; ilevel++)
    {
      amrDBL.push_back((*a_dataPtr)[ilevel]->getBoxes());
    }

  // put bogus numbers here
  Real dx = 1.0;
  Real dt = 1.0;
  Real time = 1.0;

  //extract the following info from the data
  Vector<int> refRatios(numLevels,2);
  ProblemDomain defaultDomain(Box(IntVect::Zero,7*IntVect::Unit));
  Vector<ProblemDomain> amrDomain(numLevels, defaultDomain);

  for (int ilevel = 0; ilevel < numLevels; ilevel++)
    {
      bool isDomSet = false;
      ProblemDomain domLev;
      const LevelData<EBCellFAB>& data = *((*a_dataPtr)[ilevel]);
      for (DataIterator dit = data.dataIterator(); dit.ok(); ++dit)
        {
          const EBCellFAB& fab     = data.operator[](dit);
          const EBISBox&  ebisBox = fab.getEBISBox();
          domLev = ebisBox.getDomain();
          isDomSet = true;
        }
      if (isDomSet)
      {
        amrDomain[ilevel] = domLev;
      }
    }

  for (int ilevel = 0; ilevel < numLevels-1; ilevel++)
    {
      const IntVect& coarDSize = amrDomain[ilevel].size();
      const IntVect& fineDSize = amrDomain[ilevel+1].size();
      IntVect refVect = fineDSize/coarDSize;
      CH_assert(refVect[0]==refVect[1]);
#if CH_SPACEDIM==3
      CH_assert(refVect[0]==refVect[2]);
#endif
      refRatios[ilevel] = refVect[0];
    }

  bool replaceCovered = false;
  Vector<Real> coveredValues;

  string fileString(a_filename);
  writeEBHDF5(fileString,
              amrDBL,
              *a_dataPtr,
              names,
              amrDomain[0],
              dx,
              dt,
              time,
              refRatios,
              numLevels,
              replaceCovered,
              coveredValues,
              ghostVect);
}

///
/** Writes a plotfile using the basic Chombovis format, but
    for a BaseIVFAB<Real>.  This is somewhat limited in functionality,
    since all it does is copy the BaseIVFAB to a FArrayBox in the appropriate
    cells, and then call writeFAB.  This is useful for debugging.
    *a_dataPtr is written to a file named ivfab.hdf5

*/
void
writeIVFAB(const BaseIVFAB<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = "ivfab.hdf5";
  writeIVFABname(a_dataPtr, fname);
}

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    using writeIVFAB and then running ChomboVis with a python script which
    brings up a data browser by default. The file has the same format
    as writeEBHDF5, but for a single EBCellFAB.  This is useful for
    debugging.
*/
void
viewIVFAB(const BaseIVFAB<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeIVFABname(a_dataPtr, fname);
  VisualizeFile(fname);
}


///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running browse with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single BaseIVFAB<Real>.  This is useful for debugging.
*/
void
browseIVFAB(const BaseIVFAB<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeIVFABname(a_dataPtr, fname);
  BrowseFile(fname);
}


///
/** Writes a plotfile using the same format as writeEBHDF5, but
    for a BaseIVFAB<Real>.  This is useful for debugging.
    *a_dataPtr is written to the file given by a_filename.
*/
void
writeIVFABname(const BaseIVFAB<Real>* a_dataPtr,
               const char*      a_filename)
{
  // get box from EBGraph...
  const EBGraph& ebgraph = a_dataPtr->getEBGraph();
  const IntVectSet& ivsIrreg = a_dataPtr->getIVS();

  const Box& region = ebgraph.getRegion();
  int numComp = a_dataPtr->nComp();
  FArrayBox fab(region, numComp);
  fab.setVal(0.);
  // for now, leave uninitialized, to take advantage of
  // uninitialized value in DEBUG mode. May want to revisit this later...

  // now set values in FAB to values in input.
  // this (unfortunately) will assume that there is only one interface
  // per grid cell.  However, I don't know of any other way to do this...
  VoFIterator vofit(ivsIrreg, ebgraph);
  for ( ; vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      IntVect iv = vof.gridIndex();
      for (int i=0; i<numComp; i++)
        {
          fab(iv, i) = (*a_dataPtr)(vof, i);
        }
    } // end loop over irregular vofs

  // now call normal Chombo function
  writeFABname(&fab, a_filename);
}

///
/** Write a plotfile using the same format as writeEBHDF5, but
    for a single LevelData<EBCellFAB>.  Useful for debugging.  *a_dataPtr is
    written to a file named EBLDF.hdf5.
*/
void
writeIVLevel(const LevelData<BaseIVFAB<Real> >* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = "IVLDF.hdf5";
  writeIVLevelname(a_dataPtr, fname);
}

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single LevelData<BaseIVFAB<Real> >.  This is useful for
    debugging.
*/
void
viewIVLevel(const LevelData<BaseIVFAB<Real> >* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeIVLevelname(a_dataPtr, fname);
  VisualizeFile(fname);
}

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running chombobrowser with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single LevelData<BaseIVFAB<Real> >.  This is useful for
    debugging.
*/
void
browseIVLevel(const LevelData<BaseIVFAB<Real> >* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeIVLevelname(a_dataPtr, fname);
  BrowseFile(fname);
}

///
/** Write a plotfile using the same format as writeEBHDF5, but
    for a single LevelData<BaseIVFAB<Real> >. Useful for debugging.
    *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeIVLevelname(const LevelData<BaseIVFAB<Real> >* a_dataPtr,
                    const char*                 a_filename)
{
  // first define a LevelData<FArrayBox> into which to copy
  // the values in a_dataPtr
  const DisjointBoxLayout& grids = a_dataPtr->getBoxes();
  int numComp = a_dataPtr->nComp();
  IntVect ghostVect = a_dataPtr->ghostVect();
  LevelData<FArrayBox> ldf(grids, numComp, ghostVect);
  // for now, leave uninitialized, to take advantage of
  // uninitialized value in DEBUG mode. May want to revisit this later...

  // now copy values in a_dataPtr into ldf
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const BaseIVFAB<Real>& thisIVFAB = (*a_dataPtr)[dit()];
      FArrayBox& fab = ldf[dit()];
      const EBGraph& ebgraph = thisIVFAB.getEBGraph();
      const IntVectSet& ivsIrreg = thisIVFAB.getIVS();

      // now set values in FAB to values in input.
      // this (unfortunately) will assume that there is only one interface
      // per grid cell.  However, I don't know of any other way to do this...
      VoFIterator vofit(ivsIrreg, ebgraph);
      for ( ; vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          IntVect iv = vof.gridIndex();
          for (int i=0; i<numComp; i++)
            {
              fab(iv, i) = thisIVFAB(vof, i);
            }
        } // end loop over irregular vofs

    } // end loop over grids

  writeLevelname(&ldf, a_filename);
}

//======================================================================
// new EB IO API
//======================================================================

void
createEBFile(HDF5Handle& a_handle,
           const std::string& a_filename,
           int a_numLevels,
           const Vector<int>& a_refRatios,
           const ProblemDomain& a_coarseDomain,
           const RealVect& a_coarseDx,
           const IntVect& a_ghost)
{
  RealVect origin = RealVect::Zero;
  RealVect aspect = RealVect::Unit;

  createEBFile(a_handle,
               a_filename,
               a_numLevels,
               a_refRatios,
               a_coarseDomain,
               origin,
               a_coarseDx,
               aspect,
               a_ghost);
}

void
createEBFile(HDF5Handle& a_handle,
           const std::string& a_filename,
           int a_numLevels,
           const Vector<int>& a_refRatios,
           const ProblemDomain& a_coarseDomain,
           const RealVect& a_origin,
           const RealVect& a_coarseDx,
           const RealVect& a_aspect,
           const IntVect& a_ghost)
{

  CH_assert(!a_handle.isOpen());
  CH_assert(a_numLevels > 0);

  a_handle.open(a_filename, HDF5Handle::CREATE, "Chombo_global");

  headerEBFile(a_handle,
               a_numLevels,
               a_refRatios,
               a_coarseDomain,
               a_origin,
               a_coarseDx,
               a_aspect,
               a_ghost);
}

void
headerEBFile(HDF5Handle& a_handle,
             int a_numLevels,
             const Vector<int>& a_refRatios,
             const ProblemDomain& a_coarseDomain,
             const RealVect& a_coarseDx,
             const IntVect& a_ghost)
{
  RealVect origin = RealVect::Zero;
  RealVect aspect = RealVect::Unit;

  headerEBFile(a_handle,
               a_numLevels,
               a_refRatios,
               a_coarseDomain,
               origin,
               a_coarseDx,
               aspect,
               a_ghost);
}

void
headerEBFile(HDF5Handle& a_handle,
             int a_numLevels,
             const Vector<int>& a_refRatios,
             const ProblemDomain& a_coarseDomain,
             const RealVect& a_origin,
             const RealVect& a_coarseDx,
             const RealVect& a_aspect,
             const IntVect& a_ghost)
{
  CH_assert(a_handle.isOpen());
  CH_assert(a_numLevels > 0);


  HDF5HeaderData header;

  a_handle.setGroup("/Chombo_global");
  header.m_string ["Filetype"]  = "Chombo EB File";
  header.m_int ["NumLevels"] = a_numLevels;
  header.m_box ["ProblemDomain"] = a_coarseDomain.domainBox();
  header.m_realvect ["Origin"] = a_origin;
  header.m_realvect ["Dx"] = a_coarseDx;
  header.m_realvect ["AspectRatio"] = a_aspect;
  header.m_intvect["Ghost"] = a_ghost;
  header.writeToFile(a_handle);
  header.clear();

  hid_t dataspace, dataset;

  createDataset(dataset, dataspace, a_handle, "RefRatios",
                &a_refRatios[0], a_refRatios.size());
  if (procID() == 0)
    writeDataset(dataset, dataspace, &a_refRatios[0], 0,   a_refRatios.size());

  H5Sclose(dataspace);
  H5Dclose(dataset);

  char levelName[100];
  std::string currentGroup = a_handle.getGroup();
  for (int i=0; i<a_numLevels; i++)
    {
      sprintf(levelName, "/Level%i",i);
      a_handle.setGroup(levelName);
    }

  a_handle.setGroup("/CellCenteredComponents");
  header.m_int["NumC"] = 0;
  header.writeToFile(a_handle);
  header.clear();

  a_handle.setGroup("/NodeCenteredComponents");
  header.m_int["NumN"] = 0;
  header.writeToFile(a_handle);
  header.clear();

  a_handle.setGroup("/XFaceCenteredComponents");
  header.m_int["NumX"] = 0;
  header.writeToFile(a_handle);
  header.clear();

  a_handle.setGroup("/YFaceCenteredComponents");
  header.m_int["NumY"] = 0;
  header.writeToFile(a_handle);
  header.clear();

  a_handle.setGroup("/ZFaceCenteredComponents");
  header.m_int["NumZ"] = 0;
  header.writeToFile(a_handle);
  header.clear();

  H5Fflush(a_handle.fileID(), H5F_SCOPE_LOCAL);

}

#define WRITENAMEMACRO(FUNC, GROUPPARAM, COMPSPARAM)           \
void                                                           \
write##FUNC##CenteredNames(HDF5Handle&               a_handle, \
                          const Vector<std::string>& a_names)  \
  {                                                            \
    CH_assert(a_handle.isOpen());                              \
    HDF5HeaderData header;                                     \
    a_handle.setGroup(GROUPPARAM);                             \
    header.m_int[COMPSPARAM] = a_names.size();                 \
    char cname[100];                                           \
    for (int i=0; i<a_names.size(); i++)                       \
      {                                                        \
        sprintf(cname, "Component%i",i);                       \
        header.m_string[cname] = a_names[i];                   \
      }                                                        \
    header.writeToFile(a_handle);                              \
  }

WRITENAMEMACRO(Cell,  "/CellCenteredComponents",  "NumC")
WRITENAMEMACRO(Node,  "/NodeCenteredComponents",  "NumN")
WRITENAMEMACRO(XFace, "/XFaceCenteredComponents", "NumX")
WRITENAMEMACRO(YFace, "/YFaceCenteredComponents", "NumY")
WRITENAMEMACRO(ZFace, "/ZFaceCenteredComponents", "NumZ")

class VolHDF5
{
public:
  static hid_t id;

  VolHDF5();
  VolHDF5& operator=(const VolIndex& iv);
  VolHDF5& operator=(const VolData& v);

  Real m_volFrac;
  Real m_bndryArea;
  RealVect  m_normal,  m_volCentroid,  m_bndryCentroid;
  IntVect m_cell;
  int m_ident;
};

VolHDF5::VolHDF5()
        :
    m_volFrac(0.0),
    m_bndryArea(0.0),
    m_normal(RealVect::Zero),
    m_volCentroid(RealVect::Zero),
    m_bndryCentroid(RealVect::Zero),
    m_cell(IntVect::Zero),
    m_ident(0)
{
  if (id == 0)
  {
    VolHDF5& v = *this;
    id = H5Tcreate(H5T_COMPOUND, sizeof(VolHDF5));
    H5Tinsert (id, "volFrac",   CHOFFSET(v, m_volFrac),     H5T_NATIVE_REAL);
    H5Tinsert (id, "bndryArea", CHOFFSET(v, m_bndryArea),   H5T_NATIVE_REAL);
    H5Tinsert (id, "normal",    CHOFFSET(v, m_normal),      HDF5Handle::realvect_id);
    H5Tinsert (id, "centroid",  CHOFFSET(v, m_volCentroid), HDF5Handle::realvect_id);
    H5Tinsert (id, "bndryCentroid", CHOFFSET(v, m_bndryCentroid), HDF5Handle::realvect_id);
    H5Tinsert (id, "cell",  CHOFFSET(v, m_cell),  HDF5Handle::intvect_id);
    H5Tinsert (id, "ident", CHOFFSET(v, m_ident), H5T_NATIVE_INT);
  }
}

VolHDF5& VolHDF5::operator=(const VolData& v)
{
  m_volFrac = v.m_volFrac;
  m_bndryArea = v.m_averageFace.m_bndryArea;
  m_normal = v.m_averageFace.m_normal;
  m_volCentroid = v.m_volCentroid;
  m_bndryCentroid = v.m_averageFace.m_bndryCentroid;
  return *this;
}

VolHDF5& VolHDF5::operator=(const VolIndex& iv)
{
    m_cell = iv.gridIndex();
    m_ident = iv.cellIndex();
    return *this;
}

class FaceHDF5
{
public:
  static hid_t id;

  FaceHDF5();

  Real m_areaFrac;
  RealVect m_centroid;
  IntVect m_loIv;
  int m_loVof, m_hiVof;
};

FaceHDF5::FaceHDF5()
        :
    m_areaFrac(0.0),
    m_centroid(RealVect::Zero),
    m_loIv(IntVect::Zero),
    m_loVof(0),
    m_hiVof(0)
{
  if (id == 0)
  {
    FaceHDF5& f = *this;
    id = H5Tcreate(H5T_COMPOUND, sizeof(FaceHDF5));
    H5Tinsert (id, "areaFrac", CHOFFSET(f, m_areaFrac), H5T_NATIVE_REAL);
    H5Tinsert (id, "centroid", CHOFFSET(f, m_centroid), HDF5Handle::realvect_id);
    H5Tinsert (id, "loCell",   CHOFFSET(f, m_loIv),     HDF5Handle::intvect_id);
    H5Tinsert (id, "loIdent",  CHOFFSET(f, m_loVof),    H5T_NATIVE_INT);
    H5Tinsert (id, "hiIdent",  CHOFFSET(f, m_hiVof),    H5T_NATIVE_INT);
  }
}

hid_t FaceHDF5::id = 0;
hid_t VolHDF5::id = 0;

template <class T>
void read(HDF5Handle& a_handle, const char* name, Vector<T>& v, hid_t type)
{
  hsize_t dims[1], maxdims[1];
  hid_t dataset, dataspace, memdataspace;
#ifdef H516
  dataset = H5Dopen(a_handle.groupID(), name);
#else
  dataset = H5Dopen2(a_handle.groupID(), name ,H5P_DEFAULT);
#endif
  dataspace = H5Dget_space(dataset);
  CH_assert(H5Sis_simple(dataspace));
  H5Sget_simple_extent_dims(dataspace, dims, maxdims);
  memdataspace = H5Screate_simple(1, dims, NULL);
  v.resize(dims[0]);
  H5Dread(dataset, type, memdataspace, dataspace,
          H5P_DEFAULT, &(v[0]));

  H5Sclose(memdataspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
}

void readOffsets(HDF5Handle& a_handle, int a_level,
                 Vector<Vector<long long> >& a_offsets)
{
  char levelName[100];
  sprintf(levelName, "/Level%i",a_level);
  a_handle.setGroup(levelName);

  a_offsets.resize(CH_SPACEDIM+1);
  char names[4][9] =
  {
    {
      "VOffsets"
    },
    {
      "XOffsets"
    },
    {
      "YOffsets"
    },
    {
      "ZOffsets"
    }
  };
  for (int i=0; i<1+CH_SPACEDIM; i++)
  {
    read(a_handle, names[i], a_offsets[i], H5T_NATIVE_LLONG);
  }
}

void
writeEBInfo(HDF5Handle& a_handle,
            int a_level,
            const IntVect& a_ghost,
            const LevelData<EBCellFAB>* a_data,
            Vector<Vector<long long> >& a_offsets)
{

  static VolHDF5 volType;
  static FaceHDF5 faceType;

  char levelName[100];
  sprintf(levelName, "/Level%i",a_level);
  a_handle.setGroup(levelName);

  herr_t status;
  H5E_auto_t efunc; void* edata; // turn auto error messaging off
#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  status = H5Gget_objinfo (a_handle.groupID(), "VOffsets", 0, NULL);
  H5Eset_auto(efunc, edata); //turn error messaging back on.
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  //AD : not sure if giving NULL here is OK, because this is an OUT argument
  status = H5Oget_info (a_handle.groupID(), NULL);
  H5Eset_auto2(H5E_DEFAULT,efunc, edata); //turn error messaging back on.
#endif
  if (status == 0)  //EB Data already written out.
    {
      readOffsets(a_handle, a_level, a_offsets);
      return;
    }

  OffsetBuffer boff;

  DisjointBoxLayout dbl = a_data->disjointBoxLayout();
  for (DataIterator dit = a_data->dataIterator(); dit.ok(); ++dit)
    {
      Box b = dbl.get(dit());
      const EBCellFAB& fab = a_data->operator[](dit);
      const EBGraph&   graph = fab.getEBISBox().getEBGraph();
      b.grow(a_ghost);
      boff.index.push_back(dbl.index(dit()));
      Vector<int> sizes(4,0);

      IntVectSet vofs = graph.getIrregCells(b);
      VoFIterator it(vofs, graph);
      sizes[0] = it.getVector().size();

      Vector<FaceIndex> f;
      f = graph.getIrregFaces(b, 0);
      sizes[1] = f.size();
      f.clear();

      f = graph.getIrregFaces(b, 1);
      sizes[2] = f.size();
      f.clear();

#if CH_SPACEDIM==3
      f = graph.getIrregFaces(b, 2);
      sizes[3] = f.size();
#endif

      boff.offsets.push_back(sizes);
    }

  Vector<OffsetBuffer> gathering(numProc());
  gather(gathering, boff, uniqueProc(SerialTask::compute));
  broadcast(gathering,  uniqueProc(SerialTask::compute));

  a_offsets.resize(4, dbl.size()+1);
  a_offsets[0][0] = 0;
  a_offsets[1][0] = 0;
  a_offsets[2][0] = 0;
  a_offsets[3][0] = 0;
  for (int i=0; i<numProc(); ++i)
    {
      OffsetBuffer& offbuf = gathering[i];
      for (int num=0; num<offbuf.index.size(); num++)
        {
          int index = offbuf.index[num] + 1;
          for (unsigned int j=0; j<4; ++j)
            {

              a_offsets[j][index] = offbuf.offsets[num][j];
            }
        }
    }
  for (int i=0; i<dbl.size(); i++)
    {
      for (unsigned int j=0; j<4; ++j)
        {
          a_offsets[j][i+1] += a_offsets[j][i];
        }
    }
  long long dummy = 0;
  hid_t dataset, dataspace;
  char names[4][9] =
  {
    {
      "VOffsets"
    },
    {
      "XOffsets"
    },
    {
      "YOffsets"
    },
    {
      "ZOffsets"
    }
  };
  for (int i=0; i<1+CH_SPACEDIM; i++)
  {
    createDataset(dataset, dataspace, a_handle, names[i],
                  &dummy, dbl.size()+1 );
    if (procID() == 0)
      writeDataset(dataset, dataspace, &(a_offsets[i][0]), 0,   dbl.size()+1);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  hid_t vdataset, vdataspace;
  hid_t fdataset[3], fdataspace[3];
  hsize_t nVofs = a_offsets[0].back();

  if (nVofs)
    {
    createData(vdataset, vdataspace, a_handle, std::string("VOFs"),
                  VolHDF5::id, nVofs);

    createData(fdataset[0], fdataspace[0], a_handle, "XFaces",
                  FaceHDF5::id, a_offsets[1].back());

    createData(fdataset[1], fdataspace[1], a_handle, "YFaces",
                  FaceHDF5::id, a_offsets[2].back());

#if CH_SPACEDIM == 3
    createData(fdataset[2], fdataspace[2], a_handle, "ZFaces",
               FaceHDF5::id, a_offsets[3].back());
#endif

    //  Now the actual writing of volume and face data
    for (DataIterator dit = a_data->dataIterator(); dit.ok(); ++dit)
      {
        Box b = dbl.get(dit());
        const EBCellFAB& fab = a_data->operator[](dit);
        const EBISBox& ebox = fab.getEBISBox();
        const EBGraph& graph = ebox.getEBGraph();
        b.grow(a_ghost);
        int index = dbl.index(dit());
        const BaseIVFAB<VolData>& vdata = ebox.getEBData().getVolData();
        // VoF data
        int size = a_offsets[0][index+1]-a_offsets[0][index];
        if (size > 0)
        {

          Vector<VolHDF5> buffer(size);

          IntVectSet vofs = graph.getIrregCells(b);
          VoFIterator it(vofs, graph);
          CH_assert(size == it.getVector().size());
          CH_assert(sizeof(VolData*) >= sizeof(int));

          for (int i=0; it.ok(); ++it, ++i)
            {
              VolData& v = (VolData&)vdata(it(), 0);
              buffer[i] = v;
              buffer[i] = it();

              int& volCount = (int&)(v.m_averageFace.m_volIndex);
              volCount = i;  // two phase flow hack.
                             // I know that this member isn't used.
            }
          writeDataset(vdataset, vdataspace, &(buffer[0]), a_offsets[0][index], size);

          for (int dir=0; dir<CH_SPACEDIM; ++dir)
          {
            size = a_offsets[dir+1][index+1]-a_offsets[dir+1][index];
            if (size > 0)
            {
              Vector<FaceIndex> f;
              Vector<FaceHDF5>  fbuffer(size);
              f = graph.getIrregFaces(b, dir);
              CH_assert(size == f.size());
              for (int i=0; i<size; i++)
              {
                const FaceIndex& face = f[i];
                fbuffer[i].m_areaFrac = ebox.areaFrac(face);
                fbuffer[i].m_centroid = ebox.centroid(face);
                fbuffer[i].m_loIv = face.gridIndex(Side::Lo);
                fbuffer[i].m_loVof = face.cellIndex(Side::Lo);
                fbuffer[i].m_hiVof = face.cellIndex(Side::Hi);
                // two phase flow hack, see also above.   in single phase code, these next two lines are broken.
                //fbuffer[i].m_loVof = (int)(vdata(face.getVoF(Side::Lo),0).m_averageFace.m_volIndex.cellIndex());
                //fbuffer[i].m_hiVof = (int)(vdata(face.getVoF(Side::Hi),0).m_averageFace.m_volIndex.cellIndex());
              }
              writeDataset(fdataset[dir], fdataspace[dir], &(fbuffer[0]),
                           a_offsets[dir+1][index], size);
            }
          }
        }
      }

    H5Sclose(vdataspace);
    H5Dclose(vdataset);
    for (int i=0; i<CH_SPACEDIM; i++)
      {
        H5Sclose(fdataspace[i]);
        H5Dclose(fdataset[i]);
      }
  }

  // Now generate and write Mask array
  //  0       outside problem domain
  //  1       covered
  //  2       regular
  //  3       irregular
  LevelData<BaseFab<char> > mask(dbl, 1, a_ghost);
  for (DataIterator dit = a_data->dataIterator(); dit.ok(); ++dit)
    {
      Box b = dbl.get(dit());
      const EBCellFAB& fab = a_data->operator[](dit);
      const EBISBox& ebox = fab.getEBISBox();
      const EBGraph& graph = ebox.getEBGraph();
      BaseFab<char>& m = mask[dit];
      m.setVal(0);
      graph.fillCellTypeMask(m);
    }
  mask.exchange(mask.interval());
  write(a_handle, (const BoxLayoutData<BaseFab<char> >&)mask,
        std::string("M"), a_ghost, Interval(), true);
}

void
writeCellCentered(HDF5Handle& a_handle,
                  int a_level,
                  const LevelData<EBCellFAB>* a_data,
                  Interval a_interval)
{
  CH_assert(a_handle.isOpen());

  a_handle.setGroup("/Chombo_global");
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  IntVect ghost = header.m_intvect["Ghost"];

  header.clear();

  char levelName[100];
  sprintf(levelName, "/Level%i",a_level);
  a_handle.setGroup(levelName);

  DisjointBoxLayout dbl = a_data->disjointBoxLayout();
  write(a_handle, dbl, "Boxes");

  if (a_interval.size() == 0)
    a_interval = a_data->interval();

  // write out cell centered level data
  LevelData<FArrayBox> aliasDense;
  aliasEB(aliasDense, (LevelData<EBCellFAB>&)*a_data);
  a_handle.setGroup(levelName);

  write(a_handle, (const BoxLayoutData<FArrayBox>&)aliasDense,
        std::string("C"), ghost, a_interval, true);

  a_handle.setGroup(levelName); //just to be safe

  //  OK, hardcore now, need to put out all the EB information.
  Vector<Vector<long long> >allOffsets;
  writeEBInfo(a_handle, a_level, ghost, a_data, allOffsets);

  Vector<long long>& offsets = allOffsets[0];

  //readOffsets(offsets, "VOffsets");

  hid_t dataspace, dataset;
  Real* dummy = NULL;
  int ncomps = a_interval.size();
  size_t nIrreg = offsets.back();

  if (nIrreg)
    {
    createDataset(dataset, dataspace, a_handle, "CIrregular",
                  dummy, ncomps*nIrreg);

    for (DataIterator dit = a_data->dataIterator(); dit.ok(); ++dit)
      {
        Box b = dbl.get(dit());
        const EBCellFAB& fab = a_data->operator[](dit);
        b.grow(ghost);
        int index = dbl.index(dit());
        int size = ncomps * (offsets[index+1]-offsets[index]);
        if (size > 0)
        {
          Vector<Real> buffer(size);
          int v = 0;
          const EBGraph& graph = fab.getEBISBox().getEBGraph();
          //IntVectSet ivs   = graph.getIrregCells(b);
          IntVectSet ivs = fab.getEBISBox().boundaryIVS(b);
          for (VoFIterator it(ivs, graph); it.ok(); ++it)
            {
              fab.fill(&(buffer[v]), it(), a_interval);
              v+=ncomps;
            }
          writeDataset(dataset, dataspace, &(buffer[0]),
                       ncomps*offsets[index], size);
        }
      }
    H5Sclose(dataspace);
    H5Dclose(dataset);
    }
}

void
readCellCentered(HDF5Handle& a_handle,
                 int a_level,
                 const EBIndexSpace* eb,
                 int ebghost,
                 LevelData<EBCellFAB>* a_data)
{
  a_handle.setGroup("/Chombo_global");
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  IntVect ghost = header.m_intvect["Ghost"];
  Box domain = header.m_box["ProblemDomain"];
  Vector<int> refRatios(header.m_int["NumLevels"]);
  read(a_handle, "RefRatios", refRatios, H5T_NATIVE_INT);

  for (int i=0; i<a_level; i++)
    {
      domain.refine(refRatios[i]);
    }

  int ncomp;
  a_handle.setGroup("/CellCenteredComponents");
  header.readFromFile(a_handle);
  ncomp = header.m_int["NumC"];

  Interval comps(0, ncomp-1);

  char levelName[100];
  sprintf(levelName, "/Level%i",a_level);
  a_handle.setGroup(levelName);

  Vector<Box> boxes;
  read(a_handle, boxes, "Boxes");

  Vector<int> procIDs;
  LoadBalance(procIDs, boxes);

  // pout()<<boxes<<" : "<<procIDs<<std::endl;

  DisjointBoxLayout dbl(boxes, procIDs);

  EBISLayout l;
  eb->fillEBISLayout(l, dbl, domain, ebghost);

  EBCellFactory factory(l);
  a_data->define(dbl, ncomp, ghost, factory);

  LevelData<FArrayBox> aliasDense;
  aliasEB(aliasDense, *a_data);

  {
#ifdef H516
    hid_t dataset   = H5Dopen(a_handle.groupID(), "CRegular");
#else
    hid_t dataset   = H5Dopen2(a_handle.groupID(), "CRegular" ,H5P_DEFAULT);
#endif
    hid_t dataspace = H5Dget_space(dataset);
    herr_t err;
    hsize_t count[1];
    ch_offset_t offset[1];
    Vector<Vector<char> > bufferS(1,100);
    Vector<char>& buffer = bufferS[0];
    // char dataname[9] = "COffsets";
    Vector<long long> offsets;
    read(a_handle, "COffsets", offsets, H5T_NATIVE_LLONG);
    for (DataIterator it = aliasDense.dataIterator(); it.ok(); ++it)
      {
        FArrayBox& data = aliasDense[it()];
        unsigned int index = aliasDense.boxLayout().index(it());
        Box box = aliasDense.box(it());

        offset[0] = offsets[index];
        count[0] = offsets[index+1] - offset[0];

        if (count[0] > 0)
          {
            int size = count[0] * H5Tget_size(H5T_NATIVE_REAL);
            if (size > buffer.size())
              {
                buffer.resize(size+100);
              }

            err =  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                       offset, NULL,
                                       count, NULL);
            CH_assert(err >= 0);
            hid_t memdataspace = H5Screate_simple(1, count, NULL);
            CH_assert(memdataspace >= 0);
            err = H5Dread(dataset, H5T_NATIVE_REAL, memdataspace, dataspace,
                          H5P_DEFAULT, &(buffer[0]));
            CH_assert(err >= 0);
            H5Sclose(memdataspace);

          }

        box.grow(ghost);
        read(data, bufferS,  box, comps);
      }
    H5Sclose(dataspace);
    H5Dclose(dataset);

  }

  {
    hid_t dataset, dataspace;
#ifdef H516
    dataset = H5Dopen(a_handle.groupID(), "CIrregular");
#else
    dataset = H5Dopen2(a_handle.groupID(), "CIrregular" ,H5P_DEFAULT);
#endif
    dataspace = H5Dget_space(dataset);
    CH_assert(H5Sis_simple(dataspace));

    Vector<Vector<long long> > offsets;
    readOffsets(a_handle, a_level, offsets);
    for (DataIterator dit = a_data->dataIterator(); dit.ok(); ++dit)
      {
        Box b = dbl.get(dit());
        EBCellFAB& fab = a_data->operator[](dit);
        b.grow(ghost);
        int index = dbl.index(dit());
        int size = comps.size() * (offsets[0][index+1]-offsets[0][index]);
        if (size > 0)
        {
          Vector<Real> buffer(size);
          int v = 0;
          const EBGraph& graph = fab.getEBISBox().getEBGraph();
          readDataset(dataset, dataspace, &(buffer[0]),
                      comps.size()*offsets[0][index], size);
          //IntVectSet ivs   = graph.getIrregCells(b);
          IntVectSet ivs = fab.getEBISBox().boundaryIVS(b);
          for (VoFIterator it(ivs, graph); it.ok(); ++it)
            {
              fab.assign(&(buffer[v]), it(), comps);
              v+=comps.size();
            }
        }
      }
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

}

void setWhichCellIndex(int a_whichCellIndex)
{
  g_whichCellIndex = a_whichCellIndex;
}

int getWhichCellIndex()
{
  return g_whichCellIndex;
}

void setCoveredCellValue(Real a_coveredCellValue)
{
  g_coveredCellValue = a_coveredCellValue;
}

Real getCoveredCellValue()
{
  return g_whichCellIndex;
}
#endif   // CH_USE_HDF5
#include "NamespaceFooter.H"
