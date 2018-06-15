#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "VoFIterator.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBArith.H"

#include "UsingNamespace.H"

const bool g_verbose = true;
/***************/
// define a ramp EBIS.
/***************/
int makeGeometry(Box& domain);
/***************/
int checkIO(const Box& a_domain, const Box& a_startdom);
/***************/
void dumpmemoryatexit();
/***************/
int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin forever present scoping trick
  {
    const char* in_file = "ramp.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    int eekflag = 0;
    //define the geometry object.
    Box domain;
    eekflag =  makeGeometry(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    //check that total coarse vol = total fine vol
    //check that total coarse vol centroid = total fine vol centroid
    Box startdom = domain;
    bool keepgoing = true;
    while (keepgoing)
      {
        barrier();
        eekflag = checkIO(domain, startdom);
        if (eekflag != 0)
          {
            pout() << "checkCoarse: eek = " << eekflag << endl;
            MayDay::Error("problem in checkEBISL");
          }
        startdom.coarsen(2);
        pout() << "startdom = " << startdom << endl;
        keepgoing = (startdom.size(0) >= 16);
        
        EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
        ebisPtr->clear();

        if (keepgoing)
          {
            eekflag =  makeGeometry(domain);
            if (eekflag != 0)
              {
                pout() << "non zero eek detected = " << eekflag << endl;
                MayDay::Error("problem in makeGeometry");
              }
          }
      }
    pout() << "ebisWriteAll test passed" << endl;
  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
/**********/
int makeGeometry(Box& a_domain)
{
  Real dx;
  RealVect origin;
  int upDir;
  int indepVar;
  Real startPt;
  Real slope;
  int eekflag =  0;
  //parse input file
  ParmParse pp;

  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  dx = (prob_hi-prob_lo[0])/n_cell[0];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = prob_lo[idir];
    }
  pp.get("up_dir",upDir);
  pp.get("indep_var",indepVar);
  pp.get("start_pt", startPt);
  pp.get("ramp_slope", slope);

  RealVect normal = RealVect::Zero;
  normal[upDir] = 1.0;
  normal[indepVar] = -slope;

  RealVect point = RealVect::Zero;
  point[upDir] = -slope*startPt;

  bool normalInside = true;

  PlaneIF ramp(normal,point,normalInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= dx;

  GeometryShop workshop(ramp,0,vectDx);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, origin, dx, workshop);

  return eekflag;
}
/***************/
int checkIO(const Box& a_domain, const Box& a_startdom)
{
#ifdef CH_USE_HDF5
  Real tolerance  = PolyGeom::getTolerance();
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  Vector<Box> vbox(1, a_domain);
  Vector<int>  procAssign(1, 0);
  DisjointBoxLayout dblFinest(vbox, procAssign);
  int numLevels = ebisPtr->numLevels();

  int startlev = -1;
  bool found = false;
  Box domlev = a_domain;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      if (a_startdom == domlev)
        {
          found = true;
          startlev = ilev;
          break;
        }
      domlev.coarsen(2);
    }
  if (!found)
    {
      MayDay::Error("did not find starting domain");
    }
  //create all the layouts
  Vector<DisjointBoxLayout> dblayouts(numLevels);
  Vector<Box>               domains(numLevels);
  dblayouts[0] = dblFinest;
  domains[0] = a_domain;
  for (int ilev = 1; ilev < numLevels; ilev++)
    {
      ebcoarsen(dblayouts[ilev], dblayouts[ilev-1], 2);
      domains[ilev] = ebcoarsen(domains[ilev-1], 2);
    }

  Vector<EBISLayout> ebislOutput(numLevels);
  Vector<EBISLayout> ebislInput (numLevels);
  //define a bunch of ebisls to represent ebis before output
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      ebisPtr->fillEBISLayout(ebislOutput[ilev], dblayouts[ilev], domains[ilev], 0);
    }
  //output original ebis to a file.
  const char* const  filename = "ebiswriteall.hdf5";
  HDF5Handle handleOut(filename, HDF5Handle::CREATE);
  ebisPtr->writeAllLevels(handleOut);
  handleOut.close();

  //define ebis anew from file output
  HDF5Handle handleIn(filename, HDF5Handle::OPEN_RDONLY);
  ebisPtr->readInAllLevels(handleIn, a_startdom);
  handleIn.close();

  //define a bunch of ebisls to represent ebis after input
  for (int ilev = startlev; ilev < numLevels; ilev++)
    {
      ebisPtr->fillEBISLayout(ebislInput[ilev], dblayouts[ilev], domains[ilev], 0);
    }

  //now compare the damn things.  this test might break horribly if
  //we ever do the "data base of EBISLayouts" idea
  for (int ilev = startlev; ilev < numLevels; ilev++)
    {
      const EBISLayout& ebislOld = ebislOutput[ilev];
      const EBISLayout& ebislNew = ebislInput[ilev];
      const DisjointBoxLayout& dbl = dblayouts[ilev];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBoxOld = ebislOld[dit()];
          const EBISBox& ebisBoxNew = ebislNew[dit()];
          const Box& box = dbl.get(dit());
          IntVectSet ivsOld = ebisBoxOld.getIrregIVS(box);
          IntVectSet ivsNew = ebisBoxNew.getIrregIVS(box);
          IntVectSet ivsSumOfDiffs = ((ivsOld - ivsNew) |
                                      (ivsNew - ivsOld));
          if (!ivsSumOfDiffs.isEmpty())
            {
              pout() << "irreg IVS does not match "  << endl;
              return -113;
            }
          ivsOld = ebisBoxOld.getMultiCells(box);
          ivsNew = ebisBoxNew.getMultiCells(box);
          ivsSumOfDiffs = ((ivsOld - ivsNew) |
                           (ivsNew - ivsOld));
          if (!ivsSumOfDiffs.isEmpty())
            {
              pout() << "multi IVS does not match "  << endl;
              return -115;
            }
          for (BoxIterator bit(box); bit.ok(); ++bit)
            {
              Vector<VolIndex> vofsOld = ebisBoxOld.getVoFs(bit());
              Vector<VolIndex> vofsNew = ebisBoxNew.getVoFs(bit());
              if (vofsOld.size() != vofsNew.size())
                {
                  pout() << "vof size does not match at " << bit() << endl;
                  return -1;
                }
              for (int ivof = 0; ivof < vofsOld.size(); ivof++)
                {
                  //now check all the cell-centered stuff
                  const VolIndex& vofOld = vofsOld[ivof];
                  const VolIndex& vofNew = vofsNew[ivof];
                  if (vofOld != vofNew)
                    {
                      pout() << "vofs do not match at " << bit() << endl;
                      return -2;
                    }
                  Real volFracOld = ebisBoxOld.volFrac(vofOld);
                  Real volFracNew = ebisBoxNew.volFrac(vofNew);
                  if (Abs(volFracOld - volFracNew) > tolerance)
                    {
                      pout() << "vol fracs do not match at " << bit() << endl;
                      return -3;
                    }
                  Real bndryAreaOld = ebisBoxOld.bndryArea(vofOld);
                  Real bndryAreaNew = ebisBoxNew.bndryArea(vofNew);
                  if (Abs(bndryAreaOld - bndryAreaNew) > tolerance)
                    {
                      pout() << "bndry areas do not match at " << bit() << endl;
                      return -4;
                    }
                  RealVect normalOld = ebisBoxOld.normal(vofOld);
                  RealVect normalNew = ebisBoxNew.normal(vofNew);
                  RealVect volCentroidOld = ebisBoxOld.centroid(vofOld);
                  RealVect volCentroidNew = ebisBoxNew.centroid(vofNew);
                  RealVect bndryCentroidOld = ebisBoxOld.bndryCentroid(vofOld);
                  RealVect bndryCentroidNew = ebisBoxNew.bndryCentroid(vofNew);
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      if (Abs(bndryCentroidOld[idir] - bndryCentroidNew[idir]) > tolerance)
                        {
                          pout() << "bndry centroids do not match at iv dir " << bit()
                                 << " " << idir << endl;
                          return -4;
                        }
                      if (Abs(volCentroidOld[idir] - volCentroidNew[idir]) > tolerance)
                        {
                          pout() << "vol centroids do not match at iv dir " << bit()
                                 << " " << idir << endl;
                          return -5;
                        }
                      if (Abs(normalOld[idir] - normalNew[idir]) > tolerance)
                        {
                          pout() << "normals do not match at iv dir " << bit()
                                 << " " << idir << endl;
                          return -6;
                        }
                    }

                  if (ilev > 0)
                    {
                      VolIndex coarVoFOld = ebisBoxOld.coarsen(vofOld);
                      VolIndex coarVoFNew = ebisBoxNew.coarsen(vofNew);
                      if (coarVoFOld != coarVoFNew)
                        {
                          pout() << "coarse vofs do not match at " << bit() << endl;
                          return -7;
                        }
                    }
                  if (ilev < numLevels - 1)
                    {
                      Vector<VolIndex> fineVoFsOld = ebisBoxOld.refine(vofOld);
                      Vector<VolIndex> fineVoFsNew = ebisBoxNew.refine(vofNew);
                      if (fineVoFsOld.size() != fineVoFsNew.size())
                        {
                          pout() << "fine vof sizes do not match at " << bit() << endl;
                          return -8;
                        }
                      for (int ifine = 0; ifine < fineVoFsOld.size();  ifine++)
                        {
                          const VolIndex& fineVoFOld = fineVoFsOld[ifine];
                          const VolIndex& fineVoFNew = fineVoFsNew[ifine];
                          if (fineVoFOld != fineVoFNew)
                            {
                              pout() << "fine vofs do not match at " << bit() << endl;
                              return -9;
                            }
                        }
                    }
                  //now check all the face-centered stuff
                  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
                    {
                      for (SideIterator sit; sit.ok(); ++sit)
                        {
                          Vector<FaceIndex> facesOld = ebisBoxOld.getFaces(vofOld, faceDir, sit());
                          Vector<FaceIndex> facesNew = ebisBoxNew.getFaces(vofNew, faceDir, sit());
                          if (facesOld.size() != facesNew.size())
                            {
                              pout() << "fine face sizes do not match at " << bit() << endl;
                              pout() << "facedir, side = " << faceDir << "  " << sign(sit()) << endl;
                              return -10;
                            }
                          for (int iface = 0; iface < facesOld.size(); iface++)
                            {
                              const FaceIndex& faceOld = facesOld[iface];
                              const FaceIndex& faceNew = facesNew[iface];
                              IntVect ivdeblo(D_DECL(31,12,0));
                              IntVect ivdebhi(D_DECL(32,12,0));
                              int stopHere = -1;
                              const IntVect& ivlo = faceOld.gridIndex(Side::Lo);
                              const IntVect& ivhi = faceOld.gridIndex(Side::Hi);
                              if ((ilev == 1) && (ivlo == ivdeblo) && (ivhi == ivdebhi))
                                {
                                  stopHere = 1;
                                }
                              if (faceOld != faceNew)
                                {
                                  pout() << "faces do not match at " << bit() << endl;
                                  pout() << "facedir, side = " << faceDir << "  " << sign(sit()) << endl;
                                  return -11;
                                }
                              Real areaFracOld = ebisBoxOld.areaFrac(faceOld);
                              Real areaFracNew = ebisBoxNew.areaFrac(faceNew);
                              if (Abs(areaFracOld - areaFracNew) > tolerance)
                                {
                                  pout() << "area fracs do not match at " << bit() << endl;
                                  pout() << "facedir, side = " << faceDir << "  " << sign(sit()) << endl;
                                  return -12;
                                }

                              RealVect faceCentroidOld = ebisBoxOld.centroid(faceOld);
                              RealVect faceCentroidNew = ebisBoxNew.centroid(faceNew);
                              for (int idir = 0; idir < SpaceDim; idir++)
                                {
                                  if (idir != faceDir)
                                    {
                                      if (Abs(faceCentroidOld[idir] - faceCentroidNew[idir]) > tolerance)
                                        {
                                          pout() << "face centroids do not match at " << bit() << endl;
                                          pout() << "faceCentroidOld = " << faceCentroidOld << endl;
                                          pout() << "faceCentroidNew = " << faceCentroidNew << endl;
                                          pout() << "facedir, side = " << faceDir << "  " << sign(sit()) << endl;
                                          return -13;
                                        }
                                    }
                                }
                            } //end loop over faces
                        } //end loop over sides for faces
                    } //end loop over face directions
                } //end loop over vofs at a cell
            } //end loop over cells in a box
        } //end loop over boxes
    }//end loop over refinement levels
#endif // CH_USE_HDF5
  return 0;
}
/***************/
/***************/
