#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <string>
#include "DataIterator.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "ParmParse.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "SlabService.H"
#include "BoxIterator.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "UsingNamespace.H"

#include "slab.cpp"

/***************/
int checkEBISL(const EBISLayout& a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box& a_domain,
               const Box& coveredBox);
/***************/
/***************/
int coarsenCheck(const Box& a_domain,
                 const Box& a_coveredBox);
void dumpmemoryatexit();
/***************/
/***************/
int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin scoping trick
  {

    int eekflag = 0;
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    Box domain, coveredBox;
    Real dx;
    //make the gometry.  this makes the first Chombo_EBIS
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domain,  dx, coveredBox);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    //make grids
    DisjointBoxLayout grids;
    eekflag = makeLayout(grids, domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeLayouts");
      }

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebisl;
    eekflag = makeEBISL(ebisl, grids, domain, nghost);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeEBISL");
      }

    //check everything i can think of on finest level
    eekflag = checkEBISL(ebisl, grids, domain, coveredBox);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkEBISL");
      }

    //check that coarsened things are correct as far as i can tell
    eekflag = coarsenCheck(domain, coveredBox);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in coarsenCheck");
      }

    pout() << "slab test passed" << endl;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}

/***************/
int checkEBISL(const EBISLayout& a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box& a_domain,
               const Box& a_coveredBox)
{
  int eekflag = 0;
  //check absolute basics of geometry (number of vofs, etc)
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      const Box& grid = a_grids.get(dit());
      if (!ebisBox.getRegion().contains(grid))
        {
          eekflag = 1;
          return eekflag;
        }
      if (!a_domain.contains(ebisBox.getRegion()))
        {
          eekflag = 2;
          return eekflag;
        }
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
          //check to see if it is not covered when
          //within the covered box.
          if (a_coveredBox.contains(iv))
            {
              if (!ebisBox.isCovered(iv))
                {
                  eekflag = 3;
                  return eekflag;
                }
              if (vofs.size() > 0)
                {
                  eekflag = 4;
                  return eekflag;
                }
            }
          else
            {
              if (vofs.size() != 1)
                {
                  eekflag = 5;
                  return eekflag;
                }
              const VolIndex vof = vofs[0];
              Real volFrac = ebisBox.volFrac(vof);
              if (Abs(volFrac -1.0) > 1.0e-8)
                {
                  eekflag = 6;
                  return eekflag;
                }
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      IntVect ivside = iv + sign(sit())*BASISV(idir);
                      if (a_domain.contains(ivside))
                        {
                          //                          IntVect ivdebug = a_coveredBox.smallEnd() - BASISV(1);

                          Vector<FaceIndex> faces =
                            ebisBox.getFaces(vof,idir,sit());
                          if (a_coveredBox.contains(ivside))
                            {
                              if (faces.size() != 0)
                                {
                                  eekflag = 7;
                                  return eekflag;
                                }
                           }
                          else
                            {
                              if (faces.size() != 1)
                                {
                                  eekflag = 8;
                                  return eekflag;
                                }
                              const FaceIndex& face = faces[0];
                              Real areaFrac = ebisBox.areaFrac(face);
                              if (Abs(areaFrac -1.0) > 1.0e-8)
                                {
                                  eekflag = 9;
                                  return eekflag;
                                }
                              RealVect faceCent = ebisBox.centroid(face);
                              RealVect rightCent = RealVect::Zero;
                              if (faceCent != rightCent)
                                {
                                  eekflag = 15;
                                  return eekflag;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  //check geometric derived stuff  (normals and all that).
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      const Box& grid = a_grids.get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Box sideBox =
                adjCellBox(a_coveredBox, idir, sit(), 1);
              sideBox &= grid;
              if (!sideBox.isEmpty())
                {
                  for (BoxIterator bit(sideBox); bit.ok(); ++bit)
                    {
                      const IntVect& iv = bit();
                      Vector<VolIndex> vecvof = ebisBox.getVoFs(iv);
                      CH_assert(vecvof.size() == 1);
                      const VolIndex& vof =  vecvof[0];
                      RealVect normal = ebisBox.normal(vof);
                      //normal should be unit vector in the
                      //direction outward from box
                      RealVect rightnormal = BASISREALV(idir);

                      rightnormal *= Real(sign(sit()));
                      if (normal != rightnormal)
                        {
                          eekflag = 10;
                          return eekflag;
                        }
                      //centroid should be at dead zero.
                      RealVect centroid = ebisBox.centroid(vof);
                      if (centroid != RealVect::Zero)
                        {
                          eekflag = 11;
                          return eekflag;
                        }
                      //bndryArea should be one
                      Real bndryArea = ebisBox.bndryArea(vof);
                      if (bndryArea != 1.0)
                        {
                          eekflag = 12;
                          return eekflag;
                        }
                      //bndrycentroid should be at base, in middle
                      RealVect bndryCentroid = ebisBox.bndryCentroid(vof);
                      RealVect rightCentroid = RealVect::Zero;
                      rightCentroid[idir] = 0.5*Real( sign(Side::flip(sit())));
                      if (bndryCentroid != rightCentroid)
                        {
                          eekflag = 13;
                          return eekflag;
                        }
                    } //end loop over cells in sidebox
                } //end if (!sidebox is empty)
            } //end loop over sides
        } //end loop over directions
    } //end loop over grids.
  return eekflag;
}

/***************/
/***************/
int coarsenCheck(const Box& a_domain,
                 const Box& a_coveredBox)
{
  int eekflag = 0;
  const EBIndexSpace* const ebisPtr =
    Chombo_EBIS::instance();

  int nlevels = ebisPtr->numLevels();
  //make finest layout = the layout for just the domain;
  //this will suck mightily in parallel.
  //need to do this because the DBL used elsewhere
  //might not be coarsenable to the same depth as the EBIS
  Box fineDom  = a_domain;
  DisjointBoxLayout fineDBL(Vector<Box>(1, a_domain),
                            Vector<int>(1, 0));
  EBISLayout fineEBL, finestEBL;
  ebisPtr->fillEBISLayout(fineEBL, fineDBL, fineDom, 0);
  finestEBL=fineEBL;
  int reftofinest = 1;
  for (int ilev = 0; ilev < nlevels-1; ilev++)
    {
      //create coarsened versions of layouts
      Box coarDom  = coarsen(fineDom, 2);
      DisjointBoxLayout coarDBL;
      coarsen(coarDBL, fineDBL, 2);
      EBISLayout coarEBL;
      reftofinest *= 2;
      //i do not think i need any ghost cells.
      ebisPtr->fillEBISLayout(coarEBL, coarDBL, coarDom, 0);
      //they should both work with the same iterator
      //by construction
      Real vofCoarsenFactor = 1.0;
      for (int idir = 0; idir < SpaceDim; idir++)
        vofCoarsenFactor /= 2.0;

      for (DataIterator dit = coarDBL.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& coarEBISBox = coarEBL[dit()];
          const EBISBox& fineEBISBox = fineEBL[dit()];
          const Box& fineBox = fineDBL.get(dit());
          const Box& coarBox = coarDBL.get(dit());
          if ((coarEBISBox.getRegion() != coarBox) ||
             (fineEBISBox.getRegion() != fineBox))
            {
              eekflag = 14;
              return eekflag;
            }
          if (refine(coarBox,2) != fineBox)
            {
              eekflag = 15;
              return eekflag;
            }

          for (BoxIterator bitc(coarBox); bitc.ok(); ++bitc)
            {
              const IntVect& ivc = bitc();
              //check to see that the sum of volume
              //fractions on the fine grid/factor =
              //vol frac on the coarse grid.
              Vector<VolIndex> coarVoFs= coarEBISBox.getVoFs(ivc);
              for (int ivof=0; ivof < coarVoFs.size(); ivof++)
                {
                  const VolIndex& coarVoF=coarVoFs[ivof];
                  Real coarVolFrac = coarEBISBox.volFrac(coarVoF);
                  Real finerVolFrac = 0;
                  Vector<VolIndex> finerVoFs = coarEBISBox.refine(coarVoF);
                  for (int ifine = 0; ifine < finerVoFs.size(); ifine++)
                    {
                      finerVolFrac += fineEBISBox.volFrac(finerVoFs[ifine]);
                    }
                  finerVolFrac *= vofCoarsenFactor;
                  if (coarVolFrac != finerVolFrac)
                    {
                      eekflag = 16;
                      return eekflag;
                    }
                }
              //now refine the coarse intvect to the very finest level
              //and count the number of covered and uncovered cells.
              //and come up with a new correct volume fraction.
              Real coarVolFracTot = 0.;
              for (int ivof=0; ivof < coarVoFs.size(); ivof++)
                {
                  coarVolFracTot += coarEBISBox.volFrac(coarVoFs[ivof]);
                }
              int ncovered = 0;
              int nregular = 0;
              Box finebox = refine(Box(ivc, ivc), reftofinest);
              for (BoxIterator bitf(finebox); bitf.ok(); ++bitf)
                {
                  if (a_coveredBox.contains(bitf()))
                    ncovered++;
                  else
                    nregular++;
                }
              Real fineVolFracTot = Real(nregular)/Real(ncovered+nregular);
              if (coarVolFracTot != fineVolFracTot)
                {
                  eekflag = 17;
                  return eekflag;
                }
            }
        }//end loop over grids in the layout.
      //set the finer stuff to the coarser for
      //the next go-round
      fineEBL = coarEBL;
      fineDBL = coarDBL;
      fineDom = coarDom;
    }
  return eekflag;
}
