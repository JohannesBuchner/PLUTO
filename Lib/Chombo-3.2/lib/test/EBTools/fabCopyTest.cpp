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
#include "BaseIVFAB.H"
#include "BaseIFFAB.H"
#include "BaseEBCellFAB.H"
#include "BaseEBFaceFAB.H"
#include "VoFIterator.H"
#include "UsingNamespace.H"


#include "slab.cpp"

/***************/
/***************/
int checkIrregFabCopy(const EBISLayout& a_ebisl,
                      const DisjointBoxLayout& a_grids,
                      const Box& a_domain,
                      const int& nghost,
                      const Box& coveredBox);
/***************/
/***************/
int checkRegularFabCopy(const EBISLayout& a_ebisl,
                        const DisjointBoxLayout& a_grids,
                        const Box& a_domain,
                        const int& nghost,
                        const Box& coveredBox);
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

    DisjointBoxLayout grids;
    eekflag = makeLayout(grids, domain);

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebisl;
    eekflag = makeEBISL(ebisl, grids, domain, nghost);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeEBISL");
      }

    eekflag = checkIrregFabCopy(ebisl, grids, domain, nghost,
                                coveredBox);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkIrregFabCopy");
      }

    eekflag = checkRegularFabCopy(ebisl, grids, domain, nghost,
                                  coveredBox);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkRegularFabCopy");
      }

    pout() << "fab copy test passed" << endl;
  }//end scoping trick


  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
/***************/

/***************/
int checkIrregFabCopy(const EBISLayout& a_ebisl,
                      const DisjointBoxLayout& a_grids,
                      const Box& a_domain,
                      const int& a_nghost,
                      const Box& a_coveredBox)
{
  int eekflag = 0;

  //  Box outsideBox = grow(a_coveredBox, 1);
  Box outsideBox = a_domain;
  Box interiorBox1 = grow(a_domain, -1);
  //IntVectSet ivstotal(outsideBox);
  //ivstotal -= interiorBox1;
  Box interiorBox2 = grow(a_coveredBox, 1);
  //IntVectSet innerIVS(interiorBox2);
  //innerIVS -= a_coveredBox;
  //ivstotal |= innerIVS;
  Interval interv(0, 0);
  for (DataIterator dit=a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = a_grids.get(dit());
      Box subBox = localBox & interiorBox1;
      IntVectSet subIVS(subBox);
      subIVS -= interiorBox2;
      IntVectSet localIVS(localBox);
      localIVS -= a_coveredBox;
      localIVS -= subIVS;
      //IntVectSet localIVS2 = ivstotal & localBox;
      const EBISBox& ebisBox = a_ebisl[dit()];
      BaseIVFAB<int> ivfabone(localIVS, ebisBox.getEBGraph(), 1);
      BaseIVFAB<int> ivfabtwo(localIVS, ebisBox.getEBGraph(), 1);
      ivfabone.setVal(1);
      ivfabtwo.setVal(2);
      ivfabtwo.copy(localBox,  interv, localBox, ivfabone, interv);
      for (VoFIterator vofit(localIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          if (ivfabtwo(vof, 0) != 1)
            {
              eekflag = 1;
              return eekflag;
            }
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          BaseIFFAB<int> iffabone(localIVS, ebisBox.getEBGraph(), idir, 1);
          BaseIFFAB<int> iffabtwo(localIVS, ebisBox.getEBGraph(), idir, 1);
          iffabone.setVal(1);
          iffabtwo.setVal(2);
          iffabtwo.copy(localBox,  interv, localBox, iffabone, interv);
          for (VoFIterator vofit(localIVS, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Vector<FaceIndex> faces = ebisBox.getFaces(vof, idir, sit());
                  for (int iface = 0; iface < faces.size(); ++iface)
                    {
                      const FaceIndex& face = faces[iface];
                      if (iffabtwo(face, 0) != 1)
                        {
                          eekflag = 2;
                          return eekflag;
                        }
                    }
                }
            }
        }
    }
  return eekflag;
}

/***************/
/***************/
int checkRegularFabCopy(const EBISLayout& a_ebisl,
                        const DisjointBoxLayout& a_grids,
                        const Box& a_domain,
                        const int& a_nghost,
                        const Box& a_coveredBox)
{
  int eekflag = 0;

  Interval interv(0, 0);
  //  Box interiorBox = grow(a_domain, -1);
  for (DataIterator dit=a_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box localBox = a_grids.get(dit());
      //localBox &= interiorBox;
      const EBISBox& ebisBox = a_ebisl[dit()];
      BaseEBCellFAB<int> ivfabone(ebisBox, localBox, 1);
      BaseEBCellFAB<int> ivfabtwo(ebisBox, localBox, 1);
      ivfabone.setVal(1);
      ivfabtwo.setVal(2);
      ivfabtwo.copy(localBox,  interv, localBox, ivfabone, interv);
      IntVectSet localIVS(localBox);
      for (VoFIterator vofit(localIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          if (ivfabtwo(vof, 0) != 1)
            {
              eekflag = 3;
              return eekflag;
            }
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          BaseEBFaceFAB<int> iffabone(ebisBox, localBox, idir, 1);
          BaseEBFaceFAB<int> iffabtwo(ebisBox, localBox, idir, 1);
          iffabone.setVal(1);
          iffabtwo.setVal(2);
          iffabtwo.copy(localBox,  interv, localBox, iffabone, interv);
          for (VoFIterator vofit(localIVS, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Vector<FaceIndex> faces = ebisBox.getFaces(vof, idir, sit());
                  for (int iface = 0; iface < faces.size(); ++iface)
                    {
                      const FaceIndex& face = faces[iface];
                      if (iffabtwo(face, 0) != 1)
                        {
                          eekflag = 4;
                          return eekflag;
                        }
                    }
                }
            }
        }
    }
  return eekflag;
}
