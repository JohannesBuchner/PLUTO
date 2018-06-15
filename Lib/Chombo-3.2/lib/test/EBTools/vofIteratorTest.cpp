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
#include "VoFIterator.H"
#include "UsingNamespace.H"

#include "slab.cpp"

int checkVoFIterator(const EBISLayout& a_ebisl,
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
    makeGeometry(domain,  dx, coveredBox);

    DisjointBoxLayout grids;
    makeLayout(grids, domain);

    int nghost = 2;
    EBISLayout ebisl;
    makeEBISL(ebisl, grids, domain, nghost);


    eekflag = checkVoFIterator(ebisl, grids, domain, nghost,
                               coveredBox);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkVoFIterator");
      }

    pout() << "vofiterator test passed" << endl;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}

int checkVoFIterator(const EBISLayout& a_ebisl,
                     const DisjointBoxLayout& a_grids,
                     const Box& a_domain,
                     const int& a_nghost,
                     const Box& a_coveredBox)
{
  int eekflag = 0;

  IntVectSet ivstotal(a_domain);
  Box insideBox = grow(a_domain, -1);
  ivstotal -= insideBox;
  for (DataIterator dit=a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = a_grids.get(dit());
      IntVectSet localIVS = ivstotal & localBox;
      const EBISBox& ebisBox = a_ebisl[dit()];
      BaseIVFAB<char> gotthisone(localIVS, ebisBox.getEBGraph(), 1);
      gotthisone.setVal(false);
      for (VoFIterator vofit(localIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          gotthisone(vof, 0) = true;
          if (!a_domain.contains(iv))
            {
              eekflag = 1;
              return eekflag;
            }
          if (insideBox.contains(iv))
            {
              eekflag = 2;
              return eekflag;
            }

        }
      for (IVSIterator ivsit(localIVS);
          ivsit.ok(); ++ivsit)
        {
          const IntVect& iv = ivsit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
          for (int ivof = 0; ivof < vofs.size(); ivof++)
            {
              if (!gotthisone(vofs[ivof],0))
                {
                  eekflag = 3;
                  return eekflag;
                }
            }
        }
    }

  return eekflag;
}
