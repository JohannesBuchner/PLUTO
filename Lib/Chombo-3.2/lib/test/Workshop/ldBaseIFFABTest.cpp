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
#include <cstdlib>
#include <cstring>
#include "SPACE.H"

#include "AllRegularService.H"
#include "Box.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "DataIterator.H"
#include "DebugDump.H"
#include "LoadBalance.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBFluxFAB.H"
#include "EBIndexSpace.H"

#include "Misc.H"
#include "REAL.H"
#include "SPMD.H"
#include "SlabService.H"
#include "UGIO.H"
#include "Vector.H"

#include "PlaneIF.H"
#include "GeometryShop.H"
#include "CH_Attach.H"
#include "SPMD.H"
#include "BaseIFFactory.H"
#include "EBDebugDump.H"
#include "EBArith.H"

#include "UsingNamespace.H"

bool g_diagnosticMode = false;

void do_geo(DisjointBoxLayout&              a_dbl,
            EBISLayout       &              a_ebisl,
            Box              &              a_domain,
            Real             &              a_dx )
{

  //define grids
  int len = 32;
  int mid = len/2;
  Real xdom = 1.0;
  a_dx = xdom/len;
  const IntVect hi = (len-1)*IntVect::Unit;
  a_domain=  Box(IntVect::Zero, hi);
  Vector<Box> boxes(4);
  Vector<int> procs(4);
  boxes[0] = Box(IntVect(D_DECL(0  ,  0,   0)), IntVect(D_DECL(mid-1, mid-1, len-1)));
  boxes[1] = Box(IntVect(D_DECL(0  ,mid,   0)), IntVect(D_DECL(mid-1, len-1, len-1)));
  boxes[2] = Box(IntVect(D_DECL(mid,  0,   0)), IntVect(D_DECL(len-1, mid-1, len-1)));
  boxes[3] = Box(IntVect(D_DECL(mid,mid,   0)), IntVect(D_DECL(len-1, len-1, len-1)));

  LoadBalance(procs, boxes);
  /**
  int loProc = 0;
  int hiProc = numProc() -1;
  procs[0] = loProc;
  procs[1] = hiProc;
  procs[2] = loProc;
  procs[3] = hiProc;
  **/

  a_dbl = DisjointBoxLayout(boxes, procs);

  //define geometry
  RealVect rampNormal = RealVect::Zero;
  rampNormal[0] = -0.5;
  rampNormal[1] = 0.866025404;

  Real rampAlpha = -0.0625;

  RealVect rampPoint = RealVect::Zero;
  rampPoint[0] = rampAlpha / rampNormal[0];

  bool inside = true;

  PlaneIF ramp(rampNormal,rampPoint,inside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop mygeom( ramp, 0, vectDx );
  RealVect origin = RealVect::Zero;

  EBIndexSpace *ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, origin, a_dx, mygeom);

  //fill layout
  ebisPtr->fillEBISLayout(a_ebisl, a_dbl, a_domain, 4 );
}

Real
rightAns(const FaceIndex& a_face)
{
  const IntVect& iv = a_face.gridIndex(Side::Lo);
  Real retval = 0;
  for (int  idir = 0; idir < SpaceDim; idir++)
    {
      retval += Real(iv[idir]);
    }
  return retval;
}

int
testIFFAB(const DisjointBoxLayout&              a_dbl,
          const EBISLayout       &              a_ebisl,
          const Box              &              a_domain,
          const Real             &              a_dx )
{
  int faceDir = 0;
  int nFlux = 1;
  LayoutData<IntVectSet>        irregSetsGrown;
  LevelData< BaseIFFAB<Real> > fluxInterpolant;

  EBArith::defineFluxInterpolant(fluxInterpolant,
                                 irregSetsGrown,
                                 a_dbl, a_ebisl, a_domain, nFlux, faceDir);

  //set source fab to right ans over set only on grids interior cells
  int ibox = 0;
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      BaseIFFAB<Real>& srcFab = fluxInterpolant[dit()];
      srcFab.setVal(-1.0);
      IntVectSet ivsSmall = irregSetsGrown[dit()];
      const Box& grid = a_dbl.get(dit());
      ivsSmall &=  grid;
      for (FaceIterator faceit(ivsSmall, a_ebisl[dit()].getEBGraph(), faceDir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          srcFab(faceit(), 0) = rightAns(faceit());
        }
      ibox++;
    }

  //diagnostics
  if (g_diagnosticMode)
    {
      pout() << " diagnostics for processor " << procID() << endl;
      for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
        {
          const IntVectSet& ivsGrown = irregSetsGrown[dit()];
          const Box& grid = a_dbl.get(dit());
          pout() << "============" << endl;
          pout() << " box = " << grid;
          pout() << ", full ivs = "  ;
          dumpIVS(&ivsGrown);
          pout() << "============" << endl;
          for (LayoutIterator lit = a_dbl.layoutIterator(); lit.ok(); ++lit)
            {
              const Box& grid2 = a_dbl.get(lit());
              IntVectSet ivsIntersect = ivsGrown;
              ivsIntersect &= grid2;
              pout() << "intersection with box " << grid2 << " = ";
              dumpIVS(&ivsIntersect);
              pout() << "============" << endl;
            }
        }

    }

  BaseIFFAB<Real>::setVerbose(true);
  //do the evil exchange
  Interval interv(0, nFlux-1);
  fluxInterpolant.exchange(interv);
  ibox = 0;
  //check the answer over grown set
  Real tolerance = 0.001;
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const BaseIFFAB<Real>& srcFab = fluxInterpolant[dit()];
      const IntVectSet& ivsGrown = irregSetsGrown[dit()];
      for (FaceIterator faceit(ivsGrown, a_ebisl[dit()].getEBGraph(), faceDir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          Real correct = rightAns(faceit());
          Real fabAns  =   srcFab(faceit(), 0);
          if (Abs(correct - fabAns)  > tolerance)
            {
              pout() << "iffab test failed at face "
                     << faceit().gridIndex(Side::Lo)
                     << faceit().gridIndex(Side::Hi) << endl;
              pout() << " right ans = " << correct  << endl;
              pout() << " data holds= " << fabAns << endl;

              int eekflag = -3;
              return eekflag;
            }

        }
      ibox++;
    }
  return 0;
}
int
main(int argc, char** argv)
{
  //This test is an attempt to test linearization
  //of eb data holders directly

  int eekflag;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin forever present scoping trick
  {
    // registerDebugger();

    // Set up some geometry.
    Box domain;
    Real dx;
    EBISLayout ebisl;
    DisjointBoxLayout dbl;

    do_geo( dbl, ebisl, domain, dx );


    eekflag = testIFFAB( dbl, ebisl, domain, dx );
    if (eekflag != 0)
      {
        pout() << "IFFAB exhange linearization test failed " << endl;
        return eekflag;
      }
    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
  pout() << "IFFAB exchange linearization tests passed "<< endl;
  return eekflag;
}
