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
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "IrregFAB.H"
#include "IrregFABFactory.H"
#include "EBIndexSpace.H"
#include "FArrayBox.H"
#include "LayoutIterator.H"
#include "LevelData.H"
#include "Misc.H"
#include "REAL.H"
#include "SPMD.H"
#include "SlabService.H"
#include "UGIO.H"
#include "Vector.H"

#include "GeometryShop.H"
#include "PlaneIF.H"

#include "UsingNamespace.H"

void do_geo( Box &domain, Real &dx );

int testStuff(IrregFAB& curu,
              IrregFAB& curv,
              const Box& box,
              const EBGraph&  graph)
{
  IntVectSet ivs = graph.getIrregCells(box);
  VoFIterator vofit(ivs, graph);
  Real tol = 1.0e-6;
  //first test plus 3+4 = 7
  Real rightAns = 7;
  curu.setVal(3.0);
  curv.setVal(4.0);
  curu += curv;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "fab addition operator failed at vof " << vofit().gridIndex() << endl;
          return -1;
        }
    }
  curu.setVal(3.0);
  curu += 4.0;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "real addition operator failed at vof " << vofit().gridIndex() << endl;
          return -2;
        }
    }

  //now subtraction operator 14-7 = 7
  curu.setVal(14.0);
  curv.setVal(7.0);
  curu -= curv;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "fab subtraction operator failed at vof " << vofit().gridIndex() << endl;
          return -3;
        }
    }
  curu.setVal(14.0);
  curu -= 7.0;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "real subtraction operator failed at vof " << vofit().gridIndex() << endl;
          return -4;
        }
    }

  //now multiplication. new rightans because 7 is prime (7*3 = 21)
  rightAns = 21;
  curu.setVal(7.0);
  curv.setVal(3.0);
  curu *= curv;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "fab multiplication operator failed at vof " << vofit().gridIndex() << endl;
          return -5;
        }
    }
  curu.setVal(7.0);
  curu *= 3.0;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "real multiplication operator failed at vof " << vofit().gridIndex() << endl;
          return -6;
        }
    }

  //now division (42/2 = 21)
  curu.setVal(42.0);
  curv.setVal(2.0);
  curu /= curv;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "fab division operator failed at vof " << vofit().gridIndex() << endl;
          return -7;
        }
    }
  curu.setVal(42.0);
  curu /= 2.0;
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      if (Abs(curu(vofit(), 0) - rightAns)  > tol)
        {
          pout() << "real division operator failed at vof " << vofit().gridIndex() << endl;
          return -8;
        }
    }

  return 0;
}


int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin forever present scoping trick
  {
    // Set up some geometry.
    Box domain;
    Real dx;

    int GhostThickness = 2;

    do_geo( domain, dx );

    // Make a layout for the space.
    DisjointBoxLayout dbl;
    {
      Box sbox = domain;
      Box otherbox = sbox.chop(0,5);
      pout() << "Box 1 " << sbox << endl;
      pout() << "Box 2 " << otherbox << endl;
      Vector<Box> vbox(2);
      Vector<int> vproc(2, 0);
      vbox[0] = sbox;
      vbox[1] = otherbox;

      BoxLayout bl(vbox, vproc);
      dbl.define(bl);
      dbl.close();
    }

    // Fill in the layout with the geometrical info.
    EBISLayout mylayout;
    EBIndexSpace *ebisPtr = Chombo_EBIS::instance();
    ebisPtr->fillEBISLayout( mylayout, dbl, domain, 1 );

    // Define the leveldata for my one and only level.
    IrregFABFactory ebfact(mylayout);
    LevelData<IrregFAB> Udata, Vdata;
    Udata.define(dbl, 1,  GhostThickness*IntVect::Unit, ebfact );
    Vdata.define(dbl, 1,  GhostThickness*IntVect::Unit, ebfact );

    DataIterator dit = dbl.dataIterator();
    for ( dit.begin(); dit.ok(); ++dit )
      {
        IrregFAB& curu = Udata[dit()];
        IrregFAB& curv = Vdata[dit()];
        int eekflag = testStuff(curu, curv, dbl.get(dit()), mylayout[dit()].getEBGraph());
        if (eekflag != 0)
          {
            pout() << "irregFABArith test failed with code = " << eekflag << endl;
            return eekflag;
          }
      }

    ebisPtr->clear();

  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}

void do_geo( Box &domain, Real &dx )
{
  int N = 10;
  Real xdom = 1.0;

  dx = xdom/N;
  const IntVect hi = (N-1)*IntVect::Unit;
  domain=  Box(IntVect::Zero, hi);

  int updir=1;
  int indv=0;
  Real start = 0.5;
  Real slope = tan( M_PI/6.0 );

  RealVect normal = RealVect::Zero;
  normal[updir] = 1.0;
  normal[indv] = -slope;

  RealVect point = RealVect::Zero;
  point[updir] = -slope*start;

  bool normalInside = true;

  PlaneIF myramp(normal,point,normalInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= dx;

  GeometryShop mygeom( myramp, 0, vectDx );

  RealVect origin = RealVect::Zero;

  EBIndexSpace *ebisPtr = Chombo_EBIS::instance();

  ebisPtr->define(domain, origin, dx, mygeom);

  return;
}
