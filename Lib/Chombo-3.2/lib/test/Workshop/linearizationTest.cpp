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
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "DebugDump.H"
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

#include "GeometryShop.H"
#include "PlaneIF.H"

#include "CH_Attach.H"

#include "UsingNamespace.H"

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
Real
rightAns(const VolIndex& a_vof)
{
  const IntVect& iv = a_vof.gridIndex();
  Real retval = 0;
  for (int  idir = 0; idir < SpaceDim; idir++)
    {
      retval += Real(iv[idir]);
    }
  return retval;
}

int
testIVFAB(const EBISBox& a_ebisBox, const Box& a_box)
{
  IntVectSet ivs = a_ebisBox.getIrregIVS(a_box);
  if (ivs.isEmpty()) return 0;

  Interval comps(0,0);
  BaseIVFAB<Real> srcFab(ivs, a_ebisBox.getEBGraph(), 1);
  BaseIVFAB<Real> dstFab(ivs, a_ebisBox.getEBGraph(), 1);
  //set source fab to right ans
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      srcFab(vofit(), 0) = rightAns(vofit());
    }

  //linearize the data to dst
  int sizeFab = srcFab.size(a_box, comps);
  unsigned char* buf = new unsigned char[sizeFab];
  srcFab.linearOut(buf, a_box, comps);
  dstFab.linearIn( buf, a_box, comps);
  delete[] buf;

  //check the answer
  int eekflag = 0;
  Real tolerance = 0.001;
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      Real correct  = rightAns(vofit());
      if (Abs(dstFab(vofit(), 0) - correct) > tolerance)
        {
          pout() << "ivfab test failed at vof " << vofit().gridIndex() << endl;
          eekflag = -1;
          return eekflag;
        }
    }

  return 0;
}
int
testEBCellFAB(const EBISBox& a_ebisBox, const Box& a_box)
{
  IntVectSet ivs(a_box);

  Interval comps(0,0);
  EBCellFAB srcFab(a_ebisBox, a_box, 1);
  EBCellFAB dstFab(a_ebisBox, a_box, 1);
  //set source fab to right ans
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      srcFab(vofit(), 0) = rightAns(vofit());
    }

  //linearize the data to dst
  int sizeFab = srcFab.size(a_box, comps);
  unsigned char* buf = new unsigned char[sizeFab];
  srcFab.linearOut(buf, a_box, comps);
  dstFab.linearIn( buf, a_box, comps);
  delete[] buf;

  //check the answer
  int eekflag = 0;
  Real tolerance = 0.001;
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      Real correct  = rightAns(vofit());
      if (Abs(dstFab(vofit(), 0) - correct) > tolerance)
        {
          pout() << "ivfab test failed at vof " << vofit().gridIndex() << endl;
          eekflag = -2;
          return eekflag;
        }
    }

  return 0;
}
int
testIFFAB(const EBISBox& a_ebisBox, const Box& a_box)
{
  IntVectSet ivs = a_ebisBox.getIrregIVS(a_box);
  if (ivs.isEmpty()) return 0;

  Interval comps(0,0);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      BaseIFFAB<Real> srcFab(ivs, a_ebisBox.getEBGraph(), idir, 1);
      BaseIFFAB<Real> dstFab(ivs, a_ebisBox.getEBGraph(), idir, 1);
      //set source fab to right ans
      for (FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          srcFab(faceit(), 0) = rightAns(faceit());
        }

      //linearize the data to dst
      int sizeFab = srcFab.size(a_box, comps);
      unsigned char* buf = new unsigned char[sizeFab];
      srcFab.linearOut(buf, a_box, comps);
      dstFab.linearIn( buf, a_box, comps);
      delete[] buf;

      //check the answer
      int eekflag = 0;
      Real tolerance = 0.001;
      for (FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          Real correct  = rightAns(faceit());
          if (Abs(dstFab(faceit(), 0) - correct) > tolerance)
            {
              pout() << "ivfab test failed at face "
                     << faceit().gridIndex(Side::Lo)
                     << faceit().gridIndex(Side::Hi) << endl;

              eekflag = -3;
              return eekflag;
            }
        }
    }
  return 0;
}
int
testEBFaceFAB(const EBISBox& a_ebisBox, const Box& a_box)
{
  Interval comps(0,0);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      EBFaceFAB srcFab( a_ebisBox, a_box, idir, 1);
      EBFaceFAB dstFab( a_ebisBox, a_box, idir, 1);
      //set source fab to right ans
      IntVectSet ivs(a_box);
      for (FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          srcFab(faceit(), 0) = rightAns(faceit());
        }

      //linearize the data to dst
      int sizeFab = srcFab.size(a_box, comps);
      unsigned char* buf = new unsigned char[sizeFab];
      srcFab.linearOut(buf, a_box, comps);
      dstFab.linearIn( buf, a_box, comps);
      delete[] buf;

      //check the answer
      int eekflag = 0;
      Real tolerance = 0.001;
      for (FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          Real correct  = rightAns(faceit());
          if (Abs(dstFab(faceit(), 0) - correct) > tolerance)
            {
              pout() << "ivfab test failed at face "
                     << faceit().gridIndex(Side::Lo)
                     << faceit().gridIndex(Side::Hi) << endl;

              eekflag = -4;
              return eekflag;
            }
        }
    }
  return 0;
}
int
testEBFluxFAB(const EBISBox& a_ebisBox, const Box& a_box)
{
  Interval comps(0,0);
  IntVectSet ivs(a_box);
  EBFluxFAB srcFab( a_ebisBox, a_box, 1);
  EBFluxFAB dstFab( a_ebisBox, a_box, 1);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //set source fab to right ans
      for (FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          srcFab[idir](faceit(), 0) = rightAns(faceit());
        }
    }

  //linearize the data to dst
  int sizeFab = srcFab.size(a_box, comps);
  unsigned char* buf = new unsigned char[sizeFab];
  srcFab.linearOut(buf, a_box, comps);
  dstFab.linearIn( buf, a_box, comps);
  delete[] buf;

  //check the answer
  int eekflag = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real tolerance = 0.001;
      for (FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          Real correct  = rightAns(faceit());
          if (Abs(dstFab[idir](faceit(), 0) - correct) > tolerance)
            {
              pout() << "ivfab test failed at face "
                     << faceit().gridIndex(Side::Lo)
                     << faceit().gridIndex(Side::Hi) << endl;

              eekflag = -4;
              return eekflag;
            }
        }
    }
  return 0;
}
int
main(int argc, char** argv)
{
  //This test is an attempt to test linearization
  //of eb data holders directly

  int eekflag=0;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin forever present scoping trick
  {
    //registerDebugger();

    // Set up some geometry.
    Box domain;
    Real dx;

    do_geo( domain, dx );

    Vector<Box> boxes(1, domain);
    Vector<int> procs(1, 0);
    // Make a layout for the space.
    DisjointBoxLayout dbl(boxes, procs);

    // Fill in the layout with the geometrical info.
    EBISLayout ebisl;
    EBIndexSpace *ebisPtr = Chombo_EBIS::instance();
    ebisPtr->fillEBISLayout(ebisl, dbl, domain, 1 );

    // Define the leveldata for my one and only level.
    DataIterator dit = dbl.dataIterator();
    for ( dit.begin(); dit.ok(); ++dit )
      {
        eekflag = testIVFAB(ebisl[dit()], dbl.get(dit()));
        if (eekflag != 0)
          {
            pout() << "IVFAB linearization test failed " << endl;
            return eekflag;
          }
        eekflag = testIFFAB(ebisl[dit()], dbl.get(dit()));
        if (eekflag != 0)
          {
            pout() << "IFFAB  linearization test failed " << endl;
            return eekflag;
          }
        eekflag = testEBCellFAB(ebisl[dit()], dbl.get(dit()));
        if (eekflag != 0)
          {
            pout() << "EBCellFAB  linearization test failed " << endl;
            return eekflag;
          }
        eekflag = testEBFaceFAB(ebisl[dit()], dbl.get(dit()));
        if (eekflag != 0)
          {
            pout() << "EBFaceFAB  linearization test failed " << endl;
            return eekflag;
          }

        eekflag = testEBFluxFAB(ebisl[dit()], dbl.get(dit()));
        if (eekflag != 0)
          {
            pout() << "EBFluxFAB  linearization test failed " << endl;
            return eekflag;
          }
      }

    ebisPtr->clear();

  }//end scoping trick
  pout() << "linearization tests passed "<< endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}
