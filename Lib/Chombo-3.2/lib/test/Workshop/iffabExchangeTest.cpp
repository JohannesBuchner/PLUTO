#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <stdio.h>

#include "ParmParse.H"
#include "EBArith.H"
#include "EBAMRIO.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "ParmParse.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "VoFIterator.H"
#include "TiltedCylinderIF.H"
#include  "RealVect.H"
#include "Box.H"
#include "ParmParse.H"
#include "EBIndexSpace.H"
#include "Vector.H"
#include "AllRegularService.H"
#include "SlabService.H"
#include "GeometryShop.H"
#include "TiltedCylinderIF.H"
#include "EllipsoidIF.H"
#include "TransformIF.H"
#include "BRMeshRefine.H"
#include "parstream.H"
#include "PolyGeom.H"
#include "BaseIFFactory.H"
#include "BaseIFFactory.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "UsingNamespace.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/***************/
void getDebugIVS(IntVectSet& a_ivs,
                 const bool& a_thisCalcCoar)
{
  a_ivs = IntVectSet();
  if (!a_thisCalcCoar)
    {
      IntVect ivdebugloFine(D_DECL(66,10,10));
      IntVect ivdebughiFine(D_DECL(68,11,10));
      a_ivs |= ivdebugloFine;
      a_ivs |= ivdebughiFine;
    }
}

/***************/
void
getFinestDomain(Box&       a_domain,
                Real&      a_dx)
{
  ParmParse pp;
  Vector<int> n_cell(SpaceDim, 128);
  //  pp.getarr("n_cell",n_cell,0,SpaceDim);

 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain = Box(lo, hi);

  Real prob_hi = 1.0;
  int numOpen = n_cell[0];
  a_dx = prob_hi/numOpen;
}
/***********/
void
makeGeometry(const Box&       a_domain,
             const Real&      a_dx)
{
  pout() << "ellipsoid geometry" << endl;
  bool insideCalc = false;
  RealVect ellipsoidCenter(D_DECL(0.36, 0.50, 0.5));
  RealVect ellipsoidRadii( D_DECL(0.18, 0.25, 0.3));
  RealVect ellipsoidXAxis( D_DECL(3.00, 2.00, 1.0));

  Real sum;
  PolyGeom::unifyVector(ellipsoidXAxis, sum);
  EllipsoidIF ellipsoid(ellipsoidRadii, ellipsoidCenter, insideCalc);
  RealVect origxaxis = BASISREALV(0);
  TransformIF rotazoid(ellipsoid);
  rotazoid.rotate(origxaxis, ellipsoidXAxis, ellipsoidCenter);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(rotazoid,0,vectDx);
  //this generates the new EBIS
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  RealVect origin = RealVect::Zero;
  int ebMaxSize = 16;
  int ebMaxCoarsen = 0;
  ebisPtr->define(a_domain, origin, a_dx, workshop, ebMaxSize, ebMaxCoarsen);
}
/************/
void
makeDBL(DisjointBoxLayout& a_dbl)
{
  Vector<Box> boxes;
#if CH_SPACEDIM==3
  boxes.push_back(Box(IntVect(8,48,56)  ,  IntVect(19,63,71) ));
  boxes.push_back(Box(IntVect(8,64,56)  ,  IntVect(19,79,71) ));
  boxes.push_back(Box(IntVect(8,80,24)  ,  IntVect(19,87,39) ));
  boxes.push_back(Box(IntVect(12,40,48) ,  IntVect(19,55,55) ));
  boxes.push_back(Box(IntVect(12,44,88) ,  IntVect(23,55,95) ));
  boxes.push_back(Box(IntVect(12,68,96) ,  IntVect(19,79,103)));
  boxes.push_back(Box(IntVect(16,40,56) ,  IntVect(19,47,59) ));
  boxes.push_back(Box(IntVect(16,52,100),  IntVect(19,55,103)));
  boxes.push_back(Box(IntVect(16,88,36) ,  IntVect(23,91,47) ));
  boxes.push_back(Box(IntVect(16,96,56) ,  IntVect(23,103,63)));
  boxes.push_back(Box(IntVect(20,36,52) ,  IntVect(23,39,55) ));
  boxes.push_back(Box(IntVect(20,36,88) ,  IntVect(23,43,95) ));
  boxes.push_back(Box(IntVect(20,48,64) ,  IntVect(27,51,71) ));
  boxes.push_back(Box(IntVect(20,52,64) ,  IntVect(23,59,71) ));
  boxes.push_back(Box(IntVect(20,56,56) ,  IntVect(23,63,63) ));
  boxes.push_back(Box(IntVect(20,68,36) ,  IntVect(31,79,47) ));
  boxes.push_back(Box(IntVect(20,80,24) ,  IntVect(31,87,39) ));
  boxes.push_back(Box(IntVect(20,88,80) ,  IntVect(31,103,87)));
  boxes.push_back(Box(IntVect(24,32,48) ,  IntVect(31,39,55) ));
  boxes.push_back(Box(IntVect(24,44,88) ,  IntVect(39,55,95) ));
  boxes.push_back(Box(IntVect(24,56,88) ,  IntVect(39,71,95) ));
  boxes.push_back(Box(IntVect(24,80,96) ,  IntVect(31,87,111)));
  boxes.push_back(Box(IntVect(24,88,64) ,  IntVect(31,95,79) ));
  boxes.push_back(Box(IntVect(24,96,88) ,  IntVect(27,99,91) ));
  boxes.push_back(Box(IntVect(28,88,28) ,  IntVect(31,95,35) ));
  boxes.push_back(Box(IntVect(32,24,64) ,  IntVect(39,35,79) ));
  boxes.push_back(Box(IntVect(32,32,40) ,  IntVect(39,43,47) ));
  boxes.push_back(Box(IntVect(32,40,80) ,  IntVect(39,51,87) ));
  boxes.push_back(Box(IntVect(32,48,108),  IntVect(35,55,111)));
  boxes.push_back(Box(IntVect(32,60,32) ,  IntVect(39,71,47) ));
  boxes.push_back(Box(IntVect(32,72,104),  IntVect(39,87,115)));
  boxes.push_back(Box(IntVect(32,84,24) ,  IntVect(39,95,35) ));
  boxes.push_back(Box(IntVect(32,92,80) ,  IntVect(39,107,87)));
  boxes.push_back(Box(IntVect(32,96,92) ,  IntVect(39,99,99) ));
  boxes.push_back(Box(IntVect(40,0,8)   ,  IntVect(55,15,19) ));
  boxes.push_back(Box(IntVect(40,0,52)  ,  IntVect(55,15,63) ));
  boxes.push_back(Box(IntVect(40,0,96)  ,  IntVect(55,15,107)));
  boxes.push_back(Box(IntVect(40,16,8)  ,  IntVect(55,31,19) ));
  boxes.push_back(Box(IntVect(40,16,52) ,  IntVect(55,31,63) ));
  boxes.push_back(Box(IntVect(40,16,96) ,  IntVect(55,31,107)));
  boxes.push_back(Box(IntVect(40,32,8)  ,  IntVect(55,47,19) ));
  boxes.push_back(Box(IntVect(40,32,52) ,  IntVect(55,47,63) ));
  boxes.push_back(Box(IntVect(40,32,96) ,  IntVect(55,43,107)));
  boxes.push_back(Box(IntVect(40,44,96) ,  IntVect(55,55,107)));
  boxes.push_back(Box(IntVect(40,48,8)  ,  IntVect(55,63,19) ));
  boxes.push_back(Box(IntVect(40,48,52) ,  IntVect(55,63,63) ));
  boxes.push_back(Box(IntVect(40,56,96) ,  IntVect(55,71,107)));
  boxes.push_back(Box(IntVect(40,64,8)  ,  IntVect(55,79,19) ));
  boxes.push_back(Box(IntVect(40,64,52) ,  IntVect(55,79,63) ));
  boxes.push_back(Box(IntVect(40,72,96) ,  IntVect(55,79,107)));
  boxes.push_back(Box(IntVect(40,80,8)  ,  IntVect(55,95,19) ));
  boxes.push_back(Box(IntVect(40,80,52) ,  IntVect(55,95,63) ));
  boxes.push_back(Box(IntVect(40,80,96) ,  IntVect(55,91,107)));
  boxes.push_back(Box(IntVect(40,92,96) ,  IntVect(55,107,107)));
  boxes.push_back(Box(IntVect(40,96,8)  ,  IntVect(55,111,19)));
  boxes.push_back(Box(IntVect(40,96,52) ,  IntVect(55,111,63)));
  boxes.push_back(Box(IntVect(40,108,96),  IntVect(55,115,107)));
  boxes.push_back(Box(IntVect(40,112,8) ,  IntVect(55,127,19)));
  boxes.push_back(Box(IntVect(40,112,52),  IntVect(55,127,63)));
  boxes.push_back(Box(IntVect(40,116,96),  IntVect(55,127,107)));

  boxes.push_back(Box(IntVect(8,48,72)   , IntVect(19,63,87)  ));
  boxes.push_back(Box(IntVect(8,64,72)   , IntVect(19,79,87)  ));
  boxes.push_back(Box(IntVect(8,80,40)   , IntVect(19,87,55)  ));
  boxes.push_back(Box(IntVect(12,40,60)  , IntVect(19,47,67)  ));
  boxes.push_back(Box(IntVect(12,56,36)  , IntVect(19,67,47)  ));
  boxes.push_back(Box(IntVect(12,88,48)  , IntVect(23,95,63)  ));
  boxes.push_back(Box(IntVect(16,40,88)  , IntVect(19,43,91)  ));
  boxes.push_back(Box(IntVect(16,56,48)  , IntVect(27,67,55)  ));
  boxes.push_back(Box(IntVect(16,88,88)  , IntVect(27,95,95)  ));
  boxes.push_back(Box(IntVect(20,28,56)  , IntVect(31,35,67)  ));
  boxes.push_back(Box(IntVect(20,36,56)  , IntVect(31,47,67)  ));
  boxes.push_back(Box(IntVect(20,40,48)  , IntVect(31,55,55)  ));
  boxes.push_back(Box(IntVect(20,48,72)  , IntVect(27,51,79)  ));
  boxes.push_back(Box(IntVect(20,52,72)  , IntVect(23,59,79)  ));
  boxes.push_back(Box(IntVect(20,56,80)  , IntVect(23,63,87)  ));
  boxes.push_back(Box(IntVect(20,68,80)  , IntVect(23,79,87)  ));
  boxes.push_back(Box(IntVect(20,80,40)  , IntVect(31,87,55)  ));
  boxes.push_back(Box(IntVect(20,88,96)  , IntVect(23,91,99)  ));
  boxes.push_back(Box(IntVect(24,32,88)  , IntVect(39,43,95)  ));
  boxes.push_back(Box(IntVect(24,44,96)  , IntVect(39,55,107) ));
  boxes.push_back(Box(IntVect(24,68,84)  , IntVect(27,71,87)  ));
  boxes.push_back(Box(IntVect(24,88,32)  , IntVect(27,91,35)  ));
  boxes.push_back(Box(IntVect(24,88,96)  , IntVect(27,95,103) ));
  boxes.push_back(Box(IntVect(28,40,36)  , IntVect(31,43,39)  ));
  boxes.push_back(Box(IntVect(28,88,88)  , IntVect(39,95,95)  ));
  boxes.push_back(Box(IntVect(32,28,48)  , IntVect(39,39,63)  ));
  boxes.push_back(Box(IntVect(32,36,32)  , IntVect(39,43,39)  ));
  boxes.push_back(Box(IntVect(32,44,20)  , IntVect(39,59,31)  ));
  boxes.push_back(Box(IntVect(32,56,96)  , IntVect(39,71,103) ));
  boxes.push_back(Box(IntVect(32,72,20)  , IntVect(39,83,31)  ));
  boxes.push_back(Box(IntVect(32,76,80)  , IntVect(39,91,87)  ));
  boxes.push_back(Box(IntVect(32,84,36)  , IntVect(39,95,47)  ));
  boxes.push_back(Box(IntVect(32,96,32)  , IntVect(35,99,39)  ));
  boxes.push_back(Box(IntVect(36,48,108) , IntVect(39,55,115) ));
  boxes.push_back(Box(IntVect(40,0,20)   , IntVect(55,15,31)  ));
  boxes.push_back(Box(IntVect(40,0,64)   , IntVect(55,15,75)  ));
  boxes.push_back(Box(IntVect(40,0,108)  , IntVect(55,15,115) ));
  boxes.push_back(Box(IntVect(40,16,20)  , IntVect(55,31,31)  ));
  boxes.push_back(Box(IntVect(40,16,64)  , IntVect(55,31,75)  ));
  boxes.push_back(Box(IntVect(40,16,108) , IntVect(55,31,115) ));
  boxes.push_back(Box(IntVect(40,32,20)  , IntVect(55,47,31)  ));
  boxes.push_back(Box(IntVect(40,32,64)  , IntVect(55,47,75)  ));
  boxes.push_back(Box(IntVect(40,32,108) , IntVect(55,43,115) ));
  boxes.push_back(Box(IntVect(40,44,108) , IntVect(55,55,115) ));
  boxes.push_back(Box(IntVect(40,48,20)  , IntVect(55,63,31)  ));
  boxes.push_back(Box(IntVect(40,48,64)  , IntVect(55,63,75)  ));
  boxes.push_back(Box(IntVect(40,56,108) , IntVect(55,71,115) ));
  boxes.push_back(Box(IntVect(40,64,20)  , IntVect(55,79,31)  ));
  boxes.push_back(Box(IntVect(40,64,64)  , IntVect(55,79,75)  ));
  boxes.push_back(Box(IntVect(40,72,108) , IntVect(55,79,115) ));
  boxes.push_back(Box(IntVect(40,80,20)  , IntVect(55,95,31)  ));
  boxes.push_back(Box(IntVect(40,80,64)  , IntVect(55,95,75)  ));
  boxes.push_back(Box(IntVect(40,80,108) , IntVect(55,91,115) ));
  boxes.push_back(Box(IntVect(40,92,108) , IntVect(55,107,115)));
  boxes.push_back(Box(IntVect(40,96,20)  , IntVect(55,111,31) ));
  boxes.push_back(Box(IntVect(40,96,64)  , IntVect(55,111,75) ));
  boxes.push_back(Box(IntVect(40,108,108), IntVect(55,115,115)));
  boxes.push_back(Box(IntVect(40,112,20) , IntVect(55,127,31) ));
  boxes.push_back(Box(IntVect(40,112,64) , IntVect(55,127,75) ));
  boxes.push_back(Box(IntVect(40,116,108), IntVect(55,127,115)));

  boxes.push_back(Box(IntVect(8,56,48)   , IntVect(15,67,55)  ));
  boxes.push_back(Box(IntVect(8,68,48)   , IntVect(15,79,55)  ));
  boxes.push_back(Box(IntVect(8,80,56)   , IntVect(19,87,71)  ));
  boxes.push_back(Box(IntVect(12,40,68)  , IntVect(19,47,75)  ));
  boxes.push_back(Box(IntVect(12,56,96)  , IntVect(19,67,103) ));
  boxes.push_back(Box(IntVect(12,88,64)  , IntVect(23,103,79) ));
  boxes.push_back(Box(IntVect(16,44,36)  , IntVect(31,55,47)  ));
  boxes.push_back(Box(IntVect(16,68,48)  , IntVect(27,79,55)  ));
  boxes.push_back(Box(IntVect(16,92,44)  , IntVect(19,95,47)  ));
  boxes.push_back(Box(IntVect(20,28,68)  , IntVect(31,35,75)  ));
  boxes.push_back(Box(IntVect(20,36,68)  , IntVect(31,47,75)  ));
  boxes.push_back(Box(IntVect(20,44,96)  , IntVect(23,47,99)  ));
  boxes.push_back(Box(IntVect(20,48,80)  , IntVect(31,55,87)  ));
  boxes.push_back(Box(IntVect(20,56,24)  , IntVect(31,67,35)  ));
  boxes.push_back(Box(IntVect(20,56,96)  , IntVect(31,67,111) ));
  boxes.push_back(Box(IntVect(20,68,96)  , IntVect(31,79,111) ));
  boxes.push_back(Box(IntVect(20,80,56)  , IntVect(31,87,71)  ));
  boxes.push_back(Box(IntVect(20,92,40)  , IntVect(23,99,47)  ));
  boxes.push_back(Box(IntVect(24,32,96)  , IntVect(39,43,107) ));
  boxes.push_back(Box(IntVect(24,48,28)  , IntVect(31,55,35)  ));
  boxes.push_back(Box(IntVect(24,72,80)  , IntVect(27,79,87)  ));
  boxes.push_back(Box(IntVect(24,88,36)  , IntVect(31,103,47) ));
  boxes.push_back(Box(IntVect(24,96,48)  , IntVect(31,107,63) ));
  boxes.push_back(Box(IntVect(28,44,32)  , IntVect(31,47,35)  ));
  boxes.push_back(Box(IntVect(28,88,96)  , IntVect(39,95,107) ));
  boxes.push_back(Box(IntVect(32,28,80)  , IntVect(39,39,87)  ));
  boxes.push_back(Box(IntVect(32,36,64)  , IntVect(39,47,79)  ));
  boxes.push_back(Box(IntVect(32,44,32)  , IntVect(39,59,47)  ));
  boxes.push_back(Box(IntVect(32,56,104) , IntVect(39,71,115) ));
  boxes.push_back(Box(IntVect(32,72,32)  , IntVect(39,83,47)  ));
  boxes.push_back(Box(IntVect(32,80,48)  , IntVect(39,91,63)  ));
  boxes.push_back(Box(IntVect(32,92,48)  , IntVect(39,107,63) ));
  boxes.push_back(Box(IntVect(32,96,40)  , IntVect(39,107,47) ));
  boxes.push_back(Box(IntVect(36,96,32)  , IntVect(39,103,39) ));
  boxes.push_back(Box(IntVect(40,0,32)   , IntVect(55,15,43)  ));
  boxes.push_back(Box(IntVect(40,0,76)   , IntVect(55,15,87)  ));
  boxes.push_back(Box(IntVect(40,0,116)  , IntVect(55,15,127) ));
  boxes.push_back(Box(IntVect(40,16,32)  , IntVect(55,31,43)  ));
  boxes.push_back(Box(IntVect(40,16,76)  , IntVect(55,31,87)  ));
  boxes.push_back(Box(IntVect(40,16,116) , IntVect(55,31,127) ));
  boxes.push_back(Box(IntVect(40,32,32)  , IntVect(55,47,43)  ));
  boxes.push_back(Box(IntVect(40,32,76)  , IntVect(55,47,87)  ));
  boxes.push_back(Box(IntVect(40,32,116) , IntVect(55,43,127) ));
  boxes.push_back(Box(IntVect(40,44,116) , IntVect(55,55,127) ));
  boxes.push_back(Box(IntVect(40,48,32)  , IntVect(55,63,43)  ));
  boxes.push_back(Box(IntVect(40,48,76)  , IntVect(55,63,87)  ));
  boxes.push_back(Box(IntVect(40,56,116) , IntVect(55,71,127) ));
  boxes.push_back(Box(IntVect(40,64,32)  , IntVect(55,79,43)  ));
  boxes.push_back(Box(IntVect(40,64,76)  , IntVect(55,79,87)  ));
  boxes.push_back(Box(IntVect(40,72,116) , IntVect(55,79,127) ));
  boxes.push_back(Box(IntVect(40,80,32)  , IntVect(55,95,43)  ));
  boxes.push_back(Box(IntVect(40,80,76)  , IntVect(55,95,87)  ));
  boxes.push_back(Box(IntVect(40,80,116) , IntVect(55,91,127) ));
  boxes.push_back(Box(IntVect(40,92,116) , IntVect(55,107,127)));
  boxes.push_back(Box(IntVect(40,96,32)  , IntVect(55,111,43) ));
  boxes.push_back(Box(IntVect(40,96,76)  , IntVect(55,111,87) ));
  boxes.push_back(Box(IntVect(40,108,116), IntVect(55,115,127)));
  boxes.push_back(Box(IntVect(40,112,32) , IntVect(55,127,43) ));
  boxes.push_back(Box(IntVect(40,112,76) , IntVect(55,127,87) ));
  boxes.push_back(Box(IntVect(40,116,116), IntVect(55,127,127)));

  boxes.push_back(Box(IntVect(8,56,88)   , IntVect(23,71,95)  ));
  boxes.push_back(Box(IntVect(8,72,88)   , IntVect(23,87,95)  ));
  boxes.push_back(Box(IntVect(8,80,72)   , IntVect(19,87,87)  ));
  boxes.push_back(Box(IntVect(12,40,76)  , IntVect(19,47,87)  ));
  boxes.push_back(Box(IntVect(12,68,36)  , IntVect(19,79,47)  ));
  boxes.push_back(Box(IntVect(12,88,80)  , IntVect(19,95,87)  ));
  boxes.push_back(Box(IntVect(16,48,96)  , IntVect(23,55,99)  ));
  boxes.push_back(Box(IntVect(16,80,96)  , IntVect(23,87,103) ));
  boxes.push_back(Box(IntVect(16,96,52)  , IntVect(19,99,55)  ));
  boxes.push_back(Box(IntVect(20,28,76)  , IntVect(31,35,87)  ));
  boxes.push_back(Box(IntVect(20,36,76)  , IntVect(31,47,87)  ));
  boxes.push_back(Box(IntVect(20,48,56)  , IntVect(27,55,63)  ));
  boxes.push_back(Box(IntVect(20,48,100) , IntVect(23,55,107) ));
  boxes.push_back(Box(IntVect(20,56,36)  , IntVect(31,67,47)  ));
  boxes.push_back(Box(IntVect(20,68,24)  , IntVect(31,79,35)  ));
  boxes.push_back(Box(IntVect(20,76,72)  , IntVect(27,79,79)  ));
  boxes.push_back(Box(IntVect(20,80,72)  , IntVect(31,87,87)  ));
  boxes.push_back(Box(IntVect(20,96,48)  , IntVect(23,103,55) ));
  boxes.push_back(Box(IntVect(24,36,40)  , IntVect(31,43,47)  ));
  boxes.push_back(Box(IntVect(24,48,108) , IntVect(31,55,111) ));
  boxes.push_back(Box(IntVect(24,72,88)  , IntVect(39,87,95)  ));
  boxes.push_back(Box(IntVect(24,88,48)  , IntVect(31,95,63)  ));
  boxes.push_back(Box(IntVect(24,96,64)  , IntVect(31,107,79) ));
  boxes.push_back(Box(IntVect(28,76,84)  , IntVect(31,79,87)  ));
  boxes.push_back(Box(IntVect(28,96,88)  , IntVect(31,103,95) ));
  boxes.push_back(Box(IntVect(32,28,88)  , IntVect(39,31,95)  ));
  boxes.push_back(Box(IntVect(32,40,48)  , IntVect(39,51,63)  ));
  boxes.push_back(Box(IntVect(32,44,108) , IntVect(39,47,111) ));
  boxes.push_back(Box(IntVect(32,60,20)  , IntVect(39,71,31)  ));
  boxes.push_back(Box(IntVect(32,72,96)  , IntVect(39,87,103) ));
  boxes.push_back(Box(IntVect(32,80,64)  , IntVect(39,91,79)  ));
  boxes.push_back(Box(IntVect(32,92,64)  , IntVect(39,107,79) ));
  boxes.push_back(Box(IntVect(32,96,88)  , IntVect(39,103,91) ));
  boxes.push_back(Box(IntVect(40,0,0)    , IntVect(55,15,7)   ));
  boxes.push_back(Box(IntVect(40,0,44)   , IntVect(55,15,51)  ));
  boxes.push_back(Box(IntVect(40,0,88)   , IntVect(55,15,95)  ));
  boxes.push_back(Box(IntVect(40,16,0)   , IntVect(55,31,7)   ));
  boxes.push_back(Box(IntVect(40,16,44)  , IntVect(55,31,51)  ));
  boxes.push_back(Box(IntVect(40,16,88)  , IntVect(55,31,95)  ));
  boxes.push_back(Box(IntVect(40,32,0)   , IntVect(55,47,7)   ));
  boxes.push_back(Box(IntVect(40,32,44)  , IntVect(55,47,51)  ));
  boxes.push_back(Box(IntVect(40,32,88)  , IntVect(55,43,95)  ));
  boxes.push_back(Box(IntVect(40,44,88)  , IntVect(55,55,95)  ));
  boxes.push_back(Box(IntVect(40,48,0)   , IntVect(55,63,7)   ));
  boxes.push_back(Box(IntVect(40,48,44)  , IntVect(55,63,51)  ));
  boxes.push_back(Box(IntVect(40,56,88)  , IntVect(55,71,95)  ));
  boxes.push_back(Box(IntVect(40,64,0)   , IntVect(55,79,7)   ));
  boxes.push_back(Box(IntVect(40,64,44)  , IntVect(55,79,51)  ));
  boxes.push_back(Box(IntVect(40,72,88)  , IntVect(55,79,95)  ));
  boxes.push_back(Box(IntVect(40,80,0)   , IntVect(55,95,7)   ));
  boxes.push_back(Box(IntVect(40,80,44)  , IntVect(55,95,51)  ));
  boxes.push_back(Box(IntVect(40,80,88)  , IntVect(55,91,95)  ));
  boxes.push_back(Box(IntVect(40,92,88)  , IntVect(55,107,95) ));
  boxes.push_back(Box(IntVect(40,96,0)   , IntVect(55,111,7)  ));
  boxes.push_back(Box(IntVect(40,96,44)  , IntVect(55,111,51) ));
  boxes.push_back(Box(IntVect(40,108,88) , IntVect(55,115,95) ));
  boxes.push_back(Box(IntVect(40,112,0)  , IntVect(55,127,7)  ));
  boxes.push_back(Box(IntVect(40,112,44) , IntVect(55,127,51) ));
  boxes.push_back(Box(IntVect(40,116,88) , IntVect(55,127,95) ));
#else

  boxes.push_back(Box(IntVect(40,112) , IntVect(55,127) ));
  boxes.push_back(Box(IntVect(20,116) , IntVect(39,127) ));

#endif

  Vector<int> procs(boxes.size(), 0);
  a_dbl.define(boxes, procs);
}
/************/
int
obscura(const DisjointBoxLayout&    a_dbl,
        const Box&                  a_domain,
        const Real&                 a_dx)
{

  LevelData<BaseIFFAB<Real> > fluxInterpolants[SpaceDim];

  EBISLayout ebisl;
  const CH_XD::EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int nghost = 5;
  ebisPtr->fillEBISLayout(ebisl, a_dbl, a_domain, nghost);

  LayoutData<IntVectSet>      irregSetsSmall;
  LayoutData<IntVectSet>      irregSetsGrown[SpaceDim];

  irregSetsSmall.define(a_dbl);
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      irregSetsGrown[faceDir].define(a_dbl);
    }

  int ibox = 0;
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      ibox = 0;
      for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = ebisl[dit()];

          const Box&     thisBox = a_dbl.get(dit());
          IntVectSet ivsIrreg = ebisBox.getIrregIVS(thisBox);
          irregSetsSmall[dit()] = ivsIrreg;

          IntVectSet& grownIVS = irregSetsGrown[faceDir][dit()];
          //need to do this in case there is an irregular cells
          //just outside the box.
          Box grownBox2 = grow(thisBox, 2);
          grownBox2 &= a_domain;
          IntVectSet irregivs = ebisBox.getIrregIVS(grownBox2);
          grownIVS = irregivs;
          //          Box grownBox1 = thisBox;
          //grow in all directions != faceDir
          grownIVS.grow(1);
          for (int jdir = 0; jdir < SpaceDim; jdir++)
            {
              if (jdir != faceDir)
                {
                  //grownIVS.grow(jdir, 1);
                  //                  grownBox1.grow(jdir, 1);
                }
            }
          //restrict to domain
          //grownBox1 &= a_domain;
          //grownIVS &= grownBox1;
          grownIVS &= a_domain;
          ibox++;
        }
    }

  //define flux interpolant stuff
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      BaseIFFactory<Real> faceFactory(ebisl, irregSetsGrown[faceDir], faceDir);
      fluxInterpolants[faceDir].define(a_dbl,  1,  3*IntVect::Unit, faceFactory);
    }

  //now set the values of stuff ONLY on the grids
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      ibox = 0;
      for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
        {
          const EBGraph& ebgraph = ebisl[dit()].getEBGraph();
          FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

          BaseIFFAB<Real>& interpol = fluxInterpolants[faceDir][dit()];
          interpol.setVal(-777.);
          const Box& grid = a_dbl.get(dit());
          IntVectSet ivsGrid(grid);
          ivsGrid &= interpol.getIVS();
          for (FaceIterator faceit(ivsGrid, ebgraph, faceDir, stopCrit);  faceit.ok(); ++faceit)
            {
              interpol(faceit(), 0) = ibox;
            }
          ibox ++;
        }
    }
  Interval interv(0,0);
  //exchange
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      fluxInterpolants[faceDir].exchange(interv);
    }
  IntVectSet ivsGrids;
  for (LayoutIterator lit = a_dbl.layoutIterator(); lit.ok(); ++lit)
    {
      ivsGrids |= a_dbl.get(lit());
    }

  //check everywhere that is on the grids (including ghost) has reasonable values
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      ibox = 0;
      for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
        {

          const EBGraph& ebgraph = ebisl[dit()].getEBGraph();
          FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

          BaseIFFAB<Real>& interpol = fluxInterpolants[faceDir][dit()];
          const IntVectSet& interpolIVS = interpol.getIVS();

          IntVectSet ivs = ebgraph.getIrregCells(a_dbl.get(dit()));
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir != faceDir)
                ivs.grow(idir, 1);
            }
          ivs &= a_domain;
          for (FaceIterator faceit(ivs, ebgraph, faceDir, stopCrit);  faceit.ok(); ++faceit)
            {
              //test only faces that contain cells in real grids
              bool testThisFace = false;
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  IntVect iv = faceit().gridIndex(sit());
                  if (ivsGrids.contains(iv) && interpolIVS.contains(iv))
                    {
                      testThisFace = true;
                    }

                }
              if (testThisFace && interpol(faceit(), 0) < 0)
                {
                  pout() << "failed at ibox = " << ibox << ", face = " << faceit().gridIndex(Side::Lo) << faceit().gridIndex(Side::Hi) << endl;
                  return -1;
                }
            }
          ibox++;
        }
    }

  return 0;
}
/************/
int
iffabExchangeTest()
{
  Box domain;
  Real dx;
  getFinestDomain(domain, dx);

  DisjointBoxLayout dbl;
  makeDBL(dbl);

  makeGeometry(domain,  dx);

  int retval = obscura(dbl, domain, dx);

  return retval;
}
/************/
int
main(int a_argc, char* a_argv[])
{
  int retval = 0;
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif
    // Check for an input file

    ParmParse pp(0, NULL, NULL, "cylinder.inputs");

    retval = iffabExchangeTest();

    if (retval == 0)
      {
        pout() << "iffabExchangeTest passed " << endl;
      }
    else
      {
        pout() << "iffabExchangeTest failed with code " << retval << endl;
      }
#ifdef CH_MPI
  }
  MPI_Finalize();
#endif

  return retval;
}
