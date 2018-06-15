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
#include <list>
#include <set>
using std::cout;

#include "parstream.H"

#include "DataIterator.H"
#include "Misc.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "TestCommon.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"
/////================================///////
void
setCircleTags(IntVectSet& ivs, int circleR, int thickness,
              const IntVect& center )
{
  ivs.define();
  CH_assert(circleR > 0);
  CH_assert(thickness < circleR && thickness > 0);
  for (int t = 0; t<thickness; ++t)
  {
    for (int x=-circleR-t; x<=circleR+t; ++x)
    {
#if CH_SPACEDIM == 1
      ivs |= center + IntVect(x);
#elif CH_SPACEDIM == 2
        // circle
        // variables not used, avoid warning (ndk)
        //IntVect loVect = center+IntVect(x,Y);
        //IntVect hiVect = center+IntVect(x,-Y);
      int Y = (int)sqrt((Real)abs((circleR+t)*(circleR+t) - x*x));
      // not perfect, but I'm not fussy here
        ivs |= center+IntVect(x,Y);
        ivs |= center+IntVect(x,-Y);
#elif CH_SPACEDIM == 3
        //sphere
        int Y = (int)sqrt((Real)abs((circleR+t)*(circleR+t) - x*x));
        // not perfect, but I'm not fussy here
        for (int y=-Y; y<=Y; ++y)
          {
            // must cast to Real inside sqrt()  (ndk)
            int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
            ivs |= center+IntVect(x,y,Z);
            ivs |= center+IntVect(x,y,-Z);
          }

#elif CH_SPACEDIM == 4
        // hypersphere
        int Y = (int)sqrt((Real)abs((circleR+t)*(circleR+t) - x*x));
        // not perfect, but I'm not fussy here
        for (int y=-Y; y<=Y; ++y)
          {
            // must cast to Real inside sqrt()  (ndk)
            int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
            for (int z=-Z; z<=Z; ++z)
              {
                int U = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y -z*z));
                ivs |= center+IntVect(x,y,z,U);
                ivs |= center+IntVect(x,y,z,-U);
              }
          }

#elif CH_SPACEDIM == 5
        // hypersphere
        int Y = (int)sqrt((Real)abs((circleR+t)*(circleR+t) - x*x));
        // not perfect, but I'm not fussy here
        for (int y=-Y; y<=Y; ++y)
          {
            // must cast to Real inside sqrt()  (ndk)
            int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
            for (int z=-Z; z<=Z; ++z)
              {
                int U = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y -z*z));
                for (int u=-U; u<=U; ++u)
                  {
                    int V = (int)sqrt((Real)abs(circleR*circleR-x*x-y*y-z*z-u*u));
                    ivs |= center+IntVect(x,y,z,u,V);
                    ivs |= center+IntVect(x,y,z,u,-V);
                  }
              }
          }
#elif CH_SPACEDIM == 6
        // hypersphere
        int Y = (int)sqrt((Real)abs((circleR+t)*(circleR+t) - x*x));
        // not perfect, but I'm not fussy here
        for (int y=-Y; y<=Y; ++y)
          {
            // must cast to Real inside sqrt()  (ndk)
            int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
            for (int z=-Z; z<=Z; ++z)
              {
                int U = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y -z*z));
                for (int u=-U; u<=U; ++u)
                  {
                    int V = (int)sqrt((Real)abs(circleR*circleR-x*x-y*y-z*z-u*u));
                    for (int v=-V; v<=V; ++v)
                      {
                        int W = (int)sqrt((Real)abs(circleR*circleR-x*x-y*y-z*z-u*u-v*v));
                        ivs |= center+IntVect(x,y,z,u,v,W);
                        ivs |= center+IntVect(x,y,z,u,v,-W);
                      }
                  }
              }
          }
#else
#error "setCircleTags() implemented only for 2D, 3D,4D,5D, and 6D"
#endif
    }
  }
}
/////================================///////
void
buildDisjointBoxLayout(DisjointBoxLayout& plan,const IntVectSet& tags, const Box& domain)
{
  Vector<Box> vectBox;
  Vector<int> assignments;

  IntVectSet domainivs(domain);
  int maxsize = 128;
  Vector<int> fakeNRef(2,2);
  Real fillRatio = 0.5;
  int bufferSize = 0;
  int blockFactor = 1;
  BRMeshRefine mr(domain, fakeNRef, fillRatio, blockFactor,
                  bufferSize, maxsize);

  mr.makeBoxes(vectBox, tags, domainivs, domain, maxsize, 0);

  vectBox.sort();
  int stat = LoadBalance(assignments,  vectBox);
  if ( stat != 0 )
  {
    MayDay::Error("loadBalance() FAILED.",stat);
  }

  plan.define(vectBox, assignments);

  plan.close();
}

#include "NamespaceFooter.H"
