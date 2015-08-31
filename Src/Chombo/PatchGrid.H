#ifndef PATCH_GRID_H
#define PATCH_GRID_H

extern "C" {
#include "pluto.h"
}

#include "MayDay.H"
#include "PiecewiseLinearFillPatch.H"
#include "FineInterp.H"

class PatchGrid
{
public:

 PatchGrid(){;}

  ~PatchGrid(){

    if (PGrid[0].x != NULL) { 
      FreeGrid(PGrid);  
    } 
  }

 struct GRID* getGrid()
 {
   return PGrid;
 }

protected:

 struct GRID PGrid[3];

private:


};

class PiecewiseLinearFillPluto : public PiecewiseLinearFillPatch
{
public:

  PiecewiseLinearFillPluto(){m_is_defined = false;}

 ~PiecewiseLinearFillPluto(){;}

  int get_geometry ()
  {
   #if GEOMETRY == CARTESIAN
     return 1;
   #endif

   #if GEOMETRY == CYLINDRICAL
    #if DIMENSIONS == 2
     return 2;
    #else
     MayDay::Error("Cylindrical geometry supported in Pluto/Chombo only in 2 dimensions");
    #endif
   #endif

   #if GEOMETRY == SPHERICAL
    #if CHOMBO_LOGR == NO
     return 3;
    #else
     return 4;
    #endif
   #endif

   #if GEOMETRY == POLAR
    #if CHOMBO_LOGR == NO
     return 5;
    #else
     return 6;
    #endif
   #endif

  }

};

class FineInterpPluto : public FineInterp
{
public:

  FineInterpPluto(){is_defined = false;}

 ~FineInterpPluto(){;}

  int get_geometry ()
  {
   #if GEOMETRY == CARTESIAN
     return 1;
   #endif

   #if GEOMETRY == CYLINDRICAL
    #if DIMENSIONS == 2
     return 2;
    #else
     MayDay::Error("Cylindrical geometry supported in Pluto/Chombo only in 2 dimensions");
    #endif
   #endif

   #if GEOMETRY == SPHERICAL
    #if CHOMBO_LOGR == NO
     return 3;
    #else
     return 4;
    #endif
   #endif

   #if GEOMETRY == POLAR
    #if CHOMBO_LOGR == NO
     return 5;
    #else
     return 6;
    #endif
   #endif

  }

};

#endif
