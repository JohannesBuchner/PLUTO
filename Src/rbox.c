/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Basic RBox database.

  The function SetRBox() defines an array of predefined RBox structures 
  needed to loop through different portions (or regions) of the 
  computational domain.
  
  - rbox_center: an array of RBoxes for looping over the cell-centered data.
  - rbox_x1face: an array of RBoxes for looping over the X1-staggered data.
  - rbox_x2face: an array of RBoxes for looping over the X2-staggered data.
  - rbox_x3face: an array of RBoxes for looping over the X3-staggered data.
 
  Each array of structures has 8 elements corresponding to the six sides
  of the computational domain (X1_BEG, ... , X3_END) and, in addition, 
  we also define the DOM and TOT array indices to loop over the active 
  computational zones or over the total
  computational domain (interior + ghost zones), respectively.

  The function GetRBox() can be used to retrieve a pointer to a RBox
  structure given the computational side and the variable position.

  \author A. Mignone (mignone@ph.unito.it)
  \date   May 13, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static RBox rbox_center[8], rbox_x1face[8], rbox_x2face[8], rbox_x3face[8];

/* ********************************************************************* */
void SetRBox(void)
/*! 
 *
 *
 *********************************************************************** */
{
  int s;

/* ---------------------------------------------------
    0. Set X1_BEG grid index ranges
   --------------------------------------------------- */

  s = X1_BEG; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;

  rbox_center[s].ib = IBEG-1; rbox_center[s].ie =         0;
  rbox_center[s].jb =      0; rbox_center[s].je = NX2_TOT-1;
  rbox_center[s].kb =      0; rbox_center[s].ke = NX3_TOT-1;

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(rbox_x1face[s].ib--; rbox_x1face[s].ie--;  ,
            rbox_x2face[s].jb--;                       ,
            rbox_x3face[s].kb--;)
  #endif

/* ---------------------------------------------------
    1. set X1_END grid index ranges
   --------------------------------------------------- */
  
  s = X1_END; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;

  rbox_center[s].ib = IEND+1; rbox_center[s].ie = NX1_TOT-1;
  rbox_center[s].jb =      0; rbox_center[s].je = NX2_TOT-1;
  rbox_center[s].kb =      0; rbox_center[s].ke = NX3_TOT-1;

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(                   ;   ,
            rbox_x2face[s].jb--;   ,
            rbox_x3face[s].kb--;)
  #endif

/* ---------------------------------------------------
    2. set X2_BEG grid index ranges
   --------------------------------------------------- */

  s = X2_BEG; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;

  rbox_center[s].ib =      0; rbox_center[s].ie = NX1_TOT-1;
  rbox_center[s].jb = JBEG-1; rbox_center[s].je =         0;
  rbox_center[s].kb =      0; rbox_center[s].ke = NX3_TOT-1;

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(rbox_x1face[s].ib--;                       ,
            rbox_x2face[s].jb--; rbox_x2face[s].je--;  ,
            rbox_x3face[s].kb--;)
  #endif

/* ---------------------------------------------------
    3. set X2_END grid index ranges
   --------------------------------------------------- */
  
  s = X2_END; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;

  rbox_center[s].ib =      0; rbox_center[s].ie = NX1_TOT-1;
  rbox_center[s].jb = JEND+1; rbox_center[s].je = NX2_TOT-1;
  rbox_center[s].kb =      0; rbox_center[s].ke = NX3_TOT-1;

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(rbox_x1face[s].ib--;    ,
                               ;    ,
            rbox_x3face[s].kb--;)
  #endif

/* ---------------------------------------------------
    4. set X3_BEG grid index ranges
   --------------------------------------------------- */

  s = X3_BEG; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;

  rbox_center[s].ib =      0; rbox_center[s].ie = NX1_TOT-1;
  rbox_center[s].jb =      0; rbox_center[s].je = NX2_TOT-1;
  rbox_center[s].kb = KBEG-1; rbox_center[s].ke =         0;

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(rbox_x1face[s].ib--;   ,
            rbox_x2face[s].jb--;   ,
            rbox_x3face[s].kb--; rbox_x3face[s].ke--;)
  #endif

/* ---------------------------------------------------
    5.  set X3_END grid index ranges
   --------------------------------------------------- */
  
  s = X3_END; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;

  rbox_center[s].ib =      0; rbox_center[s].ie = NX1_TOT-1; 
  rbox_center[s].jb =      0; rbox_center[s].je = NX2_TOT-1; 
  rbox_center[s].kb = KEND+1; rbox_center[s].ke = NX3_TOT-1; 

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(rbox_x1face[s].ib--;      ,
            rbox_x2face[s].jb--;      ,
                               ;)
  #endif

/* ---------------------------------------------------
    6. set DOM index ranges
   --------------------------------------------------- */
  
  s = DOM; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;
  rbox_center[s].ib = IBEG; rbox_center[s].ie = IEND; 
  rbox_center[s].jb = JBEG; rbox_center[s].je = JEND; 
  rbox_center[s].kb = KBEG; rbox_center[s].ke = KEND; 

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(rbox_x1face[s].ib--;   ,
            rbox_x2face[s].jb--;   ,
            rbox_x3face[s].kb--;)
  #endif

/* ---------------------------------------------------
    7. set TOT index ranges
   --------------------------------------------------- */
  
  s = TOT; s -= X1_BEG;

  rbox_center[s].vpos = CENTER;
  rbox_center[s].ib = 0; rbox_center[s].ie = NX1_TOT-1; 
  rbox_center[s].jb = 0; rbox_center[s].je = NX2_TOT-1; 
  rbox_center[s].kb = 0; rbox_center[s].ke = NX3_TOT-1; 

  rbox_x1face[s] = rbox_x2face[s] = rbox_x3face[s] = rbox_center[s];

  #ifndef CHOMBO /* -- useless for AMR -- */
   rbox_x1face[s].vpos = X1FACE;
   rbox_x2face[s].vpos = X2FACE;
   rbox_x3face[s].vpos = X3FACE;

   D_EXPAND(rbox_x1face[s].ib--;   ,
            rbox_x2face[s].jb--;   ,
            rbox_x3face[s].kb--;)
  #endif
}
/* ********************************************************************* */
RBox *GetRBox(int side, int vpos)
/*!
 *  Returns a pointer to a local static RBox 
 *
 *  \param[in]  side  the region of the computational domain where 
 *                    the box is required. There 8 possible values:
 *                    X1_BEG, ... , X3_END, DOM, TOT.
 *  \param[in]  vpos  the variable position inside the cell:
 *                    CENTER, X1FACE, X2FACE or X3FACE.
 *
 *********************************************************************** */
{
  if      (vpos == CENTER) return &(rbox_center[side-X1_BEG]);
  else if (vpos == X1FACE) return &(rbox_x1face[side-X1_BEG]);
  else if (vpos == X2FACE) return &(rbox_x2face[side-X1_BEG]);
  else if (vpos == X3FACE) return &(rbox_x3face[side-X1_BEG]);
}

