#include "pluto.h"

/* ************************************************************** */
void EMF_USERDEF_BOUNDARY (EMF *emf, int side, int loc, Grid *grid)
/*
 *
 *
 * PURPOSE:
 *
 *  Set boundary conditions on electric field. 
 *  Quantities are evaluated at zone face, NOT edge.
 *  
 *  The electric field components are given by 
 *  the fluxes of the induction equation:
 *
 *   +----+     
 *   |  0 |
 *   |-Ezi|  at x - faces
 *   | Eyi|
 *   +----+
 *
 *   +----+     
 *   | Ezj|
 *   | 0  |  at y - faces
 *   |-Exj|
 *   +----+
 *
 *   +----+     
 *   |-Eyk|
 *   | Exk|  at z - faces
 *   | 0  |
 *   +----+
 *
 *  Only the components required to update the normal megnetic
 *  field are required.
 *  For example, at the lower x-boundary, the necessary 
 *  electric field components to update bx are Ez and Ey (no Ex).
 *  Thus one needs to specify 
 *
 *     Ez at y-faces   (ezj)
 *     Ey at z-faces   (eyk)
 *  
 *
 * LAST MODIFIED: 
 *
 *   February 21 2008 by A. Mignone (email:mignone@to.astro.it)
 *
 *
 *
 *
 **************************************************************** */
{
  int i, j, k;

  return; /* !! comment this line if you intend to use this function !! */

  if (loc != FACE_EMF) return;

  if (side == X1_BEG) {

    X1_BEG_LOOP(k,j,i){
      D_EXPAND(                                    ,
        emf->ezj[k][j][i] = emf->ezj[k][j][IBEG];  ,
        emf->eyk[k][j][i] = emf->eyk[k][j][IBEG];
      )
    }

  }else if (side == X1_END){

    X1_END_LOOP(k,j,i){
      D_EXPAND(                                    ,
        emf->ezj[k][j][i] = emf->ezj[k][j][IEND];  ,
        emf->eyk[k][j][i] = emf->eyk[k][j][IEND];  
      )
    }

  }else if (side == X2_BEG){

    X2_BEG_LOOP(k,j,i){
      D_EXPAND(                                    ,
        emf->ezi[k][j][i] = emf->ezi[k][JBEG][i];  ,
        emf->exk[k][j][i] = emf->exk[k][JBEG][i];
      )
    }

  }else if (side == X2_END){

    X2_END_LOOP(k,j,i){
      D_EXPAND(                                    ,
        emf->ezi[k][j][i] = emf->ezi[k][JEND][i];  ,
        emf->exk[k][j][i] = emf->exk[k][JEND][i];
      )
    }

  }else if (side == X3_BEG){
 
    X3_BEG_LOOP(k,j,i){
      emf->exj[k][j][i] = emf->exj[KBEG][j][i]; 
      emf->eyi[k][j][i] = emf->eyi[KBEG][j][i];
    } 

  }else if (side == X3_END){

    X3_END_LOOP(k,j,i){
      emf->exj[k][j][i] = emf->exj[KEND][j][i];
      emf->eyi[k][j][i] = emf->eyi[KEND][j][i];
    }

  }
}


