/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Collects various functions to operate on the RBox structure.

  - RBoxDefine() is used to set the box extent in terms of six
    indices (beg..end for each direction) and the variable location
    inside the cell.

  - RBoxSetDirections() is used to set normal, tangent and binormal indices
    with respect to the specified (sweeping) direction.
    The convention adopted here is the same one used for vector indices
    permutation in PLUTO: <tt> (i,j,k) -> (j,i,k) -> (k,i,j) </tt>.
    Useful for sweeping along different directions during the time stepping
    routines.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Nov 13, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void RBoxDefine(int ib, int ie, int jb, int je, int kb, int ke, int vpos, RBox *box)
/*! 
 * 
 * \param [in]  ib    leftmost  index in the X1 direction
 * \param [in]  ie    rightmost index in the X1 direction
 * \param [in]  jb    leftmost  index in the X2 direction
 * \param [in]  je    rightmost index in the X2 direction
 * \param [in]  kb    leftmost  index in the X3 direction
 * \param [in]  ke    rightmost index in the X3 direction
 * \param [in] vpos   the variable location inside the cell
 *                    (CENTER/X1FACE/.../X3FACE)
 * \param [out] box   pointer to a RBox structure
 *
 *********************************************************************** */
{
  box->ibeg = ib;
  box->iend = ie;

  box->jbeg = jb;
  box->jend = je;

  box->kbeg = kb;
  box->kend = ke;
 
  box->vpos = vpos;   
}

/* ********************************************************************* */
void RBoxSetDirections(RBox *box, int dir)
/*!
 * Set normal, tangent and binormal directions while sweeping
 * across a box using the BOX_TRANSVERSE_LOOP macro;
 *
 * \param [in,out]  box  pointer to a RBox structure
 * \param [in]      dir  the sweeping direction giving the normal direction
 *                       and respect to which assign the tangent and binormal
 *                       indices.
 *********************************************************************** */
{
  if (dir == IDIR){
    box->nbeg = &(box->ibeg); box->nend = &(box->iend);
    box->tbeg = &(box->jbeg); box->tend = &(box->jend);
    box->bbeg = &(box->kbeg); box->bend = &(box->kend);
  }else if (dir == JDIR){
    box->nbeg = &(box->jbeg); box->nend = &(box->jend);
    box->tbeg = &(box->ibeg); box->tend = &(box->iend);
    box->bbeg = &(box->kbeg); box->bend = &(box->kend);
  }else if (dir == KDIR){
    box->nbeg = &(box->kbeg); box->nend = &(box->kend);
    box->tbeg = &(box->ibeg); box->tend = &(box->iend);
    box->bbeg = &(box->jbeg); box->bend = &(box->jend);
  }else{
    print ("! RBoxSetDirections(): invalid dir = %d\n",dir);
    QUIT_PLUTO(1);
  }

}

/* ********************************************************************* */
void RBoxShow(RBox *box)
/*
 *
 *
 *********************************************************************** */
{
  print ("===============================================================\n");
  print (" (ibeg, iend) = (%d, %d)\n",box->ibeg, box->iend);
  print (" (jbeg, jend) = (%d, %d)\n",box->jbeg, box->jend);
  print (" (kbeg, kend) = (%d, %d)\n",box->kbeg, box->kend);
  print ("===============================================================\n");
 
}
