/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief 3D wrapper for conservative/primitive conversion.

  Provide 3D wrappers to the standard 1D conversion functions
  ConsToPrim() and PrimToCons().

  \authors A. Mignone (mignone@ph.unito.it)
  \date    June 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ConsToPrim3D (Data_Arr U, Data_Arr V, unsigned char ***flag, RBox *box)
/*!
 *  Convert a 3D array of conservative variables \c U to
 *  an array of primitive variables \c V.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]     U      pointer to 3D array of conserved variables,
 *                        with array indexing <tt>[k][j][i][nv]</tt>
 * \param [out]    V      pointer to 3D array of primitive variables,
 *                        with array indexing <tt>[nv][k][j][i]</tt>
 * \param [in,out] flag   pointer to 3D array of flags.
 * \param [in]     box    pointer to RBox structure containing the domain
 *                        portion over which conversion must be performed.
 *
 *********************************************************************** */
{
  int   i, j, k, nv, err;
  int   ibeg, iend, jbeg, jend, kbeg, kend;
  int   current_dir;
  static double **v, **u;

  if (v == NULL){
    v = ARRAY_2D(NMAX_POINT, NVAR, double);
    u = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* ----------------------------------------------
    Save current sweep direction and by default,
    perform the conversion along X1 stripes
   ---------------------------------------------- */

  current_dir = g_dir; 
  g_dir = IDIR;
  
/* -----------------------------------------------
    Set (beg,end) indices in ascending order for
    proper call to ConsToPrim()
   ----------------------------------------------- */

  ibeg = (box->ib <= box->ie) ? (iend=box->ie, box->ib):(iend=box->ib, box->ie);
  jbeg = (box->jb <= box->je) ? (jend=box->je, box->jb):(jend=box->jb, box->je);
  kbeg = (box->kb <= box->ke) ? (kend=box->ke, box->kb):(kend=box->kb, box->ke);

  for (k = kbeg; k <= kend; k++){ g_k = k;
  for (j = jbeg; j <= jend; j++){ g_j = j;

#ifdef CHOMBO
    for (i = ibeg; i <= iend; i++) NVAR_LOOP(nv) u[i][nv] = U[nv][k][j][i];
    #if COOLING == MINEq || COOLING == H2_COOL
    if (g_intStage == 1) for (i = ibeg; i <= iend; i++) NormalizeIons(u[i]);
    #endif
    err = ConsToPrim (u, v, ibeg, iend, flag[k][j]);

  /* -------------------------------------------------------------
      Ensure any change to 1D conservative arrays is not lost.
      Note: Conversion must be done even when err = 0 in case
            ENTROPY_SWITCH is employed.
     ------------------------------------------------------------- */

    for (i = ibeg; i <= iend; i++) NVAR_LOOP(nv) U[nv][k][j][i] = u[i][nv];
#else
    err = ConsToPrim (U[k][j], v, ibeg, iend, flag[k][j]);
#endif
    for (i = ibeg; i <= iend; i++) NVAR_LOOP(nv) V[nv][k][j][i] = v[i][nv];
  }}
  g_dir = current_dir;

}
/* ********************************************************************* */
void PrimToCons3D (Data_Arr V, Data_Arr U, RBox *box)
/*!
 *  Convert a 3D array of primitive variables \c V  to
 *  an array of conservative variables \c U.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]    V     pointer to 3D array of primitive variables,
 *                      with array indexing <tt>[nv][k][j][i]</tt>
 * \param [out]   U     pointer to 3D array of conserved variables,
 *                      with array indexing <tt>[k][j][i][nv]</tt>
 * \param [in]    box   pointer to RBox structure containing the domain
 *                      portion over which conversion must be performed.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  int   ibeg, iend, jbeg, jend, kbeg, kend;
  int   current_dir;
  static double **v, **u;

  if (v == NULL) {
    v = ARRAY_2D(NMAX_POINT, NVAR, double);
    u = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  current_dir = g_dir; /* save current direction */
  g_dir = IDIR;

/* -----------------------------------------------
    Set (beg,end) indices in ascending order for
    proper call to ConsToPrim()
   ----------------------------------------------- */

  ibeg = (box->ib <= box->ie) ? (iend=box->ie, box->ib):(iend=box->ib, box->ie);
  jbeg = (box->jb <= box->je) ? (jend=box->je, box->jb):(jend=box->jb, box->je);
  kbeg = (box->kb <= box->ke) ? (kend=box->ke, box->kb):(kend=box->kb, box->ke);

  for (k = kbeg; k <= kend; k++){ g_k = k;
  for (j = jbeg; j <= jend; j++){ g_j = j;
    for (i = ibeg; i <= iend; i++) VAR_LOOP(nv) v[i][nv] = V[nv][k][j][i];
#ifdef CHOMBO
    PrimToCons(v, u, ibeg, iend);
    for (i = ibeg; i <= iend; i++) VAR_LOOP(nv) U[nv][k][j][i] = u[i][nv];
#else      
    PrimToCons (v, U[k][j], ibeg, iend);
#endif
  }}
  g_dir = current_dir; /* restore current direction */

}
