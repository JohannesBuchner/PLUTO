/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Miscellaneous functions for implementing the shearing-box
         boundary condition.

  The shearing-box tool file contains various functions for the 
  implementation of the shearingbox boundary conditions in serial or
  parallel mode.
  The SB_SetBoundaryVar() function applies shearing-box boundary conditions 
  to a 3D array U[k][j][i] at an X1_BEG or X1_END boundary.
  The array U[k][j][i] is defined on the RBox *box with grid indices
  <tt>
  (box->ibeg) <= i <= (box->iend), (box->jbeg) <= j <= (box->jend)
  (box->kbeg) <= k <= (box->kend),
  </tt> 
  and assumes that periodic boundary conditions have already been set.
  Indices in the y-directions are not necessary, since the 
  entire y-range is assumed.
 
  In parallel mode, we use the SB_ShiftBoundaryVar() function to perform 
  the integer shift of cells across the processors, while
  interpolation function SB_ShearingInterp() handles the fractional 
  part only.
  In serial or when there's only 1 processor along y, the 
  interpolation function does all the job (integer+fractional).
    
  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.musicanisi@cineca.it)
  \date    Nov 13, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double sb_Ly;

/* ********************************************************************* */
void SB_SetBoundaryVar(double ***U, RBox *box, int side, double t, Grid *grid)
/*!
 * Fill ghost zones using shearing-box conditions.
 *
 * \param [out] U   pointer to a 3D array (centered or staggered)
 * \param [in]  box pointer to RBox structure defining the domain 
 *                  sub-portion over which shearing-box conditions
 *                  have to be applied
 * \param [in] side the side of the X1 boundary (X1_BEG or X1_END)
 * \param [in] t    the simulation time
 * \param [in] grid  pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  int i, j, k, ngh_y;
  int nghL, nghR; 
  int coords[3], dst1, dst2;
  long int count, buf_size;
  static double *qL, *qR;

  if (side != X1_BEG && side != X1_END){
    print ("! SB_SetBoundaryVar: wrong boundary\n");
    QUIT_PLUTO(1);
  }

  if (qL == NULL){
    qL = ARRAY_1D(NMAX_POINT, double);
    qR = ARRAY_1D(NMAX_POINT, double);
  }

  sb_Ly = g_domEnd[JDIR] - g_domBeg[JDIR];

/* -- Get number of ghost zones to the left and to the right --*/

  nghL = JBEG - box->jbeg;
  nghR = box->jend - JEND;

/* -- Shift data values across parallel domains -- */

#ifdef PARALLEL
  if (grid->nproc[JDIR] > 1) SB_ShiftBoundaryVar(U, box, side, t, grid);
#endif

/* -- Exchange values between processors to fill ghost zones -- */

  SB_FillBoundaryGhost(U, box, nghL, nghR, grid);

/* ---- Perform 1D interpolation in the x 
        boundary zones along the y direction  ---- */

  if (side == X1_BEG){
    for (k = box->kbeg; k <= box->kend; k++){
    for (i = box->ibeg; i <= box->iend; i++){
      JTOT_LOOP(j) qR[j] = U[k][j][i];
      SB_ShearingInterp (qL, qR, t, X1_BEG, grid);
      JDOM_LOOP(j) U[k][j][i] = qL[j];
    }}
  }else if (side == X1_END){
    for (k = box->kbeg; k <= box->kend; k++){
    for (i = box->ibeg; i <= box->iend; i++){
      JTOT_LOOP(j) qL[j] = U[k][j][i];
      SB_ShearingInterp (qL, qR, t, X1_END, grid);
      JDOM_LOOP(j) U[k][j][i] = qR[j];
    }}
  }

/* -- Exchange values between processors to fill ghost zones -- */

  SB_FillBoundaryGhost (U, box, nghL, nghR, grid);
  return;
}

#ifdef PARALLEL
/* ********************************************************************* */
void SB_ShiftBoundaryVar(double ***q, RBox *box, int side, double t, Grid *grid)
/*!
 * Split the 3D array q[k][j][i] in two buffers and send them 
 * to processors with rank dst1 and dst2.
 * The box structure contains the original grid index ranges on
 * top of which q is defined.
 * At the same time, receive buffers from processors with rank 
 * src1,src2.
 *
 * \param [in,out] q  3D array
 * \param [in]        box the rectangular box giving the index range
 * \param [in]        side the side of the X1 boundary
 * \param [in] t      simulation time
 * \param [in] grid   pointer to an array of Grid structures
 *
 * \return This function has no return value.
 ************************************************************************* */
{
  int    i, j, k, ngh_x, ngh_y;
  int    nx_buf1, ny_buf1, nz_buf1;
  int    nx_buf2, ny_buf2, nz_buf2;
  long int count, buf1_size, buf2_size;

  int coords[3];
  int dst1, dst2, src1, src2;   
  int Delta_j, shift_coordy_mod, shift_coordy_div; 

  double scrh;
  static double *snd_buf, *rcv_buf;

  RBox buf1, buf2;
  RBox *pbuf1, *pbuf2;

  static MPI_Comm cartcomm;
  MPI_Status status;

/*  -------------------------------------------------------------------- */
/*! We allocate static memory areas for send/receive buffers and
    employ just one send buffer and one receive buffer with 
    size equals the full extent of the boundary side. 
    The two data chunks coming from q[k][j][i] are stored at 
    different positions in the send/recv buffers.

    In 3D staggered MHD we augment the buffer size by 1 point 
    in the z-direction for BZs. 
    There's no need to do the same for BYs periodic condition 
    will be applied later on.                                            */
/*  -------------------------------------------------------------------- */

  ngh_x = grid->nghost[IDIR];
  ngh_y = grid->nghost[JDIR];
  if (snd_buf == NULL){
    int nz = NX3_TOT;
    #if defined STAGGERED_MHD && DIMENSIONS == 3
     nz += 1;
    #endif
    snd_buf = ARRAY_1D(nz*NX2_TOT*ngh_x, double);
    rcv_buf = ARRAY_1D(nz*NX2_TOT*ngh_x, double);
     
    AL_Get_cart_comm(SZ, &cartcomm);
  }

/* -- compute buffer index offsets in x/z -- */

  buf1.ibeg = buf2.ibeg = box->ibeg;
  buf1.iend = buf2.iend = box->iend;
  buf1.kbeg = buf2.kbeg = box->kbeg;
  buf1.kend = buf2.kend = box->kend;
    
/* --------------------------------------------------------------------- */
/*! Depending on the boundary side, we set buffer offsets in the 
    y-direction and find the ranks of the processors to/from which we 
    must send/receive data.
    Local processor sends data to processors dst1 and dst2
    and receive data from processors src1 and src2.
    \note parallel decomposition must result in equally sized
          grids in the y-direction.
    \verbatim
    |________________|      |_____[dst1]_____||____[dst2]______|
     <----><-------->                  <---->  <-------->
      buf1    buf2                      buf1      buf2
    \endverbatim
*/
/* --------------------------------------------------------------------- */

  scrh    = (fmod(sb_vy*t, sb_Ly))/grid->dx[JDIR][JBEG]; 
  Delta_j = (int)scrh;                   

  shift_coordy_mod = Delta_j%NX2; 
  shift_coordy_div = Delta_j/NX2; 

  if (side == X1_BEG){

    ny_buf1 = NX2 - (shift_coordy_mod) + ngh_y;
    ny_buf2 = ngh_y + (shift_coordy_mod); 

    buf1.jbeg = 0; buf1.jend = ny_buf1 - 1;
    buf2.jbeg = 0; buf2.jend = ny_buf2 - 1;

    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[JDIR] += shift_coordy_div;
    MPI_Cart_rank(cartcomm, coords, &dst1);
    
    coords[JDIR] += 1;
    MPI_Cart_rank(cartcomm, coords, &dst2);
    
    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[JDIR] -= shift_coordy_div;
    MPI_Cart_rank(cartcomm, coords, &src1);

    coords[JDIR] -= 1;
    MPI_Cart_rank(cartcomm, coords, &src2);

  }else if (side == X1_END){

    ny_buf1 = shift_coordy_mod + ngh_y;
    ny_buf2 = ngh_y + (NX2 - shift_coordy_mod);

    buf1.jbeg = 0; buf1.jend = ny_buf1 - 1;
    buf2.jbeg = 0; buf2.jend = ny_buf2 - 1;

    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[JDIR] -= shift_coordy_div;
    MPI_Cart_rank(cartcomm, coords, &dst2);
    
    coords[JDIR] -= 1;
    MPI_Cart_rank(cartcomm, coords, &dst1);
    
    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[JDIR] += shift_coordy_div;
    MPI_Cart_rank(cartcomm, coords, &src2);

    coords[JDIR] += 1;
    MPI_Cart_rank(cartcomm, coords, &src1);

  }else{
    print ("! SB_ShiftBoundaryVar: wrong boundary\n");
    QUIT_PLUTO(1);
  }

/* -- Compute buffer sizes in the three directions -- */

  nx_buf1 = buf1.iend - buf1.ibeg + 1;
  nx_buf2 = buf2.iend - buf2.ibeg + 1;

  ny_buf1 = buf1.jend - buf1.jbeg + 1;
  ny_buf2 = buf2.jend - buf2.jbeg + 1;

  nz_buf1 = buf1.kend - buf1.kbeg + 1;
  nz_buf2 = buf2.kend - buf2.kbeg + 1;

/* -- Total buffer size -- */

  buf1_size = nz_buf1*ny_buf1*nx_buf1;
  buf2_size = nz_buf2*ny_buf2*nx_buf2;

/* -- Fill send buffers with values -- */

  pbuf1 = &buf1;  /* pointer to RBox are used   */
  pbuf2 = &buf2;  /* inside the BOX_LOOP macros */

  count = 0;
  BOX_LOOP(pbuf1, k, j, i) snd_buf[count++] = q[k][JBEG+j][i];
  BOX_LOOP(pbuf2, k, j, i) snd_buf[count++] = q[k][JEND-ny_buf2+1+j][i];

/* -- Send/Receive  buffers -- */

  MPI_Sendrecv(snd_buf, buf1_size, MPI_DOUBLE, dst1, 1, 
               rcv_buf, buf1_size, MPI_DOUBLE, src1, 1,
               cartcomm, &status);

  MPI_Sendrecv(snd_buf + buf1_size, buf2_size, MPI_DOUBLE, dst2, 2, 
               rcv_buf + buf1_size, buf2_size, MPI_DOUBLE, src2, 2, 
               cartcomm, &status);

/* -- Store received buffers in the correct locations -- */

  count = 0;
  BOX_LOOP(pbuf1, k, j, i) q[k][j+ny_buf2][i] = rcv_buf[count++];   
  BOX_LOOP(pbuf2, k, j, i) q[k][j][i] = rcv_buf[count++];
}
#endif /* PARALLEL */

/* ********************************************************************* */
void SB_FillBoundaryGhost(double ***U, RBox *box, 
                          int nghL, int nghR, Grid *grid)
/*!
 *  Fill ghost zones in the Y direction in the X1_BEG and 
 *  X1_END boundary regions.
 *
 * \param [in,out] U a 3D data array
 * \param [in,out] box the RBox structure containing the grid indices in 
 *                     the x and z directions. Indices in the y-directions 
 *                     are reset here for convenience. 
 * \param [in] nghL number of ghost zones on the left 
 * \param [in] nghR number of ghost zones on the right
 * \param [in] grid pointer to an 
 * \note nghL and nghR can be different for staggered fields like BY.
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  int i, j, k;
  int jb0, je0; 
#ifdef PARALLEL
  int coords[3];
  long int count, buf_size;
  static double *snd_buf1, *snd_buf2, *rcv_buf1, *rcv_buf2;
  static int dst1, dst2;
  static MPI_Comm cartcomm;
  MPI_Status status;
#endif

/* -------------------------------------------------------
    Store a copy of the original grid indices in the y-dir
   ------------------------------------------------------- */
 
  jb0 = box->jbeg;
  je0 = box->jend;

/* ------------------------------------------------------------------
    Impose periodic b.c. if there's only 1 processor along y
   ------------------------------------------------------------------ */

  if (grid->nproc[JDIR] == 1){
    box->jbeg = JBEG - nghL; 
    box->jend = JBEG - 1;
    BOX_LOOP(box, k, j, i) U[k][j][i] = U[k][j + NX2][i];

    box->jbeg = JEND + 1; 
    box->jend = JEND + nghR;
    BOX_LOOP(box, k, j, i) U[k][j][i] = U[k][j - NX2][i];

    box->jbeg = jb0; box->jend = je0;
    return;
  }
 
#ifdef PARALLEL

/* ----------------------------------------------------------------
    - Allocate memory for maximum buffer size;
    - Get destination ranks of the processors lying above and 
      below current processor. 
      This is done only once at the beginning since parallel
      decomposition is not going to change during the computation.
   ---------------------------------------------------------------- */

  if (snd_buf1 == NULL){
   
    i = grid->nghost[IDIR];
    j = grid->nghost[JDIR];
    k = NX3_TOT;
    #ifdef STAGGERED_MHD 
     j++;
     k++;
    #endif
    snd_buf1 = ARRAY_1D(i*j*k, double);
    snd_buf2 = ARRAY_1D(i*j*k, double);
    rcv_buf1 = ARRAY_1D(i*j*k, double);
    rcv_buf2 = ARRAY_1D(i*j*k, double);

    AL_Get_cart_comm(SZ, &cartcomm);

    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[JDIR] -= 1;
    MPI_Cart_rank(cartcomm, coords, &dst1);

    coords[JDIR] += 2;
    MPI_Cart_rank(cartcomm, coords, &dst2);
    
  }

/* ----------------------------------------------------------------- 
     Send/receive buffers. 
     Buffer snd_buf1 is sent to rank [dst1] while snd_buf2 is sent
     to rank [dst2] as follows:
     

     |_____[dst1]_____|      |________________|      |_____[dst2]_____|
                              <--->      <--->     
                       <----- buf1        buf2  ---->
   ----------------------------------------------------------------- */

/* -- Send buffer at JBEG -- */

  count = 0; 
  box->jbeg = JBEG; 
  box->jend = JBEG + nghR - 1;
  buf_size =  (box->kend - box->kbeg + 1)
             *(box->jend - box->jbeg + 1)
             *(box->iend - box->ibeg + 1);
  BOX_LOOP(box, k, j, i) snd_buf1[count++] = U[k][j][i];

  MPI_Sendrecv(snd_buf1, buf_size, MPI_DOUBLE, dst1, 1, 
               rcv_buf2, buf_size, MPI_DOUBLE, dst2, 1,
               cartcomm, &status);

/* -- Send buffer at JEND -- */

  count = 0; 
  box->jbeg = JEND - nghL + 1; 
  box->jend = JEND; 
  buf_size =  (box->kend - box->kbeg + 1)
             *(box->jend - box->jbeg + 1)
             *(box->iend - box->ibeg + 1);
  BOX_LOOP(box, k, j, i) snd_buf2[count++] = U[k][j][i];

  MPI_Sendrecv(snd_buf2, buf_size, MPI_DOUBLE, dst2, 2, 
               rcv_buf1, buf_size, MPI_DOUBLE, dst1, 2, 
               cartcomm, &status);

/* -- Place buffers in the correct position -- */

  count = 0; 
  box->jbeg = JEND + 1; 
  box->jend = JEND + nghR;
  BOX_LOOP(box, k, j, i) U[k][j][i] = rcv_buf2[count++];
  
  count = 0; 
  box->jbeg = JBEG - nghL; 
  box->jend = JBEG - 1;
  BOX_LOOP(box, k, j, i) U[k][j][i] = rcv_buf1[count++];

/* -- Restore original grid index in the y-dir -- */

  box->jbeg = jb0; box->jend = je0;

#endif /* PARALLEL */
}

/* ********************************************************************* */
void SB_ShearingInterp (double *qL, double *qR, double t, int side,
                         Grid *grid)
/*!
 *  Shift the 1D arrays qL or qR by an amount \f$ wt = 
 *  \Delta y\Delta J + \epsilon \f$ (if there's only one processor in y) 
 *  or just by \f$ \epsilon \f$ (otherwise) by advecting the initial
 *  profile of qR (when side == X1_BEG) or qL (when side == X1_END).
 * 
 * \param [in,out] qL 1D array containing the un-shifted values on the 
 *                 left boundary
 * \param [in,out] qR 1D array containing the un-shifted values on the
 *                 right boundary
 * \param [in] t simulation time
 * \param [in] side boundary side
 * \param [in] grid pointer to an array of Grid structures
 * 
 * \return This function has no return value.
 *********************************************************************** */
{
  int   j, Delta_j, nghost;
  int   jS, jN, jp, jm;
  int   isgn;
  double  Delta_L, Delta_y, dyfrac;
  double  scrh;
  double  *q_from, *q_to;
  static double *dq, *dql, *flux;

  if (dq == NULL){
    dq   = ARRAY_1D(NMAX_POINT, double);
    dql  = ARRAY_1D(NMAX_POINT, double);
    flux = ARRAY_1D(NMAX_POINT, double);
  }   

  Delta_y = grid->dx[JDIR][JBEG];
  nghost  = grid->nghost[JDIR];

  if (side == X1_BEG){
    q_from = qR; 
    q_to   = qL;
    isgn   = -1;
  }else if (side == X1_END){
    q_from = qL; 
    q_to   = qR;
    isgn   = 1;
  }  
  
  Delta_L = fmod(isgn*sb_vy*t, sb_Ly);
  scrh    = Delta_L/Delta_y;

/* ---- shift index and fractional remainder inside cell ---- */

  Delta_j = (int)scrh;                   
  dyfrac  = fabs(scrh) - floor(fabs(scrh));    

/* -- in parallel mode, shift has alrady been done at this point -- */

  #ifdef PARALLEL
   if (grid->nproc[JDIR] > 1) Delta_j = 0;
  #endif
 
/* -- compute limited slopes & interpolate -- */

  for (j = 1; j <= JEND + nghost; j++) {
    dq[j] = q_from[j] - q_from[j - 1];
  }

  #if SB_ORDER == 1  /* -- 1st order interpolation -- */

   for (j = 0; j <= JEND + nghost; j++) dql[j] = 0.0;

   if (side == X1_BEG)
     for (j = JBEG - 1; j <= JEND; j++) flux[j] = q_from[j];
   else if (side == X1_END) 
     for (j = JBEG - 1; j <= JEND; j++) flux[j] = q_from[j + 1];

  #elif SB_ORDER == 2  /* -------------------------------------------
                           2nd order interpolation uses the standard
                           conservative MUSCL-Hancock scheme 
                          ------------------------------------------ */
/*   vanleer_lim (dq + 1, dq, dql, JBEG - 1, JEND + 1, grid);  */
   for (j = JBEG - 1; j <= JEND + 1; j++) dql[j] = VAN_LEER(dq[j+1], dq[j]);

   if (side == X1_BEG){
     for (j = JBEG - 1; j <= JEND; j++) {
       flux[j] = q_from[j] + 0.5*dql[j]*(1.0 - dyfrac);
     }
   }else if (side == X1_END){
     for (j = JBEG - 1; j <= JEND; j++) {
       flux[j] = q_from[j + 1] - 0.5*dql[j + 1]*(1.0 - dyfrac);
     }
   }

  #elif SB_ORDER == 3  /* --------------------------------------
                            3rd order interpolation uses the 
                            traditional PPM method.
                          -------------------------------------- */
   static double *qp, *qm;
   double  ap, am, P0, P1, P2;
    
   if (qp == NULL){
     qp = ARRAY_1D(NMAX_POINT, double);
     qm = ARRAY_1D(NMAX_POINT, double);
     if ( (JBEG - 3) < 0) {
        printf ("! SB_ShearingInterp: the value of SB_ORDER (%d)\n",
                 SB_ORDER);
        printf ("! is not consistent with underlying algorithm\n");
        QUIT_PLUTO(1);
     }  
   }       
    
/*    mc_lim (dq + 1, dq, dql, JBEG - 2, JEND + 2, grid);  */
   for (j = JBEG - 2; j <= JEND + 2; j++) dql[j] = MC(dq[j+1], dq[j]);
   for (j = JBEG - 1; j <= JEND + 1; j++){
     qp[j] = 0.5*(q_from[j] + q_from[j+1]) - (dql[j+1] - dql[j])/6.0;
     qm[j] = 0.5*(q_from[j] + q_from[j-1]) - (dql[j]   - dql[j-1])/6.0;

     ap = qp[j] - q_from[j];
     am = qm[j] - q_from[j];

     if (ap*am >= 0.0) {
       ap = am = 0.0;
     }else{
       if (fabs(ap) >= 2.0*fabs(am)) ap = -2.0*am;
       if (fabs(am) >= 2.0*fabs(ap)) am = -2.0*ap;
     }
     qp[j] = q_from[j] + ap;
     qm[j] = q_from[j] + am;
   }
      
   if (side == X1_BEG){
     for (j = JBEG - 1; j <= JEND; j++){       
       P0 = 1.5*q_from[j] - 0.25*(qp[j] + qm[j]);
       P1 = qp[j] - qm[j];
       P2 = 3.0*(qp[j] - 2.0*q_from[j] + qm[j]);
       flux[j] = P0 + 0.5*P1*(1.0 - dyfrac) +
                 P2*(0.25 - 0.5*dyfrac + dyfrac*dyfrac/3.0);
     }
   }else if (side == X1_END) {
     for (j = JBEG - 1; j <= JEND; j++){       
       P0 = 1.5*q_from[j + 1] - 0.25*(qp[j + 1] + qm[j + 1]);
       P1 = qp[j + 1] - qm[j + 1];
       P2 = 3.0*(qp[j + 1] - 2.0*q_from[j + 1] + qm[j + 1]);
       flux[j] = P0 - 0.5*P1*(1.0 - dyfrac) +
                 P2*(0.25 - 0.5*dyfrac + dyfrac*dyfrac/3.0);
     }
   }       
   
  #else
  
   print ("! SB_ORDER should lie between 1 and 3\n");
   QUIT_PLUTO(1);
   
  #endif      

/* -------------------------------------------------------------
     Assign volumetric sliding averages. 
     The match point falls on the opposite side 
     between jS and jN = jS + 1.

     At t = 0:    

      - point j on the left  matches jN on the right;
      - point j on the right matches jS on the left;

    In any case, 

      (j + Delta_j) --> jN if side is X1_BEG
      (j + Delta_j) --> jS if side is X1_BEG

    and the advection operator can be defined in terms
    of jup = (j + Delta_j):

     q_to[j] = q_from[jup] -+ dyfrac*(flux[jup] - flux[jup-1]);

    with - (+) for left (right) if dyfrac is always positive.
   ------------------------------------------------------------- */

  if (side == X1_END) dyfrac *= -1.0;

  if (grid->nproc[JDIR] == 1){
    for (j = JBEG; j <= JEND; j++){
      jp = SB_JSHIFT(j + Delta_j);
      jm = SB_JSHIFT(jp - 1);
      q_to[j] = q_from[jp] - dyfrac*(flux[jp] - flux[jm]);
    }
  }else{
    for (j = JBEG; j <= JEND; j++){
      jp = j;
      jm = j-1;
      q_to[j] = q_from[jp] - dyfrac*(flux[jp] - flux[jm]);
    }
  }
}
/* ********************************************************************* */
int SB_JSHIFT (int jc)
/*!
 * Make sure jc always fall  between JBEG and JEND
 *
 *********************************************************************** */
{
  int  jj;

  jj = jc;
  if (jj > JEND) jj -= (int)NX2;
  if (jj < JBEG) jj += (int)NX2;

  if (jj < JBEG || jj > JEND){
    print (" ! j index out of range in JSHIFT  %d, %d\n", jc, jj);
    QUIT_PLUTO(1);
  }

  return(jj);
}

