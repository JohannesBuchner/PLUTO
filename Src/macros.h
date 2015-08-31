/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO header file for function-like macros.

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 31, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */

/* #####################################################################
    1. Macros that can be placed at any point in the code
   ##################################################################### */

/*! \name Spatial loop macros.
    The following macros provide a compact way to perform 1D or multi-D
    loops in selected regions of the (local) computational domain.
    The \c *BEG_LOOP and \c *END_LOOP macros are used to loop in the 
    leftmost or rightmost boundary ghost zones in the corresponding
    direction \c I, \c J or \c K.
    The \c *DOM_LOOP macros are used to loop inside the computational
    domain (boundaries excluded) while the \c *TOT_LOOP macros are
    used to loop across the entire domain (inside+boundary).
*/
/**@{ */
#define IBEG_LOOP(i)  for ((i) = IBEG; (i)--;    )  
#define JBEG_LOOP(j)  for ((j) = JBEG; (j)--;    )  
#define KBEG_LOOP(k)  for ((k) = KBEG; (k)--;    )

#define IEND_LOOP(i)  for ((i) = IEND + 1; (i) < NX1_TOT; (i)++)
#define JEND_LOOP(j)  for ((j) = JEND + 1; (j) < NX2_TOT; (j)++)
#define KEND_LOOP(k)  for ((k) = KEND + 1; (k) < NX3_TOT; (k)++)

#define IDOM_LOOP(i)  for ((i) = IBEG; (i) <= IEND; (i)++)
#define JDOM_LOOP(j)  for ((j) = JBEG; (j) <= JEND; (j)++)
#define KDOM_LOOP(k)  for ((k) = KBEG; (k) <= KEND; (k)++)

#define ITOT_LOOP(i)  for ((i) = 0; (i) < NX1_TOT; (i)++)
#define JTOT_LOOP(j)  for ((j) = 0; (j) < NX2_TOT; (j)++)
#define KTOT_LOOP(k)  for ((k) = 0; (k) < NX3_TOT; (k)++)

#define DOM_LOOP(k,j,i) KDOM_LOOP(k) JDOM_LOOP(j) IDOM_LOOP(i)

#define TOT_LOOP(k,j,i) KTOT_LOOP(k) JTOT_LOOP(j) ITOT_LOOP(i)

#define X1_BEG_LOOP(k,j,i) KTOT_LOOP(k) JTOT_LOOP(j) IBEG_LOOP(i)
#define X2_BEG_LOOP(k,j,i) KTOT_LOOP(k) JBEG_LOOP(j) ITOT_LOOP(i)
#define X3_BEG_LOOP(k,j,i) KBEG_LOOP(k) JTOT_LOOP(j) ITOT_LOOP(i)

#define X1_END_LOOP(k,j,i) KTOT_LOOP(k) JTOT_LOOP(j) IEND_LOOP(i)
#define X2_END_LOOP(k,j,i) KTOT_LOOP(k) JEND_LOOP(j) ITOT_LOOP(i)
#define X3_END_LOOP(k,j,i) KEND_LOOP(k) JTOT_LOOP(j) ITOT_LOOP(i)

#define TRANSVERSE_LOOP(indx, ip, i,j,k) \
 if (g_dir == IDIR) {ip = &i; indx.pt1 = &j; indx.pt2 = &k;}  \
 if (g_dir == JDIR) {ip = &j; indx.pt1 = &i; indx.pt2 = &k;}  \
 if (g_dir == KDIR) {ip = &k; indx.pt1 = &i; indx.pt2 = &j;}  \
 for (*(indx.pt2) = indx.t2_beg; *(indx.pt2) <= indx.t2_end; *(indx.pt2) += 1) \
 for (*(indx.pt1) = indx.t1_beg; *(indx.pt1) <= indx.t1_end; *(indx.pt1) += 1)
/**@} */

/*! The BOX_LOOP() macro implements a loop over (i,j,k) in a rectangular 
    portion of the domain with indices defined by the (pointer to)
    RBox structure B. 
    The loop increments (di,dj,dk) are members of the structure which 
    are here initialized to either 1 or -1 depending on whether the 
    lower corner index lies below or above the upper index 
    (e.g. B->ib <= B->ie or not). 
*/
#define BOX_LOOP(B,k,j,i) \
 for ((B)->dk = ((k=(B)->kb) <= (B)->ke ? 1:-1); k != (B)->ke+(B)->dk; k += (B)->dk)\
 for ((B)->dj = ((j=(B)->jb) <= (B)->je ? 1:-1); j != (B)->je+(B)->dj; j += (B)->dj)\
 for ((B)->di = ((i=(B)->ib) <= (B)->ie ? 1:-1); i != (B)->ie+(B)->di; i += (B)->di)

/*! The FOR_EACH(p, beg, intList) macro implements a loop over the
    elements of the array \c intList->indx starting at \c beg (in a
    similar way to Python lists).

    Example:
    \code
      intList list;
      list.nvar = 3;
      list.indx[0] = 2;
      list.indx[1] = 5;
      list.indx[2] = 17;
      FOR_EACH(nv, 0, list) printf ("value is = %d\n",nv);
    \endcode
*/
#define FOR_EACH(nv, beg, list)  \
  for ((list)->i = beg, nv = (list)->indx[beg]; \
       (list)->i < (list)->nvar; \
       nv = (list)->indx[++((list)->i)])


/*! Faster implementation than stdlib floor() function.
    It returns the largest integer value less than or equal to z.
*/
#define INT_FLOOR(z)   ((int)((z) + 32768.) - 32768)

/*! Return the maximum between two numbers. */
#define MAX(a,b)  ( (a) >= (b) ? (a) : (b) ) 

/*! Return the minimum between two numbers. */
#define MIN(a,b)  ( (a) <= (b) ? (a) : (b) ) 
      
/*! Return the number with the smaller absolute value. */
#define ABS_MIN(a,b)  (fabs(a) < fabs(b) ? (a):(b)) 
                         
/*! Return the sign of x. */
#define DSIGN(x)      ( (x) >= 0.0 ? (1.0) : (-1.0))

#define MINMOD(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
#define VAN_LEER(a,b) ((a)*(b) > 0.0 ? 2.0*(a)*(b)/((a)+(b)):0.0)
#define MC(a,b) (MINMOD(0.5*((a)+(b)), 2.0*MINMOD((a),(b))))
#define SWAP_VAR(x) SwapEndian(&x, sizeof(x));

/*! Exit macro. For Chombo it is defined elsewhere. */
#ifdef PARALLEL
 #define QUIT_PLUTO(e_code)   \
        {MPI_Abort(MPI_COMM_WORLD, e_code);MPI_Finalize(); exit(e_code);}
#elif (defined CH_MPI)
 #define QUIT_PLUTO(e_code)   \
        {MPI_Abort(MPI_COMM_WORLD, e_code); exit(e_code);}
#else
 #define QUIT_PLUTO(e_code)   exit(e_code);
#endif

/* #####################################################################
    2. Macros that must be placed only *AFTER* definitions.h has been 
       included
   ##################################################################### */

/* *******************************************************
    Expand dimension- or component-dependent expressions
   *******************************************************  */

/*! \def EXPAND(a,b,c)
    Allows to write component-independent expressions involving vectors 
    by evaluating as many arguments as the value of COMPONENTS. 
    The result is that only the first argument will be compiled in 1D, 
    the first two arguments in 2D and all of them in 3D. 
    As an example: 
    \code
    EXPAND( v[VX1] =  1.0;  ,
            v[VX2] =  5.0;  , 
            v[VX3] = -4.0; )
    \endcode
    becomes
    \code
     v[VX1] = 1.0;   
    \endcode
    when \c COMPONENTS is equal to 1 or
    \code
     v[VX1] = 1.0;  
     v[VX2] = 5.0;
    \endcode
    when \c COMPONENTS is equal to 2 or
    \code
     v[VX1] =  1.0;  
     v[VX2] =  5.0;
     v[VX3] = -4.0;
    \endcode
    when \c COMPONENTS is equal to 3.
*/    
    
/*! \def D_EXPAND(a,b,c)
    Similar to the EXPAND() macro but the expansion depends on DIMENSIONS.  */

/*! \def SELECT(a,b,c)
    Expand only the 1st, 2nd or 3rd argument based on the value of
    COMPONENTS.                                                       */

/*! \def D_SELECT(a,b,c)
    Expand only the 1st, 2nd or 3rd argument based on the value of
    DIMENSIONS.                                                       */

#if COMPONENTS == 1
 #define EXPAND(a,b,c) a
 #define SELECT(a,b,c) a
#endif

#if COMPONENTS == 2
 #define EXPAND(a,b,c) a b
 #define SELECT(a,b,c) b
#endif

#if COMPONENTS == 3
 #define EXPAND(a,b,c) a b c
 #define SELECT(a,b,c) c
#endif

#if DIMENSIONS == 1
 #define D_EXPAND(a,b,c)  a
 #define D_SELECT(a,b,c)  a
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 0
 #define KOFFSET 0
#endif

#if DIMENSIONS == 2
 #define D_EXPAND(a,b,c) a b
 #define D_SELECT(a,b,c) b
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 1
 #define KOFFSET 0
#endif

#if DIMENSIONS == 3
 #define D_EXPAND(a,b,c) a b c
 #define D_SELECT(a,b,c) c
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 1
 #define KOFFSET 1
#endif

#if WARNING_MESSAGES == YES
 #define WARNING(a)  a
#else
 #define WARNING(a)
#endif

#define DOT_PRODUCT(a,b)  (EXPAND((a)[0]*(b)[0], + (a)[1]*(b)[1], + (a)[2]*(b)[2]))



#define VAR_LOOP(n)   for ((n) = NVAR; (n)--;    )
#define DIM_LOOP(d)   for ((d) = 0; (d) < DIMENSIONS; (d)++)



/* -- some new macros.
   CDIFF: Central DIFFerencing
   FDIFF: Forward DIFFerencing
   
*/
 #define FDIFF_X1(b,k,j,i) (b[k][j][i+1] - b[k][j][i])
 #define FDIFF_X2(b,k,j,i) (b[k][j+1][i] - b[k][j][i])
 #define FDIFF_X3(b,k,j,i) (b[k+1][j][i] - b[k][j][i])

 #define CDIFF_X1(b,k,j,i) 0.5*(b[k][j][i+1] - b[k][j][i-1])
 #define CDIFF_X2(b,k,j,i) 0.5*(b[k][j+1][i] - b[k][j-1][i])
 #define CDIFF_X3(b,k,j,i) 0.5*(b[k+1][j][i] - b[k-1][j][i])

/* -- in 2D we set the derivative with respect to the third coordinate = 0 */

#if DIMENSIONS == 2
 #undef CDIFF_X3
 #undef FDIFF_X3

 #define CDIFF_X3(b,i,j,k)  0.0
 #define FDIFF_X3(b,i,j,k)  0.0
#endif

/* ---- define average macros ---- */

#define AVERAGE_X(q,k,j,i)   0.5*(q[k][j][i] + q[k][j][i+1])
#define AVERAGE_Y(q,k,j,i)   0.5*(q[k][j][i] + q[k][j+1][i])
#define AVERAGE_Z(q,k,j,i)   0.5*(q[k][j][i] + q[k+1][j][i])

#define AVERAGE_XY(q,k,j,i)   0.25*(  q[k][j][i]   + q[k][j][i+1] \
                                    + q[k][j+1][i] + q[k][j+1][i+1])
#define AVERAGE_XZ(q,k,j,i)   0.25*(  q[k][j][i]   + q[k][j][i+1] \
                                    + q[k+1][j][i] + q[k+1][j][i+1])
#define AVERAGE_YZ(q,k,j,i)   0.25*(  q[k][j][i]   + q[k][j+1][i] \
                                    + q[k+1][j][i] + q[k+1][j+1][i])

/* -- re-define the macros for 1 or 2 dimensions -- */

#if DIMENSIONS == 1

 #undef AVERAGE_Y
 #undef AVERAGE_Z

 #undef AVERAGE_XY
 #undef AVERAGE_XZ
 #undef AVERAGE_YZ

 #define AVERAGE_Y(q,k,j,i)    (q[0][0][i])
 #define AVERAGE_Z(q,k,j,i)    (q[0][0][i])

 #define AVERAGE_XY(q,k,j,i)   AVERAGE_X(q,0,0,i)
 #define AVERAGE_XZ(q,k,j,i)   AVERAGE_X(q,0,0,i)
 #define AVERAGE_YZ(q,k,j,i)   (q[0][0][i])

#elif DIMENSIONS == 2

 #undef AVERAGE_Z
 #undef AVERAGE_XZ
 #undef AVERAGE_YZ

 #define AVERAGE_Z(q,k,j,i)    (q[0][j][i])
 #define AVERAGE_XZ(q,k,j,i)   0.5*(q[0][j][i] + q[0][j][i+1])
 #define AVERAGE_YZ(q,k,j,i)   0.5*(q[0][j][i] + q[0][j+1][i])

#endif

