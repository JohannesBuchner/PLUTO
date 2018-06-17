/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO header file for function-like macros.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 15, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */

#ifndef DEBUG
  #define DEBUG FALSE
#endif

#if DEBUG == TRUE
  #define DEBUG_FUNC_BEG(a)                          \
    char d_funcName[64];                             \
    sprintf (d_funcName,"%s",a);                     \
    d_indent += 2;                                   \
    if (d_condition) { print("%*c", d_indent, ' ');  \
                       print (">>[%s]\n",a );  }

  #define DEBUG_FUNC_END(a)  \
    if (d_condition) {print("%*c", d_indent, ' ');  \
                     print("<<[%s]\n",a ); }        \
     d_indent -= 2;                           
  
  #define DEBUG_FUNC_NAME  d_funcName

#else
  #define DEBUG_FUNC_BEG(a)
  #define DEBUG_FUNC_END(a)
  #define DEBUG_FUNC_NAME "Not Set"
#endif
  
  

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
 for ((B)->dk = ((k=(B)->kbeg) <= (B)->kend ? 1:-1); k != (B)->kend+(B)->dk; k += (B)->dk)\
 for ((B)->dj = ((j=(B)->jbeg) <= (B)->jend ? 1:-1); j != (B)->jend+(B)->dj; j += (B)->dj)\
 for ((B)->di = ((i=(B)->ibeg) <= (B)->iend ? 1:-1); i != (B)->iend+(B)->di; i += (B)->di)


#define IBOX_LOOP(B,i) \
 for ((B)->di = ((i=(B)->ibeg) <= (B)->iend ? 1:-1); i != (B)->iend+(B)->di; i += (B)->di)

#define JBOX_LOOP(B,j) \
 for ((B)->dj = ((j=(B)->jbeg) <= (B)->jend ? 1:-1); j != (B)->jend+(B)->dj; j += (B)->dj)

#define KBOX_LOOP(B,k) \
 for ((B)->dk = ((k=(B)->kbeg) <= (B)->kend ? 1:-1); k != (B)->kend+(B)->dk; k += (B)->dk)


#define BOX_TRANSVERSE_LOOP(box,k,j,i)                                     \
  if      (g_dir == IDIR) {(box)->n = &i; (box)->t = &j; (box)->b = &k;}   \
  else if (g_dir == JDIR) {(box)->n = &j; (box)->t = &i; (box)->b = &k;}   \
  else if (g_dir == KDIR) {(box)->n = &k; (box)->t = &i; (box)->b = &j;}   \
  for (*((box)->b) = *(box)->bbeg; *((box)->b) <= *(box)->bend; *((box)->b) +=1) \
  for (*((box)->t) = *(box)->tbeg; *((box)->t) <= *(box)->tend; *((box)->t) +=1)


/*! The FOR_EACH(p, intList) macro implements a loop over the elements
    of the array \c intList->indx (in a similar way to Python lists).

    Example:
    \code
      intList list;
      list.nvar = 3;
      list.indx[0] = 2;
      list.indx[1] = 5;
      list.indx[2] = 17;
      FOR_EACH(nv, &list) printf ("value is = %d\n",nv);
    \endcode
*/
#define FOR_EACH(nv, list)  \
  for ((list)->i = 0, nv = (list)->indx[0]; \
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
/* See
http://stackoverflow.com/questions/11303135/broadcast-message-for-all-processes-to-exitmpi
*/

#ifdef PARALLEL
 #define QUIT_PLUTO(e_code)   \
        {print ("! abort\n"); MPI_Abort(MPI_COMM_WORLD, e_code); \
         MPI_Finalize(); exit(e_code);}
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

/*! \name Forward and central finite differences macros.
    The following set of macros provide a compact way to perform 
    two-point, undivided finite difference operations in a specified
    direction.
    Differences can be either \c forward or \c central.
    For instance, \c FDIFF_X2(Q,k,j,i) will compute a forward
    difference of \c Q in the \c y direction, \c (Q[j+1]-Q[j]), while
    \c CDIFF_X3(Q,k,j,i) will compute a central different approximation
    of \c Q in the \c z direction: \c (Q[k+1] - Q[k-1])/2.
*/
/**@{ */
#define FDIFF_X1(Q,k,j,i)     (Q[k][j][i+1] - Q[k][j][i])
#define CDIFF_X1(Q,k,j,i) (0.5*(Q[k][j][i+1] - Q[k][j][i-1]))

#if DIMENSIONS == 1

 #define FDIFF_X2(Q,k,j,i)  (0.0)
 #define CDIFF_X2(Q,k,j,i)  (0.0)

 #define FDIFF_X3(Q,i,j,k)  (0.0)
 #define CDIFF_X3(Q,i,j,k)  (0.0)

#elif DIMENSIONS == 2

 #define FDIFF_X2(Q,k,j,i)      (Q[k][j+1][i] - Q[k][j][i])
 #define CDIFF_X2(Q,k,j,i) (0.5*(Q[k][j+1][i] - Q[k][j-1][i]))

 #define CDIFF_X3(Q,i,j,k)  0.0
 #define FDIFF_X3(Q,i,j,k)  0.0

#elif DIMENSIONS == 3

 #define FDIFF_X2(Q,k,j,i)      (Q[k][j+1][i] - Q[k][j][i])
 #define CDIFF_X2(Q,k,j,i) (0.5*(Q[k][j+1][i] - Q[k][j-1][i]))

 #define FDIFF_X3(Q,k,j,i)      (Q[k+1][j][i] - Q[k][j][i])
 #define CDIFF_X3(Q,k,j,i) (0.5*(Q[k+1][j][i] - Q[k-1][j][i]))

#endif
/**@} */

/*! \name Spatial averages macros.
    The following set of macros provide a compact way to perform multi-D
    averages from cell centered values to interfaces.
    For instance, \C AVERAGE_X(q,k,j,i) will simply take the
    arithmetic average betwen q(i) and q(i+1) at the i+1/2 interface.
    Likewise, AVERAGE_YZ(q,k,j,i) will produce an average at the
    j+1/2 and k+1/2 edge.
*/
/**@{ */
#define AVERAGE_X(q,k,j,i)   (0.5*(q[k][j][i] + q[k][j][i+1]))

#if DIMENSIONS == 1

  #define AVERAGE_Y(q,k,j,i)    (q[0][0][i])
  #define AVERAGE_Z(q,k,j,i)    (q[0][0][i])

  #define AVERAGE_XY(q,k,j,i)   AVERAGE_X(q,0,0,i)
  #define AVERAGE_XZ(q,k,j,i)   AVERAGE_X(q,0,0,i)
  #define AVERAGE_YZ(q,k,j,i)   (q[0][0][i])
 
  #define AVERAGE_XYZ(q,k,j,i)  0.5*(q[0][0][i] + q[0][0][i+1])

#elif DIMENSIONS == 2


  #define AVERAGE_Y(q,k,j,i)   (0.5*(q[k][j][i] + q[k][j+1][i]))
  #define AVERAGE_Z(q,k,j,i)    (q[0][j][i])

  #define AVERAGE_XY(q,k,j,i)   ( 0.25*(  q[k][j][i]   + q[k][j][i+1] \
                                        + q[k][j+1][i] + q[k][j+1][i+1]) )
  #define AVERAGE_XZ(q,k,j,i)   (0.5*(q[0][j][i] + q[0][j][i+1]))
  #define AVERAGE_YZ(q,k,j,i)   (0.5*(q[0][j][i] + q[0][j+1][i]))

  #define AVERAGE_XYZ(q,k,j,i)  (0.25*(  q[0][j][i]   + q[0][j][i+1]        \
                                       + q[0][j+1][i] + q[0][j+1][i+1]))

#elif DIMENSIONS == 3

  #define AVERAGE_Y(q,k,j,i)   (0.5*(q[k][j][i] + q[k][j+1][i]))
  #define AVERAGE_Z(q,k,j,i)   (0.5*(q[k][j][i] + q[k+1][j][i]))

  #define AVERAGE_XY(q,k,j,i)  (0.25*(  q[k][j][i]   + q[k][j][i+1] \
                                       + q[k][j+1][i] + q[k][j+1][i+1]) )
  #define AVERAGE_XZ(q,k,j,i)  (0.25*(  q[k][j][i]   + q[k][j][i+1] \
                                       + q[k+1][j][i] + q[k+1][j][i+1]) )
  #define AVERAGE_YZ(q,k,j,i)  (0.25*(  q[k][j][i]   + q[k][j+1][i] \
                                       + q[k+1][j][i] + q[k+1][j+1][i]) )
  #define AVERAGE_XYZ(q,k,j,i) (0.125*(  q[k][j][i]   + q[k][j][i+1]        \
                                       + q[k][j+1][i] + q[k][j+1][i+1]      \
                                       + q[k+1][j][i]   + q[k+1][j][i+1]    \
                                       + q[k+1][j+1][i] + q[k+1][j+1][i+1]))
#endif
/**@} */

