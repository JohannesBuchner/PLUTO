/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief FARGO-MHD module header file.

  Contains contains basic definitions and declarations used by the
  FARGO-MHD module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
  \date    Sept 30, 2016
  \todo   
*/
/* ///////////////////////////////////////////////////////////////////// */

/**<  Set the order of interpolation during the linear transport step.
      Either 2 or 3. Default is 3. */
#ifndef FARGO_ORDER
 #define FARGO_ORDER    3  
#endif
        
/*! Set how often (in number of steps) the total azimuthal 
    velocity should be averaged.                           */
#ifndef FARGO_NSTEP_AVERAGE
 #define FARGO_NSTEP_AVERAGE  10   /* Default is 10 */
#endif

/*! Average background velocity is computed by averaging the orbital 
    velocity (YES) or by prescribing the velocity analytically exactly in 
    the user-supplied function FARGO_SetVelocity  */
#ifndef FARGO_AVERAGE_VELOCITY
 #ifndef SHEARINGBOX
  #define FARGO_AVERAGE_VELOCITY  YES   /* Default is YES */
 #else                             
  #define FARGO_AVERAGE_VELOCITY  NO    /* Default is NO */
 #endif
#endif
                                 
/*! Used to write total/residual velocity */
#ifndef FARGO_OUTPUT_VTOT
  #define FARGO_OUTPUT_VTOT    YES
#endif

/* ----------------------------------------------------------------
    Note: when both FARGO and SHEARINGBOX modules are used, *ALL* 
          source terms are accounted for by the body_force function 
   ---------------------------------------------------------------- */

/* -- SBEG and SEND are the initial and final 
      indices in the direction of the orbital velocity -- */

#if GEOMETRY != SPHERICAL
 #define SDIR   JDIR  /* Orbital direction. */
 #define SBEG   JBEG  /* Starting index in the orbital dir. */
 #define SEND   JEND  /* Final index in the orbital dir. */
 #define NS     NX2
 #define NS_TOT NX2_TOT
#else
 #define SDIR   KDIR
 #define SBEG   KBEG
 #define SEND   KEND
 #define NS     NX3
 #define NS_TOT NX3_TOT
#endif

/* ------------------------------------------------------------------
             Handy macros for dimensional loops
   ------------------------------------------------------------------ */

#define SDOM_LOOP(s) for ((s) = SBEG; (s) <= SEND; (s)++)
  
#if GEOMETRY == SPHERICAL
 #define FARGO_ARRAY_INDEX(A,s,k,j,i)  A[s][j][i]
/* #define FARGO_SELECT(w,k,j,i)  ((w)[j][i])  */
#else
 #define FARGO_ARRAY_INDEX(A,s,k,j,i)  A[k][s][i]
/* #define FARGO_SELECT(w,k,j,i)  ((w)[k][i])  */
#endif 

/* ----------------------------------------------------------------- 
     Prevent compilation when Fargo and AMR are given simultaneously
   ----------------------------------------------------------------- */

#ifdef CH_SPACEDIM
 #error FARGO and AMR are not compatible
#endif
 
void     FARGO_AddVelocity(const Data *, Grid *);
void     FARGO_CHECK (Data_Arr V, Data_Arr U);
void     FARGO_ComputeVelocity(const Data *, Grid *);
double **FARGO_GetVelocity(void);
void     FARGO_SubtractVelocity(const Data *, Grid *);
void     FARGO_ShiftSolution(Data_Arr, Data_Arr, Grid *);
#ifdef PARTICLES
void     FARGO_ShiftParticles(Data *, Grid *, double);
#endif
double   FARGO_SetVelocity(double, double);
void     FARGO_Source(Data_Arr, double, Grid *);
int      FARGO_TotalVelocityIsSet(void);

