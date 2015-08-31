/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO header file for structure declarations.

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 31, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */

typedef struct CMD_LINE {
  int restart;       
  int h5restart;       
  int nrestart;      
  int makegrid;      
  int write;         
  int maxsteps;
  int jet;  /* -- follow jet evolution in a given direction -- */
  int parallel_dim[3];
  int nproc[3];  /* -- user supplied number of processors -- */
  int show_dec; /* -- show domain decomposition ? -- */
  int xres; /* -- change the resolution via command line -- */
  int fill; /* useless, it makes the struct a power of 2 */ 
} Cmd_Line;
   
/* ********************************************************************* */
/*! The Data structure contains the main solution 3D arrays used by 
    the code. 
   ********************************************************************* */
typedef struct DATA{
  double ****Vc;  /**< The main four-index data array used for cell-centered
                       primitive variables. The index order is
                       <tt>Vc[nv][k][j][i]</tt> where \c nv gives the variable
                       index while \c k,\c j and \c i are the
                       locations of the cell in the \f$x_3\f$,
                       \f$x_2\f$ and \f$x_1\f$ direction. */
  double ****Uc;  /**< The main four-index data array used for cell-centered
                       conservative variables. The index order is
                       <tt>Uc[k][j][i][nv]</tt> (\c nv fast running index)
                       where \c nv gives the variable index, \c k,\c j and \c i
                       are the locations of the cell in the \f$x_3\f$,
                       \f$x_2\f$ and \f$x_1\f$ direction. */
  double ****Vs;  /**< The main four-index data array used for face-centered
                       staggered magnetic fields. 
                       The index order is <tt>Vc[nv][k][j][i]</tt>,
                       where \c nv gives the variable index, \c k,\c j and \c i
                       are the locations of the cell in the \f$x_3\f$,
                       \f$x_2\f$ and \f$x_1\f$ direction. */
  double ****Vuser; /**< Array storing user-defined supplementary variables 
                         written to disk. */ 
  double ***Ax1;  /**< Vector potential component in the \f$x_1\f$ direction.*/
  double ***Ax2;  /**< Vector potential component in the \f$x_2\f$ direction.*/
  double ***Ax3;  /**< Vector potential component in the \f$x_3\f$ direction.*/
  double ****J;   /**< Electric current defined as curl(B). */
  unsigned char ***flag; /**< Pointer to a 3D array setting useful integration
                              flags that are retrieved during integration. */
  char fill[28];  /* make the structure a power of two.  */
} Data;

/* ********************************************************************* */
/*! The PLUTO Grid structure contains information pertaining to the 
    computational mesh in a specific 1D coordinate direction. 
    Since PLUTO assumes a logically rectangular system of coordinates, 
    the whole computational domain is obtained as the cartesian product
    of 2 or 3 grid structures.\n

    In parallel, each processor owns a different portion of the domain
    and the grid structures will be different.
    For this reason, in the following member description, we use the 
    word "global" or "local" to refer the the whole computational domain
    or to the sub-domain owned by a single processor.
 
    Similarly, variables ending with a "glob" suffix are intended to be 
    global, i.e., they refer to the whole computational stencil and 
    not to the local processor sub-domain.
   ********************************************************************* */

typedef struct GRID{
  double xi, xf;        /**< Leftmost and rightmost point in the local domain. */
  double *x, *x_glob;   /**< Cell geometrical central points. */
  double *xr, *xr_glob; /**< Cell right interface. */
  double *xl, *xl_glob; /**< Cell left interface. */
  double *dx, *dx_glob; /**< Cell size.  */ 
  double *xgc;          /**< Cell volumetric centroid 
                             (!= x when geometry != CARTESIAN).  */
  double *dV;           /**< Cell volume.  */
  double *A;            /**< Right interface area, A[i] = \f$A_{i+\HALF}\f$. */
  double *r_1;          /**< Geometrical factor 1/r.  */
  double *ct;           /**< Geometrical factor cot(theta).  */
  double *inv_dx;       /**<      */
  double *inv_dxi;      /**< inverse of the distance between the center of 
                             two cells, inv_dxi = \f$\DS \frac{2}{\Delta x_i +
                             \Delta x_{i+1}}\f$.     */
  double dl_min;      /**<  minimum cell length (e.g. min[dr, r*dth,
                            r*sin(th)*dphi] (GLOBAL DOMAIN).  */
  int np_tot_glob; /**< Total number of points in the global domain 
                        (boundaries included). */
  int np_int_glob; /**< Total number of points in the global domain 
                        (boundaries excluded). */
  int np_tot;      /**< Total number of points in the local domain 
                        (boundaries included). */
  int np_int;      /**< Total number of points in the local domain 
                        (boundaries excluded). */
  int nghost;      /**< Number of ghost zones. */
  int lbound;      /**< When different from zero, it specifies the boundary
                        condition to be applied at leftmost grid side where  
                        the physical boundary is located.
                        Otherwise, it equals zero if the current 
                        processor does not touch the leftmost physical boundary. 
                        This evantuality (lbound = 0) is possible only
                        in PARALLEL mode.  */ 
  int rbound;      /**< Same as lbound, but for the right edge of the grid. */
  int gbeg;        /**< Global start index for the global array. */
  int gend;        /**< Global end   index for the global array. */
  int beg;         /**< Global start index for the local array. */
  int end;         /**< Global end   index for the local array. */
  int lbeg;        /**< Local start  index for the local array. */
  int lend;        /**< Local end    index for the local array. */
  int uniform;     /* = 1 when the grid is cartesian AND uniform everywhere  */
  int nproc;       /**< number of processors for this grid. */
  int rank_coord;  /**< Parallel coordinate in a Cartesian topology. */
  int level;       /**< The current refinement level (chombo only). */
  char fill[40];   /* useless, just to make the structure size a power of 2 */
} Grid;
   
/* ********************************************************************* */
/*! This structure contains one-dimensional vectors of conserved variables,   
    primitive variables, fluxes and so on, used during the 
    reconstruct-Solve-Average strategy. 
    It is a frequently passed to the Riemann solver routines, source and 
    flux functions, etc.
   ********************************************************************* */
typedef struct STATE_1D{
  double **v;    /**< Cell-centered primitive varables at the base time level,
                      v[i] = \f$ \vec{V}^n_i \f$ . */
  double **vL;   /**< Primitive variables to the left of the interface, 
                       \f${\rm vL[i]} \equiv \vec{V}_{i,+} = 
                           \vec{V}^L_{i+\HALF} \f$. */
  double **vR;   /**< Primitive variables to the right of the interface, 
                       \f$\mathrm{vR[i]} \equiv \vec{V}^R_{i+\HALF} \f$. */
  double **vm;   /**< prim vars at i-1/2 edge, vm[i] = vR(i-1/2)     */
  double **vp;   /**< prim vars at i+1/2 edge, vp[i] = vL(i+1/2)     */

  double **uL;   /**< same as vL, in conservative vars */
  double **uR;   /**< same as vR, in conservative vars */
  double **um;   /**< same as vm, in conservative vars */
  double **up;   /**< same as vp, in conservative vars */

  double **flux;      /**< upwind flux computed with the Riemann solver */
  double **visc_flux; /**< Viscosity flux             */
  double **visc_src;  /**< Viscosity source term      */
  double **tc_flux;   /**< Thermal conduction flux    */
  double **res_flux;  /**< Resistive flux (current)   */
                          
  double ***Lp, ***Rp; /**< Left and right primitive eigenvectors */
  double **lambda;     /**< Characteristic speed associated to Lp and Rp */
  double *lmax;   /**< Define the maximum k-characteristic speed over the domain */
  double *a2;     /**< Sound speed squared */ 
  double *h;      /**< Enthalpy. */
  double **src;     

  double **vh;      /**< Primitive    state at n+1/2 (only for one step method) */
  double **rhs;     /**< Conservative right hand side */
  double *press;    /**< Upwind pressure term computed with the Riemann solver */
  double *bn;       /**< Face magentic field, bn = bx(i+1/2) */
  double *SL;       /**< Leftmost  velocity in the Riemann fan at i+1/2 */
  double *SR;       /**< Rightmost velocity in the Riemann fan at i+1/2 */
  unsigned char *flag;
  double fill1, fill2;
} State_1D;

typedef struct TABLE2D {
  char **defined;
  int nx;  /**< Number of columns or points in the x direction */
  int ny;  /**< Number of rows    or points in the y direction */  
  int nf;
  int interpolation;   /**< LINEAR/SPLINE1  */   
  int **i; 
  int id;
  double *x;  /**< array of x-values (not uniform) */
  double *y;  /**< array of y-values (not uniform) */
  double *dx; /**< grid spacing array in the x direction (not uniform) */
  double *dy; /**< grid spacing array in the y direction (not uniform) */
  double *lnx; /**< array of log10(x) values (uniform) */
  double *lny; /**< array of log10(y) values (uniform) */
  double **f;

  double **a;  /**< Spline coefficient (x^3) */ 
  double **b;  /**< Spline coefficient (x^2) */
  double **c;  /**< Spline coefficient (x)   */
  double **d;  /**< Spline coefficiten (1)   */

  double **dfx;
  double **dfy;
  double *fmin;
  double *fmax;
  double *df;
  double lnxmin; /**< lower limit (in log10) in the x-direction */
  double lnxmax; /**< upper limit (in log10) in the x-direction */
  double lnymin; /**< lower limit (in log10) in the y-direction */
  double lnymax; /**< upper limit (in log10) in the y-direction */
  double dlnx;   /**< uniform spacing in log10(x) */
  double dlny;   /**< uniform spacing in log10(y) */
  double dlnx_1;
  double dlny_1;
} Table2D;


/* ********************************************************************* */
/*! The Time_Step structure contains essential information for 
    determining the time step.
   ********************************************************************* */
typedef struct TIME_STEP{
  double *cmax;     /**< Maximum signal velocity for hyperbolic eqns. */
  double inv_dta;   /**< Inverse of advection (hyperbolic) time step, 
                         \f$ \lambda/\Delta l\f$.*/
  double inv_dtp;   /**< Inverse of diffusion (parabolic)  time step 
                         \f$ \eta/\Delta l^2\f$. */
  double dt_cool;   /**< Cooling time step. */
  double cfl;       /**< Courant number for advection. */
  double cfl_par;   /**< Courant number for diffusion (STS only). */
  double rmax_par;  
  int    Nsts;      /**< Maximum number of substeps used in STS. */
  int    Nrkc;      /**< Maximum number of substeps used in RKC. */
  char  fill[24];   /* useless, just to make the structure size a power of 2 */
} Time_Step;


/* ********************************************************************* */
/*! The Output structure contains essential information for I/O.
   ********************************************************************* */
typedef struct OUTPUT{
  int    type;            /**< output format (DBL, FLT, ...) - one per output */
  int    nvar;            /**< tot. # of vars that can be written - same for all   */
  int    user_outs;
  int    cgs;             /**< when set to 1 saves data in c.g.s units     */
  int    nfile;           /**< current number being saved - one per output */
  int    dn;              /**< step increment between outputs - one per output */
  int    *stag_var;       /**< centered or staggered variable - same for all   */
  int    *dump_var;       /**< select vars being written      - one per output */
  char   mode[32];        /**< single or multiple files       - one per output */
  char   **var_name;      /**< variable names                 - same for all   */
  char   ext[8];          /**< output extension                            */
  char   dir[256];       /**< output directory name                        */
  double dt;           /**< time increment between outputs   - one per output */
  double dclock;       /**< time increment in clock hours     - one per output */
  double ***V[64];     /**< pointer to arrays being written   - same for all  */
  char   fill[168];    /**< useless, just to make the structure size a power of 2 */
} Output;

/* ********************************************************************* */
/*! The Runtime structure contains runtime initialization parameters
    read from pluto.ini (or equivalent). 
   ********************************************************************* */
typedef struct RUNTIME{
  int    npoint[3];           /**< Global number of zones in the interior domain */
  int    left_bound[3];       /**< Array of left boundary types */
  int    right_bound[3];      /**< Array of right boundary types */
  int    grid_is_uniform[3];  /* = 1 when grid is uniform, 0 otherwise */
  int    npatch[5];           /**< The number of grid patches  */
  int    patch_npoint[5][16]; /* number of points per patch */
  int    patch_type[5][16];             
  int    log_freq;            /**< The log frequency (\c log) */
  int    user_var;            /**< The number of additional user-variables being
                                 held in memory and written to disk */
  int    anl_dn;               /*  number of step increment for ANALYSIS */
  char   solv_type[64];         /**< The Riemann solver (\c Solver) */
  char   user_var_name[128][128];
  char   output_dir[256];         /**< The name of the output directory
                                       (\c output_dir for static PLUTO,
                                        \c Output_dir for PLUTO-Chombo)  */
  Output output[MAX_OUTPUT_TYPES];  
  double patch_left_node[5][16];  /*  self-expl. */
  double  cfl;               /**< Hyperbolic cfl number (\c CFL) */
  double  cfl_max_var;       /**< Maximum increment between consecutive time
                                  steps (\c CFL_max_var). */
  double  cfl_par;           /**< (STS) parabolic  cfl number */
  double  rmax_par;          /**< (STS) max ratio between current time
                                step and parabolic time step */
  double  tstop;           /**< The final integration time (\c tstop) */
  double  first_dt;        /**< The initial time step (\c first_dt) */
  double  anl_dt;          /**< Time step increment for Analysis()
                                ( <tt> analysis (double) </tt> )*/
  double  aux[32];         /* we keep aux inside this structure, 
                              since in parallel execution it has
                              to be comunicated to all processors  */
} Runtime;

typedef struct RESTART{
  int nstep;
  int nfile[MAX_OUTPUT_TYPES];
  double t;
  double dt;
} Restart;

typedef struct RGB{
  unsigned char r, g, b;
} RGB;

typedef struct IMAGE{
  int    nrow, ncol;    /* -- image rows and columns -- */
  int    slice_plane;   /* -- one of X12_PLANE, X13_PLANE, X23_PLANE -- */
  int    logscale;      /* -- YES/NO for log scale -- */
  char   *colormap;     /* -- colormap name -- */
  char   basename[32];  /* -- image base name (no extensions) -- */
  unsigned char r[256], g[256], b[256]; /* -- colortable saved here -- */
  RGB    **rgb;         /* -- rgb array containing image values -- */
  double max;           /* -- max image value -- */
  double min;           /* -- min image value -- */
  double slice_coord;   /* -- slice coord orthogonal to slice_plane -- */
} Image;

typedef struct FLOAT_VECT{
  float v1, v2, v3;
} Float_Vect;

typedef struct INDEX{
  int ntot, beg, end;
  int *pt1, t1, t1_beg, t1_end;
  int *pt2, t2, t2_beg, t2_end;
  char fill[20]; /* useless, just to make the structure size a power of 2 */
} Index;

/* ********************************************************************* */
/*! The List defines a collection of integer values typically used
    as argument to the FOR_EACH() macro.
   ********************************************************************* */
typedef struct INT_LIST{
  int indx[2046]; /**< Array of integers containg variables indices. */
  int nvar;       /**< Number of variables. */
  int i;          /**< Internal counter. */
} intList;

/* ********************************************************************* */
/*! The RBox (= Rectangular Box) defines a rectangular portion of the 
    domain in terms of the grid indices <tt>[ib,jb,kb]</tt> corresponding
    to the lower corner and <tt>[ie,je,ke]</tt> corresponding to the
    upper corner. 
    The integer \c vpos specifies the variable location with respect to 
    the grid (e.g. center/staggered). 

    \note The lower and upper grid indices may also be reversed 
          (e.g. box->ib > box->ie). In this case the macro ::BOX_LOOP
          automatically reset the directional increment (box->di) to -1.
   ********************************************************************* */
typedef struct RBOX{
  int ib; /**< Lower corner index in the x1 direction. */
  int ie; /**< Upper corner index in the x1 direction. */
  int jb; /**< Lower corner index in the x2 direction. */
  int je; /**< Upper corner index in the x2 direction. */
  int kb; /**< Lower corner index in the x3 direction. */
  int ke; /**< Upper corner index in the x3 direction. */
  int di; /**< Directional increment (+1 or -1) when looping over the 1st 
               dimension of the box. Automatically set by the ::BOX_LOOP macro. */
  int dj; /**< Directional increment (+1 or -1) when looping over the 2nd 
               dimension of the box. Automatically set by the ::BOX_LOOP macro. */
  int dk; /**< Directional increment (+1 or -1) when looping over the 3rd 
               dimension of the box. Automatically set by the ::BOX_LOOP macro. */
  int vpos; /**< Location of the variable inside the cell. */
} RBox;
