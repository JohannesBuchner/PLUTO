/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Math tools header file.

  This module provides a number of standard numerical routines 
  to achieve simple basic tasks such as

  - Numerical quadrature: GaussQuadrature();
  - LU decomposition functions: LUDecompose(), LUBackSubst();
  - Scalar and multi-dim. root finder methods: Brent(), Ridder(), Brodyen();
  - Vector transformation, VectorCartesianComponents();
  - Ordinary differential equation solvers: ODESolve() with RK4Solve()
    or CK45Solve();
  - Table handling: LocateIndex(), Table2DInterpolate(), GenerateTable2D(),
    InverseLookupTable2D().         
  

  \authors A. Mignone (mignone@ph.unito.it)
  \date    Sep 22, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */

/* ---- ODE solver labels ---- */
        
#define ODE_RK2   2
#define ODE_RK4   4
#define ODE_CK45  45

#define MAXITS         200
#define EPS_FD_JAC     1.0e-7
#define TOLF           EPS_FD_JAC
#define TOLX           EPS_FD_JAC
#define STPMX          100.0
#define TOLMIN         1.0e-6
#define ALF            1.0e-4
#define MAX_ROOT_EQNS  8

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES 
    Function prototyping
   *********************************************************** */

/* -- Functions containd in math_lu_decomp.c -- */

int  LUDecompose (double **, int, int *, double *);
void LUBackSubst (double **, int, int *, double b[]);

/* -- Functions containd in math_misc.c -- */

double BesselJ0(double);
double BesselJ1(double);
double BesselIO(double);
double BesselI1(double);
double BesselKO(double);
double BesselK1(double);
double BesselKn(int, double);
double RandomNumber (double, double);
void   VectorCartesianComponents(double *, double, double, double);

/* -- Functions containd in math_ode.c -- */

void ODE_Solve(double *, int, double, double, double, 
               void (*rhs)(double, double *, double *), int method);

/* -- Functions contained in math_qr_decomp.c -- */

void QRUpdate(double **, double **, int , double *, double *);
void RSolve(double **, int , double *, double *);
void QRSolve(double **, int , double *, double *, double *);
void QRDecompose (double **, int , double *, double *, int *);
void rotate(double **, double **, int , int ,  double , double );

/* -- Functions containd in math_quadrature.c -- */

double GaussQuadrature(double (*func)(double), double, double, int, int);

/* -- Functions containd in math_root_finders.c -- */

int Brent(double (*func)(double, void *), void *, double, double, 
          double, double, double *);
          
int Ridder(double (*func)(double, void *), void *, 
           double, double, double, double, double *);

void FDJacobian(int , double *, double *, double **, 
		void (*vecfunc)(int , double *, double *));
void LineSearch (int , double *, double , double *, 
		double *, double *, double *, double *, double , 
		 int *, void (*vecfunc)(int , double *, double *));
void Broyden(double *, int , int *, void (*vecfunc)(int, double *, double *));


/* -- Functions containd in math_table2D.c -- */

void InitializeTable2D (Table2D *, double, double, int, double, double, int);
void FinalizeTable2D   (Table2D *);
int  Table2DInterpolate   (Table2D *, double, double, double *);
int  InverseLookupTable2D (Table2D *, double, double, double *);
void WriteBinaryTable2D (char *, Table2D *);

/* -- Functions containd in math_interp.c -- */

void MonotoneSplineCoeffs (double *x, double *y, double *dydx, int n,
                           double *a, double *b, double *c, double *d);

void SplineCoeffs (double *x, double *f, double dfL, double dfR, int n,
                   double *a, double *b, double *c, double *d);

void TridiagonalSolve(double *, double *, double *, double *,  double *, int);

/* \endcond */
