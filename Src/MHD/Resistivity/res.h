/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Resistive MHD module header file.

  Contains prototypes for the resistive MHD module.
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date    June 20, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping
   *********************************************************** */

void ResistiveFlux (Data_Arr, Data_Arr, double **, double **, int, int, Grid *);
void Resistive_eta (double *, double, double, double, double *, double *);
void ResistiveRHS (const Data *, Data_Arr, double **,
                  double **, double, int, int, Grid *);

void ComputeStaggeredCurrent (const Data *, Grid *);
void ComputeStaggeredEta (const Data *, Grid *);
Data_Arr GetStaggeredEta();

/* \endcond */
