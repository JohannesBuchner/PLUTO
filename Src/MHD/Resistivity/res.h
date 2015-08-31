/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Resistive MHD module header file.

  Contains prototypes for the resistive MHD module.
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date    August 27, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping
   *********************************************************** */

void ResistiveFlux (Data_Arr, Data_Arr, double **, double **, int, int, Grid *);
void GetCurrent (const Data *, int, Grid *);
void Resistive_eta (double *, double, double, double, double *, double *);
#ifdef STAGGERED_MHD
 void ComputeStaggeredCurrent (const Data *, Grid *);
 void ComputeStaggeredEta (const Data *, Grid *);
 Data_Arr GetStaggeredEta();
#endif

/* \endcond */