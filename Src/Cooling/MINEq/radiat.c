#include "pluto.h"
#include "cooling_defs.h"

/* ----------------------------------------------------------------------- 
    Global variables definition. They are declared inside cooling_defs.h
   ----------------------------------------------------------------------- */

/*                           H      He      C      N      O       Ne      S      Fe  */
const double elem_ab[] = {  0.93, 0.074, 3.0e-4,  5.e-5, 4.0e-4, 7.0e-5, 1.5e-5, 2.69e-5 }; /* Number densities */
double   elem_ab_est[] = {  0.93, 0.074, 3.0e-4,  5.e-5, 4.0e-4, 7.0e-5, 1.5e-5 };  /* Number densities in Orion nebula, Esteban & al 2004 */
double   elem_ab_lod[] = {  0.92,  0.08, 2.3e-4, 6.2e-5, 4.4e-4, 7.0e-5, 1.25e-5 }; /* Number densities, Solar, Lodders 2003 ApJ */
double   elem_ab_sol[] = {  0.93, 0.074, 3.0e-4,  9.e-5,  7.e-4,  7.e-5, 1.e-5 };   /* Number fractions, Solar    */
double   elem_ab_uni[] = {  0.93, 0.072, 5.0e-4,  9.e-5,  8.e-4,  8.e-5, 2.e-5 };   /* Number fractions, Universe */
const double elem_mass[]   = { 1.007, 4.002,  12.01,  14.01,  15.99,  20.18, 32.07, 55.845 };  /*   Atomic mass, in a.m.u.   */

const int elem_part[] = {0, 1, 1 
                   C_EXPAND(2, 2, 2, 2, 2) 
                   N_EXPAND(3, 3, 3, 3, 3) 
                   O_EXPAND(4, 4, 4, 4, 4) 
                  Ne_EXPAND(5, 5, 5, 5, 5) 
                   S_EXPAND(6, 6, 6, 6, 6)
                  Fe_EXPAND(7, 7, 7) };

const double rad_rec_z[] = { 1., 1., 2. 
                       C_EXPAND(1., 2., 3., 4., 5.)
                       N_EXPAND(1., 2., 3., 4., 5.) 
                       O_EXPAND(1., 2., 3., 4., 5.)
                      Ne_EXPAND(1., 2., 3., 4., 5.) 
                       S_EXPAND(1., 2., 3., 4., 5.)
                      Fe_EXPAND(1., 2., 3.) };
const double coll_ion_dE[] = { 13.6, 24.6, 54.4 
                        C_EXPAND(11.3,    24.4,   47.9, 64.5, 392.1)
                        N_EXPAND(14.5,    29.6,   47.5, 77.5,  97.9) 
                        O_EXPAND(13.6,    35.1,   54.9, 77.4, 113.9)
                       Ne_EXPAND(21.6,    41.0,   63.5, 97.1, 126.2) 
                        S_EXPAND(10.4,    23.3,   34.8, 47.3,  72.6)
                       Fe_EXPAND(7.87, 16.1879, 30.652) };     
COOL_COEFF CoolCoeffs;
/* ********************************************************************* */
void Radiat (double *v, double *rhs)
/*! 
 *  Cooling for optically thin plasma up to about 200,000 K
 *  Plasma composition: H, HeI-II, CI-V, NI-V, OI-V, NeI-V, SI-V
 *  Assumed abundances in elem_ab
 *  Uses S : Array = Variables vector x line points
 *       rhs :  output for the system of ODE
 *       ibeg, iend : begin and end points of the current line
 *
 *
 *********************************************************************** */
{
  int  nv, j, k, cooling_nT, cooling_nNe, jj;
  int  ti1, ti2, ne1, ne2, nrt;
  double   mu, sT, scrh, tmpT, tmpNe, tt1, tt2, nn1, nn2, tf1, tf2, nf1, nf2;
  double   N, n_el, rlosst, em, cf1, cf2;
  double   T, *X, *RS;
  static double ***tab;
  static double N_rho;
  static double E_cost, Unit_Time; /* -- for dimensionalization purposes -- */
  double em2, em3;
double prs;
 
/* --------------------------------------------------------------------
            Load tables from disk at first call.
            Tables are 2-D with T and Ne being the coordinates.
   -------------------------------------------------------------------- */

  if (tab == NULL) {  

    E_cost    = UNIT_LENGTH/UNIT_DENSITY/
               (UNIT_VELOCITY*UNIT_VELOCITY*UNIT_VELOCITY);
    Unit_Time = UNIT_LENGTH/UNIT_VELOCITY;
    
    /* ---------------------------------------
           Compute cooling function tables
       --------------------------------------- */

    for (nv = 0; nv < NIONS; nv++) CoolCoeffs.dLIR_dX[nv] = 0.0;
    
    n_el = C_NeMIN;
    ne1 = 0;
    while (n_el < C_NeMAX) {
      T  = C_TMIN;
      ne2 = 0;
      while (T  < C_TMAX)  {
        tmpT = T;
        T  = T*exp(C_TSTEP);   /* should be *exp(0.02)  */
        ne2 = ne2 + 1;
      }
      tmpNe = n_el;
      n_el = n_el*exp(C_NeSTEP);   /* should be *exp(0.06)  */
      ne1 = ne1 + 1;
    }

    tab = ARRAY_3D(NIONS, ne1, ne2, double);
    Create_Losses_Tables(tab, &cooling_nT, &cooling_nNe);
    N_rho = find_N_rho();
  }    /*  end load tables   */

/* ---------------------------------------
    Force species to lie between 0 and 1 
   --------------------------------------- */

  for (nv = NFLX; nv < (NFLX + NIONS); nv++){
    v[nv] = MIN(v[nv], 1.0);
    v[nv] = MAX(v[nv], 0.0);
  }

  prs = v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }

  mu = MeanMolecularWeight(v); 
  T  = prs/v[RHO]*KELVIN*mu;

  if (mu < 0.0){
    print ("! Radiat: negative mu\n");
    QUIT_PLUTO(1);
  }

 /* ---------------------------------
     offset pointers to work more
     efficiently.
    --------------------------------- */

  X  = v + NFLX; 
  RS = rhs + NFLX;

  sT = sqrt(T);
  N  = v[RHO]*N_rho;       /* -- Total number density -- */
    
/* -----------------------------------
     compute electron number density 
   ----------------------------------- */

  CoolCoeffs.dnel_dX[0] = -N*elem_ab[el_H];   
  n_el = N*(1.0 - v[X_HI])*elem_ab[el_H]; /* -- contribution from ionized H --  */

  for (nv = 1; nv < NIONS; nv++) {
    CoolCoeffs.dnel_dX[nv] = N*(rad_rec_z[nv] - 1.0)*elem_ab[elem_part[nv]];
    n_el                  += X[nv]*CoolCoeffs.dnel_dX[nv];    
  }

  if (n_el/N < 1.e-4) n_el = N*1.e-4;  /*  OK ????? */
   
  CoolCoeffs.Ne = n_el;  
  
  Find_Rates(T, n_el, N, v);  /*  find transition rates  */

/* -----------------------------------------------------------------------------
          Compute RHS of the system of ODE.
          Evaluate the transition rates  
         (multiply the rates given by find_rates()
          by numerical densities)
   ----------------------------------------------------------------------------- */

   /* ----------------------------------------------
       Compute right hand side for all ions 
      ---------------------------------------------- */

  for (nv = 1; nv <= NIONS - 2; nv++) {  
    RS[nv] =   X[nv - 1]*CoolCoeffs.Lrate[nv] 
             - X[nv]    *CoolCoeffs.Crate[nv] 
             + X[nv + 1]*CoolCoeffs.Rrate[nv];
  }
  RS[0] = (1.0 - X[0])*CoolCoeffs.Rrate[0] - X[0]*CoolCoeffs.Crate[0];   /* separate for H */

  nv = NIONS - 1;
  RS[nv] = X[nv - 1]*CoolCoeffs.Lrate[nv] - X[nv]*CoolCoeffs.Crate[nv];

  for (nv = 0; nv < NIONS; nv++) RS[nv] *= Unit_Time;

/* --------------------------------------------------------------- 
                    check sums 
   --------------------------------------------------------------- */
/*
{ 
double sum;

sum=fabs(RS[1] + RS[2]);
if (sum > 1.e-17){
  print("> Sum(RHS) != 0 for He!!\n");
  exit(1);
}

sum = 0.0;
for (nv = 0; nv < C_IONS; nv++) sum += RS[CI-NFLX+nv];
if (fabs(sum) > 1.e-10){
  print("> Sum(RHS) != 0 for C !!  (%12.6e)\n", sum);
  exit(1);
}
}

nv = 8;
sum=fabs(RS[nv] + RS[nv+1] + RS[nv+2] + RS[nv+3] + RS[nv+4]);
if ( sum > 1.e-10){
 print("> Sum(RHS) != 0 for N!!\n");
  exit(1);
}
nv = 13;
sum=fabs(RS[nv] + RS[nv+1] + RS[nv+2] + RS[nv+3] + RS[nv+4]);
if ( sum > 1.e-10){
 print("> Sum(RHS) != 0 for O!!\n");
  exit(1);
}
nv = 18;
sum = fabs(RS[nv] + RS[nv+1] + RS[nv+2] + RS[nv+3] + RS[nv+4]);
if ( sum > 1.e-10){
 print("> Sum(RHS) != 0 for Ne!!\n");
  exit(1);
}
nv = 23;
fabs(RS[nv] + RS[nv+1] + RS[nv+2] + RS[nv+3] + RS[nv+4]);
if ( sum > 1.e-10){
 print("> Sum(RHS) != 0 for S!!\n");
  exit(1);
}
}
*/

/*
print("%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", v[HI], 
        RS[0],RS[1], RS[2], RS[3], RS[4]);
} 
exit(1);
*/

/* ********************************************************************************************************* */


    /* -------------------------------------------------------------
         Find the indexes in the tables needed for interpolation  

        Constrain both T and n_el to lie in the range provided 
        by the tables.
        This does not affect the original vector of primitive
        quantities, but it is more computationally efficient.
        In other words, when T > Tmax, the cooling function 
        will saturate at Tmax.
       ------------------------------------------------------------- */

  tmpT  = T;   /*  keep the double value for after the interpolation   */
  tmpNe = n_el;

  if (tmpT <= C_TMIN) tmpT = C_TMIN + C_TSTEP*0.5;
  if (tmpT >= C_TMAX) tmpT = C_TMAX - C_TSTEP*0.5;

  if (n_el <= C_NeMIN) tmpNe = C_NeMIN + C_NeSTEP*0.5;
  if (n_el >= C_NeMAX) tmpNe = C_NeMAX - C_NeSTEP*0.5;
 
  scrh = log(tmpT/C_TMIN)/C_TSTEP;
  ti1 = floor(scrh);
  ti2 = ceil (scrh);

  scrh = log(tmpNe/C_NeMIN)/C_NeSTEP;
  ne1  = floor(scrh);
  ne2  = ceil (scrh);
    
  tt1 = C_TMIN*exp(ti1*C_TSTEP);
  tt2 = C_TMIN*exp(ti2*C_TSTEP);

  nn1 = C_NeMIN*exp(ne1*C_NeSTEP);
  nn2 = C_NeMIN*exp(ne2*C_NeSTEP);
    
  tf1 = (tt2 - tmpT)/(tt2 - tt1);
  tf2 = (tmpT - tt1)/(tt2 - tt1);
  nf1 = (nn2 - tmpNe)/(nn2 - nn1);
  nf2 = (tmpNe - nn1)/(nn2 - nn1);

/* ------------------------------------------------------------
                 get \Delta e[k] coefficients 
   ------------------------------------------------------------ */

  for (nv = 0; nv < NIONS; nv++) {
    cf1 = tf1*tab[nv][ne1][ti1] + tf2*tab[nv][ne1][ti2];
    cf2 = tf1*tab[nv][ne2][ti1] + tf2*tab[nv][ne2][ti2];
    CoolCoeffs.de[nv] = nf1*cf1 + nf2*cf2;
  }

/* -- Save the slopes of the radiative losses coefficients from the table --  */

  scrh = 0.5/(nn2 - nn1);
  for (nv = 0; nv < NIONS; nv++) {
    CoolCoeffs.de_dne[nv] = (  tab[nv][ne2][ti1] + tab[nv][ne2][ti2]
                             - tab[nv][ne1][ti1] - tab[nv][ne1][ti2])*scrh;
  }

/* ---------------------------------------------------------------
     For HeII we impose a temperature cut-off above 9.e4 K
     since HeII --> HeIII
   --------------------------------------------------------------- */

  if (T >= 9.e4) {
    scrh = exp(-(T - 9.e+4)/4.e+4);
    CoolCoeffs.de[2]     *= scrh;
    CoolCoeffs.de_dne[2] *= scrh;
  }

  /*  ... and then continue with the radiative losses computation  */

  em = 0.0;
  for (nv = 0; nv < NIONS; nv++) {
    em += X[nv]*CoolCoeffs.de[nv]*elem_ab[elem_part[nv]];
  }
  
  if  (em < 0.0) print("> 2 negative em => possitive losses???? em = %12.6e\n",em);
     
/* -------------------------------------------------------------------
     Now finally add the ionization / recombination losses ...
   ------------------------------------------------------------------- */

 /* -----------------------------------
     thermal energy lost by ionization 
    ----------------------------------- */

  CoolCoeffs.dLIR_dX[0] = 1.08e-8*sT*elem_ab[el_H]/13.6*exp(-157890.0/T)*CONST_eV; 
  em += CoolCoeffs.dLIR_dX[0]*v[X_HI]; 

 /*
  if ( (1.08e-8*sT*v[HI]*elem_ab[el_H]/13.6*exp(-157890.0/T) *CONST_eV) / (v[HI] *CoolCoeffs.Crate[0]*elem_ab[el_H]*coll_ion_dE[0]*Unit_Time*CONST_eV) > 2.0) {
    print("too different em's: %12.6e  -   %12.6e\n",
              (1.08e-8*sT*v[HI]*elem_ab[el_H]/13.6*exp(-157890.0/T) *CONST_eV), (v[HI] *CoolCoeffs.Crate[0]*elem_ab[el_H]*coll_ion_dE[0]*Unit_Time*CONST_eV));
    
  }
  if ( (1.08e-8*sT*v[HI]*elem_ab[el_H]/13.6*exp(-157890.0/T) *CONST_eV) / (v[HI] *CoolCoeffs.Crate[0]*elem_ab[el_H]*coll_ion_dE[0]*Unit_Time*CONST_eV) < 0.5) {
    print("too different em's: %12.6e  -   %12.6e\n",
              (1.08e-8*sT*v[HI]*elem_ab[el_H]/13.6*exp(-157890.0/T) *CONST_eV), (v[HI] *CoolCoeffs.Crate[0]*elem_ab[el_H]*coll_ion_dE[0]*Unit_Time*CONST_eV));
    
  }
*/ 
/*  
  em2 = v[HI]*CoolCoeffs.Crate[0]*elem_ab[el_H]*coll_ion_dE[0]*CONST_eV*N;
  for (nv=2; nv<NIONS; nv++) {
    if (CoolCoeffs.Lrate[nv]>0) em2 += X[nv-1]*CoolCoeffs.Lrate[nv]*elem_ab[elem_part[nv]]*coll_ion_dE[nv]*CONST_eV*N;
  }
*/
      /* --------------------------------------
          thermal energy lost by recombination
          (2/3 kT per recombination).  
         -------------------------------------- */

  scrh  = 2.6e-11/sT*elem_ab[el_H]*T/11590.0*0.67*CONST_eV;
  em   += scrh*(1.0 - v[X_HI]); 
  CoolCoeffs.dLIR_dX[0] -= scrh;

/*
  if (T>4000) print("ratesHI: old %12.6e    new %12.6e    at T = %12.6e\n",2.6e-11/sT*(1.0 - v[HI])*elem_ab[el_H]*T/11590.0*0.67*CONST_eV,
                          (1.0 - v[HI])*CoolCoeffs.Rrate[0]*elem_ab[el_H]*0.67*CONST_kB*T,T);
  */
/*
  em3 = (1.0 - v[HI])*CoolCoeffs.Rrate[0]*elem_ab[el_H]*0.67*CONST_kB*T*N; 
  for (nv=1; nv<NIONS-1; nv++) {
    if (CoolCoeffs.Rrate[nv]>0) em3 += X[nv+1]*CoolCoeffs.Rrate[nv]*elem_ab[elem_part[nv]]*0.67*CONST_kB*T*N;
  }
*/
    
  /* --------------------------------------------------------
          Add cooling contribution from bremsstrahlung 
     -------------------------------------------------------- */
    
  scrh  = 1.42e-27*sT*elem_ab[el_H];
  CoolCoeffs.dLIR_dX[2]  = 1.42e-27*sT*elem_ab[el_He];
  em   += scrh*(1.0 - v[X_HI]) + CoolCoeffs.dLIR_dX[2]*v[X_HeII];

  CoolCoeffs.dLIR_dX[0] -= scrh;

  /* --------------------------------------------------------------------
          Add cooling contribution from FeII  1.6 & 25 micron
          and Mg II 2800 and Si II 35 micron
          emission lines. This is TEMPORARY until the adding 
          of FeI and FeII as ions to NEq !
     -------------------------------------------------------------------- */
  
  em +=  1.6e-12*8.63e-6*0.3*0.0495*exp(-575.0/T)/sT * 580.0*sT/(n_el + 580.0*sT) *0.00004;
  em +=  1.6e-12*8.63e-6*0.39*0.775*exp(-8980.0/T)/sT * 1130.0*sT/(n_el + 1130.0*sT) *0.00004;

  em += 1.6e-12*8.63e-6*8.0*4.43*exp(-5.13e4/T)/sT*1.e10*sT/(n_el + 1.e10*sT)*0.00002*0.7; /* Mg II 2800: Mendoza    */
  em += 1.6e-12*8.63e-6*2.85*0.0354*exp(-410.0/T)/sT*16.8*sT/(n_el + 16.8*sT)*0.000033*0.7; /* Si II 35 micron: Dufton&Kingston, MNRAS 248 */

  /* ---------------------------------------------------
      rlosst is the energy loss in units of erg/cm^3/s;
      it must be multiplied by cost_E in order to match 
      non-dimensional units.
      Source term for the neutral fraction scales with 
      Unit_Time

       Pressure source term equals zero when 
       T < g_minCoolingTemp 
     --------------------------------------------------- */

  rlosst  = em*n_el*N;
/*
  rlosst += em3 + em2;
*/

/*
  if (T > g_minCoolingTemp) rhs[PRS] = -E_cost*rlosst*(g_gamma - 1.0);
  else                rhs[PRS] = 0.0;
*/

  rhs[RHOE] = -E_cost*rlosst;
  rhs[RHOE] *= 1.0/(1.0 + exp( -(T - g_minCoolingTemp)/100.)); /* -- cut cooling function --*/

  if (rhs[RHOE] > 0.0) {
     print ("! Error: positive radiative losses!  %12.6e  %12.6e  %12.6e \n", em, n_el, N);
     QUIT_PLUTO(1);
  }  
}

/* ********************************************************************* */
void Find_Rates(double T, double Ne, double N, double *v)
/*
 *
 * PURPOSE
 *
 *  Transition rates computation routine   
 *
 *   Output: CoolCoeffs.Crate, CoolCoeffs.Rrate, CoolCoeffs.Lrate
 *
 *   dN = (   CoolCoeffs.Lrate * n_(ion-1) 
 *          + CoolCoeffs.Rrate  * n_(ion+1) - CoolCoeffs.Crate * n_ion ) * dt
 *
 *
 *********************************************************************** */
{
  double dn, lam, t4, tmprec, ft1, ft2, tmpT, scrh;
  int ti1, ti2, cindex, i, ions, nv, tindex, j, k;
  static double ***ion_data, **intData;

  if (ion_data == NULL) {  /* -- compute ionization rates tables -- */
    ion_data = ARRAY_3D(I_g_stepNumber, 7, NIONS, double);
    intData  = ARRAY_2D(7, NIONS, double);
    Create_Ion_Coeff_Tables(ion_data);   
  }  

  for (nv = 0; nv < NIONS; nv++ ) {
    CoolCoeffs.Rrate[nv] = 0.0;
    CoolCoeffs.Lrate[nv] = 0.0;
    CoolCoeffs.Crate[nv] = 0.0;
    CoolCoeffs.Ra[nv] = 0.0;
    CoolCoeffs.La[nv] = 0.0;
    CoolCoeffs.Ca[nv] = 0.0;
    CoolCoeffs.Rb[nv] = 0.0;
    CoolCoeffs.Lb[nv] = 0.0;
    CoolCoeffs.Cb[nv] = 0.0;
    CoolCoeffs.Rc[nv] = 0.0;
    CoolCoeffs.Lc[nv] = 0.0;
    CoolCoeffs.Cc[nv] = 0.0;
  }

  tmpT = T;  
  if      (tmpT > I_TEND) tmpT = I_TEND - I_TSTEP*0.5;
  else if (tmpT < I_TBEG) tmpT = I_TBEG + I_TSTEP*0.5;
  
  tindex = floor( (tmpT - I_TBEG)/I_TSTEP);

  ft1 = ( I_TBEG + I_TSTEP*((double)tindex + 1.0) - tmpT ) / I_TSTEP;
  ft2 = ( tmpT - I_TBEG - I_TSTEP*((double)tindex) )       / I_TSTEP;

  for (k = 0; k < 7; k++) {
  for (ions = 0; ions < NIONS; ions++) {
    intData[k][ions] = ion_data[tindex][k][ions]*ft1 + ion_data[tindex + 1][k][ions]*ft2;
  }}

  /* -------------------------------------------------------------
                  Computing rates for H          
     ------------------------------------------------------------- */

  lam = 157890. / T;
  CoolCoeffs.Crate[0] += Ne*intData[0][el_H]; /*  collisional ionization */
  CoolCoeffs.Rrate[0] += Ne*intData[1][el_H]; /*  total recombination  -  NIFS-DATA-54 - Kato & Asano 1999 , from Aldrovandi & Pequignot 1973 */

  CoolCoeffs.fCH = intData[0][el_H];
  CoolCoeffs.fRH = intData[1][el_H];

  /* contribution from charge - exchange (recombination and ionization of all ions) :  */
/*
  for (ions=3; ions<28; ions++) {  

        if (T/10000.<=chtrH_ion_upT[ions]) CoolCoeffs.Crate[0] += X[ions+1][ii] * N * elem_ab[elem_part[ions+1]] * 1.e-9 * chtrH_rec_a[ions] * pow( (T/10000.), chtrH_rec_b[ions] ) * ( 1. +  chtrH_rec_c[ions] * exp ( chtrH_rec_d[ions]*(T/10000.) ));
        if (T<=chtrH_rec_upT[ions]) CoolCoeffs.Rrate[0] +=  1.e-9 * chtrH_ion_a[ions] * pow( (T/10000.), chtrH_ion_b[ions] ) * ( 1. +  chtrH_ion_c[ions] * exp ( chtrH_ion_d[ions]*(T/10000.) ));
    
      tmprec = v[ions] * elem_ab[elem_part[ions]] / elem_ab[el_H] / fabs(v[HI]-1.e-13) * N * ion_data[4][ions];
      CoolCoeffs.Crate[0] += tmprec;
      CoolCoeffs.Ca[0]    += ion_data[4][ions];
  } 
*/

/* -------------------------------------------------------------------
                Do the other ions now   
   ------------------------------------------------------------------- */

  for (ions = 2; ions < NIONS - 1; ions++) {

  /*  ------------------------------------------
       collisional ionization  -  Voronov 1997  
      ------------------------------------------ */
    
    tmprec = Ne*intData[0][ions];
    CoolCoeffs.Crate[ions]   += tmprec;
    CoolCoeffs.Lrate[ions+1] += tmprec;
    CoolCoeffs.Ca[ions]   += intData[0][ions];
    CoolCoeffs.La[ions+1] += intData[0][ions];

  /* -------------------------------------------------------
      radiative recombination  -  Pequignot & al, 1991  A&A 
     ------------------------------------------------------- */
    
    tmprec = Ne*intData[1][ions-1];
    CoolCoeffs.Crate[ions]   += tmprec;
    CoolCoeffs.Rrate[ions-1] += tmprec;
    CoolCoeffs.Ca[ions]      += intData[1][ions-1];
    CoolCoeffs.Ra[ions-1]    += intData[1][ions-1];

  /*  -----------------------------------------------------------
        dielectronic recombination  -  Nussbaumer & Storey, 1983  
      ----------------------------------------------------------- */

    tmprec = Ne*intData[2][ions-1];
    CoolCoeffs.Crate[ions]   += tmprec;
    CoolCoeffs.Rrate[ions-1] += tmprec;
    CoolCoeffs.Ca[ions]      += intData[2][ions-1];
    CoolCoeffs.Ra[ions-1]    += intData[2][ions-1];

  /*  --------------------------------------------------------------------------
       charge-transfer with H+  - ionization - Kingdon & Ferland 1996, ApJSS  
      -------------------------------------------------------------------------- */

    scrh = N*elem_ab[el_H]*intData[3][ions];
    tmprec = (1.0 - v[X_HI])*scrh;
    CoolCoeffs.Crate[ions]   += tmprec;
    CoolCoeffs.Lrate[ions+1] += tmprec;
    CoolCoeffs.Cb[ions]      -= scrh;
    CoolCoeffs.Lb[ions+1]    -= scrh;

  /* ---------------------------------------------------------------------------
      charge-transfer with H   - recombination - Kingdon & Ferland 1996, ApJSS 
     --------------------------------------------------------------------------- */

    scrh = N*elem_ab[el_H]*intData[4][ions-1];
    tmprec = v[X_HI]*scrh;
    CoolCoeffs.Crate[ions]   += tmprec;
    CoolCoeffs.Rrate[ions-1] += tmprec;
    CoolCoeffs.Cb[ions]      += scrh;
    CoolCoeffs.Rb[ions-1]    += scrh;
  }

 /* -- calculation for last element (ions = NIONS - 1) -- */

  tmprec = Ne*intData[1][ions-1];
  CoolCoeffs.Crate[ions]   += tmprec;
  CoolCoeffs.Rrate[ions-1] += tmprec;
  CoolCoeffs.Ca[ions]      += intData[1][ions-1];
  CoolCoeffs.Ra[ions-1]    += intData[1][ions-1];

  tmprec = Ne*intData[2][ions-1];
  CoolCoeffs.Crate[ions]   += tmprec;
  CoolCoeffs.Rrate[ions-1] += tmprec;
  CoolCoeffs.Ca[ions]      += intData[2][ions-1];
  CoolCoeffs.Ra[ions-1]    += intData[2][ions-1];

  scrh   = N*elem_ab[el_H]*intData[4][ions-1];
  tmprec = v[X_HI]*scrh;
  CoolCoeffs.Crate[ions]   += tmprec;
  CoolCoeffs.Rrate[ions-1] += tmprec;
  CoolCoeffs.Cb[ions]      += scrh;
  CoolCoeffs.Rb[ions-1]    += scrh;

 /* -- calculation for HeI  -- */

  ions = 1;
  tmprec = Ne*intData[0][ions];
  CoolCoeffs.Crate[ions]   += tmprec;
  CoolCoeffs.Lrate[ions+1] += tmprec;
  CoolCoeffs.Ca[ions]   += intData[0][ions];
  CoolCoeffs.La[ions+1] += intData[0][ions];

  scrh   = N*elem_ab[el_H]*intData[3][ions];
  tmprec = (1.0 - v[X_HI])*scrh;
  CoolCoeffs.Crate[ions]   += tmprec;
  CoolCoeffs.Lrate[ions+1] += tmprec;
  CoolCoeffs.Cb[ions]      -= scrh;
  CoolCoeffs.Lb[ions+1]    -= scrh;

/* --------------------------------------------------------------------------
    charge-transfer with He  - recombination  
    ORNL Charge Transfer Database 
    http://www-cfadc.phy.ornl.gov/astro/ps/data/cx/helium/rates/fits.data   
   -------------------------------------------------------------------------- */

  for (ions = 4; ions < NIONS; ions++) {   /* 4 ??????????????? */
    scrh   = N*elem_ab[el_He]*intData[5][ions-1];
    tmprec = v[X_HeI]*scrh;
    CoolCoeffs.Crate[ions]   += tmprec;
    CoolCoeffs.Rrate[ions-1] += tmprec;

    CoolCoeffs.Cc[ions]   += scrh;
    CoolCoeffs.Rc[ions-1] += scrh;
  }
    
/* -------------------------------------------------------
    And now, for ions for which we only have the total 
    electron-ion recombination coefficient:   
    ( NIFS-DATA-54 - Kato & Asano 1999 )
   ------------------------------------------------------- */

  for (ions = 2; ions < NIONS; ions++) {
    tmprec = Ne*intData[6][ions-1];
    CoolCoeffs.Crate[ions]   += tmprec;
    CoolCoeffs.Rrate[ions-1] += tmprec;   

  /* -----------------------------------------------
      Now save the data for the partial derivatives
      computation
     ----------------------------------------------- */

    CoolCoeffs.Ca[ions]   += intData[6][ions-1];
    CoolCoeffs.Ra[ions-1] += intData[6][ions-1];
  } 
  
  /*  He/He+ contribution from charge-exchange processes:  */
/*
  for (ions=3; ions<28; ions++) {
    if (T<5000.) CoolCoeffs.Rrate[2] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a1[ions] * pow( (T/10000.), chtrHe_rec_b1[ions] ) * ( 1. +  chtrHe_rec_c1[ions] * exp ( chtrHe_rec_d1[ions]*(T/10000.) ));
    if ( (T>=5000.) && (T<10000.) ) CoolCoeffs.Rrate[2] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a2[ions] * pow( (T/10000.), chtrHe_rec_b2[ions] ) * ( 1. +  chtrHe_rec_c2[ions] * exp ( chtrHe_rec_d2[ions]*(T/10000.) ));
    if ( (T>=10000.) && (T<50000.) ) CoolCoeffs.Rrate[2] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a3[ions] * pow( (T/10000.), chtrHe_rec_b3[ions] ) * ( 1. +  chtrHe_rec_c3[ions] * exp ( chtrHe_rec_d3[ions]*(T/10000.) ));
    if ( (T>=50000.) && (T<100000.) ) CoolCoeffs.Rrate[2] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a4[ions] * pow( (T/10000.), chtrHe_rec_b4[ions] ) * ( 1. +  chtrHe_rec_c4[ions] * exp ( chtrHe_rec_d4[ions]*(T/10000.) ));
    if (T>=100000.) CoolCoeffs.Rrate[2] += N * v[ions] *elem_ab[elem_part[ions]] *  1.e-9 * chtrHe_rec_a5[ions] * pow( (T/10000.), chtrHe_rec_b5[ions] ) * ( 1. +  chtrHe_rec_c5[ions] * exp ( chtrHe_rec_d5[ions]*(T/10000.) ));

    if (T<5000.) CoolCoeffs.Crate[1] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a1[ions] * pow( (T/10000.), chtrHe_rec_b1[ions] ) * ( 1. +  chtrHe_rec_c1[ions] * exp ( chtrHe_rec_d1[ions]*(T/10000.) ));
    if ( (T>=5000.) && (T<10000.) ) CoolCoeffs.Crate[1] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a2[ions] * pow( (T/10000.), chtrHe_rec_b2[ions] ) * ( 1. +  chtrHe_rec_c2[ions] * exp ( chtrHe_rec_d2[ions]*(T/10000.) ));
    if ( (T>=10000.) && (T<50000.) ) CoolCoeffs.Crate[1] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a3[ions] * pow( (T/10000.), chtrHe_rec_b3[ions] ) * ( 1. +  chtrHe_rec_c3[ions] * exp ( chtrHe_rec_d3[ions]*(T/10000.) ));
    if ( (T>=50000.) && (T<100000.) ) CoolCoeffs.Crate[1] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a4[ions] * pow( (T/10000.), chtrHe_rec_b4[ions] ) * ( 1. +  chtrHe_rec_c4[ions] * exp ( chtrHe_rec_d4[ions]*(T/10000.) ));
    if (T>=100000.) CoolCoeffs.Crate[1] += N * v[ions] * elem_ab[elem_part[ions]] * 1.e-9 * chtrHe_rec_a5[ions] * pow( (T/10000.), chtrHe_rec_b5[ions] ) * ( 1. +  chtrHe_rec_c5[ions] * exp ( chtrHe_rec_d5[ions]*(T/10000.) ));
  }
*/

/* ------------------------------------------------------------- 
       C H E C K 
   ------------------------------------------------------------- */
/*
for (ions = 0; ions < NIONS; ions++) {
  CoolCoeffs.Lrate[ions] *= UNIT_LENGTH/UNIT_VELOCITY;
  CoolCoeffs.Crate[ions] *= UNIT_LENGTH/UNIT_VELOCITY;
  CoolCoeffs.Rrate[ions] *= UNIT_LENGTH/UNIT_VELOCITY;
}

print ("H:     %12.6e  %12.6e  %12.6e\n",
         CoolCoeffs.Lrate[HI-NFLX], CoolCoeffs.Crate[HI-NFLX],CoolCoeffs.Rrate[HI-NFLX]);

print ("--------------------------------------------------------------- \n");
for (nv = HeI; nv <= HeII; nv++){
print ("He(%d): %12.6e  %12.6e  %12.6e\n",nv-HeI+1,
         CoolCoeffs.Lrate[nv-NFLX], CoolCoeffs.Crate[nv-NFLX],CoolCoeffs.Rrate[nv-NFLX]);
}

print ("--------------------------------------------------------------- \n");
for (nv = CI; nv <= CV; nv++){
print ("C(%d):  %12.6e  %12.6e  %12.6e\n",nv-CI+1,
         CoolCoeffs.Lrate[nv-NFLX], CoolCoeffs.Crate[nv-NFLX],CoolCoeffs.Rrate[nv-NFLX]);
}
print ("--------------------------------------------------------------- \n");

for (nv = NI; nv <= NV; nv++){
print ("N(%d):  %12.6e  %12.6e  %12.6e\n",nv-NI+1,
         CoolCoeffs.Lrate[nv-NFLX], CoolCoeffs.Crate[nv-NFLX],CoolCoeffs.Rrate[nv-NFLX]);
}
print ("--------------------------------------------------------------- \n");

for (nv = OI; nv <= OV; nv++){
print ("O(%d):  %12.6e  %12.6e  %12.6e\n",nv-OI+1,
         CoolCoeffs.Lrate[nv-NFLX], CoolCoeffs.Crate[nv-NFLX],CoolCoeffs.Rrate[nv-NFLX]);
}
print ("--------------------------------------------------------------- \n");

for (nv = NeI; nv <= NeV; nv++){
print ("Ne(%d): %12.6e  %12.6e  %12.6e\n",nv-NeI+1,
         CoolCoeffs.Lrate[nv-NFLX], CoolCoeffs.Crate[nv-NFLX],CoolCoeffs.Rrate[nv-NFLX]);
}
print ("--------------------------------------------------------------- \n");
for (nv = SI; nv <= SV; nv++){
print ("Se(%d):  %12.6e  %12.6e  %12.6e\n",nv-SI+1,
         CoolCoeffs.Lrate[nv-NFLX], CoolCoeffs.Crate[nv-NFLX],CoolCoeffs.Rrate[nv-NFLX]);
}


exit(1);
*/
}
 
/* ****************************************************************************** */
double find_N_rho ()
/*
 *
 *  PURPOSE
 *
 *  Find the number density of atoms / ions (divided by the density so that is 
 *  a constant depending on composition only), knowing the ions ratios and density.
 *  We compute a /mu taking into account only the heavy particles (nuclei).
 *
 *
 * LAST MODIFIED: July 18, 2006 by Ovidiu Tesileanu
 *
 ******************************************************************************* */
{
  int i, j;
  double mu, mu1, mu2;
   
  mu1 = mu2 = 0.0; 
  for (i = 0; i < 7; i++) {
    mu1 += elem_ab[i]*elem_mass[i]; /*    Numerator part of mu    */
    mu2 += elem_ab[i];              /*    Denominator part of mu  */
  }
  mu = mu1 / mu2;
  return (UNIT_DENSITY/mu*CONST_NA);  /* This is N/rho, with N the total number density of atoms and ions */
}

/* ********************************************************************* */
double H_MassFrac (void)
/*! 
 * Compute the mass fraction X of Hydrogen as function of the 
 * composition of the gas.
 *
 *             f_H A_H
 *    X = ----------------   
 *         \sum_k f_k A_k
 * 
 * where 
 *
 *    f_k   : is the fractional abundance (by number), 
 *               f_k = N_k/N_tot
 *            of atomic species (no matter ionization degree).
 *
 *    A_K   : is the atomic weight
 *
 *
 *   where N_{tot} is the total number density of particles
 *
 * ARGUMENTS
 *
 *   none
 *
 *********************************************************************** */
{
  double mmw2;
  int    i, j;
  
  mmw2 = 0.0;
  for (i = 0; i < (Fe_IONS==0?7:8); i++) {
    mmw2 += elem_mass[i]*elem_ab[i];  /*    Denominator part of X  */    
  }
  return ( (elem_mass[0]*elem_ab[0]) / mmw2);
}


#ifdef CH_SPACEDIM
/* ********************************************************************* */
void NormalizeIons (double *u)
/*!
 *  This function re-normalize the set of conservative variables u in
 *  such a way that the sum of all ion fractions is equal to 1. 
 *  Since u is a vector of conserved quantities, the normalization 
 *  condition is:
 *
 *      \sum_j (\rho*X_j) = \rho
 * 
 *  where, for instance, X_j = {OI, OII, OIII, OIV, OV}.
 *  In order to fulfill the previous equation, each conservative
 *  variable is redefined according to 
 *
 *                               (\rho*X_j)
 *     (\rho*X_j) -->  \rho * -------------------
 *                             \sum_j (\rho*X_j) 
 *
 *  This function is necessary only when used with Chombo AMR since
 *  the coarse-to-fine interpolation fails to conserve the correct
 *  normalization. Re-fluxing and fine-to-coarse do not need this
 *  treatment since the normalization is preserved.
 *
 *********************************************************************** */
{
  int nv;
  double phi;

  phi = 0.0;
  for (nv = X_HeI; nv <= X_HeII; nv++) phi   += u[nv];
  for (nv = X_HeI; nv <= X_HeII; nv++) u[nv] *= u[RHO]/phi;

  phi = 0.0;
  for (nv = X_SI; nv < X_SI+S_IONS; nv++) phi   += u[nv];
  for (nv = X_SI; nv < X_SI+S_IONS; nv++) u[nv] *= u[RHO]/phi;

  phi = 0.0;
  for (nv = X_OI; nv < X_OI+O_IONS; nv++) phi   += u[nv];
  for (nv = X_OI; nv < X_OI+O_IONS; nv++) u[nv] *= u[RHO]/phi;

  phi = 0.0;
  for (nv = X_NI; nv < X_NI+N_IONS; nv++) phi   += u[nv];
  for (nv = X_NI; nv < X_NI+N_IONS; nv++) u[nv] *= u[RHO]/phi;

  phi = 0.0;
  for (nv = X_CI; nv < X_CI+C_IONS; nv++) phi   += u[nv];
  for (nv = X_CI; nv < X_CI+C_IONS; nv++) u[nv] *= u[RHO]/phi;
    
  phi = 0.0;
  for (nv = X_NeI; nv < X_NeI+Ne_IONS; nv++) phi   += u[nv];
  for (nv = X_NeI; nv < X_NeI+Ne_IONS; nv++) u[nv] *= u[RHO]/phi;

}
#endif

