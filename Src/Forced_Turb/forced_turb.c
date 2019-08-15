#include "pluto.h"

/* ************************************************************* */
void ForcedTurb_Init(ForcedTurb *Ft)
/*!
 *   Initialize the modes.
 *
 * ************************************************************* */
{
  int s, dir, ikx, iky, ikz, ikxmax=1, ikymax=1, ikzmax=1;
  double kx, ky, kz, kc, k;
  double Lx, Ly, Lz;
  double amin = 0.0;
  
  Ft->NModes = 0;
  Ft->SpectForm = 1;
  Ft->StirFreq  = FORCED_TURB_FREQ;
  Ft->StirDecay = FORCED_TURB_DECAY;
  Ft->StirEnergy = FORCED_TURB_ENERGY;
  Ft->StirKMin = FORCED_TURB_KMIN;
  Ft->StirKMax = FORCED_TURB_KMAX;
  Ft->SolveWt = FORCED_TURB_WEIGHT; /* Purely solenoidal for 1.0 and Purely dilational for 0 */
  Ft->Mode = ARRAY_2D(3, ForcedTurb_MAXMODES, double);
  Ft->aka  = ARRAY_2D(3, ForcedTurb_MAXMODES, double);
  Ft->akb  = ARRAY_2D(3, ForcedTurb_MAXMODES, double);
  Ft->OUPhases = ARRAY_1D(6*(ForcedTurb_MAXMODES), double);
  Ft->StirAmpl = ARRAY_1D(ForcedTurb_MAXMODES, double);
  Ft->Acc = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
  int a,b,c;
  TOT_LOOP(a,b,c) for (dir=0; dir<3; dir++) Ft->Acc[dir][a][b][c] = 0.0;
  
  /* Initialize them to 0.0 */
  for(s=0;s<ForcedTurb_MAXMODES;s++){
    Ft->StirAmpl[s] = 0.0;
    for(dir=0;dir<3;dir++){
      Ft->Mode[dir][s] = 0.0;
      Ft->aka[dir][s] = 0.0;
      Ft->akb[dir][s] = 0.0;
    }
  }
  
  for(s=0;s<(6*ForcedTurb_MAXMODES);s++) Ft->OUPhases[s] = 0.0;
  
  Ft->OUVar = sqrt(Ft->StirEnergy/Ft->StirDecay);
  kc = 0.5*(Ft->StirKMax + Ft->StirKMin);

  Ft->SolveWtNorm = sqrt(3.0/1.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*Ft->SolveWt
                                                       + 1.0*Ft->SolveWt*Ft->SolveWt);
  #if DIMENSIONS == 2
  Ft->SolveWtNorm = sqrt(3.0/2.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*Ft->SolveWt
                                                   + 2.0*Ft->SolveWt*Ft->SolveWt);
  #elif DIMENSIONS == 3
    Ft->SolveWtNorm = sqrt(3.0/3.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*Ft->SolveWt
                                                            + 3.0*Ft->SolveWt*Ft->SolveWt);
  #endif

  D_EXPAND(ikxmax = 20; ,
          ikymax = 20; ,
          ikzmax = 20; );

  Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  Lz = g_domEnd[KDIR] - g_domBeg[KDIR];
  
  for(ikx = 0; ikx < ikxmax; ikx++)
  {
    kx = 2.0*CONST_PI*ikx/Lx;
    
    for (iky = 0; iky < ikymax; iky++)
    {
      ky = 2.0*CONST_PI*iky/Ly;
      
      for (ikz = 0; ikz < ikzmax; ikz++){
        
        kz = 2.0*CONST_PI*ikz/Lz;
        k = sqrt(kx*kx + ky*ky + kz*kz);
        
        if ((k >= Ft->StirKMin) && (k <= Ft->StirKMax)){
          
          if (Ft->SpectForm == 0) Ft->StirAmpl[Ft->NModes]  = 1.0;
          if (Ft->SpectForm == 1) Ft->StirAmpl[Ft->NModes]  = 4.0*(amin - 1.0)/(pow((Ft->StirKMax - Ft->StirKMin),2.0))*pow((kc - k),2.0) + 1.0;
          
          Ft->Mode[IDIR][Ft->NModes] = kx;
          Ft->Mode[JDIR][Ft->NModes] = ky;
          Ft->Mode[KDIR][Ft->NModes] = kz;
          Ft->NModes += 1;
          
#if DIMENSIONS >= 2
        
          if (Ft->SpectForm == 0) Ft->StirAmpl[Ft->NModes]  = 1.0;
          if (Ft->SpectForm == 1) Ft->StirAmpl[Ft->NModes]  = 4.0*(amin - 1.0)/(pow((Ft->StirKMax - Ft->StirKMin),2.0))*pow((kc - k),2.0) + 1.0;
          
          Ft->Mode[IDIR][Ft->NModes] = kx;
          Ft->Mode[JDIR][Ft->NModes] = -ky;
          Ft->Mode[KDIR][Ft->NModes] = kz;
          Ft->NModes += 1;
#endif
          
#if DIMENSIONS == 3
          
          if (Ft->SpectForm == 0) Ft->StirAmpl[Ft->NModes]  = 1.0;
          if (Ft->SpectForm == 1) Ft->StirAmpl[Ft->NModes]  = 4.0*(amin - 1.0)/(pow((Ft->StirKMax - Ft->StirKMin),2.0))*pow((kc - k),2.0) + 1.0;
          
          Ft->Mode[IDIR][Ft->NModes] = kx;
          Ft->Mode[JDIR][Ft->NModes] = ky;
          Ft->Mode[KDIR][Ft->NModes] = -kz;
          Ft->NModes += 1;
          
          
          if (Ft->SpectForm == 0) Ft->StirAmpl[Ft->NModes]  = 1.0;
          if (Ft->SpectForm == 1) Ft->StirAmpl[Ft->NModes]  = 4.0*(amin - 1.0)/(pow((Ft->StirKMax - Ft->StirKMin),2.0))*pow((kc - k),2.0) + 1.0;
          
          Ft->Mode[IDIR][Ft->NModes] = kx;
          Ft->Mode[JDIR][Ft->NModes] = -ky;
          Ft->Mode[KDIR][Ft->NModes] = -kz;
          Ft->NModes += 1;
#endif
        }
      }
    }
  }
  
/* -- Seed random number sequence using SeedGenerator to have
      different sequences on different processors               -- */

  static int first_call = 1;
  if (first_call){
    RandomSeed(SeedGenerator(),0);
    first_call = 0;
  }
}

/* ********************************************************************* */
double ForcedTurb_GenRandNum()
/*!
 *  Generate Random Number from a Gaussian Distribution.
 *
 *********************************************************************** */
{
  double r1 = RandomNumber(0.0,1.0);
  double r2 = RandomNumber(0.0,1.0);

  return sqrt(2.0*log(1./r1))*cos(2.0*CONST_PI*r2);
}

/* ********************************************************************* */
void ForcedTurb_OUNoiseInit(double *InVec, int VecLength, double Variance)
/*!
 *
 *  Initialize the Ornstein Uhlenbeck process.
 *
 *********************************************************************** */
{
  int i;
  double GenRan;
  for(i=0;i<VecLength;i++){
    GenRan   = ForcedTurb_GenRandNum();
    InVec[i] = GenRan*Variance;
  }
  
}
/* ********************************************************************* */
void ForcedTurb_OUNoiseUpdate(double *InVec, int VecLength,
                              double Variance, double dt, double ts)
/*!
 *
 *  Update the OU process.
 *
 *********************************************************************** */
{
  int i;
  double GenRan, DampFact = exp(-dt/ts), fact;
    
  if (g_intStage == 1) DampFact = 1.0;
  else                 DampFact = exp(-dt/ts);
  
  fact = Variance*sqrt(1.0 - DampFact*DampFact);
  for(i=0; i<VecLength; i++){
    GenRan = ForcedTurb_GenRandNum();
    InVec[i] = InVec[i]*DampFact + fact*GenRan;
  }
}

/* ********************************************************************* */
void ForcedTurb_CalcPhases(ForcedTurb *Ft)
/*!
 * Calculates the phases based on modes and type of forcing.
 * SolveWt = 0 --> Purely compressive (curl F = 0)
 * SolveWt = 1 --> Purely solenoidal  (div F = 0)
 *
 * ********************************************************************* */
{
  int i,j;
  double ka, kb, kk, diva, divb, curla, curlb;
  
  for(i=0;i<Ft->NModes; i++){
    ka = 0.0;
    kb = 0.0;
    kk = 0.0;
    for(j=0;j<DIMENSIONS;j++){
      kk = kk + Ft->Mode[j][i]*Ft->Mode[j][i];
      ka = ka + Ft->Mode[j][i]*Ft->OUPhases[6*i+2*j+0+1];
      kb = kb + Ft->Mode[j][i]*Ft->OUPhases[6*i+2*j+0+0];
    }
    
    
    for(j=0;j<DIMENSIONS;j++){
      diva  = Ft->Mode[j][i]*ka/kk;
      divb  = Ft->Mode[j][i]*kb/kk;
      curla = (Ft->OUPhases[6*i+2*j+0+0] - divb);
      curlb = (Ft->OUPhases[6*i+2*j+1+0] - diva);
      
      Ft->aka[j][i] = Ft->SolveWt*curla + (1.0-Ft->SolveWt)*divb;
      Ft->akb[j][i] = Ft->SolveWt*curlb + (1.0-Ft->SolveWt)*diva;
    }
  }
}

/* ********************************************************************* */
void ForcedTurb_ComputeAcceleration(ForcedTurb *Ft, Grid *grid)
/*!
 * Compute the acceleration at each grid point in all directions.
 *
 *********************************************************************** */
{
  int m, dir,k,j,i;
  double x1, x2, x3;
  double arg_x, arg_y, arg_z;
  double coskx, cosky, coskz, sinkx, sinky, sinkz;
  static double ****realterm, ****imagterm;
  
/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (realterm == NULL){
     realterm = ARRAY_4D(Ft->NModes, NX3_TOT, NX2_TOT, NX1_TOT, double);
     imagterm = ARRAY_4D(Ft->NModes, NX3_TOT, NX2_TOT, NX1_TOT, double);
     TOT_LOOP(k,j,i){
	      x1 = grid->x[IDIR][i];
        x2 = grid->x[JDIR][j];
        x3 = grid->x[KDIR][k];
	 
      	for (m = 0; m<Ft->NModes;m++){
          arg_x = Ft->Mode[IDIR][m]*x1;
          coskx = cos(arg_x);
          sinkx = sin(arg_x);
          
          arg_y = Ft->Mode[JDIR][m]*x2;
          cosky = cos(arg_y);
          sinky = sin(arg_y);
          
        #if DIMENSIONS == 3
          arg_z = Ft->Mode[KDIR][m]*x3;
          coskz = cos(arg_z);
          sinkz = sin(arg_z);
        #else
          coskz = 1.0;
          sinkz = 0.0;
        #endif
        realterm[m][k][j][i] = (coskx*cosky - sinkx*sinky)*coskz- (sinkx*cosky + coskx*sinky)*sinkz;
        imagterm[m][k][j][i] = coskx*(cosky*sinkz+ sinky*coskz) + sinkx*(cosky*coskz - sinky*sinkz);
     }
    }
  }

/* --------------------------------------------------------
   1. Compute acceleration term
   -------------------------------------------------------- */
 
  TOT_LOOP(k,j,i){
    for(dir = 0; dir < 3; dir++) Ft->Acc[dir][k][j][i] = 0.0;
    x1 = grid->x[IDIR][i];
    x2 = grid->x[JDIR][j];
    x3 = grid->x[KDIR][k];
    for (m = 0; m<Ft->NModes;m++){
      for(dir = 0; dir < 3; dir++){
        Ft->Acc[dir][k][j][i] += 2.0*Ft->StirAmpl[m]
                                 *(  Ft->aka[dir][m]*realterm[m][k][j][i]
                                   - Ft->akb[dir][m]*imagterm[m][k][j][i]);
      }
    }

    for(dir = 0; dir < 3; dir++)  Ft->Acc[dir][k][j][i] *= Ft->SolveWtNorm;
  }
}

/* ********************************************************************* */
void ForcedTurb_CorrectRHS(const Data *d, const Sweep *sweep, int beg, 
                           int end, double dt, Grid *grid)
/*!
 *  Correct the RHS term in update stage to add the acceleration
 *  term to momentum and energy flux, sweep by sweep.
 *
 *********************************************************************** */
{
  int i, j, k;
  double at;
  double **flux, **rhs, *vg;
  const State *stateC = &(sweep->stateC);
  
  flux = sweep->flux;
  rhs  = sweep->rhs;
  ForcedTurb *Ft;
  Ft = d->Ft;

 if (g_dir == IDIR){
    j = g_j;
    k = g_k;
    for(i = beg; i <= end; i++){
      at = Ft->Acc[IDIR][k][j][i];
      vg = stateC->v[i];
      rhs[i][MX1] += dt*vg[RHO]*at;
      IF_ENERGY(rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*at;)
    }
  }
  if (g_dir == JDIR){
    i = g_i;
    k = g_k;
    for(j = beg; j <= end; j++){
        at = Ft->Acc[JDIR][k][j][i];
        vg = stateC->v[j];
        rhs[j][MX2] += dt*vg[RHO]*at;
        IF_ENERGY(rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*at;)
    }
  }
  
  if (g_dir == KDIR){
    i = g_i;
    j = g_j;
    for(k = beg; k <= end; k++){
        at = Ft->Acc[KDIR][k][j][i];
        vg = stateC->v[k];
        rhs[k][MX3] += dt*vg[RHO]*at;
        IF_ENERGY(rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*at;)
    }
  }
}

