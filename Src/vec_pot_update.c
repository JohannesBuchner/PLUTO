/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief Function to handle vector potential.
 
  Update the cell-centered vector potential by integrating
  \f[
     \pd{\vec{A}}{t} + \vec{E} = 0
  \f] 
  Depending on the time-stepping algorithm, we advance the solution
  in time using the following increments
  \f[
   \begin{array}{ll}
    \rm{RK2:}\quad &\DS A^{n+1} = A^n + \frac{\Delta t^n}{2} (E^n + E^*) 
      \\ \noalign{\medskip}
    \rm{RK3:}\quad &\DS A^{n+1} = A^n + \frac{\Delta t^n}{6} (E^n + E^* + 4E^{**}) 
      \\ \noalign{\medskip}
    \rm{CTU:}\quad &\DS A^{n+1} = A^n + \Delta t^nE^{n+\HALF}
   \end{array}        
  \f]
  The cell-centered electric field \b E is computed by averaging values 
  available at faces or edges, depending on the MHD formulation:
  - staggered_mhd: we average the electric field values 
                   available at cell edges and previously computed
                   . In this
                   case the function does the job only 
                   if called from CT;
  - Cell-cent. mhd: we average the face values using upwind fluxes. 
                    In this case the function does the job only 
                    if called from time stepping algorithms.
                    The flux components of induction equation
                    contain the emf: 
                    \f$ F_{i+\HALF} = (   0,-E_z, E_y)_{i+\HALF}\f$,
                    \f$ F_{j+\HALF} = ( E_z,   0,-E_x)_{j+\HALF}\f$,
                    \f$ F_{k+\HALF} = (-E_y, E_x,   0)_{k+\HALF}\f$,
                    
 
  Note that the vector potential is NOT used in the rest
  of the code and serves solely as a diagnostic tool.
 
  \param [in,out] d   pointer to PLUTO Data structure
  \param [in] vp      pointer to an EMF structure, used only with CT
  \param [in] sweep   pointer to Sweep structure containing the 
                      electric field components (cell-centered MHD)
  \param [in] grid    pointer to an array of Grid structures
 
  \attention In cylindrical coordinates:
  PLUTO treats cylindrical (but not polar) coordinates as 
  Cartesian so that \f$ (x_1,x_2) = (R,z) \f$ .
  Nevertheless, the correct right-handed cylindrical system is 
  defined by \f$ (x_1,x_2,x_3) = (R,\phi,z)\f$. 
  The fact that the \f$ \phi \f$ coordinate is missing causes an 
  ambiguity when computing the third component of the EMF in 2.5 D:
  \f[
     E_3 = v_2B_1 - v_1B_2
  \f]
  which is (wrongly) called \f$ E_z \f$  but is, indeed, 
  \f$ -E_\phi \f$ .
  This is not a problem in the update of the magnetic field, since
  the contribution coming from \f$ E_3 \f$ is also reversed:
  \f[
  \Delta B_1 = -\Delta t dE_3/dx_2 \quad\Longrightarrow\quad
  \Delta B_r =  \Delta t dE_\phi/dz
  \f]
  as it should be.
  However, it is a problem in the definition of the vector potential
  which is consistently computed using its physical definition.
  Therefore, to correctly update the phi component of the vector
  potential we change sign to \f$ E_z \f$ (\f$ = -E_\phi \f$).
 
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date   Sep 24, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if (PHYSICS == MHD || PHYSICS == RMHD) && (UPDATE_VECTOR_POTENTIAL == YES)
/* ********************************************************************* */
void VectorPotentialUpdate (const Data *d, const void *vp, 
                            const Sweep *sweep, const Grid *grid)
/*!
 *
 *********************************************************************** */
{
  int    i,j,k;
  double scrh, dt, sEz = 0.25;
  double Ex, Ey, Ez;

/* ---------------------------------------------------------
    When this function is called from Time stepping 
    routines and CT is used, do nothing.
   --------------------------------------------------------- */

  #ifdef STAGGERED_MHD
   if (vp == NULL) return;
  #endif

  #if GEOMETRY == CYLINDRICAL 
   sEz = -0.25;
  #endif

/* ------------------------------------------------------
     Determine dt according to the integration step
   ------------------------------------------------------ */

  dt = g_dt;
  #if TIME_STEPPING == RK2
   dt *= 0.5;
  #elif TIME_STEPPING == RK3
   dt /= 6.0;
   if (g_intStage == 3) dt *= 4.0;
  #endif

  #ifdef CTU  
   if (g_intStage == 1) return;
  #endif

  #ifdef STAGGERED_MHD
/* -----------------------------------------------------
     If staggered fields are used, we average the 
     EMF available at cell edges:

     A-----------B
     |           |
     |           |
     |     o     |    o = (A+B+C+D)/4
     |           |
     |           |
     D-----------C
   ----------------------------------------------------- */
  {
    EMF *emf;
    emf = (EMF *) vp;

    DOM_LOOP(k,j,i){
      Ez = sEz*(emf->ez[k][j][i]     + emf->ez[k][j-1][i] + 
                emf->ez[k][j-1][i-1] + emf->ez[k][j][i-1]);

      #if DIMENSIONS == 3
       Ex = 0.25*(emf->ex[k][j][i]     + emf->ex[k][j-1][i] + 
                  emf->ex[k-1][j-1][i] + emf->ex[k-1][j][i]);

       Ey = 0.25*(emf->ey[k][j][i]     + emf->ey[k][j][i-1] + 
                  emf->ey[k-1][j][i-1] + emf->ey[k-1][j][i]);
      #endif

      d->Ax3[k][j][i] -= dt*Ez; 
      #if DIMENSIONS == 3
       d->Ax1[k][j][i] -= dt*Ex; 
       d->Ax2[k][j][i] -= dt*Ey; 
      #endif

    }   
  }
  #else 
/* ------------------------------------------------------------
     If cell-center is used, we average the 
     EMF available at cell faces:

     +-----B-----+
     |           |
     |           |
     A     o     C    o = (A+B+C+D)/4
     |           |
     |           |
     +-----D-----+
     
     The components of the flux in the induction equation give 
     the emf:
     
     when g_dir == IDIR we have  F(i+1/2) = (0,-Ez, Ey)
     when g_dir == JDIR we have  F(j+1/2) = (Ez,0,-Ex)
     when g_dir == KDIR we have  F(k+1/2) = (-Ey,Ex,0)
     
   ------------------------------------------------------------ */
  {
    double **f;

    f = sweep->flux;
    if (g_dir == IDIR){

      IDOM_LOOP(i){
        d->Ax3[g_k][g_j][i]  -= -sEz*dt*(f[i][BX2] + f[i-1][BX2]); /* -Ez */
        #if DIMENSIONS == 3
         d->Ax2[g_k][g_j][i] -= 0.25*dt*(f[i][BX3] + f[i-1][BX3]); /* Ey */
        #endif
      }
    }else if (g_dir == JDIR){

      JDOM_LOOP(j) {
        #if DIMENSIONS == 3
         d->Ax1[g_k][j][g_i] -= -0.25*dt*(f[j][BX3] + f[j-1][BX3]); /* -Ex */
        #endif
        d->Ax3[g_k][j][g_i]  -= sEz*dt*(f[j][BX1] + f[j-1][BX1]);  /* Ez */
      }

    }else if (g_dir == KDIR){
      KDOM_LOOP(k){
        d->Ax1[k][g_j][g_i] -=  0.25*dt*(f[k][BX2] + f[k-1][BX2]);
        d->Ax2[k][g_j][g_i] -= -0.25*dt*(f[k][BX1] + f[k-1][BX1]);
      }
    }
  }
  #endif /* STAGGERED_MHD */
}
#endif  /* UPDATE_VECTOR_POTENTIAL == YES */
