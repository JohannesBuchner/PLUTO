/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the dust module.

  Contains variable names and prototypes for the dust module

  \author A. Mignone (mignone@ph.unito.it)
  \date   April, 2, 2015
  
  TODO:
  
   - Implicit implementation of drag coupling term.
   - More careful characteristic analysis (Miniati 2010)
*/
/* ///////////////////////////////////////////////////////////////////// */

#define  RHO_D   (NDUST_FLUID_BEG)
#define  MX1_D   (RHO_D + 1)
#define  MX2_D   (COMPONENTS >= 2 ? (RHO_D + 2):255)
#define  MX3_D   (COMPONENTS == 3 ? (RHO_D + 3):255)
#define  VX1_D   MX1_D
#define  VX2_D   MX2_D
#define  VX3_D   MX3_D

#if GEOMETRY == POLAR
  #define iMPHI_D  MX2_D
#elif GEOMETRY == SPHERICAL
  #define iMPHI_D  MX3_D
#endif

#define NDUST_FLUID (1 + COMPONENTS)


void Dust_Solver (const Sweep *, int, int, double *, Grid *);
void Dust_DragForce(const Sweep *, int, int, double, Grid *);
