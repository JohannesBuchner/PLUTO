#include "pluto.h"
#include <sys/stat.h>
#include <sys/types.h>

static void Particles_Log (Data *, timeStep *, Grid *);

/* ********************************************************************* */
void OutputLogPre(Data *data, timeStep *Dts, Runtime *ini, Grid *grid)
/*!
 *  Provide log-file information *before* integration starts
 *********************************************************************** */
{
  char *indent = IndentString();
  print ("step:%d; t = %10.4e; dt = %10.4e; %3.1f %%\n",
           g_stepNumber, g_time, g_dt, 100.0*g_time/ini->tstop);
}

/* ********************************************************************* */
void OutputLogPost(Data *data, timeStep *Dts, Runtime *ini, Grid *grid)
/*!
 *  Provide log-file information *after* integration ends
 *
 *********************************************************************** */
{
  char *indent = IndentString();
  
  print ("%s [Mach = %f", indent, g_maxMach);
  if (g_maxRiemannIter > 0){
    print (", NRiemann = %d", g_maxRiemannIter);
  }
  
  #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
/*  print ("%s [Nsts = %d]\n",indent,Dts->Nsts); */
  print (", Nsts = %d",Dts->Nsts);
  #endif
  #if (PARABOLIC_FLUX & RK_LEGENDRE)
/*  print ("%s [Nrkl = %d]\n",indent,Dts->Nrkl); */
  print (", Nrkl = %d",Dts->Nrkl);
  #endif
  #if (PHYSICS == RMHD) && (RESISTIVITY == IMEX)
  print (", Nimex = %d",g_maxIMEXIter);
  #endif
  print ("]\n");
  #ifdef PARTICLES
  Particles_Log (data, Dts, grid);
  #endif
}

/* ********************************************************************* */
void Particles_Log (Data *data, timeStep *Dts, Grid *grid)
/*!
 *  Global MPI reduction operations for PARTICLES Diagnostics
 *
 *********************************************************************** */
{
#ifdef PARTICLES
  int n;
  long int np_glob;
  particleNode *CurNode;
  Particle *pp;

  CurNode = data->PHead;
  double kin, kin_glob;
  #if PARTICLES_LP_SPECTRA == YES
   double sEmin, sEmin_glob, sEmax, sEmax_glob;
   sEmin  = 0.0;
   sEmax  = 0.0;
  #endif

  kin = 0.0;
  while(CurNode != NULL){
    pp   = &(CurNode->p);
    kin += 0.5*(EXPAND(  pp->speed[IDIR]*pp->speed[IDIR],  
                       + pp->speed[JDIR]*pp->speed[JDIR],
                       + pp->speed[KDIR]*pp->speed[KDIR]));
    #if PARTICLES_LP_SPECTRA == YES
    sEmin  += pp->eng[0];
    sEmax  += pp->eng[PARTICLES_LP_NEBINS-1];
    #endif

    CurNode = CurNode->next;
  }

#ifdef PARALLEL
  MPI_Allreduce(&p_nparticles, &np_glob, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&kin, &kin_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
  kin = kin_glob/(np_glob+1.e-6);  /* Avoid division by zero when
                                      there're no particles */
  #if PARTICLES_LP_SPECTRA == YES
  MPI_Allreduce(&sEmin, &sEmin_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sEmin = sEmin_glob/(np_glob+1.e-6);  /* Avoid division by zero when there're no particles */
  MPI_Allreduce(&sEmax, &sEmax_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sEmax = sEmax_glob/(np_glob+1.e-6);  /* Avoid division by zero when there're no particles */
  #endif
  
#else
  np_glob = p_nparticles;
  kin    /= p_nparticles + 1.e-6;  /* Avoid division by zero when
                                      there're no particles */
   #if PARTICLES_LP_SPECTRA == YES
   sEmin  /= p_nparticles + 1.0e-6;
   sEmax  /= p_nparticles + 1.0e-6;
   #endif
#endif

/* -- print both local and global number of particles -- */
  print ("%s [Nparticles/tot: %ld / %ld; Nsub = %d; <Ek> = %10.4e]\n",
           IndentString(), p_nparticles, np_glob, Dts->Nsub_particles, kin);

#if PARTICLES_LP_SPECTRA == YES
  print ("%s [<SpecE_min> = %10.4e, <SpecE_max> = %10.4e]\n",
          IndentString(), sEmin, sEmax);
#endif
/* -- print just global number of particles --  */

//  print ("%s [Nparticles (tot) = %ld; Nsub = %d; <Ek> = %10.4e]\n",
//             IndentString(), np_glob, Dts->Nsub_particles, kin);

#endif
}

/* /////////////////////////////////////////////////////////////////////
    The next set of functions provides basic functionalities to
     
     - set the log file
     - formatted output to the log file through the print() and print()
       functions
   ///////////////////////////////////////////////////////////////////// */

static char log_file_name[512];

/* ********************************************************************* */
int SetLogFile(char *log_dir, Cmd_Line *cmd)
/*!
 * Set log file name in parallel mode.
 * Each processor has its own log file "pluto.prank.log" 
 *
 * \param [in] output_dir  the name of the output directory
 * \param [in] cmd         pointer to cmd line option structure.
 *                           
 *********************************************************************** */
{
#ifdef PARALLEL
  FILE *fl;

  //if (LOG_FILE_FOLDER ) {
  //  mkdir("Log_Files",0777);
  //}
 
  sprintf (log_file_name, "%s/pluto.%d.log",log_dir,prank);

  if (cmd->restart == NO && cmd->h5restart == NO){
    fl = fopen(log_file_name,"w");
  }else{
    fl = fopen(log_file_name,"aw");
  } 

  /* -- check that we have a valid directory name -- */

  if (fl == NULL){
    printf ("! SetLogFile(): %s cannot be written.\n",log_file_name);
    printf ("  Using current directory instead.\n");
    sprintf (log_file_name, "./pluto.%d.log",prank);
    if (cmd->restart == NO && cmd->h5restart == NO){
      fl = fopen(log_file_name,"w");
    }else{
      fl = fopen(log_file_name,"aw");
    } 
  }
  fprintf(fl,"\n");
  fclose(fl);
#endif
}

#ifndef CHOMBO
/* ********************************************************************* */
void print (const char *fmt, ...)
/*!
 * Define print function for the static grid version
 * of PLUTO. The Chombo version is defined in Chombo/amrPLUTO.cpp
 *
 *********************************************************************** */
{
  FILE *fl;

  va_list args;
  va_start(args, fmt);

#ifdef PARALLEL
  fl = fopen(log_file_name,"a");
  vfprintf(fl, fmt, args);
  fclose (fl);
#else
  vprintf(fmt, args);
#endif

  va_end(args);
}
#endif

/* ********************************************************************* */
char *IndentString()
/*
 *
 *********************************************************************** */
{ 
  static char str[64];

  if      (g_stepNumber < 10)     sprintf (str,"%7s"," "); 
  else if (g_stepNumber < 100)    sprintf (str,"%8s"," ");
  else if (g_stepNumber < 1000)   sprintf (str,"%9s"," ");
  else if (g_stepNumber < 10000)  sprintf (str,"%10s"," ");
  else if (g_stepNumber < 100000) sprintf (str,"%11s"," ");
  else                            sprintf (str,"%12s"," ");

  return str;
}

