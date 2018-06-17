/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Main output driver for fluid variables.

  WriteData() is the main driver for writing data arrays in any of
  the available formats (binary, VTK, HDF5, etc...).  

  - For .dbl, .flt or .vtk file formats, access to binary files is 
    provided by the functions in bin_io.c.
  - HDF5 files are handled by hdf5_io.c.
  - image files are handled by write_img.c
  - tabulated ascii files are handled by write_tab.c

  This function also updates the corresponding .out file associated 
  with the output data format.

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)

  \date   Sep 30, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void BOV_Header(Output *output, char *fdata);

/* ********************************************************************* */
void WriteData (const Data *d, Output *output, Grid *grid)
/*!
 * Write data to disk using any of the available formats.
 *
 * \param [in] d      pointer to PLUTO Data structre 
 * \param [in] output the output structure corresponding to a given
 *                    format
 * \param [in] grid   pointer to an array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nv;
  int    single_file;
  size_t dsize;
  char   filename[512], sline[512];
  static int last_computed_var = -1;
  double units[MAX_OUTPUT_VARS]; 
  float ***Vpt3;
  void *Vpt;
  FILE *fout, *fbin;
  time_t tbeg, tend;
  long long offset;

/* -----------------------------------------------------------
   0. Increment the file number and initialize units
   ----------------------------------------------------------- */

  output->nfile++;

  print ("> Writing file #%d (%s) to disk...", output->nfile, output->ext);

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (prank == 0) time(&tbeg);
#endif

  for (nv = 0; nv < MAX_OUTPUT_VARS; nv++) units[nv] = 1.0;
  if (output->cgs) GetCGSUnits(units);

/* --------------------------------------------------------
   1. Compute user-defined output arrays if necessary 
   -------------------------------------------------------- */

  if (last_computed_var != g_stepNumber && d->Vuser != NULL) {
    ComputeUserVar (d, grid);
    last_computed_var = g_stepNumber;
  }

/* --------------------------------------------------------
   2. With FARGO, we can output total or residual velocity 
   -------------------------------------------------------- */

#if (defined FARGO)
  #if  (FARGO_OUTPUT_VTOT == YES)
  if (FARGO_TotalVelocityIsSet() == 0) FARGO_AddVelocity (d,grid); 
  #else 
  if (FARGO_TotalVelocityIsSet() == 1) FARGO_SubtractVelocity (d,grid); 
  #endif
#endif

/* --------------------------------------------------------
   3. Select the output type 
   -------------------------------------------------------- */

  if (output->type == DBL_OUTPUT) {

  /* ------------------------------------------------------------------- */
  /*! - \b DBL output:
        Double-precision data files can be written using single or
        multiple file mode. 
        - for single file, serial: we open the file just once before
          the main variable loop, dump variables and then close.
        - for single file, parallel the distributed array descriptor sz is
          different for cell-centered or staggered data type and we
          thus have to open and close the file after each variable
          has been dumped.
        - when writing multiple files we open, write to and close the
          file one each loop cycle.
        \note In all cases, the pointer to the data array that has to be 
              written must be cast into (void *) and the starting index of 
              the array must be zero.
  */
  /* ------------------------------------------------------------------- */

    int sz;
    single_file = strcmp(output->mode,"single_file") == 0;
    dsize = sizeof(double);

    if (single_file){  /* -- single output file -- */

      sprintf (filename, "%s/data.%04d.%s", output->dir,output->nfile, 
                                            output->ext);
 
      FileDelete (filename);  /* Avoid partial fill of pre-existing files */
      offset = 0;
      #ifndef PARALLEL
       fbin = FileOpen (filename, 0, "w");
      #endif
      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;

        if      (output->stag_var[nv] == -1) {  /* -- cell-centered data -- */
          sz = SZ;
          Vpt = (void *)output->V[nv][0][0];
        } else if (output->stag_var[nv] == 0) { /* -- x-staggered data -- */
          sz  = SZ_stagx;
          Vpt = (void *)(output->V[nv][0][0]-1);
        } else if (output->stag_var[nv] == 1) { /* -- y-staggered data -- */
          sz = SZ_stagy;
          Vpt = (void *)output->V[nv][0][-1];
        } else if (output->stag_var[nv] == 2) { /* -- z-staggered data -- */
           sz = SZ_stagz;
           Vpt = (void *)output->V[nv][-1][0];
        }
        #ifdef PARALLEL
         fbin = FileOpen (filename, sz, "w");
         AL_Set_offset(sz, offset);
        #endif
        FileWriteData (Vpt, dsize, sz, fbin, output->stag_var[nv]);
        #ifdef PARALLEL
         offset = AL_Get_offset(sz);
         FileClose(fbin, sz);
        #endif
      }
      #ifndef PARALLEL
       FileClose(fbin, sz);
      #endif

    }else{              /* -- multiple files -- */

      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;
        sprintf (filename, "%s/%s.%04d.%s", output->dir, output->var_name[nv], 
                                            output->nfile, output->ext);

        FileDelete (filename);  /* Avoid partial fill of pre-existing files */
        if      (output->stag_var[nv] == -1) {  /* -- cell-centered data -- */
          sz = SZ;
          Vpt = (void *)output->V[nv][0][0];
        } else if (output->stag_var[nv] == 0) { /* -- x-staggered data -- */
          sz  = SZ_stagx;
          Vpt = (void *)(output->V[nv][0][0]-1);
        } else if (output->stag_var[nv] == 1) { /* -- y-staggered data -- */
          sz = SZ_stagy;
          Vpt = (void *)output->V[nv][0][-1];
        } else if (output->stag_var[nv] == 2) { /* -- z-staggered data -- */
           sz = SZ_stagz;
           Vpt = (void *)output->V[nv][-1][0];
        }
        fbin = FileOpen (filename, sz, "w");
        FileWriteData (Vpt, dsize, sz, fbin, output->stag_var[nv]);
        FileClose (fbin, sz);
      }
    }

  } else if (output->type == FLT_OUTPUT) {

  /* ----------------------------------------------------------
                 FLT output for cell-centered data
     ---------------------------------------------------------- */

    single_file = strcmp(output->mode,"single_file") == 0;

    if (single_file){  /* -- single output file -- */
      sprintf (filename, "%s/data.%04d.%s", output->dir, output->nfile, 
                                            output->ext);
      FileDelete (filename);  /* Avoid partial fill of pre-existing files */
      fbin = FileOpen (filename, SZ_float, "w");
      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;
        Vpt3 = Convert_dbl2flt(output->V[nv], units[nv],0);
        Vpt = (void *)Vpt3[0][0];
        FileWriteData (Vpt, sizeof(float), SZ_float, fbin, 
                          output->stag_var[nv]);
      }
      FileClose(fbin, SZ_float);
/*
BOV_Header(output, filename);
*/
    }else{              /* -- multiple files -- */

      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;
        sprintf (filename, "%s/%s.%04d.%s", output->dir, output->var_name[nv], 
                                            output->nfile, output->ext);

        FileDelete (filename);  /* Avoid partial fill of pre-existing files */
        fbin = FileOpen (filename, SZ_float, "w");
        Vpt3 = Convert_dbl2flt(output->V[nv], units[nv],0);
        Vpt = (void *)Vpt3[0][0];
        FileWriteData (Vpt, sizeof(float), SZ_float, fbin, 
                          output->stag_var[nv]);
        FileClose (fbin, SZ_float);
      }
    }

  }else if (output->type == DBL_H5_OUTPUT || output->type == FLT_H5_OUTPUT){

  /* ------------------------------------------------------
       HDF5 (static grid) output (single/double precision)
     ------------------------------------------------------ */

    #ifdef USE_HDF5 
     single_file = YES;
     WriteHDF5 (output, grid);
    #else
     print ("! WriteData: HDF5 library not available\n");
     return;
    #endif

  }else if (output->type == VTK_OUTPUT) { 

  /* ------------------------------------------------------------------- */
  /*! - \b VTK output:  
      in order to enable parallel writing, files must be closed and
      opened again for scalars, since the distributed array descriptors 
      used by ArrayLib (Float_Vect) and (float) are different. This is
      done using the AL_Get_offset() and AL_Set_offset() functions.      */
  /* ------------------------------------------------------------------- */
    
    single_file = strcmp(output->mode,"single_file") == 0;
    sprintf (filename, "%s/data.%04d.%s", output->dir, output->nfile,
                                          output->ext);

    if (single_file){  /* -- single output file -- */

      FileDelete (filename);  /* Avoid partial fill of pre-existing files */
      fbin  = FileOpen(filename, SZ_Float_Vect, "w");
      WriteVTK_Header(fbin, grid);
      for (nv = 0; nv < output->nvar; nv++) {  /* -- write vectors -- */
        if (output->dump_var[nv] != VTK_VECTOR) continue;
        WriteVTK_Vector (fbin, output->V + nv, units[nv],
                         output->var_name[nv], grid);
      }

      #ifdef PARALLEL
       offset = AL_Get_offset(SZ_Float_Vect);
       FileClose(fbin, SZ_Float_Vect);
       fbin  = FileOpen(filename, SZ_float, "w");
       AL_Set_offset(SZ_float, offset);
      #endif
      
      for (nv = 0; nv < output->nvar; nv++) { /* -- write scalars -- */
        if (output->dump_var[nv] != YES) continue;
        WriteVTK_Scalar (fbin, output->V[nv], units[nv],
                         output->var_name[nv], grid);
      }
      FileClose(fbin, SZ_float);

    }else{          /* -- multiple output files -- */

      for (nv = 0; nv < output->nvar; nv++) { /* -- write vectors -- */
        if (output->dump_var[nv] != VTK_VECTOR) continue;
        if (strcmp(output->var_name[nv],"vx1") == 0) {
          sprintf (filename, "%s/vfield.%04d.%s", output->dir, output->nfile, 
                                                  output->ext);
        }else if (strcmp(output->var_name[nv],"bx1") == 0) {
          sprintf (filename, "%s/bfield.%04d.%s", output->dir, output->nfile, 
                                                  output->ext);
        }else{
          print ("! WriteData: unknown vector type in VTK output\n"); 
          QUIT_PLUTO(1);
        }

        FileDelete (filename);  /* Avoid partial fill of pre-existing files */
        fbin = FileOpen(filename, SZ_Float_Vect, "w");
        WriteVTK_Header(fbin, grid);
        WriteVTK_Vector(fbin, output->V + nv, units[nv],
                        output->var_name[nv], grid);
        FileClose(fbin, SZ_Float_Vect);
      }

      for (nv = 0; nv < output->nvar; nv++) {  /* -- write scalars -- */
        if (output->dump_var[nv] != YES) continue;
        sprintf (filename, "%s/%s.%04d.%s", output->dir, output->var_name[nv], 
                                            output->nfile,  output->ext);
        FileDelete (filename);  /* Avoid partial fill of pre-existing files */
        fbin = FileOpen(filename, SZ_Float_Vect, "w");
        WriteVTK_Header(fbin, grid);
      #ifdef PARALLEL
        offset = AL_Get_offset(SZ_Float_Vect);
        FileClose(fbin, SZ_Float_Vect);
        fbin  = FileOpen(filename, SZ_float, "w");
        AL_Set_offset(SZ_float, offset);
      #endif
        WriteVTK_Scalar(fbin, output->V[nv], units[nv],
                        output->var_name[nv], grid);
        FileClose (fbin, SZ_float);
      }
    }

  }else if (output->type == TAB_OUTPUT) { 

  /* ------------------------------------------------------
               Tabulated (ASCII) output
     ------------------------------------------------------ */

    single_file = YES;
    sprintf (filename,"%s/data.%04d.%s", output->dir, output->nfile,
                                         output->ext);
    WriteTabArray (output, filename, grid);

  }else if (output->type == PPM_OUTPUT) { 

  /* ------------------------------------------------------
                   PPM output
     ------------------------------------------------------ */

    single_file = NO;
    for (nv = 0; nv < output->nvar; nv++) {
      if (!output->dump_var[nv]) continue;
      sprintf (filename, "%s/%s.%04d.%s", output->dir, output->var_name[nv], 
                                          output->nfile, output->ext);
      WritePPM (output->V[nv], output->var_name[nv], filename, grid);
    }

  }else if (output->type == PNG_OUTPUT) { 

  /* ------------------------------------------------------
                   PNG  output
     ------------------------------------------------------ */

    #ifdef USE_PNG
     single_file = NO;
     for (nv = 0; nv < output->nvar; nv++) {
       if (!output->dump_var[nv]) continue;
       sprintf (filename, "%s/%s.%04d.%s", output->dir, output->var_name[nv], 
                                           output->nfile, output->ext);
       WritePNG (output->V[nv], output->var_name[nv], filename, grid);
     }
    #else
     print ("! PNG library not available\n");
     return;
    #endif
  }       

/* -------------------------------------------------------------
   3. Update corresponding ".out" file
   ------------------------------------------------------------- */

  sprintf (filename,"%s/%s.out",output->dir, output->ext);

  if (prank == 0) {
    if (output->nfile == 0) {
      fout = fopen (filename, "w");
    }else {
      fout = fopen (filename, "r+");
      for (nv = 0; nv < output->nfile; nv++) fgets (sline, 512, fout);
      fseek (fout, ftell(fout), SEEK_SET);
    }

  /* -- write a multi-column file -- */

    fprintf (fout, "%d %12.6e %12.6e %ld ",
             output->nfile, g_time, g_dt, g_stepNumber);

    if (single_file) fprintf (fout,"single_file ");
    else             fprintf (fout,"multiple_files ");

    if (IsLittleEndian()) fprintf (fout, "little ");
    else                  fprintf (fout, "big ");

    for (nv = 0; nv < output->nvar; nv++) { 
      if (output->dump_var[nv]) fprintf (fout, "%s ", output->var_name[nv]);
    }

    fprintf (fout,"\n");
    fclose (fout);
  }

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (prank == 0){
    time(&tend);
    print (" [%5.2f sec]",difftime(tend,tbeg));
  }
#endif
  print ("\n");

}

/* ********************************************************************* */
void GetCGSUnits (double *u)
/*!
 * Compute an array of c.g.s units
 * 
 * \param [in]  u   an array containing the c.g.s units of the primitive
 *                  variables, e.g., <tt> u[RHO] </tt> will be in gr/cm^3, 
 *                  <tt> u[VX1] </tt> will be in cm/s, etc...
 *
 *********************************************************************** */
{
  int nv;
  double unit_velocity = UNIT_VELOCITY;
  double unit_pressure = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
  double unit_mag      = UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY);
  
  #if PHYSICS == RHD || PHYSICS == RMHD
   unit_velocity = CONST_c;
  #endif

  u[RHO] = UNIT_DENSITY;
  EXPAND(u[VX1] = unit_velocity;  ,
         u[VX2] = unit_velocity;  ,
         u[VX3] = unit_velocity;)
  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(u[BX1] = unit_mag;  ,
          u[BX2] = unit_mag;  ,
          u[BX3] = unit_mag;)
  #endif
  #if HAVE_ENERGY
   u[PRS] = unit_pressure;
  #endif
}

void BOV_Header(Output *output, char *fdata)
{
  int  nv, nvar;
  char filename[256];
  FILE *fp;

  if (prank != 0) return;

/* -- count number of variables -- */

  nvar = 0;
  for (nv = 0; nv < output->nvar; nv++){
    if (!output->dump_var[nv]) continue;
    nvar++;
  }
  
  sprintf (filename,"data.%04d.bov",output->nfile);
  fp = fopen (filename, "w");
  fprintf (fp,"TIME: %f\n",g_time);
  fprintf (fp,"DATA_FILE: %s\n",fdata);
  fprintf (fp,"DATA_SIZE: %ld %ld %ld\n", NX1, NX2, NX3);
  fprintf (fp,"DATA_FORMAT: FLOAT\n");
  fprintf (fp,"VARIABLE: Density\n");
  fprintf (fp,"DATA_ENDIAN: LITTLE\n");
  fprintf (fp,"CENETRING: zonal\n");
  fprintf (fp,"BRICK_ORIGIN: %f %f %f\n",
           g_domBeg[IDIR], g_domBeg[JDIR], g_domBeg[KDIR]);
  fprintf (fp, "BRICK_SIZE: %f %f %f\n", 
           g_domEnd[IDIR]-g_domBeg[IDIR],
           g_domEnd[JDIR]-g_domBeg[JDIR],
           g_domEnd[KDIR]-g_domBeg[KDIR]);

  fprintf (fp,"DATA_COMPONENTS: %d\n", 2);
  
  fclose(fp);
}
