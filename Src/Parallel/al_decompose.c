/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ArrayLib functions to decompose the domain among the MPI processes. 

  ArrayLib functions to decompose the domain among the MPI processes. 

  \authors A. Malagoli (University of Chicago) 
  \authors G. Muscianisi (g.muscianisi@cineca.it)
  \authors A. Mignone (mignone@ph.unito.it)

  \date   Aug 24, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"  /*I "al_hidden.h" I*/

/* 
  If this is set, then each processor
  gathers the entire index set for a
  distributed array from all the group
*/
#define AL_GATHER_INDEXES

/* 
   The SZ structure stack is defined and maintained
   in al_szptr_.c
   Here we include an external reference to it in
   order to be able to make internal references to it.
*/
extern SZ *sz_stack[AL_MAX_ARRAYS];
extern int stack_ptr[AL_MAX_ARRAYS];

/* PROTOTYPES */
int AL_Find_decomp_(int sz_ptr, int mode, int *procs);
int AL_Global_to_local_(int sz_ptr);
int AL_Decomp1d_(int gdim, int lproc, int lloc, int *start, int *end);

/* ********************************************************************** */
int AL_Decompose(int sz_ptr, int *procs, int mode)
/*!  Create a distributed array descriptor and compile it
 *
 * \param [in]     sz_ptr integer pointer to the distributed array descriptor 
 * \param [in,out] procs  array with the processor decomposition 
 * \param [in]     mode   available value:\n 
                          AL_AUTO_DECOMP (internal decomposition) [Only for powers of two dimensions];\n
                          AL_MPI_DECOMP (MPI decomposition); \n
                          AL_USER_DECOMP (user defined)
 ************************************************************************ */
{
  int myrank, nproc;
  SZ *s;

  register int i, nd, nb;

  int ndim;
  int gdims[AL_MAX_DIM], ldims[AL_MAX_DIM], starts[AL_MAX_DIM];
  int count[AL_MAX_DIM], blocklen[AL_MAX_DIM], periods[AL_MAX_DIM];
  int reord, nleft, nright;
  int stagpt, type_size;
  int stride[AL_MAX_DIM];
  MPI_Datatype itype,ivector;
  MPI_Comm comm, cart_comm;

  int Find_decomp = AL_TRUE; /* By default, we use the internal decomposition */

  /* DIAGNOSTICS
    Check that sz_ptr points to an allocated SZ
  */
  if( stack_ptr[sz_ptr] == AL_STACK_FREE){
    printf("AL_Decompose: wrong SZ pointer\n");
  }

  s = sz_stack[sz_ptr];

  myrank = s->rank;
  nproc = s->size;

  ndim = s->ndim;

  /*
    Find if we need to find our own decomposition
  */

  if( mode != AL_USER_DECOMP && mode != AL_AUTO_DECOMP && mode != AL_MPI_DECOMP ) {
    printf("AL_Decompose: wrong mode. Using default AL_MPI_DECOMP\n");
    mode = AL_MPI_DECOMP;
  }

  if( mode == AL_USER_DECOMP ){ Find_decomp = AL_FALSE; }

  /* 
     If Find_decomp is TRUE, we proceed to find the
     decomposition
  */
  if( Find_decomp == AL_TRUE ){
    AL_Find_decomp_(sz_ptr, mode, procs);
  } else {
    for(nd=0;nd<ndim;nd++){
      s->lsize[nd] = procs[nd];
    }
  }
  
  /*
    Now we set the MPI datatypes for the ghost points
    exchange
  */
  comm = s->comm;
  
  /*
    First, we create the cartesian communicator
  */
  reord = 0;

  for(nd=0;nd<ndim;nd++){
    ldims[nd] = s->lsize[nd];
    periods[nd] = s->isperiodic[nd];
  }

#ifdef DEBUG
  if(sz_ptr==4){  /* -- print SZ_stagx -- */
     if(myrank==0) printf("MPI_Cart_create: %d x %d x %d, ndim %d, nproc %d\n",ldims[0],ldims[1],ldims[2],ndim,nproc);
  }
#endif

  MPI_Cart_create(comm, ndim, ldims, periods, reord, &cart_comm);
  
  s->cart_comm = cart_comm;

  /*
    Then we generate the cartesian topology 
    of the processors
  */
  MPI_Cart_coords(cart_comm, myrank, ndim, ldims);

  for(nd=0;nd<ndim;nd++){
    s->lrank[nd] = ldims[nd];
    MPI_Cart_shift(cart_comm, nd, 1, &nleft, &nright);
    s->left[nd] = nleft;
    s->right[nd] = nright;
  }

  /* 
     We also generate a set of 1D communicators that
     are useful when doing global operations along
     single directions.
  */
  {
    int rem_dims[AL_MAX_DIM];

    for(nd=0;nd<ndim;nd++){
      for(i=0;i<ndim;i++){ 
	if( i==nd ) { 
	  rem_dims[i] = AL_TRUE;
	} else {
	  rem_dims[i] = AL_FALSE;
	}
      }
      MPI_Cart_sub(cart_comm, rem_dims, &(s->oned_comm[nd]));
    }
  }
  
  /*
    Now we create the strided data types for the 
    exchange of the ghost points
  */

  /*
    Get the local array indexes (in global addressing mode)
  */
  AL_Global_to_local_(sz_ptr);

  /*
    Set the size of the local buffer for this array 
  */
  s->buffsize = 1;
  for(nd=0;nd<ndim;nd++){
    s->buffsize *=
      ( (s->end[nd]) - (s->beg[nd]) + (s->bg[nd]) + (s->eg[nd]) + 1 );
  }

#ifdef AL_GATHER_INDEXES

  /*
    If it has not yet been done, allocate the array of begs
    and ends, then gather the indexing information from all
    processors. We check that the array of 'begs' was not
    previously allocated, in order to avoid problems with
    duplicated descriptors.
  */
/*   if( (s->begs) )  free(s->begs); */

  if( !(s->begs = (int *)malloc(sizeof(int)*nproc*ndim*2))) {
    printf("AL_Decompose: could not allocate s->begs\n");
  }
  s->ends = s->begs + nproc*ndim;

  /* Gather the 'beg' index information from all nodes */
  MPI_Allgather(s->beg, ndim, MPI_INT, s->begs, ndim, MPI_INT, comm);

  /* Gather the 'end' index information from all nodes */
  MPI_Allgather(s->end, ndim, MPI_INT, s->ends, ndim, MPI_INT, comm);

#endif /* End #ifdef AL_GATHER_INDEXES */

  /* -- Set the type size of the array -- */

  MPI_Type_size(s->type, &(s->type_size));

  /*--------------------- Begin creation of strided data types -------*/

  /*
    Set the count, blocklen, and stride elements for
    the MPI_Vector calls.
    This is VERY, VERY cryptic, but it seems to work ..
  */

  /* Right->Left exchange */
  for(nd=0;nd<ndim;nd++){
    count[nd] = 1;
    blocklen[nd] = s->bg[nd];
    stride[nd] = s->larrdim_gp[nd]; 
    for(nb=0;nb<nd;nb++){
      blocklen[nd] *= s->larrdim_gp[nb];
      stride[nd]   *= s->larrdim_gp[nb];
    }
    for(nb=nd+1;nb<ndim;nb++){
      count[nd] *= s->larrdim_gp[nb];
    }
  }

  itype = s->type;
  for(nd=0;nd<ndim;nd++){
    if( count[nd] != 1){
      MPI_Type_vector( count[nd], blocklen[nd], stride[nd], itype, &ivector);
    } else {
      MPI_Type_contiguous( blocklen[nd], itype, &ivector);
    }
    MPI_Type_commit(&ivector);
    s->strided[nd] = ivector;
    s->type_rl[nd] = ivector;
  }

  /* Left->Right exchange */
  for(nd=0;nd<ndim;nd++){
    count[nd] = 1;
    /* If the array is staggered, we give 
       unique ownership of the overlapping
       point to the LEFT processor */
    if( s->isstaggered[nd] == AL_TRUE ){
      blocklen[nd] = s->bg[nd]+1;
    } else {
       blocklen[nd] = s->bg[nd];
    }

    stride[nd] = s->larrdim_gp[nd];
    for(nb=0;nb<nd;nb++){
      blocklen[nd] *= s->larrdim_gp[nb];
      stride[nd]   *= s->larrdim_gp[nb];
    }
    for(nb=nd+1;nb<ndim;nb++){
      count[nd] *= s->larrdim_gp[nb];
    }
  }

  itype = s->type;
  for(nd=0;nd<ndim;nd++){
    if( count[nd] != 1){
      MPI_Type_vector( count[nd], blocklen[nd], stride[nd], itype, &ivector);
    } else {
      MPI_Type_contiguous( blocklen[nd], itype, &ivector);
    }
    MPI_Type_commit(&ivector);
    s->type_lr[nd] = ivector;
  }

  /*--------------------- End creation of strided data types -------*/

  /* 
     Create the buffer pointers to the array
     for the SendRecv operations
  */
  MPI_Type_size( s->type, &type_size);

  for(nd=0;nd<ndim;nd++){
    s->tag1[nd] = nd*100;
    s->tag2[nd] = nd*100+1;
    /* This is the correction for a staggered mesh */
    if( s->isstaggered[nd] == AL_TRUE ){
      stagpt = 1;
    } else { 
      stagpt = 0; 
    }
    
    s->sendb1[nd] = s->lbeg[nd]             + stagpt;
    s->recvb1[nd] = s->lend[nd]+1;
    s->sendb2[nd] = s->lend[nd]-s->bg[nd]+1 - stagpt;
    s->recvb2[nd] = 0;

    for(nb=0;nb<nd;nb++){
      s->sendb1[nd] *= s->larrdim_gp[nb];
      s->recvb1[nd] *= s->larrdim_gp[nb];
      s->sendb2[nd] *= s->larrdim_gp[nb];
      s->recvb2[nd] *= s->larrdim_gp[nb];
    }

    /* 
       The buffer will be considered as char, so we
       have to correct for the size of the data  type
    */
    s->sendb1[nd] *= type_size;
    s->recvb1[nd] *= type_size;
    s->sendb2[nd] *= type_size;
    s->recvb2[nd] *= type_size;

  }

  /*
    Create the subarrays for the MPI-IO operations
  */
  {
    int istag;
 
    MPI_Datatype igsubarr, ilsubarr;

    MPI_Datatype igsubarr_stag[AL_MAX_DIM];
    MPI_Datatype ilsubarr_stag[AL_MAX_DIM];

    /* -----------------------------------------------------
        Create the global subarray for MPI_Set_file_view 
       ----------------------------------------------------- */

    for(nd=0;nd<ndim;nd++){
      gdims[nd]  = s->arrdim[nd];
      ldims[nd]  = s->larrdim[nd];
      starts[nd] = s->beg[nd] - s->bg[nd];
    }
#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
printf("%d, gsubarr: gdims[0:2] %d,%d,%d, ldims[0:2] %d,%d,%d, starts[0:2] %d,%d,%d\n", s->rank, gdims[0], gdims[1], gdims[2], ldims[0], ldims[1], ldims[2], starts[0], starts[1], starts[2]);
}
#endif

    AL_Type_create_subarray(ndim, gdims, ldims, starts,
                            AL_ORDER_FORTRAN, s->type, &igsubarr); 
    MPI_Type_commit(&igsubarr);
    s->gsubarr = igsubarr;

    /* -----------------------------------------------------
        Create the global staggered subarrays
        for MPI_Set_file_view 
       ----------------------------------------------------- */

    for (istag = 0; istag < ndim; istag++){
      for(nd = 0; nd < ndim; nd++){
        gdims[nd]  = s->arrdim[nd];
        ldims[nd]  = s->larrdim[nd];
        starts[nd] = s->beg[nd] - s->bg[nd];
      }
      /* --------------------------------------------------------------- */
      /*! <b> Bug fixed on Aug 26, 2012:</b>  
          When the global view of the file for a staggered
          variable is computed, the gdims[nd] has to be computed
          and then passed to the function AL_Type_create_subarray. 
          
          Since gdims[nd]=s->arrdim[nd], and in s->arrdim[nd] is 
          taking into account that the variable is staggered, we comment 
          the line "gdims[istag]++;" because it is no more needed. */
      /* --------------------------------------------------------------- */
        /* gdims[istag]++; */



      /* --------------------------------------------------------------- */
      /*! <b> Bug fixed on Aug 26, 2012: </b>
          The following "if" has been modified:

          OLD version:
          if (s->beg[istag] == s->bg[istag]) {
             ldims[istag]++; 
          }else{ 
             starts[istag]++; 
          }
  
          NEW version: 
          if (s->beg[istag] == s->bg[istag]) {
          }else{ 
             starts[istag]++; 
             ldims[istag]--; 
          }
  
          \li "ldims[istag]++;" has been cancelled from the firt part
              of the if for the same motivation explaned before; 
          \li "ldims[istag]--;" has been added in the second part of the
              if because in PLUTO the index of the staggered variables 
              start from -1, while in the ArrayLib they start from 0. */
      /* --------------------------------------------------------------- */

      if (s->beg[istag] == s->bg[istag]) {
      }else{ 
        starts[istag]++; 
        ldims[istag]--; 
      }

#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
    printf("%d, gsubarr_stag[%d]: gdims[0:2] %d,%d,%d, ldims[0:2] %d,%d,%d, starts[0:2] %d,%d,%d \n", s->rank, istag, gdims[0], gdims[1], gdims[2], ldims[0], ldims[1], ldims[2], starts[0], starts[1], starts[2]);
}
#endif

      AL_Type_create_subarray(ndim, gdims, ldims, starts,
                              AL_ORDER_FORTRAN, s->type, igsubarr_stag + istag); 
      MPI_Type_commit(igsubarr_stag + istag);
      s->gsubarr_stag[istag] = igsubarr_stag[istag];
    }

    /* -----------------------------------------------------
        Create the local subarray for the MPI_Set_write_all 
       ----------------------------------------------------- */

    for(nd=0;nd<ndim;nd++){
      gdims[nd]  = s->larrdim_gp[nd];
      ldims[nd]  = s->larrdim[nd];
      starts[nd] = s->lbeg[nd];
    }
#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
printf("%d, lsubarr: gdims[0:2] %d,%d,%d, ldims[0:2] %d,%d,%d, starts[0:2] %d,%d,%d \n", s->rank, gdims[0], gdims[1], gdims[2], ldims[0], ldims[1], ldims[2], starts[0], starts[1], starts[2]);
}
#endif
    AL_Type_create_subarray(ndim, gdims, ldims, starts,
			     AL_ORDER_FORTRAN, s->type, &ilsubarr);
    MPI_Type_commit(&ilsubarr);
    s->lsubarr = ilsubarr;


    /* -----------------------------------------------------
        Create the local staggered subarrays
        for MPI_Set_write_all
       ----------------------------------------------------- */

    for (istag = 0; istag < ndim; istag++){
      for(nd = 0; nd < ndim; nd++){
        gdims[nd]  = s->larrdim_gp[nd];
        ldims[nd]  = s->larrdim[nd];
        starts[nd] = s->lbeg[nd];
      }
      /* --------------------------------------------------------------- */
      /*! <b> Bugs fixed on Aug 26, 2012: </b>
          The following if has been modified: 

          OLD version: 
          if (s->beg[istag] == s->bg[istag]) {
             ldims[istag]++;  
             starts[istag]--; 
          }

          NEW version:
          if (s->beg[istag] == s->bg[istag]) {
          }else{
            ldims[istag]--; 
            starts[istag]++; 
          }

          \li "ldims[istag]++;" has been cancelled from the firt part of 
              the if because when the local subarray for the MPI_Set_write_all 
              for a staggered variable is computed, the ldims[nd] has to be 
              computed and then passed to the function AL_Type_create_subarray. 
              Since ldims[nd]=s->larrdim[nd], and in s->larrdim[nd] is 
              taking into account that the variable is staggered, we remouved 
              the line "ldims[istag]++;" because it is no more needed; 
          \li "starts[istag]--;" has been cancelled from the first part of 
              the if because in PLUTO the indexes for the staggered variables 
              start locally from -1, while in ArrayLib they start from 0; 
          \li the "else" has been added to take into account the motivation
              explaned in the point before. */
      /* --------------------------------------------------------------- */

      if (s->beg[istag] == s->bg[istag]) {
      }else{
        ldims[istag]--; 
        starts[istag]++; 
      }


#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
printf("%d, lsubarr_stag[%d]: gdims[0:2] %d,%d,%d, ldims[0:2] %d,%d,%d, starts[0:2] %d,%d,%d \n", s->rank, istag, gdims[0], gdims[1], gdims[2], ldims[0], ldims[1], ldims[2], starts[0], starts[1], starts[2]);
}
#endif
      AL_Type_create_subarray(ndim, gdims, ldims, starts,
                              AL_ORDER_FORTRAN, s->type, ilsubarr_stag + istag); 
      MPI_Type_commit(ilsubarr_stag + istag);
      s->lsubarr_stag[istag] = ilsubarr_stag[istag];
    }

  }

  /* 
     Now array is officially compiled
  */
  s->compiled = AL_TRUE;

  /* DIAGNOSTICS */
#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
  printf("AL_Decompose: Decomposition successful\n");
  for(nd=0;nd<ndim;nd++){
    printf("AL_Decompose: %d %d\n", s->beg[nd], s->end[nd]);
  }
}
#endif

  return (int) AL_SUCCESS;
}

/* -------------- INTERNAL ROUTINES ------------------*/

/* ********************************************************************* */
int AL_Find_decomp_(int sz_ptr, int mode, int *procs)
/*!
  Find the decomposition for a distributed array given its 
  global size and the number of nodes
 *
 * \param [in]    sz_ptr integer pointer to the SZ structure
 * \param [in]    mode   mode parameter: AL_AUTO_DECOMP, AL_MPI_DECOMP, AL_USER_DECOMP
 * \param [in,out] procs  array of processor decomposition
 *********************************************************************** */
{
  register int i, nd;

  SZ *s;

  int myid, nproc, ntdim, npdim;
  int ldims[AL_MAX_DIM], ipdims[AL_MAX_DIM];
  int s_inds[AL_MAX_DIM], t_ldims[AL_MAX_DIM];
  int gdims[AL_MAX_DIM];

  s = sz_stack[sz_ptr];

  myid = s->rank;
  nproc = s->size;

  ntdim = s->ndim;

#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
    if(myid==0) printf("AL_Find_decomp_: Using Mode - %d\n",mode);
}
#endif

  /*
    If mode is AL_USER_DECOMP, we get the user defined decomposition

  */
  if( mode == AL_USER_DECOMP ){
#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
    if(myid==0) printf("AL_Find_decomp_: Using User Decomp Mode - %d\n",ntdim);
}
#endif
    for(nd=0;nd<ntdim;nd++){
      if( ldims[nd] != 0 ){
	s->lsize[nd] = ldims[nd];
      } else {
	s->lsize[nd] = 1; /* This is a safety patch, a more
                             sophisticated tratement will be needed */
      }
    }
    return 0;
  }

  /*
    Determine the parallel directions and
    save them in ipdims[]
  */
  npdim = 0;
  for(nd=0;nd<ntdim;nd++){
    s->lsize[nd] = 1;
    if( s->isparallel[nd] == AL_TRUE ){
      ipdims[npdim] = nd;
      /* Correct for the overlapping points in a staggered mesh */
      if( s->isstaggered[nd] == AL_TRUE ) {
	gdims[npdim] = s->arrdim[nd]-AL_STAGGERED_OVERLAP;
      } else {
	gdims[npdim] = s->arrdim[nd];
      }
      npdim = npdim+1;
    } else {
      /* The dimension is NOT parallel */
      procs[nd] = 1;
    }
  }

  /* 
     Introduce special test for power of two 
     dimensions and number of processors.

     We check whether both the size of the communicator
     AND each of the dimensions are powers of two.
     If this is the case, then the AL_AUTO_DECOMP is used.
  */
  {
    int isp;
    isp = 0;
    if(s->ndim <4){
      if( AL_POWEROF2(s->size) ){
	for(nd=0;nd<s->ndim;nd++){
	  if(AL_POWEROF2(s->arrdim[nd])) isp++;
	}
      }
      if(isp==s->ndim) {
	mode = AL_AUTO_DECOMP;
#ifdef DEBUG
      if(myid==0) printf("Using mode AL_AUTO_DECOMP\n");
#endif
      }
    }
  }

  /* 
     If mode is AL_AUTO_DECOMP, AL finds the decomposition 
     We try to find a "maximally cubic" decomposition, which
     is useful for the application of the multigrid algorithm.
  */
  if( mode == AL_AUTO_DECOMP ){
    if( AL_Auto_Decomp_(nproc,npdim,ldims,gdims) == AL_FAILURE) {
      mode = AL_MPI_DECOMP;
#ifdef DEBUG
      if( myid == 0 ) printf("Auto Decomp failed: using MPI Decomp Mode\n");
#endif
    }
  }

  /* 
    If mode is AL_MPI_DECOMP, MPI finds the decomposition
  */
#ifdef DEBUG
  if(myid==0) printf("Auto Decomp: Mode is: %d %d %d %d\n",mode, AL_MPI_DECOMP, nproc, npdim);
#endif

  if( mode == AL_MPI_DECOMP ){
    for(i=0;i<npdim;i++){ ldims[i]=0 ;}
    MPI_Dims_create(nproc, npdim, ldims);
    /* This may not get released */
    /*     AL_Dims_create(nproc, npdim, ldims); */
    {
    /* Make an attempt at aligning the decomposition with the
       array dimensions. We try to match the longest dimensions
       of the decomposition with those of the array */
      AL_Sort_(npdim,ldims,s_inds);
      for(i=0;i<npdim;i++) t_ldims[i] = ldims[s_inds[i]];
      AL_Sort_(npdim,gdims,s_inds);
      for(i=0;i<npdim;i++) ldims[s_inds[i]] = t_ldims[i];
    }
#ifdef DEBUG
    if(myid==0) printf("MPI Decomp: %d %d %d %d %d\n", ldims[0],ldims[1],ldims[2],nproc,npdim);
#endif
  }

  /*
     Fill the decomposition information back
     int the SZ structure
  */
#ifdef DEBUG
  if(myid==0) printf("End of AL_Fing_decomp_: ");
#endif
  for(i=0;i<npdim;i++){ 
    s->lsize[ipdims[i]]=ldims[i];
    procs[ipdims[i]] = ldims[i];
#ifdef DEBUG
    if(myid==0) printf("%d ",s->lsize[ipdims[i]]);
#endif
  }
#ifdef DEBUG 
  if(myid==0) printf("\n");
#endif 
  return 0;
}


/* ********************************************************************* */
int AL_Global_to_local_(int sz_ptr)
/*!
 * Compute addresses of local portions of array from 
 * global dimensions and processor decomposition
 *
 * \param [in] sz_ptr integer pointer to SZ structure
 *********************************************************************** */
{
  register int i, j;
  SZ *s;

  int myrank, nproc, ndim;
  int gdim, lproc, lloc;
  int start, end;

  s = sz_stack[sz_ptr];

  myrank = s->rank;
  nproc = s->size;
  ndim = s->ndim;

  for(i=0;i<ndim;i++){
    gdim = s->arrdim[i];
    lproc = s->lsize[i];
    lloc = s->lrank[i];

    /* We apply the following trick if the array is staggered */
    if( s->isstaggered[i] == AL_TRUE ){ gdim = gdim-1;}

    AL_Decomp1d_(gdim, lproc, lloc, &start, &end);

    s->beg[i] = start+s->bg[i];
    s->end[i] = end+s->bg[i];
    if( s->isstaggered[i] == AL_TRUE ){ s->end[i] = end+1+s->bg[i];}
    s->lbeg[i] = s->bg[i];
    /* -------------------------------------------------------------------- */
    /*! <b> Bug fixed on Dec 7, 2011: </b>
        The lines: 

        if( s->isstaggered[i] == AL_TRUE ){ s->lend[i] = (end-start+1+s->lbeg[i]);}

        if( s->isstaggered[i] == AL_TRUE ){ s->larrdim_gp[i] = (end-start+1+s->bg[i]+s->eg[i]+1);}

        if( s->isstaggered[i] == AL_TRUE ){ s->larrdim[i] = (end-start+1+1);}
      
        have been added to manage in the right way the correspondence among 
        global and local indeces for staggered variables. */
    /* -------------------------------------------------------------------- */
    s->lend[i] = (end-start+s->lbeg[i]);
    if( s->isstaggered[i] == AL_TRUE ){ s->lend[i] = (end-start+1+s->lbeg[i]);}
    s->larrdim_gp[i] = (end-start+s->bg[i]+s->eg[i]+1);
    if( s->isstaggered[i] == AL_TRUE ){ s->larrdim_gp[i] = (end-start+1+s->bg[i]+s->eg[i]+1);}
#ifdef DEBUG
  if(sz_ptr==4){  /* print SZ_stagx */
    printf("%d AL_Global_to_local, dim %d : s->beg %d, s->end %d, lbeg %d, lend %d, bg %d, larrdim_gp %d\n",s->rank,i, s->beg[i],s->end[i], s->lbeg[i], s->lend[i], s->bg[i], s->larrdim_gp[i]);
}
#endif
    s->larrdim[i] = end-start+1;
    if( s->isstaggered[i] == AL_TRUE ){ s->larrdim[i] = (end-start+1+1);}
	      
  }
  
  /* 
     These are the offsets for the DO loop index computations.

     We will gradually replace "offset" with "stride", so that
     the definition is more in line with the MPI strided data
     types definitions.
  */
  for(i=0;i<ndim;i++){
    s->offset[i] = 1; 
    s->stride[i] = 1; 
    for(j=0;j<=i-1;j++){
      s->offset[i] *= s->larrdim_gp[j];  
      s->stride[i] *= s->larrdim_gp[j];  
    }
#ifdef DEBUG
  if(sz_ptr==4){  /* print solo SZ_stagx */
    printf("%d Strides: %d %d\n",s->rank,i,s->stride[i]);
}
#endif
  }

  return (int) AL_SUCCESS;
}

/* ********************************************************************** */
/* --  NOTE: the following routine, AL_Decomp1d_, has been
             modified from the MPE_Decomp1d_ routine, which
             is part of the MPICH distribution.       -- */

int AL_Decomp1d_(int gdim, int lproc, int lloc, int *start, int *end)
/*!
 * Decompose a 1D array, given the number of processors along the direction
 *
 * \param [in]  gdim  integer size of the global dimension
 * \param [in]  lproc integer size of the number of processors along the dimension
 * \param [in]  lloc  integer location of this node along the dimension
 * \param [out] start integer pointer to start address for the array (C-convention)
 * \param [out] end   integer pointer to end address for the array (C-convention)
  *********************************************************************** */
{
  int nlocal, deficit, itmp;

  nlocal = gdim/lproc;
  *start = lloc * nlocal;
  deficit = gdim % lproc;
  itmp = lloc < deficit ? lloc : deficit;
  *start = *start + itmp;
  if( lloc < deficit ){ nlocal = nlocal+1; }
  *end = *start + nlocal -1;

#ifdef DEBUG
  printf("AL_Decomp1_: %d %d %d %d - %d %d  \n", lloc, deficit, *start, *end, nlocal, lproc);
#endif
  
  if( *end >= gdim || lloc == lproc-1 ) {*end = gdim-1;}
  
  return (int) AL_SUCCESS;
}

