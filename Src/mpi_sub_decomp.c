#include "pluto.h"


/* real sb_A = -0.75;
real sb_vy, sb_Lx, sb_Ly; */

double sb_Omega = 1.0;
double sb_q    = 1.5;

/* ************************************************************** */
void Init (double *us, double x1, double x2, double x3)
/* 
 *
 * 
 *
 *
 **************************************************************** */
{
  static int first_call = 1;
  double b0, rnd;
  double x, y, z;
  int  zeroflux;	

  #ifndef SHEARINGBOX
   print1 ("! ShearingBox module has not been included.\n");
   print1 ("! Cannot continue.\n");
   QUIT_PLUTO(1);
  #endif
  
	zeroflux = (int) g_inputParam[ZEROFLUX];

  x = x1; y = x2; z = x3;

/* --------------------------------------
    seeds the random number differently
    according to processor.
   -------------------------------------- */

  if (first_call == 1){
    srand(time(NULL) + prank);
    first_call = 0;
  }
  rnd = (double)(rand())/((double)RAND_MAX + 1.0);

  us[RHO] = exp(-0.5*z*z/g_inputParam[TBOUND]);
  us[VX] = 0.0;
  us[VY] = 2.0*sb_A*x + 0.01*rnd;
  us[VZ] = 0.0;
  us[PRS] = g_inputParam[TBOUND]*us[DN];
   
  us[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD
   b0     = sqrt(2.0/g_inputParam[BETA]);
   us[BX] = 0.;
   us[BY] = 0.;
	
   if (zeroflux == 1) 
	   us[BZ] = b0*sin(2.*CONST_PI*x);
   else
	   us[BZ] = b0;
   


   us[AX] = 0.;
	
  if (zeroflux == 1)
   us[AY] = -b0*cos(2.*CONST_PI*x)/2./CONST_PI;
  else
   us[AY] = b0*x;
	
   us[AZ] = 0.0;
  #endif
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
#define aMAXWELL   0
#define aREYNOLDS  1
#define aRHO       2
#define aT1        3
#define aT4        4
#define aEM        5
#define aKN        6
#define aBX        7
#define aBY        8
#define aBZ        9
#define aFC       10
#define aFT       11
#define aEM_FLUCT 12

#define NAVERAGES 13
{
  int n, i, j, k, k1, kk, nv, nx, ny, nz, ierr, nxl, nyl, nzl;
        int nprx, npry, nprz, xp, yp, zp, midx, midy, midz;
	int midxl, midyl, midzl;
	int surface;
	int dims[3] = {1,1,0};
	static int gs[2];
	static int str[2], ldm[2]; 
	static int myrank_y, nptxy, nptyz, nptxz;
	static int first_call=1;
	static int nfile1 = 0;
	double rho, vx, vy, vz, prs, Bx, By, Bz;
	double T, gradT;
	double *dx, *dy, *dz;
	double *x, *y, *z;
	static double **av, **avg, **recv_buf;
	static float **section2dxy, **section2dyz, **section2dxz;
	static float **avg2dy;
	static MPI_Comm XYcomm, XZcomm, YZcomm, Xcomm, Ycomm, Zcomm;
	MPI_Comm cart_comm;
	static MPI_Datatype subarr2dxy, subarr2dyz, subarr2dxz;
	MPI_File ifp, ifp1;
	MPI_Status status;
	MPI_Offset offset;
	char fname[128];
    
	
	/* -- get total number of points in x and y (global domain) -- */
	nx = grid[IDIR].np_int_glob;
	ny = grid[JDIR].np_int_glob;
	nz = grid[KDIR].np_int_glob;
	nxl = grid[IDIR].np_int;
	nyl = grid[JDIR].np_int;
	nzl = grid[KDIR].np_int;
	nprx = grid[IDIR].nproc;
	npry = grid[JDIR].nproc;
	nprz = grid[KDIR].nproc;
	xp   = grid[IDIR].rank_coord;
	yp   = grid[JDIR].rank_coord;
	zp   = grid[KDIR].rank_coord;
	midx = nprx/2;
	midy = npry/2;
	midz = nprz/2;
	midxl = (IBEG+IEND)/2 - IEND;
	midyl = (JBEG+JEND)/2 - JBEG;
	midzl = (KBEG+KEND)/2 - KBEG;
	
	
	
	/* -- allocate memory and set pointer shortcuts -- */
	
	if (av == NULL){
		av       = ARRAY_2D(NX3, NAVERAGES, double);
		avg      = ARRAY_2D(NX3, NAVERAGES, double);
		recv_buf = ARRAY_2D(nz,  NAVERAGES, double);
		section2dxy = ARRAY_2D(nyl,nxl,float);
		section2dyz = ARRAY_2D(nzl,nyl,float);
		section2dxz = ARRAY_2D(nzl,nxl,float);
		avg2dy      = ARRAY_2D(nzl, nxl, float);
		for (k = 0; k < nz; k++){			 
			for (n = 0; n < NAVERAGES; n++){
				recv_buf[k][n] = -1.0;
			}}
	}
	
	dx = grid[IDIR].dx; dy = grid[JDIR].dx; dz = grid[KDIR].dx;
	x = grid[IDIR].x;   y = grid[JDIR].x;   z = grid[KDIR].x;
	
#ifdef PARALLEL	
	if (first_call == 1) {
		
		JDOM_LOOP(j) IDOM_LOOP(i) section2dxy[j-JBEG][i-IBEG] = 0.;
		KDOM_LOOP(k) IDOM_LOOP(i) section2dxz[k-KBEG][i-IBEG] = 0.;
		KDOM_LOOP(k) JDOM_LOOP(i) section2dyz[k-KBEG][j-JBEG] = 0.;
		
			
		
 	
		/* -- get different sub-communicators corresponding
		 to 2D horizontal  planes. -- */
		AL_Get_cart_comm(SZ, &cart_comm);
		k = MPI_Cart_sub(cart_comm, dims, &XYcomm);
		dims[0] = 1; dims[1] = 0; dims[2] = 1;
		k = MPI_Cart_sub(cart_comm, dims, &XZcomm);
		dims[0] = 0; dims[1] = 1; dims[2] = 1;
		k = MPI_Cart_sub(cart_comm, dims, &YZcomm);
		/* -- create sub-communicators corresponding to vertical columns -- */
		
		dims[0] = 1; dims[0] = 0; dims[2] = 0;
		MPI_Cart_sub(cart_comm, dims, &Xcomm);
		dims[0] = 0; dims[0] = 1; dims[2] = 0;
		MPI_Cart_sub(cart_comm, dims, &Ycomm);		
		dims[0] = 0; dims[1] = 0; dims[2] = 1;
		MPI_Cart_sub(cart_comm, dims, &Zcomm);
		
		MPI_Comm_rank(Ycomm, &myrank_y);


		gs[0]  = nx;
		gs[1]  = ny;
		str[0] = grid[IDIR].rank_coord * nxl;
		str[1] = grid[JDIR].rank_coord * nyl;
		ldm[0]=nxl;
		ldm[1]=nyl;                     
		nptxy = nxl*nyl;
		MPI_Type_create_subarray(2,gs,ldm,str,MPI_ORDER_FORTRAN, MPI_FLOAT, &subarr2dxy);
		MPI_Type_commit(&subarr2dxy);
		
		
		gs[0] = ny;
		gs[1] = nz;
		str[0] = grid[JDIR].rank_coord * nyl;
		str[1] = grid[KDIR].rank_coord * nzl;
		ldm[0]=nyl;
		ldm[1]=nzl;                     
		nptyz = nzl*nyl;
		MPI_Type_create_subarray(2,gs,ldm,str,MPI_ORDER_FORTRAN, MPI_FLOAT, &subarr2dyz);
		MPI_Type_commit(&subarr2dyz);
		
		gs[0] = nx;
		gs[1] = nz;
		str[0] = grid[IDIR].rank_coord * nxl;
		str[1] = grid[KDIR].rank_coord * nzl;
		ldm[0]=nxl;
		ldm[1]=nzl;                     
		nptxz = nzl*nxl;
		MPI_Type_create_subarray(2,gs,ldm,str,MPI_ORDER_FORTRAN, MPI_FLOAT, &subarr2dxz);
		MPI_Type_commit(&subarr2dxz);		
		
		first_call = 0;	
	}
	
	/* Compute y averages */
	KDOM_LOOP(k) IDOM_LOOP(i) {
		JDOM_LOOP(j) section2dxz[k-KBEG][i-IBEG] += (float) d->Vc[BX2][k][j][i];
	
	}
	MPI_Reduce(section2dxz[0], avg2dy[0], nptxz, MPI_FLOAT, MPI_SUM, 0, Ycomm);
	if (myrank_y == 0) {
		printf(" %i   %i\n",grid[KDIR].rank_coord, grid[IDIR].rank_coord);
		/* sprintf (fname,"avgy.%05d.dat",nfile1);
		ierr = MPI_File_open(XZcomm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR | 
							 MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &ifp);*/		
		T = 1./(float)ny;
		KDOM_LOOP(k) IDOM_LOOP(i) {
			avg2dy[k-KBEG][i-IBEG] /= T;
		}
		//ierr  = MPI_File_set_view(ifp, offset, MPI_FLOAT, subarr2dxz,"native", MPI_INFO_NULL); 
		//ierr  = MPI_File_write_all(ifp, avg2dy[0], nptxz, MPI_FLOAT, &status);
		//MPI_File_close(&ifp);

	}
	
      
	/* Write section at the bottom boundary */
	 if (grid[KDIR].rank_coord == 0) {
		offset = 0; 
		sprintf (fname,"2d_x3bg.%05d.dat",nfile1);
		ierr = MPI_File_open(XYcomm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR | 
							 MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &ifp);
		for (nv=0; nv < NVAR; nv++) {
		  JDOM_LOOP(j) IDOM_LOOP(i) { 
					section2dxy[j-JBEG][i-IBEG] = (float) d->Vc[nv][KBEG][j][i];		      
		    }
		    ierr  = MPI_File_set_view(ifp, offset, MPI_FLOAT, subarr2dxy,"native", MPI_INFO_NULL); 
			ierr  = MPI_File_write_all(ifp, section2dxy[0], nptxy, MPI_FLOAT, &status);
			offset = offset + (long)(long)nx*(long)(long)ny *(long)(long)4.;
		}
				ierr  = MPI_File_close(&ifp);

				}   
	  if (zp == midz ) {
		  offset = 0; 
		  sprintf (fname,"2d_x3mid.%05d.dat",nfile1);
		  ierr = MPI_File_open(XYcomm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR | 
							   MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &ifp);
		  for (nv=0; nv < NVAR; nv++) {
			  JDOM_LOOP(j) IDOM_LOOP(i) { 
				  section2dxy[j-JBEG][i-IBEG] = (float) d->Vc[nv][midzl][j][i];		      
			  }
			  ierr  = MPI_File_set_view(ifp, offset, MPI_FLOAT, subarr2dxy,"native", MPI_INFO_NULL); 
			  ierr  = MPI_File_write_all(ifp, section2dxy[0], nptxy, MPI_FLOAT, &status);
			  offset = offset + (long)(long)nx*(long)(long)ny *(long)(long)4.;
		  }
		  ierr  = MPI_File_close(&ifp);
		  		  
	  }
	if (xp == midx ) {
		offset = 0; 
		sprintf (fname,"2d_x1mid.%05d.dat",nfile1);
		ierr = MPI_File_open(YZcomm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR | 
							 MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &ifp);
		for (nv=0; nv < NVAR; nv++) {
			KDOM_LOOP(k) JDOM_LOOP(j) { 
				section2dyz[k-KBEG][j-JBEG] = (float) d->Vc[nv][k][j][midxl];		      
			}
			ierr  = MPI_File_set_view(ifp, offset, MPI_FLOAT, subarr2dyz,"native", MPI_INFO_NULL); 
			ierr  = MPI_File_write_all(ifp, section2dyz[0], nptyz, MPI_FLOAT, &status);
			offset = offset + (long)(long)ny*(long)(long)nz *(long)(long)4.;
		}
		ierr  = MPI_File_close(&ifp);
		
	}
	if (yp == midy ) {
		offset = 0; 
		sprintf (fname,"2d_x2mid.%05d.dat",nfile1);
		ierr = MPI_File_open(XZcomm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR | 
							 MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &ifp);
		for (nv=0; nv < NVAR; nv++) {
			KDOM_LOOP(k) IDOM_LOOP(i) { 
				section2dxz[k-KBEG][i-IBEG] = (float) d->Vc[nv][k][midyl][i];		      
			}
			ierr  = MPI_File_set_view(ifp, offset, MPI_FLOAT, subarr2dxz,"native", MPI_INFO_NULL); 
			ierr  = MPI_File_write_all(ifp, section2dxz[0], nptxz, MPI_FLOAT, &status);
			offset = offset + (long)(long)nx*(long)(long)nz *(long)(long)4.;
		}
		ierr  = MPI_File_close(&ifp);
		
	} 
//	MPI_Barrier(cart_comm);
         nfile1++;                     
		
#endif	
	/* -- Loop on k index and get averages in xy planes -- */
	
	Boundary (d, ALL_DIR, grid);
	KDOM_LOOP(k){
		k1 = k - KBEG;
		for (n = 0; n < NAVERAGES; n++) av[k1][n] = 0.0;
		JDOM_LOOP(j) IDOM_LOOP(i){
			rho = d->Vc[RHO][k][j][i];
			
			vx  = d->Vc[VX1][k][j][i];
			vy  = d->Vc[VX2][k][j][i] - 2.0*sb_A*x[i];
			vz  = d->Vc[VX3][k][j][i];
			
			Bx  = d->Vc[BX1][k][j][i];
			By  = d->Vc[BX2][k][j][i];
			Bz  = d->Vc[BX3][k][j][i];
			
			prs = d->Vc[PRS][k][j][i];
			
			av[k1][aMAXWELL]  += Bx*By;
	        av[k1][aREYNOLDS] += rho*vx*vy;
			
			T        = prs/rho;
	        av[k1][aT1]  += T;
			av[k1][aT4]  += T*T*T*T;
			av[k1][aRHO] += rho;
			
			av[k1][aEM]  += 0.5*(Bx*Bx + By*By + Bz*Bz);
			
			av[k1][aKN]  += 0.5*rho*(vx*vx + vy*vy + vz*vz);
			
			av[k1][aBX] += Bx;
			av[k1][aBY] += By;
			av[k1][aBZ] += Bz;
			
			gradT   = d->Vc[PRS][k+1][j][i]/d->Vc[RHO][k+1][j][i];
			gradT  -= d->Vc[PRS][k-1][j][i]/d->Vc[RHO][k-1][j][i];
			gradT  /= 2.0*dz[k];      
			av[k1][aFC] += -g_inputParam[KAPPA]*rho*gradT;  /* conductive flux */            
		}    
	}
	
#ifdef PARALLEL
	MPI_Allreduce (av[0], avg[0], NX3*NAVERAGES, MPI_DOUBLE, MPI_SUM, XYcomm);
	MPI_Barrier (MPI_COMM_WORLD);
	
#endif
    
	T = 1.0/((double)nx*(double)ny);
	KDOM_LOOP(k){
	    k1 = k - KBEG;
	    JDOM_LOOP(j) IDOM_LOOP(i) {
			av[k1][aFT] += 2.5*d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*(d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]-avg[k1][aT1]*T) ;
			av[k1][aEM_FLUCT] += 0.5*((d->Vc[BX1][k][j][i] - avg[k1][aBX]*T)*(d->Vc[BX1][k][j][i] - avg[k1][aBX]*T)               + (d->Vc[BX2][k][j][i] - avg[k1][aBY]*T)*(d->Vc[BX2][k][j][i] - avg[k1][aBY]*T) +
									  (d->Vc[BX3][k][j][i] - avg[k1][aBZ]*T)*(d->Vc[BX3][k][j][i] - avg[k1][aBZ]*T));
			
		}
		//	printf("%12.6e\n", av[k1][aFT]);	
	}
	
#ifdef PARALLEL
	MPI_Allreduce (av[0], avg[0], NX3*NAVERAGES, MPI_DOUBLE, MPI_SUM, XYcomm);
	MPI_Barrier (MPI_COMM_WORLD);
	
#endif	
	
	
#ifdef PARALLEL	
	
	
	MPI_Gather (avg[0], NX3*NAVERAGES, MPI_DOUBLE, recv_buf[0], NX3*NAVERAGES, 
				MPI_DOUBLE, 0, Zcomm);
	
#endif
	
	/* --------------------------------------
	 write data to disk
	 -------------------------------------- */
	
	if (prank == 0){
		static int nfile = -1;
		char fname[64];
		FILE *fp;
		
		
		/* -- set file number:
		 if we start from the beginning, nfile = 0;
		 if we restart, scan all the average.nnnnn.dat files
		 until the current time level is reached. Use 0 otherwise -- */
		
		if (nfile == -1){  /* -- nfile is not defined -- */
			if (g_stepNumber == 0) nfile = 0;
			else {
				for (nfile = 0; nfile < 64000; nfile++){  /* -- read at most 64000 files -- */
					sprintf (fname,"average.%05d.dat",nfile);
					fp = fopen(fname,"r");
					if (fp == NULL){
						print1 ("! Analysis: could not find file number. Setting nfile=0\n");
						nfile = 0;
						break;
					}
					fread (&T, sizeof(double), 1, fp);
					printf (">> Analysis: scanning %s: t = %12.6e\n", fname, T);
					if (T >= g_time) break;
					fclose(fp);
				}
			}
		}
		
		sprintf (fname,"average.%05d.dat",nfile);
		printf (">> Analysis: writing file %s\n",fname);
		
		fp = fopen(fname,"wb");
		
		/* -- write time and grid size -- */
		
		fwrite (&g_time, sizeof(double), 1, fp);
		
		T = (double)NAVERAGES;
		fwrite (&T, sizeof(double), 1, fp);
		
		T = (double)nz;
		fwrite (&T, sizeof(double), 1, fp);
		
		/* -- normalize data -- */
		
		T = 1.0/((double)nx*(double)ny);
		for (k = 0; k < nz; k++) for (n = 0; n < NAVERAGES; n++ ) {
			recv_buf[k][n] *= T;
		}
		
		
		fwrite (recv_buf[0], sizeof(double), nz*NAVERAGES, fp);
		fclose(fp);
		nfile++;
	}
	
}


/* ************************************************************** */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *
 **************************************************************** */
{
  int   i, j, k, nv;
  int fixedT;	
  real  *x1, *x2, *x3;
  real ***q;
  real Tkbeg, Tkend, T0, Tk1, Tks, dTdz, dpdz, alp;
	real zb, dzb, Sigma, kcond;


	fixedT = (int) g_inputParam[FIXEDT];
	T0 = g_inputParam[TBOUND];
	Sigma = g_inputParam[SIGMA];
  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < 1.e-5) d->Vc[RHO][k][j][i] = 1.e-5;
    }
  }
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    X1_BEG_LOOP(k,j,i){}
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    X1_END_LOOP(k,j,i){}
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    X2_BEG_LOOP(k,j,i){}
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    X2_END_LOOP(k,j,i){}
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
	          
	  if (box->vpos == CENTER){
		  zb = grid[KDIR].xl[KBEG];
		  dzb = grid[KDIR].dx[KBEG];
         BOX_LOOP(box,k,j,i){
      for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = d->Vc[nv][2*KBEG-k-1][j][i];
		if (fixedT == 1) { 
			alp = (1.+zb*dzb/2./T0)/(1.-zb*dzb/2./T0);
			 Tk1 = T0;
			 }
		 else {
			 kcond = g_inputParam[KAPPA]*d->Vc[RHO][KBEG][j][i];
			 Tkbeg = d->Vc[PRS][KBEG][j][i]/d->Vc[RHO][KBEG][j][i];
			 dTdz    = Sigma*pow(Tkbeg, 3.)*dzb/2./kcond;
			 dTdz    = dTdz/(1.+4.*dTdz);
			 Tk1   = Tkbeg * (1. - 2.*dTdz);
			 Tks   = Tkbeg * (1. - dTdz);
			 alp = (1.+zb*dzb/2./Tks)/(1.-zb*dzb/2./Tks); 
		 }
         d->Vc[VZ][k][j][i] *= -1.0;
         d->Vc[BX][k][j][i] *= -1.0;
         d->Vc[BY][k][j][i] *= -1.0;
	     d->Vc[PRS][k][j][i] = d->Vc[PRS][KBEG][j][i]*alp;
		 d->Vc[RHO][k][j][i] = d->Vc[PRS][k][j][i]/Tk1;
		 }
	  }else if (box->vpos == X1FACE) {
        #ifdef STAGGERED_MHD
		  BOX_LOOP(box, k, j, i){
		  d->Vs[BXs][k][j][i] = -d->Vs[BXs][2*KBEG-k-1][j][i];
		  }
        #endif  
	  } else if (box->vpos == X2FACE) {
           #ifdef STAGGERED_MHD
		  BOX_LOOP(box, k, j, i){  
		  d->Vs[BYs][k][j][i] = -d->Vs[BYs][2*KBEG-k-1][j][i];
		  }
         #endif
	  }
    }

 
	
	if (side == X3_END) {  /* -- X3_END boundary -- */
		
		if (box->vpos == CENTER){
			zb = grid[KDIR].xl[KEND];
			dzb = grid[KDIR].dx[KEND];
			BOX_LOOP(box,k,j,i){
				for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = d->Vc[nv][2*KEND-k+1][j][i];
				if (fixedT == 1) {
					alp = (1.-zb*dzb/2./T0)/(1.+zb*dzb/2./T0);
					Tk1 = T0;
				}
				else {
					kcond = g_inputParam[KAPPA]*d->Vc[RHO][KEND][j][i];
					Tkend = d->Vc[PRS][KEND][j][i]/d->Vc[RHO][KEND][j][i];
					dTdz  = Sigma*pow(Tkend, 3.)*dzb/2./kcond;
					dTdz  = dTdz/(1.+4.*dTdz);
					Tk1   = Tkend * (1. - 2.*dTdz);
					Tks   = Tkend * (1. - dTdz);
					alp = (1.-zb*dzb/2./Tks)/(1.+zb*dzb/2./Tks);
				}
				d->Vc[VZ][k][j][i] *= -1.0;
				d->Vc[BX][k][j][i] *= -1.0;
				d->Vc[BY][k][j][i] *= -1.0;
				d->Vc[PRS][k][j][i] = d->Vc[PRS][KEND][j][i]*alp;
				d->Vc[RHO][k][j][i] = d->Vc[PRS][k][j][i]/Tk1;
			}
		}else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
			BOX_LOOP(box,k,j,i) {
				d->Vs[BXs][k][j][i] = -d->Vs[BXs][2*KEND-k+1][j][i];
			}
       #endif
		}else if (box->vpos == X2FACE){
	   #ifdef STAGGERED_MHD
			BOX_LOOP(box,k,j,i) {
				d->Vs[BYs][k][j][i] = -d->Vs[BYs][2*KEND-k+1][j][i];
			}
       #endif
		}
	}
	
}

/* ****************************************************************** */
void BodyForceVector (double *v, 
                 double *g,
                 double x1, double x2, double x3)
/*
 *
 *  PURPOSE
 *
 *   Include body force terms. These are:
 *
 *    - gravitational acceleration vector
 *    - gravitational potential
 *    - Coriolis vector
 *
 *  
 ******************************************************************* */
{
  g[IDIR] =  3.0*x1;
  g[JDIR] =  0.0;
  g[KDIR] = -x3;

  #ifdef CTU
  if (x3 > g_domEnd[2] ) g[KDIR] += 2.0*g_domEnd[2];
  if (x3 < g_domBeg[2] ) g[KDIR] += 2.0*g_domBeg[2];
  #endif
}

