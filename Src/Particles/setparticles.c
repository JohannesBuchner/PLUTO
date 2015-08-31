#include "pluto.h"


extern int nop;

void INTERPOL_VELOCITY (real *v_interpol,struct PARTICLES *pl,
                        real ***uu[], struct GRID gxyz[]);


void INTERPOL (real *a_interpol,struct PARTICLES *pl,
                        real ***uu[], struct GRID gxyz[], int VAR);

void DIV ( real *diver_a, int VAR, real ***uu[], struct GRID gxyz[],
                                       int indici[] );

void LOCATE_PART ( struct PARTICLES *pl, struct GRID gxyz[]); 

void SAVE_STEP ( HEAD *list, real dt, int nop );

/* We define a linked list for PARTICLES */

void ADD_CELL(HEAD *list,struct PARTICLES *pl){

    CELL  *new_cell;

    /* we always add cells at the head: it's easy and cheap */
    /* it means that they are sorted in the exact opposite order of input */
    
    new_cell =(CELL *)malloc (sizeof(CELL));
    new_cell->part= *pl ;
    new_cell->next=list->first;
    list->first=new_cell;

}

/* At end this part should be parametrizable (in order to choose
   the type of initialization in space) and should stand in init.c         */


void SET_PARTICLES (HEAD *list, HEAD *temp_list, struct GRID gxyz[],
                    real ***uu[] )
{
    
    int       i , j , k ; 
    int coords[DIMENSIONS];
    /* variables for interpolation */

    struct PARTICLES pl, pl_temp;
    struct PARTICLES *p_l;
    real *v_interpol,*lx;
    real *b_interpol,*pr_interpol,*rho_interpol;   /* bocchi */
    real *div_v;
    real dxA,dyA;  /* variables for test1 */
    int indici[3];

#ifdef PARALLEL
    int TAG=100;
    int rank_part;
    int cart_comm;
    MPI_Datatype type_PARTICLES;
    MPI_Datatype *type_ptr;
#endif

    lx=array_1D(3);
    v_interpol = array_1D(3);
    b_interpol = array_1D(3);  /* bocchi */
    pr_interpol = array_1D(1);
    rho_interpol = array_1D(1);
    div_v = array_1D(1);

    lx[0]=1./gxyz[0].dx[0];
    lx[1]=1./gxyz[1].dx[0];
    lx[2]=1./gxyz[2].dx[0];
    
    print1("   > Particles intialization \n");

    dxA = (0.8*(gxyz[0].xf-gxyz[0].xi))/(real)NPARTICLES;
    dyA = (gxyz[1].xf-gxyz[1].xi)/3.;

    for(j=0;j<NPARTICLES;++j){

#ifdef PARALLEL

        /*initialize the cartesian communicator of ArrayLib */
        AL_Get_cart_comm(SZ, &cart_comm);

        if(prank==0){
            type_ptr=&(type_PARTICLES);
            pl.identity=j;
            /* the random sequence is always the SAME since the seed is unchanged */
            for(i=0;i<DIMENSIONS;++i){
                pl.coor[i]=gxyz[i].xi+ (gxyz[i].xf-gxyz[i].xi)*(rand()/(1.+RAND_MAX) ) ;
                pl.cell[i]=(int)( ( (pl.coor[i]) - gxyz[i].xi)*lx[i] ) + gxyz[i].nghost;
                coords[i]=pl.cell[i];
                p_l=&pl;
            }

            BUILD_PARTICLES_TYPE(&pl, type_ptr);
            MPI_Bcast( &pl, 1, type_PARTICLES, 0, MPI_COMM_WORLD );

            /* find the processor associated to particle position */
            MPI_Cart_rank(cart_comm, coords, &rank_part);       

            /* send the particle on its processor: rank_part */
            MPI_SSend( pl,1,type_PARTICLES , rank_part, TAG, cart_comm  );

        }
        /* receive the particle on the processor rank_part */
        if(prank==rank_part){
            MPI_RECV( p_l,1,type_PARTICLES , 0, TAG,cart_comm  );
            INTERPOL ( v_interpol, &pl, uu,  gxyz, VX);
            for(i=0;i<DIMENSIONS;++i){pl.speed[i] = v_interpol[i] ;}
            ADD_CELL(list, &pl);
        }


#else 

        pl.identity=j;

        for(i=0;i<DIMENSIONS;++i){
          /*        pl.coor[i]=gxyz[i].xi+ (gxyz[i].xf-gxyz[i].xi)*(rand()/(1.+RAND_MAX) ) ;*/
           /* pl.coor[i]=-1.+1.*(rand()/(1.+RAND_MAX) ) ;*/
            /*      pl.cell[i]=(int)( ( (pl.coor[i]) - gxyz[i].xi)*lx[i] ) + gxyz[i].nghost;        */
        }

          /* bocchi - special initial conditions */
          /* all particles on the jet : r=1.0    */
        /*        pl.coor[0]=gxyz[0].xi + (rand()/(1.+RAND_MAX));
                  pl.coor[1]=gxyz[1].xi+ (gxyz[1].xf-gxyz[1].xi)*(rand()/(1.+RAND_MAX) ) ;*/

          /* bocchi - initial conditions for test1 */

        if(j<(NPARTICLES/2)){
        
          pl.coor[0]=gxyz[0].xi + 0.1*(gxyz[0].xf-gxyz[0].xi) + j*dxA;
          pl.coor[1]=gxyz[1].xi + dyA;

        }else{

          pl.coor[0]=gxyz[0].xi + 0.1*(gxyz[0].xf-gxyz[0].xi) + 
            (j-NPARTICLES/2)*dxA;
          pl.coor[1]=gxyz[1].xi + 2.*dyA;

        }

        /* bocchi */

        /*      printf("\nix_p:%d iy_p:%d\n", pl.cell[0], pl.cell[1]);*/
        LOCATE_PART(&pl, gxyz);
        /*      printf("\nix_d:%d iy_d:%d\n", pl.cell[0], pl.cell[1]);*/

        INTERPOL ( v_interpol, &pl, uu,  gxyz, VX);    
        INTERPOL ( b_interpol, &pl, uu,  gxyz, BX);   
        INTERPOL ( pr_interpol, &pl, uu, gxyz, PR);
        INTERPOL ( rho_interpol,&pl, uu, gxyz, DN);     

      /*        printf("\nic_prima:%d\n", pl.cell[0]);  */

        for(i=0;i<DIMENSIONS;++i){
        indici[i] = pl.cell[i];
        }

        DIV      ( div_v, VX, uu, gxyz, indici );

        for(i=0;i<DIMENSIONS;++i){
          pl.speed[i] = v_interpol[i];
          pl.mag[i] = b_interpol[i]; 
        }    
        pl.pression = *pr_interpol;                            
        pl.density = *rho_interpol;
        pl.divv = *div_v;

        /* bocchi - end */

        ADD_CELL(list, &pl);    
        ADD_CELL(temp_list, &pl_temp);

    }

    nop=NPARTICLES;

    SAVE_STEP(list, 0.0, nop);

#endif

    free_array_1D(v_interpol);
    free_array_1D(b_interpol);   /* bocchi */
    free_array_1D(pr_interpol);
    free_array_1D(rho_interpol);
    free_array_1D(div_v);
    free_array_1D(lx);


}



void ADVANCE_PARTICLES_predictor( HEAD *list, HEAD *temp_list, real ***uu[], 
                                  struct GRID gxyz[], real dt, real striction){

    int i ;

    /* variable for interpolation*/
    real *v_interpol;

    /* bocchi - RK2 - variable for interpolation */
    /*    real *v_int_old;
          struct PARTICLES *pl_dummy; */
    /* bocchi - end */
  
    /* variables for RK2*/
    real *a,*v1,*x1;
    real *lx;
    real *gravity;

    /* Variables to handle data structure*/
    struct PARTICLES *pl,*pl_temp;
    CELL  *actual_cell, *temp_actual_cell;    
    CELL  *previous, *deleted;  /* bocchi */ 
    CELL  *temp_previous, *temp_deleted;  

    gravity= array_1D(3); 
    v_interpol = array_1D(3);
    /*    v_int_old = array_1D(3); */ /* bocchi */
    a  = array_1D(3);
    v1 = array_1D(3);
    x1 = array_1D(3);
    lx = array_1D(3);  
    lx[0]=1./gxyz[0].dx[0];
    lx[1]=1./gxyz[1].dx[0];
    lx[2]=1./gxyz[2].dx[0];
  
    if(list->first != NULL ){
    
        actual_cell=list->first;
        temp_actual_cell=temp_list->first;

        previous=actual_cell;
        temp_previous=temp_actual_cell;

        while(actual_cell !=NULL){
      
            pl=&(actual_cell->part);
            pl_temp=&(temp_actual_cell->part);
            INTERPOL ( v_interpol, pl, uu,  gxyz, VX);
/*
                  GRAVITY_FORCE (pl->coor, gravity);
*/		  
      
            /*-----------------------------------*/
            /*                                   */
            /*   RK 2 scheme to move particles   */
            /*                                   */
            /*-----------------------------------*/     
      
            for(i=0;i<2;++i){
              /*                a[i]  = - striction*( pl->speed[i] - v_interpol[i] )+ gravity[i];       
                v1[i] = pl->speed[i] + dt*a[i];
                x1[i] = pl->coor[i]  + dt*pl->speed[i];*/   /*rem by bocchi*/

              /* bocchi - euler method */
              /*                v1[i] = v_interpol[i];
                                x1[i] = pl->coor[i]  + dt*pl->speed[i];  */
              /* bocchi - end */

              /* bocchi - RK2 */
              /*                v1[i] = v_interpol[i];
                                x1[i] = pl->coor[i] + 0.5*dt*v_interpol[i]; */
              /* bocchi - end */

                /* bocchi - RK2 - alternative */
                v1[i] = v_interpol[i];
                x1[i] = pl->coor[i] + dt*pl->speed[i];
                /*bocchi -RK2 - alternative - end */

            }


            /*---------------------------------------------------*/
            /*                                                   */
            /*          bocchi - Boundary Conditions             */
            /*        xleft: mirror - xright: outflow            */
            /*  yleft(bottom): periodic - yright(top): periodic  */
            /*                                                   */ 
            /*---------------------------------------------------*/

            /* check x outflow first */
                    
            if(x1[0]>gxyz[0].xf){
              /*  DESTROY ACTUAL PARTICLE  */

              deleted=actual_cell;
              temp_deleted=temp_actual_cell;

              previous->next=actual_cell->next;
              temp_previous->next=temp_actual_cell->next;

              if(actual_cell==list->first  ){
                list->first=list->first->next;
                temp_list->first=temp_list->first->next;
              }

              if(actual_cell->next == NULL && actual_cell==list->first ){
                actual_cell=NULL;
                temp_actual_cell=NULL;
              }else{
                actual_cell=actual_cell->next;
                temp_actual_cell=temp_actual_cell->next;
                   }

              free(deleted);
              free(temp_deleted);

              --nop; /* correct the number of particles */

              continue;  /* jump to the end of the loop */
            }

            else{

              if(x1[0]<gxyz[0].xi){ 
                 x1[0] = 2*gxyz[0].xi - x1[0];
                 v1[0] = -v1[0];                
              }
                   
            
              /* check y axis */
            
              if(x1[1]>gxyz[1].xf){ x1[1] = gxyz[1].xi + x1[1] - gxyz[1].xf; }

              if(x1[1]<gxyz[1].xi){ x1[1] = gxyz[1].xf + x1[1] - gxyz[1].xi; }

            }
            
            /*-----------------------*/  
            /* bocchi - Boundary end */
            /*-----------------------*/



            for(i=0;i<DIMENSIONS;++i)
            {
                pl_temp->coor[i]=x1[i] ;
                /*              pl_temp.cell[i]=(int)( (x1[i]-gxyz[i].xi)*lx[i] ) + gxyz[i].nghost;*/
                pl_temp->speed[i]=v1[i];
            }

            LOCATE_PART (pl_temp, gxyz);  /* bocchi */

            /* bocchi - RK2 - fluid velocity in x1 at time t */
            /*
            pl_dummy = &pl_temp;
            INTERPOL ( v_int_old, pl_dummy, uu, gxyz, VX);
            for(i=0;i<DIMENSIONS;++i)
            {
              pl_temp.v_old[i]=v_int_old[i];
            }
            */
            /* bocchi - end */





# if( GEOMETRY==SPHERICAL  )

# if( DIMENSIONS>1)
            /* special conditions for theta-boundary crossing*/
            if( x1[1]>gxyz[1].xf ){ x1[1] = x1[1]-gxyz[1].xf; }
            if( x1[1]<gxyz[1].xi ){ x1[1] = gxyz[1].xi-x1[1]; }

            pl_temp->coor[1]=x1[1];
            pl_temp->cell[1]=(int)( (x1[1]-gxyz[1].xi)*lx[1] ) + gxyz[1].nghost;
            pl_temp->speed[1]=v1[1] ;

#endif
# if( DIMENSIONS>2)
            /* special conditions for phi-boundary  crossing*/
            if( x1[2]>gxyz[2].xf ){ x1[2] = x1[2]-gxyz[2].xf; }
            if( x1[2]<gxyz[2].xi ){ x1[2] = gxyz[2].xi-x1[2]; }
            pl_temp->coor[2]=x1[2];
            pl_temp->cell[2]=(int)( (x1[2]-gxyz[2].xi)*lx[2] ) + gxyz[2].nghost;
            pl_temp->speed[2]=v1[2] ;
#endif
#endif

# if( GEOMETRY==POLAR )
# if( DIMENSIONS>1)
            /* special conditions for theta-boundary crossing*/
            if( x1[1]>gxyz[1].xf ){ x1[1] = x1[1]-gxyz[1].xf; }
            if( x1[1]<gxyz[1].xi ){ x1[1] = gxyz[1].xi-x1[1]; }

            pl_temp->coor[1]=x1[1];
            pl_temp->cell[1]=(int)( (x1[1]-gxyz[1].xi)*lx[1] ) + gxyz[1].nghost;
            pl_temp->speed[1]=v1[1] ;

                    
#endif
#endif
# if( GEOMETRY== CYLINDRICAL )
# if( DIMENSIONS==3)
            /* special conditions for theta-boundary crossing*/
            if( x1[2]>gxyz[2].xf ){ x1[2] = x1[2]-gxyz[2].xf; }
            if( x1[2]<gxyz[2].xi ){ x1[2] = gxyz[2].xi-x1[2]; }
                    
            pl_temp->coor[2]=x1[2] ;
            pl_temp->cell[2]=(int)( (x1[2]-gxyz[2].xi)*lx[2] ) + gxyz[2].nghost;
            pl_temp->speed[2]=v1[2];
#endif
#endif
/*
             ADD_CELL(temp_list, &pl_temp);
*/	     
            previous=actual_cell;  
            actual_cell=actual_cell->next;

            temp_previous=temp_actual_cell;  
            temp_actual_cell=temp_actual_cell->next;
        
        }
    }

    free_array_1D(gravity);
    free_array_1D(a);
    free_array_1D(v1);
    free_array_1D(x1);
    free_array_1D(lx);
    free_array_1D(v_interpol);
    /*    free_array_1D(v_int_old); */   /* bocchi */

}



void ADVANCE_PARTICLES_corrector( HEAD *list, HEAD *temp_list, real ***uu[], struct GRID gxyz[]
                                  , real dt, real striction)
{
    int i ;

    /* variable for interpolation*/
    real *v_interpol;
    real *b_interpol,*rho_interpol,*pr_interpol;  /* bocchi */
    real *div_v;
    int indici[3];

    /* variables for RK2*/
    real *a,*v,*x;
    real *lx;
    real *gravity;

    /* Variables to handle data structure*/
    struct PARTICLES *pl,*pl_temp;

    CELL  *actual_cell, *previous, *deleted;    
    CELL  *temp_actual_cell, *temp_previous, *temp_deleted; 
    
    gravity= array_1D(3); 
    v_interpol = array_1D(3);
    b_interpol = array_1D(3);   /* bocchi */
    rho_interpol = array_1D(1);
    pr_interpol = array_1D(1);
    div_v = array_1D(1);
    a  = array_1D(3);
    v  = array_1D(3);
    x  = array_1D(3);
    lx = array_1D(3);

    lx[0]=1./gxyz[0].dx[0];
    lx[1]=1./gxyz[1].dx[0];
    lx[2]=1./gxyz[2].dx[0];

    nop = 0; /* bocchi */
    
    if(list->first != NULL ){

        actual_cell=list->first;
        temp_actual_cell=temp_list->first;
        
        previous=actual_cell;
        temp_previous=temp_actual_cell;

        while(actual_cell !=NULL){
            pl=&(actual_cell->part);
            pl_temp=&(temp_actual_cell->part);
            
            ++nop; /* bocchi */

            INTERPOL ( v_interpol, pl_temp, uu,  gxyz, VX);

            /*      print1("A u1:%f,u2:%f\n", v_interpol[0],v_interpol[1]);*/
/*
                  GRAVITY_FORCE (pl_temp->coor, gravity);
*/		  

            /*-----------------------------------*/
            /*                                   */
            /*   RK 2 scheme to move particles   */
            /*                                   */
            /*-----------------------------------*/     

            for(i=0;i<2;++i){    
              /*                a[i]  = - striction*( pl_temp->speed[i] - v_interpol[i] )+ gravity[i];
                v[i]  = 0.5*( pl_temp->speed[i] + pl->speed[i] + dt*a[i] );
                x[i]  = pl_temp->coor[i]  + dt*0.5*(pl_temp->speed[i] + pl->speed[i] );*/    /*rem by bocchi*/ 

              /* bocchi - euler method */
              /*                v[i] = v_interpol[i];
                                x[i] = pl_temp->coor[i];*/
              /* bocchi - end */

              /* bocchi - RK2 */
              /*                x[i] = pl->coor[i] + dt*0.5*( v_interpol[i] + pl_temp->v_old[i] );  
                                v[i] = v_interpol[i]; */   /* dummy value */
              /* bocchi - RK2 - end */

                /* bocchi - RK2 - alternative */
                v[i] = v_interpol[i]; /* dummy value */
                x[i] = pl->coor[i] + dt*0.5*(pl->speed[i] + v_interpol[i]);
                /* bocchi - RK2 - alternative - end */

            }

            /*
                if(pl->identity==10){
                  printf("\npart10:%lf,%lf",pl->speed[1],v_interpol[1]);
                }
                if(pl->identity==100){
                  printf("\npart100:%lf,%lf",pl->speed[1],v_interpol[1]);
                }
            */



            /*---------------------------------------------------*/
            /*                                                   */
            /*          bocchi - Boundary Conditions             */
            /*        xleft: mirror - xright: outflow            */
            /*  yleft(bottom): periodic - yright(top): periodic  */
            /*                                                   */ 
            /*---------------------------------------------------*/

            /* check x outflow first */
                    
            if(x[0]>gxyz[0].xf){
              /*  DESTROY ACTUAL PARTICLE  */

              deleted=actual_cell;
              temp_deleted=temp_actual_cell;

              previous->next=actual_cell->next;
              temp_previous->next=temp_actual_cell->next;

              if(actual_cell==list->first  ){
                list->first=list->first->next;
                temp_list->first=temp_list->first->next;
              }

              if(actual_cell->next == NULL && actual_cell==list->first ){
                actual_cell=NULL;
                temp_actual_cell=NULL;
              }else{
                actual_cell=actual_cell->next;
                temp_actual_cell=temp_actual_cell->next;
                   }

              free(deleted);
              free(temp_deleted);

              --nop; /* correct the number of particles */

              continue;  /* jump to the end of the loop */
            }

            else{

              if(x[0]<gxyz[0].xi){ 
                 x[0] = 2*gxyz[0].xi - x[0];
                 /* speed correction not needed: calculated after this part */                
              }
                   
            
              /* check y axis */
            
              if(x[1]>gxyz[1].xf){ x[1] = gxyz[1].xi + x[1] - gxyz[1].xf; }

              if(x[1]<gxyz[1].xi){ x[1] = gxyz[1].xf + x[1] - gxyz[1].xi; }

            }
            
            /*-----------------------*/  
            /* bocchi - Boundary end */
            /*-----------------------*/




#if( GEOMETRY==CARTESIAN)
# if( DIMENSIONS==1)
            if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf )
#endif
# if( DIMENSIONS==2)
                if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf &&
                   x[1]>gxyz[1].xi && x[1]<gxyz[1].xf )
#endif
# if( DIMENSIONS==3)
                    if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf &&
                       x[1]>gxyz[1].xi && x[1]<gxyz[1].xf &&
                       x[2]>gxyz[2].xi && x[2]<gxyz[2].xf)
#endif
                    {
                        for(i=0;i<DIMENSIONS;++i)
                        {
                            pl->coor[i]=x[i] ;
                            pl->cell[i]=(int)( (x[i]-gxyz[i].xi)*lx[i] ) + gxyz[i].nghost;
                            pl->speed[i]=v[i];
                        }




                        /*print1("   >          particle : %d, x1:%f, x2:%f \n",
                               actual_cell->part.identity,
                               actual_cell->part.coor[0],
                               actual_cell->part.coor[1]   );*/


                        previous=actual_cell;
                        actual_cell=actual_cell->next;
                        
                        temp_previous=temp_actual_cell;
                        temp_actual_cell=temp_actual_cell->next;
                    }
                    
#endif


# if( GEOMETRY==SPHERICAL  )

            if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf){
                pl->coor[0]=x[0] ;
                pl->cell[0]=(int)( (x[0]-gxyz[0].xi)*lx[0] ) + gxyz[0].nghost;
                pl->speed[0]=v[0];
# if( DIMENSIONS>1)
                /* special conditions for theta-boundary crossing*/
                if( x[1]>gxyz[1].xf ){ x[1] = x[1]-gxyz[1].xf; }
                if( x[1]<gxyz[1].xi ){ x[1] = gxyz[1].xi-x[1]; }
                pl->coor[1]=x[1];
                pl->cell[1]=(int)( (x[1]-gxyz[1].xi)*lx[1] ) + gxyz[1].nghost;
                pl->speed[1]=v[1] ;

#endif
# if( DIMENSIONS>2)
                /* special conditions for phi-boundary  crossing*/
                if( x[2]>gxyz[2].xf ){ x[2] = x[2]-gxyz[2].xf; }
                if( x[2]<gxyz[2].xi ){ x[2] = gxyz[2].xi-x[2]; }
                pl->coor[2]=x[2];
                pl->cell[2]=(int)( (x[2]-gxyz[2].xi)*lx[2] ) + gxyz[2].nghost;
                pl->speed[2]=v[2] ;
#endif
                previous=actual_cell;
                actual_cell=actual_cell->next;
                temp_previous=temp_actual_cell;
                temp_actual_cell=temp_actual_cell->next;
            }       

#endif


# if( GEOMETRY==POLAR )

# if( DIMENSIONS<3)
            if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf)
#endif

# if( DIMENSIONS==3)
                if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf &&
                   x[2]>gxyz[2].xi && x[2]<gxyz[2].xf  )
#endif
                {
                    pl->coor[0]=x[0] ;
                    pl->cell[0]=(int)( (x[0]-gxyz[0].xi)*lx[0] ) + gxyz[0].nghost;
                    pl->speed[0]=v[0];
                    
# if( DIMENSIONS==3)
                    pl->coor[2]=x[2] ;
                    pl->cell[2]=(int)( (x[2]-gxyz[2].xi)*lx[2] ) + gxyz[2].nghost;
                    pl->speed[2]=v[2];
#endif

# if( DIMENSIONS>1)
                    /* special conditions for theta-boundary crossing*/
                    if( x[1]>gxyz[1].xf ){ x[1] = x[1]-gxyz[1].xf; }
                    if( x[1]<gxyz[1].xi ){ x[1] = gxyz[1].xi-x[1]; }

                    pl->coor[1]=x[1];
                    pl->cell[1]=(int)( (x[1]-gxyz[1].xi)*lx[1] ) + gxyz[1].nghost;
                    pl->speed[1]=v[1] ;

                    /*print1("   >          particle : %d , r:%f, theta:%f \n",
                           actual_cell->part.identity,
                           actual_cell->part.coor[0],
                           57.29578*actual_cell->part.coor[1]   );*/

#endif
                    previous=actual_cell;
                    actual_cell=actual_cell->next;

                    temp_previous=temp_actual_cell;
                    temp_actual_cell=temp_actual_cell->next;
                }
#endif


# if( GEOMETRY==CYLINDRICAL )

            if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf)

# if( DIMENSIONS==3)
                if(x[0]>gxyz[0].xi && x[0]<gxyz[0].xf &&
                   x[1]>gxyz[1].xi && x[1]<gxyz[1].xf  )
#endif
                {
                    pl->coor[0]=x[0] ;
                    /*              pl->cell[0]=(int)( (x[0]-gxyz[0].xi)*lx[0] ) + gxyz[0].nghost;*/
                    pl->speed[0]=v[0];

# if( DIMENSIONS>1)
                    pl->coor[1]=x[1] ;
                    /*              pl->cell[1]=(int)( (x[1]-gxyz[1].xi)*lx[1] ) + gxyz[1].nghost;*/
                    pl->speed[1]=v[1];
#endif
# if( DIMENSIONS>2)
                    /* special conditions for theta-boundary crossing*/
                    if( x[2]>gxyz[2].xf ){ x[2] = x[2]-gxyz[2].xf; }
                    if( x[2]<gxyz[2].xi ){ x[2] = gxyz[2].xi-x[2]; }
                    pl->coor[2]=x[2];
                    /*              pl->cell[2]=(int)( (x[2]-gxyz[2].xi)*lx[2] ) + gxyz[2].nghost;*/
                    pl->speed[2]=v[2] ;

#endif

                        /* bocchi - RK2 - adjust velocity */

                        LOCATE_PART ( pl, gxyz );

                        INTERPOL ( v_interpol, pl, uu, gxyz, VX);
                        INTERPOL ( b_interpol, pl, uu, gxyz, BX);
                        INTERPOL ( rho_interpol,pl,uu, gxyz, DN);
                        INTERPOL ( pr_interpol, pl,uu, gxyz, PR);

                        for(i=0;i<DIMENSIONS;++i){
                          indici[i] = pl->cell[i];
                        }

                        DIV ( div_v, VX, uu, gxyz, indici );

                        for(i=0;i<DIMENSIONS;++i)
                        {
                          pl->speed[i]=v_interpol[i];
                          pl->mag[i]=b_interpol[i];
                        }

                        pl->pression=*pr_interpol;
                        pl->density=*rho_interpol;
                        pl->divv = *div_v;

                
                        /* bocchi - RK2 - end */
        
                /*                      print1("B u1:%f,u2:%f\n", v_interpol[0],v_interpol[1]);*/
                        /*                      print1("\n part:%d, v1:%f, v2:%f, b2:%lf, divv:%f\n",
                               actual_cell->part.identity,
                               actual_cell->part.speed[0],
                               actual_cell->part.speed[1],
                               actual_cell->part.mag[1],
                               actual_cell->part.divv); */ /* bocchi */



                        /*                          print1("   >          particule: %d, rayon: %f\n",
                           actual_cell->part.identity,
                           actual_cell->part.coor[0]  );*/
                    
                    previous=actual_cell;
                    actual_cell=actual_cell->next;
                    
                    temp_previous=temp_actual_cell;
                    temp_actual_cell=temp_actual_cell->next;
                }

#endif
                else{
                    /* enter here if only particle is out the domain*/
                    /*print1("   > Canceled particle : %d, r: %f, theta:%f\n",
                      actual_cell->part.identity, 
                      x[0],
                      actual_cell->part.coor[1]);*/
                    
                    deleted=actual_cell;
                    temp_deleted=temp_actual_cell;

                    previous->next=actual_cell->next;
                    temp_previous->next=temp_actual_cell->next;

                    if(actual_cell==list->first  ){
                        list->first=list->first->next;
                        temp_list->first=temp_list->first->next;
                    }

                    if(actual_cell->next == NULL && actual_cell==list->first ){
                        actual_cell=NULL;
                        temp_actual_cell=NULL;
                    }else{
                        actual_cell=actual_cell->next;
                        temp_actual_cell=temp_actual_cell->next;
                    }

                    free(deleted);
                    free(temp_deleted);
                }
        }
    }


    /* last thing to do: destruct all cells of the temporary list*/
    /*    if(temp_list->first != NULL ){

        temp_actual_cell=temp_list->first;
        temp_previous=temp_actual_cell;

        while(temp_list->first !=NULL){
            deleted=temp_list->first;
            temp_list->first=temp_list->first->next;
            free(deleted);
        }

        }*/     /* now the temp_list is mem-allocked with the list so 
                   no need to destroy it now */

    free_array_1D(gravity);
    free_array_1D(a);
    free_array_1D(v);
    free_array_1D(x);
    free_array_1D(lx);
    free_array_1D(v_interpol);
    free_array_1D(b_interpol);  /* bocchi */
    free_array_1D(rho_interpol);
    free_array_1D(pr_interpol);
    free_array_1D(div_v);


}
void SAVE_PARTICLES( HEAD *list , int nt, real t, real striction,
                     real ***uu[] ){

    FILE      *pout ,*verif;
    int j,i0,i1,count=0 ;
    real u0, u1, v0, v1, r,theta;
    struct PARTICLES *pl;
    CELL  *actual_cell;
    
    if(list->first!= NULL){
        actual_cell=list->first;        
        if(nt==0){
            print("WRITING Particles_0\n" );
            pout = fopen("0.out", "w");
        }
        if(nt==50){
            print("WRITING Particles_1\n" );
            pout = fopen("100.out", "w");
        }
        if(nt==100){
            print("WRITING Particles_2\n" );
            pout = fopen("200.out", "w");
        }
        if(nt==150){
            print("WRITING Particles_3\n" );
            pout = fopen("300.out", "w");
        }
        while(actual_cell!= NULL) {

            count+=1;   
            pl=&(actual_cell->part);
        
            /*r     = pl->coor[0];
              theta = pl->coor[1];*/

            fprintf(pout," %f  %f",pl->coor[0],pl->coor[1] );
            fprintf(pout,"\n" );
            actual_cell=actual_cell->next;
        }   

        print(">number of particles written:%d \n",count );
    
        fclose(pout);
    
    }
}

void INTERPOL_VELOCITY (real *v_interpol,struct PARTICLES *pl, real ***uu[]
                        , struct GRID gxyz[]){


    /* variables for interpolation*/
    int i,il, ir, jr, jl ,kr, kl, ind0, ind1, ind2;
    real  ix, iy, iz;
    
    ind0=pl->cell[0];
    ind1=pl->cell[1];

    /*-------------------------------------------------------*/
    /*                                                       */
    /*   This Perform Lagrangian interpolation of velocity   */
    /*            base of polynoms: < 1, x, y, xy >          */
    /*                                                       */           
    /*-------------------------------------------------------*/

    if(gxyz[0].x[ind0] >pl->coor[0]){
        il = ind0-1;
        ir = ind0;
    }else{
        il = ind0;
        ir = ind0+1;
    }
    ix = ( pl->coor[0] - gxyz[0].x[il] )/( gxyz[0].x[ir] - gxyz[0].x[il] );    

#if ( DIMENSIONS == 2 ||DIMENSIONS ==3 ) 

    if(gxyz[1].x[ind1] >pl->coor[1]){
        jl = ind1-1;
        jr = ind1 ;
    }else{
        jl = ind1;
        jr = ind1+1 ;
    }
    iy = ( pl->coor[1] - gxyz[1].x[jl] )/( gxyz[1].x[jr] - gxyz[1].x[jl] );
        
#endif

#if DIMENSIONS ==3
    if(gxyz[2].x[ind2] >pl->coor[2]){
        kl = ind2-1;
        kr = ind2 ;
    }else{
        kl = ind2;
        kr = ind2+1 ;
    }
    iz = ( pl->coor[2] - gxyz[2].x[kl] )/( gxyz[2].x[kr] - gxyz[2].x[kl] );
#endif

#if DIMENSIONS == 1     
    v_interpol[0] = uu[VX1][0][0][il]
        + ix*(uu[VX1][0][0][ir] - uu[VX1][0][0][il] ) ;
#endif

#if DIMENSIONS == 2     

    for(i=0;i<DIMENSIONS;++i){  
        v_interpol[i] = 
            (1.-ix)*( uu[VX+i][0][jl][il]+iy*(uu[VX+i][0][jr][il]-uu[VX+i][0][jl][il]))
            + ix*(uu[VX+i][0][jl][ir]+iy*(uu[VX+i][0][jr][ir] -uu[VX+i][0][jl][ir]));
    }

#endif

#if DIMENSIONS == 3

    /* We use a formulation which asks for 18 calls to variables and 8 products*/
    for(i=0;i<DIMENSIONS;++i){
        v_interpol[i] =
            (1.-ix)*( uu[VX+i][kl][jl][il]
                      + iz*(uu[VX+i][kr][jl][il] - uu[VX+i][kl][jl][il] )
                      + iy *( (uu[VX+i][kl][jr][il] - uu[VX+i][kl][jl][il] )
                              +iz*(  uu[VX+i][kr][jr][il] - uu[VX+i][kl][jr][il]
                                     - uu[VX+i][kr][jl][il] + uu[VX+i][kl][jl][il])
                          )
                )
            +   ix*( uu[VX+i][kl][jl][ir]
                     + iz*(uu[VX+i][kr][jl][ir] - uu[VX+i][kl][jl][ir] )
                     + iy *( (uu[VX+i][kl][jr][ir] - uu[VX+i][kl][jl][ir] )
                             +iz*(    uu[VX+i][kr][jr][ir] - uu[VX+i][kl][jr][ir]
                                      - uu[VX+i][kr][jl][ir] + uu[VX+i][kl][jl][ir])
                         )
                );

    }
#endif

}




void INTERPOL (real *a_interpol,struct PARTICLES *pl, real ***uu[]
                        , struct GRID gxyz[], int VAR){


    /* variables for interpolation*/
    int i,il, ir, jr, jl ,kr, kl, ind0, ind1, ind2;
    real  ix, iy, iz;
    
    ind0=pl->cell[0];
    ind1=pl->cell[1];


    /*-------------------------------------------------------*/
    /*                                                       */
    /*   This Perform Lagrangian interpolation of velocity   */
    /*            base of polynoms: < 1, x, y, xy >          */
    /*                                                       */           
    /*-------------------------------------------------------*/

    if(gxyz[0].x[ind0] >pl->coor[0]){
        il = ind0-1;
        ir = ind0;
    }else{
        il = ind0;
        ir = ind0+1;
    }
    ix = ( pl->coor[0] - gxyz[0].x[il] )/( gxyz[0].x[ir] - gxyz[0].x[il] );    

#if ( DIMENSIONS == 2 ||DIMENSIONS ==3 ) 

    if(gxyz[1].x[ind1] >pl->coor[1]){
        jl = ind1-1;
        jr = ind1 ;
    }else{
        jl = ind1;
        jr = ind1+1 ;
    }
    iy = ( pl->coor[1] - gxyz[1].x[jl] )/( gxyz[1].x[jr] - gxyz[1].x[jl] );
        
#endif

#if DIMENSIONS ==3
    if(gxyz[2].x[ind2] >pl->coor[2]){
        kl = ind2-1;
        kr = ind2 ;
    }else{
        kl = ind2;
        kr = ind2+1 ;
    }
    iz = ( pl->coor[2] - gxyz[2].x[kl] )/( gxyz[2].x[kr] - gxyz[2].x[kl] );
#endif


    if(VAR==VX||VAR==BX){


#if DIMENSIONS == 1     
    a_interpol[0] = uu[var][0][0][il]
        + ix*(uu[VAR][0][0][ir] - uu[VAR][0][0][il] ) ;
#endif

#if DIMENSIONS == 2     

    for(i=0;i<DIMENSIONS;++i){  
        a_interpol[i] = 
            (1.-ix)*( uu[VAR+i][0][jl][il]+iy*(uu[VAR+i][0][jr][il]-uu[VAR+i][0][jl][il]))
            + ix*(uu[VAR+i][0][jl][ir]+iy*(uu[VAR+i][0][jr][ir] -uu[VAR+i][0][jl][ir]));
    }

#endif

#if DIMENSIONS == 3

    /* We use a formulation which asks for 18 calls to variables and 8 products*/
    for(i=0;i<DIMENSIONS;++i){
        a_interpol[i] =
            (1.-ix)*( uu[VAR+i][kl][jl][il]
                      + iz*(uu[VAR+i][kr][jl][il] - uu[VAR+i][kl][jl][il] )
                      + iy *( (uu[VAR+i][kl][jr][il] - uu[VAR+i][kl][jl][il] )
                              +iz*(  uu[VAR+i][kr][jr][il] - uu[VAR+i][kl][jr][il]
                                     - uu[VAR+i][kr][jl][il] + uu[VAR+i][kl][jl][il])
                          )
                )
            +   ix*( uu[VAR+i][kl][jl][ir]
                     + iz*(uu[VAR+i][kr][jl][ir] - uu[VAR+i][kl][jl][ir] )
                     + iy *( (uu[VAR+i][kl][jr][ir] - uu[VAR+i][kl][jl][ir] )
                             +iz*(    uu[VAR+i][kr][jr][ir] - uu[VAR+i][kl][jr][ir]
                                      - uu[VAR+i][kr][jl][ir] + uu[VAR+i][kl][jl][ir])
                         )
                );

    }
#endif
    }
    else{


#if DIMENSIONS == 1     
    *a_interpol = uu[VAR][0][0][il]
        + ix*(uu[VAR][0][0][ir] - uu[VAR][0][0][il] ) ;
#endif

#if DIMENSIONS == 2     

        
        *a_interpol = 
            (1.-ix)*( uu[VAR][0][jl][il]+iy*(uu[VAR][0][jr][il]-uu[VAR][0][jl][il]))
            + ix*(uu[VAR][0][jl][ir]+iy*(uu[VAR][0][jr][ir] -uu[VAR][0][jl][ir]));
    

#endif

#if DIMENSIONS == 3

    /* We use a formulation which asks for 18 calls to variables and 8 products*/
    
        *a_interpol =
            (1.-ix)*( uu[VAR][kl][jl][il]
                      + iz*(uu[VAR][kr][jl][il] - uu[VAR][kl][jl][il] )
                      + iy *( (uu[VAR][kl][jr][il] - uu[VAR][kl][jl][il] )
                              +iz*(  uu[VAR][kr][jr][il] - uu[VAR][kl][jr][il]
                                     - uu[VAR][kr][jl][il] + uu[VAR][kl][jl][il])
                          )
                )
            +   ix*( uu[VAR][kl][jl][ir]
                     + iz*(uu[VAR][kr][jl][ir] - uu[VAR][kl][jl][ir] )
                     + iy *( (uu[VAR][kl][jr][ir] - uu[VAR][kl][jl][ir] )
                             +iz*(    uu[VAR][kr][jr][ir] - uu[VAR][kl][jr][ir]
                                      - uu[VAR][kr][jl][ir] + uu[VAR][kl][jl][ir])
                         )
                );
#endif

    }

}





void RENAME_PARTICLES( HEAD *list ){
    

    /* This routine simply rename particles, starting from the root*/
    
    int count=0;
    struct PARTICLES *pl;
    CELL  *actual_cell;


    if(list->first!= NULL){
        actual_cell=list->first;
    
        while(actual_cell!= NULL) {
        
            pl=&(actual_cell->part);
            pl->identity=count;
            count+=1;
            actual_cell=actual_cell->next;
        }   
    
        print("   > particles re_named: %d \n",count );
    }
}

/*void SEND_PARTICLES(PARTICLES *pl, int proc_ini, int proc_end ){*/

/*   This Function sends a particle:  */
/*   from  processor proc_ini         */
/*   to    processor proc_end         */




/* bocchi */
void DIV ( real *diver, int VAR, real ***uu[], struct GRID gxyz[],
            int indici[] ){

  /* perform divergence of the array_1D VAR: es VX */

  int il, ic, ir, jl, jc, jr, kl, kc, kr;
  real dx1, dx2, dy1, dy2, dz1, dz2;
  real ax, bx, cx, ay, by, cy, az, bz, cz;

#if ( GEOMETRY == CYLINDRICAL )

# if DIMENSIONS == 1 

  ic = indici[0];
  il = ic - 1;
  ir = ic + 1;

  dx1 = gxyz[0].x[ic] - gxyz[0].x[il];
  dx2 = gxyz[0].x[ir] - gxyz[0].x[ic];

  ax = (dx2*dx2-dx1*dx1)/(dx1*dx2*(dx1+dx2));
  bx = -dx2/(dx1*(dx1+dx2));
  cx =  dx1/(dx2*(dx1+dx2));

  *diver = ax*uu[VAR][0][0][ic] + bx*uu[VAR][0][0][il] + cx*uu[VAR][0][0][ir]
           + uu[VAR][0][0][ic]/gxyz[0].x[ic];

#endif


# if DIMENSIONS == 2

  ic = indici[0];
  il = ic - 1;
  ir = ic + 1;

  jc = indici[1];
  jl = jc - 1;
  jr = jc + 1;

  dx1 = gxyz[0].x[ic] - gxyz[0].x[il];
  dx2 = gxyz[0].x[ir] - gxyz[0].x[ic];

  dy1 = gxyz[1].x[jc] - gxyz[1].x[jl];
  dy2 = gxyz[1].x[jr] - gxyz[1].x[jc];


  ax = (dx2*dx2-dx1*dx1)/(dx1*dx2*(dx1+dx2));
  bx = -dx2/(dx1*(dx1+dx2));
  cx =  dx1/(dx2*(dx1+dx2));

  ay = (dy2*dy2-dy1*dy1)/(dy1*dy2*(dy1+dy2));
  by = -dy2/(dy1*(dy1+dy2));
  cy =  dy1/(dy2*(dy1+dy2));

  *diver = ax*uu[VAR][0][jc][ic] + bx*uu[VAR][0][jc][il] 
           + cx*uu[VAR][0][jc][ir]
           + uu[VAR][0][jc][ic]/gxyz[0].x[ic]
           + ay*uu[VAR+1][0][jc][ic] + by*uu[VAR+1][0][jl][ic] 
           + cy*uu[VAR+1][0][jr][ic];

  

#endif

# if DIMENSIONS == 3 

  ic = pl->cell[0];
  il = ic - 1;
  ir = ic + 1;

  jc = pl->cell[1];
  jl = jc - 1;
  jr = jc + 1;

  kc = pl->cell[2];
  kl = kc - 1;
  kr = kc + 1;

  dx1 = gxyz[0].x[ic] - gxyz[0].x[il];
  dx2 = gxyz[0].x[ir] - gxyz[0].x[ic];

  dy1 = gxyz[1].x[jc] - gxyz[1].x[jl];
  dy2 = gxyz[1].x[jr] - gxyz[1].x[jc];

  dz1 = gxyz[2].x[kc] - gxyz[2].x[kl];
  dz2 = gxyz[2].x[kr] - gxyz[2].x[kc];

  ax = (dx2*dx2-dx1*dx1)/(dx1*dx2*(dx1+dx2));
  bx = -dx2/(dx1*(dx1+dx2));
  cx =  dx1/(dx2*(dx1+dx2));

  ay = (dy2*dy2-dy1*dy1)/(dy1*dy2*(dy1+dy2));
  by = -dy2/(dy1*(dy1+dy2));
  cy =  dy1/(dy2*(dy1+dy2));

  az = (dz2*dz2-dz1*dz1)/(dz1*dz2*(dz1+dz2));
  bz = -dz2/(dz1*(dz1+dz2));
  cz =  dz1/(dz2*(dz1+dz2));

  *diver = ax*uu[VAR][kc][jc][ic] + bx*uu[VAR][kc][jc][il] 
           + cx*uu[VAR][kc][jc][ir]
           + uu[VAR][kc][jc][ic]/gxyz[0].x[ic]
           + ay*uu[VAR+1][kc][jc][ic] + by*uu[VAR+1][kc][jl][ic] 
           + cy*uu[VAR+1][kc][jr][ic]
           + (az*uu[VAR+2][kc][jc][ic] + bz*uu[VAR+2][kl][jc][ic]
           + cz*uu[VAR+2][kr][jc][ic])/gxyz[0].x[ic];

#endif  

#endif




}
/* bocchi end */

#ifdef PARALLEL
void BUILD_PARTICLES_TYPE(struct PARTICLES* pl,MPI_Datatype* type_ptr){

    /* This function create a MPI type for the particles to ease     */
    /* SEND and RECV procedures.                                     */ 

    /* For information see "A User Guide to MPI" by Peter S. Pachero */
    /* (his responsible for my mistakes...) page 28-30               */


    int block_lengths[4];
    MPI_Aint displacements[4];
    MPI_Aint addresses[5];
    MPI_Datatype typelist[4];

    /* specify the types */
    typelist[0] = MPI_INT    ;
    typelist[1] = MPI_DOUBLE ;
    typelist[2] = MPI_DOUBLE ;
    typelist[3] = MPI_INT    ;

    /* specify the number of elements of each type */

    block_lengths[0] = 1;
    block_lengths[1] = DIMENSIONS;
    block_lengths[2] = DIMENSIONS;
    block_lengths[3] = DIMENSIONS;

    /* Calculate the displacements of the members */
  
    MPI_Address(pl             , &addresses[0]);
    MPI_Address(&(pl->identity), &addresses[1]);
    MPI_Address(&(pl->coor)    , &addresses[2]);
    MPI_Address(&(pl->speed)   , &addresses[3]);
    MPI_Address(&(pl->cell)    , &addresses[4]);

    displacements[0]=addresses[1]-addresses[0];
    displacements[1]=addresses[2]-addresses[0];
    displacements[2]=addresses[3]-addresses[0];
    displacements[3]=addresses[4]-addresses[0];

    /* Create the derived type */
  
    MPI_Type_struct(4,block_lengths,displacements,typelist,type_ptr);

    /* Commit it so it can be used */
  
    MPI_Type_commit(type_ptr);

}  

#endif



/* bocchi */

void SAVE_STEP ( HEAD *list, real dt, int nop ){

  /* save structure particle, with dt and number of particles */

  FILE *stream;
  struct PARTICLES *pl;
  CELL *actual_cell;

  stream = fopen("part_data.bin.out", "a");

  fwrite(&nop, sizeof(int), 1, stream);
  fwrite(&dt, sizeof(real), 1, stream);

  if (list->first != NULL){
  
    actual_cell = list->first;

    while (actual_cell != NULL){

      pl = &(actual_cell->part);

      fwrite(pl, sizeof(particle), 1, stream);

      actual_cell = actual_cell->next;

    }
  }

  fclose(stream);

}



void LOCATE_PART ( struct PARTICLES *pl, struct GRID gxyz[] ){


  /* determine the cell of the particle, even in stretched grids */
  

  int l_ind, r_ind;
  int m_ind, i, j;
  int true = 1;
  
    for(i=0;i<DIMENSIONS;++i){
   

    l_ind = gxyz[i].nghost;
    r_ind = gxyz[i].np_int_glob + gxyz[i].nghost - 1;


    while ( true ){

      m_ind = l_ind + (r_ind - l_ind)/2;

      if (pl->coor[i] <= gxyz[i].xr[m_ind]){
        r_ind = m_ind;
      }
      else{
        l_ind = m_ind + 1;
      }

      if( l_ind == r_ind ) break;

      }      


  /*now ind_left == ind_right == ind_part */

    pl->cell[i] = r_ind;

  }

}

/* bocchi - end */
