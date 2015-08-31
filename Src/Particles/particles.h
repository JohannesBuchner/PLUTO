/*  **************************************************************
            Header file for particle module.
    **************************************************************  */

#define NPARTICLES 256

typedef struct PARTICLES
{  
    int    identity             ;  /*   identifiant of  particle  */
    int    cell [DIMENSIONS]    ;  /*   particle'   cell index    */
    real coor [DIMENSIONS]    ;  /*   coordinates of particle   */
    real speed[DIMENSIONS]    ;  /*   speed       of particle   */
    real mag  [DIMENSIONS]    ;  /*   magnetic field            */
    real density              ;  /*   density in part.coor      */
    real pression             ;  /*   pression in part.coor     */
    real divv                 ;  /*   speed divergence          */

}particle;

typedef struct head {
    struct cell *first;
}HEAD;

typedef struct cell{
    struct PARTICLES part;
    struct cell *next;
}CELL;

/*  ----  Prototyping goes here  ----  */

void SET_PARTICLES (HEAD *list, HEAD *temp_list, struct GRID gxyz[],
		    real ***uu[] );
void ADVANCE_PARTICLES_predictor (HEAD *list,HEAD *temp_list, real ***uu[], 
                                  struct GRID gxyz[], real dt, real striction);
void ADVANCE_PARTICLES_corrector (HEAD *list,HEAD *temp_list, real ***uu[], 
                                  struct GRID gxyz[], real dt, real striction);
void SAVE_PARTICLES(HEAD *list, int nt, real t, 
                    real striction, real ***uu[]);
void SAVE_STEP ( HEAD *list, real dt, int nop );
void RENAME_PARTICLES( HEAD *list);

/* ROUTINES FOR LINKED LIST */

void ADD_CELL(HEAD *list, struct PARTICLES *pl);
void DELETE_CELL(HEAD *list, int i);
void DISPLAY_LIST(HEAD *list);




